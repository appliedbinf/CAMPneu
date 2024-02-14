#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_dir = "/scicomp/home-pure/ubt4/mycoplasma/test_Raw/*_{1,2}.fq"
// params.input_dir = "/scicomp/home-pure/ubt4/mycoplasma/allRawData/*_{1,2}.fastq"
params.reference_dir = "/scicomp/home-pure/ubt4/mycoplasma/nextflow/updated/references/*.fna"
params.ref23S = "/scicomp/home-pure/ubt4/mycoplasma/nextflow/updated/23S_reference_positions.csv"
params.snp23S = "/scicomp/home-pure/ubt4/mycoplasma/nextflow/updated/23sSNPS.bed"
params.snippyInput = "/scicomp/home-pure/ubt4/mycoplasma/nextflow/updated/inputSnippy.tab"

process assembly {

    publishDir 'assemblies'

    input:
    tuple val(sampleID), path(reads)

    output:
    path("*.fasta"), emit: genomes

    script:
    """
    spades.py -1 ${reads[0]} -2 ${reads[1]} -o ${sampleID} 
    mv ./${sampleID}/contigs.fasta ./${sampleID}.fasta
    """
}

process fastANI{

    publishDir 'fastANI'

    input:
    tuple path(assembly), val(sample), val(ref_label), path(reference)

    output: 
    tuple val(sample), path("${sample}_${reference.baseName}.out"), val(ref_label)

    script:
    """
    fastANI -q ${assembly} -r ${reference} -o ${sample}_${reference.baseName}.out
    """

}

process bestRef {

    publishDir 'bestReference'

    input:
    tuple val(sample), path(ani_res), val(ref_label)

    output:
    tuple val(sample), env(ref)

    script:
    """
    cat ${ani_res} > ${sample}_allRef.txt
    ref=\$(sort -n -k 3 ${sample}_allRef.txt | tail -n 1 | cut -f2 | sed 's/.fna//')
    """
}

process minimap2 {

    publishDir 'minimap2'

    input:
    tuple val(sample), val(ref_label), path(isolate), path(reference)

    output:
    tuple val(sample), val(ref_label), path("${isolate[0].baseName}.sam"), path("${reference}")

    script:
    """
    minimap2 -ax sr -o ${isolate[0].baseName}.sam ${reference} ${isolate[0]} ${isolate[1]}
    """
}

process samtools {

    publishDir 'samtools'

    input:
    tuple val(sample), val(ref_label), path(minimapOut), path(reference)

    output:
    tuple val(sample), val(ref_label), path("${minimapOut.baseName}_${reference}.sorted.bam"), path("${reference}")

    script:
    """
    samtools view -b ${minimapOut} > ${minimapOut.baseName}_${reference}.bam
    samtools sort ${minimapOut.baseName}_${reference}.bam > ${minimapOut.baseName}_${reference}.sorted.bam
    """
}

process freebayes {

    publishDir 'freebayes'

    input:
    tuple val(sample), val(ref_label), path(bamFile), path(reference)

    output:
    tuple val(ref_label), path(reference), path("${bamFile.baseName}.vcf")

    script:
    """
    freebayes -f ${reference} ${bamFile} > ${bamFile.baseName}.vcf
    """
}

process vcf_subset {
    publishDir 'final_vcf'

    input:
    tuple val(ref_label), path(reference), path(vcf), val(start), val(end)
    path(snps)

    output:
    path("*")

    script:
    """
    bcftools view ${vcf} -Oz -o ${vcf}.gz
    bcftools index ${vcf}.gz
    python3 $baseDir/23SsnpAnalysis.py ${vcf}.gz ${reference} ${start} ${end} ${snps}
    """

}

// additional SNP analysis

process amrfinder {
    publishDir 'amrfinderplus'

    input:
    path(genomes)

    output:
    path("*.out")

    script:
    """
    amrfinder.pl -n ${genomes} -o ${genomes.baseName}.amr.out --plus
    """
}

// Phylogenetic analysis 

process snippy {

    publishDir 'snippyOut'

    input:
    tuple val(ref_label), path(reference)
    path(snippyInput)

    output:
    path("*")

    script:
    """
    snippy-multi ${snippyInput} --ref ${reference[0]} --cpus 16 > runme.sh
    sh runme.sh
    """
}

workflow {

    Channel.fromFilePairs(params.input_dir)
           .set {paired_reads}           

    genomes = assembly(paired_reads)

    references = Channel.fromPath(params.reference_dir)
        .map(file -> tuple(file.baseName.replaceAll(/.fna/,""), file))

    // create tuple of genome combinations with the references
    genomes
        .map(file -> tuple(file, file.simpleName.replaceFirst(/.fasta/,"")))
        .combine(references)
        .set {inputSet}
    
    // getting results of fastANI and grouping based on the samples
    results = fastANI(inputSet)
    grouped = results.groupTuple()

    // getting the best reference for each isolate/sample 
    out = bestRef(grouped)
        .combine(paired_reads, by:0)
        .map { tuple( it[1], *it ) }
        .combine(references, by:0)
        .map {it[1..-1]}

    // running minimap2 to get sam alignment file of isolates to the best reference
    // running samtools to convert sam files to bam and then sorting the bam files
    samOut = minimap2(out)
    bamOut = samtools(samOut)
    freebayesOut = freebayes(bamOut)
        
    // snps bed file with known snp regions
    snp_regions = channel.fromPath(params.snp23S)

    // as SNPs are present in the 23SrRNA, extract SNPS in that region based on 23SrRNA position in the reference genome
    ref23S_path = Channel.fromPath(params.ref23S)
                .splitCsv(header:true)
                .map {row -> 
                    def newRefName = row.reference.minus(~/\.fna$/)
                    return tuple(newRefName,row.start,row.end)
               }
    bcfInput = freebayesOut.combine(ref23S_path, by:0)
    vcf_subset(bcfInput, snp_regions)
    
    // additional SNP analysis
    amrfinder(genomes)
    snippy(references,params.snippyInput)

}