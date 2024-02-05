#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.input_dir = "/scicomp/home-pure/ubt4/mycoplasma/nextflow/simReads/*_{R1,R2}.fq"
//params.input_dir = "/scicomp/home-pure/ubt4/mycoplasma/allRawData/*_{1,2}.fastq"
params.reference_dir = "/scicomp/home-pure/ubt4/mycoplasma/nextflow/updated/references/*.fna"
params.ref23S = "/scicomp/home-pure/ubt4/mycoplasma/nextflow/updated/23S_reference_positions.csv"
params.snp23S = "/scicomp/home-pure/ubt4/mycoplasma/nextflow/updated/23sSNPS.bed"

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
    tuple val(sample), path(genome), path(reference)

    output: 
    tuple val(sample), path("${sample}_${reference.baseName}.out"), emit: ani_results

    script:
    """
    fastANI -q ${genome} -r ${reference} -o ${sample}_${reference.baseName}.out
    """

}

process bestRef {

    publishDir 'out'

    input:
    tuple val(sample), path(ani_res)

    output:
    tuple val(sample), env(ref)

    script:
    """
    cat ${ani_res} > ${sample}_allRef.txt
    ref=\$(sort -n -k 3 ${sample}_allRef.txt | tail -n 1 | cut -f2 | sed 's/.fna//')
    """
}

process bestRef_track {

    publishDir "out"

    input:
    tuple val(sample), path(ani_res)
    
    output:
    path("*_bestRef.txt")

    script:
    """
    cat ${ani_res} > ${sample}_allRef.txt
    sort -n -k 3 ${sample}_allRef.txt | tail -n 1 | cut -f1-2  > ${sample}_bestRef.txt
    """
}

// Phylogeny

// process snippy {

//     publishDir 'snippyOut'

//     input:
//     path(reads) 
//     path(reference)

//     output:
//     path("*")

//     script:
//     """
//     snippy --cpus 16 --outdir ${reads.baseName} --ref ${reference[0]} --R1 ${reads[0]} --R2 ${reads[1]} --force
//     """
// }

// SNP analysis

process minimap2 {

    publishDir 'minimap2'

    input:
    tuple path(isolate), path(reference)

    output:
    tuple path("${isolate[0].baseName}.sam"), path("${reference}")

    script:
    """
    minimap2 -ax sr -o ${isolate[0].baseName}.sam ${reference} ${isolate[0]} ${isolate[1]}
    """
}

process samtools {

    publishDir 'samtools'

    input:
    tuple path(minimapOut), path(reference)

    output:
    tuple path("${minimapOut.baseName}_${reference}.sorted.bam"), path("${reference}")

    script:
    """
    samtools view -b ${minimapOut} > ${minimapOut.baseName}_${reference}.bam
    samtools sort ${minimapOut.baseName}_${reference}.bam > ${minimapOut.baseName}_${reference}.sorted.bam
    """
}


process freebayes {

    publishDir 'freebayes'

    input:
    tuple path(bamFile), path(reference)

    output:
    path("*")

    script:
    """
    freebayes -f ${reference} ${bamFile} > ${bamFile.baseName}.vcf
    bcftools view ${bamFile.baseName}.vcf -Oz -o ${bamFile.baseName}.vcf.gz
    bcftools index ${bamFile.baseName}.vcf.gz
    """
}

process vcf_subset {
    publishDir 'final_vcf'

    input:
    path(vcf)
    tuple path(reference), val(start), val(end)
    path(snps)

    output:
    path("*")

    script:
    """
    23SsnpAnalysis.py ${vcf} ${reference} ${start} ${end} ${snps}
    """

}

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
workflow {

    Channel.fromFilePairs(params.input_dir)
           .set {paired_reads}           

    genomes = assembly(paired_reads)

    // channels for genomes and references
    //genomes = Channel.fromPath(params.input_dir)

    references = Channel.fromPath(params.reference_dir)
    ref_label = references
        .map(file -> tuple(file.baseName.replaceAll(/.fna/,""), file))


    // create tuple of genome combinations with the references
    genomes
        .map(file -> tuple(file.simpleName.replaceFirst(/.fasta/,""), file))
        .combine(references)
        .set {inputSet}

    // classification based on fastANI output
    results = fastANI(inputSet)
    grouped = results.groupTuple()

    // find the best reference and creating a tuple of paired reads and reference
    out = bestRef(grouped)
        .combine(paired_reads, by:0)
        .map { tuple( it[1], *it ) }
        .combine(ref_label, by:0)
        .map {it[1..-1]}
        .map {it[1..-1]}
        .map {it[1..-1]}

    bestRef_track(grouped)
    
    // Phylogenetic analysis
    // snippyOut = snippy(paired_reads, references)

    // SNP analysis
    minimapOut = minimap2(out)
    bamOut = samtools(minimapOut)
    vcfout = freebayes(bamOut)
    amrfinder(genomes)

    // snp file
    snp_regions = channel.fromPath(params.snp23S)

    // as SNPs are present in the 23SrRNA, extract SNPS in that region based on 23SrRNA position in the reference genome
    ref23S_path = Channel
                .fromPath(params.ref23S)
                .splitCsv(header:true)
                // .map { row -> tuple(row.reference, row.start, row.end)}
                .map {row -> 
                    def newRefName = row.reference.minus(~/\.fna$/)
                    return tuple(newRefName,row.start,row.end)
               }

    //channel with reference path and start and end locations of the 23S rRNA region
    combined = ref_label.combine(ref23S_path, by:0)
                .map {it[1..-1]}


    vcf_subset(vcfout, combined, snp_regions)
}