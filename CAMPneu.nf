#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// params.input_dir = "/scicomp/home-pure/ubt4/mycoplasma/nextflow/simulation_tests/type2reads_vs_type1ref_lowQual/reads/*_{1,2}.fq"
// params.reference_dir = "/scicomp/home-pure/ubt4/mycoplasma/nextflow/simulation_tests/type1reads_vs_type1ref/GCF_000027345.1_ASM2734v1_genomic.fna"
// params.ref23S = "/scicomp/home-pure/ubt4/mycoplasma/nextflow/updated/23S_reference_positions.csv"

params.input_dir = "/home/wengland7/campneu_testfiles/*_{1,2}.fastq"
params.reference_dir = "${baseDir}/references/"
params.ref23S = "${baseDir}/23S_reference_positions.csv"
params.reference_type2="${baseDir}/references/GCF_001272835.1_ASM127283v1_genomic.fna"
params.known_snps="${baseDir}/knownSNPs.txt"

process assembly {

    publishDir 'assemblies'

    input:
    tuple val(sampleID), path(reads)

    output:
    path("*.fasta"), emit: genomes // emit the sample ID as well

    script:
    """
    unicycler -1 ${reads[0]} -2 ${reads[1]} -o ${sampleID} 
    mv ./${sampleID}/assembly.fasta ./${sampleID}.fasta
    """
}

process snpCheck {

    publishDir 'snpCheck'

    input:
    tuple path(assembly), val(sample), path(t2reference), path(snps)

    output:
    path("${sample}.known_snps")

    script:
    """
    nucmer -p ${sample} ${t2reference} ${assembly}
    dnadiff -d ${sample}.delta -p ${sample}
    show-snps -ClrT ${sample}.delta > ${sample}.snpcheck
    bash ${baseDir}/snpcheck.sh ${sample}.snpcheck ${snps} > ${sample}.known_snps
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
    tuple val(ref_label), path(reference), path(vcf), val(start), val(end), path(snps)

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
    amrfinder -n ${genomes} -o ${genomes.baseName}.amr.out --plus
    """
}

workflow {

    Channel.fromFilePairs(params.input_dir)
           .set {paired_reads}           

    genomes = assembly(paired_reads)

    Channel.fromPath(params.reference_type2)
        .set {t2ref}

    Channel.fromPath(params.known_snps)
        .set {known_snps}

    // align genomes to Type 2 reference to identify known SNPs between Types 1 & 2
    genomes
        .map(file -> tuple(file, file.simpleName.replaceFirst(/.fasta/,"")))
        .combine(t2ref)
        .combine(known_snps)
        .set {snpSet}

    snpCheck(snpSet)

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

    // as SNPs are present in the 23SrRNA, extract SNPS in that region based on 23SrRNA position in the reference genome
    // along with path to snps bed file with known snp regions 
    ref23S_path = Channel.fromPath(params.ref23S)
                .splitCsv(header:true)
                .map{row -> tuple(row.reference,row.start,row.end,row.snpsFile)}

    ref23S_path.view()

    bcfInput = freebayesOut.combine(ref23S_path, by:0)
    vcf_subset(bcfInput)
    
    // additional SNP analysis
    amrfinder(genomes)
}
