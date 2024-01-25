#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_dir = "/scicomp/home-pure/ubt4/mycoplasma/allRawData/*_{1,2}.fastq"
params.reference_dir = "/scicomp/home-pure/ubt4/mycoplasma/nextflow/updated/references/*.fna"

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

    // output:
    // path("*_bestRef.txt")

    output:
    tuple val(sample), env(ref)

    script:
    """
    cat ${ani_res} > ${sample}_allRef.txt
    ref=\$(sort -n -k 3 ${sample}_allRef.txt | tail -n 1 | cut -f2 | sed 's/.fna//')
    """

    // script:
    // """
    // cat ${ani_res} > ${sample}_allRef.txt
    // sort -n -k 3 ${sample}_allRef.txt | tail -n 1 | cut -f1-2  > ${sample}_bestRef.txt
    // """    

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
//     tuple path(isolate), path(reference)

//     output:
//     path("*")

//     script:
//     """
//     snippy --cpus 16 --outdir ${isolate[0].baseName} --ref ${reference} --R1 ${isolate[0]} --R2 ${isolate[1]} --force
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
    path("*.vcf")

    script:
    """
    freebayes -f ${reference} ${bamFile} > ${bamFile.baseName}.vcf
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

    results = fastANI(inputSet)
    grouped = results.groupTuple()

    out = bestRef(grouped)
        .combine(paired_reads, by:0)
        .map { tuple( it[1], *it ) }
        .combine(ref_label, by:0)
        .map {it[1..-1]}
        .map {it[1..-1]}
        .map {it[1..-1]}

    bestRef_track(grouped)
    
    // Phylogenetic analysis
    //snippyOut = snippy(out)

    // SNP analysis
    minimapOut = minimap2(out)
    minimapOut.view()
    bamOut = samtools(minimapOut)
    freebayes(bamOut)
    amrfinder(genomes)
}