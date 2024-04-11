#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_dir = "/home/mkadam7/CAMPneu/Simulated_Reads/*_{1,2}.fq"
params.reference_dir = "/home/mkadam7/CAMPneu/References/*.fna"
params.ref23S = "/home/mkadam7/CAMPneu/23S_reference_positions.csv"
params.krakendb = "/home/mkadam7/CAMPneu/minikraken2_v2_8GB_201904_UPDATE/"

process kraken {

    publishDir 'Kraken'

    input:
    tuple val(sampleID), path(reads)
    path(db)   

    output:
    tuple path("${sampleID}.Kraken.out"), path("${sampleID}.report")
    
    script:
    """
    kraken2 -db ${db} \
    --threads 16 \
    --report ${sampleID}.report \
    --paired ${reads[0]} ${reads[1]} > ${sampleID}.Kraken.out
    """
}

process fastp {

    publishDir 'fastp'

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), path("${reads[0].baseName}_qc.fq"), path("${reads[1].baseName}_qc.fq"), env(qc)
    

    shell:
    """
    fastp \
      --in1 ${reads[0]} \
      --in2 ${reads[1]} \
      --out1 ${reads[0].baseName}_qc.fq \
      --out2 ${reads[1].baseName}_qc.fq \
      --average_qual 30 \
      --json ${reads[0].baseName}.json

    qc=\$(if [ "\$(jq '.summary.before_filtering.q30_bases > 0' ${reads[0].baseName}.json)" = true ]; then echo "PASS"; else echo "FAIL"; fi)
    """

    emit:
    file("${reads[0].baseName}.json")

}

process coverage_check {

    publishDir 'Coverage_check'

    input:
    tuple val(sampleID), path(qc_read1), path(qc_read2), val(ref), path(reference)

    output:
    tuple val(sampleID), path(qc_read1), path(qc_read2), env(coverage)

    script:
    """
    minimap2 -ax sr -o ${sampleID}.sam ${reference} ${qc_read1} ${qc_read2}
    samtools view -b ${sampleID}.sam > ${sampleID}.bam
    samtools sort ${sampleID}.bam > ${sampleID}.sorted.bam
    samtools coverage ${sampleID}.sorted.bam > ${sampleID}.cov.txt
    coverage=\$(awk 'NR==2 { if (\$7 > 90) print "PASS"; else print "FAIL" }' ${sampleID}.cov.txt)
    """

    emit:
    file("${sampleID}.cov.txt")

}

process combine{
    
    input:
    file(fastp_report)   // Input file for Fastp report
    file(coverage_report)

    output:
    path(${combined_report})

    script:
    """
    cat ${fastp_report} ${coverage_report} > ${combined_report}
    """
}

process assembly {

    publishDir 'assemblies'

    input:
    tuple val(sampleID), path(read1), path(read2)

    output:
    path("*.fasta"), emit: genomes // emit the sample ID as well

    script:
    """
    unicycler -1 ${read1} -2 ${read2} -o ${sampleID} 
    mv ./${sampleID}/assembly.fasta ./${sampleID}.fasta
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

// process minimap2 {

//     publishDir 'minimap2'

//     input:
//     tuple val(sample), val(ref_label), path(isolate), path(reference)

//     output:
//     tuple val(sample), val(ref_label), path("${isolate[0].baseName}.sam"), path("${reference}")

//     script:
//     """
//     minimap2 -ax sr -o ${isolate[0].baseName}.sam ${reference} ${isolate[0]} ${isolate[1]}
//     """
// }

// process samtools {

//     publishDir 'samtools'

//     input:
//     tuple val(sample), val(ref_label), path(minimapOut), path(reference)

//     output:
//     tuple val(sample), val(ref_label), path("${minimapOut.baseName}_${reference}.sorted.bam"), path("${reference}")

//     script:
//     """
//     samtools view -b ${minimapOut} > ${minimapOut.baseName}_${reference}.bam
//     samtools sort ${minimapOut.baseName}_${reference}.bam > ${minimapOut.baseName}_${reference}.sorted.bam
//     """
// }

// process freebayes {

//     publishDir 'freebayes'

//     input:
//     tuple val(sample), val(ref_label), path(bamFile), path(reference)

//     output:
//     tuple val(ref_label), path(reference), path("${bamFile.baseName}.vcf")

//     script:
//     """
//     freebayes -f ${reference} ${bamFile} > ${bamFile.baseName}.vcf
//     """
// }

// process vcf_subset {
//     publishDir 'final_vcf'

//     input:
//     tuple val(ref_label), path(reference), path(vcf), val(start), val(end), path(snps)

//     output:
//     path("*")

//     script:
//     """
//     bcftools view ${vcf} -Oz -o ${vcf}.gz
//     bcftools index ${vcf}.gz
//     python3 $baseDir/23SsnpAnalysis.py ${vcf}.gz ${reference} ${start} ${end} ${snps}
//     """
// }

// // additional SNP analysis

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

    kraken(paired_reads, params.krakendb)
    qual_check_out = fastp(paired_reads)  


    filteredQual = qual_check_out
                    .filter { tuple -> tuple[-1] == "PASS" }
                    .map { tuple -> tuple.take(tuple.size() - 1)}
    
    references = Channel.fromPath(params.reference_dir)
        .map(file -> tuple(file.baseName.replaceAll(/.fna/,""), file))

    
    input = filteredQual
                .combine(references.first())
    
    cov_check_out = coverage_check(input)

    filteredCov = cov_check_out
                    .filter { tuple -> tuple[-1] == "PASS" }
                    .map { tuple -> tuple.take(tuple.size() - 1)}
    filteredCov.view()


    // assembling the QC and Coverage threshold passed reads
    genomes = assembly(filteredCov)

    // create tuple of genome combinations with the references
    genomes
        .map(file -> tuple(file, file.simpleName.replaceFirst(/.fasta/,"")))
        .combine(references)
        .set {inputSet}
    
    // getting results of fastANI and grouping based on the samples
    results = fastANI(inputSet)
    grouped = results.groupTuple()

    // // getting the best reference for each isolate/sample 
    out = bestRef(grouped)
        .combine(paired_reads, by:0)
        .map { tuple( it[1], *it ) }
        .combine(references, by:0)
        .map {it[1..-1]}

    out.view()

    // // running minimap2 to get sam alignment file of isolates to the best reference
    // // running samtools to convert sam files to bam and then sorting the bam files
    // samOut = minimap2(out)
    // bamOut = samtools(samOut)
    // freebayesOut = freebayes(bamOut)

    // // as SNPs are present in the 23SrRNA, extract SNPS in that region based on 23SrRNA position in the reference genome
    // // along with path to snps bed file with known snp regions 
    // ref23S_path = Channel.fromPath(params.ref23S)
    //             .splitCsv(header:true)
    //             .map{row -> tuple(row.reference,row.start,row.end,row.snpsFile)}

    // ref23S_path.view()

    // bcfInput = freebayesOut.combine(ref23S_path, by:0)
    // vcf_subset(bcfInput)
    
    // // additional SNP analysis
    amrfinder(genomes)
}
