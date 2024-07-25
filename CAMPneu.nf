#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_dir = null
params.reference_dir = null
params.krakendb =null
params.snpsBed = null
params.help = false

if (params.help) {
        help = """nextflow run CAMPneu.nf --input_dir <fastq_reads_dir> --reference_dir <reference_genome_dir>
              |
              |
              |Required arguments:  
              |  --input_dir     Location of the input directory with the Paired Fastq Reads  
              |  --reference_dir Location of directory containing fna files for Mycoplasma Pneumoniae References 
              |                  Type 1 - GCF_000027345.1_ASM2734v1_genomic.fna
              |                  Type 2 - GCF_001272835.1_ASM127283v1_genomic.fna
              |  --krakendb      Path to the Kraken database for Taxonomic Classification 
              |                  Database can be found at : https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240112.tar.gz
              |                  Database ".tar.gz" can be unzipped using: tar -xvzf k2_standard_08gb_20240112.tar.gz
              |  --bed           A bed file containing the positions of Macrolide Resistant snps which is available along with the pipeline
              |              
              |Optional arguments:  
              |  --help           Print this message and exit""".stripMargin()

    println(help)
    exit(0)
}

if (!params.input_dir) {
    println "ERROR: Missing required parameter: --input_dir"
    exit(1)
}

if (!params.reference_dir) {
    println "ERROR: Missing required parameter: --reference_dir"
    exit(2)
}

process kraken {

    publishDir 'Kraken'

    input:
    tuple val(sampleID), path(reads)
    path(db)   

    output:
    tuple path("${sampleID}.Kraken.out"), path("${sampleID}.report"), emit: kraken_out
    tuple val(sampleID), path("${sampleID}_Kraken.tsv"), emit: kraken_report
    tuple val(sampleID), env(percent), env(sp), emit: kraken_summary

    script:
    """
    kraken2 -db ${db} \
    --threads 16 \
    --report ${sampleID}.report \
    --paired ${reads[0]} ${reads[1]} > ${sampleID}.Kraken.out

    awk -F'\\t' '{if (\$4 == "G") {\$4 = "Genus"} else if (\$4 == "S") {\$4 = "Species"}; if ((\$4 == "Genus" || \$4 == "Species") && \$1 > 5) {gsub(/[[:space:]]+\$/, "", \$NF); print \$1 "\t" \$4 "\t" \$6}}' ${sampleID}.report > ${sampleID}.tsv
    echo -e "Percent_Reads_Covered\tRank\tScientific_name" > header.tsv
    cat header.tsv ${sampleID}.tsv > ${sampleID}_Kraken_1.tsv
    sp=\$(grep -w 'Species' ${sampleID}_Kraken_1.tsv | cut -f3 | sed 's/^[ \t]*//' | sed 's/ /_/g')
    percent=\$(grep -w 'Species' ${sampleID}_Kraken_1.tsv | cut -f1 | sed 's/^[ \t]*//')
    column -t ${sampleID}_Kraken_1.tsv > ${sampleID}_Kraken.tsv
    """
}

process fastp {

    publishDir 'fastp'

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), path("${reads[0].baseName}_qc.fq"), path("${reads[1].baseName}_qc.fq"), env(qc), emit: fastp_out
    tuple val(sampleID), path("${reads[0].baseName}_fastpQC.tsv"), emit: fastp_report
    tuple val(sampleID), env(rate), env(qc), emit: fastp_summary


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
    jq -r '.summary | [.before_filtering.total_reads, .after_filtering.total_reads, .before_filtering.q30_rate] | @csv' ${reads[0].baseName}.json | awk -F ',' '{print \$1 "\\t" \$2 "\\t" \$3}' > ${reads[0].baseName}.tsv
    echo -e "Total_reads_before_filtering\tTotal_reads_after_filtering\tQ30_rate" | cat - ${reads[0].baseName}.tsv > temp &&  mv temp ${reads[0].baseName}.tsv
    rate=\$(cut -f 3 ${reads[0].baseName}.tsv | grep '^[0-9].*')
    column -t ${reads[0].baseName}.tsv > ${reads[0].baseName}_fastpQC.tsv
    """

}

process coverage_check {

    publishDir 'Coverage_check'

    input:
    tuple val(sampleID), path(qc_read1), path(qc_read2), val(qc), val(ref), path(reference)

    output:
    tuple val(sampleID), path(qc_read1), path(qc_read2), env(qc), emit: cov_out
    tuple val(sampleID), path("${sampleID}.cov.tsv"), emit: cov_report
    tuple val(sampleID), env(percent), env(coverage), env(qc), emit: cov_summary

    script:
    """
    if [ "${qc}" == "FAIL" ]; then
        echo "Sample failed coverage check" > ${sampleID}.cov.txt
        coverage="FAIL"
        percent=0
        qc="FAIL"
    else
        minimap2 -ax sr -o ${sampleID}.sam ${reference} ${qc_read1} ${qc_read2}
        samtools view -b ${sampleID}.sam > ${sampleID}.bam
        samtools sort ${sampleID}.bam > ${sampleID}.sorted.bam
        samtools coverage ${sampleID}.sorted.bam > ${sampleID}.cov.txt
        coverage=\$(awk 'NR==2 { if (\$7 < 10) print "FAIL"; else if (\$7 <=30 )print "PASS-Coverage<30x"; else print "PASS-Coverage>30x"}' ${sampleID}.cov.txt)
        percent=\$(cut -f 7 ${sampleID}.cov.txt | grep '^[0-9].*')
        column -t ${sampleID}.cov.txt > ${sampleID}.cov.tsv
        qc="PASS"
    fi
    """

}

process assembly {

    publishDir 'assemblies'

    input:
    tuple val(sampleID), path(read1), path(read2), val(qc)

    output:
    tuple val(sampleID), path("${sampleID}.fasta"), val(qc), emit: genomes 

    script:
    """
    if [ "${qc}" == "PASS" ]; then
        unicycler -1 ${read1} -2 ${read2} -o ${sampleID} --min_fasta_length 500 
        mv ./${sampleID}/assembly.fasta ./${sampleID}.fasta
    else
        touch ${sampleID}.fasta
        echo ">${sampleID}" > ${sampleID}.fasta
        echo "Skipping assembly for ${sampleID} due to QC failure" >> ${sampleID}.fasta
    fi
    """
}

process fastANI{

    publishDir 'fastANI'

    input:
    tuple val(sample), path(assembly), val(qc), val(ref_label), val(type), path(reference)

    output: 
    tuple val(sample), path("${sample}_${ref_label}_fastANI.out"), val(ref_label), path(reference), val(type), val(qc), emit: fastANI_out
    tuple val(sample), path("${sample}_${ref_label}_fastANI.out"), emit: fastANI_report

    script:
    """
    if [ "${qc}" == "PASS" ]; then
        fastANI -q ${assembly} -r ${reference} -o ${sample}_fastANI_1.out
        awk -F'\t' 'BEGIN {OFS="\t"} {print \$0, "${type}"}' ${sample}_fastANI_1.out > ${sample}_fastANI_2.out
        cut -f1-3,6 ${sample}_fastANI_2.out > ${sample}_${ref_label}_fastANI.out
    else 
        touch ${sample}_fastANI.out
        echo ">${sample}" > ${sample}_fastANI.out
        echo "Sample skipped due to QC failure" >> ${sample}_${ref_label}_fastANI.out
    fi
    """

}

process bestRef {

    publishDir 'bestReference'

    input:
    tuple val(sample), path(ani_res), val(ref_label), path(reference), val(type), val(qc)

    output:
    tuple val(sample), env(ref), env(type_new), env(qc_new), emit: bestRef_out
    tuple val(sample), path("${sample}_bestRef.tsv"), emit:bestRef_report
    tuple val(sample), env(type_new), emit: bestRef_summary 

    script:
    """
    if [ "${qc}" == "[PASS, PASS]" ]; then
        cat ${ani_res} > ${sample}_allRef.txt
        sort -n -k 3 ${sample}_allRef.txt | tail -n 1 > ${sample}_bestRef.txt
        ref=\$(less ${sample}_bestRef.txt | cut -f2 | sed 's/.fna//')
        type_new=\$(less ${sample}_bestRef.txt | cut -f4 | cut -d' ' -f2)
        echo -e "sampleID\tReference\tANI_value\ttype" | cat - ${sample}_bestRef.txt | column -t > temp &&  mv temp  ${sample}_bestRef.tsv
        qc_new="PASS"

    else
        touch ${sample}_bestRef.txt
        echo "No best reference because sample failed QC" >> ${sample}_bestRef.tsv
        ref="NO_BEST_REF"
        qc_new="FAIL"
        type_new="no_type"
    fi
    """
}

// Individual reports for all samples are generated here
process combine_reports{

    publishDir "sample_reports"
    
    input:
    tuple val(sample), path(kraken), path(fastp), path(coverage), path(bestRef)

    output:
    path("${sample}_report.out")
    
    script:
    """
    touch ${sample}_report.out
    echo "Kraken Classfication\n" >> ${sample}_report.out
    cat ${kraken} >> ${sample}_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> ${sample}_report.out
    echo "Quality Filtering using Fastp\nSamples that have Qscore < 30 are marked as FAILED\n" >> ${sample}_report.out
    cat ${fastp} >> ${sample}_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> ${sample}_report.out
    echo "Coverage Filtering using Samtools Coverage\nSamples below 10x are marked as FAILED\n" >> ${sample}_report.out
    cat ${coverage} >> ${sample}_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> ${sample}_report.out
    echo "FASTani to select the best reference\n" >> ${sample}_report.out
    cat ${bestRef} >> ${sample}_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> ${sample}_report.out
    """
}

// create a summary for everything done so far

process summary1 {

    publishDir "summary"

    input:
    val(sampleSummaries)

    output:
    path("final_summary.txt"), emit: summary1

    script:
    def header = "SampleID\tKraken_Percent\tKraken_Class\tFastp_Score_Rate\tCoverage_X\tType\tQC"
    def row = sampleSummaries.join("\n")
    """
    touch final_summary_1.txt
    echo ${header} >> final_summary_1.txt
    echo '${row}' >> final_summary_1.txt
    column -t final_summary_1.txt > final_summary.txt
    """
}

// The previous steps are done to QC the input reads and subtype them into type1 or type2
// Once reads have passed both Quality and Coverage threshold, all the passed reads are then aligned to type1 reference only
// Freebayes is then used find the SNPs and compare them to known list of SNPS

process minimap2 {

    publishDir 'minimap2'

    input:
    tuple val(sample), path(read1), path(read2), val(qc), val(type), path(reference)

    output:
    tuple val(sample), path("${reference}"), path("${read1.baseName}.sam")

    script:
    """
    minimap2 -ax sr -o ${read1.baseName}.sam ${reference} ${read1} ${read2}
    """
}

process samtools {

    publishDir 'samtools'

    input:
    tuple val(sample), path(reference), path(minimapOut)

    output:
    tuple val(sample), path("${reference}"), path("${minimapOut.baseName}_${reference}.sorted.bam")

    script:
    """
    samtools view -b ${minimapOut} > ${minimapOut.baseName}_${reference}.bam
    samtools sort ${minimapOut.baseName}_${reference}.bam > ${minimapOut.baseName}_${reference}.sorted.bam
    samtools index ${minimapOut.baseName}_${reference}.sorted.bam
    """
}

process freebayes {

    publishDir 'freebayes'

    input:
    tuple val(sample), path(reference), path(bamFile)

    output:
    tuple val(sample), path(reference), path("${bamFile.baseName}.vcf")

    script:
    """
    freebayes -f ${reference} --ploidy 1 ${bamFile} > ${bamFile.baseName}.vcf
    """
}

process vcf_subset {
    publishDir 'final_vcf'

    input:
    tuple val(sample), path(reference), path(vcf), val(start), val(end), path(snps)

    output:
    tuple val(sample), path("*identified.snps.txt"), path("*all23S.subset.txt")

    script:
    """
    bcftools view -i 'QUAL>=30' ${vcf} -Oz -o ${vcf}.gz
    bcftools index ${vcf}.gz
    python3 $baseDir/23SsnpAnalysis.py ${vcf}.gz ${reference} ${start} ${end} ${snps}
    """
}

process snp_summary {
    publishDir 'snp_summary'

    input:
    tuple val(sample), path(mrSnps)

    output:
    path("${sample}_snps.txt")

    script:
    """
    touch ${sample}_snps.txt
    """

}

process amrfinder {
    publishDir 'amrfinderplus'

    input:
    tuple val(sample), path(fasta), val(qc)

    output:
    path("*.out")

    script:
    """
    if [ "${qc}" == "PASS" ]; then
        amrfinder -n ${fasta} -o ${fasta.baseName}.amr.out
    else
        touch ${fasta.baseName}.amr.out
        echo "FAILED SAMPLE" >> ${fasta.baseName}.amr.out
    fi
    """
}

workflow {

    Channel.fromFilePairs("${params.input_dir}/*_{1,2,R1,R2,r1,r2}.{fastq,fq,FASTQ,FQ}")
           .ifEmpty{ error "NO {reads}.fastq/fq files found in the specified directory: ${params.input_dir}"}
           .set {paired_reads}   


    // Run Kraken and generate Kraken classification and report
    kraken_out = kraken(paired_reads, params.krakendb)
    kraken_summary = kraken_out.kraken_summary
              .map{ line -> line.collect {it.replaceAll(' ','|')}}

    // Run fastp to assign PASS or FAIL check to reads based on Q-scores
    qual_check_out = fastp(paired_reads)

    // References
    references = Channel.fromPath("${params.reference_dir}/*.fna")
        .ifEmpty{ error "NO {reference}.fna files found in the specified directory: ${params.reference_dir}"}
        .map(file -> tuple(file.baseName.replaceAll(/.fna/,""), file))

    // Create input for coverage check where PASS reads with go through coverage filter
    input = qual_check_out.fastp_out
                .combine(references.first())

    // Run the coverage check
    cov_check = coverage_check(input)

    // assembling the QC and Coverage threshold passed reads
    genomes = assembly(cov_check.cov_out)

    ref_type = Channel.of(
        ["GCF_000027345.1_ASM2734v1_genomic", "type1"],
        ["GCF_001272835.1_ASM127283v1_genomic", "type2"]
    ).combine(references, by:0)

    // create tuple of genome combinations with the references
    genomes
        .combine(ref_type)
        .set {inputSet}
    
    // getting results of fastANI for all samples compared to both references and grouping based on the samples
    results = fastANI(inputSet)
    grouped = results.fastANI_out.groupTuple()

    // getting the best reference for each isolate/sample 
    out = bestRef(grouped)

    // Combined report for each sample 
    combined = kraken_out.kraken_report
            .combine(qual_check_out.fastp_report, by:0)
            .combine(cov_check.cov_report, by:0)  
            .combine(out.bestRef_report, by:0)

    combine_reports(combined)

    // summary
    //summary_values = [sampleID ,kraken_percent, kraken_class, fastp_score_rate, coverage_x, cov_qc, type]
    summary_values = kraken_summary
                    .combine(qual_check_out.fastp_summary, by:0)
                    .map {it[0..3]}
                    .combine(cov_check.cov_summary, by:0)
                    .map {it[0..5]}
                    .combine(out.bestRef_summary, by:0)
                    .collect{ summary -> summary.join('\t')}
                    .set { sampleSummaries }

    summary1(sampleSummaries)

    // Once the combined reports are generated, the passed samples are moved ahead for further processing

    ref_type = ref_type.map {it[1..2]}
    passed_samples = cov_check.cov_out
                    .filter { tuple -> tuple[-1] == "PASS" }
                    .combine(ref_type)
                    .filter { tuple -> tuple[-2] == "type1"}

    minimapOut = minimap2(passed_samples)
    samOut = samtools(minimapOut)
    freebayesOut = freebayes(samOut)

    inputVcf = Channel.of(["120057", "122961", params.snpsBed])
    inputVcf = freebayesOut.combine(inputVcf)
    vcf_out = vcf_subset(inputVcf)
    vcf_out.view()

    amrfinder(genomes)

}
