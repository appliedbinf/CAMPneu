#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_dir = "/home/mkadam7/CAMPneu/Simulated_Reads/*_{1,2}.fq"
params.reference_dir = "/home/mkadam7/CAMPneu/References/*.fna"
params.krakendb = "/home/mkadam7/CAMPneu/minikraken2_v2_8GB_201904_UPDATE/"
params.refType = "/home/mkadam7/CAMPneu/References/Reference_type.csv"
params.snpsBed = "/home/mkadam7/CAMPneu/test/snps_ref1.bed"


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
    cat header.tsv ${sampleID}.tsv > ${sampleID}_Kraken.tsv
    sp=\$(grep -w 'Species' ${sampleID}_Kraken.tsv | cut -f3 | sed 's/^[ \t]*//' | sed 's/ /_/g')
    percent=\$(grep -w 'Species' ${sampleID}_Kraken.tsv | cut -f1 | sed 's/^[ \t]*//')
    """
}

process fastp {

    publishDir 'fastp'

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), path("${reads[0].baseName}_qc.fq"), path("${reads[1].baseName}_qc.fq"), env(qc), emit: fastp_out
    tuple val(sampleID), path("${reads[0].baseName}.tsv"), emit: fastp_report
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
    """

}

process coverage_check {

    publishDir 'Coverage_check'

    input:
    tuple val(sampleID), path(qc_read1), path(qc_read2), val(qc), val(ref), path(reference)

    output:
    tuple val(sampleID), path(qc_read1), path(qc_read2), env(coverage), emit: cov_out
    tuple val(sampleID), path("${sampleID}.cov.txt"), emit: cov_report
    tuple val(sampleID), env(percent), env(coverage), emit: cov_summary

    script:
    """
    if [ "${qc}" == "FAIL" ]; then
        echo "Sample failed coverage check" > ${sampleID}.cov.txt
        coverage="FAIL"
        percent=0
    else
        minimap2 -ax sr -o ${sampleID}.sam ${reference} ${qc_read1} ${qc_read2}
        samtools view -b ${sampleID}.sam > ${sampleID}.bam
        samtools sort ${sampleID}.bam > ${sampleID}.sorted.bam
        samtools coverage ${sampleID}.sorted.bam > ${sampleID}.cov.txt
        coverage=\$(awk 'NR==2 { if (\$7 > 90) print "PASS"; else print "FAIL" }' ${sampleID}.cov.txt)
        percent=\$(cut -f 7 ${sampleID}.cov.txt | grep '^[0-9].*')
    fi
    """

}

process assembly {

    publishDir 'assemblies'

    input:
    tuple val(sampleID), path(read1), path(read2), val(qc)

    output:
    tuple val(sampleID), path("${sampleID}.fasta"), val(qc), emit: genomes // emit the sample ID as well

    script:
    """
    if [ "${qc}" == "PASS" ]; then
        unicycler -1 ${read1} -2 ${read2} -o ${sampleID} 
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
        awk '{print \$0, "${type}"}' ${sample}_fastANI_1.out > ${sample}_${ref_label}_fastANI.out
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
    tuple val(sample), path("${sample}_bestRef.txt"), emit:bestRef_report
    tuple val(sample), env(type_new), env(qc_new), emit: bestRef_summary 

    script:
    """
    if [ "${qc}" == "[PASS, PASS]" ]; then
        cat ${ani_res} > ${sample}_allRef.txt
        sort -n -k 3 ${sample}_allRef.txt | tail -n 1 > ${sample}_bestRef.txt
        ref=\$(less ${sample}_bestRef.txt | cut -f2 | sed 's/.fna//')
        type_new=\$(less ${sample}_bestRef.txt | cut -f5 | cut -d' ' -f2)
        qc_new="PASS"

    else
        touch ${sample}_bestRef.txt
        echo "No best reference because sample failed QC" >> ${sample}_bestRef.txt
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
    echo "Coverage Filtering using Samtools Coverage\nSamples below 100x are marked as FAILED\n" >> ${sample}_report.out
    cat ${coverage} >> ${sample}_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> ${sample}_report.out
    echo "FASTani to select the best reference\n" >> ${sample}_report.out
    cat ${bestRef} >> ${sample}_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> ${sample}_report.out
    """
}

// create a summary for everything done so far

process writeToFile {

    publishDir "summary"

    input:
    val(sampleSummaries)

    output:
    path("final_summary_1.txt"), emit: summary1

    script:
    """
    echo "SampleID Kraken_Percent Kraken_Class Fastp_Score_Rate Coverage_X Type QC" > summary.txt
    for summary in "${sampleSummaries.join('\t')}"; do
        echo "\$summary" | sed 's/,/\\t/g' >> summary.txt
    done
    mv summary.txt final_summary.txt
    sed -i 's/sample/\\n&/g' final_summary.txt
    header=\$(head -n 1 final_summary.txt)
    data=\$(tail -n +2 final_summary.txt | sort -nk1)
    echo "Summary of Kraken classification and QC thresholds\nReads below a qscore of 30 are marked as FAILED and dropped\nReads with Coverage below 100x are marked as FAILED and dropped\n" > final_summary_1.txt
    echo "\$header\n\$data" | column -t >> final_summary_1.txt

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
    path("*")

    script:
    """
    bcftools view ${vcf} -Oz -o ${vcf}.gz
    bcftools index ${vcf}.gz
    python3 $baseDir/23SsnpAnalysis.py ${vcf}.gz ${reference} ${start} ${end} ${snps}
    """
}

workflow {

    Channel.fromFilePairs(params.input_dir)
           .set {paired_reads}         

    // Run Kraken and generate Kraken classification and report
    kraken_out = kraken(paired_reads, params.krakendb)

    // Run fastp to assign PASS or FAIL check to reads based on Q-scores
    qual_check_out = fastp(paired_reads)

    // References
    references = Channel.fromPath(params.reference_dir)
        .map(file -> tuple(file.baseName.replaceAll(/.fna/,""), file))

    // Create input for coverage check where PASS reads with go through coverage filter
    input = qual_check_out.fastp_out
                .combine(references.first())

    // Run the coverage check
    cov_check = coverage_check(input)

    // // assembling the QC and Coverage threshold passed reads
    genomes = assembly(cov_check.cov_out)

    // Read a csv file that has the type for the references
    ref_type = Channel.fromPath(params.refType)
        .splitCsv(header:true)
        .map{row  -> tuple(row.Reference,row.Type)}
        .combine(references, by:0)

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
    //kraken_out.kraken_summary.view()
    //qual_check_out.fastp_summary.view()
    //cov_check.cov_summary.view()
    //out.bestRef_summary.view()

    //summary_values = [sampleID ,kraken_percent, kraken_class, fastp_score_rate, coverage_x, type, qc]
    summary_values = kraken_out.kraken_summary
                    .combine(qual_check_out.fastp_summary, by:0)
                    .map {it[0..3]}
                    .combine(cov_check.cov_summary, by:0)
                    .map {it[0..4]}
                    .combine(out.bestRef_summary, by:0)
                    .collect()
                    .set { sampleSummaries }

    writeToFile(sampleSummaries)


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
    vcf_subset(inputVcf)

}
