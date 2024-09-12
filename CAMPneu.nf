#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input = ''
params.output = ''
params.snpFile = ''
params.help = false

///// HELP MESSAGE /////

if (params.help) {
        help = """
              |Usage: 
              |CAMPneu.nf --input <fastq_reads_dir> --output <output_dir>
              |
              |Custom Run: Run all the processes in pipeline using the downloaded database and references with a user specified SNP bed file
              |Command: nextflow run CAMPneu.nf --input <fastq_reads_dir> --output <output_dir> --snpFile <snp_bed_file>
              |         
              |Required arguments:     
              |  --input     Path to the Paired Fastq Reads directory  
              |  --output    Directory where process outputs are saved          
              |Optional arguments:  
              |  --snpFile   Path to the custom SNP bed file
              |  --help      Print this message and exit""".stripMargin()

    println(help)
    exit(0)
}

process downloadKrakenDB {
    publishDir "${CONDA_PREFIX}/bin/data/krakendb/"

    output:
    path("minikraken_8GB_202003")
 
    script:
    """
    mkdir -p minikraken_8GB_202003
    wget -O minikraken_8GB_202003.tgz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240112.tar.gz
    tar -xvzf minikraken_8GB_202003.tgz -C minikraken_8GB_202003
    rm minikraken_8GB_202003.tgz
    """
}

process download_refs {
    publishDir "${CONDA_PREFIX}/bin/data/References"

    output:
    tuple path('GCF_000027345.1_ASM2734v1_genomic.fna'), path('GCF_001272835.1_ASM127283v1_genomic.fna')
 
    script:
    """
    wget -O GCF_000027345.1_ASM2734v1_genomic.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/345/GCF_000027345.1_ASM2734v1/GCF_000027345.1_ASM2734v1_genomic.fna.gz 
    wget -O GCF_001272835.1_ASM127283v1_genomic.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/272/835/GCF_001272835.1_ASM127283v1/GCF_001272835.1_ASM127283v1_genomic.fna.gz
    gunzip GCF_000027345.1_ASM2734v1_genomic.fna.gz
    gunzip GCF_001272835.1_ASM127283v1_genomic.fna.gz
    """
}

process createBedFile {
    output:
    path("snp_ref1.bed"), emit: bed

    script:
    """
    cat > snp_ref1.bed <<EOF
    NC_000912.1 120272	120273
    NC_000912.1	121167	121168
    NC_000912.1	122118	122119
    NC_000912.1	122119	122120
    NC_000912.1	122486	122487
    NC_000912.1	122666	122667
    EOF
    """
}
 
process gunzip_reads {

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), path("${reads[0].simpleName}_unzip.fastq"), path("${reads[1].simpleName}_unzip.fastq")

    script:
    """
    if [[ "${reads[0]}" == *.gz ]]; then
        gunzip -c "${reads[0]}" > "${reads[0].simpleName}_unzip.fastq"
    else
        mv ${reads[0]} ${reads[0].simpleName}_unzip.fastq
    fi
    if [[ "${reads[1]}" == *.gz ]]; then
        gunzip -c "${reads[1]}" > "${reads[1].simpleName}_unzip.fastq"
    else
        mv ${reads[1]} ${reads[1].simpleName}_unzip.fastq
    fi  
    """
}
 
process kraken {



    publishDir "${params.output}/Kraken"

    input:
    tuple val(sampleID), path(read1), path(read2), path(db)   

    output:
    tuple val(sampleID), path(read1), path(read2), env(qc), emit: kraken_out
    tuple val(sampleID), path("${sampleID}_Kraken.tsv"), emit: kraken_report
    tuple val(sampleID), env(percent), env(sp), emit: kraken_summary


    script:
    """
    kraken2 -db ${db} \
    --threads 16 \
    --report ${sampleID}.report \
    --paired ${read1} ${read2} > ${sampleID}.Kraken.out

    awk -F'\\t' '{if (\$4 == "G") {\$4 = "Genus"} else if (\$4 == "S") {\$4 = "Species"}; if ((\$4 == "Genus" || \$4 == "Species") && \$1 > 5) {gsub(/[[:space:]]+\$/, "", \$NF); print \$1 "\t" \$4 "\t" \$6}}' ${sampleID}.report > ${sampleID}.tsv
    echo -e "Percent_Reads_Covered\tRank\tScientific_name" > header.tsv
    cat header.tsv ${sampleID}.tsv > ${sampleID}_Kraken_1.tsv
    sp=\$(grep -w 'Species' ${sampleID}_Kraken_1.tsv | cut -f3 | sed 's/^[ \t]*//' | sed 's/ /_/g')
    percent=\$(grep -w 'Species' ${sampleID}_Kraken_1.tsv | cut -f1 | sed 's/^[ \t]*//')
    qc="FAIL"
    if [ \$(echo "\${percent} >= 90" | bc) -eq 1 ]; then
        qc="PASS"
        column -t ${sampleID}_Kraken_1.tsv > ${sampleID}_Kraken.tsv
    else
        touch ${sampleID}_Kraken.tsv
        echo "Sample failed classfication check." > ${sampleID}_Kraken.tsv       
    fi
    """
}

process fastp {

    publishDir "${params.output}/fastp"

    input:
    tuple val(sampleID), path(read1), path(read2), val(qc)

    output:
    tuple val(sampleID), path("${read1.baseName}_qc.fq"), path("${read2.baseName}_qc.fq"), env(fastp_qc), emit: fastp_out
    tuple val(sampleID), path("${read1.baseName}_fastpQC.tsv"), emit: fastp_report
    tuple val(sampleID), env(rate), env(fastp_qc), emit: fastp_summary


    shell:
    """
    if [ "${qc}" == "PASS" ]; then
        fastp \
        --in1 ${read1} \
        --in2 ${read2} \
        --out1 ${read1.baseName}_qc.fq \
        --out2 ${read2.baseName}_qc.fq \
        --average_qual 30 \
        --json ${read1.baseName}.json

        fastp_qc=\$(if [ "\$(jq '.summary.before_filtering.q30_bases > 0' ${read1.baseName}.json)" = true ]; then echo "PASS"; else echo "FAIL"; fi)
        jq -r '.summary | [.before_filtering.total_reads, .after_filtering.total_reads, .before_filtering.q30_rate] | @csv' ${read1.baseName}.json | awk -F ',' '{print \$1 "\\t" \$2 "\\t" \$3}' > ${read1.baseName}.tsv
        echo -e "Total_reads_before_filtering\tTotal_reads_after_filtering\tQ30_rate" | cat - ${read1.baseName}.tsv > temp &&  mv temp ${read1.baseName}.tsv
        rate=\$(cut -f 3 ${read1.baseName}.tsv | grep '^[0-9].*')
        column -t ${read1.baseName}.tsv > ${read1.baseName}_fastpQC.tsv
    fi
    """
}

process coverage_check {

    publishDir "${params.output}/Coverage_check"

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

    publishDir "${params.output}/assemblies"

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

    publishDir "${params.output}/fastANI"

    input:
    tuple val(sample), path(assembly), val(qc), val(ref_label), val(type), path(reference)

    output: 
    tuple val(sample), path("${sample}_${ref_label}_fastANI.out"), val(ref_label), path(reference), val(type), val(qc), emit: fastANI_out
    tuple val(sample), path("${sample}_${ref_label}_fastANI.out"), emit: fastANI_report

    script:
    """
    if [ "${qc}" == "PASS" ]; then
        fastANI -q ${assembly} -r ${reference} --minFraction 0.5 -o ${sample}_fastANI_1.out
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

    publishDir "${params.output}/bestReference"

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

process minimap2 {

    publishDir "${params.output}/minimap2"

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

    publishDir "${params.output}/samtools"

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

    publishDir "${params.output}/freebayes"

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
    publishDir "${params.output}/final_vcf"

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

process amrfinder {
    publishDir "${params.output}/amrfinderplus"

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

///// SUMMARIES

process summary1 {

    input:
    val(sampleSummaries)

    output:
    path("final_summary.txt")

    script:
    def header = "SampleID\tKraken_Percent\tKraken_Class\tFastp_Score_Rate\tCoverage_X\tType\tQC"
    def row = sampleSummaries.join("\n")
    """
    touch final_summary_1.txt
    newlines="Summary of Kraken Classification and QC thresholds\nReads below a qscore of 30 are marked as FAILED and dropped\nReads with coverage below 10X are marked as failed and dropped\n"
    echo ${header} >> final_summary_1.txt
    echo '${row}' >> final_summary_1.txt
    column -t final_summary_1.txt > final_summary.txt
    echo -e "\$newlines" | cat - final_summary.txt > temp.txt && mv temp.txt final_summary.txt
    """
}

process each_snp_summary {

    input:
    tuple val(sample), path(mrSnps), path(allSnps)

    output:
    tuple val(sample), path("${sample}_snps.txt")

    script:
    """
    if [ -s "${mrSnps}" ]; then
        touch ${sample}_snps.txt
        cut -f1-2,4-5 ${mrSnps} >> ${sample}_snps.txt
        awk '{ new_col = \$2 - 120055; print "${sample}", \$0, \$3 new_col \$4 }' ${sample}_snps.txt > ${sample}_snps_out.txt
        awk '{\$2=""; print \$0}' ${sample}_snps_out.txt > ${sample}_snps.txt
    else
        touch ${sample}_snps.txt
        echo "No macrolide resistant SNPs observed" >> ${sample}_snps_1.txt
        awk '{ print "${sample}", \$0 }' ${sample}_snps_1.txt > ${sample}_snps.txt
    fi
    """

}

process combine_snp_summary {
    publishDir "${params.output}/summary"

    input:
    path(snp_files)
    path(mid_summary)

    output:
    path("summary_stats.txt")

    script:
    """
    echo -e "Sample\tPos\tALT\tREF\tSNP" >> header.tsv
    cat header.tsv ${snp_files} | column -t > snp_summary.txt
    echo -e "\nMacrolide resistant SNP analysis\n" | cat - snp_summary.txt > temp.txt && mv temp.txt snp_summary.txt
    cat ${mid_summary} snp_summary.txt > summary_stats.txt
    """
}

process combine_reports{

    publishDir "${params.output}/sample_reports"
    
    input:
    tuple val(sample), path(kraken), path(fastp), path(coverage), path(bestRef), path(mrSnps)

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
    echo "Identification of Macrolide Resistant SNPs using Freebayes and bcftools" >> ${sample}_report.out
    echo -e "Sample\tPos\tALT\tREF\tSNP" | cat - ${mrSnps} > temp && mv temp ${mrSnps}
    cat ${mrSnps} | column -t >> ${sample}_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> ${sample}_report.out
    """
}


 
workflow {

    ///// KRAKEN DB /////
    // Check if minikraken database exists, if not, download to $CONDA_PREFIX
    def kraken_db_dir = file("${CONDA_PREFIX}/bin/data/krakendb/minikraken_8GB_202003")
    if (kraken_db_dir.exists()) {
        println "Minikraken database exists, skipping download."
        kraken_db = Channel.value(kraken_db_dir)
    } else {
        println "Minikraken database downloading now..."
        kraken_db = downloadKrakenDB()        
    }

    ///// REFERENCE FILES - TYPE 1 AND TYPE 2 /////
    def ref1 = file("${CONDA_PREFIX}/bin/data/References/GCF_000027345.1_ASM2734v1_genomic.fna")
    def ref2 = file("${CONDA_PREFIX}/bin/data/References/GCF_001272835.1_ASM127283v1_genomic.fna")
    if (ref1.exists() && ref2.exists()) {
        println "Type1 and Type2 Reference files also exist, skipping download"
        references = Channel.fromPath([ref1, ref2])
        // references = tuple(ref1, ref2)
        // references = Channel.of(references)
    } else {
        println "Reference files downloading now..."
        references = download_refs()
    }
    
    ///// UNZIP ZIPPED READS /////
    if (!params.input) {
        error "ERROR: Missing required input parameter. Please specify the input directory using '--input'."
    }
    if (!params.output) {
        error "ERROR: Missing required output parameter. Please specify the output directory using '--output'."
    }

    Channel.fromFilePairs("${params.input}/*_{1,2,R1,R2,r1,r2}.{fastq,fq,FASTQ,FQ,fastq.gz,fq.gz,FASTQ.GZ,FQ.GZ}")
           .ifEmpty{ error "NO {reads}.fastq/fq files found in the specified directory: ${params.input}"}
           .set {paired_reads}
    
    unzipped_reads = gunzip_reads(paired_reads)

    ///// CREATE SNP FILE IF NOT INPUT /////
    if (!params.snpFile) {
        snpFile = createBedFile()
    } else {
        snpFile = params.snpFile
    }

    // Run Kraken and generate Kraken classification and report 
    kraken_input = unzipped_reads.combine(kraken_db)
    kraken_run = kraken(kraken_input)
    kraken_summary = kraken_run.kraken_summary
                     .map{ line -> line.collect {it.replaceAll(' ','|')}}
    
    // Run fastp to filter reads based on Q-scores and assign QC value of PASS or FAIL
    qual_check_out = fastp(kraken_run.kraken_out)
    references_ch = references.map { file ->
        def id = file.getBaseName().split('\\.1')[0]  // Extract the identifier 
        [id, file]  // Return a tuple of [id, file]
    }

    // Run the coverage check
    cov_input = qual_check_out.fastp_out
                .combine(references_ch.first())
    cov_check = coverage_check(cov_input)
    
    // assembling the QC and Coverage threshold passed reads
    assemblies = assembly(cov_check.cov_out)

    ref_type = Channel.of(
        ["GCF_000027345", "type1"],
        ["GCF_001272835", "type2"]
    ).combine(references_ch, by:0)
    
    // create tuple of genome combinations with the references
    assemblies
        .combine(ref_type)
        .set {inputSet}
    
    // getting results of fastANI for all samples compared to both references and grouping based on the samples
    results = fastANI(inputSet)
    grouped = results.fastANI_out.groupTuple()

    // getting the best reference for each isolate/sample 
    out = bestRef(grouped)

    ref_type = ref_type.map {it[1..2]}
    passed_samples = cov_check.cov_out
                    .filter { tuple -> tuple[-1] == "PASS" }
                    .combine(ref_type)
                    .filter { tuple -> tuple[-2] == "type1"}
    
    minimapOut = minimap2(passed_samples)
    samOut = samtools(minimapOut)
    freebayesOut = freebayes(samOut)

    inputVcf = Channel.of(["120057", "122961", snpFile.bed])
    inputVcf = freebayesOut.combine(inputVcf)
    
    cleanedChannel = inputVcf.map { tuple ->
        def path = tuple[5].get() // Assuming the DataflowVariable is the 6th element
        tuple[0..4] + [path]      // Return the tuple with the path instead of the DataflowVariable
    }

    vcf_out = vcf_subset(cleanedChannel)
    // amrfinder
    amrfinder(assemblies)

    // SUMMARIES

    snpSumOut = each_snp_summary(vcf_out)
    snpSumOut.map { it[1] }
             .collect()
             .set { snpSumOut1 }

    combined = kraken_run.kraken_report
            .combine(qual_check_out.fastp_report, by:0)
            .combine(cov_check.cov_report, by:0)  
            .combine(out.bestRef_report, by:0)
            .combine(snpSumOut, by:0)

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

    summary1 = summary1(sampleSummaries)

    // summarising results from SNP analysis and combining with previous results to create one single run summary

    combine_snp_summary(snpSumOut1, summary1)
}
