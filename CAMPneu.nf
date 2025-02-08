#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.version = '1.1.1'
params.input = ''
params.output = '/scicomp/groups-pure/OID/NCIRD/DBD/RDB/PRS/PRS_ABiL_URDO/mycoplasma/CAMPneu/LOD_methodValidation'
params.snpFile = ''
params.help = false

///// HELP MESSAGE /////

if (params.help) {
        help = """
              |Usage: 
              |CAMPneu.nf --input <fastq_reads_dir> --output <output_dir>
              |
              |Optional run:
              |CAMPneu.nf --input <fastq_reads_dir> --output <output_dir> --snpFile [path/to/snps.bed]
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
    
    //publishDir "${HOME}/CAMPneu/db/krakendb", mode: 'copy'
    publishDir "${params.kraken_db_dir}", mode: 'copy'

    output:
    path("minikraken_8GB_202003")
 
    script:
    """
    mkdir -p minikraken_8GB_202003
    wget -O minikraken_8GB_202003.tgz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240904.tar.gz
    tar -xvzf minikraken_8GB_202003.tgz -C minikraken_8GB_202003
    rm minikraken_8GB_202003.tgz
    """
}

process download_refs {

    //publishDir "${HOME}/CAMPneu/db/References", mode: 'copy'
    publishDir "${params.reference_dir}", mode: 'copy'

    output:
    //tuple path('GCF_000027345.1_ASM2734v1_genomic.fna'), path('GCF_001272835.1_ASM127283v1_genomic.fna')
    path('GCF_000027345.1_ASM2734v1_genomic.fna'), emit: ref1
    path('GCF_001272835.1_ASM127283v1_genomic.fna'), emit: ref2

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
    printf "NC_000912.1\\t120272\\t120273\\nNC_000912.1\\t121167\\t121168\\nNC_000912.1\\t122118\\t122119\\nNC_000912.1\\t122119\\t122120\\nNC_000912.1\\t122486\\t122487\\nNC_000912.1\\t122666\\t122667\\nNC_000912.1\\t122672\\t122673" > snp_ref1.bed
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

    publishDir "${params.output}/Kraken", mode: 'copy', pattern: '*tsv'

    input:
    tuple val(sampleID), path(read1), path(read2), path(db)   

    output:
    tuple val(sampleID), path(read1), path(read2), env(qc), emit: kraken_out
    tuple val(sampleID), path("${sampleID}_Kraken.tsv"), emit: kraken_report
    tuple val(sampleID), env(percent), env(sp), env(qc), emit: kraken_summary


    script:
    """
    kraken2 -db ${db} \
    --threads 16 \
    --report ${sampleID}.report \
    --paired ${read1} ${read2} > ${sampleID}.Kraken.out

    grep -w "S" ${sampleID}.report | head -n 1 | awk '{printf "%-10s%s_%s\\n", \$1, \$6, \$7}' > ${sampleID}_Kraken.tsv
    sp=\$(awk '{print \$2}' ${sampleID}_Kraken.tsv)
    percent=\$(awk '{print \$1}' ${sampleID}_Kraken.tsv)
    percent_int=\${percent%.*}

    if [[ \${percent_int} -ge 90 && \${sp} == "Mycoplasmoides_pneumoniae" ]]; then
        qc="PASS"
        sp="Mycoplasma_pneumoniae"
    else
        qc="FAIL"
        sp="NA"
        percent="NA"
    fi
    """
}

process fastp {
    publishDir "${params.output}/qc_reads", mode: 'copy', pattern: '*.fq'
    
    input:
    tuple val(sampleID), path(read1), path(read2), val(qc)

    output:
    tuple val(sampleID), path("${read1.baseName}_qc.fq"), path("${read2.baseName}_qc.fq"), path("${read1.baseName}.json"), val(qc), emit: fastp_out

    shell:
    """
    fastp \
    --in1 ${read1} \
    --in2 ${read2} \
    --out1 ${read1.baseName}_qc.fq \
    --out2 ${read2.baseName}_qc.fq \
    --average_qual 30 \
    --json ${read1.baseName}.json
    """
}

process fastp_jq {
    publishDir "${params.output}/fastp", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(sampleID), path(read1), path(read2), path(json), val(qc)

    output:
    tuple val(sampleID), path(read1), path(read2), env(fastp_qc_new), emit: fastp_out
    tuple val(sampleID), path("${read1.baseName}_fastpQC.tsv"), emit: fastp_report
    tuple val(sampleID), env(rate), env(avg_qscore), env(fastp_qc_new), emit: fastp_summary


    shell:
    """
    fastp_qc=\$(if [ "\$(jq '.summary.after_filtering.q30_bases > 0' ${json})" = true ]; then echo "PASS"; else echo "FAIL"; fi)
    jq -r '.summary | [.before_filtering.total_reads, .after_filtering.total_reads, .after_filtering.q30_rate] | @csv' ${json} | awk -F ',' '{print \$1 "\\t" \$2 "\\t" \$3}' > ${read1.baseName}.tsv
    rate=\$(cut -f 3 ${read1.baseName}.tsv | grep '^[0-9].*')

    if [ "${qc}" == "PASS" ] && [ "\${fastp_qc}" == PASS ]; then
        avg_q1=\$(jq '(.read1_before_filtering.quality_curves.mean | add / length )' ${json})
        avg_q2=\$(jq '(.read2_before_filtering.quality_curves.mean | add / length )' ${json})
        avg_qscore=\$(awk "BEGIN {print (\$avg_q1 + \$avg_q2)/2}")
        awk -v avg_qscore=\$avg_qscore '{print \$0 "\\t" avg_qscore}' ${read1.baseName}.tsv > temp && mv temp ${read1.baseName}.tsv
        echo -e "Total_reads_before_filtering\tTotal_reads_after_filtering\tQ30_rate\tAvg_QScore" | cat - ${read1.baseName}.tsv > temp && mv temp ${read1.baseName}.tsv
        awk '{printf "%-30s\\t%-30s\\t%-20s\\t%-20s\\n", \$1, \$2, \$3, \$4}' ${read1.baseName}.tsv > ${read1.baseName}_fastpQC.tsv
        fastp_qc_new="PASS"
    elif [ "\${fastp_qc}" == PASS ] && [ "${qc}" == "FAIL" ]; then
        > ${read1.baseName}_fastpQC.tsv
        echo "sample failed quality check" > ${read1.baseName}_fastpQC.tsv
        fastp_qc_new="FAIL"
        avg_qscore="NA"
        rate="NA"
    fi
    """
}

process coverage_check {

    publishDir "${params.output}/Coverage_check", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(sampleID), path(qc_read1), path(qc_read2), val(qc), val(ref), path(reference)

    output:
    tuple val(sampleID), path(qc_read1), path(qc_read2), env(qc), emit: cov_out
    tuple val(sampleID), path("${sampleID}.cov.tsv"), emit: cov_report
    tuple val(sampleID), env(coverage), env(qc_cov), env(qc), emit: cov_summary

    script:
    """
    if [ "${qc}" == "PASS" ]; then
        minimap2 -ax sr -o ${sampleID}.sam ${reference} ${qc_read1} ${qc_read2}
        samtools view -b ${sampleID}.sam > ${sampleID}.bam
        samtools sort ${sampleID}.bam > ${sampleID}.sorted.bam
        samtools coverage ${sampleID}.sorted.bam > ${sampleID}.tsv
        awk '{printf "%-10s\\t%-10s\\t%-10s\\t%-10s\\t%-10s\\t%-10s\\t%-10s\\t%-10s\\n", \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8}' ${sampleID}.tsv > ${sampleID}.cov.tsv
        coverage=\$(cut -f 7 ${sampleID}.cov.tsv | grep '^[0-9].*')
        coverage=\${coverage%.*}

        if [ "\${coverage}" -gt 30 ]; then
            qc_cov="PASS-Coverage>30x"
            qc="PASS"
        elif [ "\${coverage}" -ge 10 ]; then
            qc_cov="PASS-Coverage<30x"
            qc="PASS"
        else
            qc_cov="FAIL-Coverage<10x"
            qc="FAIL"
        fi
    else
        echo "Sample failed quality check" > ${sampleID}.cov.tsv
        qc_cov="FAIL"
        qc="FAIL"
        coverage=0
    fi
    """

}

process assembly {

    publishDir "${params.output}/assemblies", mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(sampleID), path(read1), path(read2), val(qc)

    output:
    tuple val(sampleID), path("${sampleID}.fasta"), env(qc_new), emit: genomes 
    tuple val(sampleID), path(read1), path(read2), env(qc_new), emit: assembly_out 

    script:
    """
    if [ "${qc}" == "PASS" ]; then
        unicycler -1 ${read1} -2 ${read2} -o ${sampleID} --min_fasta_length 500 -t $task.cpus
        mv ./${sampleID}/assembly.fasta ./${sampleID}.fasta
        qc_new="PASS"
        if [ ! -s ./${sampleID}.fasta ]; then
            qc_new="FAIL"
            touch ${sampleID}.fasta
            echo ">${sampleID}" > ${sampleID}.fasta
            echo "Empty assmebly file generated.\nPossible reasons could be low quality of input reads, incorrect/incomplete data or contamination/misclassified data" >> ${sampleID}.fasta
        fi
    else
        qc_new="FAIL"
        touch ${sampleID}.fasta
        echo ">${sampleID}" > ${sampleID}.fasta
        echo "Skipping assembly for ${sampleID} due to QC failure" >> ${sampleID}.fasta
    fi
    """
}

// Perform mlst 

process mlst {

    publishDir "${params.output}/mlst", mode: 'copy'

    input:
    tuple val(sample), path(assembly), val(qc)

    output:
    tuple val(sample), path("${sample}.mlst.out"), emit: mlst_report
    tuple val(sample), env(st), env(profile), env(qc), emit: mlst_summary

    script:
    """
    if [ "${qc}" == "PASS" ]; then
        mlst ${assembly} --legacy --scheme mpneumoniae > ${sample}.mlst.out
        profile=\$(awk 'NR==2 {print \$4","\$5","\$6","\$7","\$8","\$9","\$10","\$11}' ${sample}.mlst.out)
        if [[ "\${profile}" == *~* ]]; then 
            st="Novel_Allele"
        else
            st=\$(awk '{if (NR==2) print \$3}' ${sample}.mlst.out)
        fi
        qc="PASS"
    else
        touch ${sample}.mlst.out
        echo "no sequence typing data for failed sample" > ${sample}.mlst.out
        st="NA"
        profile="NA"
        qc="FAIL"
    fi
    """
}

// best ref

process fastANI{

    publishDir "${params.output}/fastANI", mode: 'copy', pattern: '*.out'

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

    publishDir "${params.output}/bestReference", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(sample), path(ani_res), val(ref_label), path(reference), val(type), val(qc)

    output:
    tuple val(sample), env(ref), env(type_new), env(qc_new), emit: bestRef_out
    tuple val(sample), path("${sample}_bestRef.tsv"), emit:bestRef_report
    tuple val(sample), env(type_new), env(ani), emit: bestRef_summary 

    script:
    """
    if [ "${qc}" == "[PASS, PASS]" ]; then
        cat ${ani_res} > ${sample}_allRef.txt
        sort -n -k 3 ${sample}_allRef.txt | tail -n 1 > ${sample}_bestRef.txt
        ref=\$(less ${sample}_bestRef.txt | cut -f2 | sed 's/.fna//')
        type_new=\$(less ${sample}_bestRef.txt | cut -f4 | cut -d' ' -f2)
        ani=\$(less ${sample}_bestRef.txt | cut -f3 | cut -d' ' -f2)
        # Print the header with specific column widths
        awk 'BEGIN {printf "%-30s\\t%-40s\\t%-20s\\t%-20s\\n", "sampleID", "Reference", "ANI_value", "type"}' > ${sample}_bestRef.tsv
        # Print the data rows using the same column widths
        awk '{printf "%-30s\\t%-40s\\t%-20s\\t%-20s\\n", \$1, \$2, \$3, \$4}' ${sample}_bestRef.txt >> ${sample}_bestRef.tsv
        qc_new="PASS"
    else
        touch ${sample}_bestRef.txt
        echo "No best reference because sample failed QC" >> ${sample}_bestRef.tsv
        ref="NO_BEST_REF"
        qc_new="FAIL"
        type_new="NA"
        ani="NA"
    fi
    """
}

/// Post assembly processing to identify SNPs from passed samples

process minimap2 {

    publishDir "${params.output}/minimap2", mode: 'copy', pattern: '*.sam'

    input:
    tuple val(sample), path(read1), path(read2), val(qc), val(type), path(reference)

    output:
    tuple val(sample), path("${reference}"), path("${read1.baseName}.sam"), val(qc)

    script:
    """
    minimap2 -ax sr -o ${read1.baseName}.sam ${reference} ${read1} ${read2}
    """
}

process samtools {

    publishDir "${params.output}/samtools", mode: 'copy', pattern: '*.bam'

    input:
    tuple val(sample), path(reference), path(minimapOut), val(qc)

    output:
    tuple val(sample), path("${reference}"), path("${minimapOut.baseName}.sorted.bam"), val(qc)

    script:
    """
    samtools view -b ${minimapOut} > ${minimapOut.baseName}.bam
    samtools sort ${minimapOut.baseName}.bam > ${minimapOut.baseName}.sorted.bam
    samtools index ${minimapOut.baseName}.sorted.bam
    """
}

process freebayes {

    publishDir "${params.output}/freebayes", mode: 'copy', pattern: '*.vcf'

    input:
    tuple val(sample), path(reference), path(bamFile), val(qc)

    output:
    tuple val(sample), path(reference), path("${bamFile.baseName}.vcf"), val(qc)

    script:
    """
    freebayes -f ${reference} --ploidy 1 ${bamFile} > ${bamFile.baseName}.vcf
    """
}

process vcf_subset {
    publishDir "${params.output}/final_vcf", mode: 'copy', pattern: '*.vcf'

    input:
    tuple val(sample), path(reference), path(vcf), val(qc), val(start), val(end), path(snps)

    output:
    tuple val(sample), path("${vcf.simpleName}_all23S.subset.vcf"), path("${vcf.simpleName}_identified.snps.vcf"), val(qc)

    script:
    """
    bcftools view -i 'QUAL>=30' ${vcf} -Oz -o ${vcf}.gz
    bcftools index ${vcf}.gz
    chrom=\$(bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' ${vcf}.gz | head -1 | cut -f 1)
    bcftools view -r \${chrom}:${start}-${end} --no-header ${vcf}.gz > ${vcf.simpleName}_all23S.subset.vcf
    bcftools view -R ${snps} --no-header ${vcf}.gz > ${vcf.simpleName}_identified.snps.vcf
    """
}

/// AMR Gene identification

process amrfinder {
    publishDir "${params.output}/amrfinderplus", mode: 'copy', pattern: '*.out'

    input:
    tuple val(sample), path(fasta), val(qc)

    output:
    path("*.out")

    script:
    """
    if [ "${qc}" == "PASS" ]; then
        amrfinder -n ${fasta} -o ${fasta.baseName}.amr.out

        #check if amr genes are identified
        #header is created so we know the output will have atleast one line, checking that
        num_lines=\$(wc -l < ${fasta.baseName}.amr.out)
        if [ "\${num_lines}" -le 1 ]; then
            >  ${fasta.baseName}.amr.out
            echo "No AMR genes were identified" > ${fasta.baseName}.amr.out
        fi
    else
        touch ${fasta.baseName}.amr.out
        echo "FAILED SAMPLE" >> ${fasta.baseName}.amr.out
    fi
    """
}


///// SUMMARIES

process qc_summary {

    input:
    val(sampleSummaries)

    output:
    path("final_summary.txt")

    script:
    def header = "SampleID\tKraken_Percent\tKraken_Class\tFastp_Score_Rate\tAverage_QC\tCoverage_X\tQC\tP1_Type\tANI\tST\tST_Profile"
    def row = sampleSummaries.join("\n")
    """
    touch final_summary_1.txt
    newlines="Summary of Kraken Classification and QC thresholds\nReads below a qscore of 30 are marked as FAILED and dropped\nReads with coverage below 10X are marked as failed and dropped\n"
    echo ${header} >> final_summary_1.txt
    echo '${row}' >> final_summary_1.txt
    awk '{printf "%-30s\\t%-35s\\t%-40s\\t%-20s\\t%-15s\\t%-15s\\t%-20s\\t%-15s\\t%-25s\\t%-15s\\t%-15s\\t%-15s\\n", \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12}' final_summary_1.txt > final_summary.txt
    echo -e "\$newlines" | cat - final_summary.txt > temp.txt && mv temp.txt final_summary.txt
    """
}

process snp_summary {

    input:
    tuple val(sample), path(allSnps), path(mrSnps), val(qc)

    output:
    tuple val(sample), path("${sample}_snps.txt")

    script:
    """
    if [ "${qc}" == "PASS" ]; then
        if [ -s "${mrSnps}" ]; then
            touch ${sample}_snps.txt
            cut -f1-2,4-5 ${mrSnps} >> ${sample}_snps.txt
            awk '{ new_col = \$2 - 120056; print "${sample}", \$0, \$3 new_col \$4 }' ${sample}_snps.txt > ${sample}_snps_out.txt
            awk '{\$2=""; print \$0}' ${sample}_snps_out.txt | awk '{print \$0, "Resistant"}' > ${sample}_snps.txt
        else
            touch ${sample}_snps.txt
            echo "NA NA NA NA" | awk '{print \$0, "Sensitive"}' >> ${sample}_snps_1.txt
            awk '{ print "${sample}", \$0 }' ${sample}_snps_1.txt > ${sample}_snps.txt
        fi
    else 
        touch ${sample}_snps.txt
    fi
    """

}

process combine_snp_summary {
    publishDir "${params.output}/summary", mode: 'copy'

    input:
    path(snp_files)
    path(mid_summary)

    output:
    path("summary_stats.txt")

    script:
    """
    echo -e "Sample\tPos\tREF\tALT\tSNP\tMacrolide_Susceptibility(Sensitive/Resistant)" >> header.tsv
    cat header.tsv ${snp_files} > snp_summary_1.txt
    awk '{printf "%-30s\\t%-10s\\t%-10s\\t%-10s\\t%-10s\\t%-10s\\t%-10s\\n", \$1, \$2, \$3, \$4, \$5, \$6, \$7}' snp_summary_1.txt > snp_summary.txt
    echo -e "\nMacrolide resistant SNP analysis\nOnly passed samples are tested for presence of SNPs\n\n" | cat - snp_summary.txt > temp.txt && mv temp.txt snp_summary.txt
    cat ${mid_summary} snp_summary.txt > summary_stats.txt
    """
}

process combine_reports{

    publishDir "${params.output}/sample_reports", mode: 'copy'
    
    input:
    tuple val(sample), path(kraken), path(fastp), path(coverage), path(mlst), path(bestRef), path(mrSnps)

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
    echo "Sequence typing using MLST\n" >> ${sample}_report.out
    cat ${mlst} >> ${sample}_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> ${sample}_report.out
    echo "FASTani to select the best reference\n" >> ${sample}_report.out
    cat ${bestRef} >> ${sample}_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> ${sample}_report.out
    echo "Identification of Macrolide Resistant SNPs using Freebayes and bcftools" >> ${sample}_report.out
    echo -e "Sample\tPos\tALT\tREF\tSNP\tType" | cat - ${mrSnps} > temp && mv temp ${mrSnps}
    cat ${mrSnps} >> ${sample}_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> ${sample}_report.out
    """
}
 
workflow {

    ///// KRAKEN DB /////
    // Check if minikraken database exists, if not, download to $CONDA_PREFIX
    def kraken_db_dir = file("${params.kraken_db_dir}/minikraken_8GB_202003")
    if (kraken_db_dir.exists()) {
        println "Minikraken database exists, skipping download."
        kraken_db = Channel.value(kraken_db_dir)
    } else {
        println "Minikraken database downloading now..."
        kraken_db = downloadKrakenDB()        
    }

    /// REFERENCE FILES - TYPE 1 AND TYPE 2 /////
    def ref1 = file("${params.reference_dir}/GCF_000027345.1_ASM2734v1_genomic.fna")
    def ref2 = file("${params.reference_dir}/GCF_001272835.1_ASM127283v1_genomic.fna")    
    if (ref1.exists() && ref2.exists()) {
        println "Type1 and Type2 Reference files exist, skipping download"
        references = Channel.fromPath([ref1, ref2])
    } else {
        println "Reference files downloading now..."
        download_refs()
        references = download_refs.out.ref1
                        .concat(download_refs.out.ref2)
    }

    references_ch = references
        .map { 
            file ->
            def id = file.getName().split('\\.1')[0]  // Extract the identifier 
            [id, file]  // Return a tuple of [id, file]
        }
        .view()
}
/*
    ///// UNZIP ZIPPED READS /////
    if (!params.input) {
        error "ERROR: Missing required input parameter. Please specify the input directory using '--input'."
    }
    if (!params.output) {
        error "ERROR: Missing required output parameter. Please specify the output directory using '--output'."
    }

    Channel.fromFilePairs("${params.input}/*_{1,2,R1,R2,r1,r2}*.{fastq,fq,FASTQ,FQ,fastq.gz,fq.gz,FASTQ.GZ,FQ.GZ}")
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

    // Run fastp to filter reads based on Q-scores and assign QC value of PASS or FAIL
    fastp_json = fastp(kraken_run.kraken_out)
    
    fastp_run = fastp_jq(fastp_json)

    //Run the coverage check
    cov_input = fastp_run.fastp_out
                .combine(references_ch.first())

    cov_check = coverage_check(cov_input)

    // assembling the QC and Coverage threshold passed reads
    assemblies = assembly(cov_check.cov_out)

    //mlst
    mlst = mlst(assemblies.genomes)
    
    passed_samples = assemblies.genomes
                    .filter { tuple -> tuple[-1] == "PASS" }

    //passed_samples.ifEmpty { error "All the input samples failed QC and assembly. Terminating..."}

    passed_samples.map{ it[0] }
                  .collect().subscribe {samples ->
                        println samples.isEmpty() ? "Error: All the input reads failed QC. Terminating..." : "${samples.size()} samples passed QC."    
                }
    
    ref_type = Channel.of(
        ["GCF_000027345", "type1"],
        ["GCF_001272835", "type2"]
    ).combine(references_ch, by:0)

    inputSet = assemblies.genomes.combine(ref_type)
    

    // getting results of fastANI for all samples compared to both references and grouping based on the samples
    results = fastANI(inputSet)
    grouped = results.fastANI_out.groupTuple()

    // getting the best reference for each isolate/sample 
    out = bestRef(grouped)

    ref_type = ref_type.map {it[1..2]}
    
    passed_samples = assemblies.assembly_out
                    .combine(ref_type)
                    .filter { tuple -> tuple[-2] == "type1"}
    

    minimapOut = minimap2(passed_samples)
    samOut = samtools(minimapOut)
    freebayesOut = freebayes(samOut)

    inputVcf = Channel.of(["120057", "122961", snpFile.bed])
    inputVcf = freebayesOut.combine(inputVcf)

    cleanedChannel = inputVcf.map { tuple ->
        def path = tuple[6].get() // Assuming the DataflowVariable is the 6th element
        tuple[0..5] + [path]      // Return the tuple with the path instead of the DataflowVariable
    }
    
    vcf_out = vcf_subset(cleanedChannel)
    
    // amrfinder
    amrfinder(assemblies.genomes)

    // SUMMARIES

    snpSumOut = snp_summary(vcf_out)
    snpSumOut.map { it[1] }
             .collect()
             .set { snpSumOut1 }

    combined = kraken_run.kraken_report
            .combine(fastp_run.fastp_report, by:0)
            .combine(cov_check.cov_report, by:0) 
            .combine(mlst.mlst_report, by:0)
            .combine(out.bestRef_report, by:0)
            .combine(snpSumOut, by:0)
            
    combine_reports(combined)

    // summary
    //summary_values = [sampleID ,kraken_percent, kraken_class, fastp_score_rate, avgQscore, coverage_x, cov_qc, type, ani]

    summary_values = kraken_summary
                    .combine(fastp_run.fastp_summary, by:0)
                    .map {it[0..2] + it[4..5]}
                    .combine(cov_check.cov_summary, by:0)
                    .map {it[0..6]}
                    .combine(out.bestRef_summary, by:0)
                    .combine(mlst.mlst_summary, by:0)
                    .map {it[0..10]}
                    .collect{ summary -> summary.join('\t')}
                    .set { sampleSummaries }
    
    qc_summary = qc_summary(sampleSummaries)

    // summarising results from SNP analysis and combining with previous results to create one single run summary

    combine_snp_summary(snpSumOut1, qc_summary)
}
*/