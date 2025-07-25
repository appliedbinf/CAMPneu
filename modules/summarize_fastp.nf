process summarize_fastp {
    publishDir "${params.output}/fastp", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(sampleID), path(read1), path(read2), path(json), val(qc)

    output:
    tuple val(sampleID), path(read1), path(read2), env(fastp_qc_new), emit: out
    tuple val(sampleID), path("${read1.baseName}_fastpQC.tsv"), emit: report
    tuple val(sampleID), env(rate), env(avg_qscore), env(fastp_qc_new), emit: summary


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
        fastp_qc_new="Failed_QC"
        avg_qscore="NA"
        rate="NA"
    fi
    """
}