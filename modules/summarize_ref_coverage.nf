process summarize_reference_coverage {
    cpus params.max_cpus
    publishDir "${params.output}/Coverage_check", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(sampleID), path(qc_read1), path(qc_read2), val(qc), val(ref), path(reference)

    output:
    tuple val(sampleID), path(qc_read1), path(qc_read2), env(qc), emit: out
    tuple val(sampleID), path("${sampleID}.cov.tsv"), emit: report
    tuple val(sampleID), env(coverage), env(qc_cov), env(qc), emit: summary

    script:
    """
    if [ "${qc}" == "PASS" ]; then
        minimap2 -t $task.cpus -ax sr ${reference} ${qc_read1} ${qc_read2} | samtools sort -@ $task.cpus -o ${sampleID}.bam
        samtools coverage ${sampleID}.bam > ${sampleID}.tsv
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