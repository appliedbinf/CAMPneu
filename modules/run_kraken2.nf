process kraken {
    cpus params.max_cpus
    maxForks 1
    publishDir "${params.output}/Kraken", mode: 'copy', pattern: '*tsv'

    input:
    tuple val(sampleID), path(read1), path(read2), path(db)   

    output:
    tuple val(sampleID), path(read1), path(read2), env(qc), emit: out
    tuple val(sampleID), path("${sampleID}_Kraken.tsv"), emit: report
    tuple val(sampleID), env(percent), env(sp), env(qc), emit: summary


    script:
    """
    kraken2 -db ${db} \
    --threads $task.cpus \
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