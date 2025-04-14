process fastANI{
    cpus params.max_cpus
    publishDir "${params.output}/fastANI", mode: 'copy', pattern: '*.out'

    input:
    tuple val(sample), path(assembly), val(qc), val(ref_label), val(type), path(reference)

    output: 
    tuple val(sample), path("${sample}_${ref_label}_fastANI.out"), val(ref_label), path(reference), val(type), val(qc), emit: fastANI_out
    tuple val(sample), path("${sample}_${ref_label}_fastANI.out"), emit: fastANI_report

    script:
    """
    if [ "${qc}" == "PASS" ]; then
        fastANI -t $task.cpus -q ${assembly} -r ${reference} --minFraction 0.5 -o ${sample}_fastANI_1.out
        awk -F'\t' 'BEGIN {OFS="\t"} {print \$0, "${type}"}' ${sample}_fastANI_1.out > ${sample}_fastANI_2.out
        cut -f1-3,6 ${sample}_fastANI_2.out > ${sample}_${ref_label}_fastANI.out
    else 
        touch ${sample}_fastANI.out
        echo ">${sample}" > ${sample}_fastANI.out
        echo "Sample skipped due to QC failure" >> ${sample}_${ref_label}_fastANI.out
    fi
    """
}