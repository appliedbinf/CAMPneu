process bestRef {

    publishDir "${params.output}/bestReference", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(sample), path(ani_res), val(ref_label), path(reference), val(type), val(qc)

    output:
    tuple val(sample), env(ref), env(type_new), env(qc_new), emit: out
    tuple val(sample), path("${sample}_bestRef.tsv"), emit: report
    tuple val(sample), env(type_new), env(ani), emit: summary

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