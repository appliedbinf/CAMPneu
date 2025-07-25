process mlst {

    publishDir "${params.output}/mlst", mode: 'copy'

    input:
    tuple val(sample), path(assembly), val(qc)

    output:
    tuple val(sample), path("${sample}.mlst.out"), emit: report
    tuple val(sample), env(st), env(profile), env(qc), emit: summary

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
        echo "No sequence typing data for failed sample" > ${sample}.mlst.out
        st="NA"
        profile="NA"
        qc="FAIL"
    fi
    """
}