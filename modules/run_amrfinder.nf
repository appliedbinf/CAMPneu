process amrfinder {
    cpus params.max_cpus
    publishDir "${params.output}/amrfinderplus", mode: 'copy', pattern: '*.out'

    input:
    tuple val(sample), path(fasta), val(qc)
    path(db)

    output:
    tuple val(sample), path("*.out"), emit: report
    tuple val(sample), env(amrfinder_gene_list), emit: summary

    script:
    """
    if [ "${qc}" == "PASS" ]; then
        amrfinder --threads $task.cpus --database ${db} -n ${fasta} -o ${fasta.baseName}.amr.out

        #check if amr genes are identified
        #header is created so we know the output will have atleast one line, checking that
        num_lines=\$(wc -l < ${fasta.baseName}.amr.out)
        if [ "\${num_lines}" -le 1 ]; then
            >  ${fasta.baseName}.amr.out
            echo "No AMR genes were identified" > ${fasta.baseName}.amr.out
            amrfinder_gene_list="None"
        else
            amrfinder_gene_list=\$(awk '{print \$1}' | paste -d, -s)
        fi
    else
        touch ${fasta.baseName}.amr.out
        echo "FAILED SAMPLE" >> ${fasta.baseName}.amr.out
        amrfinder_gene_list="FAILED SAMPLE"
    fi
    """
}