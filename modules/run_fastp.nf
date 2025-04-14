process fastp {
    cpus params.max_cpus
    publishDir "${params.output}/qc_reads", mode: 'copy', pattern: '*.fq'
    
    input:
    tuple val(sampleID), path(read1), path(read2), val(qc)

    output:
    tuple val(sampleID), path("${read1.baseName}_qc.fq"), path("${read2.baseName}_qc.fq"), path("${read1.baseName}.json"), val(qc), emit: fastp_out

    shell:
    """
    fastp \
    --thread ${task.cpus} \
    --in1 ${read1} \
    --in2 ${read2} \
    --out1 ${read1.baseName}_qc.fq \
    --out2 ${read2.baseName}_qc.fq \
    --average_qual 30 \
    --json ${read1.baseName}.json
    """
}