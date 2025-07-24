process minimap2 {
    cpus params.max_cpus
    publishDir "${params.output}/minimap2", mode: 'copy', pattern: '*.sam'

    input:
    tuple val(sample), path(read1), path(read2), val(qc), val(type), path(reference)

    output:
    tuple val(sample), path("${reference}"), path("${read1.baseName}.bam"), val(qc)

    script:
    """
    minimap2 -t $task.cpus -ax sr ${reference} ${read1} ${read2} | samtools sort -@ task.cpus -o ${read1.baseName}.bam
    samtools index ${read1.baseName}.bam
    """
}