process minimap2 {
    cpus params.max_cpus
    publishDir "${params.output}/minimap2", mode: 'copy', pattern: '*.sam'

    input:
    tuple val(sample), path(read1), path(read2), val(qc), val(type), path(reference)

    output:
    tuple val(sample), path("${reference}"), path("${minimapOut.baseName}.bam"), val(qc)

    script:
    """
    minimap2 -t $task.cpus -ax sr -o ${read1.baseName}.sam ${reference} ${read1} ${read2}
    samtools view -h -@ $task.cpus ${minimapOut} | samtools sort -@ $task.cpus -o ${minimapOut.baseName}.bam
    samtools index ${minimapOut.baseName}.bam
    """
}