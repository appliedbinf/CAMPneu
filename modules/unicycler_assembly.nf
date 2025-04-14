process assembly {
    cpus params.max_cpus
    publishDir "${params.output}/assemblies", mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(sampleID), path(read1), path(read2), val(qc)

    output:
    tuple val(sampleID), path("${sampleID}.fasta"), env(qc_new), emit: genomes 
    tuple val(sampleID), path(read1), path(read2), env(qc_new), emit: assembly_out 

    script:
    """
    if [ "${qc}" == "PASS" ]; then
        unicycler -1 ${read1} -2 ${read2} -o ${sampleID} --min_fasta_length 500 -t $task.cpus
        mv ./${sampleID}/assembly.fasta ./${sampleID}.fasta
        qc_new="PASS"
        if [ ! -s ./${sampleID}.fasta ]; then
            qc_new="FAIL"
            touch ${sampleID}.fasta
            echo ">${sampleID}" > ${sampleID}.fasta
            echo "Empty assmebly file generated.\nPossible reasons could be low quality of input reads, incorrect/incomplete data or contamination/misclassified data" >> ${sampleID}.fasta
        fi
    else
        qc_new="FAIL"
        touch ${sampleID}.fasta
        echo ">${sampleID}" > ${sampleID}.fasta
        echo "Skipping assembly for ${sampleID} due to QC failure" >> ${sampleID}.fasta
    fi
    """
}