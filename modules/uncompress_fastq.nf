process uncompress_reads {

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), path("${reads[0].simpleName}_unzip.fastq"), path("${reads[1].simpleName}_unzip.fastq")

    script:
    """
    if [[ "${reads[0]}" == *.gz ]]; then
        gunzip -c "${reads[0]}" > "${reads[0].simpleName}_unzip.fastq"
    else
        mv ${reads[0]} ${reads[0].simpleName}_unzip.fastq
    fi
    if [[ "${reads[1]}" == *.gz ]]; then
        gunzip -c "${reads[1]}" > "${reads[1].simpleName}_unzip.fastq"
    else
        mv ${reads[1]} ${reads[1].simpleName}_unzip.fastq
    fi  
    """
}