process freebayes {

    publishDir "${params.output}/freebayes", mode: 'copy', pattern: '*.vcf'

    input:
    tuple val(sample), path(reference), path(bamFile), val(qc)

    output:
    tuple val(sample), path(reference), path("${bamFile.baseName}.vcf"), val(qc)

    script:
    """
    freebayes -f ${reference} --ploidy 1 ${bamFile} > ${bamFile.baseName}.vcf
    """
}