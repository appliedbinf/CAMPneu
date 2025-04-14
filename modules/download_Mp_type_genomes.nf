process download_refs {

    publishDir "${params.reference_dir}", mode: 'copy'

    output:
    path('GCF_000027345.1_ASM2734v1_genomic.fna'), emit: ref1
    path('GCF_001272835.1_ASM127283v1_genomic.fna'), emit: ref2

    script:
    """
    wget -O GCF_000027345.1_ASM2734v1_genomic.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/345/GCF_000027345.1_ASM2734v1/GCF_000027345.1_ASM2734v1_genomic.fna.gz 
    wget -O GCF_001272835.1_ASM127283v1_genomic.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/272/835/GCF_001272835.1_ASM127283v1/GCF_001272835.1_ASM127283v1_genomic.fna.gz
    gunzip GCF_000027345.1_ASM2734v1_genomic.fna.gz
    gunzip GCF_001272835.1_ASM127283v1_genomic.fna.gz
    """
}