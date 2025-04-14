process downloadKrakenDB {
    
    publishDir "${params.kraken_db_dir}", mode: 'copy'

    output:
    path("minikraken_8GB_202003")
 
    script:
    """
    mkdir -p minikraken_8GB_202003
    wget -O minikraken_8GB_202003.tgz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240904.tar.gz
    tar -xvzf minikraken_8GB_202003.tgz -C minikraken_8GB_202003
    rm minikraken_8GB_202003.tgz
    """
}