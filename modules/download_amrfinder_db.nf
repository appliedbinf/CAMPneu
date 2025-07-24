process download_amrfinder_db {
    
    publishDir "${params.amrfinder_db}", mode: 'copy'

    output:
    path("latest")
 
    script:
    """
    amrfinder_update -d ./

    cp -r ./latest/ ./latest-backup
    rm -r latest
    mv latest-backup latest
    """
}