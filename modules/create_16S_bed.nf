process create_16S_bed {
    output:
    path("snp_16S_tet.bed"), emit: bed_16S

    script:
    """
    printf "NC_000912.1\\t119505\\t119506\\nNC_000912.1\\t119280\\t119281\\n" > snp_16S_tet.bed
    """ 
}