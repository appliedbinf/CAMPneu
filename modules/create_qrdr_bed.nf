process create_quinolone_amr_locations {
    output:
    path("quinolones_ref.tsv"), emit: tsv

    script:
    """
    #gene	nucl_change	aa_change	amr_class
    printf "gyrA\tG295A\tAsp99Xaa\tquinolones\ngyrB\tG1327A\tAsp443Xaa\tquinolones\ngyrB\tG1391A\tArg464Lys\tquinolones\ngyrB\tA1448G\tGlu483Gly\tquinolones\nparC\tG241T\tGly81Cys\tquinolones\nparC\tC248T\tAla83Val\tquinolones\nparC\tG259A\tAsp87Xaa\tquinolones\nparE\tC1345T\tPro449Ser\tquinolones\ngyrA\tA141C\tPro47Pro\tTEST" > quinolones_ref.tsv
    """
}