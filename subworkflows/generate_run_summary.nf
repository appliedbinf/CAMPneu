include { generate_run_report } from '../modules/concat_all_sample_summaries.nf'

workflow RUN_REPORT {
    take:
    kraken_summary
    bestref_summary
    mlst_summary
    fastp_summary
    coverage_summary
    amrfinder_summary
    macrolide_summary
    tet_summary
    quinolone_summary

    main:
    summary_values = kraken_summary.map { it -> tuple(it[0], it[2], it[1])} //ID, K_Classification, K_percent
                    .combine(bestref_summary, by:0) //ID, K_Classification, K_percent, Type, TypeANI
                    .combine(mlst_summary.map {it[0..1]}, by:0) //ID, K_Classification, K_percent, Type, TypeANI, ST, (removed Alleles, it[2])
                    .combine(fastp_summary.map{it[0..2]}, by:0) //ID, K_Classification, K_percent, Type, TypeANI, ST, Percent>Q30, AvgQ
                    .combine(coverage_summary.map{it -> tuple(it[0], it[1],it[3])}, by:0) //ID, K_Classification, K_percent, Type, TypeANI, ST, Percent>Q30, AvgQ, covX, QC(pass/fail)
                    .combine(amrfinder_summary, by:0) //ID, K_Classification, K_percent, Type, TypeANI, ST, Percent>Q30, AvgQ, covX, QC(pass/fail), AMRFinder_Gene_list
                    .combine(macrolide_summary.map{it[0,2]}, by:0)
                    .combine(tet_summary.map{it[0,2]}, by:0)
                    .combine(quinolone_summary.map{it[0,2]}, by:0)
                    .map{it[0,9,1,2,3,4,5,6,7,8,11,12,13,10]}
                    .collect{ summary -> summary.join('\t')}
                    .set { sampleSummaries }
    run_report = generate_run_report(sampleSummaries)
}