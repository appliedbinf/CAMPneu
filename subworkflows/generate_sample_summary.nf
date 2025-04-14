include { generate_sample_report } from '../modules/concat_all_reports.nf'

workflow SAMPLE_SUMMARY {
    take:
    kraken_report
    fastp_report
    cov_report
    mlst_report
    best_ref
    amrfinder_report
    macrolide_report
    tet_report
    quinolone_report

    main:
    // PER SAMPLE SUMMARIES
    combined_sample_outputs = kraken_report
            .combine(fastp_report, by:0)
            .combine(cov_report, by:0) 
            .combine(mlst_report, by:0)
            .combine(best_ref, by:0)
            .combine(amrfinder_report, by:0)
            .combine(macrolide_report.map{it[0,3]}, by:0)
            .combine(tet_report.map{it[0,3]}, by:0)
            .combine(quinolone_report.map{it[0,3]}, by:0)
  
    generate_sample_report(combined_sample_outputs)
}