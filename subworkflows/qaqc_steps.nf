include { kraken } from '../modules/run_kraken2.nf'
include { fastp } from '../modules/run_fastp.nf'
include { summarize_fastp } from '../modules/summarize_fastp.nf'
include { summarize_reference_coverage } from '../modules/summarize_ref_coverage.nf'

workflow QAQC {
    take:
    unzipped_reads
    kraken_db
    references

    main:
    // Run Kraken and generate Kraken classification and report 
    kraken_input = unzipped_reads.combine(kraken_db)
    kraken_run = kraken(kraken_input)

    // Run fastp to filter reads based on Q-scores and assign QC value of PASS or FAIL
    fastp_json = fastp(kraken_run.out)
    fastp_run = summarize_fastp(fastp_json)

    //Run the coverage check
    cov_input = fastp_run.out
                .combine(references.first())
    cov_check = summarize_reference_coverage(cov_input)

    emit:
    kraken_summary = kraken_run.summary
    kraken_report = kraken_run.report
    fastp_report = fastp_run.report
    fastp_summary = fastp_run.summary
    coverage_out = cov_check.out
    coverage_report = cov_check.report
    coverage_summary = cov_check.summary
}