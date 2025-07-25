#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.version = '1.3.0'
params.input = ''
params.output = ''
params.help = false
params.max_cpus = '1'

///// HELP MESSAGE /////
if (params.help) {
        help = """
              |Usage: 
              |CAMPneu.nf --input <fastq_reads_dir> --output <output_dir>
              |       
              |Required arguments:     
              |  --input     Path to the Paired Fastq Reads directory  
              |  --output    Directory where process outputs are saved          
              |Optional arguments:  
              |  --help      Print this message and exit""".stripMargin()


    println(help)
    exit(0)
}

//IMPORTS
include { INITIALIZE_PIPELINE } from './subworkflows/initialize_pipeline.nf'
include { QAQC } from './subworkflows/qaqc_steps.nf'
include { ASSEMBLY_BASED_ANALYSIS } from './subworkflows/assembly_annotation.nf'
include { DETECT_SNPS } from './subworkflows/detect_amr_snps.nf'
include { SAMPLE_SUMMARY } from './subworkflows/generate_sample_summary.nf'
include { RUN_REPORT } from './subworkflows/generate_run_summary.nf'

workflow {

    if (!params.input) {
        error "ERROR: Missing required input parameter. Please specify the input directory using '--input'."
    }
    if (!params.output) {
        error "ERROR: Missing required output parameter. Please specify the output directory using '--output'."
    }

    Channel.fromFilePairs("${params.input}/*_{1,2,R1,R2,r1,r2}*.{fastq,fq,FASTQ,FQ,fastq.gz,fq.gz,FASTQ.GZ,FQ.GZ}")
           .ifEmpty{ error "NO {reads}.fastq/fq files found in the specified directory: ${params.input}"}
           .set {paired_reads}

    //download necessary databases & reference files
    INITIALIZE_PIPELINE( paired_reads )

    //run fastq, kraken, and coverage check vs typ1 reference
    QAQC( 
        INITIALIZE_PIPELINE.out.unzipped_reads,
        INITIALIZE_PIPELINE.out.kraken_db, 
        INITIALIZE_PIPELINE.out.references_ch
        )

    ASSEMBLY_BASED_ANALYSIS(
        QAQC.out.coverage_out,
        INITIALIZE_PIPELINE.out.references_ch,
        INITIALIZE_PIPELINE.out.amrfinder_db
    )

    DETECT_SNPS(
        ASSEMBLY_BASED_ANALYSIS.out.samples,
        channel.fromPath("$projectDir/data/macrolide_resistance_snps_23S.bed"),
        channel.fromPath("$projectDir/data/tetracycline_resistance_snps_16S.bed"),
        channel.fromPath("$projectDir/data/quinolone_resistance_snps_qrdr.bed")
    )

    SAMPLE_SUMMARY(
        QAQC.out.kraken_report,
        QAQC.out.fastp_report,
        QAQC.out.coverage_report,
        ASSEMBLY_BASED_ANALYSIS.out.mlst_report,
        ASSEMBLY_BASED_ANALYSIS.out.bestref_report,
        ASSEMBLY_BASED_ANALYSIS.out.amrfinder_report,
        DETECT_SNPS.out.macrolide_report,
        DETECT_SNPS.out.tet_report,
        DETECT_SNPS.out.quinolone_report
    )

    RUN_REPORT(
        QAQC.out.kraken_summary,
        ASSEMBLY_BASED_ANALYSIS.out.bestref_summary,
        ASSEMBLY_BASED_ANALYSIS.out.mlst_summary,
        QAQC.out.fastp_summary,
        QAQC.out.coverage_summary,
        ASSEMBLY_BASED_ANALYSIS.out.amrfinder_summary,
        DETECT_SNPS.out.macrolide_summary,
        DETECT_SNPS.out.tet_summary,
        DETECT_SNPS.out.quinolone_summary
    )
}
