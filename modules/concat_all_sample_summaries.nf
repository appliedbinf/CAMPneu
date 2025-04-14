process generate_run_report {
    publishDir "${params.output}", mode: 'copy'
    input:
    val(sampleSummaries) // tuple( ID, QC(pass/fail), K_Classification, K_percent, Type, TypeANI, ST, Percent>Q30, AvgQ, covX, AMRFinder_Gene_list )

    output:
    path("CAMPneu_report_*.csv")

    script:
    //script info
    def scriptName = workflow.scriptName.replace('.nf','')
    def version = params.version
    def runDate = new Date().format('yyyy-MM-dd')
    // header and data
    def row = sampleSummaries.join("\n")
    def csv_header = "SampleID,QC Status,Kraken Classification,Kraken Percent,Closest Mp type strain,ANI to Mp type strain,MLST ST,Fraction Bases Q30+,Avg QScore,Coverage,Predicted Macrolide Resistance,Predicted Tetracycline Resistance,Predicted Quinolone Resistance,AMRFinder Detected Genes"
    def csv_row = row.replace("\t", ",")

    """  
    echo -e "${scriptName}\nComprehensive Analysis of Mycoplasma Pneumoniae\nVersion,${version}\nDate,${runDate}\n" > final_summary.csv
    echo ${csv_header} >> final_summary.csv
    echo '${csv_row}' >> final_summary.csv
    echo -e "\n*Reads below a qscore of 30 are marked as FAILED and dropped. Reads with coverage below 10X are marked as failed and dropped\n" >> final_summary.csv
    mv final_summary.csv "CAMPneu_report_${runDate}.csv"
    """
}