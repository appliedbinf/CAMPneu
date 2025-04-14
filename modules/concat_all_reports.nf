process generate_sample_report{

    publishDir "${params.output}/sample_reports", mode: 'copy'
    
    input:
    tuple val(sample), path(kraken), path(fastp), path(coverage), path(mlst), path(bestRef), path(amrfinder), path(snps_23S), path(snps_16S), path(quin_res_vcf)

    output:
    path("${sample}_report.out")
    
    script:
    """
    touch ${sample}_report.out
    echo "Kraken Classfication\n" >> ${sample}_report.out
    cat ${kraken} >> ${sample}_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> ${sample}_report.out
    echo "Quality Filtering using Fastp\nSamples that have Qscore < 30 are marked as FAILED\n" >> ${sample}_report.out
    cat ${fastp} >> ${sample}_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> ${sample}_report.out
    echo "Coverage Filtering using Samtools Coverage\nSamples below 10x are marked as FAILED\n" >> ${sample}_report.out
    cat ${coverage} >> ${sample}_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> ${sample}_report.out
    echo "Sequence typing using MLST\n" >> ${sample}_report.out
    cat ${mlst} >> ${sample}_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> ${sample}_report.out
    echo "FASTani to select the best reference\n" >> ${sample}_report.out
    cat ${amrfinder} >> ${sample}_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> ${sample}_report.out
    echo "AMRFinder results\n" >> ${sample}_report.out
    cat ${bestRef} >> ${sample}_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> ${sample}_report.out
    echo "Identification of Macrolide Resistant SNPs using Freebayes and bcftools" >> ${sample}_report.out
    echo -e "Sample\tPos\tALT\tREF\tSNP\tType" | cat - ${snps_23S} > temp && mv temp ${snps_23S}
    cat ${snps_23S} >> ${sample}_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> ${sample}_report.out
    echo "Experimental Identification of Tetracycline Resistant SNPs using Freebayes and bcftools" >> ${sample}_report.out
    echo -e "Sample\tPos\tALT\tREF\tSNP\tType" | cat - ${snps_16S} > temp && mv temp ${snps_16S}
    cat ${snps_16S} >> ${sample}_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> ${sample}_report.out
    echo "Experimental Identification of Quinolone Resistant SNPs using Freebayes, snpEff, and bcftools" >> ${sample}_report.out
    echo -e "Sample\tGene\tNucleotide\tAminoAcid\tType" | cat - ${quin_res_vcf} > temp && mv temp ${quin_res_vcf}
    cat ${quin_res_vcf} >> ${sample}_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> ${sample}_report.out
    """
}