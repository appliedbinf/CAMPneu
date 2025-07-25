process vcf_subset_16S {
    label 'bcftools'
    publishDir "${params.output}/final_vcf", mode: 'copy', pattern: '*.vcf'

    input:
    tuple val(sample), path(reference), path(vcf), val(qc), val(start), val(end), path(snps_16S)

    output:
    tuple val(sample), path("${vcf.simpleName}_all16S.subset.vcf"), path("${vcf.simpleName}_16S_Tet_identified.snps.vcf"), path("${sample}_16Ssnps.txt"), val(qc), emit: report
    tuple val(sample), val(qc), env(tetracycline_resistance), emit: summary

    script:
    """
    if [ "${qc}" == "PASS" ]; then
        bcftools view -i 'QUAL>=30' ${vcf} -Oz -o ${vcf}.gz
        bcftools index ${vcf}.gz
        chrom=\$(bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' ${vcf}.gz | head -1 | cut -f 1)
        bcftools view -r \${chrom}:${start}-${end} --no-header ${vcf}.gz > ${vcf.simpleName}_all16S.subset.vcf
        bcftools view -R ${snps_16S} --no-header ${vcf}.gz > ${vcf.simpleName}_16S_Tet_identified.snps.vcf

        if [ -s "${vcf.simpleName}_16S_Tet_identified.snps.vcf" ]; then
            touch ${sample}_snps.txt
            cut -f1-2,4-5 "${vcf.simpleName}_16S_Tet_identified.snps.vcf" >> ${sample}_snps.txt
            awk '{ new_col = \$2 - 118313; print "${sample}", \$0, \$3 new_col \$4 }' ${sample}_snps.txt > ${sample}_snps_out.txt
            awk '{\$2=""; print \$0}' ${sample}_snps_out.txt | awk '{print \$0, "Resistant"}' > ${sample}_16Ssnps.txt
            tetracycline_resistance="Resistant"
        else
            touch ${sample}_snps.txt
            echo "NA NA NA NA" | awk '{print \$0, "Sensitive"}' >> ${sample}_snps_1.txt
            awk '{ print "${sample}", \$0 }' ${sample}_snps_1.txt > ${sample}_16Ssnps.txt
            tetracycline_resistance="Susceptible"
        fi
    else 
        echo "sampled failed QC" > ${sample}_16Ssnps.txt
        tetracycline_resistance="Fail_QC"

        ##dummy files
        touch "${vcf.simpleName}_all16S.subset.vcf"
        touch "${vcf.simpleName}_16S_Tet_identified.snps.vcf"
    fi
    """
}