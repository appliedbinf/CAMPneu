process vcf_subset_23S {
    label 'bcftools'
    publishDir "${params.output}/final_vcf", mode: 'copy', pattern: '*.vcf'

    input:
    tuple val(sample), path(reference), path(vcf), val(qc), val(start), val(end), path(snps_23S)

    output:
    tuple val(sample), path("${vcf.simpleName}_all23S.subset.vcf"), path("${vcf.simpleName}_identified.snps.vcf"), path("${sample}_23Ssnps.txt"), val(qc), emit: report
    tuple val(sample), val(qc), env(macrolide_resistance), emit: summary

    script:
    """
    if [ "${qc}" == "PASS" ]; then
        bcftools view -i 'QUAL>=30' ${vcf} -Oz -o ${vcf}.gz
        bcftools index ${vcf}.gz
        chrom=\$(bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' ${vcf}.gz | head -1 | cut -f 1)
        bcftools view -r \${chrom}:${start}-${end} --no-header ${vcf}.gz > ${vcf.simpleName}_all23S.subset.vcf
        bcftools view -R ${snps_23S} --no-header ${vcf}.gz > ${vcf.simpleName}_identified.snps.vcf
    
        if [ -s "${vcf.simpleName}_identified.snps.vcf" ]; then
            touch ${sample}_snps.txt
            cut -f1-2,4-5 "${vcf.simpleName}_identified.snps.vcf" >> ${sample}_snps.txt
            awk '{ new_col = \$2 - 120056; print "${sample}", \$0, \$3 new_col \$4 }' ${sample}_snps.txt > ${sample}_snps_out.txt
            awk '{\$2=""; print \$0}' ${sample}_snps_out.txt | awk '{print \$0, "Resistant"}' > ${sample}_23Ssnps.txt
            macrolide_resistance="Resistant"
        else
            touch ${sample}_snps.txt
            echo "NA NA NA NA" | awk '{print \$0, "Sensitive"}' >> ${sample}_snps_1.txt
            awk '{ print "${sample}", \$0 }' ${sample}_snps_1.txt > ${sample}_23Ssnps.txt
            macrolide_resistance="Susceptible"
        fi
    else
        echo "sampled failed QC" > ${sample}_23Ssnps.txt
        macrolide_resistance="Failed_QC"

        ##dummy files
        touch "${vcf.simpleName}_all23S.subset.vcf"
        touch "${vcf.simpleName}_identified.snps.vcf"
    fi
    """

}