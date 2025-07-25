process vcf_subset_qrdr {
    publishDir "${params.output}/final_vcf", mode: 'copy', pattern: '*.vcf'

    input:
    tuple val(sample), path(vcf), val(qc), path(quinolone_res)

    output:
    tuple val(sample), path("${vcf.simpleName}_ann.vcf"), path("${vcf.simpleName}_quinolone_res_features.vcf"), path("${sample}_qrdr.txt"), val(qc), emit:report
    tuple val(sample), val(qc), env(quinolone_resistance), emit: summary
    
    script:
    """
    ##quinolone aa mutations
    if [ "${qc}" == "PASS" ]; then
        # run snpEff, need to replace reference chromosome name with "Chromosome" used by snpeff
        bcftools view -i 'QUAL>=30' ${vcf} -Ov | sed "s/NC_000912.1/Chromosome/" | snpEff "Mycoplasma_pneumoniae_m129" > ${vcf.simpleName}_ann.vcf
        
        # search the annotated vcf file for specific changes listed in the aa_res file
        while read -r line; do
            gene=\$(echo \$line | awk '{print \$1}')
            aa_change=\$(echo \$line | awk '{print \$3}' | sed -E "s/Xaa//g")
            regex="\${gene}.*p.\${aa_change}"
            grep -E "\${regex}" "${vcf.simpleName}_ann.vcf" >> "${vcf.simpleName}_quinolone_res_features.vcf" || true
        done < ${quinolone_res}
    
        if [ -s "${vcf.simpleName}_quinolone_res_features.vcf" ]; then
            ## Ugly bash code, grabs the INFO column which has the snpeff annotation
            #Splits the snpeff ANN section (42) out
            #Grabs the first part of the annotation & reformats it
            ## \$4 == GENE, \$10 = nucleotide change, \$11 = aa change
            cut -f8 "${vcf.simpleName}_quinolone_res_features.vcf" | cut -d";" -f42 | cut -d"," -f1 | awk -v sample=${sample} 'BEGIN {FS="|"} {if(\$0 ~ /^ANN/ && \$0 ~ /missense_variant/) {print sample, \$4, gensub(/c\\.([0-9]+)([ACTG])>([ATCG])/, "\\\\2\\\\1\\\\3", "g", \$10), gensub(/p\\./,"","g",\$11), "Resistant"} }' > "${sample}_qrdr.txt"
            quinolone_resistance="Resistant"
        else
            echo -e "${sample}\tNA\tNA\tNA\tQuinolone_Sensitive\n" > ${sample}_qrdr.txt
            quinolone_resistance="Susceptible"
        fi
    else 
        echo "sampled failed QC" > ${sample}_qrdr.txt
        quinolone_resistance="Failed_QC"

        #dummy files
        touch "${vcf.simpleName}_ann.vcf"
        touch "${vcf.simpleName}_quinolone_res_features.vcf"
    fi
    """
}