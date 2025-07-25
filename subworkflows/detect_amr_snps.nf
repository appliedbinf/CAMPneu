include { minimap2 } from '../modules/run_minimap2.nf'
include { freebayes } from '../modules/run_freebayes.nf'
include { vcf_subset_23S } from '../modules/macrolide_23S_snps.nf'
include { vcf_subset_16S } from '../modules/tetracycline_16S_snps.nf'
include { vcf_subset_qrdr } from '../modules/quinolone_qrdr_snps.nf'

workflow DETECT_SNPS {
    take:
    samples
    bed_23S //from file in data
    bed_16S //from file in data
    quinFile //from file in data
    
    main:
    minimapOut = minimap2(samples)
    freebayesOut = freebayes(minimapOut)

    //23S variants conferring macrolide resistance
    inputVcf_23S = Channel.of(["120057", "122961"])
                .combine(bed_23S)
    inputVcf_23S = freebayesOut.combine(inputVcf_23S)
    macrolide_res_23S = vcf_subset_23S(inputVcf_23S)

    // 16S variants conferring tetracycline resistance
    inputVcf_16S = Channel.of(["118314", "119829"])
                .combine(bed_16S)
    inputVcf_16S = freebayesOut.combine(inputVcf_16S)
    tetracycline_res_16S = vcf_subset_16S(inputVcf_16S)

    // QRDR conferring variants
    quin_search_in = inputVcf_23S.map {
        it -> tuple it[0], it[2], it[3] //sample, reference, vcf_in, qc
        }
        .combine(quinFile)
    quinolone_res = vcf_subset_qrdr(quin_search_in)

    emit:
    macrolide_report = macrolide_res_23S.report
    macrolide_summary = macrolide_res_23S.summary
    tet_report = tetracycline_res_16S.report
    tet_summary = tetracycline_res_16S.summary
    quinolone_report = quinolone_res.report
    quinolone_summary = quinolone_res.summary
}