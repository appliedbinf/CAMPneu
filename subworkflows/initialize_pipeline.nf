include { uncompress_reads } from '../modules/uncompress_fastq.nf'
include { downloadKrakenDB } from "../modules/download_kraken_db.nf"
include { create_23S_bed } from "../modules/create_23S_bed.nf"
include { create_16S_bed } from "../modules/create_16S_bed.nf"
include { create_quinolone_amr_locations } from "../modules/create_qrdr_bed.nf"

workflow INITIALIZE_PIPELINE {
    take:
    paired_reads

    main:
    // Check if minikraken database exists, if not, download to $CONDA_PREFIX
    def kraken_db_dir = file("${params.kraken_db_dir}/minikraken_8GB_202003")
    if (kraken_db_dir.exists()) {
        println "Minikraken database exists, skipping download."
        kraken_db = Channel.value(kraken_db_dir)
    } else {
        println "Minikraken database downloading now..."
        kraken_db = downloadKrakenDB()        
    }

    /// REFERENCE FILES - TYPE 1 AND TYPE 2 /////
    def ref1 = file("${params.reference_dir}/GCF_000027345.1_ASM2734v1_genomic.fna")
    def ref2 = file("${params.reference_dir}/GCF_001272835.1_ASM127283v1_genomic.fna")    
    if (ref1.exists() && ref2.exists()) {
        println "Type1 and Type2 Reference files exist, skipping download"
        references = Channel.fromPath([ref1, ref2])
    } else {
        println "Reference files downloading now..."
        download_refs()
        references = download_refs.out.ref1
                        .concat(download_refs.out.ref2)
    }

    references_ch = references
        .map { 
            file ->
            def id = file.getName().split('\\.1')[0]  // Extract the identifier 
            [id, file]  // Return a tuple of [id, file]
        }
    ///// CREATE SNP FILES /////
    macrolide_file = create_23S_bed()
    tet_file = create_16S_bed()
    quinFile = create_quinolone_amr_locations()

    //Inflate gzipped fastq
    unzipped_reads = uncompress_reads(paired_reads)

    emit:
    unzipped_reads
    references_ch
    kraken_db
    macrolide_file
    tet_file
    quinFile
}