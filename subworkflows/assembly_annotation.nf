include { assembly } from '../modules/unicycler_assembly.nf'
include { mlst } from '../modules/run_mlst.nf'
include { fastANI } from '../modules/run_fastani.nf'
include { bestRef } from '../modules/closest_type_reference.nf'
include { amrfinder } from '../modules/run_amrfinder.nf'

workflow ASSEMBLY_BASED_ANALYSIS {
    take:
    cov_out
    references_ch
    amrfinder_db

    main:
    // assembling the QC and Coverage threshold passed reads
    assemblies = assembly(cov_out)

    //mlst
    mlst = mlst(assemblies.genomes)
    
    //prompt the user how many have passed
    assemblies.genomes
        .filter { tuple -> tuple[-1] == "PASS" }
        .count().subscribe { count ->
            if (count == 0) {
                println "Error: All the input reads failed QC. CAMPNeu will step through the rest of the pipeline using blank files and generate the reports."
            } else {
                println "${count} samples passed QC"
            }
        }

    ref_type = Channel.of(
        ["GCF_000027345", "type1"],
        ["GCF_001272835", "type2"]
    ).combine(references_ch, by:0)

    //run on both passing and failing samples, fastani process will skip failed samples
    inputSet = assemblies.genomes.combine(ref_type)
    
    // getting results of fastANI for all samples compared to both references and grouping based on the samples
    results = fastANI(inputSet)
    grouped = results.fastANI_out.groupTuple()

    // getting the best reference for each isolate/sample 
    bestRef(grouped)

    //give it all samples, including failed so we generate the reports correctly
    ref_type = ref_type.map {it[1..2]}
    samples = assemblies.assembly_out
                    .combine(ref_type.first())
                    //.filter { tuple -> tuple[2] == "PASS"}
                    //.filter { tuple -> tuple[-2] == "type1"}

    // amrfinder
    amrfinder(assemblies.genomes, amrfinder_db)

    emit:
    samples
    amrfinder_report = amrfinder.out.report
    amrfinder_summary = amrfinder.out.summary
    mlst_report = mlst.report
    mlst_summary = mlst.summary
    best_reference = bestRef.out.out
    bestref_report = bestRef.out.report
    bestref_summary = bestRef.out.summary
}