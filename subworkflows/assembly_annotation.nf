include { assembly } from '../modules/unicycler_assembly.nf'
include { mlst } from '../modules/run_mlst.nf'
include { fastANI } from '../modules/run_fastani.nf'
include { bestRef } from '../modules/closest_type_reference.nf'
include { amrfinder } from '../modules/run_amrfinder.nf'

workflow ASSEMBLY_BASED_ANALYSIS {
    take:
    cov_out
    references_ch

    main:
    // assembling the QC and Coverage threshold passed reads
    assemblies = assembly(cov_out)

    //mlst
    mlst = mlst(assemblies.genomes)
    
    passed_samples = assemblies.genomes
                    .filter { tuple -> tuple[-1] == "PASS" }

    passed_samples.map{ it[0] }
                  .collect().subscribe {samples ->
                        println samples.isEmpty() ? "Error: All the input reads failed QC. Terminating..." : "${samples.size()} samples passed QC."    
                }
    
    ref_type = Channel.of(
        ["GCF_000027345", "type1"],
        ["GCF_001272835", "type2"]
    ).combine(references_ch, by:0)

    inputSet = assemblies.genomes.combine(ref_type)
    

    // getting results of fastANI for all samples compared to both references and grouping based on the samples
    results = fastANI(inputSet)
    grouped = results.fastANI_out.groupTuple()

    // getting the best reference for each isolate/sample 
    bestRef(grouped)

    ref_type = ref_type.map {it[1..2]}
    
    passed_samples = assemblies.assembly_out
                    .combine(ref_type)
                    .filter { tuple -> tuple[-2] == "type1"}

    // amrfinder
    amrfinder(assemblies.genomes)

    emit:
    passed_samples
    amrfinder_report = amrfinder.out.report
    amrfinder_summary = amrfinder.out.summary
    mlst_report = mlst.report
    mlst_summary = mlst.summary
    best_reference = bestRef.out.out
    bestref_report = bestRef.out.report
    bestref_summary = bestRef.out.summary
}