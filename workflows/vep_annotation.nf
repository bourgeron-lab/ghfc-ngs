/*
 * VEP Annotation Workflow
 * Annotates family VCF files using Ensembl VEP with chunking for parallel processing
 */

// Include modules
include { generateSplits } from '../modules/generate_splits'
include { splitVCF } from '../modules/split_vcf'
include { runVEP } from '../modules/vep_annotation'
include { mergeVCF } from '../modules/merge_vcf'

workflow VEP_ANNOTATION {

    take:
    family_vcfs    // channel: [fid, vcf, tbi]

    main:

    // Prepare VEP input channel
    vep_input = family_vcfs
        .map { fid, vcf, tbi ->
            def meta = [:]
            meta.index_type = "tbi"  // Assuming TBI index
            [meta, fid, vcf, tbi, params.vep_config]
        }

    // Run VEP with chunking
    // Generate split files that each contain bin_size number of variants
    generateSplits(vep_input)

    // Split VCF using split files
    generateSplits.out
        .transpose()
        .set { split_input }

    splitVCF(split_input)

    // Run VEP for each split VCF file
    vep_input_split = splitVCF.out
        .map { meta, original, vcfs, indices, vep_config ->
            // Create individual tuples for each split VCF
            def result = []
            vcfs.eachWithIndex { vcf, i ->
                def split_meta = meta.clone()
                split_meta.split_index = i
                result << [split_meta, original, vcf, indices[i], vep_config]
            }
            result
        }
        .flatten()
        .collate(5)

    runVEP(vep_input_split)

    // Merge split VCF files back into single annotated VCF
    runVEP.out.files
        .groupTuple(by: [0, 1, 4])
        .set { merge_input }

    mergeVCF(merge_input)

    emit:
    annotated_vcfs = mergeVCF.out
}