/*
 * WisecondorX Workflow
 * Converts CRAM files to NPZ format and runs WisecondorX predict for CNV analysis
 * Merges aberration BED files at family and cohort level
 */

// Include modules
include { NPZ_CONVERT } from '../modules/wisecondorx/npz_convert'
include { PREDICT } from '../modules/wisecondorx/predict'
include { MERGE_FAMILY_ABERRATIONS } from '../modules/wisecondorx/merge_family'
include { MERGE_COHORT_ABERRATIONS } from '../modules/wisecondorx/merge_cohort'

workflow WISECONDORX {

    take:
    cram_files         // channel: [barcode, cram, crai]
    existing_npz_files // channel: [barcode, npz]
    need_predict       // list of barcodes that need predict
    family_members     // map: [barcode: fid] - mapping of barcodes to family IDs
    need_family_merge  // map: [fid: boolean] - whether family merge is needed
    need_cohort_merge  // boolean: whether cohort merge is needed

    main:
    
    // Add reference genome to each input for NPZ conversion
    cram_with_ref = cram_files
        .map { barcode, cram, crai -> 
            tuple(barcode, cram, crai, file(params.ref))
        }
    
    // Run NPZ conversion for samples that need it
    NPZ_CONVERT(cram_with_ref)
    
    // Mix existing NPZ files with newly converted ones
    all_npz_files = existing_npz_files
        .mix(NPZ_CONVERT.out)
    
    // Filter NPZ files for samples that need predict
    npz_for_predict = all_npz_files
        .filter { barcode, npz ->
            barcode in need_predict
        }
        .map { barcode, npz ->
            tuple(barcode, npz, file(params.wisecondorx_reference))
        }
    
    // Run WisecondorX predict
    PREDICT(npz_for_predict)
    
    // Extract aberrations BED files from predict output
    aberrations_beds = PREDICT.out
        .map { barcode, bed, stats, segments, bins, plots ->
            tuple(barcode, bed)
        }
    
    // Merge family aberrations if needed
    if (need_family_merge && !need_family_merge.isEmpty()) {
        // Create channel for existing individual aberrations
        existing_aberrations = channel
            .fromPath("${params.data}/samples/*/svs/wisecondorx/*_aberrations.bed")
            .map { bed ->
                def barcode = bed.name.replaceAll(/_aberrations\.bed$/, '')
                tuple(barcode, bed)
            }
        
        // Mix existing with newly created aberrations
        all_aberrations = existing_aberrations.mix(aberrations_beds)
        
        // Group by family ID
        family_aberrations = all_aberrations
            .map { barcode, bed ->
                def fid = family_members[barcode]
                tuple(fid, barcode, bed)
            }
            .groupTuple(by: 0)
            .filter { fid, barcode_list, bed_list ->
                // Only process families that need merging
                need_family_merge[fid] == true
            }
            .map { fid, barcode_list, bed_list ->
                tuple(fid, barcode_list, bed_list)
            }
        
        MERGE_FAMILY_ABERRATIONS(family_aberrations)
        family_output = MERGE_FAMILY_ABERRATIONS.out.family_aberrations_bed
    } else {
        family_output = channel.empty()
    }
    
    // Merge cohort aberrations if needed
    if (need_cohort_merge) {
        // Create channel for existing family aberrations
        existing_family_aberrations = channel
            .fromPath("${params.data}/families/*/svs/wisecondorx/*_aberrations.bed")
            .collect()
        
        // Mix with newly created family aberrations if any
        if (need_family_merge && !need_family_merge.isEmpty()) {
            all_family_aberrations = existing_family_aberrations
                .mix(family_output.map { fid, bed -> bed }.collect())
                .flatten()
                .collect()
        } else {
            all_family_aberrations = existing_family_aberrations
        }
        
        // Create cohort merge input
        cohort_input = all_family_aberrations
            .map { bed_files -> tuple(params.cohort_name, bed_files) }
        
        MERGE_COHORT_ABERRATIONS(cohort_input)
        cohort_output = MERGE_COHORT_ABERRATIONS.out.cohort_aberrations_bed
    } else {
        cohort_output = channel.empty()
    }
    
    emit:
    npz_files = NPZ_CONVERT.out
    predict_results = PREDICT.out
    family_aberrations = family_output
    cohort_aberrations = cohort_output
}
