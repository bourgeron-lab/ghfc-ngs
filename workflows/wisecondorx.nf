/*
 * WisecondorX Workflow
 * Converts CRAM files to NPZ format and runs WisecondorX predict for CNV analysis
 * Merges aberration BED files at family and cohort level
 * Annotates family aberrations with gencode gene and exon information
 */

// Include modules
include { NPZ_CONVERT } from '../modules/wisecondorx/npz_convert'
include { PREDICT } from '../modules/wisecondorx/predict'
include { REFORMAT_CHR } from '../modules/wisecondorx/reformat_chr'
include { MERGE_FAMILY_ABERRATIONS } from '../modules/wisecondorx/merge_family'
include { ANNOTATE_ABERRATIONS } from '../modules/wisecondorx/annotate_aberrations'
include { MERGE_COHORT_ABERRATIONS } from '../modules/wisecondorx/merge_cohort'

workflow WISECONDORX {

    take:
    cram_files         // channel: [barcode, cram, crai]
    existing_npz_files // channel: [barcode, npz]
    need_predict       // list of barcodes that need predict
    family_members     // map: [barcode: fid] - mapping of barcodes to family IDs
    need_family_merge  // map: [fid: boolean] - whether family merge is needed
    need_family_annotate // map: [fid: boolean] - whether family annotation is needed
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
        .map { barcode, bed, _stats, _segments, _bins, _plots ->
            tuple(barcode, bed)
        }
    
    // Reformat chromosome names to add chr prefix
    REFORMAT_CHR(aberrations_beds)
    
    // Merge family aberrations if needed
    if (need_family_merge && !need_family_merge.isEmpty()) {
        // Create channel for existing individual chr-prefixed aberrations
        existing_chr_aberrations = channel
            .fromPath("${params.data}/samples/*/svs/wisecondorx/*_aberrations.chr.bed")
            .map { bed ->
                def barcode = bed.name.replaceAll(/_aberrations\.chr\.bed$/, '')
                tuple(barcode, bed)
            }
        
        // Mix existing with newly created chr-prefixed aberrations
        all_chr_aberrations = existing_chr_aberrations.mix(REFORMAT_CHR.out.chr_aberrations)
        
        // Group by family ID
        family_aberrations = all_chr_aberrations
            .map { barcode, bed ->
                def fid = family_members[barcode]
                tuple(fid, barcode, bed)
            }
            .groupTuple(by: 0)
            .filter { fid, _barcode_list, _bed_list ->
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
    
    // Annotate family aberrations if needed
    if (need_family_annotate && !need_family_annotate.isEmpty()) {
        // Create channel for existing family aberrations
        existing_family_aberrations_for_annot = channel
            .fromPath("${params.data}/families/*/svs/wisecondorx/*_aberrations.bed")
            .map { bed ->
                def fid = bed.name.replaceAll(/_aberrations\.bed$/, '')
                tuple(fid, bed)
            }
        
        // Mix with newly created family aberrations
        all_family_aberrations_for_annot = existing_family_aberrations_for_annot.mix(family_output)
        
        // Filter for families that need annotation
        families_to_annotate = all_family_aberrations_for_annot
            .filter { fid, bed ->
                need_family_annotate[fid] == true
            }
            .map { fid, bed ->
                tuple(fid, bed, params.annotation_annotation_path, params.annotation_gencode)
            }
        
        ANNOTATE_ABERRATIONS(families_to_annotate)
        annotated_output = ANNOTATE_ABERRATIONS.out.annotated_aberrations
    } else {
        annotated_output = channel.empty()
    }
    
    // Merge cohort aberrations if needed
    if (need_cohort_merge) {
        // Create channel for existing annotated family aberrations
        existing_annotated_family_aberrations = channel
            .fromPath("${params.data}/families/*/svs/wisecondorx/*_aberrations.annotated.bed")
            .collect()
        
        // Mix with newly annotated family aberrations if any
        if (need_family_annotate && !need_family_annotate.isEmpty()) {
            all_annotated_family_aberrations = existing_annotated_family_aberrations
                .mix(annotated_output.map { fid, bed -> bed }.collect())
                .flatten()
                .collect()
        } else {
            all_annotated_family_aberrations = existing_annotated_family_aberrations
        }
        
        // Create cohort merge input
        cohort_input = all_annotated_family_aberrations
            .map { bed_files -> tuple(params.cohort_name, bed_files) }
        
        MERGE_COHORT_ABERRATIONS(cohort_input)
        cohort_output = MERGE_COHORT_ABERRATIONS.out.cohort_aberrations_bed
    } else {
        cohort_output = channel.empty()
    }
    
    emit:
    npz_files = NPZ_CONVERT.out
    predict_results = PREDICT.out
    chr_aberrations = REFORMAT_CHR.out.chr_aberrations
    family_aberrations = family_output
    annotated_family_aberrations = annotated_output
    cohort_aberrations = cohort_output
}
