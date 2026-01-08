/*
 * Extractor Workflow
 * Extracts and validates variants from user-provided TSV files
 * against family BCFs, wombat outputs, and individual gVCFs
 */

// Include modules
include { REFACTOR } from '../modules/extractor/refactor'
include { PROCESS_FAM_TSV } from '../modules/extractor/process_fam_tsv'
include { PROCESS_FAM_BCF } from '../modules/extractor/process_fam_bcf'
include { PROCESS_IND_GVCF } from '../modules/extractor/process_ind_gvcf'
include { FAM_AGGREGATE } from '../modules/extractor/fam_aggregate'
include { AGGREGATE } from '../modules/extractor/aggregate'

workflow EXTRACTOR {

    take:
    tsv_files           // channel: [original_filename, tsv_path]
    pedigree            // path to pedigree file
    chain_file          // path to liftover chain file
    norm_bcfs           // channel: [fid, bcf, csi]
    wombat_tsvs         // channel: [fid, tsv.gz]
    gvcfs               // channel: [barcode, gvcf, tbi]

    main:
    
    // Step 1: REFACTOR - split input TSVs into family/sample files
    refactor_input = tsv_files
        .map { original_filename, tsv_path ->
            tuple(original_filename, file(tsv_path), file(pedigree), file(chain_file))
        }
    
    REFACTOR(refactor_input)
    
    // Parse refactor outputs to get family/sample channels
    // Family files: families/{FID}/extractor/{FID}.{original_filename}.tsv
    family_extractor_files = REFACTOR.out.family_list
        .splitText()
        .map { it.trim() }
        .combine(REFACTOR.out.refactored_files)
        .flatMap { fid, original_filename, family_files, sample_files ->
            def family_tsv = family_files.find { it.name.contains("/${fid}.") && it.name.endsWith('.tsv') }
            if (family_tsv) {
                return [[fid, original_filename, family_tsv]]
            }
            return []
        }
    
    // Sample files: samples/{barcode}/extractor/{barcode}.{original_filename}.tsv
    sample_extractor_files = REFACTOR.out.sample_list
        .splitText()
        .map { it.trim() }
        .combine(REFACTOR.out.refactored_files)
        .flatMap { barcode, original_filename, family_files, sample_files ->
            def sample_tsv = sample_files.find { it.name.contains("/${barcode}.") && it.name.endsWith('.tsv') }
            if (sample_tsv) {
                return [[barcode, original_filename, sample_tsv]]
            }
            return []
        }
    
    // Step 2: PROCESS_FAM_TSV - match with wombat bcf2tsv output
    fam_tsv_input = family_extractor_files
        .map { fid, original_filename, extractor_tsv ->
            tuple(fid, original_filename, extractor_tsv)
        }
        .combine(wombat_tsvs, by: 0)  // Join on fid
        .map { fid, original_filename, extractor_tsv, wombat_tsv ->
            tuple(fid, original_filename, extractor_tsv, wombat_tsv)
        }
    
    PROCESS_FAM_TSV(fam_tsv_input)
    
    // Step 3: PROCESS_FAM_BCF - match with normalized BCF
    fam_bcf_input = family_extractor_files
        .map { fid, original_filename, extractor_tsv ->
            tuple(fid, original_filename, extractor_tsv)
        }
        .combine(norm_bcfs, by: 0)  // Join on fid
        .map { fid, original_filename, extractor_tsv, bcf, csi ->
            tuple(fid, original_filename, extractor_tsv, bcf, csi)
        }
    
    PROCESS_FAM_BCF(fam_bcf_input)
    
    // Step 4: PROCESS_IND_GVCF - extract from individual gVCFs
    ind_gvcf_input = sample_extractor_files
        .map { barcode, original_filename, extractor_tsv ->
            tuple(barcode, original_filename, extractor_tsv)
        }
        .combine(gvcfs, by: 0)  // Join on barcode
        .map { barcode, original_filename, extractor_tsv, gvcf, tbi ->
            tuple(barcode, original_filename, extractor_tsv, gvcf, tbi)
        }
    
    PROCESS_IND_GVCF(ind_gvcf_input)
    
    // Step 5: FAM_AGGREGATE - combine results per family
    // Need to join family-level results with sample-level results
    fam_aggregate_input = family_extractor_files
        .map { fid, original_filename, extractor_tsv ->
            tuple(fid, original_filename, extractor_tsv)
        }
        .join(PROCESS_FAM_TSV.out, by: [0, 1])  // fid, original_filename
        .join(PROCESS_FAM_BCF.out, by: [0, 1])  // fid, original_filename
        .map { fid, original_filename, extractor_tsv, fam_tsv_found, fam_tsv_notfound, fam_bcf_found, fam_bcf_notfound ->
            tuple(fid, original_filename, extractor_tsv, fam_tsv_found, fam_tsv_notfound, fam_bcf_found, fam_bcf_notfound)
        }
        .combine(
            PROCESS_IND_GVCF.out
                .map { barcode, original_filename, ind_gvcf ->
                    tuple(original_filename, ind_gvcf)
                }
                .groupTuple(by: 0),
            by: 1  // original_filename
        )
        .map { original_filename, fid, extractor_tsv, fam_tsv_found, fam_tsv_notfound, fam_bcf_found, fam_bcf_notfound, ind_gvcf_files ->
            tuple(fid, original_filename, extractor_tsv, fam_tsv_found, fam_tsv_notfound, fam_bcf_found, fam_bcf_notfound, ind_gvcf_files)
        }
    
    FAM_AGGREGATE(fam_aggregate_input)
    
    // Step 6: AGGREGATE - combine all family summaries into cohort summary
    aggregate_input = FAM_AGGREGATE.out
        .map { fid, original_filename, summary ->
            tuple(original_filename, summary)
        }
        .groupTuple(by: 0)
    
    AGGREGATE(aggregate_input)
    
    emit:
    family_summaries = FAM_AGGREGATE.out
    cohort_summaries = AGGREGATE.out
}
