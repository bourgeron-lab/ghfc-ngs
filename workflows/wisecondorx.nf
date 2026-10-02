/*
 * WisecondorX Workflow
 * Converts CRAM files to NPZ format and runs WisecondorX predict for CNV analysis
 * Merges aberration BED files at family and cohort level
 * Annotates family aberrations with gencode gene and exon information
 * 
 * Skips processing if output files already exist - only runs necessary modules
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
    cram_files         // channel: [barcode, cram, crai] - only CRAMs that need NPZ conversion
    existing_npz_files // channel: [barcode, npz] - existing NPZ files
    need_predict       // list of barcodes that need predict
    family_members     // map: [barcode: fid] - mapping of barcodes to family IDs
    need_family_merge  // map: [fid: boolean] - whether family merge is needed
    need_family_annotate // map: [fid: boolean] - whether family annotation is needed
    need_cohort_merge  // boolean: whether cohort merge is needed
    families           // collection: family IDs in the current pedigree - scopes the cohort merge

    main:
    
    // Add reference genome to each input for NPZ conversion
    cram_with_ref = cram_files
        .map { barcode, cram, crai -> 
            tuple(barcode, cram, crai, file(params.ref))
        }
    
    // Run NPZ conversion for samples that need it (cram_files already filtered by main.nf)
    NPZ_CONVERT(cram_with_ref)
    
    // Mix existing NPZ files with newly converted ones
    all_npz_files = existing_npz_files
        .mix(NPZ_CONVERT.out)
    
    // Check which samples need predict vs have existing results
    npz_with_status = all_npz_files
        .filter { barcode, _npz ->
            barcode in need_predict
        }
        .map { barcode, npz ->
            def smp_dir = Sharding.getSampleDir(params.data, barcode)
            def predict_bed = file("${smp_dir}/svs/wisecondorx/${barcode}_aberrations.bed")
            def chr_bed = file("${smp_dir}/svs/wisecondorx/${barcode}_aberrations.chr.bed")

            def has_predict = predict_bed.exists()
            def has_chr = chr_bed.exists()
            
            [barcode: barcode, npz: npz, 
             predict_bed: predict_bed, chr_bed: chr_bed,
             has_predict: has_predict, has_chr: has_chr]
        }
    
    // Samples needing PREDICT
    npz_for_predict = npz_with_status
        .filter { rec -> !rec.has_predict }
        .map { rec -> tuple(rec.barcode, rec.npz, file(params.wisecondorx_reference.replace('{binsize}', params.wisecondorx_binsize.toString()))) }
    
    // Existing predict results
    existing_predict_beds = npz_with_status
        .filter { rec -> rec.has_predict }
        .map { rec -> tuple(rec.barcode, rec.predict_bed) }
    
    // Run WisecondorX predict
    PREDICT(npz_for_predict)
    
    // Extract aberrations BED files from predict output
    new_aberrations_beds = PREDICT.out
        .map { barcode, bed, _stats, _segments, _bins ->
            tuple(barcode, bed)
        }
    
    // Combine new + existing predict results
    all_aberrations_beds = new_aberrations_beds.mix(existing_predict_beds)
    
    // Check which samples need chr reformat
    aberrations_with_chr_status = all_aberrations_beds
        .map { barcode, bed ->
            def chr_bed = file("${Sharding.getSampleDir(params.data, barcode)}/svs/wisecondorx/${barcode}_aberrations.chr.bed")
            def has_chr = chr_bed.exists()
            [barcode: barcode, bed: bed, chr_bed: chr_bed, has_chr: has_chr]
        }
    
    // Samples needing chr reformat
    beds_for_chr_reformat = aberrations_with_chr_status
        .filter { rec -> !rec.has_chr }
        .map { rec -> tuple(rec.barcode, rec.bed) }
    
    // Existing chr-reformatted beds
    existing_chr_beds = aberrations_with_chr_status
        .filter { rec -> rec.has_chr }
        .map { rec -> tuple(rec.barcode, rec.chr_bed) }
    
    // Reformat chromosome names to add chr prefix
    REFORMAT_CHR(beds_for_chr_reformat)
    
    // Combine new + existing chr-reformatted results
    all_chr_aberrations = REFORMAT_CHR.out.chr_aberrations.mix(existing_chr_beds)
    
    // Merge family aberrations if needed
    if (need_family_merge && !need_family_merge.isEmpty()) {
        // Every existing chr-prefixed aberrations file of a member of a family to be merged, not
        // just the newly created ones. Looked up at the sharded path of each such member rather
        // than globbed across the whole samples/ tree: the merge below keeps only those families
        // anyway, so the glob found nothing more, just slowly.
        all_existing_chr_aberrations = channel.fromList(
            family_members
                .findAll { _barcode, fid -> need_family_merge[fid] == true }
                .collect { barcode, _fid ->
                    tuple(barcode, file("${Sharding.getSampleDir(params.data, barcode)}/svs/wisecondorx/${barcode}_aberrations.chr.bed"))
                }
                .findAll { _barcode, bed -> bed.exists() })
        
        // Mix ALL existing (from disk) with newly created chr-prefixed aberrations.
        // A sample reformatted in this run also matches the disk glob once published, so dedupe by
        // barcode - staging both copies would make merge_family's `ls` miss the renamed duplicate.
        combined_chr_aberrations = REFORMAT_CHR.out.chr_aberrations
            .mix(all_existing_chr_aberrations)
            .unique { barcode, _bed -> barcode }

        // Group by family ID
        family_aberrations = combined_chr_aberrations
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
        // Existing family aberrations of the families to annotate, by path rather than by glob
        existing_family_aberrations_for_annot = channel.fromList(
            need_family_annotate
                .findAll { _fid, needed -> needed == true }
                .collect { fid, _needed ->
                    tuple(fid, file("${Sharding.getFamilyDir(params.data, fid)}/svs/wisecondorx/${fid}_aberrations.bed"))
                }
                .findAll { _fid, bed -> bed.exists() })
        
        // Mix with newly created family aberrations
        all_family_aberrations_for_annot = existing_family_aberrations_for_annot.mix(family_output)
        
        // Filter for families that need annotation
        families_to_annotate = all_family_aberrations_for_annot
            .filter { fid, _bed ->
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
        // Existing annotated family aberrations of the families in the current pedigree, and
        // only those: every family ever processed under params.data must never be merged into
        // this cohort. Looked up per family rather than globbed and filtered.
        existing_annotated_family_aberrations = channel.fromList(
            families
                .collect { fid -> file("${Sharding.getFamilyDir(params.data, fid)}/svs/wisecondorx/${fid}_aberrations.annotated.bed") }
                .findAll { bed -> bed.exists() })

        // Concat rather than collect-then-mix so the barrier genuinely waits for this run's own
        // ANNOTATE_ABERRATIONS output before the cohort is merged. Dedupe by family in case a
        // freshly annotated family has already been published to disk.
        cohort_input = annotated_output
            .map { _fid, bed -> bed }
            .concat(existing_annotated_family_aberrations)
            .unique { bed -> bed.name }
            .collect()
            .filter { bed_files ->
                // Never publish a cohort built from an incomplete set of families
                def expected = families.findAll { fid ->
                    need_family_annotate[fid] == true ||
                    file("${Sharding.getFamilyDir(params.data, fid)}/svs/wisecondorx/${fid}_aberrations.annotated.bed").exists()
                }
                def found = bed_files.collect { bed -> bed.name.replaceAll(/_aberrations\.annotated\.bed$/, '') }
                def missing = expected - found
                if (missing) {
                    log.warn "Skipping cohort aberrations merge: annotated aberrations missing for ${missing.join(', ')}"
                    return false
                }
                return true
            }
            .map { bed_files -> tuple(params.cohort_name, bed_files) }

        MERGE_COHORT_ABERRATIONS(cohort_input)
        cohort_output = MERGE_COHORT_ABERRATIONS.out.cohort_aberrations_bed
    } else {
        cohort_output = channel.empty()
    }
    
    emit:
    npz_files = NPZ_CONVERT.out
    predict_results = PREDICT.out
    chr_aberrations = all_chr_aberrations
    family_aberrations = family_output
    annotated_family_aberrations = annotated_output
    cohort_aberrations = cohort_output
}
