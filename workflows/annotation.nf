/*
 * Annotation Workflow
 * Integrated workflow for variant annotation and filtering:
 * 1. Annotate with gnomAD frequencies
 * 2. Filter into rare and common variants
 * 3. Filter common variants (keep only GT)
 * 4. Annotate rare variants with VEP
 */

// Include modules
include { GNOMAD_FREQ_ANNOT } from '../modules/annotation/gnomad_freq_annot'
include { GNOMAD_FREQ_FILTER } from '../modules/annotation/gnomad_freq_filter'
include { COMMON_FILTERS } from '../modules/annotation/common_filters'
include { VEP_ANNOTATION } from '../modules/annotation/vep_annotation'
include { BCFTOOLS_ANNOTATE } from '../modules/annotation/bcftools_annotate'
include { DNM_EXTRACTION } from '../modules/annotation/dnm_extraction'

workflow ANNOTATION {

    take:
    normalized_bcfs    // channel: [fid, bcf, csi]
    family_pedigrees   // channel: [fid, pedigree_file]

    main:

    // Step 1: Annotate with gnomAD frequencies
    gnomad_input = normalized_bcfs
        .map { fid, bcf, csi ->
            [fid, bcf, csi, params.gnomad_file]
        }

    GNOMAD_FREQ_ANNOT(gnomad_input)

    // Step 2: Filter into rare and common variants
    filter_input = GNOMAD_FREQ_ANNOT.out.annotated_bcfs
        .map { fid, bcf, csi ->
            [fid, bcf, csi, params.gnomad_filter_field, params.gnomad_filter_threshold]
        }

    GNOMAD_FREQ_FILTER(filter_input)

    // Step 3: Filter common variants (keep only GT field)
    COMMON_FILTERS(GNOMAD_FREQ_FILTER.out.common_bcfs)

    // Step 4: Annotate rare variants with VEP
    // First, separate rare VCFs into those needing VEP vs those with existing VEP results
    rare_vcfs_with_check = GNOMAD_FREQ_FILTER.out.rare_vcfs
        .map { fid, vcf, tbi ->
            def vep_vcf = file("${params.data}/families/${fid}/vcfs/${fid}.rare.${params.vep_config_name}.vcf.gz")
            def vep_tbi = file("${params.data}/families/${fid}/vcfs/${fid}.rare.${params.vep_config_name}.vcf.gz.tbi")
            
            // Check if VEP results already exist
            if (vep_vcf.exists() && vep_tbi.exists()) {
                return [fid: fid, vcf: vcf, tbi: tbi, vep_vcf: vep_vcf, vep_tbi: vep_tbi, has_vep: true]
            } else {
                return [fid: fid, vcf: vcf, tbi: tbi, vep_vcf: null, vep_tbi: null, has_vep: false]
            }
        }
    
    // Split into two channels: those needing VEP and those with existing VEP
    vep_needed = rare_vcfs_with_check
        .filter { it.has_vep == false }
        .map { [it.fid, it.vcf, it.tbi, params.vep_config] }
    
    vep_existing = rare_vcfs_with_check
        .filter { it.has_vep == true }
        .map { [it.fid, it.vep_vcf, it.vep_tbi] }

    // Run VEP only for families that need it
    VEP_ANNOTATION(vep_needed)

    // Combine new VEP outputs with existing VEP files
    all_vep_vcfs = VEP_ANNOTATION.out.annotated_vcfs.mix(vep_existing)
    
    // Step 5: Add additional annotations from BCF files
    // Check for existing bcftools annotated files
    vep_with_bcftools_check = all_vep_vcfs
        .map { fid, vcf, tbi ->
            def annotated_bcf = file("${params.data}/families/${fid}/vcfs/${fid}.rare.${params.vep_config_name}.annotated.bcf")
            def annotated_csi = file("${params.data}/families/${fid}/vcfs/${fid}.rare.${params.vep_config_name}.annotated.bcf.csi")
            
            // Check if bcftools annotated results already exist
            if (annotated_bcf.exists() && annotated_csi.exists()) {
                return [fid: fid, vcf: vcf, tbi: tbi, annotated_bcf: annotated_bcf, annotated_csi: annotated_csi, has_bcftools: true]
            } else {
                return [fid: fid, vcf: vcf, tbi: tbi, annotated_bcf: null, annotated_csi: null, has_bcftools: false]
            }
        }
    
    // Split into two channels: those needing bcftools annotation and those with existing results
    bcftools_needed = vep_with_bcftools_check
        .filter { it.has_bcftools == false }
        .map { [it.fid, it.vcf, it.tbi, params.annotation_annotation_path, params.annotation_annotation_list] }
    
    bcftools_existing = vep_with_bcftools_check
        .filter { it.has_bcftools == true }
        .map { [it.fid, it.annotated_bcf, it.annotated_csi] }

    // Run bcftools annotate only for families that need it
    BCFTOOLS_ANNOTATE(bcftools_needed)

    // Combine new bcftools outputs with existing annotated BCFs
    all_annotated_bcfs = BCFTOOLS_ANNOTATE.out.annotated_bcfs.mix(bcftools_existing)

    // Step 6: Extract de novo mutations
    dnm_input = all_annotated_bcfs
        .join(family_pedigrees)
        .map { fid, bcf, csi, pedigree ->
            [fid, bcf, csi, pedigree, 
             params.annotation_dnm_min_callrate, 
             params.annotation_dnm_min_DP, 
             params.annotation_dnm_min_GQ, 
             params.annotation_dnm_min_VAF,
             params.ref_par1_start,
             params.ref_par1_end,
             params.ref_par2_start,
             params.ref_par2_end,
             params.annotation_annotation_path]
        }

    DNM_EXTRACTION(dnm_input)

    emit:
    gnomad_bcfs = GNOMAD_FREQ_ANNOT.out.annotated_bcfs
    rare_vcfs = GNOMAD_FREQ_FILTER.out.rare_vcfs
    common_bcfs = GNOMAD_FREQ_FILTER.out.common_bcfs
    filtered_common_bcfs = COMMON_FILTERS.out.filtered_common_bcfs
    vep_annotated_vcfs = all_vep_vcfs
    fully_annotated_bcfs = all_annotated_bcfs
    dnm_bcfs = DNM_EXTRACTION.out.dnm_results
}
