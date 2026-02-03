/*
 * Annotation Workflow
 * Integrated workflow for variant annotation and filtering:
 * 1. Annotate with gnomAD frequencies
 * 2. Filter into rare and common variants
 * 3. Filter common variants (keep only GT)
 * 4. Annotate rare variants with VEP
 * 5. Add additional annotations with bcftools
 * 
 * Skips processing if output files already exist - only runs necessary modules
 */

// Include modules
include { GNOMAD_FREQ_ANNOT } from '../modules/annotation/gnomad_freq_annot'
include { GNOMAD_FREQ_FILTER } from '../modules/annotation/gnomad_freq_filter'
include { COMMON_FILTERS } from '../modules/annotation/common_filters'
include { VEP_ANNOTATION } from '../modules/annotation/vep_annotation'
include { BCFTOOLS_ANNOTATE } from '../modules/annotation/bcftools_annotate'

workflow ANNOTATION {

    take:
    normalized_bcfs    // channel: [fid, bcf, csi]
    _family_pedigrees  // channel: [fid, pedigree_file] - kept for API compatibility

    main:

    // Check existing outputs for each family to determine what needs to run
    bcfs_with_status = normalized_bcfs
        .map { fid, bcf, csi ->
            def rare_vcf = file("${params.data}/families/${fid}/vcfs/${fid}.rare.vcf.gz")
            def rare_tbi = file("${params.data}/families/${fid}/vcfs/${fid}.rare.vcf.gz.tbi")
            def common_bcf = file("${params.data}/families/${fid}/vcfs/${fid}.common.bcf")
            def common_csi = file("${params.data}/families/${fid}/vcfs/${fid}.common.bcf.csi")
            def common_gt_bcf = file("${params.data}/families/${fid}/vcfs/${fid}.common_gt.bcf")
            def common_gt_csi = file("${params.data}/families/${fid}/vcfs/${fid}.common_gt.bcf.csi")
            def vep_vcf = file("${params.data}/families/${fid}/vcfs/${fid}.rare.${params.vep_config_name}.vcf.gz")
            def vep_tbi = file("${params.data}/families/${fid}/vcfs/${fid}.rare.${params.vep_config_name}.vcf.gz.tbi")
            def annotated_bcf = file("${params.data}/families/${fid}/vcfs/${fid}.rare.${params.vep_config_name}.annotated.bcf")
            def annotated_csi = file("${params.data}/families/${fid}/vcfs/${fid}.rare.${params.vep_config_name}.annotated.bcf.csi")
            
            def has_rare_common = rare_vcf.exists() && rare_tbi.exists() && 
                                  common_bcf.exists() && common_csi.exists()
            def has_common_gt = common_gt_bcf.exists() && common_gt_csi.exists()
            def has_vep = vep_vcf.exists() && vep_tbi.exists()
            def has_annotated = annotated_bcf.exists() && annotated_csi.exists()
            
            [fid: fid, bcf: bcf, csi: csi,
             rare_vcf: rare_vcf, rare_tbi: rare_tbi,
             common_bcf: common_bcf, common_csi: common_csi,
             common_gt_bcf: common_gt_bcf, common_gt_csi: common_gt_csi,
             vep_vcf: vep_vcf, vep_tbi: vep_tbi,
             annotated_bcf: annotated_bcf, annotated_csi: annotated_csi,
             has_rare_common: has_rare_common, has_common_gt: has_common_gt,
             has_vep: has_vep, has_annotated: has_annotated]
        }

    // =====================================================================
    // Step 1 & 2: gnomAD annotation and filtering (rare/common split)
    // =====================================================================
    
    // Families needing gnomAD annotation + filtering (no rare/common files exist)
    bcfs_need_gnomad = bcfs_with_status
        .filter { rec -> !rec.has_rare_common }
        .map { rec -> [rec.fid, rec.bcf, rec.csi, params.gnomad_file] }
    
    // Families with existing rare/common files
    existing_rare_vcfs = bcfs_with_status
        .filter { rec -> rec.has_rare_common }
        .map { rec -> [rec.fid, rec.rare_vcf, rec.rare_tbi] }
    
    existing_common_bcfs = bcfs_with_status
        .filter { rec -> rec.has_rare_common }
        .map { rec -> [rec.fid, rec.common_bcf, rec.common_csi] }

    // Run gnomAD annotation for families that need it
    GNOMAD_FREQ_ANNOT(bcfs_need_gnomad)

    // Filter into rare and common variants
    filter_input = GNOMAD_FREQ_ANNOT.out.annotated_bcfs
        .map { fid, bcf, csi ->
            [fid, bcf, csi, params.gnomad_filter_field, params.gnomad_filter_threshold]
        }

    GNOMAD_FREQ_FILTER(filter_input)

    // Combine newly created + existing rare/common files
    all_rare_vcfs = GNOMAD_FREQ_FILTER.out.rare_vcfs.mix(existing_rare_vcfs)
    all_common_bcfs = GNOMAD_FREQ_FILTER.out.common_bcfs.mix(existing_common_bcfs)

    // =====================================================================
    // Step 3: Filter common variants (keep only GT field)
    // =====================================================================
    
    // Families needing common_gt filtering
    common_bcfs_need_gt_filter = bcfs_with_status
        .filter { rec -> rec.has_rare_common && !rec.has_common_gt }
        .map { rec -> [rec.fid, rec.common_bcf, rec.common_csi] }
    
    // Existing common_gt files
    existing_common_gt_bcfs = bcfs_with_status
        .filter { rec -> rec.has_common_gt }
        .map { rec -> [rec.fid, rec.common_gt_bcf, rec.common_gt_csi] }
    
    // Run common filters on newly created + those needing GT filtering
    common_bcfs_for_filtering = GNOMAD_FREQ_FILTER.out.common_bcfs.mix(common_bcfs_need_gt_filter)
    COMMON_FILTERS(common_bcfs_for_filtering)
    
    all_common_gt_bcfs = COMMON_FILTERS.out.filtered_common_bcfs.mix(existing_common_gt_bcfs)

    // =====================================================================
    // Step 4: VEP annotation on rare variants
    // =====================================================================
    
    // Families needing VEP annotation
    rare_vcfs_need_vep = bcfs_with_status
        .filter { rec -> rec.has_rare_common && !rec.has_vep }
        .map { rec -> [rec.fid, rec.rare_vcf, rec.rare_tbi, params.vep_config] }
    
    // Existing VEP annotated files
    existing_vep_vcfs = bcfs_with_status
        .filter { rec -> rec.has_vep }
        .map { rec -> [rec.fid, rec.vep_vcf, rec.vep_tbi] }
    
    // Run VEP on newly created rare VCFs + those needing VEP
    newly_created_rare_for_vep = GNOMAD_FREQ_FILTER.out.rare_vcfs
        .map { fid, vcf, tbi -> [fid, vcf, tbi, params.vep_config] }
    
    all_vcfs_for_vep = newly_created_rare_for_vep.mix(rare_vcfs_need_vep)
    VEP_ANNOTATION(all_vcfs_for_vep)

    // Combine new VEP outputs with existing VEP files
    all_vep_vcfs = VEP_ANNOTATION.out.annotated_vcfs.mix(existing_vep_vcfs)
    
    // =====================================================================
    // Step 5: bcftools annotation
    // =====================================================================
    
    // Families needing bcftools annotation
    vep_vcfs_need_bcftools = bcfs_with_status
        .filter { rec -> rec.has_vep && !rec.has_annotated }
        .map { rec -> [rec.fid, rec.vep_vcf, rec.vep_tbi, params.annotation_annotation_path, params.annotation_annotation_list] }
    
    // Existing fully annotated BCFs
    existing_annotated_bcfs = bcfs_with_status
        .filter { rec -> rec.has_annotated }
        .map { rec -> [rec.fid, rec.annotated_bcf, rec.annotated_csi] }
    
    // Run bcftools annotate on newly created VEP files + those needing annotation
    newly_created_vep_for_bcftools = VEP_ANNOTATION.out.annotated_vcfs
        .map { fid, vcf, tbi -> [fid, vcf, tbi, params.annotation_annotation_path, params.annotation_annotation_list] }
    
    all_vcfs_for_bcftools = newly_created_vep_for_bcftools.mix(vep_vcfs_need_bcftools)
    BCFTOOLS_ANNOTATE(all_vcfs_for_bcftools)

    // Combine new bcftools outputs with existing annotated BCFs
    all_annotated_bcfs = BCFTOOLS_ANNOTATE.out.annotated_bcfs.mix(existing_annotated_bcfs)

    emit:
    gnomad_bcfs = GNOMAD_FREQ_ANNOT.out.annotated_bcfs
    rare_vcfs = all_rare_vcfs
    common_bcfs = all_common_bcfs
    filtered_common_bcfs = all_common_gt_bcfs
    vep_annotated_vcfs = all_vep_vcfs
    fully_annotated_bcfs = all_annotated_bcfs
}
