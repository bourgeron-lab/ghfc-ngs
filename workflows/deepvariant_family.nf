#!/usr/bin/env nextflow

/*
 * DeepVariant Family Workflow
 * Integrates GLnexus family calling, VCF normalization, and pedigree extraction
 * Skips processing if output files already exist - only runs necessary modules
 */

// Include modules
include { GLNEXUS_FAMILY } from '../modules/deepvariant_family/glnexus_family'
include { NORMALIZE } from '../modules/deepvariant_family/normalize_vcf'
include { FAMILIAL_PEDIGREE } from '../modules/deepvariant_family/familial_pedigree'

workflow DEEPVARIANT_FAMILY {

    take:
    family_gvcfs    // channel: [fid, [barcodes], [gvcfs], [tbis]]
    pedigree_file   // path: cohort pedigree file

    main:

    // Check existing outputs for each family
    family_with_status = family_gvcfs
        .map { fid, barcodes, gvcfs, tbis ->
            def norm_bcf = file("${params.data}/families/${fid}/vcfs/${fid}.norm.bcf")
            def norm_csi = file("${params.data}/families/${fid}/vcfs/${fid}.norm.bcf.csi")
            def pedigree = file("${params.data}/families/${fid}/${fid}.pedigree.tsv")
            
            def has_norm_bcf = norm_bcf.exists() && norm_csi.exists()
            def has_pedigree = pedigree.exists()
            
            [fid: fid, barcodes: barcodes, gvcfs: gvcfs, tbis: tbis,
             norm_bcf: norm_bcf, norm_csi: norm_csi, pedigree: pedigree,
             has_norm_bcf: has_norm_bcf, has_pedigree: has_pedigree]
        }
    
    // Families needing GLnexus + Normalize (no normalized BCF exists)
    families_need_glnexus = family_with_status
        .filter { rec -> !rec.has_norm_bcf }
        .map { rec -> [rec.fid, rec.barcodes, rec.gvcfs, rec.tbis] }
    
    // Families with existing normalized BCF (for passthrough)
    families_with_norm_bcf = family_with_status
        .filter { rec -> rec.has_norm_bcf }
        .map { rec -> [rec.fid, rec.norm_bcf, rec.norm_csi] }
    
    // Families needing pedigree extraction only
    families_need_pedigree_only = family_with_status
        .filter { rec -> rec.has_norm_bcf && !rec.has_pedigree }
        .map { rec -> [rec.fid, pedigree_file, params.cohort_name] }
    
    // Families with existing pedigree (for passthrough)
    families_with_pedigree = family_with_status
        .filter { rec -> rec.has_pedigree }
        .map { rec -> [rec.fid, rec.pedigree] }

    // Run GLnexus for families that need it
    GLNEXUS_FAMILY(families_need_glnexus)

    // Normalize family VCF files (only for newly created ones)
    NORMALIZE(GLNEXUS_FAMILY.out.family_vcf)

    // Extract family-specific pedigree files
    // Include families that ran through GLnexus (need pedigree) + families that only need pedigree
    families_needing_pedigree_from_glnexus = GLNEXUS_FAMILY.out.family_vcf
        .map { fid, _vcf, _tbi ->
            // Check if pedigree already exists for this family
            def pedigree = file("${params.data}/families/${fid}/${fid}.pedigree.tsv")
            [fid, pedigree.exists()]
        }
        .filter { _fid, has_pedigree -> !has_pedigree }
        .map { fid, _has_pedigree ->
            [fid, pedigree_file, params.cohort_name]
        }
    
    all_families_needing_pedigree = families_needing_pedigree_from_glnexus.mix(families_need_pedigree_only)
    FAMILIAL_PEDIGREE(all_families_needing_pedigree)

    // Combine newly created outputs with existing ones
    all_normalized_bcfs = NORMALIZE.out.normalized_bcf.mix(families_with_norm_bcf)
    all_family_pedigrees = FAMILIAL_PEDIGREE.out.family_pedigree.mix(families_with_pedigree)

    emit:
    family_vcfs = GLNEXUS_FAMILY.out.family_vcf
    normalized_bcfs = all_normalized_bcfs
    family_pedigrees = all_family_pedigrees
}
