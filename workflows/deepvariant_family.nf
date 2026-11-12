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
    family_gvcfs        // channel: [fid, [barcodes], [gvcfs], [tbis]]
    pedigree_file       // path: cohort pedigree file
    outdated_pedigrees  // list: families whose {FID}.pedigree.tsv does not hold their cohort rows,
                        // taken as missing so that it is written again
    pedigree_rewrites   // list: those of them not in family_gvcfs, whose pedigree alone is written

    main:

    // Check existing outputs for each family
    family_with_status = family_gvcfs
        .map { fid, barcodes, gvcfs, tbis ->
            def fam_dir = Sharding.getFamilyDir(params.data, fid)
            def norm_bcf = file("${fam_dir}/vcfs/${fid}.norm.bcf")
            def norm_csi = file("${fam_dir}/vcfs/${fid}.norm.bcf.csi")
            def pedigree = file("${fam_dir}/${fid}.pedigree.tsv")
            
            def has_norm_bcf = norm_bcf.exists() && norm_csi.exists()
            // An empty pedigree must not count as present, or it is passed through unchanged forever;
            // nor an out-of-date one, which would be
            def has_pedigree = pedigree.exists() && pedigree.size() > 0 && !(fid in outdated_pedigrees)
            
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
            def pedigree = file("${Sharding.getFamilyDir(params.data, fid)}/${fid}.pedigree.tsv")
            [fid, pedigree.exists() && pedigree.size() > 0 && !(fid in outdated_pedigrees)]
        }
        .filter { _fid, has_pedigree -> !has_pedigree }
        .map { fid, _has_pedigree ->
            [fid, pedigree_file, params.cohort_name]
        }
    
    // Called families whose pedigree file is out of date, and that are not called again here
    families_pedigree_rewrite = Channel.fromList(pedigree_rewrites)
        .map { fid -> [fid, pedigree_file, params.cohort_name] }

    all_families_needing_pedigree = families_needing_pedigree_from_glnexus
        .mix(families_need_pedigree_only)
        .mix(families_pedigree_rewrite)
    FAMILIAL_PEDIGREE(all_families_needing_pedigree)

    // Combine newly created outputs with existing ones
    all_normalized_bcfs = NORMALIZE.out.normalized_bcf.mix(families_with_norm_bcf)
    all_family_pedigrees = FAMILIAL_PEDIGREE.out.family_pedigree.mix(families_with_pedigree)

    emit:
    family_vcfs = GLNEXUS_FAMILY.out.family_vcf
    normalized_bcfs = all_normalized_bcfs
    family_pedigrees = all_family_pedigrees
}
