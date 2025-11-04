/*
 * Familial Pedigree Workflow
 * Extracts family-specific pedigree files from cohort pedigree file
 */

// Include modules
include { familial_pedigree } from '../modules/familial_pedigree'

workflow FAMILIAL_PEDIGREE {

    take:
    family_data    // channel: [fid, pedigree_file, project]

    main:

    // Run pedigree extraction for each family
    familial_pedigree(family_data)

    emit:
    family_pedigrees = familial_pedigree.out.family_pedigree
}
