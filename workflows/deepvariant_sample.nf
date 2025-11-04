/*
 * DeepVariant Sample Variant Calling Workflow
 * Handles variant calling from CRAM files using DeepVariant and generates VAF bedgraphs
 */

// Include modules
include { DV_MAKE_EXAMPLES } from '../modules/dv_make_examples'
include { DV_CALL_VARIANTS } from '../modules/dv_call_variants'
include { DV_POSTPROCESS_VARIANTS } from '../modules/dv_postprocess_variants'
include { DV_VCF2BEDGRAPH_VAF } from '../modules/dv_vcf2bedgraph_vaf'

workflow DEEPVARIANT_SAMPLE {
    
    take:
    cram_files      // channel: [barcode, cram_file, crai_file]
    
    main:
    
    // Make examples from CRAM files
    DV_MAKE_EXAMPLES(
        cram_files,
        params.ref
    )
    
    // Call variants using the examples
    DV_CALL_VARIANTS(
        DV_MAKE_EXAMPLES.out.examples,
        DV_MAKE_EXAMPLES.out.examples_info
    )
    
    // Postprocess variants to create final VCF and gVCF files
    DV_POSTPROCESS_VARIANTS(
        DV_CALL_VARIANTS.out.called_variants
            .join(DV_MAKE_EXAMPLES.out.gvcf_records)
            .join(DV_MAKE_EXAMPLES.out.call_variant_outputs),
        params.ref
    )
    
    // Generate VAF bedgraphs from VCF files
    DV_VCF2BEDGRAPH_VAF(
        DV_POSTPROCESS_VARIANTS.out.vcf
    )
    
    emit:
    vcf = DV_POSTPROCESS_VARIANTS.out.vcf
    gvcf = DV_POSTPROCESS_VARIANTS.out.gvcf
    vaf_bedgraph = DV_VCF2BEDGRAPH_VAF.out.vaf_bedgraph
}
