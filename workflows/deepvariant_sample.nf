/*
 * DeepVariant Sample Variant Calling Workflow
 * Handles variant calling from CRAM files using DeepVariant and generates VAF bedgraphs
 */

// Include modules
include { MAKE_EXAMPLES } from '../modules/deepvariant_sample/make_examples'
include { CALL_VARIANTS } from '../modules/deepvariant_sample/call_variants'
include { POSTPROCESS_VARIANTS } from '../modules/deepvariant_sample/postprocess_variants'
include { VCF2BEDGRAPH_VAF } from '../modules/deepvariant_sample/vcf2bedgraph_vaf'

workflow DEEPVARIANT_SAMPLE {
    
    take:
    cram_files      // channel: [barcode, cram_file, crai_file]
    
    main:
    
    // Make examples from CRAM files
    MAKE_EXAMPLES(
        cram_files,
        params.ref
    )
    
    // Call variants using the examples
    CALL_VARIANTS(
        MAKE_EXAMPLES.out.examples,
        MAKE_EXAMPLES.out.examples_info
    )
    
    // Postprocess variants to create final VCF and gVCF files
    POSTPROCESS_VARIANTS(
        CALL_VARIANTS.out.called_variants
            .join(MAKE_EXAMPLES.out.gvcf_records)
            .join(MAKE_EXAMPLES.out.call_variant_outputs),
        params.ref
    )
    
    // Generate VAF bedgraphs from VCF files
    VCF2BEDGRAPH_VAF(
        POSTPROCESS_VARIANTS.out.vcf
    )
    
    emit:
    vcf = POSTPROCESS_VARIANTS.out.vcf
    gvcf = POSTPROCESS_VARIANTS.out.gvcf
    vaf_bedgraph = VCF2BEDGRAPH_VAF.out.vaf_bedgraph
}
