/*
 * DeepVariant Sample Variant Calling Workflow
 * Handles variant calling from CRAM files using DeepVariant and generates VAF bedgraphs
 * Skips processing if output files already exist
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
    
    // Check which samples need full DeepVariant processing vs just VAF bedgraph
    cram_with_status = cram_files
        .map { barcode, cram, crai ->
            def gvcf_path = file("${params.data}/samples/${barcode}/deepvariant/${barcode}.g.vcf.gz")
            def gvcf_tbi_path = file("${params.data}/samples/${barcode}/deepvariant/${barcode}.g.vcf.gz.tbi")
            def vcf_path = file("${params.data}/samples/${barcode}/deepvariant/${barcode}.vcf.gz")
            def vcf_tbi_path = file("${params.data}/samples/${barcode}/deepvariant/${barcode}.vcf.gz.tbi")
            def vaf_path = file("${params.data}/samples/${barcode}/sequences/${barcode}.vaf.bedgraph.gz")
            def vaf_tbi_path = file("${params.data}/samples/${barcode}/sequences/${barcode}.vaf.bedgraph.gz.tbi")
            
            def has_gvcf = gvcf_path.exists() && gvcf_tbi_path.exists()
            def has_vcf = vcf_path.exists() && vcf_tbi_path.exists()
            def has_vaf = vaf_path.exists() && vaf_tbi_path.exists()
            
            [barcode: barcode, cram: cram, crai: crai, 
             gvcf_path: gvcf_path, gvcf_tbi_path: gvcf_tbi_path,
             vcf_path: vcf_path, vcf_tbi_path: vcf_tbi_path,
             has_gvcf: has_gvcf, has_vcf: has_vcf, has_vaf: has_vaf]
        }
    
    // Samples needing full DeepVariant processing (no gVCF exists)
    crams_need_deepvariant = cram_with_status
        .filter { rec -> !rec.has_gvcf }
        .map { rec -> [rec.barcode, rec.cram, rec.crai] }
    
    // Samples with existing gVCF/VCF but missing VAF bedgraph
    samples_need_vaf_only = cram_with_status
        .filter { rec -> rec.has_gvcf && rec.has_vcf && !rec.has_vaf }
        .map { rec -> [rec.barcode, rec.vcf_path, file("${rec.vcf_path}.tbi")] }
    
    // Existing gVCFs that can be passed through
    existing_gvcfs = cram_with_status
        .filter { rec -> rec.has_gvcf }
        .map { rec -> [rec.barcode, rec.gvcf_path, rec.gvcf_tbi_path] }
    
    // Existing VCFs that can be passed through
    existing_vcfs = cram_with_status
        .filter { rec -> rec.has_vcf }
        .map { rec -> [rec.barcode, rec.vcf_path, file("${rec.vcf_path}.tbi")] }
    
    // Make examples from CRAM files (only for samples needing full processing)
    MAKE_EXAMPLES(
        crams_need_deepvariant,
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
    
    // Generate VAF bedgraphs from VCF files (newly created + existing without VAF)
    vcfs_for_vaf = POSTPROCESS_VARIANTS.out.vcf.mix(samples_need_vaf_only)
    VCF2BEDGRAPH_VAF(vcfs_for_vaf)
    
    // Combine newly created outputs with existing ones
    all_vcfs = POSTPROCESS_VARIANTS.out.vcf.mix(existing_vcfs)
    all_gvcfs = POSTPROCESS_VARIANTS.out.gvcf.mix(existing_gvcfs)
    
    emit:
    vcf = all_vcfs
    gvcf = all_gvcfs
    vaf_bedgraph = VCF2BEDGRAPH_VAF.out.vaf_bedgraph
}
