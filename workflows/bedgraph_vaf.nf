/*
 * VAF Bedgraph Generation Workflow
 * Converts VCF files to VAF (Variant Allele Frequency) bedgraph format
 */

// Include modules
include { VCF2BEDGRAPH_VAF } from '../modules/vcf2bedgraph_vaf'

workflow BEDGRAPH_VAF {

    take:
    vcf_files    // channel: [barcode, vcf, tbi]

    main:

    // Run VCF to VAF bedgraph conversion
    VCF2BEDGRAPH_VAF(vcf_files)

    emit:
    vaf_bedgraphs = VCF2BEDGRAPH_VAF.out.vaf_bedgraph
}
