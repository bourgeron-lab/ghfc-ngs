#!/usr/bin/env nextflow

/*
 * GLnexus family calling module for DeepVariant family workflow
 */

nextflow.enable.dsl=2

process GLNEXUS_FAMILY {

    label 'slow'
    tag "$fid"
    
    // No publishDir - this is an intermediate file kept only in work directory
    
    input:
    tuple val(fid), val(barcodes), path(gvcfs), path(tbis)
    
    output:
    tuple val(fid), path("${fid}.vcf.gz"), path("${fid}.vcf.gz.tbi"), emit: family_vcf
    
    script:
    def gvcf_list = gvcfs.join(' ')
    // Read in place rather than staged, like params.ref in NORMALIZE. A whole line or nothing,
    // so that without it the rendered script - and so the task hash - is exactly what it was.
    // Only a string counts: a bare --glnexus_bed on the command line arrives as Boolean true.
    def bed = params.glnexus_bed instanceof CharSequence ? params.glnexus_bed.toString().trim() : ''
    def bed_line = bed ? "        --bed ${bed} \\\n" : ''

    """
    # Create temporary directory for GLnexus
    rm -rf tmp_glnexus_${fid}
    
    # Run GLnexus
    glnexus_cli \\
        --config ${params.glnexus_config ?: 'DeepVariant_unfiltered'} \\
${bed_line}        --threads ${task.cpus} \\
        --mem-gbytes ${task.memory.toGiga()} \\
        --dir tmp_glnexus_${fid} \\
        ${gvcf_list} \\
        | bcftools view -Oz -o ${fid}.vcf.gz
    
    # Index the output VCF
    tabix -p vcf ${fid}.vcf.gz
    
    # Clean up temporary directory
    rm -rf tmp_glnexus_${fid}
    """
    
    stub:
    """
    touch ${fid}.vcf.gz
    touch ${fid}.vcf.gz.tbi
    """
}
