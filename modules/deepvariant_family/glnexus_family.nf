#!/usr/bin/env nextflow

/*
 * GLnexus family calling module for DeepVariant family workflow
 */

nextflow.enable.dsl=2

process GLNEXUS_FAMILY {
    
    tag "$fid"
    
    publishDir "${params.data}/families/${fid}/vcfs", mode: 'copy'
    
    input:
    tuple val(fid), val(barcodes), path(gvcfs), path(tbis)
    
    output:
    tuple val(fid), path("${fid}.vcf.gz"), path("${fid}.vcf.gz.tbi"), emit: family_vcf
    
    script:
    def gvcf_list = gvcfs.join(' ')
    
    """
    # Create temporary directory for GLnexus
    rm -rf tmp_glnexus_${fid}
    
    # Run GLnexus
    glnexus_cli \\
        --config ${params.glnexus_config ?: 'DeepVariant_unfiltered'} \\
        --threads ${task.cpus} \\
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
