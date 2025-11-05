#!/usr/bin/env nextflow

/*
 * Migration module: Convert legacy normalized VCF.gz files to BCF format
 */

nextflow.enable.dsl=2

process CONVERT_NORM_VCF_TO_BCF {
    /*
    Convert legacy {FID}.norm.vcf.gz files to {FID}.norm.bcf format
    
    Parameters
    ----------
    fid : str
        Family ID
    vcf : path
        Legacy normalized VCF.gz file
    tbi : path
        VCF.gz index file (.tbi)
    
    Returns
    -------
    Tuple of family ID, new BCF file, and its index
    */
    
    tag "$fid"
    
    publishDir "${params.data}/families/${fid}/vcfs",
        mode: 'copy',
        pattern: "${fid}.norm.bcf*"
    
    input:
    tuple val(fid), path(vcf), path(tbi)
    
    output:
    tuple val(fid), path("${fid}.norm.bcf"), path("${fid}.norm.bcf.csi"), emit: converted_bcf
    
    script:
    """
    echo "Converting ${vcf} to BCF format for family ${fid}"
    
    # Convert VCF.gz to BCF
    bcftools view \\
        -O b \\
        -o ${fid}.norm.bcf \\
        ${vcf}
    
    # Index the new BCF file
    bcftools index ${fid}.norm.bcf
    
    echo "Conversion completed: ${fid}.norm.bcf created"
    """
    
    stub:
    """
    touch ${fid}.norm.bcf
    touch ${fid}.norm.bcf.csi
    """
}
