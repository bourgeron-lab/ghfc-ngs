#!/usr/bin/env nextflow

/*
 * Migration module: Convert legacy common VCF.gz files to BCF format
 */

nextflow.enable.dsl=2

process CONVERT_COMMON_VCF_TO_BCF {
    /*
    Convert legacy {FID}.common.vcf.gz files to {FID}.common.bcf format
    
    Parameters
    ----------
    fid : str
        Family ID
    vcf : path
        Legacy common variants VCF.gz file
    tbi : path
        VCF.gz index file (.tbi)
    
    Returns
    -------
    Tuple of family ID, new BCF file, and its index
    */
    
    tag "$fid"
    
    publishDir "${params.data}/families/${fid}/vcfs",
        mode: 'copy',
        pattern: "${fid}.common.bcf*"
    
    input:
    tuple val(fid), path(vcf), path(tbi)
    
    output:
    tuple val(fid), path("${fid}.common.bcf"), path("${fid}.common.bcf.csi"), emit: converted_bcf
    
    script:
    """
    echo "Converting ${vcf} to BCF format for family ${fid}"
    
    # Convert VCF.gz to BCF
    bcftools view \\
        -O b \\
        -o ${fid}.common.bcf \\
        ${vcf}
    
    # Index the new BCF file
    bcftools index ${fid}.common.bcf
    
    echo "Conversion completed: ${fid}.common.bcf created"
    """
    
    stub:
    """
    touch ${fid}.common.bcf
    touch ${fid}.common.bcf.csi
    """
}
