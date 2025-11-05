#!/usr/bin/env nextflow

/*
 * Migration module: Clean up legacy intermediate files
 */

nextflow.enable.dsl=2

process CLEANUP_LEGACY_FILES {
    /*
    Remove obsolete legacy files after successful migration
    
    Parameters
    ----------
    fid : str
        Family ID
    data_dir : val
        Data directory path
    files_to_remove : val
        List of legacy file patterns to remove
    
    Returns
    -------
    Tuple of family ID and cleanup status
    */
    
    tag "$fid"
    
    container 'docker://staphb/bcftools:latest'
    
    input:
    tuple val(fid), val(data_dir), val(files_to_remove)
    
    output:
    tuple val(fid), val("cleaned"), emit: status
    
    script:
    def vcfs_dir = "${data_dir}/families/${fid}/vcfs"
    """
    echo "Cleaning up legacy files for family ${fid}"
    
    # Remove legacy normalized VCF.gz if it exists
    if [ -f "${vcfs_dir}/${fid}.norm.vcf.gz" ]; then
        echo "  Removing ${fid}.norm.vcf.gz and its index"
        rm -f "${vcfs_dir}/${fid}.norm.vcf.gz"
        rm -f "${vcfs_dir}/${fid}.norm.vcf.gz.tbi"
    fi
    
    # Remove legacy common VCF.gz if it exists
    if [ -f "${vcfs_dir}/${fid}.common.vcf.gz" ]; then
        echo "  Removing ${fid}.common.vcf.gz and its index"
        rm -f "${vcfs_dir}/${fid}.common.vcf.gz"
        rm -f "${vcfs_dir}/${fid}.common.vcf.gz.tbi"
    fi
    
    # Remove obsolete intermediate files
    if [ -f "${vcfs_dir}/${fid}.vcf.gz" ]; then
        echo "  Removing obsolete ${fid}.vcf.gz (GLnexus intermediate)"
        rm -f "${vcfs_dir}/${fid}.vcf.gz"
        rm -f "${vcfs_dir}/${fid}.vcf.gz.tbi"
    fi
    
    if [ -f "${vcfs_dir}/${fid}.gnomad.vcf.gz" ]; then
        echo "  Removing obsolete ${fid}.gnomad.vcf.gz (gnomAD intermediate)"
        rm -f "${vcfs_dir}/${fid}.gnomad.vcf.gz"
        rm -f "${vcfs_dir}/${fid}.gnomad.vcf.gz.tbi"
    fi
    
    # Also check for BCF versions of obsolete intermediates
    if [ -f "${vcfs_dir}/${fid}.gnomad.bcf" ]; then
        echo "  Removing obsolete ${fid}.gnomad.bcf (gnomAD intermediate)"
        rm -f "${vcfs_dir}/${fid}.gnomad.bcf"
        rm -f "${vcfs_dir}/${fid}.gnomad.bcf.csi"
    fi
    
    echo "Cleanup completed for family ${fid}"
    """
    
    stub:
    """
    echo "Stub: Would clean up legacy files for ${fid}"
    """
}
