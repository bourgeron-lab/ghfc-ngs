#!/usr/bin/env nextflow

/*
 * VCF normalization module for DeepVariant family workflow
 * Outputs BCF format instead of VCF
 */

nextflow.enable.dsl=2

process GL_NORMALIZE {

    tag "$fid"

    publishDir "${params.data}/families/${fid}/vcfs", mode: 'copy'

    input:
    tuple val(fid), path(vcf), path(tbi)

    output:
    tuple val(fid), path("${fid}.norm.bcf"), path("${fid}.norm.bcf.csi"), emit: normalized_bcf

    script:
    """
    # Normalize the VCF file and output as BCF
    zcat ${vcf} | \\
        sed -e 's/ID=AD,Number=\\./ID=AD,Number=R/' | \\
        bcftools norm -m - -w 1000 -c s -f ${params.ref} --threads ${task.cpus} -Ob -o ${fid}.norm.bcf

    # Index the normalized BCF
    bcftools index ${fid}.norm.bcf
    """

    stub:
    """
    touch ${fid}.norm.bcf
    touch ${fid}.norm.bcf.csi
    """
}
