#!/usr/bin/env nextflow

/*
 * SNVs frequency annotation module using gnomAD
 */

nextflow.enable.dsl=2

process snvs_gnomad_freq_annot {
  /*
  Annotate normalized VCF files with gnomAD frequency data

  Parameters
  ----------
  fid : str
    Family ID
  vcf : path
    Normalized VCF file
  vcf_index : path
    VCF index file (.tbi or .csi)
  gnomad_file : val
    Path to gnomAD annotation file

  Returns
  -------
  Tuple of family ID, gnomAD annotated VCF file, and its index
  */

  tag "$fid"

  publishDir "${params.data}/families/${fid}/vcfs",
    mode: 'copy',
    pattern: "${fid}.gnomad.vcf.gz*"

  container 'docker://staphb/bcftools:latest'
  
  label 'gnomad_annotation'

  input:
  tuple val(fid), path(vcf), path(vcf_index), val(gnomad_file)

  output:
  tuple val(fid), path("${output_vcf}.gz"), path("${output_vcf}.gz.tbi"), emit: annotated_vcfs

  script:
  output_vcf = "${fid}.gnomad.vcf"

  """
  # Annotate VCF with gnomAD frequencies
  bcftools annotate \\
      -a ${gnomad_file} \\
      --threads ${task.cpus} \\
      -c INFO \\
      -O z \\
      -o ${output_vcf}.gz \\
      ${vcf}

  # Index the annotated VCF
  bcftools index --tbi ${output_vcf}.gz
  """

  stub:
  output_vcf = "${fid}.gnomad.vcf"
  """
  touch ${output_vcf}.gz
  touch ${output_vcf}.gz.tbi
  """
}