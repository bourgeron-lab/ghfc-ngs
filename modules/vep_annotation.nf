#!/usr/bin/env nextflow

/*
 * Simple VEP annotation module for family VCF files
 */

nextflow.enable.dsl=2

process runVEP {
  /*
  Run VEP annotation on VCF files, preserving original filename structure

  Parameters
  ----------
  fid : str
    Family ID
  vcf : path
    Input VCF file (e.g., FID.rare.vcf.gz)
  vcf_index : path
    VCF index file (.tbi or .csi)
  vep_config : val
    Path to VEP configuration file

  Returns
  -------
  Tuple of family ID, annotated VCF file with preserved name structure (e.g., FID.rare.vep_config_name.vcf.gz), and its index
  */

  tag "$fid"

  publishDir "${params.data}/families/${fid}/vcfs",
    mode: 'copy',
    pattern: "*.${params.vep_config_name}.vcf.gz*"

  container 'ensemblorg/ensembl-vep:release_110.1'
  
  label 'vep'

  input:
  tuple val(fid), path(vcf), path(vcf_index), val(vep_config)

  output:
  tuple val(fid), path("*.${params.vep_config_name}.vcf.gz"), path("*.${params.vep_config_name}.vcf.gz.tbi"), emit: annotated_vcfs

  script:
  // Extract base filename without .vcf.gz extension and add VEP config name
  base_name = vcf.baseName.replaceAll(/\.vcf$/, '')
  output_vcf = "${base_name}.${params.vep_config_name}.vcf"

  """
  # Run VEP annotation
  vep \\
      --input_file ${vcf} \\
      --output_file ${output_vcf}.gz \\
      --config ${vep_config} \\
      --compress_output bgzip \\
      --fork ${task.cpus} \\
      --format vcf \\
      --vcf \\
      --force_overwrite

  # Index the annotated VCF
  tabix -f -p vcf ${output_vcf}.gz
  """

  stub:
  // Extract base filename without .vcf.gz extension and add VEP config name
  base_name = vcf.baseName.replaceAll(/\.vcf$/, '')
  output_vcf = "${base_name}.${params.vep_config_name}.vcf"
  """
  touch ${output_vcf}.gz
  touch ${output_vcf}.gz.tbi
  """
}