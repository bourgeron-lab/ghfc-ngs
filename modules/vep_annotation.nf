#!/usr/bin/env nextflow

/*
 * Simple VEP annotation module for family VCF files
 */

nextflow.enable.dsl=2

process VEP_ANNOTATION {
  /*
  Run VEP annotation on family VCF files

  Parameters
  ----------
  fid : str
    Family ID
  vcf : path
    Input VCF file
  vcf_index : path
    VCF index file (.tbi or .csi)
  vep_config : val
    Path to VEP configuration file

  Returns
  -------
  Tuple of family ID, annotated VCF file, and its index
  */

  tag "$fid"

  publishDir "${params.output_dir}/families/${fid}/vcfs",
    mode: 'copy',
    pattern: "${fid}.${params.vep_config_name}.vcf.gz*"

  container 'ensemblorg/ensembl-vep:release_110.1'
  
  label 'vep'

  input:
  tuple val(fid), path(vcf), path(vcf_index), val(vep_config)

  output:
  tuple val(fid), path("${output_vcf}.gz"), path("${output_vcf}.gz.tbi"), emit: annotated_vcfs

  script:
  output_vcf = "${fid}.${params.vep_config_name}.vcf"

  """
  # Run VEP annotation
  vep \\
      --input_file ${vcf} \\
      --output_file ${output_vcf} \\
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
  output_vcf = "${fid}.${params.vep_config_name}.vcf"
  """
  touch ${output_vcf}.gz
  touch ${output_vcf}.gz.tbi
  """
}