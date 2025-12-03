#!/usr/bin/env nextflow

/*
 * Annotation: VEP annotation module for family VCF files
 */

nextflow.enable.dsl=2

process VEP_ANNOTATION {
  /*
  Run VEP annotation on VCF.gz files, preserving original filename structure and outputting VCF.gz

  Parameters
  ----------
  fid : str
    Family ID
  vcf : path
    Input VCF.gz file (e.g., FID.rare.vcf.gz)
  vcf_index : path
    VCF index file (.tbi)
  vep_config : val
    Path to VEP configuration file

  Returns
  -------
  Tuple of family ID, annotated VCF.gz file with preserved name structure (e.g., FID.rare.vep_config_name.vcf.gz), and its index
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
  tuple val(fid), path("${output_vcf}"), path("${output_vcf}.tbi"), emit: annotated_vcfs

  script:
  // Extract base filename without .vcf.gz extension and add VEP config name
  base_name = vcf.baseName.replaceAll(/\.vcf$/, '')
  output_vcf = "${base_name}.${params.vep_config_name}.vcf.gz"

  """
  # Run VEP annotation (VEP outputs VCF)
  vep \\
      --input_file ${vcf} \\
      --output_file ${output_vcf} \\
      --config ${vep_config} \\
      --compress_output bgzip \\
      --fork ${(1.3 * task.cpus).round()} \\
      --format vcf \\
      --vcf \\
      --force_overwrite

  # Index the VCF.gz with tabix
  tabix -p vcf ${output_vcf}
  """

  stub:
  // Extract base filename without .vcf.gz extension and add VEP config name
  base_name = vcf.baseName.replaceAll(/\.vcf$/, '')
  output_vcf = "${base_name}.${params.vep_config_name}.vcf.gz"
  """
  touch ${output_vcf}
  touch ${output_vcf}.tbi
  """
}
