#!/usr/bin/env nextflow

/*
 * Wombat: VEP annotation module for family BCF files
 */

nextflow.enable.dsl=2

process W_VEP_ANNOTATION {
  /*
  Run VEP annotation on BCF files, preserving original filename structure and outputting BCF

  Parameters
  ----------
  fid : str
    Family ID
  bcf : path
    Input BCF file (e.g., FID.rare.bcf)
  bcf_index : path
    BCF index file (.csi)
  vep_config : val
    Path to VEP configuration file

  Returns
  -------
  Tuple of family ID, annotated BCF file with preserved name structure (e.g., FID.rare.vep_config_name.bcf), and its index
  */

  tag "$fid"

  publishDir "${params.data}/families/${fid}/vcfs",
    mode: 'copy',
    pattern: "*.${params.vep_config_name}.bcf*"

  container 'ensemblorg/ensembl-vep:release_110.1'
  
  label 'vep'

  input:
  tuple val(fid), path(bcf), path(bcf_index), val(vep_config)

  output:
  tuple val(fid), path("${output_bcf}"), path("${output_bcf}.csi"), emit: annotated_bcfs

  script:
  // Extract base filename without .bcf extension and add VEP config name
  base_name = bcf.baseName.replaceAll(/\.bcf$/, '')
  output_vcf = "${base_name}.${params.vep_config_name}.vcf"
  output_bcf = "${base_name}.${params.vep_config_name}.bcf"

  """
  # Run VEP annotation (VEP outputs VCF)
  vep \\
      --input_file ${bcf} \\
      --output_file ${output_vcf}.gz \\
      --config ${vep_config} \\
      --compress_output bgzip \\
      --fork ${task.cpus} \\
      --format vcf \\
      --vcf \\
      --force_overwrite

  # Convert VEP output to BCF
  bcftools view \\
      -O b \\
      -o ${output_bcf} \\
      ${output_vcf}.gz

  # Index the BCF
  bcftools index ${output_bcf}
  """

  stub:
  // Extract base filename without .bcf extension and add VEP config name
  base_name = bcf.baseName.replaceAll(/\.bcf$/, '')
  output_bcf = "${base_name}.${params.vep_config_name}.bcf"
  """
  touch ${output_bcf}
  touch ${output_bcf}.csi
  """
}
