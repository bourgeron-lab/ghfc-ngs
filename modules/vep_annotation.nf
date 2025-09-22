#!/usr/bin/env nextflow

/*
 * Script to run VEP on split VCF files
 */

nextflow.enable.dsl=2

process runVEP {
  /*
  Run VEP on VCF files

  Returns
  -------
  Tuple of original VCF, split VCF file after running VEP, tabix index of that file, vep config file
  */

  publishDir "${params.data}/families/${fid}/vcfs",
    pattern: "${fid}.${params.vep_config_name}.*.vcf.gz*",
    mode:'copy'

  label 'vep'

  input:
  tuple val(meta), val(fid), path(input), path(index), path(vep_config)

  output:
  tuple val(meta), val(fid), path("${out}{.gz,}"), path("${out}{.gz,}.{tbi,csi}"), val("${vep_config}"), emit: files

  script:
  index_type = meta.index_type
  out = "${fid}.${params.vep_config_name}.${input.baseName}.vcf"

  """
  # Run VEP annotation
  vep \\
      --input_file ${input} \\
      --output_file ${out} \\
      --config ${vep_config} \\
      --compress_output bgzip \\
      --fork ${task.cpus} \\
      --format vcf \\
      --vcf \\
      --force_overwrite

  # Index the annotated VCF
  tabix -f -p vcf ${out}.gz
  """

  stub:
  """
  touch ${out}.gz
  touch ${out}.gz.tbi
  """
}