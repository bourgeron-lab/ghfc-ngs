#!/usr/bin/env nextflow

/*
 * Script to merge VCF files into a single file
 */

nextflow.enable.dsl=2

process mergeVCF {
  /*
  Merge VCF files into a single file
  */

  cpus params.cpus
  label 'bcftools'
  cache 'lenient'

  input:
  tuple val(meta), val(original_file), path(vcf_files), path(index_files), val(vep_config)

  output:
  tuple val(meta), path("${original_file.baseName}_VEP.vcf.gz"), path("${original_file.baseName}_VEP.vcf.gz.{tbi,csi}")

  script:
  index_type = meta.index_type
  index_flag = index_type == "tbi" ? "-t" : "-c"
  merged_vcf = "${original_file.baseName}_VEP.vcf.gz"

  """
  sorted_vcfs=\$(echo ${vcf_files} | xargs -n1 | sort | xargs)
  bcftools concat --no-version --naive \${sorted_vcfs} -Oz -o ${merged_vcf}
  bcftools index ${index_flag} ${merged_vcf}
  """
}