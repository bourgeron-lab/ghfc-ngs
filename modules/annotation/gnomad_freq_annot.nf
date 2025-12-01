#!/usr/bin/env nextflow

/*
 * Annotation: SNVs frequency annotation module using gnomAD
 */

nextflow.enable.dsl=2

process GNOMAD_FREQ_ANNOT {
  /*
  Annotate normalized BCF files with gnomAD frequency data, output as BCF

  Parameters
  ----------
  fid : str
    Family ID
  bcf : path
    Normalized BCF file
  bcf_index : path
    BCF index file (.csi)
  gnomad_file : val
    Path to gnomAD annotation file

  Returns
  -------
  Tuple of family ID, gnomAD annotated BCF file, and its index
  */

  tag "$fid"

  // No publishDir - this is an intermediate file kept only in work directory

  container 'docker://staphb/bcftools:latest'
  
  label 'gnomad_annotation'

  input:
  tuple val(fid), path(bcf), path(bcf_index), val(gnomad_file)

  output:
  tuple val(fid), path("${output_bcf}"), path("${output_bcf}.csi"), emit: annotated_bcfs

  script:
  output_bcf = "${fid}.gnomad.bcf"

  """
  # Annotate BCF with gnomAD frequencies
  bcftools annotate \\
      -a ${gnomad_file} \\
      --threads ${task.cpus} \\
      -c INFO \\
      -O b \\
      -o ${output_bcf} \\
      ${bcf}

  # Index the annotated BCF
  bcftools index ${output_bcf}
  """

  stub:
  output_bcf = "${fid}.gnomad.bcf"
  """
  touch ${output_bcf}
  touch ${output_bcf}.csi
  """
}
