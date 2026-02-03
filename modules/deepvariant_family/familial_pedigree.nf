#!/usr/bin/env nextflow

/*
 * Familial pedigree extraction module for DeepVariant family workflow
 * Extracts family-specific pedigree information from cohort pedigree file
 */

nextflow.enable.dsl=2

process FAMILIAL_PEDIGREE {
  /*
  Extract family-specific pedigree subset from cohort pedigree file

  Parameters
  ----------
  fid : str
    Family ID
  pedigree_file : path
    Cohort pedigree TSV file
  project : val
    Project/cohort name

  Returns
  -------
  Tuple of family ID and family-specific pedigree TSV file
  */

  tag "$fid"

  publishDir "${params.data}/families/${fid}",
    mode: 'copy',
    pattern: "${fid}.pedigree.tsv"

  label 'pedigree_extraction'

  input:
  tuple val(fid), path(pedigree_file), val(project)

  output:
  tuple val(fid), path("${fid}.pedigree.tsv"), emit: family_pedigree

  script:
  """
  #!/bin/bash
  set -euo pipefail
  
  # Check if the first line is a header (starts with "FID" followed by tab or space)
  first_field=\$(head -n 1 ${pedigree_file} | awk '{print \$1}')
  
  if [ "\${first_field}" = "FID" ]; then
    # Has header - include it and extract family rows (skip header in awk)
    head -n 1 ${pedigree_file} > ${fid}.pedigree.tsv
    awk -F'\\t' 'NR > 1 && \$1 == "${fid}"' ${pedigree_file} >> ${fid}.pedigree.tsv
  else
    # No header - extract all matching rows directly
    awk -F'\\t' '\$1 == "${fid}"' ${pedigree_file} > ${fid}.pedigree.tsv
  fi
  """

  stub:
  """
  touch ${fid}.pedigree.tsv
  """
}
