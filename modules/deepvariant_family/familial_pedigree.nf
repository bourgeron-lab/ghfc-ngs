#!/usr/bin/env nextflow

/*
 * Familial pedigree extraction module for DeepVariant family workflow
 * Extracts family-specific pedigree information from cohort pedigree file
 */

nextflow.enable.dsl=2

process FAMILIAL_PEDIGREE {

  label 'fast'
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

  publishDir {
    def (s1, s2) = Sharding.getShards(fid)
    "${params.data}/families/${s1}/${s2}/${fid}"
  },
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

  # Collect this family's rows first so they can be counted before anything is written
  if [ "\${first_field}" = "FID" ]; then
    awk -F'\\t' 'NR > 1 && \$1 == "${fid}"' ${pedigree_file} > family_rows.tmp
  else
    awk -F'\\t' '\$1 == "${fid}"' ${pedigree_file} > family_rows.tmp
  fi

  n_rows=\$(wc -l < family_rows.tmp | tr -d '[:space:]')

  # Never publish an empty pedigree: a 0-byte file would be indistinguishable from a completed one
  if [ "\${n_rows}" -eq 0 ]; then
    echo "ERROR: no rows matched family ID '${fid}' in ${pedigree_file}" >&2
    echo "The pedigree must be TAB-separated with the family ID in column 1." >&2
    echo "Family IDs found in column 1 (first 20):" >&2
    awk -F'\\t' '\$1 != "FID" {print "  [" \$1 "]"}' ${pedigree_file} | sort -u | head -20 >&2
    exit 1
  fi

  if [ "\${first_field}" = "FID" ]; then
    head -n 1 ${pedigree_file} > ${fid}.pedigree.tsv
  else
    : > ${fid}.pedigree.tsv
  fi
  cat family_rows.tmp >> ${fid}.pedigree.tsv

  echo "Extracted \${n_rows} individuals for family ${fid}"
  """

  stub:
  """
  # Must be non-empty: publishDir copies this into the data tree, where an empty
  # pedigree would later be mistaken for a completed one
  printf '%s\\tSTUB_SAMPLE\\t0\\t0\\t0\\t-9\\n' "${fid}" > ${fid}.pedigree.tsv
  """
}
