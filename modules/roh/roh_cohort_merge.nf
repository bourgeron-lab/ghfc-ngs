#!/usr/bin/env nextflow

/*
 * ROH: cohort table merge module
 * Concatenates the family ROH tables into cohort tables and counts cohort ROH islands
 */

nextflow.enable.dsl=2

process ROH_COHORT_MERGE {
  /*
  Concatenate the family ROH tables across the cohort, and count ROH islands

  One task for every table kind, unlike the ancestry merge, because two kinds are
  legitimately absent for some families (a singleton has no kinship table) and the
  window count needs the per-sample and segment tables together. See the script's
  docstring for the window table.

  Parameters
  ----------
  cohort_name : val
    Name of the cohort
  family_tables : list
    Every family table to merge, of any kind
  panel_name : val
    Label identifying this panel/threshold combination in output file names
  roh_cohort_merge : path
    The merge script, staged in so it is reachable inside the container

  Returns
  -------
  Tuple of cohort name and the merged tables
  */

  tag "${cohort_name}"

  publishDir "${params.data}/cohorts/${cohort_name}/roh",
    mode: 'copy',
    pattern: "${cohort_name}.${panel_name}.{froh,roh,roh_genes,kinship,roh_windows}.tsv"

  label 'roh_cohort_merge'

  input:
  tuple val(cohort_name), path(family_tables, stageAs: 'input_tables/*'), val(panel_name), path(roh_cohort_merge)

  output:
  tuple val(cohort_name), path("${cohort_name}.${panel_name}.*.tsv"), emit: cohort_tables

  script:
  """
  set -euo pipefail

  python3 ${roh_cohort_merge} \\
      --input-dir input_tables \\
      --cohort ${cohort_name} \\
      --panel ${panel_name}
  """

  stub:
  """
  # Like the script: a kind is written only when some family staged a table of it
  kinds="froh roh roh_windows"
  for kind in roh_genes kinship; do
    if ls input_tables/*.\${kind}.tsv >/dev/null 2>&1; then kinds="\${kinds} \${kind}"; fi
  done
  for kind in \${kinds}; do
    printf 'IID\\tfamily_id\\nS1\\tFAM1\\n' > ${cohort_name}.${panel_name}.\${kind}.tsv
  done
  """
}
