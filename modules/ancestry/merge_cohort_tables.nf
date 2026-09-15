#!/usr/bin/env nextflow

/*
 * Ancestry: cohort table merge module
 * Concatenates one kind of family ancestry/PGS table into a cohort-level table
 */

nextflow.enable.dsl=2

process ANCESTRY_COHORT_MERGE {
  /*
  Concatenate a family-level ancestry or PGS table across the cohort

  Concatenation is exact here, not an approximation. Every value in these tables is
  a per-sample projection against the reference bundle: no cohort-derived quantity
  enters, so a sample's row is identical whether its family was scored alone or as
  part of a cohort run. A family_id column is appended so rows stay traceable.

  Parameters
  ----------
  cohort_name : val
    Name of the cohort
  table_kind : val
    Which table is being merged (pcs, ancestry, Q, pgs_raw, pgs_adjusted, pgs_zscore)
  family_tables : list
    The family-level tables to concatenate
  panel_name : val
    Label identifying this panel/threshold combination in output file names

  Returns
  -------
  Tuple of cohort name, table kind, and the merged cohort table
  */

  tag "${cohort_name}-${table_kind}"

  publishDir "${params.data}/cohorts/${cohort_name}/ancestry",
    mode: 'copy',
    pattern: "${cohort_name}.${panel_name}.${table_kind}.tsv"

  label 'ancestry_cohort_merge'

  input:
  tuple val(cohort_name), val(table_kind), path(family_tables, stageAs: 'input_tables/*'), val(panel_name)

  output:
  tuple val(cohort_name), val(table_kind), path("${output_tsv}"), emit: cohort_table

  script:
  output_tsv = "${cohort_name}.${panel_name}.${table_kind}.tsv"

  """
  set -euo pipefail

  # Tolerate a non-matching glob so the explicit check below reports it clearly
  ls -1 input_tables/*.${table_kind}.tsv 2>/dev/null | sort > file_list.txt || true

  if [ ! -s file_list.txt ]; then
    echo "ERROR: no family ${table_kind} tables were staged for cohort ${cohort_name}" >&2
    exit 1
  fi

  # Take the header from the first file and append a family_id column. Every file
  # comes from the same catalog and the same bundle, so the columns must agree;
  # a mismatch means tables from different runs have been mixed and the merged
  # table would silently misalign values under the wrong trait names.
  first_file=\$(head -n 1 file_list.txt)
  header=\$(head -n 1 "\${first_file}")
  printf '%s\\tfamily_id\\n' "\${header}" > ${output_tsv}

  while IFS= read -r file; do
    this_header=\$(head -n 1 "\${file}")
    if [ "\${this_header}" != "\${header}" ]; then
      echo "ERROR: \${file} has different columns from \${first_file}" >&2
      exit 1
    fi
    # Family ID is everything before the panel label in the file name. Literal
    # suffix removal, not sed: a panel label containing regex metacharacters
    # would otherwise strip the wrong thing.
    fid=\$(basename "\${file}")
    fid="\${fid%.${panel_name}.${table_kind}.tsv}"
    tail -n +2 "\${file}" | awk -v fid="\${fid}" 'NF{print \$0 "\\t" fid}' >> ${output_tsv}
  done < file_list.txt

  n_families=\$(wc -l < file_list.txt | tr -d '[:space:]')
  n_rows=\$(( \$(wc -l < ${output_tsv} | tr -d '[:space:]') - 1 ))
  echo "Cohort ${cohort_name} ${table_kind}: \${n_rows} rows from \${n_families} families"
  if [ "\${n_rows}" -eq 0 ]; then
    echo "ERROR: merged ${table_kind} table has no data rows" >&2
    exit 1
  fi
  """

  stub:
  output_tsv = "${cohort_name}.${panel_name}.${table_kind}.tsv"
  """
  printf 'IID\\tfamily_id\\nS1\\tFAM1\\n' > ${output_tsv}
  """
}
