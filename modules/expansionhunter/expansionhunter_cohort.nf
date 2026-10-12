#!/usr/bin/env nextflow

/*
 * ExpansionHunter: cohort merge module
 * Merges the family ExpansionHunter VCFs and tables into cohort-level files
 */

nextflow.enable.dsl=2

process EXPANSIONHUNTER_COHORT {
  /*
  Merge every family's ExpansionHunter VCF and genotype table across the cohort

  The VCF is a plain bcftools merge of the family VCFs, one record per locus. The
  tables come from eh_cohort_merge: the concatenated genotypes with each member's
  percentile among the cohort's founders, a per-locus summary, and the flagged
  genotypes. See the script's docstring.

  Parameters
  ----------
  cohort_name : val
    Name of the cohort
  tables : path
    Every family genotype table
  vcfs : path
    Every family VCF
  tbis : path
    Their indexes
  catalog_name : val
    Label identifying the catalog in output file names
  cohort_merge : path
    The merge script, staged in so it is reachable inside the container

  Returns
  -------
  Tuple of cohort name, merged VCF and its index, and the three tables
  */

  tag "${cohort_name}"

  publishDir "${params.data}/cohorts/${cohort_name}/str/expansionhunter",
    mode: 'copy',
    pattern: "${cohort_name}.${catalog_name}{.vcf.gz,.vcf.gz.tbi,.genotypes.tsv,.loci.tsv,.flags.tsv}"

  label 'expansionhunter_cohort'

  input:
  // Inputs live in their own directories: a cohort named after its only family would
  // otherwise stage a file under the very name of an output
  tuple val(cohort_name), path(tables, stageAs: 'input_tables/*'), path(vcfs, stageAs: 'input_vcfs/*'), path(tbis, stageAs: 'input_vcfs/*'), val(catalog_name), path(cohort_merge)

  output:
  tuple val(cohort_name), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf
  tuple val(cohort_name), path("${prefix}.{genotypes,loci,flags}.tsv"), emit: tables

  script:
  prefix = "${cohort_name}.${catalog_name}"

  """
  set -euo pipefail

  export TMPDIR=\$(pwd)
  # bcftools merge holds every input open at once
  ulimit -n "\$(ulimit -Hn)" 2>/dev/null || true

  # find, not ls: a large cohort expands the glob past the kernel argument limit
  find input_vcfs -maxdepth 1 -name '*.vcf.gz' | sort > vcf_list.txt
  if [ ! -s vcf_list.txt ]; then
    echo "ERROR: no family ExpansionHunter VCFs were staged for cohort ${cohort_name}" >&2
    exit 1
  fi

  if [ "\$(wc -l < vcf_list.txt)" -eq 1 ]; then
    cp "\$(cat vcf_list.txt)" ${prefix}.vcf.gz
  else
    bcftools merge --threads ${task.cpus} --merge all --file-list vcf_list.txt -Oz -o ${prefix}.vcf.gz
  fi
  bcftools index -t ${prefix}.vcf.gz

  n_out=\$(bcftools index -n ${prefix}.vcf.gz)
  echo "Cohort ${cohort_name}: \$(wc -l < vcf_list.txt | tr -d '[:space:]') families -> \${n_out} merged records"
  if [ "\${n_out}" -eq 0 ]; then
    echo "ERROR: the cohort merge produced no records" >&2
    exit 1
  fi

  python3 ${cohort_merge} \\
      --input-dir input_tables \\
      --cohort ${cohort_name} \\
      --label ${catalog_name}
  """

  stub:
  prefix = "${cohort_name}.${catalog_name}"
  """
  echo | gzip > ${prefix}.vcf.gz
  touch ${prefix}.vcf.gz.tbi
  for kind in genotypes loci flags; do
    printf 'locus\\tfamily_id\\nHTT\\tFAM1\\n' > ${prefix}.\${kind}.tsv
  done
  """
}
