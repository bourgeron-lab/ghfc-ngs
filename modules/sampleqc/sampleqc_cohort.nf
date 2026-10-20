#!/usr/bin/env nextflow

/*
 * Sample QC: cohort module
 * The family tables in one, with the relatedness between families
 */

nextflow.enable.dsl=2

process SAMPLEQC_COHORT {
  /*
  Gather the families' QC tables, and find the relatives the pedigree puts in other families

  somalier relate of every extract of the cohort, with the whole pedigree; it compares all
  pairs, about a minute for 4,500 samples, and its HTML is kept for browsing them. Its
  pairs.tsv holds every pair the pedigree relates and every pair of relatedness above 0.05,
  plus, past 200,000 pairs, a random sample of the others, which sampleqc_classify cohort
  drops as unrelated. sampleqc_classify cohort then writes the family tables in one, the
  pairs of different families related at the third degree or closer added, with the flags
  only the cohort can raise: CRYPTIC_RELATED, DUPLICATE_OR_MZ across families,
  SAMPLE_SWAP_CANDIDATE and RELATEDNESS_UNRELIABLE. relate runs with sampleqc_hom_cutoff as
  its homozygous cutoff, as the families' does; the unrelated pairs of its pairs.tsv are the
  IBS0 yardstick a second-degree or closer pair across families is checked against.

  The snapshot is the cohort's, all its families' rows: an added, removed or edited family
  changes it, and the plan then rebuilds these tables.

  Parameters
  ----------
  cohort : val
    Cohort name
  extracts : path
    Every .somalier file of the cohort
  tables : path
    Every family's samples, pairs and flags tables
  snapshot : val
    The cohort's inputs snapshot, as text
  sites : path
    somalier's sites VCF
  label : val
    sampleqc_label, identifying the QC in output file names
  classify : path
    sampleqc_classify, staged in so it is reachable inside the container

  Returns
  -------
  Tuple of cohort name, its samples, pairs and flags tables, the HTML and the snapshot
  */

  tag "${cohort}"

  publishDir "${params.data}/cohorts/${cohort}/qc",
    mode: 'copy',
    pattern: "${cohort}.${label}.{samples.tsv,pairs.tsv,flags.tsv,inputs.tsv,html}"

  label 'sampleqc_cohort'

  input:
  tuple val(cohort), path(extracts, stageAs: 'extracts/*'), path(tables, stageAs: 'families/*'),
        val(snapshot), path(sites), val(label), path(classify)

  output:
  tuple val(cohort), path("${prefix}.samples.tsv"), path("${prefix}.pairs.tsv"), path("${prefix}.flags.tsv"),
        path("${prefix}.html"), path("${prefix}.inputs.tsv"), emit: tables

  script:
  prefix = "${cohort}.${label}"

  """
  set -euo pipefail

  ${InputsSnapshot.decodeCommand(snapshot, "${prefix}.inputs.tsv")}

  awk -F'\\t' -v OFS='\\t' '!/^#/ && \$6 == "yes" { print \$1, \$2, \$3, \$4, \$5, -9 }' ${prefix}.inputs.tsv > cohort.ped
  # A list rather than the paths themselves, which a large cohort would take past ARG_MAX
  find extracts -name '*.somalier' | sort > extracts.list
  echo "${cohort}: somalier relate of \$(wc -l < extracts.list | tr -d ' ') extracts started at \$(date)"
  SOMALIER_AB_HOM_CUTOFF=${params.sampleqc_hom_cutoff} somalier relate --sites ${sites} --ped cohort.ped -o relate extracts.list
  mv relate.html ${prefix}.html

  python3 ${classify} cohort \\
      --cohort ${cohort} \\
      --snapshot ${prefix}.inputs.tsv \\
      --somalier-pairs relate.pairs.tsv \\
      --family-dir families \\
      --label ${label} \\
      --out-prefix ${prefix}
  """

  stub:
  prefix = "${cohort}.${label}"
  """
  ${InputsSnapshot.decodeCommand(snapshot, "${prefix}.inputs.tsv")}
  printf 'family_id\\tsample\\tflags\\n' > ${prefix}.samples.tsv
  printf 'family_id_a\\tsample_a\\tfamily_id_b\\tsample_b\\tstatus\\n' > ${prefix}.pairs.tsv
  printf 'level\\tfamily_id\\tsample_a\\tsample_b\\tflag\\tseverity\\tevidence\\n' > ${prefix}.flags.tsv
  echo '<html></html>' > ${prefix}.html
  """
}
