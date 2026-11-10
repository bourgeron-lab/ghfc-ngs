#!/usr/bin/env nextflow

/*
 * UPD: cohort module
 * The family tables in one, with the events found in other families too
 */

nextflow.enable.dsl=2

process UPD_COHORT {
  /*
  Gather the families' UPD tables

  upd_segments cohort concatenates every family's segments, events and members, and marks an
  event recurrent when events of three other families or more overlap it: a region where the
  calls are artefacts more often than UPD.

  The snapshot is the cohort's, all its families' rows: an added, removed or edited family
  changes it, and the plan then rebuilds these tables.

  Parameters
  ----------
  cohort : val
    Cohort name
  tables : path
    Every family's segments, events and members tables
  snapshot : val
    The cohort's inputs snapshot, as text
  label : val
    upd_label, identifying the analysis in output file names
  segments : path
    upd_segments, staged in so it is reachable inside the container

  Returns
  -------
  Tuple of cohort name, its segments, events and members tables, and the snapshot
  */

  tag "${cohort}"

  publishDir "${params.data}/cohorts/${cohort}/upd",
    mode: 'copy',
    pattern: "${cohort}.${label}.{segments,events,members,inputs}.tsv"

  label 'upd_cohort'

  input:
  tuple val(cohort), path(tables, stageAs: 'families/*'), val(snapshot), val(label), path(segments)

  output:
  tuple val(cohort), path("${prefix}.segments.tsv"), path("${prefix}.events.tsv"), path("${prefix}.members.tsv"),
        path("${prefix}.inputs.tsv"), emit: tables

  script:
  prefix = "${cohort}.${label}"

  """
  set -euo pipefail

  ${InputsSnapshot.decodeCommand(snapshot, "${prefix}.inputs.tsv")}

  python3 ${segments} cohort \\
      --cohort ${cohort} \\
      --snapshot ${prefix}.inputs.tsv \\
      --family-dir families \\
      --label ${label} \\
      --out-prefix ${prefix}
  """

  stub:
  prefix = "${cohort}.${label}"
  """
  ${InputsSnapshot.decodeCommand(snapshot, "${prefix}.inputs.tsv")}
  printf 'family_id\\tsample\\tchrom\\tstart\\tend\\n' > ${prefix}.segments.tsv
  printf 'family_id\\tsample\\tchrom\\tstart\\tend\\trecurrent\\n' > ${prefix}.events.tsv
  printf 'family_id\\tsample\\tmode\\tstatus\\n' > ${prefix}.members.tsv
  """
}
