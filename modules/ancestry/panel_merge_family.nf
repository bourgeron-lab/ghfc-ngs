#!/usr/bin/env nextflow

/*
 * Ancestry: family panel genotype merge module
 * Merges per-sample panel genotype BCFs into a family-level BCF
 */

nextflow.enable.dsl=2

process PANEL_MERGE_FAMILY {
  /*
  Merge a family's per-sample panel genotype BCFs into one family BCF

  Every sample already carries a record at every panel site, so this merge must
  NOT use --missing-to-ref. That flag exists for merging files with different site
  sets and it would overwrite the honest ./. of an uncallable site with 0/0 -
  discarding exactly the coverage information the gVCF extraction was done to
  recover.

  Parameters
  ----------
  fid : val
    Family ID
  sample_bcfs : list
    Per-sample panel genotype BCF files
  sample_csis : list
    Their index files
  panel_name : val
    Label identifying this panel/threshold combination in output file names
  expected_samples : val
    Number of samples the pedigree expects in this family

  Returns
  -------
  Tuple of family ID, merged panel genotype BCF, and its index
  */

  tag "$fid"

  publishDir {
    def (s1, s2) = Sharding.getShards(fid)
    "${params.data}/families/${s1}/${s2}/${fid}/ancestry"
  },
    mode: 'copy',
    pattern: "${fid}.panel_gt.${panel_name}.bcf*"

  label 'panel_merge_family'

  input:
  // Staged into a subdirectory so an input can never collide with the output file
  // name, for the same reason documented in snvs_cohort/snvs_cohort_merge.nf: a
  // family whose ID matches one of its own sample barcodes would otherwise stage an
  // input over the output and have bcftools truncate the published file.
  tuple val(fid), path(sample_bcfs, stageAs: 'input_bcfs/*'), path(sample_csis, stageAs: 'input_bcfs/*'), val(panel_name), val(expected_samples)

  output:
  tuple val(fid), path("${output_bcf}"), path("${output_bcf}.csi"), emit: family_panel_gt

  script:
  output_bcf = "${fid}.panel_gt.${panel_name}.bcf"

  """
  set -euo pipefail

  ls input_bcfs/*.bcf > bcf_file_list.txt

  if [ ! -s bcf_file_list.txt ]; then
    echo "ERROR: no per-sample panel BCFs were staged for family ${fid}" >&2
    exit 1
  fi

  n_in=\$(wc -l < bcf_file_list.txt | tr -d '[:space:]')
  if [ "\${n_in}" -ne "${expected_samples}" ]; then
    echo "ERROR: family ${fid} expects ${expected_samples} samples but \${n_in} panel BCFs were staged" >&2
    exit 1
  fi

  # Every input holds the same sites, so --merge none is a straight column join.
  # --missing-to-ref is deliberately absent; see the note above.
  bcftools merge \\
      --threads ${task.cpus} \\
      --merge none \\
      --file-list bcf_file_list.txt \\
      --output-type b \\
      --output ${output_bcf}

  bcftools index ${output_bcf}

  # A merge that silently dropped sites or samples would degrade every downstream
  # score, and bcftools can exit 0 having written only a header.
  n_sites_in=\$(bcftools index -n \$(head -n 1 bcf_file_list.txt) 2>/dev/null)
  n_sites_out=\$(bcftools index -n ${output_bcf} 2>/dev/null)
  n_samples_out=\$(bcftools query -l ${output_bcf} | wc -l | tr -d '[:space:]')
  echo "Family ${fid}: \${n_samples_out} samples x \${n_sites_out} sites"
  if [ "\${n_sites_out}" -ne "\${n_sites_in}" ]; then
    echo "ERROR: merged BCF has \${n_sites_out} sites, inputs have \${n_sites_in}" >&2
    exit 1
  fi
  if [ "\${n_samples_out}" -ne "${expected_samples}" ]; then
    echo "ERROR: merged BCF has \${n_samples_out} samples, expected ${expected_samples}" >&2
    exit 1
  fi
  """

  stub:
  output_bcf = "${fid}.panel_gt.${panel_name}.bcf"
  """
  touch ${output_bcf}
  touch ${output_bcf}.csi
  """
}
