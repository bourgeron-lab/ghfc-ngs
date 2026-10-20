#!/usr/bin/env nextflow

/*
 * ExpansionHunter: per-family merge module
 * Merges the members' ExpansionHunter VCFs and tabulates their genotypes
 */

nextflow.enable.dsl=2

process EXPANSIONHUNTER_FAMILY {
  /*
  Merge a family's ExpansionHunter VCFs into one, and write its genotype table

  ExpansionHunter's VCFs carry no ##contig lines, which bcftools needs to index and
  merge them, so each one gets those of the contigs it uses, from the reference index,
  before it is compressed. The merge keeps one record per locus (--merge all): the
  members' alleles are symbolic (<STR14>, <STR23>...) and differ from one member to
  the next.

  The table is written by eh_family_table: one row per locus and member, with each
  allele's class against the thresholds table when one is given, its confidence from its
  reads and the member's read length (from its JSON), and the row's severity.

  The snapshot, built by the plan (lib/InputsSnapshot.groovy), is written as
  <prefix>.inputs.tsv: the family's pedigree rows and the checksums of the thresholds table
  and of eh_family_table. The plan rebuilds the family when the one it builds for a run
  differs from the one on disk, so an edited threshold or a new classification reaches every
  family without its members being called again.

  Parameters
  ----------
  fid : val
    Family ID
  vcfs : path
    The members' ExpansionHunter VCFs, one sample each
  jsons : path
    The members' ExpansionHunter JSONs, for their read length
  pedigree : path
    Cohort pedigree file; only the family's rows are read
  ref_fai : path
    Index of the reference the samples were genotyped against
  thresholds : path
    Per-locus thresholds table, or an empty list to leave alleles unclassified
  catalog_name : val
    Label identifying the catalog in output file names
  family_table : path
    The table script, staged in so it is reachable inside the container
  snapshot : val
    The family's inputs snapshot, as text

  Returns
  -------
  Tuple of family ID, merged VCF and its index, and the genotype table; the snapshot is
  published beside them
  */

  tag "${fid}"

  publishDir {
    def (s1, s2) = Sharding.getShards(fid)
    "${params.data}/families/${s1}/${s2}/${fid}/str/expansionhunter"
  },
    mode: 'copy',
    pattern: "${fid}.${catalog_name}{.vcf.gz,.vcf.gz.tbi,.genotypes.tsv,.inputs.tsv}"

  label 'expansionhunter_family'

  input:
  // Staged into a subdirectory, so a family named after its only member cannot have its
  // input and its output under the same name
  tuple val(fid), path(vcfs, stageAs: 'input_vcfs/*'), path(jsons, stageAs: 'input_jsons/*'), path(pedigree),
        path(ref_fai), path(thresholds), val(catalog_name), path(family_table), val(snapshot)

  output:
  tuple val(fid), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf
  tuple val(fid), path("${prefix}.genotypes.tsv"), emit: table
  tuple val(fid), path("${prefix}.inputs.tsv"), emit: snapshot

  script:
  prefix = "${fid}.${catalog_name}"
  def thresholds_arg = thresholds ? "--thresholds ${thresholds}" : ""

  """
  set -euo pipefail

  # bcftools sort spills to TMPDIR; the host value is not mounted in the container
  export TMPDIR=\$(pwd)

  ${InputsSnapshot.decodeCommand(snapshot, "${prefix}.inputs.tsv")}

  find input_vcfs -maxdepth 1 -name '*.vcf' | sort > vcf_list.txt
  if [ ! -s vcf_list.txt ]; then
    echo "ERROR: no ExpansionHunter VCFs were staged for family ${fid}" >&2
    exit 1
  fi

  mkdir bgz
  : > bgz_list.txt
  while read -r vcf; do
    name=\$(basename "\${vcf}" .vcf)
    if grep -q '^##contig=' "\${vcf}"; then
      : > contigs.txt
    else
      awk -F'\\t' 'NR == FNR { if (\$0 !~ /^#/) used[\$1] = 1; next }
                  (\$1 in used) { print "##contig=<ID=" \$1 ",length=" \$2 ">" }' "\${vcf}" ${ref_fai} > contigs.txt
    fi
    { grep '^##' "\${vcf}"; cat contigs.txt; grep -v '^##' "\${vcf}"; } \\
      | bcftools sort -Oz -o "bgz/\${name}.vcf.gz"
    bcftools index -t "bgz/\${name}.vcf.gz"
    echo "bgz/\${name}.vcf.gz" >> bgz_list.txt
  done < vcf_list.txt

  # bcftools merge refuses a single input
  if [ "\$(wc -l < bgz_list.txt)" -eq 1 ]; then
    cp "\$(cat bgz_list.txt)" ${prefix}.vcf.gz
  else
    bcftools merge --merge all --file-list bgz_list.txt -Oz -o ${prefix}.vcf.gz
  fi
  bcftools index -t ${prefix}.vcf.gz

  # bcftools can exit 0 having written only a header
  n_out=\$(bcftools index -n ${prefix}.vcf.gz)
  echo "Family ${fid}: \$(wc -l < vcf_list.txt | tr -d '[:space:]') members -> \${n_out} merged records"
  if [ "\${n_out}" -eq 0 ]; then
    echo "ERROR: the family merge produced no records" >&2
    exit 1
  fi

  json_args=()
  while read -r json; do
    json_args+=(--json "\${json}")
  done < <(find input_jsons -maxdepth 1 -name '*.json' | sort)

  python3 ${family_table} \\
      --family ${fid} \\
      --pedigree ${pedigree} \\
      ${thresholds_arg} \\
      \${json_args[@]+"\${json_args[@]}"} \\
      --out ${prefix}.genotypes.tsv \\
      \$(cat vcf_list.txt)
  """

  stub:
  prefix = "${fid}.${catalog_name}"
  """
  ${InputsSnapshot.decodeCommand(snapshot, "${prefix}.inputs.tsv")}
  echo | gzip > ${prefix}.vcf.gz
  touch ${prefix}.vcf.gz.tbi
  printf 'family_id\\tsample\\tfather\\tmother\\tlocus\\tallele1\\tallele2\\tstatus\\tseverity\\tfilter\\n${fid}\\t${fid}_s1\\t.\\t.\\tHTT\\t19\\t19\\tnormal\\t.\\tPASS\\n' > ${prefix}.genotypes.tsv
  """
}
