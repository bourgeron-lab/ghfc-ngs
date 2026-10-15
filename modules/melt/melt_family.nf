#!/usr/bin/env nextflow

/*
 * MELT: per-family merge module
 * Merges the members' MELT VCFs, all three classes, into one family VCF
 * Tabulates each child's insertions against its parents' calls
 */

nextflow.enable.dsl=2

process MELT_FAMILY {
  /*
  Merge a family's MELT calls into one VCF, and write its transmission table

  Each member's ALU, LINE1 and SVA VCFs are sorted, compressed and concatenated into one,
  then the members are merged with --merge none: an insertion is one record for every
  member called at the same position with the same class, and a separate record
  otherwise. MELT calls each sample on its own, so the same insertion is often placed a
  few bases apart in two members and then makes two records; a member without a record
  has no call there, which is not proof that it lacks the insertion.

  The transmission table is written by melt_family_transmission from the members' own
  VCFs, not from the merge: for every PASS insertion of a child with a genotyped parent,
  each parent's nearest call of the same class within 50 bp, PASS or not, and whether the
  insertion is inherited, absent in both parents (a de novo candidate), or cannot tell
  because a parent is not genotyped.

  Parameters
  ----------
  fid : val
    Family ID
  vcfs : path
    The members' MELT VCFs, three per member
  pedigree : path
    Cohort pedigree file; only the family's rows are read
  label : val
    melt_label, identifying the calls in output file names
  transmission : path
    The transmission table script, staged in so it is reachable inside the container

  Returns
  -------
  Tuple of family ID, merged VCF and its index, and the transmission table
  */

  tag "${fid}"

  publishDir {
    def (s1, s2) = Sharding.getShards(fid)
    "${params.data}/families/${s1}/${s2}/${fid}/mei/melt"
  },
    mode: 'copy',
    pattern: "${fid}.${label}{.vcf.gz,.vcf.gz.tbi,.transmission.tsv}"

  label 'melt_family'

  input:
  // Staged into a subdirectory, so a family named after its only member cannot have its
  // input and its output under the same name
  tuple val(fid), path(vcfs, stageAs: 'input_vcfs/*'), path(pedigree), val(label), path(transmission)

  output:
  tuple val(fid), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf
  tuple val(fid), path("${prefix}.transmission.tsv"), emit: transmission

  script:
  prefix = "${fid}.${label}"

  """
  set -euo pipefail

  # bcftools sort spills to TMPDIR; the host value is not mounted in the container
  export TMPDIR=\$(pwd)

  find input_vcfs -maxdepth 1 -name '*.${label}.ALU.vcf' | sort > alu_list.txt
  if [ ! -s alu_list.txt ]; then
    echo "ERROR: no MELT VCFs were staged for family ${fid}" >&2
    exit 1
  fi

  mkdir bgz
  : > member_list.txt
  while read -r alu; do
    member=\$(basename "\${alu}" .${label}.ALU.vcf)
    for cls in ALU LINE1 SVA; do
      vcf="input_vcfs/\${member}.${label}.\${cls}.vcf"
      if [ ! -s "\${vcf}" ]; then
        echo "ERROR: \${member} has no ${label} \${cls} VCF" >&2
        exit 1
      fi
      bcftools sort -Oz -o "bgz/\${member}.\${cls}.vcf.gz" "\${vcf}"
      bcftools index -t "bgz/\${member}.\${cls}.vcf.gz"
    done
    bcftools concat --allow-overlaps -Oz -o "bgz/\${member}.vcf.gz" \\
      "bgz/\${member}.ALU.vcf.gz" "bgz/\${member}.LINE1.vcf.gz" "bgz/\${member}.SVA.vcf.gz"
    bcftools index -t "bgz/\${member}.vcf.gz"
    echo "bgz/\${member}.vcf.gz" >> member_list.txt
  done < alu_list.txt

  # bcftools merge refuses a single input
  if [ "\$(wc -l < member_list.txt)" -eq 1 ]; then
    cp "\$(cat member_list.txt)" ${prefix}.vcf.gz
  else
    bcftools merge --merge none --file-list member_list.txt -Oz -o ${prefix}.vcf.gz
  fi
  bcftools index -t ${prefix}.vcf.gz

  # bcftools can exit 0 having written only a header
  n_out=\$(bcftools index -n ${prefix}.vcf.gz)
  echo "Family ${fid}: \$(wc -l < member_list.txt | tr -d '[:space:]') members -> \${n_out} merged records"
  if [ "\${n_out}" -eq 0 ]; then
    echo "ERROR: the family merge produced no records" >&2
    exit 1
  fi

  python3 ${transmission} \\
      --family ${fid} \\
      --pedigree ${pedigree} \\
      --input-dir input_vcfs \\
      --label ${label} \\
      --out ${prefix}.transmission.tsv
  """

  stub:
  prefix = "${fid}.${label}"
  """
  echo | gzip > ${prefix}.vcf.gz
  touch ${prefix}.vcf.gz.tbi
  printf 'family_id\\tchild\\tfather\\tmother\\tmei_class\\tchrom\\tpos\\ttransmission\\n${fid}\\t${fid}_s1\\t${fid}_s2\\t${fid}_s3\\tALU\\tchr1\\t1000\\tinherited_father\\n' > ${prefix}.transmission.tsv
  """
}
