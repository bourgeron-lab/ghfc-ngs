#!/usr/bin/env nextflow

/*
 * Ancestry: per-sample panel genotype extraction module
 * Genotypes the reference panel sites directly from a DeepVariant gVCF
 */

nextflow.enable.dsl=2

process PANEL_EXTRACT {
  /*
  Extract the panel sites from a sample's gVCF into a panel-shaped BCF

  A gVCF carries reference blocks, so a site with coverage and no variant can be
  reported as 0/0 while a site with no coverage stays ./.. That distinction is what
  makes a single family scorable against the panel: the family's own common_gt.bcf
  holds only sites where the family carries an alt allele, which for a trio is
  around 57% of the panel - well under the 90% per-sample coverage that admixture
  requires, and enough to distort the projected PCs.

  Parameters
  ----------
  barcode : val
    Sample barcode
  gvcf : path
    DeepVariant gVCF file
  gvcf_index : path
    gVCF index file (.tbi)
  sites_file : path
    Union site list from PANEL_SITES
  regions_file : path
    bcftools targets file from PANEL_SITES
  reference : val
    Path to the mounted ancestry-pgs reference bundle
  panel_name : val
    Label identifying this panel/threshold combination in output file names

  Returns
  -------
  Tuple of barcode, panel genotype BCF, its index, and a per-sample stats TSV
  */

  tag "$barcode"

  publishDir {
    def (s1, s2) = Sharding.getShards(barcode)
    "${params.data}/samples/${s1}/${s2}/${barcode}/ancestry"
  },
    mode: 'copy',
    pattern: "${barcode}.panel_gt.${panel_name}.{bcf,bcf.csi,stats.tsv}"

  label 'panel_extract'

  input:
  tuple val(barcode), path(gvcf), path(gvcf_index), path(sites_file), path(regions_file), val(reference), val(panel_name)

  output:
  tuple val(barcode), path("${output_bcf}"), path("${output_bcf}.csi"), emit: panel_gt
  tuple val(barcode), path("${output_stats}"), emit: panel_stats

  script:
  output_bcf = "${barcode}.panel_gt.${panel_name}.bcf"
  output_stats = "${barcode}.panel_gt.${panel_name}.stats.tsv"
  panel_genotype = "${projectDir}/modules/ancestry/scripts/panel_genotype"

  """
  set -euo pipefail

  # Contig lines come from the gVCF so the emitted VCF declares real contig
  # lengths: the downstream assembly check reads them, and a header without
  # lengths falls back to a weaker heuristic.
  bcftools view -h ${gvcf} | grep '^##contig=' > contigs.txt

  # Reference blocks report MIN_DP, not DP. Without it every block would be
  # treated as failing the depth threshold and the whole point of reading the
  # gVCF would be lost, so refuse rather than silently emit an all-missing file.
  if ! bcftools view -h ${gvcf} | grep -q '^##FORMAT=<ID=MIN_DP'; then
    echo "ERROR: ${gvcf} declares no FORMAT/MIN_DP - is this a DeepVariant gVCF?" >&2
    exit 1
  fi

  bundle_version=\$(python3 -c "import json;print(json.load(open('${reference}/manifest.json')).get('bundle_version',''))")

  # --targets-overlap 1 is REQUIRED. The bcftools default for -T is 0 ("POS in the
  # region"), which drops every spanning reference block without an error and
  # leaves only the variant records.
  bcftools query \\
      -T ${regions_file} \\
      --targets-overlap 1 \\
      -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/END\\t[%GT\\t%GQ\\t%DP\\t%MIN_DP]\\n' \\
      ${gvcf} \\
  | ${panel_genotype} \\
      --sites ${sites_file} \\
      --records - \\
      --contigs contigs.txt \\
      --sample ${barcode} \\
      --out ${barcode}.panel_gt.vcf \\
      --stats ${output_stats} \\
      --min-dp ${params.ancestry_min_dp} \\
      --min-gq ${params.ancestry_min_gq} \\
      --panel-name ${panel_name} \\
      --bundle-version "\${bundle_version}"

  bcftools view --threads ${task.cpus} -O b -o ${output_bcf} ${barcode}.panel_gt.vcf
  bcftools index ${output_bcf}
  rm -f ${barcode}.panel_gt.vcf

  # Every panel site must be present as a record, hom-ref ones included, or the
  # allele-set match downstream will drop them.
  # gzip -cd, not zcat: BSD zcat only understands .Z and fails on a .gz
  n_sites=\$(gzip -cd ${sites_file} | tail -n +2 | wc -l | tr -d '[:space:]')
  n_out=\$(bcftools index -n ${output_bcf} 2>/dev/null)
  echo "Sample ${barcode}: \${n_out} records emitted for \${n_sites} panel sites"
  if [ "\${n_out}" -ne "\${n_sites}" ]; then
    echo "ERROR: expected \${n_sites} records, got \${n_out}" >&2
    exit 1
  fi
  """

  stub:
  output_bcf = "${barcode}.panel_gt.${panel_name}.bcf"
  output_stats = "${barcode}.panel_gt.${panel_name}.stats.tsv"
  """
  touch ${output_bcf}
  touch ${output_bcf}.csi
  printf 'sample\\tpanel_name\\tn_panel_sites\\tn_called\\tcall_rate\\n${barcode}\\t${panel_name}\\t1\\t1\\t1.0\\n' > ${output_stats}
  """
}
