#!/usr/bin/env nextflow

/*
 * ExpansionHunter: per-family REViewer module
 * Draws the read pileups of the members' flagged repeat genotypes
 */

nextflow.enable.dsl=2

process EXPANSIONHUNTER_REVIEWER {
  /*
  Draw a REViewer read pileup for each locus where a member's genotype is flagged

  One task per family rather than per sample: most members have nothing flagged, and a
  pileup takes about a second, so a task per sample would be mostly scheduling.

  A genotype is flagged when its severity in the family table is high or medium: the rows
  of the cohort's flags.tsv. Intermediate alleles, carriers and risk alleles are not drawn. The table names variants (VariantId),
  and REViewer takes loci (LocusId). They differ only at loci with several variants, which
  variant_loci maps: FXN_A and FXN are both drawn as FXN, in one pileup.

  REViewer works from the realigned reads ExpansionHunter keeps, which must be sorted and
  indexed, and from the catalog and reference the call was made with. Its SVGs declare a
  width and a height but no viewBox, so a browser shows them at full size, a corner of the
  pileup filling the window; the viewBox is added, and the width set to 100%.

  Each member's files are written in a directory named after it, which publishing maps to
  the member's sample directory. Its metrics table is written even when nothing is flagged,
  with its header alone: it is what tells a later run that the member has been reviewed.

  Parameters
  ----------
  fid : val
    Family ID
  barcodes : val
    The members to review
  vcfs : path
    Their ExpansionHunter VCFs, <barcode>.<catalog_name>.vcf
  bams : path
    The realigned reads ExpansionHunter wrote, sorted, <barcode>.<catalog_name>_realigned.bam
  bais : path
    Their indexes
  family_table : path
    The family's genotype table
  ref : path
    The reference the members were genotyped against
  ref_fai : path
    Its index, staged beside it
  catalog : path
    The variant catalog the members were genotyped against
  catalog_name : val
    Label identifying the catalog in output file names
  variant_loci : val
    VariantId=LocusId pairs, comma separated, of the catalog's variants named otherwise
    than their locus

  Returns
  -------
  Tuple of family ID and the members' SVGs, of family ID and their metrics tables, and of
  family ID and their tables of the haplotype pairs REViewer weighed; a member with nothing
  flagged has no SVG and no phasing table
  */

  tag "${fid}"

  publishDir "${params.data}/samples",
    mode: 'copy',
    saveAs: { name ->
      def (barcode, file_name) = name.tokenize('/')
      def (s1, s2) = Sharding.getShards(barcode)
      "${s1}/${s2}/${barcode}/str/expansionhunter/reviewer/${file_name}"
    }

  label 'expansionhunter_reviewer'

  input:
  // Staged into a subdirectory, apart from the member directories the outputs go to
  tuple val(fid), val(barcodes), path(vcfs, stageAs: 'input/*'), path(bams, stageAs: 'input/*'), path(bais, stageAs: 'input/*'), path(family_table), path(ref), path(ref_fai), path(catalog), val(catalog_name), val(variant_loci)

  output:
  tuple val(fid), path("*/*.svg"), emit: svg, optional: true
  tuple val(fid), path("*/*.metrics.tsv"), emit: metrics
  tuple val(fid), path("*/*.phasing.tsv"), emit: phasing, optional: true

  script:
  """
  set -euo pipefail

  for barcode in ${barcodes.join(' ')}; do
    prefix=\${barcode}/\${barcode}.${catalog_name}
    mkdir -p \${barcode}

    loci=\$(awk -F'\\t' -v sample=\${barcode} -v variant_loci='${variant_loci}' '
      BEGIN {
        n = split(variant_loci, pairs, ",")
        for (i = 1; i <= n; i++) { split(pairs[i], kv, "="); locus_of[kv[1]] = kv[2] }
      }
      FNR == 1 {
        for (i = 1; i <= NF; i++) col[\$i] = i
        if (!("sample" in col) || !("locus" in col) || !("severity" in col)) {
          print "ERROR: ${family_table} has no sample, locus or severity column" > "/dev/stderr"
          bad = 1
          exit 1
        }
        next
      }
      \$col["sample"] == sample && \$col["severity"] ~ /^(high|medium)\$/ {
        locus = (\$col["locus"] in locus_of) ? locus_of[\$col["locus"]] : \$col["locus"]
        if (!(locus in seen)) { seen[locus] = 1; out = out (out == "" ? "" : ",") locus }
      }
      END {
        if (bad) exit 1
        print out
      }' ${family_table})

    if [ -z "\${loci}" ]; then
      echo "\${barcode}: no flagged genotype to review"
      printf 'VariantId\\tGenotype\\tAlleleDepth\\n' > \${prefix}.metrics.tsv
      continue
    fi
    echo "\${barcode}: REViewer on \${loci}"

    REViewer \\
        --reads input/\${barcode}.${catalog_name}_realigned.bam \\
        --vcf input/\${barcode}.${catalog_name}.vcf \\
        --reference ${ref} \\
        --catalog ${catalog} \\
        --locus "\${loci}" \\
        --output-prefix \${prefix}

    for locus in \$(echo "\${loci}" | tr ',' ' '); do
      svg=\${prefix}.\${locus}.svg
      if [ ! -s "\${svg}" ]; then
        echo "ERROR: REViewer did not draw \${svg}" >&2
        exit 1
      fi
      sed -E -i '1s/^<svg width="([0-9.]+)" height="([0-9.]+)"/<svg width="100%" viewBox="0 0 \\1 \\2" preserveAspectRatio="xMinYMin meet"/' "\${svg}"
    done
    if [ ! -s \${prefix}.metrics.tsv ]; then
      echo "ERROR: REViewer did not write \${prefix}.metrics.tsv" >&2
      exit 1
    fi
  done
  """

  stub:
  """
  for barcode in ${barcodes.join(' ')}; do
    mkdir -p \${barcode}
    printf 'VariantId\\tGenotype\\tAlleleDepth\\n' > \${barcode}/\${barcode}.${catalog_name}.metrics.tsv
  done
  """
}
