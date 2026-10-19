#!/usr/bin/env nextflow

/*
 * Sample QC: per-sample extraction module
 * The sample's allele counts at somalier's sites, from its CRAM
 */

nextflow.enable.dsl=2

process SOMALIER_EXTRACT {
  /*
  Count the reference and alternate reads of one CRAM at somalier's sites

  somalier extract reads the CRAM at every site of the sites VCF, 17,384 autosomal SNPs of
  allele frequency near 0.5, 365 on chrX and 17 on chrY in hg38's, and writes a
  few kB of counts that every somalier relate of the family and of the cohort reads. The
  counts are independent of variant calling, and of the gVCF the family calls come from.

  somalier names the sample after the CRAM's read group SM, which an old CRAM may have
  set to another name. The name stored in the file is set to the barcode, and the change
  logged, so that relate and the pedigree agree on every sample's name.

  Parameters
  ----------
  barcode : val
    Sample ID
  cram : path
    Sample CRAM file
  crai : path
    CRAM index file
  ref : path
    The reference the CRAM was aligned to
  ref_fai : path
    Its index, staged beside it
  sites : path
    somalier's sites VCF for the reference (somalier_sites)
  label : val
    sampleqc_label, identifying the extraction in output file names

  Returns
  -------
  Tuple of barcode and its .somalier file
  */

  tag "${barcode}"

  publishDir {
    def (s1, s2) = Sharding.getShards(barcode)
    "${params.data}/samples/${s1}/${s2}/${barcode}/qc/somalier"
  },
    mode: 'copy',
    pattern: "${barcode}.${label}.somalier"

  label 'somalier_extract'

  input:
  tuple val(barcode), path(cram), path(crai), path(ref), path(ref_fai), path(sites), val(label)

  output:
  tuple val(barcode), path("${barcode}.${label}.somalier"), emit: extract

  script:
  """
  set -euo pipefail

  somalier extract -d extracted --sites ${sites} -f ${ref} ${cram}

  n=\$(find extracted -name '*.somalier' | wc -l | tr -d ' ')
  if [ "\${n}" -ne 1 ]; then
    echo "ERROR: somalier wrote \${n} extracts for ${barcode}'s CRAM, not one: its read groups name \${n} samples" >&2
    exit 1
  fi

  # A .somalier file starts with its format version, the length of the sample name, and the name
  python3 - extracted/*.somalier ${barcode} ${barcode}.${label}.somalier <<'PY'
import sys
src, barcode, dst = sys.argv[1:]
data = open(src, "rb").read()
if data[0] != 2:
    sys.exit("ERROR: {} is somalier format {}, not 2".format(src, data[0]))
name = data[2:2 + data[1]].decode()
if name != barcode:
    print("{}: the CRAM names its sample {}, stored as {}".format(barcode, name, barcode))
new = barcode.encode()
open(dst, "wb").write(bytes([2, len(new)]) + new + data[2 + data[1]:])
PY
  """

  stub:
  """
  echo stub > ${barcode}.${label}.somalier
  """
}
