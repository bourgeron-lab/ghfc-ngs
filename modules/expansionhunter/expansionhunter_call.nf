#!/usr/bin/env nextflow

/*
 * ExpansionHunter: per-sample short tandem repeat genotyping module
 * Genotypes every locus of the variant catalog from the sample's CRAM
 */

nextflow.enable.dsl=2

process EXPANSIONHUNTER_CALL {
  /*
  Genotype the catalog's repeats in one sample with ExpansionHunter v5.0.0

  ExpansionHunter names the sample after the reads file, without its last extension, so
  the CRAM is linked as <barcode>.cram first: read as <barcode>.GRCh38_GIABv3.cram the
  sample would be called "<barcode>.GRCh38_GIABv3" in the VCF and the JSON, and every merge
  and every join with the pedigree would need it cleaned up.

  CRAM 3.1 is read directly. The image's ExpansionHunter is linked against htslib 1.22;
  Illumina's release binary predates htslib 1.14 and cannot decode CRAM 3.1 slices.

  Parameters
  ----------
  barcode : val
    Sample ID
  cram : path
    Sample CRAM file
  crai : path
    CRAM index file
  sex : val
    'male', 'female' or empty. Empty omits --sex, and ExpansionHunter then genotypes the
    sample as female: two chrX alleles. Only the chrX loci depend on it.
  ref : path
    The reference the CRAM was aligned to, needed to decode it
  ref_fai : path
    Its index, staged beside it
  catalog : path
    ExpansionHunter variant catalog (JSON)
  catalog_name : val
    Label identifying the catalog in output file names

  Returns
  -------
  Tuple of barcode, VCF and JSON, and the realigned reads as a sorted, indexed BAM
  */

  tag "${barcode}"

  publishDir {
    def (s1, s2) = Sharding.getShards(barcode)
    "${params.data}/samples/${s1}/${s2}/${barcode}/str/expansionhunter"
  },
    mode: 'copy',
    pattern: "${barcode}.${catalog_name}{.vcf,.json,_realigned.bam,_realigned.bam.bai}"

  label 'expansionhunter_call'

  input:
  tuple val(barcode), path(cram), path(crai), val(sex), path(ref), path(ref_fai), path(catalog), val(catalog_name)

  output:
  tuple val(barcode), path("${prefix}.vcf"), emit: vcf
  tuple val(barcode), path("${prefix}.json"), emit: json
  tuple val(barcode), path("${prefix}_realigned.bam"), path("${prefix}_realigned.bam.bai"), emit: bam

  script:
  prefix = "${barcode}.${catalog_name}"
  def sex_arg = sex ? "--sex ${sex}" : ""

  """
  set -euo pipefail

  # The sample name comes from this file name (see above)
  if [ "${cram}" != "${barcode}.cram" ]; then
    ln -s "${cram}" "${barcode}.cram"
    ln -s "${crai}" "${barcode}.cram.crai"
  fi

  ExpansionHunter \\
      --reads ${barcode}.cram \\
      --reference ${ref} \\
      --variant-catalog ${catalog} \\
      --output-prefix ${prefix} \\
      --threads ${task.cpus} ${sex_arg}

  for f in ${prefix}.vcf ${prefix}.json ${prefix}_realigned.bam; do
    if [ ! -s "\$f" ]; then
      echo "ERROR: ExpansionHunter did not write \$f" >&2
      exit 1
    fi
  done

  sample=\$(awk -F'\\t' '/^#CHROM/ { print \$10; exit }' ${prefix}.vcf)
  if [ "\$sample" != "${barcode}" ]; then
    echo "ERROR: the VCF names its sample '\$sample', not ${barcode}" >&2
    exit 1
  fi
  echo "${barcode}: \$(grep -vc '^#' ${prefix}.vcf) catalog variants genotyped"

  # ExpansionHunter writes the realigned reads unsorted; sorted and indexed they open in IGV
  # (or REViewer) as they are
  samtools sort -@ ${task.cpus} -T sort_tmp -o ${prefix}_realigned.sorted.bam ${prefix}_realigned.bam
  mv ${prefix}_realigned.sorted.bam ${prefix}_realigned.bam
  samtools index ${prefix}_realigned.bam
  """

  stub:
  prefix = "${barcode}.${catalog_name}"
  """
  printf '##fileformat=VCFv4.1\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${barcode}\\n' > ${prefix}.vcf
  printf 'chr4\\t3074876\\t.\\tC\\t.\\t.\\tPASS\\tEND=3074933;REF=19;RL=57;RU=CAG;VARID=HTT;REPID=HTT\\tGT:SO:REPCN:REPCI:ADSP:ADFL:ADIR:LC\\t0/0:SPANNING/SPANNING:19/19:19-19/19-19:20/20:10/10:0/0:30.0\\n' >> ${prefix}.vcf
  echo '{}' > ${prefix}.json
  touch ${prefix}_realigned.bam ${prefix}_realigned.bam.bai
  """
}
