#!/usr/bin/env nextflow

/*
 * MELT: per-sample mobile element insertion calling module
 * Calls ALU, LINE1 and SVA insertions in the sample's CRAM with MELT Single
 */

nextflow.enable.dsl=2

process MELT_CALL {
  /*
  Call one sample's ALU, LINE1 and SVA insertions with MELT v2.2.2 Single

  MELT is not in the image (its licence forbids redistributing it; see
  containers/melt/Dockerfile): the distribution is staged in from params.melt_dir, and
  this runs its MELT.jar, transposon references and gene annotation.

  MELT names the sample after the reads file, so the CRAM is linked as <barcode>.cram, as
  for ExpansionHunter. MELT writes its intermediates, about 20 GB, beside that file, which
  is why the process runs in node-local scratch.

  The htsjdk bundled in MELT 2.2.2 cannot decode CRAM 3.1, so a 3.1 CRAM is transcoded
  to 3.0 first, in the task directory; a 2.x or 3.0 one is read as it is.

  The CRAM's contigs must be the reference's. A CRAM aligned to another naming (1 rather
  than chr1) fails at once with exit status 64, which is never retried, rather than after
  MELT has scanned it. Java runs with -XX:+ExitOnOutOfMemoryError, so an exhausted heap
  ends the task with status 3, which is retried with more memory, instead of leaving
  MELT's worker threads to hang; the heap is 80% of the task's memory.

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
  melt_dir : path
    The MELT v2.2.2 distribution: MELT.jar, me_refs/Hg38/, add_bed_files/Hg38/
  label : val
    melt_label, identifying the calls in output file names

  Returns
  -------
  Tuple of barcode and its ALU, LINE1 and SVA VCFs, as MELT writes them
  */

  tag "${barcode}"

  publishDir {
    def (s1, s2) = Sharding.getShards(barcode)
    "${params.data}/samples/${s1}/${s2}/${barcode}/mei/melt"
  },
    mode: 'copy',
    pattern: "${barcode}.${label}.{ALU,LINE1,SVA}.vcf"

  label 'melt_call'

  input:
  tuple val(barcode), path(cram), path(crai), path(ref), path(ref_fai), path(melt_dir), val(label)

  output:
  tuple val(barcode), path("${prefix}.ALU.vcf"), path("${prefix}.LINE1.vcf"), path("${prefix}.SVA.vcf"), emit: vcfs

  script:
  prefix = "${barcode}.${label}"
  def heap_mb = Math.max(1024, (long) (task.memory.toMega() * 0.8))

  """
  set -euo pipefail

  # Every CRAM contig must be one of the reference's
  samtools view -H ${cram} | awk -F'\\t' '\$1 == "@SQ" { for (i = 2; i <= NF; i++) if (\$i ~ /^SN:/) print substr(\$i, 4) }' \\
    | sort -u > cram_contigs.txt
  cut -f1 ${ref_fai} | sort -u > ref_contigs.txt
  absent=\$(comm -23 cram_contigs.txt ref_contigs.txt | head -5 | paste -sd, -)
  if [ ! -s cram_contigs.txt ] || [ -n "\${absent}" ]; then
    echo "ERROR: ${barcode}'s CRAM is not aligned to ${ref}: contigs \${absent:-none} are not in it" >&2
    exit 64
  fi

  # CRAM major and minor version, bytes 5 and 6 of the file
  version=\$(od -An -tu1 -j4 -N2 ${cram} | awk '{ print \$1 "." \$2 }')
  case "\${version}" in
    2.0|2.1|3.0)
      if [ "${cram}" != "${barcode}.cram" ]; then
        ln -s "${cram}" "${barcode}.cram"
        ln -s "${crai}" "${barcode}.cram.crai"
      fi
      ;;
    *)
      echo "${barcode}: CRAM \${version}, transcoded to 3.0 for MELT's htsjdk"
      samtools view -@ ${task.cpus} -T ${ref} -C --output-fmt-option version=3.0 -o ${barcode}.cram ${cram}
      samtools index -@ ${task.cpus} ${barcode}.cram
      ;;
  esac

  # MELT's -r: the longest read among the first thousand primary ones, 150 if none is read.
  # head closes the pipe early, which pipefail would take for a failure
  readlen=\$( (samtools view -F 0x900 -T ${ref} ${barcode}.cram 2>/dev/null || true) \\
    | head -n 1000 | awk '{ if (length(\$10) > m) m = length(\$10) } END { print m + 0 }' || true)
  [ "\${readlen:-0}" -gt 0 ] || readlen=150
  echo "${barcode}: read length \${readlen}"

  for cls in ALU LINE1 SVA; do
    echo "\$(pwd)/${melt_dir}/me_refs/Hg38/\${cls}_MELT.zip"
  done > mei_list.txt

  java -Xmx${heap_mb}m -XX:+ExitOnOutOfMemoryError -jar ${melt_dir}/MELT.jar Single \\
      -bamfile "\$(pwd)/${barcode}.cram" \\
      -h ${ref} \\
      -n ${melt_dir}/add_bed_files/Hg38/Hg38.genes.bed \\
      -t mei_list.txt \\
      -r "\${readlen}" \\
      -w "\$(pwd)" \\
      -samtools "\$(command -v samtools)"

  for cls in ALU LINE1 SVA; do
    if [ ! -s "\${cls}.final_comp.vcf" ]; then
      echo "ERROR: MELT did not write \${cls}.final_comp.vcf" >&2
      exit 1
    fi
    sample=\$(awk -F'\\t' '/^#CHROM/ { print \$10; exit }' "\${cls}.final_comp.vcf")
    if [ "\${sample}" != "${barcode}" ]; then
      echo "ERROR: MELT's \${cls} VCF names its sample '\${sample}', not ${barcode}" >&2
      exit 1
    fi
    mv "\${cls}.final_comp.vcf" "${prefix}.\${cls}.vcf"
    echo "${barcode}: \$(grep -vc '^#' "${prefix}.\${cls}.vcf" || true) \${cls} insertions, \$(awk -F'\\t' '!/^#/ && \$7 == "PASS"' "${prefix}.\${cls}.vcf" | wc -l | tr -d ' ') PASS"
  done
  """

  stub:
  prefix = "${barcode}.${label}"
  """
  for cls in ALU LINE1 SVA; do
    printf '##fileformat=VCFv4.2\\n##source=MELTv2.2.2\\n##contig=<ID=chr1,length=248956422>\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${barcode}\\n' > ${prefix}.\${cls}.vcf
    printf 'chr1\\t1000\\t.\\tA\\t<INS:ME:%s>\\t.\\tPASS\\tASSESS=5;INTERNAL=null,null;SVTYPE=%s\\tGT\\t0/1\\n' "\${cls}" "\${cls}" >> ${prefix}.\${cls}.vcf
  done
  """
}
