#!/usr/bin/env nextflow

/*
 * Wombat: De novo mutation extraction module using slivar
 */

nextflow.enable.dsl=2

process DNM_EXTRACTION {
  /*
  Extract de novo mutations from family VCF files using slivar
  Filters variants based on call rate, depth, genotype quality, and variant allele frequency

  Parameters
  ----------
  fid : str
    Family ID
  bcf : path
    Input BCF file (e.g., FID.rare.vep_config_name.annotated.bcf)
  bcf_index : path
    BCF index file (.csi)
  pedigree : path
    Family pedigree file (FID.pedigree.tsv)
  min_callrate : val
    Minimum call rate threshold (default: 0.9)
  min_dp : val
    Minimum depth (DP) threshold (default: 10)
  min_gq : val
    Minimum genotype quality (GQ) threshold (default: 20)
  min_vaf : val
    Minimum variant allele frequency (VAF) threshold (default: 0.2)

  Returns
  -------
  Tuple of family ID, de novo BCF file (FID.rare.vep_config_name.annotated.dnm.bcf), 
  its index, and TSV file (FID.rare.vep_config_name.annotated.dnm.tsv)
  */

  tag "$fid"

  publishDir "${params.data}/families/${fid}/vcfs",
    mode: 'copy',
    pattern: "*.dnm.{bcf,bcf.csi,tsv}"

  container 'docker://brentp/slivar:latest'
  
  label 'dnm_extraction'

  input:
  tuple val(fid), path(bcf), path(bcf_index), path(pedigree), val(min_callrate), val(min_dp), val(min_gq), val(min_vaf)

  output:
  tuple val(fid), path("${output_bcf}"), path("${output_bcf}.csi"), path("${output_tsv}"), emit: dnm_results

  script:
  // Extract base filename without .bcf extension and add .dnm
  base_name = bcf.baseName
  output_bcf = "${base_name}.dnm.bcf"
  output_tsv = "${base_name}.dnm.tsv"

  """
  # Check if pedigree file has a header (first column contains "FID")
  # If it does, remove the header and create a temp file for slivar
  if head -n 1 ${pedigree} | cut -f1 | grep -q "^FID\$"; then
    echo "Pedigree file has header, removing it for slivar"
    tail -n +2 ${pedigree} > temp_pedigree.tsv
    PED_FILE="temp_pedigree.tsv"
  else
    echo "Pedigree file has no header, using as is"
    PED_FILE="${pedigree}"
  fi

  # Extract de novo variants using slivar
  # The slivar expr command applies filters for de novo variants
  # Filters:
  # - Call rate >= ${min_callrate}
  # - DP >= ${min_dp} in child
  # - GQ >= ${min_gq} in child
  # - VAF (AB - allele balance) >= ${min_vaf} in child
  #
  # De novo patterns:
  # 1. Autosomal het: kid.het && mom.hom_ref && dad.hom_ref (classic de novo)
  # 2. Autosomal hom_alt (hemizygous): kid.hom_alt && mom.hom_ref && dad.hom_ref (e.g., deletion)
  # 3. X chromosome in males: kid.hom_alt && mom.hom_ref (kid.sex == 1)
  # 4. Y chromosome in males: kid.hom_alt && dad.hom_ref (kid.sex == 1)
  # Note: For hom_ref parents, we also enforce AB <= 0.02 to ensure true homozygous reference
  
  slivar expr \\
    --vcf ${bcf} \\
    --ped \$PED_FILE \\
    --pass-only \\
    --info "variant.call_rate >= ${min_callrate}" \\
    --out-vcf temp_output.bcf \\
    --trio "denovo_auto_het:(kid.het && mom.hom_ref && mom.AB <= 0.02 && dad.hom_ref && dad.AB <= 0.02 && kid.DP >= ${min_dp} && kid.GQ >= ${min_gq} && kid.AB >= ${min_vaf})" \\
    --trio "denovo_auto_hom:(kid.hom_alt && kid.AB >= 0.98 && mom.hom_ref && mom.AB <= 0.02 && dad.hom_ref && dad.AB <= 0.02 && kid.DP >= ${min_dp} && kid.GQ >= ${min_gq})" \\
    --trio "denovo_x_male:(variant.CHROM == 'chrX' && kid.sex == 1 && kid.hom_alt && mom.hom_ref && mom.AB <= 0.02 && kid.DP >= ${min_dp} && kid.GQ >= ${min_gq})" \\
    --trio "denovo_y_male:(variant.CHROM == 'chrY' && kid.sex == 1 && kid.hom_alt && dad.hom_ref && dad.AB <= 0.02 && kid.DP >= ${min_dp} && kid.GQ >= ${min_gq})"

  # Convert output to BCF format
  bcftools view -O b -o ${output_bcf} temp_output.bcf

  # Index the output BCF
  bcftools index ${output_bcf}

  # Generate TSV output from the BCF file
  slivar tsv \\
    --ped \$PED_FILE \\
    --sample-field denovo_auto_het \\
    --sample-field denovo_auto_hom \\
    --sample-field denovo_x_male \\
    --sample-field denovo_y_male \\
    --info-field denovo_auto_het \\
    --info-field denovo_auto_hom \\
    --info-field denovo_x_male \\
    --info-field denovo_y_male \\
    --csq-field CSQ \\
    ${output_bcf} > ${output_tsv}
  """

  stub:
  base_name = bcf.baseName
  output_bcf = "${base_name}.dnm.bcf"
  output_tsv = "${base_name}.dnm.tsv"
  """
  touch ${output_bcf}
  touch ${output_bcf}.csi
  touch ${output_tsv}
  """
}
