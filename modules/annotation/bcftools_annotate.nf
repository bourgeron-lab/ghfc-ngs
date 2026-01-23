#!/usr/bin/env nextflow

/*
 * Annotation: BCFtools annotate module for adding additional annotations from BCF files
 */

nextflow.enable.dsl=2

process BCFTOOLS_ANNOTATE {
  /*
  Annotate VCF.gz files with additional annotations from BCF files using bcftools annotate
  Transfers complete INFO fields from annotation BCF files

  Parameters
  ----------
  fid : str
    Family ID
  vcf : path
    Input VCF.gz file (e.g., FID.rare.vep_config_name.vcf.gz)
  vcf_index : path
    VCF index file (.tbi)
  annotation_path : val
    Directory path containing annotation BCF files
  annotation_list : val
    List of BCF filenames to use for annotation

  Returns
  -------
  Tuple of family ID, annotated BCF file (e.g., FID.rare.vep_config_name.annotated.bcf), and its index
  */

  tag "$fid"

  publishDir "${params.data}/families/${fid}/vcfs",
    mode: 'copy',
    pattern: "*.annotated.bcf*"

  container 'docker://staphb/bcftools:latest'
  
  label 'bcftools_annotate'

  input:
  tuple val(fid), path(vcf), path(vcf_index), val(annotation_path), val(annotation_list)

  output:
  tuple val(fid), path("${output_bcf}"), path("${output_bcf}.csi"), emit: annotated_bcfs

  script:
  // Extract base filename without .vcf.gz extension and add .annotated
  base_name = vcf.baseName.replaceAll(/\.vcf$/, '')
  output_bcf = "${base_name}.annotated.bcf"
  
  // Build the annotation command chain
  def annotation_commands = annotation_list.collect { annot_file ->
    // Check if file is bed.gz or bcf
    def is_bed = annot_file.endsWith('.bed.gz')
    
    if (is_bed) {
      // Handle bed.gz files (e.g., LCR regions)
      def annotation_name = annot_file.replaceAll(/\.grch38\.bed\.gz$/, '').replaceAll(/\.bed\.gz$/, '').toUpperCase()
      """
      echo "Adding ${annotation_name} annotation"
      echo '##INFO=<ID=${annotation_name},Number=0,Type=Flag,Description="Variant falls in ${annotation_name} region">' > ${annotation_name}_header.txt
      
      bcftools annotate \\
        --threads ${task.cpus} \\
        -a ${annotation_path}/${annot_file} \\
        -h ${annotation_name}_header.txt \\
        -c CHROM,FROM,TO \\
        --mark-sites +${annotation_name} \\
        -O b \\
        -o tmp_annotated.bcf \\
        input.bcf
      
      mv tmp_annotated.bcf input.bcf
      bcftools index input.bcf
      """
    } else {
      // Handle BCF files
      """
      bcftools annotate \\
        --threads ${task.cpus} \\
        -a ${annotation_path}/${annot_file} \\
        -c INFO \\
        -O b \\
        -o tmp_annotated.bcf \\
        input.bcf
      
      mv tmp_annotated.bcf input.bcf
      bcftools index input.bcf
      """
    }
  }.join('\n\n')

  """
  # Convert input VCF.gz to BCF for processing
  bcftools view \\
    -O b \\
    -o input.bcf \\
    ${vcf}
  
  bcftools index input.bcf

  # Apply annotations from each BCF file sequentially
  ${annotation_commands}

  # Rename final output
  mv input.bcf ${output_bcf}
  mv input.bcf.csi ${output_bcf}.csi
  """

  stub:
  base_name = vcf.baseName.replaceAll(/\.vcf$/, '')
  output_bcf = "${base_name}.annotated.bcf"
  """
  touch ${output_bcf}
  touch ${output_bcf}.csi
  """
}
