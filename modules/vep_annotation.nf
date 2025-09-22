process VEP_ANNOTATION {

    tag "$fid"

    publishDir "${params.data}/families/${fid}/vcfs", mode: 'copy'

    input:
    tuple val(fid), path(vcf), path(tbi)

    output:
    tuple val(fid), path("${fid}.${params.vep_config_name}.vcf.gz"), path("${fid}.${params.vep_config_name}.vcf.gz.tbi"), emit: annotated_vcf

    script:
    """
    # Run VEP annotation
    vep \\
        --input_file ${vcf} \\
        --output_file ${fid}.${params.vep_config_name}.vcf \\
        --config ${params.vep_config} \\
        --compress_output bgzip \\
        --fork 4 \\
        --format vcf \\
        --vcf \\
        --force_overwrite

    # Index the annotated VCF
    tabix -f -p vcf ${fid}.${params.vep_config_name}.vcf.gz
    """

    stub:
    """
    touch ${fid}.${params.vep_config_name}.vcf.gz
    touch ${fid}.${params.vep_config_name}.vcf.gz.tbi
    """
}