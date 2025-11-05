process POSTPROCESS_VARIANTS {
    
    tag "$barcode"
    
    publishDir "${params.data}/samples/${barcode}/deepvariant", pattern: "${barcode}.vcf.gz*", mode: 'copy'
    publishDir "${params.data}/samples/${barcode}/deepvariant", pattern: "${barcode}.g.vcf.gz*", mode: 'copy'
    
    input:
    tuple val(barcode), path(called_variants), path(gvcf_records), path(call_variant_outputs)
    val ref
    
    output:
    tuple val(barcode), path("${barcode}.vcf.gz"), path("${barcode}.vcf.gz.tbi"), emit: vcf
    tuple val(barcode), path("${barcode}.g.vcf.gz"), path("${barcode}.g.vcf.gz.tbi"), emit: gvcf
    
    script:
    num_threads = params.deepvariant_threads
    """
    mkdir -p ${params.data}/samples/${barcode}/deepvariant

    ulimit -v unlimited

    echo "Listing files in current directory:"
    ls -al
    echo "============================="
    
    /opt/deepvariant/bin/postprocess_variants \\
        --ref ${ref} \\
        --infile call_variants_output@1.tfrecord.gz \\
        --outfile ${barcode}.vcf.gz \\
        --cpus "${num_threads}" \\
        --small_model_cvo_records make_examples_call_variant_outputs.tfrecord@${num_threads}.gz \\
        --gvcf_outfile ${barcode}.g.vcf.gz \\
        --nonvariant_site_tfrecord_path gvcf.tfrecord@${num_threads}.gz
    
    # Index the output files only if index doesn't exist
    if [[ ! -f ${barcode}.vcf.gz.tbi ]]; then
        tabix -p vcf ${barcode}.vcf.gz
    fi
    
    if [[ ! -f ${barcode}.g.vcf.gz.tbi ]]; then
        tabix -p vcf ${barcode}.g.vcf.gz
    fi
    """
    
    stub:
    """
    touch ${barcode}.vcf.gz
    touch ${barcode}.vcf.gz.tbi
    touch ${barcode}.g.vcf.gz
    touch ${barcode}.g.vcf.gz.tbi
    """
}