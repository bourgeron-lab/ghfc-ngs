process CALL_VARIANTS {
    
    tag "$barcode"
    
    input:
    tuple val(barcode), path(example_files)
    tuple val(barcode), path(example_info)
    
    output:
    tuple val(barcode), path("call_variants_output-00000-of-00001.tfrecord.gz"), emit: called_variants
    
    script:
    num_threads = params.deepvariant_threads
    """
    /opt/deepvariant/bin/call_variants \\
        --outfile call_variants_output.tfrecord.gz \\
        --examples make_examples.tfrecord@${num_threads}.gz \\
        --checkpoint /opt/models/wgs
    """
    
    stub:
    """
    touch call_variants_output-00000-of-00001.tfrecord.gz
    """
}
