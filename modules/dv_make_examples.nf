process DV_MAKE_EXAMPLES {
    
    tag "$barcode"
    
    input:
    tuple val(barcode), path(cram_file), path(crai_file)
    val ref
    
    output:
    tuple val(barcode), path("make_examples.tfrecord-*-of-00096.gz"), emit: examples
    tuple val(barcode), path("make_examples.tfrecord-*-of-00096.gz.example_info.json"), emit: examples_info
    tuple val(barcode), path("gvcf.tfrecord-*-of-00096.gz"), emit: gvcf_records
    tuple val(barcode), path("make_examples_call_variant_outputs.tfrecord-*-of-00096.gz"), emit: call_variant_outputs
    
    script:
    """
    # Run all 96 tasks in parallel using background processes
    for i in {0..95}; do
        task_id=\$(printf "%05d" \$i)
        (
            /opt/deepvariant/bin/make_examples \\
                --mode calling \\
                --ref ${ref} \\
                --reads ${cram_file} \\
                --examples make_examples.tfrecord@96.gz \\
                --checkpoint "/opt/models/wgs" \\
                --call_small_model_examples \\
                --gvcf gvcf.tfrecord@96.gz \\
                --small_model_indel_gq_threshold "28" \\
                --small_model_snp_gq_threshold "20" \\
                --small_model_vaf_context_window_size "51" \\
                --track_ref_reads \\
                --trained_small_model_path "/opt/smallmodels/wgs" \\
                --task \$i
        ) &
    done
    
    # Wait for all background processes to complete
    wait
    
    # Verify all files were created
    for i in {0..95}; do
        task_id=\$(printf "%05d" \$i)
        if [[ ! -f make_examples.tfrecord-\${task_id}-of-00096.gz ]]; then
            echo "ERROR: Task \$i failed to create output file"
            exit 1
        fi
    done
    """
    
    stub:
    """
    for i in {0..95}; do
        task_id=\$(printf "%05d" \$i)
        touch make_examples.tfrecord-\${task_id}-of-00096.gz
        touch gvcf.tfrecord-\${task_id}-of-00096.gz
    done
    """
}

// previous command
// if [[ ! -f make_examples.tfrecord-\${task_id}-of-00096.gz ]]; then
//     (
//         make_examples \\
//             --mode calling \\
//             --ref ${ref} \\
//             --reads ${cram_file} \\
//             --examples make_examples.tfrecord@96.gz \\
//             --gvcf gvcf.tfrecord@96.gz \\
//             --channel_list=read_base,base_quality,mapping_quality,strand,read_supports_variant,base_differs_from_ref,insert_size \\
//             --task \$i
//     ) &
// fi