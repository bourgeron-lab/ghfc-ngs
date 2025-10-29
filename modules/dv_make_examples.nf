process DV_MAKE_EXAMPLES {
    
    tag "$barcode"
    
    input:
    tuple val(barcode), path(cram_file), path(crai_file)
    val ref
    
    output:
    tuple val(barcode), path("make_examples.tfrecord-*-of-*.gz"), emit: examples
    tuple val(barcode), path("make_examples.tfrecord-*-of-*.gz.example_info.json"), emit: examples_info
    tuple val(barcode), path("gvcf.tfrecord-*-of-*.gz"), emit: gvcf_records
    tuple val(barcode), path("make_examples_call_variant_outputs.tfrecord-*-of-*.gz"), emit: call_variant_outputs
    
    script:
    num_threads = params.deepvariant_threads
    num_threads_str = String.format("%05d", num_threads)
    max_task_id = num_threads - 1
    """
    # Run all ${num_threads} tasks in parallel using background processes
    for i in {0..${max_task_id}}; do
        task_id=\$(printf "%05d" \$i)
        (
            /opt/deepvariant/bin/make_examples \\
                --mode calling \\
                --ref ${ref} \\
                --reads ${cram_file} \\
                --examples make_examples.tfrecord@${num_threads}.gz \\
                --checkpoint "/opt/models/wgs" \\
                --call_small_model_examples \\
                --gvcf gvcf.tfrecord@${num_threads}.gz \\
                --small_model_indel_gq_threshold "28" \\
                --small_model_snp_gq_threshold "20" \\
                --small_model_vaf_context_window_size "51" \\
                --nosmall_model_call_multiallelics \\
                --track_ref_reads \\
                --trained_small_model_path "/opt/smallmodels/wgs" \\
                --task \$i
        ) &
    done
    
    # Wait for all background processes to complete
    wait
    
    # Verify all files were created
    for i in {0..${max_task_id}}; do
        task_id=\$(printf "%05d" \$i)
        if [[ ! -f make_examples.tfrecord-\${task_id}-of-${num_threads_str}.gz ]]; then
            echo "ERROR: Task \$i failed to create output file"
            exit 1
        fi
    done
    """
    
    stub:
    num_threads = params.deepvariant_threads
    num_threads_str = String.format("%05d", num_threads)
    max_task_id = num_threads - 1
    """
    for i in {0..${max_task_id}}; do
        task_id=\$(printf "%05d" \$i)
        touch make_examples.tfrecord-\${task_id}-of-${num_threads_str}.gz
        touch gvcf.tfrecord-\${task_id}-of-${num_threads_str}.gz
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