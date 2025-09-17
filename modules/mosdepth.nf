process MOSDEPTH {

    tag "$barcode"
    
    publishDir "${data_dir}/bedgraphs", mode: 'copy', pattern: "*.bedgraph.gz"

    input:
    tuple val(barcode), path(cram), path(crai)
    val ref
    val bin
    val data_dir

    output:
    tuple val(barcode), path("${barcode}.by${bin}.bedgraph.gz"), emit: bedgraph

    script:
    """
    echo "Starting mosdepth for sample ${barcode} at \$(date)"
    
    mosdepth -t ${task.cpus} -b ${bin} -n -f ${ref} -x -Q 40 -m ${barcode}.by${bin} ${cram}
    
    echo "Mosdepth completed for sample ${barcode} at \$(date)"

    sleep 5

    echo "Renaming output file for sample ${barcode} at \$(date)"

    # Rename the output file to match expected naming
    mv ${barcode}.by${bin}.regions.bed.gz ${barcode}.by${bin}.bedgraph.gz
    
    echo "Renaming completed for sample ${barcode} at \$(date)"
    """

    stub:
    """
    touch ${barcode}.by${bin}.bedgraph.gz
    """
}
