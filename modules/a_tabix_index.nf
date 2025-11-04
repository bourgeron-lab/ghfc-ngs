process A_TABIX_INDEX {

    tag "$barcode"
    
    publishDir "${data_dir}/samples/${barcode}/sequences", mode: 'copy', pattern: "*.tbi"

    input:
    tuple val(barcode), path(bedgraph)
    val data_dir

    output:
    tuple val(barcode), path(bedgraph), path("${bedgraph}.tbi"), emit: indexed_bedgraph

    script:
    """
    echo "Starting tabix indexing for sample ${barcode} at \$(date)"
    
    tabix -p bed ${bedgraph}
    
    echo "Tabix indexing completed for sample ${barcode} at \$(date)"
    """

    stub:
    """
    touch ${bedgraph}.tbi
    """
}
