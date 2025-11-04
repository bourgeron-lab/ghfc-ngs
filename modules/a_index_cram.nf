process A_INDEX_CRAM {
    
    tag "$barcode:$unit"
    
    publishDir "${params.scratch}/cram", mode: 'copy'
    
    input:
    tuple val(barcode), val(unit), path(cram)
    
    output:
    tuple val(barcode), val(unit), path("${cram}.crai"), emit: crai
    
    script:
    """
    ${params.samtools} index ${cram}
    """
    
    stub:
    """
    touch ${cram}.crai
    """
}
