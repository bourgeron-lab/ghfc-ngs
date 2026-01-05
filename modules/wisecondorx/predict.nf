process PREDICT {
    tag "${barcode}"
    
    publishDir "${params.data}/samples/${barcode}/svs/wisecondorx", mode: 'copy'
    
    input:
    tuple val(barcode), path(npz), path(reference)
    
    output:
    tuple val(barcode), path("${barcode}_aberrations.bed"), path("${barcode}_chr_statistics.txt"), path("${barcode}_chr_plots.pdf")
    
    script:
    """
    WisecondorX predict ${npz} ${reference} ${barcode} --plot --bed
    """
}
