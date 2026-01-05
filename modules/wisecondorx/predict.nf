process PREDICT {
    tag "${barcode}"
    
    publishDir "${params.data}/samples/${barcode}/svs/wisecondorx", mode: 'copy'
    
    input:
    tuple val(barcode), path(npz), path(reference)
    
    output:
    tuple val(barcode), 
          path("${barcode}_aberrations.bed"), 
          path("${barcode}_statistics.txt"), 
          path("${barcode}_segments.bed"), 
          path("${barcode}_bins.bed"), 
          path("${barcode}.plots/*.png")
    
    script:
    """
    WisecondorX predict ${npz} ${reference} ${barcode} --plot --bed
    """
}
