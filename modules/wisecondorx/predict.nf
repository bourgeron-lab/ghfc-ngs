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
          path("${barcode}.plots/genome_wide.png"),
          path("${barcode}.plots/chr1.png"),
          path("${barcode}.plots/chr2.png"),
          path("${barcode}.plots/chr3.png"),
          path("${barcode}.plots/chr4.png"),
          path("${barcode}.plots/chr5.png"),
          path("${barcode}.plots/chr6.png"),
          path("${barcode}.plots/chr7.png"),
          path("${barcode}.plots/chr8.png"),
          path("${barcode}.plots/chr9.png"),
          path("${barcode}.plots/chr10.png"),
          path("${barcode}.plots/chr11.png"),
          path("${barcode}.plots/chr12.png"),
          path("${barcode}.plots/chr13.png"),
          path("${barcode}.plots/chr14.png"),
          path("${barcode}.plots/chr15.png"),
          path("${barcode}.plots/chr16.png"),
          path("${barcode}.plots/chr17.png"),
          path("${barcode}.plots/chr18.png"),
          path("${barcode}.plots/chr19.png"),
          path("${barcode}.plots/chr20.png"),
          path("${barcode}.plots/chr21.png"),
          path("${barcode}.plots/chr22.png"),
          path("${barcode}.plots/chrX.png"),
          path("${barcode}.plots/chrY.png")
    
    script:
    """
    WisecondorX predict ${npz} ${reference} ${barcode} --plot --bed
    """
}
