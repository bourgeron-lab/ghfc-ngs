process PREDICT {
    tag "${barcode}"
    
    publishDir {
        def (s1, s2) = Sharding.getShards(barcode)
        "${params.data}/samples/${s1}/${s2}/${barcode}/svs/wisecondorx"
    }, mode: 'copy'
    
    input:
    tuple val(barcode), path(npz), path(reference)
    
    output:
    tuple val(barcode),
          path("${barcode}_aberrations.bed"),
          path("${barcode}_statistics.txt"),
          path("${barcode}_segments.bed"),
          path("${barcode}_bins.bed")

    script:
    """
    wisecondorx predict ${npz} ${reference} ${barcode} --bed ${params.wisecondorx_predict_args}
    """
}
