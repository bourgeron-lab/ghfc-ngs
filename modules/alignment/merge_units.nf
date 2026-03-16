process MERGE_UNITS {
    
    tag "$barcode"
    
    publishDir {
        def (s1, s2) = Sharding.getShards(barcode)
        "${data_dir}/samples/${s1}/${s2}/${barcode}/sequences"
    }, mode: 'copy'
    
    input:
    tuple val(barcode), path(crams), path(crais)
    val samtools_path
    val data_dir
    val ref_name
    
    output:
    tuple val(barcode), path("${barcode}.${ref_name}.cram"), path("${barcode}.${ref_name}.cram.crai"), emit: cram
    
    script:
    if (crams.size() > 1) {
        """
        ${samtools_path} merge ${barcode}.${ref_name}.cram ${crams.join(' ')}
        ${samtools_path} index ${barcode}.${ref_name}.cram
        """
    } else {
        """
        cp ${crams[0]} ${barcode}.${ref_name}.cram
        cp ${crais[0]} ${barcode}.${ref_name}.cram.crai
        """
    }
    
    stub:
    """
    touch ${barcode}.${ref_name}.cram
    touch ${barcode}.${ref_name}.cram.crai
    """
}
