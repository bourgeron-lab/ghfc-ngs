/*
 * WisecondorX Workflow
 * Converts CRAM files to NPZ format for WisecondorX CNV analysis
 */

// Include modules
include { NPZ_CONVERT } from '../modules/wisecondorx/npz_convert'

workflow WISECONDORX {

    take:
    cram_files    // channel: [barcode, cram, crai]

    main:
    
    // Add reference genome to each input
    cram_with_ref = cram_files
        .map { barcode, cram, crai -> 
            tuple(barcode, cram, crai, file(params.ref))
        }
    
    // Run NPZ conversion
    NPZ_CONVERT(cram_with_ref)
    
    emit:
    npz_files = NPZ_CONVERT.out
}
