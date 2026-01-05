/*
 * WisecondorX Workflow
 * Converts CRAM files to NPZ format and runs WisecondorX predict for CNV analysis
 */

// Include modules
include { NPZ_CONVERT } from '../modules/wisecondorx/npz_convert'
include { PREDICT } from '../modules/wisecondorx/predict'

workflow WISECONDORX {

    take:
    cram_files         // channel: [barcode, cram, crai]
    existing_npz_files // channel: [barcode, npz]
    need_predict       // list of barcodes that need predict

    main:
    
    // Add reference genome to each input for NPZ conversion
    cram_with_ref = cram_files
        .map { barcode, cram, crai -> 
            tuple(barcode, cram, crai, file(params.ref))
        }
    
    // Run NPZ conversion for samples that need it
    NPZ_CONVERT(cram_with_ref)
    
    // Mix existing NPZ files with newly converted ones
    all_npz_files = existing_npz_files
        .mix(NPZ_CONVERT.out)
    
    // Filter NPZ files for samples that need predict
    npz_for_predict = all_npz_files
        .filter { barcode, npz ->
            barcode in need_predict
        }
        .map { barcode, npz ->
            tuple(barcode, npz, file(params.wisecondorx_reference))
        }
    
    // Run WisecondorX predict
    PREDICT(npz_for_predict)
    
    emit:
    npz_files = NPZ_CONVERT.out
    predict_results = PREDICT.out
}
