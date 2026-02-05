/*
 * Wombat Workflow
 * Converts annotated BCF files to Parquet and runs PyWombat analysis
 */

// Include modules
include { BCF2PARQUET } from '../modules/wombat/bcf2parquet'
include { PYWOMBAT } from '../modules/wombat/pywombat'

workflow WOMBAT {

    take:
    annotated_bcfs     // channel: [fid, bcf, csi]
    family_pedigrees   // channel: [fid, pedigree_file]
    need_bcf2parquet   // map: [fid: boolean] - whether BCF2PARQUET is needed for each family

    main:

    // Filter families that need BCF2PARQUET conversion
    bcf2parquet_input = annotated_bcfs
        .filter { fid, bcf, csi ->
            need_bcf2parquet[fid] == true
        }
        .map { fid, bcf, csi ->
            tuple(fid, bcf, csi, params.vep_config_name)
        }

    // Run BCF2PARQUET conversion
    BCF2PARQUET(bcf2parquet_input)

    // Create channel for existing Parquet files
    existing_parquet = annotated_bcfs
        .filter { fid, bcf, csi ->
            need_bcf2parquet[fid] == false
        }
        .map { fid, bcf, csi ->
            def parquet_path = file("${params.data}/families/${fid}/wombat/${fid}.rare.${params.vep_config_name}.annotated.parquet")
            tuple(fid, parquet_path)
        }

    // Combine all available Parquet files (existing + newly created)
    all_parquet_files = existing_parquet.mix(BCF2PARQUET.out)
    
    // Check if wombat_config_list is defined and not empty
    if (!params.wombat_config_list || params.wombat_config_list.isEmpty()) {
        log.warn "WARNING: wombat_config_list is empty. No PyWombat analysis will be performed."
        wombat_outputs_channel = Channel.empty()
        merged_outputs_channel = Channel.empty()
    } else {
        // Create input channel for PyWombat by combining Parquet files with pedigrees and wombat configs
        // For each Parquet file, create entries for each wombat config
        pywombat_input = all_parquet_files
            .combine(family_pedigrees, by: 0)  // Combine by fid to get pedigree
            .map { fid, parquet, pedigree ->
                // Create a list of tuples, one for each wombat config
                params.wombat_config_list.collect { config_file ->
                    def config_name = config_file.replaceAll(/\.ya?ml$/, '')
                    def config_path = file("${params.wombat_config_path}/${config_file}")
                    tuple(fid, parquet, pedigree, config_path, config_name, params.vep_config_name)
                }
            }
            .flatMap()  // Flatten the list of tuples into individual tuples

        // Run PyWombat
        PYWOMBAT(pywombat_input)
        wombat_outputs_channel = PYWOMBAT.out
    }

    emit:
    parquet_files = all_parquet_files
    wombat_outputs = wombat_outputs_channel
}
