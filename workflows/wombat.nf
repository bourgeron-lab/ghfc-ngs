/*
 * Wombat Workflow
 * Converts annotated BCF files to TSV and runs PyWombat analysis
 */

// Include modules
include { BCF2TSV } from '../modules/wombat/bcf2tsv'
include { PYWOMBAT } from '../modules/wombat/pywombat'

workflow WOMBAT {

    take:
    annotated_bcfs     // channel: [fid, bcf, csi]
    family_pedigrees   // channel: [fid, pedigree_file]
    need_bcf2tsv       // map: [fid: boolean] - whether BCF2TSV is needed for each family

    main:

    // Filter families that need BCF2TSV conversion
    bcf2tsv_input = annotated_bcfs
        .filter { fid, bcf, csi -> 
            need_bcf2tsv[fid] == true 
        }
        .map { fid, bcf, csi -> 
            tuple(fid, bcf, csi, params.vep_config_name)
        }

    // Run BCF2TSV conversion
    BCF2TSV(bcf2tsv_input)
    
    // Create channel for existing TSV files
    existing_tsv = annotated_bcfs
        .filter { fid, bcf, csi -> 
            need_bcf2tsv[fid] == false 
        }
        .map { fid, bcf, csi ->
            def tsv_path = file("${params.data}/families/${fid}/wombat/${fid}.rare.${params.vep_config_name}.annotated.tsv.gz")
            tuple(fid, tsv_path)
        }
    
    // Combine all available TSV files (existing + newly created)
    all_tsv_files = existing_tsv.mix(BCF2TSV.out)
    
    // Check if wombat_config_list is defined and not empty
    if (!params.wombat_config_list || params.wombat_config_list.isEmpty()) {
        log.warn "WARNING: wombat_config_list is empty. No PyWombat analysis will be performed."
        wombat_outputs_channel = Channel.empty()
        merged_outputs_channel = Channel.empty()
    } else {
        // Create input channel for PyWombat by combining TSV files with pedigrees and wombat configs
        // For each TSV file, create entries for each wombat config
        pywombat_input = all_tsv_files
            .combine(family_pedigrees, by: 0)  // Combine by fid to get pedigree
            .map { fid, tsv, pedigree ->
                // Create a list of tuples, one for each wombat config
                params.wombat_config_list.collect { config_file ->
                    def config_name = config_file.replaceAll(/\.ya?ml$/, '')
                    def config_path = file("${params.wombat_config_path}/${config_file}")
                    tuple(fid, tsv, pedigree, config_path, config_name, params.vep_config_name)
                }
            }
            .flatMap()  // Flatten the list of tuples into individual tuples
        
        // Run PyWombat
        PYWOMBAT(pywombat_input)
        wombat_outputs_channel = PYWOMBAT.out
    }
    
    emit:
    tsv_files = all_tsv_files
    wombat_outputs = wombat_outputs_channel
}
