/*
 * Wombat Workflow
 * Converts annotated BCF files to TSV and runs PyWombat analysis
 */

// Include modules
include { BCF2TSV } from '../modules/wombat/bcf2tsv'
include { PYWOMBAT } from '../modules/wombat/pywombat'
include { MERGE_WOMBAT } from '../modules/wombat/merge_wombat'

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
    
    // Group PyWombat outputs by wombat_config_name and output_name for merging
    // First, we need to read the wombat config files to get the output names
    // For now, we'll collect all outputs and merge by pattern matching
    
    // Collect all wombat output files grouped by config name
    wombat_outputs = PYWOMBAT.out
        .groupTuple(by: 1)  // Group by wombat_config_name
    
    if (!params.wombat_config_list || params.wombat_config_list.isEmpty()) {
        merged_outputs_channel = Channel.empty()
    } else {
        // Collect all PyWombat outputs for merging
        // Group files by config name and then create separate channels for each output type
        // We'll merge all files matching specific patterns
        merge_jobs = wombat_outputs_channel
            .map { fid, wombat_config_name, files ->
                // files is a list of output files from one family for one config
                // Return tuples for each file with its config name
                files.collect { file ->
                    tuple(wombat_config_name, file)
                }
            }
            .flatMap()  // Flatten to individual tuples
            .groupTuple(by: 0)  // Group all files by wombat_config_name
            .flatMap { wombat_config_name, file_list ->
                // Extract unique output names from file names
                // File pattern: {FID}.rare.{vep_config_name}.{wombat_config_name}.{output_name}.tsv.gz
                def output_names = file_list.collect { file ->
                    // Extract output_name from filename
                    def parts = file.name.tokenize('.')
                    // Find the index of the wombat_config_name and get the next element
                    def idx = parts.indexOf(wombat_config_name)
                    if (idx >= 0 && idx < parts.size() - 2) {
                        parts[idx + 1]  // output_name is after wombat_config_name
                    } else {
                        null
                    }
                }.findAll { it != null }.unique()
                
                // Create a merge job for each unique output name
                output_names.collect { output_name ->
                    def matching_files = file_list.findAll { file ->
                        file.name.contains(".${wombat_config_name}.${output_name}.")
                    }
                    tuple(params.cohort_name, matching_files, params.vep_config_name, wombat_config_name, output_name)
                }
            }
        
        // Run merge for each output name
        MERGE_WOMBAT(
            merge_jobs.map { it[0] },  // cohort_name
            merge_jobs.map { it[1] },  // files
            merge_jobs.map { it[2] },  // vep_config_name
            merge_jobs.map { it[3] },  // wombat_config_name
            merge_jobs.map { it[4] }   // output_name
        )
        
        merged_outputs_channel = MERGE_WOMBAT.out.cohort_wombat_tsv
    }
    
    emit:
    tsv_files = all_tsv_files
    wombat_outputs = wombat_outputs_channel
    merged_outputs = merged_outputs_channel
}
