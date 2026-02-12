/*
 * Wombat Workflow
 * Converts annotated BCF files to Parquet and runs PyWombat analysis
 *
 * Supports two paths per wombat config:
 *   - Standard: PYWOMBAT (wombat filter → published TSV)
 *   - External MNV annotation: PYWOMBAT_FILTER → VEP_MNV_ANNOTATION → WOMBAT_MERGE_ANNOTATIONS
 */

import org.yaml.snakeyaml.Yaml

// Include modules
include { BCF2PARQUET } from '../modules/wombat/bcf2parquet'
include { PYWOMBAT } from '../modules/wombat/pywombat'
include { PYWOMBAT_FILTER } from '../modules/wombat/pywombat_filter'
include { VEP_MNV_ANNOTATION } from '../modules/wombat/vep_mnv_annotation'
include { WOMBAT_MERGE_ANNOTATIONS } from '../modules/wombat/wombat_merge_annotations'

workflow WOMBAT {

    take:
    annotated_bcfs     // channel: [fid, bcf, csi]
    family_pedigrees   // channel: [fid, pedigree_file]
    normalized_bcfs    // channel: [fid, bcf, csi]
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
    } else {
        // Parse each wombat config YAML to check for external MNV annotation
        def yaml_parser = new Yaml()

        // Create input channel for PyWombat by combining Parquet files with pedigrees and wombat configs
        // For each Parquet file, create entries for each wombat config
        all_wombat_input = all_parquet_files
            .combine(family_pedigrees, by: 0)  // Combine by fid to get pedigree
            .combine(normalized_bcfs, by: 0)   // Combine by fid to get normalized BCF
            .map { fid, parquet, pedigree, norm_bcf, norm_csi ->
                // Create a list of tuples, one for each wombat config
                params.wombat_config_list.collect { config_file ->
                    def config_name = config_file.replaceAll(/\.ya?ml$/, '')
                    def config_path = file("${params.wombat_config_path}/${config_file}")
                    def config_content = yaml_parser.load(config_path.text)
                    def is_external = config_content?.mnv?.annotate == 'external'
                    tuple(fid, parquet, pedigree, config_path, config_name, params.vep_config_name, norm_bcf, norm_csi, is_external)
                }
            }
            .flatMap()  // Flatten the list of tuples into individual tuples
            .branch {
                external: it[8] == true
                standard: true
            }

        // Standard path: existing PYWOMBAT process (strip the is_external flag)
        pywombat_input = all_wombat_input.standard
            .map { it[0..7] }  // Remove is_external flag

        PYWOMBAT(pywombat_input)

        // External MNV annotation path: 3-step pipeline
        pywombat_filter_input = all_wombat_input.external
            .map { it[0..7] }  // Remove is_external flag

        PYWOMBAT_FILTER(pywombat_filter_input)

        // VEP annotation on the .to_annotate.vcf (extract fid, config_name, vep_config_name, vcf)
        vep_mnv_input = PYWOMBAT_FILTER.out
            .map { fid, wombat_config_name, vep_config_name, filtered_tsv, to_annotate_vcf ->
                tuple(fid, wombat_config_name, vep_config_name, to_annotate_vcf)
            }

        VEP_MNV_ANNOTATION(vep_mnv_input)

        // Merge annotations: combine filtered TSV with VEP annotations
        merge_input = PYWOMBAT_FILTER.out
            .map { fid, wombat_config_name, vep_config_name, filtered_tsv, to_annotate_vcf ->
                tuple(fid, wombat_config_name, vep_config_name, filtered_tsv)
            }
            .join(VEP_MNV_ANNOTATION.out, by: [0, 1, 2])  // Join by fid, wombat_config_name, vep_config_name
            .map { fid, wombat_config_name, vep_config_name, filtered_tsv, vep_annotations ->
                tuple(fid, wombat_config_name, vep_config_name, filtered_tsv, vep_annotations)
            }

        WOMBAT_MERGE_ANNOTATIONS(merge_input)

        // Combine outputs from both paths (both emit [fid, wombat_config_name, tsv])
        wombat_outputs_channel = PYWOMBAT.out.mix(WOMBAT_MERGE_ANNOTATIONS.out)
    }

    emit:
    parquet_files = all_parquet_files
    wombat_outputs = wombat_outputs_channel
}
