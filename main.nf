#!/usr/bin/env nextflow

/*
========================================================================================
    GHFC WGS Family-based Variant Calling Pipeline
========================================================================================
    Github : https://github.com/your-repo/ghfc-ngs
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// Include workflows
include { ALIGNMENT } from './workflows/alignment'
include { DEEPVARIANT_SAMPLE } from './workflows/deepvariant_sample'
include { DEEPVARIANT_FAMILY } from './workflows/deepvariant_family'
include { ANNOTATION } from './workflows/annotation'
include { SNVS_COHORT } from './workflows/snvs_cohort'
include { WISECONDORX } from './workflows/wisecondorx'
include { WOMBAT } from './workflows/wombat'
include { EXTRACTOR } from './workflows/extractor'

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
if (!params.data) {
    exit 1, "ERROR: --data parameter is required"
}

if (!params.steps || params.steps.isEmpty()) {
    exit 1, "ERROR: --steps parameter is required. Available steps: alignment, deepvariant_sample, deepvariant_family, snvs_freq_annot, snvs_freq_filter, bedgraphs, vep_annotation"
}

// Validate steps
def valid_steps = ['alignment', 'deepvariant_sample', 'deepvariant_family', 'annotation', 'snvs_cohort', 'wisecondorx', 'wombat', 'extractor']
def invalid_steps = params.steps - valid_steps
if (invalid_steps) {
    exit 1, "ERROR: Invalid steps specified: ${invalid_steps.join(', ')}. Valid steps are: ${valid_steps.join(', ')}"
}

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    
    // Read and validate pedigree file
    def pedigree_file = params.pedigree ?: "${params.data}/pedigree.tsv"
    
    if (!new File(pedigree_file).exists()) {
        exit 1, "ERROR: Pedigree file not found: ${pedigree_file}"
    }
    
    log.info """
    ========================================================================================
                            GHFC WGS Family-based Pipeline
    ========================================================================================
    Pedigree file    : ${pedigree_file}
    Data directory   : ${params.data}
    Steps to run     : ${params.steps.join(', ')}
    Reference        : ${params.ref_name}
    Work directory   : ${workflow.workDir}
    VEP Config Name  : ${params.vep_config_name} (${params.vep_config})
    ========================================================================================
    """
    
    // Parse pedigree file
    def pedigree_data = parsePedigreeFile(pedigree_file)
    def families = pedigree_data.families
    def individuals = pedigree_data.individuals
    def family_members = pedigree_data.family_members
    
    log.info "Found ${families.size()} families with ${individuals.size()} individuals total"
    
    // Check existing files and determine what needs to be done
    def analysis_plan = createAnalysisPlan(families, individuals, family_members)
    
    // Display analysis summary
    displayAnalysisSummary(analysis_plan)
    
    // Validate that required steps are available
    validateStepsAvailability(analysis_plan)
    
    // Create channels for different steps
    def channels = createChannels(analysis_plan)
    
    // Run alignment if needed and allowed
    if (analysis_plan.alignment.needed.size() > 0 && 'alignment' in params.steps) {
        log.info "Running alignment for ${analysis_plan.alignment.needed.size()} individuals..."
        
        ALIGNMENT(
            channels.fastq_files,
            channels.cram_37_files,
            channels.cram_38_files
        )
        
        aligned_crams = ALIGNMENT.out.crams
    } else {
        aligned_crams = Channel.empty()
    }
    
    // Collect all available CRAM files (existing + newly aligned)
    all_available_crams = channels.existing_crams.mix(aligned_crams ?: Channel.empty())
    
    // Run DeepVariant sample workflow if needed and allowed
    if (analysis_plan.deepvariant_sample.needed.size() > 0 && 'deepvariant_sample' in params.steps) {
        log.info "Running DeepVariant sample workflow for ${analysis_plan.deepvariant_sample.needed.size()} individuals..."
        
        // Filter CRAM files for individuals that need DeepVariant
        deepvariant_crams = all_available_crams
            .filter { barcode, cram, crai -> 
                barcode in analysis_plan.deepvariant_sample.needed 
            }
        
        DEEPVARIANT_SAMPLE(deepvariant_crams)
        
        deepvariant_gvcfs = DEEPVARIANT_SAMPLE.out.gvcf
    } else {
        deepvariant_gvcfs = Channel.empty()
    }
    
    // Collect all available gVCF files (existing + newly created)
    all_available_gvcfs = channels.existing_gvcfs.mix(deepvariant_gvcfs ?: Channel.empty())
    
    // Run family calling, normalization, and pedigree extraction if needed and allowed
    family_vcfs_output = Channel.empty()
    normalized_bcfs_output = Channel.empty()
    family_pedigrees_output = Channel.empty()
    if (analysis_plan.deepvariant_family.needed.size() > 0 && 'deepvariant_family' in params.steps) {
        log.info "Running family calling, normalization, and pedigree extraction for ${analysis_plan.deepvariant_family.needed.size()} families..."
        
        // Group gVCF files by family
        family_gvcfs = all_available_gvcfs
            .map { barcode, gvcf, tbi -> 
                def fid = family_members[barcode]
                [fid, barcode, gvcf, tbi]
            }
            .filter { fid, barcode, gvcf, tbi -> 
                fid in analysis_plan.deepvariant_family.needed 
            }
            .groupTuple(by: 0)
            .map { fid, barcodes, gvcfs, tbis ->
                [fid, barcodes, gvcfs, tbis]
            }
        
        DEEPVARIANT_FAMILY(family_gvcfs, pedigree_file)
        family_vcfs_output = DEEPVARIANT_FAMILY.out.family_vcfs
        normalized_bcfs_output = DEEPVARIANT_FAMILY.out.normalized_bcfs
        family_pedigrees_output = DEEPVARIANT_FAMILY.out.family_pedigrees
    }
    
    // Run annotation (gnomAD annotation, filtering, VEP) if needed and allowed
    annotation_common_bcfs_output = Channel.empty()
    if (analysis_plan.annotation.needed.size() > 0 && 'annotation' in params.steps) {
        log.info "Running annotation for ${analysis_plan.annotation.needed.size()} families..."
        
        // Get all available normalized family BCFs (existing + newly created)
        all_available_normalized_bcfs = channels.existing_normalized_bcfs.mix(normalized_bcfs_output ?: Channel.empty())
        
        // Filter for families that need annotation
        annotation_bcfs = all_available_normalized_bcfs
            .filter { fid, bcf, csi -> 
                fid in analysis_plan.annotation.needed 
            }
        
        // Get all available family pedigrees (existing + newly created)
        all_available_family_pedigrees = channels.existing_family_pedigrees.mix(family_pedigrees_output ?: Channel.empty())
        
        // Filter for families that need annotation
        annotation_pedigrees = all_available_family_pedigrees
            .filter { fid, pedigree -> 
                fid in analysis_plan.annotation.needed 
            }
        
        ANNOTATION(annotation_bcfs, annotation_pedigrees)
        annotation_common_bcfs_output = ANNOTATION.out.filtered_common_bcfs
        annotation_annotated_bcfs_output = ANNOTATION.out.fully_annotated_bcfs
    } else {
        annotation_annotated_bcfs_output = Channel.empty()
    }
    
    // Run Wombat analysis if needed and allowed
    wombat_output = Channel.empty()
    if (analysis_plan.wombat.needed.size() > 0 && 'wombat' in params.steps) {
        log.info "Running Wombat analysis for ${analysis_plan.wombat.needed.size()} families..."
        
        // Get all available annotated BCFs (existing + newly created)
        all_available_annotated_bcfs = channels.existing_annotated_bcfs.mix(annotation_annotated_bcfs_output ?: Channel.empty())
        
        // Filter for families that need Wombat
        wombat_bcfs = all_available_annotated_bcfs
            .filter { fid, bcf, csi -> 
                fid in analysis_plan.wombat.needed 
            }
        
        // Get all available family pedigrees (existing + newly created)
        all_available_family_pedigrees_wombat = channels.existing_family_pedigrees.mix(family_pedigrees_output ?: Channel.empty())
        
        // Filter for families that need Wombat
        wombat_pedigrees = all_available_family_pedigrees_wombat
            .filter { fid, pedigree -> 
                fid in analysis_plan.wombat.needed 
            }
        
        WOMBAT(wombat_bcfs, wombat_pedigrees, analysis_plan.wombat.need_bcf2parquet)
        wombat_output = WOMBAT.out.wombat_outputs
    }
    
    // Run cohort common variants merge if needed and allowed
    if (analysis_plan.snvs_cohort.needed.size() > 0 && 'snvs_cohort' in params.steps) {
        def merge_tasks = []
        if (analysis_plan.snvs_cohort.need_bcf_merge) merge_tasks.add("BCF merge")
        def wombat_merge_count = analysis_plan.snvs_cohort.need_wombat_merges?.count { k, v -> v == true } ?: 0
        if (wombat_merge_count > 0) merge_tasks.add("Wombat merge (${wombat_merge_count} configs)")
        log.info "Running cohort: ${merge_tasks.join(' and ')}..."
        
        // Get all available common filtered BCFs (existing + newly created)
        if (analysis_plan.annotation.needed.size() > 0 && 'annotation' in params.steps) {
            // Mix existing BCFs with newly created ones
            all_available_common_bcfs = channels.existing_common_filtered_bcfs.mix(annotation_common_bcfs_output)
        } else {
            // Use only existing BCFs if no new ones were created
            all_available_common_bcfs = channels.existing_common_filtered_bcfs
        }
        
        // Get all available Wombat files (existing + newly created)
        all_available_wombat_files = channels.existing_wombat_files.mix(wombat_output ?: Channel.empty())
        
        SNVS_COHORT(all_available_common_bcfs, all_available_wombat_files,
                    analysis_plan.snvs_cohort.need_bcf_merge, 
                    analysis_plan.snvs_cohort.need_wombat_merges)
    }
    
    // Run WisecondorX predict if needed and allowed
    if ((analysis_plan.wisecondorx.needed.size() > 0 || 
         analysis_plan.wisecondorx.need_family_merge.any { k, v -> v == true } || 
         analysis_plan.wisecondorx.need_family_annotate.any { k, v -> v == true } || 
         analysis_plan.wisecondorx.need_cohort_merge) && 'wisecondorx' in params.steps) {
        
        def tasks = []
        if (analysis_plan.wisecondorx.needed.size() > 0) {
            tasks.add("predict for ${analysis_plan.wisecondorx.needed.size()} individuals")
        }
        def family_merge_count = analysis_plan.wisecondorx.need_family_merge.count { k, v -> v == true } ?: 0
        if (family_merge_count > 0) {
            tasks.add("family merge for ${family_merge_count} families")
        }
        def family_annotate_count = analysis_plan.wisecondorx.need_family_annotate.count { k, v -> v == true } ?: 0
        if (family_annotate_count > 0) {
            tasks.add("family annotation for ${family_annotate_count} families")
        }
        if (analysis_plan.wisecondorx.need_cohort_merge) {
            tasks.add("cohort merge")
        }
        log.info "Running WisecondorX: ${tasks.join(', ')}..."
        if (analysis_plan.wisecondorx.needed.size() > 0) {
            log.info "  - NPZ conversion needed for ${analysis_plan.wisecondorx.need_npz.size()} individuals"
            log.info "  - Predict needed for ${analysis_plan.wisecondorx.need_predict.size()} individuals"
        }
        
        // Filter CRAM files for individuals that need NPZ conversion
        wisecondorx_crams = all_available_crams
            .filter { barcode, cram, crai -> 
                barcode in analysis_plan.wisecondorx.need_npz
            }
        
        // Run WisecondorX workflow with existing NPZ files and CRAMs that need conversion
        WISECONDORX(
            wisecondorx_crams,
            channels.existing_npz_files,
            analysis_plan.wisecondorx.need_predict,
            pedigree_data.family_members,
            analysis_plan.wisecondorx.need_family_merge,
            analysis_plan.wisecondorx.need_family_annotate,
            analysis_plan.wisecondorx.need_cohort_merge
        )
    }
    
    // Run Extractor if TSV files are provided and step is allowed
    if (analysis_plan.extractor.tsv_count > 0 && 'extractor' in params.steps) {
        log.info "Running Extractor for ${analysis_plan.extractor.tsv_count} TSV files on ${analysis_plan.extractor.families.size()} families / ${analysis_plan.extractor.samples.size()} samples..."
        
        // Create channel from TSV list
        extractor_tsvs = Channel.fromList(params.extractor_tsvs_list)
            .map { tsv_path ->
                def original_filename = file(tsv_path).baseName.replaceAll(/\.(?i)(tsv|txt|csv)$/, '')
                tuple(original_filename, tsv_path)
            }
        
        // Create channel for normalized BCFs (for PROCESS_FAM_BCF)
        norm_bcfs_for_extractor = channels.existing_normalized_bcfs
            .map { fid, bcf, csi ->
                tuple(fid, bcf, csi)
            }
        
        // Create channel for wombat bcf2parquet outputs (for PROCESS_FAM_TSV)
        wombat_parquets_for_extractor = Channel
            .fromPath("${params.data}/families/*/wombat/*.rare.${params.vep_config_name}.annotated.parquet")
            .map { parquet ->
                def fid = parquet.parent.parent.name
                tuple(fid, parquet)
            }
        
        // Create channel for gVCFs (for PROCESS_IND_GVCF)
        gvcfs_for_extractor = channels.existing_gvcfs
            .map { barcode, gvcf, tbi ->
                tuple(barcode, gvcf, tbi)
            }
        
        EXTRACTOR(
            extractor_tsvs,
            pedigree_file,
            params.liftover_chain,
            norm_bcfs_for_extractor,
            wombat_parquets_for_extractor,
            gvcfs_for_extractor
        )
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

def parsePedigreeFile(pedigree_file) {
    def families = [] as Set
    def individuals = [] as Set
    def family_members = [:]
    
    file(pedigree_file).readLines().eachWithIndex { line, index ->
        // Skip comments and empty lines
        if (line.startsWith('#') || line.trim().isEmpty()) return
        
        def cols = line.split('\t')
        
        // Skip header row if first column is "FID"
        if (index == 0 && cols[0] == 'FID') {
            log.info "Skipping pedigree header row"
            return
        }
        
        if (cols.size() < 6) {
            exit 1, "ERROR: Pedigree file must have 6 columns (FID, barcode, father, mother, sex, phenotype). Found ${cols.size()} columns in line: ${line}"
        }
        
        def fid = cols[0]
        def barcode = cols[1]
        
        families.add(fid)
        individuals.add(barcode)
        family_members[barcode] = fid
    }
    
    return [
        families: families,
        individuals: individuals,
        family_members: family_members
    ]
}

def createAnalysisPlan(families, individuals, family_members) {
    def plan = [
        alignment: [needed: [], existing: []],
        deepvariant_sample: [needed: [], existing: []],
        deepvariant_family: [needed: [], existing: []],
        annotation: [needed: [], existing: []],
        snvs_cohort: [needed: [], existing: []],
        wisecondorx: [needed: [], existing: [], need_npz: [], need_predict: [], need_family_merge: [:], need_family_annotate: [:], need_cohort_merge: false],
        wombat: [needed: [], existing: [], need_bcf2parquet: [:]],
        extractor: [tsv_count: 0, families: [] as Set, samples: [] as Set]
    ]
    
    // Check extractor TSV list
    if (params.extractor_tsvs_list && !params.extractor_tsvs_list.isEmpty()) {
        plan.extractor.tsv_count = params.extractor_tsvs_list.size()
        // All families and individuals with data will be processed
        plan.extractor.families = families
        plan.extractor.samples = individuals
    }
    
    // Check existing family outputs (normalized BCF and pedigree - all part of deepvariant_family)
    // Note: {FID}.vcf.gz is intermediate (GLnexus output) and not published
    families.each { fid ->
        def norm_bcf_path = "${params.data}/families/${fid}/vcfs/${fid}.norm.bcf"
        def norm_csi_path = "${params.data}/families/${fid}/vcfs/${fid}.norm.bcf.csi"
        def pedigree_path = "${params.data}/families/${fid}/${fid}.pedigree.tsv"
        
        if (new File(norm_bcf_path).exists() && new File(norm_csi_path).exists() &&
            new File(pedigree_path).exists()) {
            plan.deepvariant_family.existing.add(fid)
        } else {
            plan.deepvariant_family.needed.add(fid)
        }
    }
    
    // Check existing annotation outputs (rare/common VCF.gz/BCFs, common_gt BCF, VEP VCF.gz, and final annotated BCF - all part of annotation)
    // Note: {FID}.gnomad.bcf is intermediate and not published
    families.each { fid ->
        def rare_vcf_path = "${params.data}/families/${fid}/vcfs/${fid}.rare.vcf.gz"
        def rare_tbi_path = "${params.data}/families/${fid}/vcfs/${fid}.rare.vcf.gz.tbi"
        def common_bcf_path = "${params.data}/families/${fid}/vcfs/${fid}.common.bcf"
        def common_csi_path = "${params.data}/families/${fid}/vcfs/${fid}.common.bcf.csi"
        def common_gt_bcf_path = "${params.data}/families/${fid}/vcfs/${fid}.common_gt.bcf"
        def common_gt_csi_path = "${params.data}/families/${fid}/vcfs/${fid}.common_gt.bcf.csi"
        def vep_vcf_path = "${params.data}/families/${fid}/vcfs/${fid}.rare.${params.vep_config_name}.vcf.gz"
        def vep_tbi_path = "${params.data}/families/${fid}/vcfs/${fid}.rare.${params.vep_config_name}.vcf.gz.tbi"
        def annotated_bcf_path = "${params.data}/families/${fid}/vcfs/${fid}.rare.${params.vep_config_name}.annotated.bcf"
        def annotated_csi_path = "${params.data}/families/${fid}/vcfs/${fid}.rare.${params.vep_config_name}.annotated.bcf.csi"
        
        if (new File(rare_vcf_path).exists() && new File(rare_tbi_path).exists() &&
            new File(common_bcf_path).exists() && new File(common_csi_path).exists() &&
            new File(common_gt_bcf_path).exists() && new File(common_gt_csi_path).exists() &&
            new File(vep_vcf_path).exists() && new File(vep_tbi_path).exists() &&
            new File(annotated_bcf_path).exists() && new File(annotated_csi_path).exists()) {
            plan.annotation.existing.add(fid)
        } else {
            // Need annotation if normalized BCFs exist or will be created
            if (fid in plan.deepvariant_family.existing || fid in plan.deepvariant_family.needed) {
                plan.annotation.needed.add(fid)
            }
        }
    }

    // Check existing cohort files (common BCF and Wombat TSVs)
    def cohort_bcf_path = "${params.data}/cohorts/${params.cohort_name}/vcfs/${params.cohort_name}.common_gt.bcf"
    def cohort_csi_path = "${params.data}/cohorts/${params.cohort_name}/vcfs/${params.cohort_name}.common_gt.bcf.csi"
    
    def cohort_bcf_exists = new File(cohort_bcf_path).exists() && new File(cohort_csi_path).exists()
    
    // Check wombat cohort files for each config
    plan.snvs_cohort.need_wombat_merges = [:]
    if (params.wombat_config_list && !params.wombat_config_list.isEmpty()) {
        params.wombat_config_list.each { config_file ->
            def config_name = config_file.replaceAll(/\.ya?ml$/, '')
            def cohort_wombat_path = "${params.data}/cohorts/${params.cohort_name}/wombat/${params.cohort_name}.rare.${params.vep_config_name}.annotated.${config_name}.results.tsv"
            plan.snvs_cohort.need_wombat_merges[config_name] = !new File(cohort_wombat_path).exists()
        }
    }
    
    if (cohort_bcf_exists && !plan.snvs_cohort.need_wombat_merges.any { k, v -> v == true }) {
        plan.snvs_cohort.existing.add('cohort')  // Single cohort entry
    } else {
        // Need cohort merge if any families have common filtered BCFs or will create them
        def families_with_common_bcfs = plan.annotation.existing + plan.annotation.needed
        if (families_with_common_bcfs.size() > 0) {
            plan.snvs_cohort.needed.add('cohort')  // Single cohort entry
            // Track what needs to be generated
            plan.snvs_cohort.need_bcf_merge = !cohort_bcf_exists
        }
    }
    
    // Check existing individual gVCF files and VAF bedgraphs (both outputs of deepvariant_sample)
    individuals.each { barcode ->
        def gvcf_path = "${params.data}/samples/${barcode}/deepvariant/${barcode}.g.vcf.gz"
        def gvcf_tbi_path = "${params.data}/samples/${barcode}/deepvariant/${barcode}.g.vcf.gz.tbi"
        def vaf_bedgraph_path = "${params.data}/samples/${barcode}/sequences/${barcode}.vaf.bedgraph.gz"
        def vaf_bedgraph_tbi_path = "${params.data}/samples/${barcode}/sequences/${barcode}.vaf.bedgraph.gz.tbi"
        
        if (new File(gvcf_path).exists() && new File(gvcf_tbi_path).exists() && 
            new File(vaf_bedgraph_path).exists() && new File(vaf_bedgraph_tbi_path).exists()) {
            plan.deepvariant_sample.existing.add(barcode)
        } else {
            // Only need DeepVariant if the individual's family needs deepvariant_family
            def fid = family_members[barcode]
            if (fid in plan.deepvariant_family.needed) {
                plan.deepvariant_sample.needed.add(barcode)
            }
        }
    }
    
    // Check existing CRAM and bedgraph files (both produced by alignment workflow)
    individuals.each { barcode ->
        def cram_path = "${params.data}/samples/${barcode}/sequences/${barcode}.${params.ref_name}.cram"
        def crai_path = "${params.data}/samples/${barcode}/sequences/${barcode}.${params.ref_name}.cram.crai"
        def bedgraph_path = "${params.data}/samples/${barcode}/sequences/${barcode}.by${params.bin}.bedgraph.gz"
        def bedgraph_tbi_path = "${params.data}/samples/${barcode}/sequences/${barcode}.by${params.bin}.bedgraph.gz.tbi"
        
        if (new File(cram_path).exists() && new File(crai_path).exists() &&
            new File(bedgraph_path).exists() && new File(bedgraph_tbi_path).exists()) {
            plan.alignment.existing.add(barcode)
        } else {
            // Only need alignment if the individual needs DeepVariant
            if (barcode in plan.deepvariant_sample.needed) {
                plan.alignment.needed.add(barcode)
            }
        }
    }
    
    // Check existing WisecondorX NPZ and predict files
    individuals.each { barcode ->
        def npz_path = "${params.data}/samples/${barcode}/svs/wisecondorx/${barcode}.npz"
        def predict_bed_path = "${params.data}/samples/${barcode}/svs/wisecondorx/${barcode}_aberrations.bed"
        def chr_bed_path = "${params.data}/samples/${barcode}/svs/wisecondorx/${barcode}_aberrations.chr.bed"
        
        def npz_exists = new File(npz_path).exists()
        def predict_exists = new File(predict_bed_path).exists()
        def chr_exists = new File(chr_bed_path).exists()
        
        if (predict_exists && chr_exists) {
            // Both predict and chr reformat are done
            plan.wisecondorx.existing.add(barcode)
        } else {
            // Need predict if we have or will have CRAM files
            if (barcode in plan.alignment.existing || barcode in plan.alignment.needed) {
                plan.wisecondorx.needed.add(barcode)
                plan.wisecondorx.need_predict.add(barcode)
                
                // If NPZ doesn't exist, also need NPZ conversion
                if (!npz_exists) {
                    plan.wisecondorx.need_npz.add(barcode)
                }
            }
        }
    }
    
    // Check existing WisecondorX family aberrations files
    families.each { fid ->
        def family_aberrations_path = "${params.data}/families/${fid}/svs/wisecondorx/${fid}_aberrations.bed"
        def family_aberrations_exists = new File(family_aberrations_path).exists()
        
        if (!family_aberrations_exists) {
            // Check if any family members have or will have individual aberrations
            def family_has_aberrations = individuals.any { barcode ->
                family_members[barcode] == fid && 
                (barcode in plan.wisecondorx.existing || barcode in plan.wisecondorx.needed)
            }
            
            if (family_has_aberrations) {
                plan.wisecondorx.need_family_merge[fid] = true
            }
        }
    }
    
    // Check existing WisecondorX annotated family aberrations files
    families.each { fid ->
        def annotated_aberrations_path = "${params.data}/families/${fid}/svs/wisecondorx/${fid}_aberrations.annotated.bed"
        def annotated_aberrations_exists = new File(annotated_aberrations_path).exists()
        
        if (!annotated_aberrations_exists) {
            // Check if family has or will have merged aberrations
            def family_aberrations_path = "${params.data}/families/${fid}/svs/wisecondorx/${fid}_aberrations.bed"
            def has_family_aberrations = new File(family_aberrations_path).exists() || plan.wisecondorx.need_family_merge[fid] == true
            
            if (has_family_aberrations) {
                plan.wisecondorx.need_family_annotate[fid] = true
            }
        }
    }
    
    // Check existing WisecondorX cohort aberrations file
    def cohort_aberrations_path = "${params.data}/cohorts/${params.cohort_name}/svs/wisecondorx/${params.cohort_name}_aberrations.bed"
    def cohort_aberrations_exists = new File(cohort_aberrations_path).exists()
    
    if (!cohort_aberrations_exists) {
        // Check if any families have or will have annotated family aberrations
        def cohort_has_families = families.any { fid ->
            def annotated_family_aberrations_path = "${params.data}/families/${fid}/svs/wisecondorx/${fid}_aberrations.annotated.bed"
            new File(annotated_family_aberrations_path).exists() || plan.wisecondorx.need_family_annotate[fid] == true
        }
        
        if (cohort_has_families) {
            plan.wisecondorx.need_cohort_merge = true
        }
    }
    
    // Check existing Wombat files
    families.each { fid ->
        // Check if BCF2PARQUET output exists
        def bcf2parquet_output_path = "${params.data}/families/${fid}/wombat/${fid}.rare.${params.vep_config_name}.annotated.parquet"
        def bcf2parquet_exists = new File(bcf2parquet_output_path).exists()

        // Check if PYWOMBAT outputs exist (need to check for all configs if defined)
        def pywombat_complete = false
        if (params.wombat_config_list && !params.wombat_config_list.isEmpty() && bcf2parquet_exists) {
            // Check if all PYWOMBAT outputs exist
            pywombat_complete = params.wombat_config_list.every { config_file ->
                def config_name = config_file.replaceAll(/\.ya?ml$/, '')
                def pywombat_output_path = "${params.data}/families/${fid}/wombat/${fid}.rare.${params.vep_config_name}.annotated.${config_name}.tsv"
                new File(pywombat_output_path).exists()
            }
        }

        // Track BCF2PARQUET status
        plan.wombat.need_bcf2parquet[fid] = !bcf2parquet_exists
        
        // Determine if family needs Wombat processing
        if (fid in plan.annotation.existing || fid in plan.annotation.needed) {
            if (pywombat_complete) {
                plan.wombat.existing.add(fid)
            } else {
                plan.wombat.needed.add(fid)
            }
        }
    }
    
    return plan
}

def displayAnalysisSummary(analysis_plan) {
    log.info """
    ========================================================================================
                                    ANALYSIS SUMMARY
    ========================================================================================
    ALIGNMENT: ${analysis_plan.alignment.existing.size()} individuals done and ${analysis_plan.alignment.needed.size()} to do
    == SNVs/INDELs Calling ==
    DEEPVARIANT_SAMPLE: ${analysis_plan.deepvariant_sample.existing.size()} individuals done and ${analysis_plan.deepvariant_sample.needed.size()} to do
    DEEPVARIANT_FAMILY: ${analysis_plan.deepvariant_family.existing.size()} families done and ${analysis_plan.deepvariant_family.needed.size()} to do
    ANNOTATION: ${analysis_plan.annotation.existing.size()} families done and ${analysis_plan.annotation.needed.size()} to do
    WOMBAT: ${analysis_plan.wombat.existing.size()} families done and ${analysis_plan.wombat.needed.size()} to do
    == Common Variants ==
    SNVS_COHORT: common variants cohort bcf is needed: ${analysis_plan.snvs_cohort.needed.size() > 0 ? 'Yes' : 'No'}
    == SVs Calling ==
    WISECONDORX PREDICT: ${analysis_plan.wisecondorx.existing.size()} individuals done and ${analysis_plan.wisecondorx.needed.size()} to do
    == Other ==
    EXTRACTOR: ${analysis_plan.extractor.tsv_count > 0 ? "${analysis_plan.extractor.tsv_count} TSV files to process on ${analysis_plan.extractor.families.size()} families / ${analysis_plan.extractor.samples.size()} samples" : 'Skipped (no TSV files provided)'}
    ========================================================================================
    """
    
    if (analysis_plan.alignment.needed) {
        log.info "Individuals needing alignment: ${analysis_plan.alignment.needed.join(', ')}"
    }
    if (analysis_plan.deepvariant_sample.needed) {
        log.info "Individuals needing variant calling: ${analysis_plan.deepvariant_sample.needed.join(', ')}"
    }
    if (analysis_plan.deepvariant_family.needed) {
        log.info "Families needing family calling, normalization, and pedigree: ${analysis_plan.deepvariant_family.needed.join(', ')}"
    }
    if (analysis_plan.annotation.needed) {
        log.info "Families needing annotation (gnomAD annotation, filtering, and VEP): ${analysis_plan.annotation.needed.join(', ')}"
    }
    if (analysis_plan.snvs_cohort.needed) {
        log.info "Cohort needing common variant merge: Yes"
    }
}

def validateStepsAvailability(analysis_plan) {
    def errors = []
    
    // Check if alignment is needed but not available
    if (analysis_plan.alignment.needed.size() > 0 && !('alignment' in params.steps)) {
        errors.add("Alignment step is required for ${analysis_plan.alignment.needed.size()} individuals but not included in steps parameter")
    }
    
    // Check if deepvariant_sample is needed but not available  
    if (analysis_plan.deepvariant_sample.needed.size() > 0 && !('deepvariant_sample' in params.steps)) {
        errors.add("DeepVariant sample step is required for ${analysis_plan.deepvariant_sample.needed.size()} individuals but not included in steps parameter")
    }
    
    // Check if deepvariant_family is needed but not available
    if (analysis_plan.deepvariant_family.needed.size() > 0 && !('deepvariant_family' in params.steps)) {
        errors.add("DeepVariant family step is required for ${analysis_plan.deepvariant_family.needed.size()} families but not included in steps parameter")
    }
    
    // Check if annotation is needed but not available
    if (analysis_plan.annotation.needed.size() > 0 && !('annotation' in params.steps)) {
        errors.add("Annotation step is required for ${analysis_plan.annotation.needed.size()} families but not included in steps parameter")
    }
    
    if (errors) {
        log.error """
        ========================================================================================
                                        ERRORS DETECTED
        ========================================================================================
        ${errors.join('\n        ')}
        
        Please add the required steps to your parameters or ensure all required files exist.
        Available steps: alignment, deepvariant_sample, deepvariant_family, annotation, snvs_cohort
        ========================================================================================
        """
        exit 1, "Pipeline stopped due to missing required steps"
    }
}

def createChannels(analysis_plan) {
    def channels = [:]
    
    // Create FASTQ channel for alignment
    if (params.fastq_pattern && analysis_plan.alignment.needed.size() > 0) {
        channels.fastq_files = Channel
            .fromFilePairs("${params.data}/fastq/${params.fastq_pattern}", size: 2)
            .map { sample_id, reads ->
                def parts = reads[0].name.tokenize('_')
                def barcode = parts[4]
                def project = parts[2][-3..-1]
                def flowcell = reads[0].name.tokenize('.')[0].tokenize('_')[-1]
                def dual = reads[0].name.tokenize('.')[1]
                def lane = parts[5]
                def unit = "${barcode}_${flowcell}_${lane}"
                
                [barcode, unit, reads[0], reads[1], project, flowcell, dual, lane]
            }
            .filter { barcode, unit, r1, r2, project, flowcell, dual, lane -> 
                barcode in analysis_plan.alignment.needed 
            }
    } else {
        channels.fastq_files = Channel.empty()
    }
    
    // Create CRAM channel for GRCh37 realignment
    if (params.old_cram_37 && analysis_plan.alignment.needed.size() > 0) {
        channels.cram_37_files = Channel
            .fromPath("${params.old_cram_37}/*.cram")
            .map { cram ->
                def barcode = cram.name.tokenize('.')[0]
                def crai_path = "${cram}.crai"
                [barcode, cram, crai_path]
            }
            .filter { barcode, cram, crai_path -> 
                barcode in analysis_plan.alignment.needed && new File(crai_path).exists()
            }
            .map { barcode, cram, crai_path ->
                [barcode, cram, file(crai_path)]
            }
    } else {
        channels.cram_37_files = Channel.empty()
    }
    
    // Create CRAM channel for GRCh38 realignment
    if (params.old_cram_38 && analysis_plan.alignment.needed.size() > 0) {
        channels.cram_38_files = Channel
            .fromPath("${params.old_cram_38}/*.cram")
            .map { cram ->
                def barcode = cram.name.tokenize('.')[0]
                def crai_path = "${cram}.crai"
                [barcode, cram, crai_path]
            }
            .filter { barcode, cram, crai_path -> 
                barcode in analysis_plan.alignment.needed && new File(crai_path).exists()
            }
            .map { barcode, cram, crai_path ->
                [barcode, cram, file(crai_path)]
            }
    } else {
        channels.cram_38_files = Channel.empty()
    }
    
    // Create channel for existing CRAM files
    channels.existing_crams = Channel
        .fromPath("${params.data}/samples/*/sequences/*.${params.ref_name}.cram")
        .map { cram ->
            def barcode = cram.name.tokenize('.')[0]
            def crai_path = "${cram}.crai"
            [barcode, cram, crai_path]
        }
        .filter { barcode, cram, crai_path -> 
            barcode in analysis_plan.alignment.existing && new File(crai_path).exists()
        }
        .map { barcode, cram, crai_path ->
            [barcode, cram, file(crai_path)]
        }
    
    // Create channel for existing gVCF files
    channels.existing_gvcfs = Channel
        .fromPath("${params.data}/samples/*/deepvariant/*.g.vcf.gz")
        .map { gvcf ->
            def barcode = gvcf.name.tokenize('.')[0]
            def tbi_path = "${gvcf}.tbi"
            [barcode, gvcf, tbi_path]
        }
        .filter { barcode, gvcf, tbi_path -> 
            barcode in analysis_plan.deepvariant_sample.existing && new File(tbi_path).exists()
        }
        .map { barcode, gvcf, tbi_path ->
            [barcode, gvcf, file(tbi_path)]
        }

    // Create channel for existing individual VCF files (not gVCF)
    channels.existing_vcfs = Channel
        .fromPath("${params.data}/samples/*/deepvariant/*.vcf.gz")
        .filter { vcf -> !vcf.name.contains('.g.vcf.gz') }  // Exclude gVCF files
        .map { vcf ->
            def barcode = vcf.name.tokenize('.')[0]
            def tbi_path = "${vcf}.tbi"
            [barcode, vcf, tbi_path]
        }
        .filter { barcode, vcf, tbi_path -> 
            barcode in analysis_plan.deepvariant_sample.existing && new File(tbi_path).exists()
        }
        .map { barcode, vcf, tbi_path ->
            [barcode, vcf, file(tbi_path)]
        }
    
    // Create channel for existing family VCF files
    channels.existing_family_vcfs = Channel
        .fromPath("${params.data}/families/*/vcfs/*.vcf.gz")
        .filter { vcf -> !vcf.name.contains(params.vep_config_name) && !vcf.name.contains('.norm.') }  // Exclude VEP annotated and normalized files
        .map { vcf ->
            def fid = vcf.parent.parent.name  // Get family ID from path
            def tbi_path = "${vcf}.tbi"
            [fid, vcf, tbi_path]
        }
        .filter { fid, vcf, tbi_path -> 
            fid in analysis_plan.deepvariant_family.existing && new File(tbi_path).exists()
        }
        .map { fid, vcf, tbi_path ->
            [fid, vcf, file(tbi_path)]
        }
    
    // Create channel for existing normalized BCF files
    channels.existing_normalized_bcfs = Channel
        .fromPath("${params.data}/families/*/vcfs/*.norm.bcf")
        .map { bcf ->
            def fid = bcf.parent.parent.name  // Get family ID from path
            def csi_path = "${bcf}.csi"
            [fid, bcf, csi_path]
        }
        .filter { fid, bcf, csi_path -> 
            fid in analysis_plan.deepvariant_family.existing && new File(csi_path).exists()
        }
        .map { fid, bcf, csi_path ->
            [fid, bcf, file(csi_path)]
        }
    
    // Create channel for existing family pedigree files
    channels.existing_family_pedigrees = Channel
        .fromPath("${params.data}/families/*/*.pedigree.tsv")
        .map { pedigree ->
            def fid = pedigree.parent.name  // Get family ID from path
            [fid, pedigree]
        }
        .filter { fid, pedigree -> 
            fid in analysis_plan.deepvariant_family.existing
        }
    
    // Create channel for existing common filtered BCF files (output from annotation, used by snvs_cohort)
    channels.existing_common_filtered_bcfs = Channel
        .fromPath("${params.data}/families/*/vcfs/*.common_gt.bcf")
        .map { bcf ->
            def fid = bcf.parent.parent.name  // Get family ID from path
            def csi_path = "${bcf}.csi"
            [fid, bcf, csi_path]
        }
        .filter { fid, bcf, csi_path -> 
            fid in analysis_plan.annotation.existing && new File(csi_path).exists()
        }
        .map { fid, bcf, csi_path ->
            [fid, bcf, file(csi_path)]
        }
    
    // Create channel for existing annotated BCF files (output from annotation, used by wombat)
    channels.existing_annotated_bcfs = Channel
        .fromPath("${params.data}/families/*/vcfs/*.rare.${params.vep_config_name}.annotated.bcf")
        .map { bcf ->
            def fid = bcf.parent.parent.name  // Get family ID from path
            def csi_path = "${bcf}.csi"
            [fid, bcf, csi_path]
        }
        .filter { fid, bcf, csi_path -> 
            fid in analysis_plan.annotation.existing && new File(csi_path).exists()
        }
        .map { fid, bcf, csi_path ->
            [fid, bcf, file(csi_path)]
        }
    
    // Create channel for existing Wombat files (output from wombat, used by snvs_cohort)
    if (params.wombat_config_list && !params.wombat_config_list.isEmpty()) {
        channels.existing_wombat_files = Channel.fromList(params.wombat_config_list)
            .map { config_file ->
                def config_name = config_file.replaceAll(/\.ya?ml$/, '')
                // For each config, create entries for each family that has it
                analysis_plan.wombat.existing.collect { fid ->
                    def wombat_file = file("${params.data}/families/${fid}/wombat/${fid}.rare.${params.vep_config_name}.annotated.${config_name}.tsv")
                    if (wombat_file.exists()) {
                        tuple(fid, config_name, wombat_file)
                    } else {
                        null
                    }
                }.findAll { it != null }
            }
            .flatMap()
    } else {
        channels.existing_wombat_files = Channel.empty()
    }
    
    // Create channel for existing NPZ files (output from NPZ_CONVERT, used by predict)
    channels.existing_npz_files = Channel
        .fromPath("${params.data}/samples/*/svs/wisecondorx/*.npz")
        .map { npz ->
            def barcode = npz.parent.parent.parent.name  // Get barcode from path
            [barcode, npz]
        }
        .filter { barcode, npz ->
            // Only include NPZ files that exist but still need predict
            barcode in analysis_plan.wisecondorx.need_predict && !(barcode in analysis_plan.wisecondorx.need_npz)
        }
    
    return channels
}
