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
include { ANCESTRY } from './workflows/ancestry'

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
    exit 1, "ERROR: --steps parameter is required. Available steps: alignment, deepvariant_sample, deepvariant_family, annotation, snvs_cohort, wisecondorx, wombat, extractor, ancestry"
}

// Validate steps
def valid_steps = ['alignment', 'deepvariant_sample', 'deepvariant_family', 'annotation', 'snvs_cohort', 'wisecondorx', 'wombat', 'extractor', 'ancestry']
def invalid_steps = params.steps - valid_steps
if (invalid_steps) {
    exit 1, "ERROR: Invalid steps specified: ${invalid_steps.join(', ')}. Valid steps are: ${valid_steps.join(', ')}"
}

// The ancestry step reads its panel and weights from paths that have no sensible
// default, and every one of its processes would fail on an empty string.
if ('ancestry' in params.steps) {
    def missing_ancestry_params = ['ancestry_reference', 'ancestry_catalog', 'ancestry_panel_name']
        .findAll { key -> !params[key] }
    if (missing_ancestry_params) {
        exit 1, "ERROR: the 'ancestry' step requires ${missing_ancestry_params.join(', ')} to be set"
    }
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
    
    // Stale family outputs are wrong results rather than missing ones, so the run can be
    // stopped on them. Warning is the default deliberately: a cohort whose pedigree has
    // drifted should not silently become unrunnable without the operator asking for that.
    if (params.pedigree_strict && analysis_plan.pedigree_drift) {
        def n = analysis_plan.pedigree_drift.size()
        exit 1, "ERROR: ${n} ${n == 1 ? 'family has' : 'families have'} stale outputs and pedigree_strict is set - see the STALE FAMILY OUTPUTS summary above"
    }

    // Resolve the input files behind everything the plan wants to align, once, so that the
    // guard below and the channels built afterwards agree on what work is actually possible
    def resolved_inputs = resolveAlignmentInputs(analysis_plan)

    // Validate that required steps are available
    validateStepsAvailability(analysis_plan, resolved_inputs, family_members)
    
    // Create channels for different steps
    def channels = createChannels(analysis_plan, resolved_inputs)
    
    // Run alignment if needed and allowed
    // Also entered when CRAMs exist but their coverage bedgraphs are missing, in which case only
    // MOSDEPTH/TABIX_INDEX run - no realignment is triggered.
    if ((analysis_plan.alignment.needed.size() > 0 || analysis_plan.alignment.need_bedgraph.size() > 0) && 'alignment' in params.steps) {
        if (analysis_plan.alignment.needed.size() > 0) {
            log.info "Running alignment for ${analysis_plan.alignment.needed.size()} individuals..."
        }
        if (analysis_plan.alignment.need_bedgraph.size() > 0) {
            log.info "Generating coverage bedgraphs for ${analysis_plan.alignment.need_bedgraph.size()} existing CRAMs..."
        }

        ALIGNMENT(
            channels.fastq_files,
            channels.cram_37_files,
            channels.cram_38_files,
            channels.bedgraph_only_crams
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

        // Get all available normalized BCFs (existing + newly created)
        all_available_normalized_bcfs_wombat = channels.existing_normalized_bcfs.mix(normalized_bcfs_output ?: Channel.empty())

        // Filter for families that need Wombat
        wombat_normalized_bcfs = all_available_normalized_bcfs_wombat
            .filter { fid, bcf, csi ->
                fid in analysis_plan.wombat.needed
            }

        WOMBAT(wombat_bcfs, wombat_pedigrees, wombat_normalized_bcfs, analysis_plan.wombat.need_bcf2parquet)
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
                    analysis_plan.snvs_cohort.need_wombat_merges,
                    families,
                    analysis_plan.annotation.needed,
                    analysis_plan.wombat.needed)
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
            analysis_plan.wisecondorx.need_cohort_merge,
            pedigree_data.families
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
            .fromPath("${params.data}/families/*/*/*/wombat/*.rare.${params.vep_config_name}.annotated.parquet")
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

    // Run ancestry and PGS if needed and allowed.
    //
    // This reads gVCFs directly rather than any merged call set: the panel sites are
    // genotyped per sample from the gVCF's reference blocks, so a site with coverage
    // and no variant is 0/0 rather than absent. A family's own common_gt.bcf carries
    // only its variant sites - around 57% of the panel for a trio - which is below
    // what admixture accepts and enough to distort the projected PCs.
    if ((analysis_plan.ancestry.needed.size() > 0 ||
         analysis_plan.ancestry.need_cohort_merge) && 'ancestry' in params.steps) {

        def ancestry_tasks = []
        if (analysis_plan.ancestry.need_extract.size() > 0) {
            ancestry_tasks.add("panel extraction for ${analysis_plan.ancestry.need_extract.size()} samples")
        }
        def ancestry_merge_count = analysis_plan.ancestry.need_family_merge.count { _fid, needed -> needed == true } ?: 0
        if (ancestry_merge_count > 0) {
            ancestry_tasks.add("family panel merge for ${ancestry_merge_count} families")
        }
        def ancestry_score_count = analysis_plan.ancestry.need_family_score.count { _fid, needed -> needed == true } ?: 0
        if (ancestry_score_count > 0) {
            ancestry_tasks.add("scoring for ${ancestry_score_count} families")
        }
        if (analysis_plan.ancestry.need_cohort_merge) {
            ancestry_tasks.add("cohort table merge")
        }
        log.info "Running ancestry/PGS: ${ancestry_tasks.join(', ')}..."

        ANCESTRY(
            all_available_gvcfs,
            pedigree_data.family_members,
            pedigree_data.families,
            analysis_plan.ancestry.need_extract,
            analysis_plan.ancestry.need_family_merge,
            analysis_plan.ancestry.need_family_score,
            analysis_plan.ancestry.need_cohort_merge
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
        alignment: [needed: [], existing: [], need_bedgraph: []],
        deepvariant_sample: [needed: [], existing: []],
        deepvariant_family: [needed: [], existing: []],
        annotation: [needed: [], existing: []],
        snvs_cohort: [needed: [], existing: [], need_bcf_merge: false, need_wombat_merges: [:]],
        wisecondorx: [needed: [], existing: [], need_npz: [], need_predict: [], need_family_merge: [:], need_family_annotate: [:], need_cohort_merge: false],
        wombat: [needed: [], existing: [], need_bcf2parquet: [:]],
        extractor: [tsv_count: 0, families: [] as Set, samples: [] as Set],
        ancestry: [needed: [], existing: [], need_extract: [], need_family_merge: [:], need_family_score: [:], need_cohort_merge: false],
        pedigree_drift: [:]
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
        def fam_dir = Sharding.getFamilyDir(params.data, fid)
        def norm_bcf_path = "${fam_dir}/vcfs/${fid}.norm.bcf"
        def norm_csi_path = "${fam_dir}/vcfs/${fid}.norm.bcf.csi"
        def pedigree_path = "${fam_dir}/${fid}.pedigree.tsv"

        // length() > 0 covers both "missing" and "present but empty": a 0-byte pedigree must not
        // count as done, otherwise FAMILIAL_PEDIGREE is skipped forever and it never self-heals
        def pedigree_ok = new File(pedigree_path).length() > 0
        if (new File(pedigree_path).exists() && !pedigree_ok) {
            log.warn "Family ${fid} has an empty pedigree at ${pedigree_path} - it will be regenerated"
        }

        if (new File(norm_bcf_path).exists() && new File(norm_csi_path).exists() && pedigree_ok) {
            plan.deepvariant_family.existing.add(fid)

            // The family has been called, but the pedigree may have gained members since.
            // A member with no gVCF cannot possibly be in the joint call, so the family's
            // outputs are stale - and nothing will fix that on its own, because the very
            // existence of norm.bcf is what stops the missing members being scheduled.
            //
            // Membership is inferred from which members have a gVCF rather than read from
            // the BCF: bcftools is not on the launching node's PATH, and the family's own
            // {FID}.pedigree.tsv is no help because it is rewritten from the current
            // pedigree and so already lists the members the call is missing.
            def members = family_members.findAll { _barcode, member_fid -> member_fid == fid }.keySet()
            def without_gvcf = members.findAll { barcode ->
                !new File("${Sharding.getSampleDir(params.data, barcode)}/deepvariant/${barcode}.g.vcf.gz").exists()
            }.sort()
            if (without_gvcf) {
                plan.pedigree_drift[fid] = [members: members.size(), missing: without_gvcf]
            }
        } else {
            plan.deepvariant_family.needed.add(fid)
        }
    }
    
    // Check existing annotation outputs (rare/common VCF.gz/BCFs, common_gt BCF, VEP VCF.gz, and final annotated BCF - all part of annotation)
    // Note: {FID}.gnomad.bcf is intermediate and not published
    families.each { fid ->
        def fam_dir = Sharding.getFamilyDir(params.data, fid)
        def rare_vcf_path = "${fam_dir}/vcfs/${fid}.rare.vcf.gz"
        def rare_tbi_path = "${fam_dir}/vcfs/${fid}.rare.vcf.gz.tbi"
        def common_bcf_path = "${fam_dir}/vcfs/${fid}.common.bcf"
        def common_csi_path = "${fam_dir}/vcfs/${fid}.common.bcf.csi"
        def common_gt_bcf_path = "${fam_dir}/vcfs/${fid}.common_gt.bcf"
        def common_gt_csi_path = "${fam_dir}/vcfs/${fid}.common_gt.bcf.csi"
        def vep_vcf_path = "${fam_dir}/vcfs/${fid}.rare.${params.vep_config_name}.vcf.gz"
        def vep_tbi_path = "${fam_dir}/vcfs/${fid}.rare.${params.vep_config_name}.vcf.gz.tbi"
        def annotated_bcf_path = "${fam_dir}/vcfs/${fid}.rare.${params.vep_config_name}.annotated.bcf"
        def annotated_csi_path = "${fam_dir}/vcfs/${fid}.rare.${params.vep_config_name}.annotated.bcf.csi"
        
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
    
    // Check existing ancestry outputs (per-sample panel genotypes, family panel
    // genotypes, family tables, cohort tables - all part of ancestry)
    //
    // The panel label is part of every file name on purpose: the depth/quality
    // thresholds and the reference bundle are not recorded anywhere this check can
    // see, so bumping ancestry_panel_name is what invalidates the old extractions.
    if ('ancestry' in params.steps) {
        def panel_name = params.ancestry_panel_name
        def table_kinds = ['pcs', 'ancestry', 'Q', 'pgs_raw', 'pgs_adjusted', 'pgs_zscore']

        individuals.each { barcode ->
            def smp_dir = Sharding.getSampleDir(params.data, barcode)
            def panel_bcf_path = "${smp_dir}/ancestry/${barcode}.panel_gt.${panel_name}.bcf"
            def panel_csi_path = "${panel_bcf_path}.csi"

            if (!(new File(panel_bcf_path).exists() && new File(panel_csi_path).exists())) {
                plan.ancestry.need_extract.add(barcode)
            }
        }

        families.each { fid ->
            def fam_dir = Sharding.getFamilyDir(params.data, fid)
            def fam_bcf_path = "${fam_dir}/ancestry/${fid}.panel_gt.${panel_name}.bcf"
            def fam_csi_path = "${fam_bcf_path}.csi"
            def fam_panel_exists = new File(fam_bcf_path).exists() && new File(fam_csi_path).exists()

            // A family whose panel BCF predates one of its samples' extractions must
            // be rebuilt, or it would keep a sample that is no longer current
            def members = family_members.findAll { _barcode, member_fid -> member_fid == fid }.keySet()
            def member_needs_extract = members.any { barcode -> barcode in plan.ancestry.need_extract }

            plan.ancestry.need_family_merge[fid] = !fam_panel_exists || member_needs_extract

            def tables_exist = table_kinds.every { kind ->
                new File("${fam_dir}/ancestry/${fid}.${panel_name}.${kind}.tsv").exists()
            }
            plan.ancestry.need_family_score[fid] = !tables_exist || plan.ancestry.need_family_merge[fid]

            if (plan.ancestry.need_family_merge[fid] || plan.ancestry.need_family_score[fid]) {
                plan.ancestry.needed.add(fid)
            } else {
                plan.ancestry.existing.add(fid)
            }
        }

        def cohort_tables_exist = table_kinds.every { kind ->
            new File("${params.data}/cohorts/${params.cohort_name}/ancestry/${params.cohort_name}.${panel_name}.${kind}.tsv").exists()
        }
        plan.ancestry.need_cohort_merge = !cohort_tables_exist ||
            plan.ancestry.need_family_score.any { _fid, needed -> needed == true }
    }

    // Check existing individual gVCF files and VAF bedgraphs (both outputs of deepvariant_sample)
    individuals.each { barcode ->
        def smp_dir = Sharding.getSampleDir(params.data, barcode)
        def gvcf_path = "${smp_dir}/deepvariant/${barcode}.g.vcf.gz"
        def gvcf_tbi_path = "${smp_dir}/deepvariant/${barcode}.g.vcf.gz.tbi"
        def vaf_bedgraph_path = "${smp_dir}/sequences/${barcode}.vaf.bedgraph.gz"
        def vaf_bedgraph_tbi_path = "${smp_dir}/sequences/${barcode}.vaf.bedgraph.gz.tbi"
        
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
    // CRAM availability and bedgraph availability are tracked separately: a sample with a usable
    // CRAM but no coverage bedgraph only needs MOSDEPTH/TABIX_INDEX, never a full realignment.
    individuals.each { barcode ->
        def smp_dir = Sharding.getSampleDir(params.data, barcode)
        def cram_path = "${smp_dir}/sequences/${barcode}.${params.ref_name}.cram"
        def crai_path = "${smp_dir}/sequences/${barcode}.${params.ref_name}.cram.crai"
        def bedgraph_path = "${smp_dir}/sequences/${barcode}.by${params.bin}.bedgraph.gz"
        def bedgraph_tbi_path = "${smp_dir}/sequences/${barcode}.by${params.bin}.bedgraph.gz.tbi"

        def has_cram = new File(cram_path).exists() && new File(crai_path).exists()
        def has_bedgraph = new File(bedgraph_path).exists() && new File(bedgraph_tbi_path).exists()

        if (has_cram) {
            // The CRAM is usable, so every downstream step can consume it
            plan.alignment.existing.add(barcode)
            if (!has_bedgraph) {
                plan.alignment.need_bedgraph.add(barcode)
            }
        } else {
            // Only need alignment if the individual needs DeepVariant
            if (barcode in plan.deepvariant_sample.needed) {
                plan.alignment.needed.add(barcode)
            } else if (new File("${smp_dir}/deepvariant/${barcode}.g.vcf.gz").exists()) {
                // Benign: nothing downstream reads the CRAM. Family calling consumes the
                // gVCF and the ancestry panel extraction reads it directly, so a missing
                // CRAM here costs nothing but the ability to regenerate coverage bedgraphs.
                def bedgraph_note = has_bedgraph ? "" : " (its coverage bedgraph is also absent and cannot be regenerated without the CRAM)"
                log.info "Individual ${barcode} has no ${params.ref_name} CRAM, but its gVCF is present, so nothing downstream needs it${bedgraph_note}"
            } else {
                // Genuinely stuck: no CRAM, no gVCF, and the family is already called, so
                // the plan will never schedule anything for this individual and the family's
                // outputs cannot contain it. See the STALE FAMILY OUTPUTS summary.
                log.warn "Individual ${barcode} has neither a ${params.ref_name} CRAM nor a gVCF, and family ${family_members[barcode]} is already called - nothing will be scheduled for it and the family's outputs cannot include it. Provide input data for it, or remove it from the pedigree"
            }
        }
    }
    
    // Check existing WisecondorX NPZ and predict files
    individuals.each { barcode ->
        def smp_dir = Sharding.getSampleDir(params.data, barcode)
        def npz_path = "${smp_dir}/svs/wisecondorx/${barcode}.${params.wisecondorx_binsize}.npz"
        def predict_bed_path = "${smp_dir}/svs/wisecondorx/${barcode}_aberrations.bed"
        def chr_bed_path = "${smp_dir}/svs/wisecondorx/${barcode}_aberrations.chr.bed"
        
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
        def family_aberrations_path = "${Sharding.getFamilyDir(params.data, fid)}/svs/wisecondorx/${fid}_aberrations.bed"
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
        def fam_dir = Sharding.getFamilyDir(params.data, fid)
        def annotated_aberrations_path = "${fam_dir}/svs/wisecondorx/${fid}_aberrations.annotated.bed"
        def annotated_aberrations_exists = new File(annotated_aberrations_path).exists()

        if (!annotated_aberrations_exists) {
            // Check if family has or will have merged aberrations
            def family_aberrations_path = "${fam_dir}/svs/wisecondorx/${fid}_aberrations.bed"
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
            def annotated_family_aberrations_path = "${Sharding.getFamilyDir(params.data, fid)}/svs/wisecondorx/${fid}_aberrations.annotated.bed"
            new File(annotated_family_aberrations_path).exists() || plan.wisecondorx.need_family_annotate[fid] == true
        }
        
        if (cohort_has_families) {
            plan.wisecondorx.need_cohort_merge = true
        }
    }
    
    // Check existing Wombat files
    families.each { fid ->
        def fam_dir = Sharding.getFamilyDir(params.data, fid)
        // Check if BCF2PARQUET output exists
        def bcf2parquet_output_path = "${fam_dir}/wombat/${fid}.rare.${params.vep_config_name}.annotated.parquet"
        def bcf2parquet_exists = new File(bcf2parquet_output_path).exists()

        // Check if PYWOMBAT outputs exist (need to check for all configs if defined)
        def pywombat_complete = false
        if (params.wombat_config_list && !params.wombat_config_list.isEmpty() && bcf2parquet_exists) {
            // Check if all PYWOMBAT outputs exist
            pywombat_complete = params.wombat_config_list.every { config_file ->
                def config_name = config_file.replaceAll(/\.ya?ml$/, '')
                def pywombat_output_path = "${fam_dir}/wombat/${fid}.rare.${params.vep_config_name}.annotated.${config_name}.tsv"
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
    ALIGNMENT: ${analysis_plan.alignment.existing.size()} individuals done, ${analysis_plan.alignment.needed.size()} to align, ${analysis_plan.alignment.need_bedgraph.size()} needing bedgraph only
    == SNVs/INDELs Calling ==
    DEEPVARIANT_SAMPLE: ${analysis_plan.deepvariant_sample.existing.size()} individuals done and ${analysis_plan.deepvariant_sample.needed.size()} to do
    DEEPVARIANT_FAMILY: ${analysis_plan.deepvariant_family.existing.size()} families done and ${analysis_plan.deepvariant_family.needed.size()} to do
    ANNOTATION: ${analysis_plan.annotation.existing.size()} families done and ${analysis_plan.annotation.needed.size()} to do
    WOMBAT: ${analysis_plan.wombat.existing.size()} families done and ${analysis_plan.wombat.needed.size()} to do
    == Common Variants ==
    SNVS_COHORT: common variants cohort bcf merge: ${analysis_plan.snvs_cohort.need_bcf_merge ? 'Yes' : 'No'} - wombat cohort merges due: ${analysis_plan.snvs_cohort.need_wombat_merges?.count { _k, v -> v == true } ?: 0}
    == SVs Calling ==
    WISECONDORX PREDICT: ${analysis_plan.wisecondorx.existing.size()} individuals done and ${analysis_plan.wisecondorx.needed.size()} to do
    == Ancestry / PGS ==
    ANCESTRY: ${'ancestry' in params.steps ? "${analysis_plan.ancestry.existing.size()} families done and ${analysis_plan.ancestry.needed.size()} to do (${analysis_plan.ancestry.need_extract.size()} samples needing panel extraction)" : 'Skipped (step not requested)'}
    == Other ==
    EXTRACTOR: ${analysis_plan.extractor.tsv_count > 0 ? "${analysis_plan.extractor.tsv_count} TSV files to process on ${analysis_plan.extractor.families.size()} families / ${analysis_plan.extractor.samples.size()} samples" : 'Skipped (no TSV files provided)'}
    ========================================================================================
    """
    
    if (analysis_plan.alignment.needed) {
        log.info "Individuals needing alignment: ${analysis_plan.alignment.needed.join(', ')}"
    }
    if (analysis_plan.alignment.need_bedgraph) {
        log.info "Individuals needing coverage bedgraph only: ${analysis_plan.alignment.need_bedgraph.join(', ')}"
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
    if (analysis_plan.snvs_cohort.need_bcf_merge) {
        log.info "Cohort needing common variant merge: Yes"
    }
    def due_wombat_merges = analysis_plan.snvs_cohort.need_wombat_merges?.findAll { _k, v -> v == true }?.keySet()
    if (due_wombat_merges) {
        log.info "Cohort wombat merges needed for: ${due_wombat_merges.join(', ')}"
    }
    if (analysis_plan.ancestry.needed) {
        log.info "Families needing ancestry/PGS: ${analysis_plan.ancestry.needed.join(', ')}"
    }
    if (analysis_plan.pedigree_drift) {
        def n_drift = analysis_plan.pedigree_drift.size()
        def drift_noun = n_drift == 1 ? "FAMILY" : "FAMILIES"
        def drift_lines = analysis_plan.pedigree_drift.collect { fid, drift ->
            "${fid}: ${drift.missing.size()} of ${drift.members} members have no gVCF (${drift.missing.join(', ')})"
        }
        log.warn "STALE FAMILY OUTPUTS: ${n_drift} ${drift_noun.toLowerCase()} with already-called outputs that cannot contain every member the pedigree lists - details below"
        log.info """
    ========================================================================================
                    STALE FAMILY OUTPUTS: ${n_drift} ${drift_noun}
    ========================================================================================
    These families have already been called, but the pedigree lists members that have no
    gVCF - so the existing family outputs cannot contain them. This does not self-heal: the
    presence of the family's norm.bcf is exactly what stops the missing members from being
    scheduled for alignment or variant calling.

    ${drift_lines.join('\n    ')}

    Confirm which samples a family actually contains with:
      bcftools query -l <data>/families/{S1}/{S2}/<FID>/vcfs/<FID>.norm.bcf

    Then either
      - provide the missing input data and delete that family's norm.bcf, its .csi and the
        annotation/wombat/ancestry outputs built from it, so it is re-called in full; or
      - remove the members that have no data from the pedigree, if they were never sequenced.

    Set pedigree_strict: true to stop the run on this instead of warning.
    ========================================================================================
    """
    }
}

def validateStepsAvailability(analysis_plan, resolved_inputs, family_members) {
    def errors = []
    
    // Check if alignment is needed but not available
    if (analysis_plan.alignment.needed.size() > 0 && !('alignment' in params.steps)) {
        errors.add("Alignment step is required for ${analysis_plan.alignment.needed.size()} individuals but not included in steps parameter")
    }

    // Check that individuals needing alignment actually have an input to align from, otherwise the
    // ALIGNMENT workflow receives empty channels and the run silently does nothing
    if (analysis_plan.alignment.needed.size() > 0 && !params.fastq_pattern && !params.old_cram_37 && !params.old_cram_38) {
        errors.add("Alignment is required for ${analysis_plan.alignment.needed.size()} individuals (${analysis_plan.alignment.needed.join(', ')}) but no input source is configured - set one of fastq_pattern, old_cram_37 or old_cram_38")
    } else {
        // Configured is not the same as present: check that each individual resolves to a file
        errors.addAll(reconcilePlanWithInputs(analysis_plan, resolved_inputs, family_members))
    }

    // Check if deepvariant_sample is needed but not available  
    if (analysis_plan.deepvariant_sample.needed.size() > 0 && !('deepvariant_sample' in params.steps)) {
        errors.add("DeepVariant sample step is required for ${analysis_plan.deepvariant_sample.needed.size()} individuals but not included in steps parameter")
    }
    
    // A missing prerequisite is only an error when something in this run would
    // consume it. The ancestry step reads gVCFs straight from deepvariant_sample and
    // never touches the normalized or annotated call sets, so `steps: ["ancestry"]`
    // must not be blocked by families that have no norm.bcf yet. Every steps list
    // that includes one of the consumers below behaves exactly as before.
    def family_consumers = ['deepvariant_family', 'annotation', 'wombat', 'snvs_cohort', 'extractor']
    def annotation_consumers = ['annotation', 'wombat', 'snvs_cohort']

    // Check if deepvariant_family is needed but not available
    if (analysis_plan.deepvariant_family.needed.size() > 0 && !('deepvariant_family' in params.steps) &&
        family_consumers.any { step -> step in params.steps }) {
        errors.add("DeepVariant family step is required for ${analysis_plan.deepvariant_family.needed.size()} families but not included in steps parameter")
    }
    
    // Check if annotation is needed but not available
    if (analysis_plan.annotation.needed.size() > 0 && !('annotation' in params.steps) &&
        annotation_consumers.any { step -> step in params.steps }) {
        errors.add("Annotation step is required for ${analysis_plan.annotation.needed.size()} families but not included in steps parameter")
    }
    
    if (errors) {
        log.error """
        ========================================================================================
                                        ERRORS DETECTED
        ========================================================================================
        ${errors.join('\n        ')}
        
        Please add the required steps to your parameters or ensure all required files exist.
        Available steps: alignment, deepvariant_sample, deepvariant_family, annotation, snvs_cohort, wisecondorx, wombat, extractor, ancestry
        ========================================================================================
        """
        exit 1, "Pipeline stopped due to missing required steps"
    }
}

// Barcode derivation lives in one place so the input scan and the channels can never
// disagree about which file belongs to which individual. That drift is the whole bug:
// createAnalysisPlan asks whether an OUTPUT cram is absent, the channels glob INPUT
// directories, and nothing used to compare the two.
def barcodeFromCramName(cram_name) {
    return cram_name.tokenize('.')[0]
}

def barcodeFromFastqName(fastq_name) {
    return fastq_name.tokenize('_')[4]
}

def unitFromFastqName(fastq_name) {
    def parts = fastq_name.tokenize('_')
    def flowcell = fastq_name.tokenize('.')[0].tokenize('_')[-1]
    return "${parts[4]}_${flowcell}_${parts[5]}"
}

// Resolve the actual input files behind every individual the plan wants to align.
//
// alignment.needed means "the output CRAM is missing and something downstream wants it" -
// it says nothing about an input existing. Resolving inputs here, eagerly and once, gives
// validateStepsAvailability something real to check and gives createChannels its rows, so
// the plan and the channels describe the same set of work.
def resolveAlignmentInputs(analysis_plan) {
    def resolved = [
        cram_37       : [],
        cram_38       : [],
        resolvable    : [] as Set,
        fastq_barcodes: [] as Set,
        index_missing : [:],
        unpaired_fastq: [:],
        searched      : []
    ]

    if (analysis_plan.alignment.needed.size() == 0) {
        return resolved
    }

    def needed = analysis_plan.alignment.needed as Set

    // Old CRAMs: flat directories, one file per individual
    [['cram_37', params.old_cram_37], ['cram_38', params.old_cram_38]].each { key, dir ->
        if (!dir) {
            return
        }
        def pattern = "${dir}/*.cram"
        resolved.searched.add(pattern)
        def matches = file(pattern)
        if (!(matches instanceof List)) {
            matches = matches ? [matches] : []
        }
        matches.each { cram ->
            def barcode = barcodeFromCramName(cram.name)
            if (!(barcode in needed)) {
                return
            }
            def crai_path = "${cram}.crai"
            if (new File(crai_path).exists()) {
                resolved[key].add([barcode, cram, file(crai_path)])
                resolved.resolvable.add(barcode)
            } else {
                // Recorded rather than dropped: the data is there, only the index is not,
                // and that needs a different fix than missing input
                resolved.index_missing[barcode] = cram.toString()
            }
        }
    }

    // FASTQ: the channel still uses fromFilePairs for the pairing itself, so this scan
    // mirrors it - a barcode is resolvable once some unit of it has both mates on disk,
    // which is the condition for fromFilePairs(size: 2) to emit that unit. createChannels
    // re-checks what was actually emitted, so an approximation here cannot go unnoticed.
    if (params.fastq_pattern) {
        def pattern = "${params.data}/fastq/${params.fastq_pattern}"
        resolved.searched.add(pattern)
        def matches = file(pattern)
        if (!(matches instanceof List)) {
            matches = matches ? [matches] : []
        }
        def unit_count = [:]
        def unit_barcode = [:]
        matches.each { fq ->
            def barcode
            def unit
            try {
                barcode = barcodeFromFastqName(fq.name)
                unit = unitFromFastqName(fq.name)
            } catch (Exception _e) {
                log.warn "Ignoring FASTQ ${fq.name}: its name does not follow the expected underscore-separated layout, so no barcode could be read from it"
                return
            }
            if (!(barcode in needed)) {
                return
            }
            unit_count[unit] = (unit_count[unit] ?: 0) + 1
            unit_barcode[unit] = barcode
        }
        unit_count.each { unit, n ->
            if (n == 2) {
                resolved.resolvable.add(unit_barcode[unit])
                resolved.fastq_barcodes.add(unit_barcode[unit])
            } else {
                resolved.unpaired_fastq[unit] = n
            }
        }
        resolved.unpaired_fastq.each { unit, n ->
            log.warn "FASTQ unit ${unit} has ${n} read file(s) instead of 2 and will be ignored by fromFilePairs"
        }
    }

    return resolved
}

// Walk the plan the way the channels will, so an individual with no usable input is named
// instead of vanishing. Nothing downstream can rescue it: an individual that needs
// alignment is by construction absent from alignment.existing, so it never enters
// all_available_crams and every later stage collapses with it, silently and with exit 0.
def reconcilePlanWithInputs(analysis_plan, resolved_inputs, family_members) {
    def errors = []

    if (analysis_plan.alignment.needed.size() == 0) {
        return errors
    }

    def unresolved = analysis_plan.alignment.needed.findAll { barcode -> !(barcode in resolved_inputs.resolvable) }
    if (!unresolved) {
        return errors
    }

    def searched = resolved_inputs.searched ?: ['nothing - no input source is configured']
    // Only hint at the naming rule of the sources actually searched
    def hints = []
    if (params.old_cram_37 || params.old_cram_38) {
        hints.add("A CRAM is matched to an individual by the part of its name before the first '.', so C0XYVSO.cram matches but C0XYVSO_hs38DH.cram does not, and its index must be <file>.cram.crai.")
    }
    if (params.fastq_pattern) {
        hints.add("A FASTQ is matched to an individual by the 5th underscore-separated field of its name, and both mates of a unit must be present.")
    }
    errors.add("No input data found for ${unresolved.size()} of the ${analysis_plan.alignment.needed.size()} individuals that need alignment: ${unresolved.join(', ')}\n" +
               "        Searched: ${searched.join(', ')}" +
               (hints ? "\n        " + hints.join("\n        ") : ""))

    // A CRAM that is present but unindexed is a different problem with a different fix
    def index_missing = resolved_inputs.index_missing.findAll { barcode, _cram -> barcode in unresolved }
    if (index_missing) {
        errors.add("Input CRAM found but unusable for ${index_missing.size()} individual(s) - the index must exist and be named <file>.cram.crai, run 'samtools index' on each:\n        " +
                   index_missing.collect { barcode, cram -> "${barcode}: ${cram}" }.join("\n        "))
    }

    // Name the downstream work that cannot happen either, so this report covers everything
    // the ANALYSIS SUMMARY above promised
    def blocked = unresolved as Set

    def blocked_samples = analysis_plan.deepvariant_sample.needed.findAll { barcode -> barcode in blocked }
    if (blocked_samples) {
        errors.add("Blocked by the above - variant calling for ${blocked_samples.size()} individual(s): ${blocked_samples.join(', ')}")
    }

    def members_by_family = [:]
    family_members.each { barcode, fid ->
        members_by_family[fid] = (members_by_family[fid] ?: []) + [barcode]
    }

    def fully_blocked = []
    def partly_blocked = []
    analysis_plan.deepvariant_family.needed.each { fid ->
        def members = members_by_family[fid] ?: []
        def blocked_members = members.findAll { barcode -> barcode in blocked }
        if (!blocked_members) {
            return
        }
        if (blocked_members.size() == members.size()) {
            fully_blocked.add(fid)
        } else {
            partly_blocked.add("${fid} (no data for ${blocked_members.join(', ')} of ${members.size()} members)")
        }
    }
    if (fully_blocked) {
        errors.add("Blocked by the above - family calling, annotation and wombat for ${fully_blocked.size()} family/families with no usable member: ${fully_blocked.join(', ')}")
    }
    if (partly_blocked) {
        // Worth its own line: this is the set that would otherwise be joint-called from an
        // incomplete pedigree rather than not called at all
        errors.add("Blocked by the above - ${partly_blocked.size()} family/families would be joint-called from an incomplete set of members: ${partly_blocked.join('; ')}")
    }

    return errors
}

def createChannels(analysis_plan, resolved_inputs) {
    def channels = [:]
    
    // Create FASTQ channel for alignment
    //
    // fromFilePairs owns the pairing, so this cannot be built from the eager scan the way
    // the CRAM channels are. Instead the scan's verdict is re-checked against what the
    // channel actually emitted: toList/flatMap keeps that to a single consumption and
    // still aborts before any alignment task is submitted.
    if (params.fastq_pattern && analysis_plan.alignment.needed.size() > 0) {
        def fastq_expected = resolved_inputs.fastq_barcodes
        channels.fastq_files = Channel
            .fromFilePairs("${params.data}/fastq/${params.fastq_pattern}", size: 2)
            .map { sample_id, reads ->
                def parts = reads[0].name.tokenize('_')
                def barcode = barcodeFromFastqName(reads[0].name)
                def project = parts[2][-3..-1]
                def flowcell = reads[0].name.tokenize('.')[0].tokenize('_')[-1]
                def dual = reads[0].name.tokenize('.')[1]
                def lane = parts[5]
                def unit = unitFromFastqName(reads[0].name)
                
                [barcode, unit, reads[0], reads[1], project, flowcell, dual, lane]
            }
            .filter { barcode, unit, r1, r2, project, flowcell, dual, lane -> 
                barcode in analysis_plan.alignment.needed 
            }
            .toList()
            .flatMap { rows ->
                def emitted = rows.collect { row -> row[0] } as Set
                def dropped = fastq_expected.findAll { barcode -> !(barcode in emitted) }
                if (dropped) {
                    error "FASTQ pairing dropped ${dropped.size()} individual(s) that the input scan found data for: ${dropped.join(', ')} - check that every read file has its mate and that '${params.fastq_pattern}' pairs them"
                }
                rows
            }
    } else {
        channels.fastq_files = Channel.empty()
    }
    
    // Create CRAM channels for realignment, from the single resolution pass rather than a
    // second glob - re-globbing here is what let the channels and the plan drift apart
    channels.cram_37_files = Channel.fromList(resolved_inputs.cram_37)
    channels.cram_38_files = Channel.fromList(resolved_inputs.cram_38)
    
    // Create channel for existing CRAM files
    channels.existing_crams = Channel
        .fromPath("${params.data}/samples/*/*/*/sequences/*.${params.ref_name}.cram")
        .map { cram ->
            def barcode = barcodeFromCramName(cram.name)
            def crai_path = "${cram}.crai"
            [barcode, cram, crai_path]
        }
        .filter { barcode, cram, crai_path ->
            barcode in analysis_plan.alignment.existing && new File(crai_path).exists()
        }
        .map { barcode, cram, crai_path ->
            [barcode, cram, file(crai_path)]
        }

    // Create channel for existing CRAMs that only need a coverage bedgraph (no realignment)
    channels.bedgraph_only_crams = Channel
        .fromPath("${params.data}/samples/*/*/*/sequences/*.${params.ref_name}.cram")
        .map { cram ->
            def barcode = barcodeFromCramName(cram.name)
            def crai_path = "${cram}.crai"
            [barcode, cram, crai_path]
        }
        .filter { barcode, _cram, crai_path ->
            barcode in analysis_plan.alignment.need_bedgraph && new File(crai_path).exists()
        }
        .map { barcode, cram, crai_path ->
            [barcode, cram, file(crai_path)]
        }

    // Create channel for existing gVCF files
    channels.existing_gvcfs = Channel
        .fromPath("${params.data}/samples/*/*/*/deepvariant/*.g.vcf.gz")
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
        .fromPath("${params.data}/samples/*/*/*/deepvariant/*.vcf.gz")
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
        .fromPath("${params.data}/families/*/*/*/vcfs/*.vcf.gz")
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
        .fromPath("${params.data}/families/*/*/*/vcfs/*.norm.bcf")
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
        .fromPath("${params.data}/families/*/*/*/*.pedigree.tsv")
        .map { pedigree ->
            def fid = pedigree.parent.name  // Get family ID from path
            [fid, pedigree]
        }
        .filter { fid, pedigree -> 
            fid in analysis_plan.deepvariant_family.existing
        }
    
    // Create channel for existing common filtered BCF files (output from annotation, used by snvs_cohort)
    channels.existing_common_filtered_bcfs = Channel
        .fromPath("${params.data}/families/*/*/*/vcfs/*.common_gt.bcf")
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
        .fromPath("${params.data}/families/*/*/*/vcfs/*.rare.${params.vep_config_name}.annotated.bcf")
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
                    def wombat_file = file("${Sharding.getFamilyDir(params.data, fid)}/wombat/${fid}.rare.${params.vep_config_name}.annotated.${config_name}.tsv")
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
        .fromPath("${params.data}/samples/*/*/*/svs/wisecondorx/*.${params.wisecondorx_binsize}.npz")
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
