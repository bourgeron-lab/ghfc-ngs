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
include { DEEPVARIANT } from './workflows/deepvariant'
include { FAMILY_CALLING } from './workflows/family_calling'
include { NORM } from './workflows/norm'
include { SNVS_FREQ_ANNOT } from './workflows/snvs_freq_annot'
include { SNVS_FREQ_FILTER } from './workflows/snvs_freq_filter'
include { SNVS_COMMON_FILTERS } from './workflows/snvs_common_filters'
include { SNVS_COMMON_COHORT } from './workflows/snvs_common_cohort'
include { BEDGRAPHS } from './workflows/bedgraphs'
include { VEP_ANNOTATION } from './workflows/vep_annotation'

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
    exit 1, "ERROR: --steps parameter is required. Available steps: alignment, deepvariant, family_calling, norm, snvs_freq_annot, snvs_freq_filter, bedgraphs, vep_annotation"
}

// Validate steps
def valid_steps = ['alignment', 'deepvariant', 'family_calling', 'norm', 'snvs_freq_annot', 'snvs_freq_filter', 'snvs_common_filters', 'snvs_common_cohort', 'bedgraphs', 'vep_annotation']
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
    
    // Run DeepVariant if needed and allowed
    if (analysis_plan.deepvariant.needed.size() > 0 && 'deepvariant' in params.steps) {
        log.info "Running DeepVariant for ${analysis_plan.deepvariant.needed.size()} individuals..."
        
        // Filter CRAM files for individuals that need DeepVariant
        deepvariant_crams = all_available_crams
            .filter { barcode, cram, crai -> 
                barcode in analysis_plan.deepvariant.needed 
            }
        
        DEEPVARIANT(deepvariant_crams)
        
        deepvariant_gvcfs = DEEPVARIANT.out.gvcf
    } else {
        deepvariant_gvcfs = Channel.empty()
    }
    
    // Collect all available gVCF files (existing + newly created)
    all_available_gvcfs = channels.existing_gvcfs.mix(deepvariant_gvcfs ?: Channel.empty())
    
    // Run family calling if needed and allowed
    family_vcfs_output = Channel.empty()
    if (analysis_plan.family_calling.needed.size() > 0 && 'family_calling' in params.steps) {
        log.info "Running family calling for ${analysis_plan.family_calling.needed.size()} families..."
        
        // Group gVCF files by family
        family_gvcfs = all_available_gvcfs
            .map { barcode, gvcf, tbi -> 
                def fid = family_members[barcode]
                [fid, barcode, gvcf, tbi]
            }
            .filter { fid, barcode, gvcf, tbi -> 
                fid in analysis_plan.family_calling.needed 
            }
            .groupTuple(by: 0)
            .map { fid, barcodes, gvcfs, tbis ->
                [fid, barcodes, gvcfs, tbis]
            }
        
        FAMILY_CALLING(family_gvcfs)
        family_vcfs_output = FAMILY_CALLING.out.family_vcfs
    }
    
    // Run VCF normalization if needed and allowed
    normalized_vcfs_output = Channel.empty()
    if (analysis_plan.norm.needed.size() > 0 && 'norm' in params.steps) {
        log.info "Running VCF normalization for ${analysis_plan.norm.needed.size()} families..."
        
        // Get all available family VCFs (existing + newly created)
        all_available_family_vcfs = channels.existing_family_vcfs.mix(family_vcfs_output ?: Channel.empty())
        
        // Filter for families that need normalization
        norm_vcfs = all_available_family_vcfs
            .filter { fid, vcf, tbi -> 
                fid in analysis_plan.norm.needed 
            }
        
        NORM(norm_vcfs)
        normalized_vcfs_output = NORM.out.normalized_vcfs
    }
    
    // Run gnomAD frequency annotation if needed and allowed
    gnomad_vcfs_output = Channel.empty()
    if (analysis_plan.snvs_freq_annot.needed.size() > 0 && 'snvs_freq_annot' in params.steps) {
        log.info "Running gnomAD frequency annotation for ${analysis_plan.snvs_freq_annot.needed.size()} families..."
        
        // Get all available normalized family VCFs (existing + newly created)
        all_available_normalized_vcfs = channels.existing_normalized_vcfs.mix(normalized_vcfs_output ?: Channel.empty())
        
        // Filter for families that need gnomAD annotation
        gnomad_vcfs = all_available_normalized_vcfs
            .filter { fid, vcf, tbi -> 
                fid in analysis_plan.snvs_freq_annot.needed 
            }
        
        SNVS_FREQ_ANNOT(gnomad_vcfs)
        gnomad_vcfs_output = SNVS_FREQ_ANNOT.out.annotated_vcfs
    }
    
    // Run gnomAD frequency filtering if needed and allowed
    rare_vcfs_output = Channel.empty()
    common_vcfs_output = Channel.empty()
    if (analysis_plan.snvs_freq_filter.needed.size() > 0 && 'snvs_freq_filter' in params.steps) {
        log.info "Running rare/common variant separation for ${analysis_plan.snvs_freq_filter.needed.size()} families..."
        
        // Get all available gnomAD annotated family VCFs (existing + newly created)
        all_available_gnomad_vcfs = channels.existing_gnomad_vcfs.mix(gnomad_vcfs_output ?: Channel.empty())
        
        // Filter for families that need frequency filtering
        filter_vcfs = all_available_gnomad_vcfs
            .filter { fid, vcf, tbi -> 
                fid in analysis_plan.snvs_freq_filter.needed 
            }
        
        SNVS_FREQ_FILTER(filter_vcfs)
        rare_vcfs_output = SNVS_FREQ_FILTER.out.rare_vcfs
        common_vcfs_output = SNVS_FREQ_FILTER.out.common_vcfs
    }
    
    // Run common variants filtering if needed and allowed
    if (analysis_plan.snvs_common_filters.needed.size() > 0 && 'snvs_common_filters' in params.steps) {
        log.info "Running common variants filtering for ${analysis_plan.snvs_common_filters.needed.size()} families..."
        
        // Get all available common VCFs (existing + newly created)
        all_available_common_vcfs = channels.existing_common_vcfs.mix(common_vcfs_output ?: Channel.empty())
        
        // Filter for families that need common filtering
        common_filter_vcfs = all_available_common_vcfs
            .filter { fid, vcf, tbi -> 
                fid in analysis_plan.snvs_common_filters.needed 
            }
        
        SNVS_COMMON_FILTERS(common_filter_vcfs)
    }
    
    // Run cohort common variants merge if needed and allowed
    if (analysis_plan.snvs_common_cohort.needed.size() > 0 && 'snvs_common_cohort' in params.steps) {
        log.info "Running cohort common variants merge..."
        
        // Get all available common filtered BCFs (existing + newly created)
        if (analysis_plan.snvs_common_filters.needed.size() > 0 && 'snvs_common_filters' in params.steps) {
            // Mix existing BCFs with newly created ones
            all_available_common_bcfs = channels.existing_common_filtered_bcfs.mix(SNVS_COMMON_FILTERS.out.filtered_common_bcfs)
        } else {
            // Use only existing BCFs if no new ones were created
            all_available_common_bcfs = channels.existing_common_filtered_bcfs
        }
        
        SNVS_COMMON_COHORT(all_available_common_bcfs)
    }
    
    // Run VEP annotation if needed and allowed
    if (analysis_plan.vep_annotation.needed.size() > 0 && 'vep_annotation' in params.steps) {
        log.info "Running VEP annotation for ${analysis_plan.vep_annotation.needed.size()} families..."
        
        // Get all available rare family VCFs (existing + newly created)
        all_available_filtered_vcfs = channels.existing_filtered_vcfs.mix(rare_vcfs_output ?: Channel.empty())
        
        // Filter for families that need VEP annotation
        vep_vcfs = all_available_filtered_vcfs
            .filter { fid, vcf, tbi -> 
                fid in analysis_plan.vep_annotation.needed 
            }
        
        VEP_ANNOTATION(vep_vcfs)
    }
    
    // Run bedgraphs generation if needed and allowed
    if (analysis_plan.bedgraphs.needed.size() > 0 && 'bedgraphs' in params.steps) {
        log.info "Running bedgraph generation for ${analysis_plan.bedgraphs.needed.size()} individuals..."
        
        // Filter CRAM files for individuals that need bedgraphs
        bedgraphs_crams = all_available_crams
            .filter { barcode, cram, crai -> 
                barcode in analysis_plan.bedgraphs.needed 
            }
        
        BEDGRAPHS(bedgraphs_crams)
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
    
    file(pedigree_file).readLines().each { line ->
        if (line.startsWith('#') || line.trim().isEmpty()) return
        
        def cols = line.split('\t')
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
        deepvariant: [needed: [], existing: []],
        family_calling: [needed: [], existing: []],
        norm: [needed: [], existing: []],
        snvs_freq_annot: [needed: [], existing: []],
        snvs_freq_filter: [needed: [], existing: []],
        snvs_common_filters: [needed: [], existing: []],
        snvs_common_cohort: [needed: [], existing: []],
        bedgraphs: [needed: [], existing: []],
        vep_annotation: [needed: [], existing: []]
    ]
    
    // Check existing family VCF files 
    families.each { fid ->
        def family_vcf_path = "${params.data}/families/${fid}/vcfs/${fid}.vcf.gz"
        def family_tbi_path = "${params.data}/families/${fid}/vcfs/${fid}.vcf.gz.tbi"
        
        if (new File(family_vcf_path).exists() && new File(family_tbi_path).exists()) {
            plan.family_calling.existing.add(fid)
        } else {
            plan.family_calling.needed.add(fid)
        }
    }
    
    // Check existing normalized VCF files
    families.each { fid ->
        def norm_vcf_path = "${params.data}/families/${fid}/vcfs/${fid}.norm.vcf.gz"
        def norm_tbi_path = "${params.data}/families/${fid}/vcfs/${fid}.norm.vcf.gz.tbi"
        
        if (new File(norm_vcf_path).exists() && new File(norm_tbi_path).exists()) {
            plan.norm.existing.add(fid)
        } else {
            // Need normalization if family VCF exists or will be created
            if (fid in plan.family_calling.existing || fid in plan.family_calling.needed) {
                plan.norm.needed.add(fid)
            }
        }
    }
    
    // Check existing gnomAD annotated files
    families.each { fid ->
        def gnomad_vcf_path = "${params.data}/families/${fid}/vcfs/${fid}.gnomad.vcf.gz"
        def gnomad_tbi_path = "${params.data}/families/${fid}/vcfs/${fid}.gnomad.vcf.gz.tbi"
        
        if (new File(gnomad_vcf_path).exists() && new File(gnomad_tbi_path).exists()) {
            plan.snvs_freq_annot.existing.add(fid)
        } else {
            // Need gnomAD annotation if normalized VCFs exist or will be created
            if (fid in plan.norm.existing || fid in plan.norm.needed) {
                plan.snvs_freq_annot.needed.add(fid)
            }
        }
    }
    
    // Check existing gnomAD filtered files (rare variants)
    families.each { fid ->
        def rare_vcf_path = "${params.data}/families/${fid}/vcfs/${fid}.rare.vcf.gz"
        def rare_tbi_path = "${params.data}/families/${fid}/vcfs/${fid}.rare.vcf.gz.tbi"
        def common_vcf_path = "${params.data}/families/${fid}/vcfs/${fid}.common.vcf.gz"
        def common_tbi_path = "${params.data}/families/${fid}/vcfs/${fid}.common.vcf.gz.tbi"
        
        if (new File(rare_vcf_path).exists() && new File(rare_tbi_path).exists() && 
            new File(common_vcf_path).exists() && new File(common_tbi_path).exists()) {
            plan.snvs_freq_filter.existing.add(fid)
        } else {
            // Need filtering if gnomAD annotated VCFs exist or will be created
            if (fid in plan.snvs_freq_annot.existing || fid in plan.snvs_freq_annot.needed) {
                plan.snvs_freq_filter.needed.add(fid)
            }
        }
    }

    // Check existing common filtered files
    families.each { fid ->
        def common_bcf_path = "${params.data}/families/${fid}/vcfs/${fid}.common_gt.bcf"
        def common_csi_path = "${params.data}/families/${fid}/vcfs/${fid}.common_gt.bcf.csi"
        
        if (new File(common_bcf_path).exists() && new File(common_csi_path).exists()) {
            plan.snvs_common_filters.existing.add(fid)
        } else {
            // Need common filtering if freq_filter VCFs exist or will be created (common variants available)
            if (fid in plan.snvs_freq_filter.existing || fid in plan.snvs_freq_filter.needed) {
                plan.snvs_common_filters.needed.add(fid)
            }
        }
    }

    // Check existing cohort common BCF file
    def cohort_bcf_path = "${params.data}/cohorts/${params.cohort_name}/vcfs/${params.cohort_name}.common_gt.bcf"
    def cohort_csi_path = "${params.data}/cohorts/${params.cohort_name}/vcfs/${params.cohort_name}.common_gt.bcf.csi"
    
    if (new File(cohort_bcf_path).exists() && new File(cohort_csi_path).exists()) {
        plan.snvs_common_cohort.existing.add('cohort')  // Single cohort entry
    } else {
        // Need cohort merge if any families have common filtered BCFs or will create them
        def families_with_common_bcfs = plan.snvs_common_filters.existing + plan.snvs_common_filters.needed
        if (families_with_common_bcfs.size() > 0) {
            plan.snvs_common_cohort.needed.add('cohort')  // Single cohort entry
        }
    }

    // Check existing VEP annotated files (rare variants)
    families.each { fid ->
        def vep_vcf_path = "${params.data}/families/${fid}/vcfs/${fid}.rare.${params.vep_config_name}.vcf.gz"
        def vep_tbi_path = "${params.data}/families/${fid}/vcfs/${fid}.rare.${params.vep_config_name}.vcf.gz.tbi"

        if (new File(vep_vcf_path).exists() && new File(vep_tbi_path).exists()) {
            plan.vep_annotation.existing.add(fid)
        } else {
            // Need VEP if filtered VCFs exist or will be created
            if (fid in plan.snvs_freq_filter.existing || fid in plan.snvs_freq_filter.needed) {
                plan.vep_annotation.needed.add(fid)
            }
        }
    }
    
    // Check existing individual gVCF files
    individuals.each { barcode ->
        def gvcf_path = "${params.data}/samples/${barcode}/deepvariant/${barcode}.g.vcf.gz"
        def gvcf_tbi_path = "${params.data}/samples/${barcode}/deepvariant/${barcode}.g.vcf.gz.tbi"
        
        if (new File(gvcf_path).exists() && new File(gvcf_tbi_path).exists()) {
            plan.deepvariant.existing.add(barcode)
        } else {
            // Only need DeepVariant if the individual's family needs family calling
            def fid = family_members[barcode]
            if (fid in plan.family_calling.needed) {
                plan.deepvariant.needed.add(barcode)
            }
        }
    }
    
    // Check existing CRAM files
    individuals.each { barcode ->
        def cram_path = "${params.data}/samples/${barcode}/sequences/${barcode}.${params.ref_name}.cram"
        def crai_path = "${params.data}/samples/${barcode}/sequences/${barcode}.${params.ref_name}.cram.crai"
        
        if (new File(cram_path).exists() && new File(crai_path).exists()) {
            plan.alignment.existing.add(barcode)
        } else {
            // Only need alignment if the individual needs DeepVariant
            if (barcode in plan.deepvariant.needed) {
                plan.alignment.needed.add(barcode)
            }
        }
    }
    
    // Check existing bedgraph files
    individuals.each { barcode ->
        def bedgraph_path = "${params.data}/samples/${barcode}/sequences/${barcode}.by${params.bin}.bedgraph.gz"
        def bedgraph_tbi_path = "${params.data}/samples/${barcode}/sequences/${barcode}.by${params.bin}.bedgraph.gz.tbi"
        
        if (new File(bedgraph_path).exists() && new File(bedgraph_tbi_path).exists()) {
            plan.bedgraphs.existing.add(barcode)
        } else {
            // Only need bedgraphs if the step is requested and CRAM files exist or will be created
            if ('bedgraphs' in params.steps && (barcode in plan.alignment.existing || barcode in plan.alignment.needed)) {
                plan.bedgraphs.needed.add(barcode)
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
    ALIGNMENT:
        - Existing CRAM files: ${analysis_plan.alignment.existing.size()} individuals
        - Need alignment: ${analysis_plan.alignment.needed.size()} individuals
    
    BEDGRAPHS:
        - Existing bedgraph files: ${analysis_plan.bedgraphs.existing.size()} individuals
        - Need bedgraph generation: ${analysis_plan.bedgraphs.needed.size()} individuals

    DEEPVARIANT:
        - Existing gVCF files: ${analysis_plan.deepvariant.existing.size()} individuals  
        - Need variant calling: ${analysis_plan.deepvariant.needed.size()} individuals
    
    FAMILY_CALLING:
        - Existing family VCFs: ${analysis_plan.family_calling.existing.size()} families
        - Need family calling: ${analysis_plan.family_calling.needed.size()} families
    
    NORM:
        - Existing normalized VCFs: ${analysis_plan.norm.existing.size()} families
        - Need VCF normalization: ${analysis_plan.norm.needed.size()} families
    
    SNVS_FREQ_ANNOT:
        - Existing gnomAD annotated VCFs: ${analysis_plan.snvs_freq_annot.existing.size()} families
        - Need gnomAD frequency annotation: ${analysis_plan.snvs_freq_annot.needed.size()} families
    
    SNVS_FREQ_FILTER:
        - Existing rare/common VCF pairs: ${analysis_plan.snvs_freq_filter.existing.size()} families
        - Need rare/common variant separation: ${analysis_plan.snvs_freq_filter.needed.size()} families
    
    SNVS_COMMON_FILTERS:
        - Existing filtered common BCFs: ${analysis_plan.snvs_common_filters.existing.size()} families
        - Need common variant filtering: ${analysis_plan.snvs_common_filters.needed.size()} families
    
    SNVS_COMMON_COHORT:
        - Existing cohort common BCF: ${analysis_plan.snvs_common_cohort.existing.size() > 0 ? 'Yes' : 'No'}
        - Need cohort common merge: ${analysis_plan.snvs_common_cohort.needed.size() > 0 ? 'Yes' : 'No'}
    
    VEP_ANNOTATION:
        - Existing VEP annotated files: ${analysis_plan.vep_annotation.existing.size()} families
        - Need VEP annotation: ${analysis_plan.vep_annotation.needed.size()} families
    ========================================================================================
    """
    
    if (analysis_plan.alignment.needed) {
        log.info "Individuals needing alignment: ${analysis_plan.alignment.needed.join(', ')}"
    }
    if (analysis_plan.deepvariant.needed) {
        log.info "Individuals needing variant calling: ${analysis_plan.deepvariant.needed.join(', ')}"
    }
    if (analysis_plan.family_calling.needed) {
        log.info "Families needing family calling: ${analysis_plan.family_calling.needed.join(', ')}"
    }
    if (analysis_plan.norm.needed) {
        log.info "Families needing VCF normalization: ${analysis_plan.norm.needed.join(', ')}"
    }
    if (analysis_plan.snvs_freq_annot.needed) {
        log.info "Families needing gnomAD frequency annotation: ${analysis_plan.snvs_freq_annot.needed.join(', ')}"
    }
    if (analysis_plan.snvs_freq_filter.needed) {
        log.info "Families needing rare/common variant separation: ${analysis_plan.snvs_freq_filter.needed.join(', ')}"
    }
    if (analysis_plan.snvs_common_filters.needed) {
        log.info "Families needing common variant filtering: ${analysis_plan.snvs_common_filters.needed.join(', ')}"
    }
    if (analysis_plan.snvs_common_cohort.needed) {
        log.info "Cohort needing common variant merge: Yes"
    }
    if (analysis_plan.vep_annotation.needed) {
        log.info "Families needing VEP annotation: ${analysis_plan.vep_annotation.needed.join(', ')}"
    }
    if (analysis_plan.bedgraphs.needed) {
        log.info "Individuals needing bedgraph generation: ${analysis_plan.bedgraphs.needed.join(', ')}"
    }
}

def validateStepsAvailability(analysis_plan) {
    def errors = []
    
    // Check if alignment is needed but not available
    if (analysis_plan.alignment.needed.size() > 0 && !('alignment' in params.steps)) {
        errors.add("Alignment step is required for ${analysis_plan.alignment.needed.size()} individuals but not included in steps parameter")
    }
    
    // Check if deepvariant is needed but not available  
    if (analysis_plan.deepvariant.needed.size() > 0 && !('deepvariant' in params.steps)) {
        errors.add("DeepVariant step is required for ${analysis_plan.deepvariant.needed.size()} individuals but not included in steps parameter")
    }
    
    // Check if family_calling is needed but not available
    if (analysis_plan.family_calling.needed.size() > 0 && !('family_calling' in params.steps)) {
        errors.add("Family calling step is required for ${analysis_plan.family_calling.needed.size()} families but not included in steps parameter")
    }
    
    // Check if bedgraphs is needed but not available
    if (analysis_plan.bedgraphs.needed.size() > 0 && !('bedgraphs' in params.steps)) {
        errors.add("Bedgraphs step is required for ${analysis_plan.bedgraphs.needed.size()} individuals but not included in steps parameter")
    }
    
    if (errors) {
        log.error """
        ========================================================================================
                                        ERRORS DETECTED
        ========================================================================================
        ${errors.join('\n        ')}
        
        Please add the required steps to your parameters or ensure all required files exist.
        Available steps: alignment, deepvariant, family_calling, bedgraphs
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
            barcode in analysis_plan.deepvariant.existing && new File(tbi_path).exists()
        }
        .map { barcode, gvcf, tbi_path ->
            [barcode, gvcf, file(tbi_path)]
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
            fid in analysis_plan.family_calling.existing && new File(tbi_path).exists()
        }
        .map { fid, vcf, tbi_path ->
            [fid, vcf, file(tbi_path)]
        }
    
    // Create channel for existing normalized VCF files
    channels.existing_normalized_vcfs = Channel
        .fromPath("${params.data}/families/*/vcfs/*.norm.vcf.gz")
        .map { vcf ->
            def fid = vcf.parent.parent.name  // Get family ID from path
            def tbi_path = "${vcf}.tbi"
            [fid, vcf, tbi_path]
        }
        .filter { fid, vcf, tbi_path -> 
            fid in analysis_plan.norm.existing && new File(tbi_path).exists()
        }
        .map { fid, vcf, tbi_path ->
            [fid, vcf, file(tbi_path)]
        }
    
    // Create channel for existing gnomAD annotated VCF files
    channels.existing_gnomad_vcfs = Channel
        .fromPath("${params.data}/families/*/vcfs/*.gnomad.vcf.gz")
        .map { vcf ->
            def fid = vcf.parent.parent.name  // Get family ID from path
            def tbi_path = "${vcf}.tbi"
            [fid, vcf, tbi_path]
        }
        .filter { fid, vcf, tbi_path -> 
            fid in analysis_plan.snvs_freq_annot.existing && new File(tbi_path).exists()
        }
        .map { fid, vcf, tbi_path ->
            [fid, vcf, file(tbi_path)]
        }

    // Create channel for existing rare VCF files (output from gnomAD filtering)
    channels.existing_filtered_vcfs = Channel
        .fromPath("${params.data}/families/*/vcfs/*.rare.vcf.gz")
        .map { vcf ->
            def fid = vcf.parent.parent.name  // Get family ID from path
            def tbi_path = "${vcf}.tbi"
            [fid, vcf, tbi_path]
        }
        .filter { fid, vcf, tbi_path -> 
            fid in analysis_plan.snvs_freq_filter.existing && new File(tbi_path).exists()
        }
        .map { fid, vcf, tbi_path ->
            [fid, vcf, file(tbi_path)]
        }

    // Create channel for existing common VCF files (output from gnomAD filtering)
    channels.existing_common_vcfs = Channel
        .fromPath("${params.data}/families/*/vcfs/*.common.vcf.gz")
        .map { vcf ->
            def fid = vcf.parent.parent.name  // Get family ID from path
            def tbi_path = "${vcf}.tbi"
            [fid, vcf, tbi_path]
        }
        .filter { fid, vcf, tbi_path -> 
            fid in analysis_plan.snvs_freq_filter.existing && new File(tbi_path).exists()
        }
        .map { fid, vcf, tbi_path ->
            [fid, vcf, file(tbi_path)]
        }

    // Create channel for existing common filtered BCF files (output from common filtering)
    channels.existing_common_filtered_bcfs = Channel
        .fromPath("${params.data}/families/*/vcfs/*.common_gt.bcf")
        .map { bcf ->
            def fid = bcf.parent.parent.name  // Get family ID from path
            def csi_path = "${bcf}.csi"
            [fid, bcf, csi_path]
        }
        .filter { fid, bcf, csi_path -> 
            fid in analysis_plan.snvs_common_filters.existing && new File(csi_path).exists()
        }
        .map { fid, bcf, csi_path ->
            [fid, bcf, file(csi_path)]
        }
    
    return channels
}
