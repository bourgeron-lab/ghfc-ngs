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
include { WOMBAT_COHORT } from './workflows/wombat_cohort'
include { WISECONDORX } from './workflows/wisecondorx'
include { WOMBAT } from './workflows/wombat'
include { EXTRACTOR } from './workflows/extractor'
include { ANCESTRY } from './workflows/ancestry'
include { ROH } from './workflows/roh'

/*
========================================================================================
    COHORT RUN STATE
========================================================================================
    Everything that writes .ghfc-ngs.state.json, the per-cohort record of when this cohort
    was last run, whether it finished, which pedigree and parameters produced its outputs,
    and how complete each step is. See COHORT_STATE.md for the schema.

    The mechanics (checksums, atomic write, history) live in lib/CohortState.groovy, which
    has no Nextflow dependency and is unit testable on its own. What a "step" is, and all
    logging, stays here.
*/

// A mutable holder rather than plain variables, and this is not a style choice: a bare
// assignment made inside `workflow {}` lands in a per-run binding that workflow.onComplete
// cannot see. Mutating a script-level container works from both scopes; reassigning it
// does not.
ghfc_run_state = [terminal_written: false, warned_no_cohort: false,
                  pedigree_file: null, pedigree_data: null, plan: null,
                  alignment_scan: null, stale_family_clean: null, progress: null]

// The cohort directory holds this cohort's params file, its pedigree and its outputs, so it
// is also where its state file belongs. Null when we cannot know it - params.cohort_name has
// no default in nextflow.config, and writing to .../cohorts/null/ would be worse than nothing.
def cohortStateDir() {
    if (!params.data || !params.cohort_name) return null
    return java.nio.file.Paths.get("${params.data}/cohorts/${params.cohort_name}")
}

// The -params-file path is not exposed on the workflow object in any Nextflow version, so the
// command line is the only place it can be recovered from.
def resolveParamsFilePath() {
    def matcher = (workflow.commandLine =~ /-params-file\s+(?:'([^']+)'|"([^"]+)"|(\S+))/)
    return matcher.find() ? (matcher.group(1) ?: matcher.group(2) ?: matcher.group(3)) : null
}

// workflow.start is an OffsetDateTime on Nextflow >= 22.04 and a java.util.Date before it.
// Both are written as ISO-8601 so a consumer never has to care which one produced the file.
def isoTimestamp(value) {
    if (value == null) return null
    if (value instanceof Date) return value.format("yyyy-MM-dd'T'HH:mm:ssXXX")
    return value.toString()
}

// A params file checksum cannot see --steps or --data given on the command line, which change
// the results without changing the file. Hashing the resolved map catches those too.
def effectiveParamsSha() {
    try {
        return CohortState.sha256OfString(groovy.json.JsonOutput.toJson(new TreeMap(params)))
    }
    catch (Exception e) {
        return null
    }
}

// Per-step completion, measured against the pedigree and never against the plan: `existing`
// and `needed` are not a partition of the cohort, because `needed` is gated on upstream
// prerequisites, so their sum understates the denominator.
//
// A null step is one that was not measured. That is not the same as one that is 0% done, and
// conflating the two would be the easiest way for this file to mislead someone.
def buildCompletionBlock(plan, Integer n_families, Integer n_individuals) {
    if (!plan) return null
    return [
        alignment:          CohortState.progress(plan.alignment.existing.size(), n_individuals),
        deepvariant_sample: CohortState.progress(plan.deepvariant_sample.existing.size(), n_individuals),
        deepvariant_family: CohortState.progress(plan.deepvariant_family.existing.size(), n_families),
        annotation:         CohortState.progress(plan.annotation.existing.size(), n_families),
        wombat:             CohortState.progress(plan.wombat.existing.size(), n_families),
        // Like ancestry below, planned only when requested - empty lists mean "unmeasured"
        wisecondorx:        ('wisecondorx' in pipeline_steps)
                                ? CohortState.progress(plan.wisecondorx.existing.size(), n_individuals)
                                : null,
        // The whole ancestry block of the plan is skipped when the step is not requested, so
        // its empty lists mean "unmeasured", not "nothing done"
        ancestry:           ('ancestry' in pipeline_steps)
                                ? CohortState.progress(plan.ancestry.existing.size(), n_families)
                                : null,
        // Same rule: planned only when requested
        roh:                ('roh' in pipeline_steps)
                                ? CohortState.progress(plan.roh.existing.size(), n_families)
                                : null,
        // One sentinel entry for the whole cohort, not a per-entity count
        snvs_cohort:        CohortState.progress(plan.snvs_cohort.existing.contains('cohort') ? 1 : 0, 1),
        wombat_cohort:      CohortState.progress(plan.wombat_cohort.existing.contains('cohort') ? 1 : 0, 1),
        // The plan does no existence check for extractor at all, so there is no honest number
        extractor:          null
    ]
}

// Which individuals have no CRAM, and whether anything on disk could align them.
//
// This is the question the run log could only ever half-answer: the plan's three CRAM-less
// buckets went to log.info/log.warn, and the input scan used to skip the stuck ones entirely.
// Every fact here was already gathered by createAnalysisPlan and resolveAlignmentInputs, so
// this does no filesystem work of its own and the completion handler can rebuild it against
// the re-scanned plan for free.
//
// Null when the scan never ran - unmeasured, not "none" - exactly as a null step in
// `completion` means unmeasured rather than 0%. A count of 0 means everyone has a CRAM.
def buildSamplesWithoutCramBlock(plan, scan, Map family_members, int limit = 200) {
    if (!plan || scan == null) return null

    def no_cram = plan.alignment.no_cram ?: []
    def with_gvcf = (plan.alignment.no_cram_with_gvcf ?: []) as Set
    def needed = (plan.alignment.needed ?: []) as Set

    def rows = no_cram.collect { barcode ->
        def source = scan.by_barcode?.get(barcode)
        def needs_alignment = barcode in needed
        [
            barcode        : barcode,
            family_id      : family_members?.get(barcode),
            has_gvcf       : barcode in with_gvcf,
            // Null rather than 'none' when the scan has no opinion at all, which happens when
            // no input source is configured - "we did not look" is not "we looked and found nothing"
            input_source   : source ? source.source : (scan.searched ? 'none' : null),
            input_path     : source?.path,
            // Whether this run schedules it, which is a different question from whether
            // anything could: a stuck member of an already-called family is scheduled for
            // nothing precisely because the stale family outputs make it look unnecessary
            needs_alignment: needs_alignment,
            // The one field worth grepping for: no gVCF, no CRAM, and nothing on disk that
            // could produce either. Deliberately not gated on needs_alignment - the samples
            // behind a stale family are scheduled for nothing at all, and they are exactly
            // the ones an operator has to go and find data for.
            blocked        : !(barcode in with_gvcf) && !(source?.usable)
        ]
    }

    def counts = [
        gvcf_only : rows.count { it.has_gvcf },
        will_align: rows.count { !it.has_gvcf && !it.blocked },
        no_input  : rows.count { it.blocked }
    ]

    // Actionable first, so the entries that matter are the ones that survive the cap. A
    // never-aligned cohort puts every individual in here, and at ~12 records per file that
    // is megabytes of JSON rewritten atomically on every run, for rows that are all alike.
    def sorted = rows.sort(false) { row ->
        [row.blocked ? 0 : (row.has_gvcf ? 2 : 1), row.barcode]
    }

    def block = [
        count    : rows.size(),
        blocked  : counts.no_input,
        by_status: counts,
        searched : scan.searched ?: []
    ]
    if (sorted.size() > limit) {
        block.truncated = true
        block.samples = sorted[0..<limit]
    } else {
        block.samples = sorted
    }
    return block
}

// Where this run was launched, and from which Slurm job if any. The job ID is what gives a
// 'running' record the liveness check it otherwise lacks: `squeue -j` either knows it or not.
def launchInfo() {
    def host = null
    try { host = java.net.InetAddress.localHost.hostName } catch (Exception ignored) { }
    return [
        dir:          workflow.launchDir?.toString(),
        slurm_job_id: System.getenv('SLURM_JOB_ID'),
        host:         host
    ]
}

// Start .ghfc-ngs.progress.json: Nextflow's per-process counters every 30 s, and the step
// completion re-measured against the pedigree every 10 min. See RunProgress and COHORT_STATE.md.
// Like the state file, never for a rehearsal, and never allowed to take the run down.
def startRunProgress(plan) {
    try {
        if (isRehearsalRun()) return
        def dir = cohortStateDir()
        if (!dir) return
        def pedigree_data = ghfc_run_state.pedigree_data
        def n_families    = pedigree_data?.families?.size()
        def n_individuals = pedigree_data?.individuals?.size()
        def launch = launchInfo()
        ghfc_run_state.progress = new RunProgress(
            dir: dir,
            header: [
                run_id:       workflow.sessionId?.toString(),
                run_name:     workflow.runName,
                cohort_name:  params.cohort_name,
                started_at:   isoTimestamp(workflow.start),
                launch_dir:   launch.dir,
                slurm_job_id: launch.slurm_job_id,
                host:         launch.host
            ],
            // The pedigree the run started with, not a re-read one: an edit made mid-run must
            // not move the denominator under a run that is still working to the old one
            rescan: {
                def p = createAnalysisPlan(pedigree_data.families, pedigree_data.individuals,
                                           pedigree_data.family_members, true)
                buildCompletionBlock(p, n_families, n_individuals)
            },
            warn: { String msg -> log.warn msg }
        ).start(buildCompletionBlock(plan, n_families, n_individuals))
    }
    catch (Throwable t) {
        log.warn "Could not start ${RunProgress.FILE_NAME}: ${t}"
    }
}

def buildStateRecord(Map opts) {
    // Fall back to whatever the run stashed, so most call sites only have to say what happened
    def pedigree_file = opts.containsKey('pedigree_file') ? opts.pedigree_file : ghfc_run_state.pedigree_file
    def pedigree_data = opts.containsKey('pedigree_data') ? opts.pedigree_data : ghfc_run_state.pedigree_data
    def plan          = opts.containsKey('plan')          ? opts.plan          : ghfc_run_state.plan
    def scan          = opts.containsKey('alignment_scan') ? opts.alignment_scan : ghfc_run_state.alignment_scan

    def n_families    = pedigree_data?.families?.size()
    def n_individuals = pedigree_data?.individuals?.size()
    def params_file   = resolveParamsFilePath()

    def record = [
        run_id:      workflow.sessionId?.toString(),
        status:      opts.status,
        started_at:  isoTimestamp(workflow.start),
        finished_at: opts.status == 'running' ? null : isoTimestamp(new Date()),
        cohort_name: params.cohort_name,
        pipeline: [
            version:    workflow.manifest?.version,
            revision:   workflow.revision,
            commit_id:  workflow.commitId,
            repository: workflow.repository
        ],
        pedigree: [
            path:        pedigree_file,
            sha256:      CohortState.sha256(pedigree_file as String),
            families:    n_families,
            individuals: n_individuals
        ],
        params_file: params_file ? [path: params_file, sha256: CohortState.sha256(params_file)] : null,
        // Where the run lives, so a tool can go from a cohort to its log and its Slurm job
        launch: launchInfo(),
        params_effective_sha256: effectiveParamsSha(),
        steps_requested: pipeline_steps,
        completion_measured: opts.measured,
        outputs_may_be_incomplete: opts.incomplete ? true : false,
        completion: buildCompletionBlock(plan, n_families, n_individuals),
        samples_without_cram: buildSamplesWithoutCramBlock(plan, scan, pedigree_data?.family_members),
        // Null when no clean was asked for, or when the run died before it could run - the
        // same "this was not measured" convention the blocks above use
        stale_family_clean: ghfc_run_state.stale_family_clean
    ]
    if (opts.reason) record.reason = opts.reason
    return record
}

// A rehearsal run: -stub-run fabricates real files with fake content, -preview builds the DAG
// without executing anything. Both still run the whole `workflow {}` body, so anything with a
// side effect out in the data tree - writing the state file, deleting stale outputs - has to
// ask. Assigned at script level so both the workflow and the completion handler can see it.
def isRehearsalRun() {
    if (workflow.stubRun) return true
    return workflow.hasProperty('preview')
        ? workflow.preview
        : (workflow.commandLine =~ /(^|\s)-preview(\s|\$)/).find()
}

def recordRunState(Map opts) {
    try {
        // A state file reporting a stub run's fabricated outputs as complete would be worse
        // than no file at all, and a preview did not do the work it would be describing.
        if (isRehearsalRun()) return

        if (opts.status != 'running') {
            // More than one exit path can reach here; the first terminal record is the true one
            if (ghfc_run_state.terminal_written) return
            ghfc_run_state.terminal_written = true
        }

        def dir = cohortStateDir()
        if (!dir) {
            // Once per run: this is called from the running marker, from every exit site and
            // from the completion handler, and repeating it would just be noise
            if (!ghfc_run_state.warned_no_cohort) {
                ghfc_run_state.warned_no_cohort = true
                log.warn "cohort_name is not set - skipping ${CohortState.FILE_NAME}"
            }
            return
        }

        CohortState.write(dir, CohortState.merge(CohortState.read(dir), buildStateRecord(opts), CohortState.HISTORY_LIMIT))
    }
    catch (Throwable t) {
        // Provenance bookkeeping must never take down a run that may have been going for days
        log.warn "Could not write ${CohortState.FILE_NAME}: ${t}"
    }
}

// Record a failure that is about to exit. workflow.onComplete never fires for these: Nextflow's
// `exit` is a System.exit with no shutdown hook, so this is the only chance to leave a trace.
def recordFailedRun(String reason) {
    recordRunState(status: 'failed', measured: 'before', reason: reason)
}

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Nextflow types params.steps from wherever it came: a YAML list arrives as a List, while
// --steps "a,b" on the command line arrives as a String. Every use below tests membership or
// joins, and `'alignment' in "alignment,wombat"` is quietly false rather than an error, so
// both forms are reconciled into one list here and nothing downstream reads params.steps.
// Assigned at script level without def, so the onComplete handler and the functions it calls
// can see it - the same reason ghfc_run_state is written that way.
def normaliseSteps(value) {
    if (value == null) return []
    def items
    if (value instanceof CharSequence) {
        items = value.toString().split(',') as List
    } else if (value instanceof Collection) {
        items = value as List
    } else {
        items = [value]
    }
    return items.collect { it?.toString()?.trim() }.findAll { it }
}

pipeline_steps = normaliseSteps(params.steps)

// Groovy truth makes the non-empty string "false" true, and `--flag false` on the command line
// arrives as exactly that string - measured on 26.04, for both `--flag false` and `--flag=false`.
// Only a YAML params file yields a real Boolean. So a bare truth test on a flag reads the
// operator's "off" as "on", which for a flag that deletes call sets is not survivable.
def asBoolean(value) {
    if (value == null) return false
    if (value instanceof Boolean) return value
    return !(value.toString().trim().toLowerCase() in ['', 'false', 'no', '0', 'null'])
}

// Validate input parameters
if (!params.data) {
    exit 1, "ERROR: --data parameter is required"
}

// Every cohort-level output path, the default pedigree location and the state file itself are
// built from cohort_name. Without it the run writes .../cohorts/null/null.* and records nothing
// about itself. Checked here rather than at first use, and with a bare exit like data above,
// because recordFailedRun has nowhere to write a state file until both of them are known.
if (!params.cohort_name) {
    exit 1, "ERROR: cohort_name parameter is required"
}

if (pipeline_steps.isEmpty()) {
    recordFailedRun("no steps requested")
    exit 1, "ERROR: --steps parameter is required. Available steps: alignment, deepvariant_sample, deepvariant_family, annotation, wombat, wombat_cohort, snvs_cohort, wisecondorx, extractor, ancestry, roh"
}

// Validate steps
def valid_steps = ['alignment', 'deepvariant_sample', 'deepvariant_family', 'annotation', 'wombat', 'wombat_cohort', 'snvs_cohort', 'wisecondorx', 'extractor', 'ancestry', 'roh']
def invalid_steps = pipeline_steps - valid_steps
if (invalid_steps) {
    recordFailedRun("invalid steps requested: ${invalid_steps.join(', ')}")
    exit 1, "ERROR: Invalid steps specified: ${invalid_steps.join(', ')}. Valid steps are: ${valid_steps.join(', ')}"
}

// snvs_cohort used to merge the wombat tables too. A parameters file written then still lists
// only snvs_cohort, and would stop refreshing its cohort wombat tables without a word.
if ('snvs_cohort' in pipeline_steps && !('wombat_cohort' in pipeline_steps) && params.wombat_config_list) {
    log.warn "snvs_cohort no longer merges the wombat tables into cohort tables: add wombat_cohort to steps for that"
}

// The ancestry step reads its panel and weights from paths that have no sensible
// default, and every one of its processes would fail on an empty string.
if ('ancestry' in pipeline_steps) {
    def missing_ancestry_params = ['ancestry_reference', 'ancestry_catalog', 'ancestry_panel_name']
        .findAll { key -> !params[key] }
    if (missing_ancestry_params) {
        recordFailedRun("ancestry step is missing ${missing_ancestry_params.join(', ')}")
        exit 1, "ERROR: the 'ancestry' step requires ${missing_ancestry_params.join(', ')} to be set"
    }
}

// roh reads the ancestry step's family panel and labels, named by the panel label, and calls
// against the same bundle - which must carry the roh/ tier. The PGS catalog is not needed.
if ('roh' in pipeline_steps) {
    def missing_roh_params = ['ancestry_reference', 'ancestry_panel_name']
        .findAll { key -> !params[key] }
    if (missing_roh_params) {
        recordFailedRun("roh step is missing ${missing_roh_params.join(', ')}")
        exit 1, "ERROR: the 'roh' step requires ${missing_roh_params.join(', ')} to be set"
    }
    // Otherwise every family's ROH_CALL would fail on it, one task at a time
    def roh_sites = "${params.ancestry_reference}/roh/sites.tsv.gz"
    if (!new File(roh_sites).exists()) {
        recordFailedRun("ancestry_reference has no roh/ tier")
        exit 1, "ERROR: the 'roh' step needs the reference bundle's roh/ tier, and ${roh_sites} does not exist. Add it with ancestry-pgs bundle/add_roh_tier.py (bundle >= 1.1.0)"
    }
}

// The cohort's `resources:` overrides are only read by sized() in nextflow.config, task by task,
// so a typo there would surface hours into the run, or - for a key that matches nothing - never.
// containsKey first: reading an unset params key, even with get(), logs a warning every run.
def resource_check = ResourceOverrides.check(params.containsKey('resources') ? params.resources : null,
    ResourceOverrides.processNames(projectDir.toFile()),
    [memory: params.max_memory, cpus: params.max_cpus, time: params.max_time])
resource_check.warnings.each { log.warn(it) }
if (resource_check.errors) {
    recordFailedRun("invalid resources: ${resource_check.errors.join('; ')}")
    exit 1, "ERROR: invalid resources in the parameters file:\n  ${resource_check.errors.join('\n  ')}"
}

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    
    // Read and validate pedigree file. The pedigree lives in the cohort directory, next to the
    // parameters file and the cohort's outputs, so that convention is the default and only a
    // pedigree kept somewhere else needs the parameter set at all.
    def pedigree_file = params.pedigree ?: "${params.data}/cohorts/${params.cohort_name}/${params.cohort_name}.pedigree.tsv"
    // Stashed before the existence check so a failure record can still name the path it wanted
    ghfc_run_state.pedigree_file = pedigree_file

    if (!new File(pedigree_file).exists()) {
        recordFailedRun("pedigree file not found: ${pedigree_file}")
        // Which of the two it is changes what the operator has to fix, so say so
        def hint = params.pedigree
            ? " (from the pedigree parameter)"
            : " (default location for cohort '${params.cohort_name}' - set pedigree: to read it from elsewhere)"
        exit 1, "ERROR: Pedigree file not found: ${pedigree_file}${hint}"
    }
    
    log.info """
    ========================================================================================
                            GHFC WGS Family-based Pipeline
    ========================================================================================
    Pedigree file    : ${pedigree_file}
    Data directory   : ${params.data}
    Steps to run     : ${pipeline_steps.join(', ')}
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
    ghfc_run_state.pedigree_data = pedigree_data

    log.info "Found ${families.size()} families with ${individuals.size()} individuals total"
    
    // --clean-stale-families reaches Nextflow as clean_stale_families; cleanStaleFamilies is
    // the camelCase spelling Nextflow produces when the hyphenated flag bypasses the wrapper
    def clean_requested = asBoolean(params.clean_stale_families) || asBoolean(params.cleanStaleFamilies)

    // Check existing files and determine what needs to be done.
    // Quiet on this pass when a clean is coming: its per-sample messages describe a tree that
    // is about to change, and the rebuild below re-emits them against the one the run uses.
    def analysis_plan = createAnalysisPlan(families, individuals, family_members, clean_requested)

    // Resolve the input files behind every individual with no CRAM, once. The guard below and
    // the channels built afterwards need it, the run-state report is built from it, and the
    // clean cannot judge a drifted family without it.
    //
    // Before the clean, and reused after it. That is sound only because the CRAM-less set
    // this scan is keyed on cannot change when family VCFs are deleted - so the verdicts stay
    // valid for the rebuilt plan, and the directory globs happen once.
    def resolved_inputs = resolveAlignmentInputs(analysis_plan)
    ghfc_run_state.alignment_scan = resolved_inputs

    if (clean_requested) {
        if (!params.vep_config_name) {
            // Every annotation and wombat path is named after it. Unset, they would all be
            // built as '...rare.null...', deleteIfExists would quietly match nothing, and the
            // clean would report success having removed only the deepvariant_family files.
            recordFailedRun("clean_stale_families requires vep_config_name")
            exit 1, "ERROR: clean_stale_families needs vep_config_name to know which annotation and wombat outputs belong to a family - set it in the parameters file"
        }
        if (analysis_plan.pedigree_drift) {
            ghfc_run_state.stale_family_clean = cleanStaleFamilies(analysis_plan, resolved_inputs, family_members)
            // A rehearsal deleted nothing, so re-planning would have the run proceed as
            // though the outputs were gone when every one of them is still there
            if (ghfc_run_state.stale_family_clean.families_cleaned &&
                !ghfc_run_state.stale_family_clean.rehearsal) {
                analysis_plan = createAnalysisPlan(families, individuals, family_members)
            }
        } else {
            log.info "clean_stale_families is set, but no family has stale outputs - nothing to clean"
        }
    }
    ghfc_run_state.plan = analysis_plan

    // Display analysis summary. After a clean this describes the tree the run will actually
    // work on, and any family still listed as drifted is one the clean refused.
    displayAnalysisSummary(analysis_plan, resolved_inputs, ghfc_run_state.stale_family_clean)

    // Stale family outputs are wrong results rather than missing ones, so the run can be
    // stopped on them. Warning is the default deliberately: a cohort whose pedigree has
    // drifted should not silently become unrunnable without the operator asking for that.
    // Evaluated after the clean, so the two together mean "fix what you can, stop on the rest".
    if (asBoolean(params.pedigree_strict) && analysis_plan.pedigree_drift) {
        def n = analysis_plan.pedigree_drift.size()
        recordFailedRun("${n} ${n == 1 ? 'family has' : 'families have'} stale outputs and pedigree_strict is set")
        exit 1, "ERROR: ${n} ${n == 1 ? 'family has' : 'families have'} stale outputs and pedigree_strict is set - see the STALE FAMILY OUTPUTS summary above"
    }

    // Validate that required steps are available
    validateStepsAvailability(analysis_plan, resolved_inputs, family_members)

    // Validation passed and real work is about to start. Nothing fires on a SLURM walltime
    // kill or a Ctrl-C, so this marker is what makes such a death visible afterwards: the
    // record stays 'running' forever, and the next run turns it into 'interrupted'.
    recordRunState(status: 'running', measured: 'before')
    startRunProgress(analysis_plan)
    
    // Create channels for different steps
    def channels = createChannels(analysis_plan, resolved_inputs)
    
    // Run alignment if needed and allowed
    // Also entered when CRAMs exist but their coverage bedgraphs are missing, in which case only
    // MOSDEPTH/TABIX_INDEX run - no realignment is triggered.
    if ((analysis_plan.alignment.needed.size() > 0 || analysis_plan.alignment.need_bedgraph.size() > 0) && 'alignment' in pipeline_steps) {
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
    if (analysis_plan.deepvariant_sample.needed.size() > 0 && 'deepvariant_sample' in pipeline_steps) {
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
    if (analysis_plan.deepvariant_family.needed.size() > 0 && 'deepvariant_family' in pipeline_steps) {
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
    if (analysis_plan.annotation.needed.size() > 0 && 'annotation' in pipeline_steps) {
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
    if (analysis_plan.wombat.needed.size() > 0 && 'wombat' in pipeline_steps) {
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
    if (analysis_plan.snvs_cohort.needed.size() > 0 && 'snvs_cohort' in pipeline_steps) {
        log.info "Running cohort: common variants BCF merge..."
        
        // Get all available common filtered BCFs (existing + newly created)
        if (analysis_plan.annotation.needed.size() > 0 && 'annotation' in pipeline_steps) {
            // Mix existing BCFs with newly created ones
            all_available_common_bcfs = channels.existing_common_filtered_bcfs.mix(annotation_common_bcfs_output)
        } else {
            // Use only existing BCFs if no new ones were created
            all_available_common_bcfs = channels.existing_common_filtered_bcfs
        }
        
        SNVS_COHORT(all_available_common_bcfs,
                    analysis_plan.snvs_cohort.need_bcf_merge,
                    families,
                    analysis_plan.annotation.needed)
    }
    
    // Run the cohort wombat merges if needed and allowed
    if (analysis_plan.wombat_cohort.needed.size() > 0 && 'wombat_cohort' in pipeline_steps) {
        def wombat_merge_count = analysis_plan.wombat_cohort.need_merges.count { k, v -> v == true }
        log.info "Running cohort: wombat merge (${wombat_merge_count} configs)..."
        
        // Get all available Wombat files (existing + newly created)
        all_available_wombat_files = channels.existing_wombat_files.mix(wombat_output ?: Channel.empty())
        
        WOMBAT_COHORT(all_available_wombat_files,
                      analysis_plan.wombat_cohort.need_merges,
                      families,
                      analysis_plan.wombat.needed)
    }
    
    // Run WisecondorX predict if needed and allowed
    if ((analysis_plan.wisecondorx.needed.size() > 0 || 
         analysis_plan.wisecondorx.need_family_merge.any { k, v -> v == true } || 
         analysis_plan.wisecondorx.need_family_annotate.any { k, v -> v == true } || 
         analysis_plan.wisecondorx.need_cohort_merge) && 'wisecondorx' in pipeline_steps) {
        
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
    if (analysis_plan.extractor.tsv_count > 0 && 'extractor' in pipeline_steps) {
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
        
        // Create channel for wombat bcf2parquet outputs (for PROCESS_FAM_TSV). The planner
        // records, for every family in the pedigree, whether this file is missing.
        wombat_parquets_for_extractor = Channel.fromList(analysis_plan.wombat.need_bcf2parquet
            .findAll { _fid, missing -> missing == false }
            .collect { fid, _missing ->
                tuple(fid, file("${Sharding.getFamilyDir(params.data, fid)}/wombat/${fid}.rare.${params.vep_config_name}.annotated.parquet"))
            })
        
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
    ancestry_family_panels = Channel.empty()
    ancestry_family_labels = Channel.empty()
    if ((analysis_plan.ancestry.needed.size() > 0 ||
         analysis_plan.ancestry.need_cohort_merge) && 'ancestry' in pipeline_steps) {

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
        ancestry_family_panels = ANCESTRY.out.family_panel_gt
        ancestry_family_labels = ANCESTRY.out.family_ancestry
    }

    // Run ROH / inbreeding if needed and allowed.
    //
    // Calls on the ancestry step's family panel, for the same reason ancestry reads gVCFs: ROH
    // calling reads homozygosity from hom-ref sites, which common_gt.bcf does not carry. A
    // family the ancestry step rebuilds in this run takes that run's panel and label; every
    // other family reads its own from disk.
    if ((analysis_plan.roh.needed.size() > 0 ||
         analysis_plan.roh.need_cohort_merge) && 'roh' in pipeline_steps) {

        def roh_tasks = []
        if (analysis_plan.roh.needed.size() > 0) {
            roh_tasks.add("calling for ${analysis_plan.roh.needed.size()} families")
        }
        if (analysis_plan.roh.need_cohort_merge) {
            roh_tasks.add("cohort table merge")
        }
        log.info "Running ROH/inbreeding: ${roh_tasks.join(', ')}..."

        def roh_panel_name = params.ancestry_panel_name
        def roh_from_ancestry = analysis_plan.roh.from_ancestry
        roh_disk_inputs = Channel.fromList(
            analysis_plan.roh.needed
                .findAll { fid -> !(fid in roh_from_ancestry) }
                .collect { fid ->
                    def dir = "${Sharding.getFamilyDir(params.data, fid)}/ancestry"
                    def bcf = file("${dir}/${fid}.panel_gt.${roh_panel_name}.bcf")
                    tuple(fid, bcf, file("${bcf}.csi"), file("${dir}/${fid}.${roh_panel_name}.ancestry.tsv"))
                })
        roh_ancestry_inputs = ancestry_family_panels
            .filter { fid, _bcf, _csi -> fid in roh_from_ancestry }
            .join(ancestry_family_labels)

        ROH(
            roh_disk_inputs.mix(roh_ancestry_inputs),
            pedigree_file,
            pedigree_data.families,
            analysis_plan.roh.needed,
            analysis_plan.roh.need_cohort_merge,
            rohGeneBeds()
        )
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

// `quiet` suppresses the operator-facing messages only. The completion handler re-runs this
// against the finished output tree, and every one of those lines has already been printed once
// at the start of the run.
def parsePedigreeFile(pedigree_file, quiet = false) {
    def families = [] as Set
    def individuals = [] as Set
    def family_members = [:]
    
    file(pedigree_file).readLines().eachWithIndex { line, index ->
        // Skip comments and empty lines
        if (line.startsWith('#') || line.trim().isEmpty()) return
        
        def cols = line.split('\t')
        
        // Skip header row if first column is "FID"
        if (index == 0 && cols[0] == 'FID') {
            if (!quiet) log.info "Skipping pedigree header row"
            return
        }
        
        if (cols.size() < 6) {
            recordFailedRun("pedigree row with ${cols.size()} columns instead of 6")
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

// family ID -> its members' barcodes, in pedigree order. Built once, because asking
// family_members "who is in this family?" once per family is a scan of every individual per
// family - invisible on a few hundred families, hours of planning on a cohort the size of SPARK.
def membersByFamily(family_members) {
    def members_by_family = [:]
    family_members.each { barcode, fid ->
        def members = members_by_family[fid]
        if (members == null) {
            members = []
            members_by_family[fid] = members
        }
        members.add(barcode)
    }
    return members_by_family
}

// An ordered set rather than a list for every bucket of the plan. Each bucket is asked
// "is this ID in you?" once per individual or family, both here and in the channel filters,
// and a list answers that by scanning itself. Insertion order is kept, so every log line and
// report that joins a bucket reads exactly as it did when these were lists.
def planBucket() {
    return new LinkedHashSet()
}

// The gene BEDs the ROH step intersects with, from the annotation's gencode settings, as
// [symbol BED, ID BED]. Both empty when gencode is not configured, and then no gene table is
// written or expected. containsKey first: reading an unset params key logs a warning.
def rohGeneBeds() {
    def gencode = params.containsKey('annotation_gencode') ? params.annotation_gencode : null
    def path = params.containsKey('annotation_annotation_path') ? params.annotation_annotation_path : null
    if (!gencode || !path) return ['', '']
    return ["${path}/gencode/${gencode}.symbol.genes.bed.gz".toString(),
            "${path}/gencode/${gencode}.ensg.genes.bed.gz".toString()]
}

// The per-family ROH tables whose presence means a family is done. The kinship table is left
// out: a singleton has none.
def rohFamilyTableKinds() {
    return ['froh', 'roh'] + (rohGeneBeds()[0] ? ['roh_genes'] : [])
}

def createAnalysisPlan(families, individuals, family_members, quiet = false) {
    def plan = [
        // no_cram is every individual with no usable CRAM, whatever the reason and whether or
        // not anything will be scheduled for it. It is the set the input scan probes, which is
        // what lets the run report on samples that nothing would otherwise look at.
        alignment: [needed: planBucket(), existing: planBucket(), need_bedgraph: planBucket(), no_cram: planBucket(), no_cram_with_gvcf: planBucket()],
        deepvariant_sample: [needed: planBucket(), existing: planBucket()],
        deepvariant_family: [needed: planBucket(), existing: planBucket()],
        annotation: [needed: planBucket(), existing: planBucket()],
        snvs_cohort: [needed: planBucket(), existing: planBucket(), need_bcf_merge: false],
        wombat_cohort: [needed: planBucket(), existing: planBucket(), need_merges: [:]],
        wisecondorx: [needed: planBucket(), existing: planBucket(), need_npz: planBucket(), need_predict: planBucket(), need_family_merge: [:], need_family_annotate: [:], need_cohort_merge: false],
        wombat: [needed: planBucket(), existing: planBucket(), need_bcf2parquet: [:]],
        extractor: [tsv_count: 0, families: [] as Set, samples: [] as Set],
        ancestry: [needed: planBucket(), existing: planBucket(), need_extract: planBucket(), need_family_merge: [:], need_family_score: [:], need_cohort_merge: false],
        // from_ancestry: families whose panel or labels the ancestry step rebuilds in this run, so
        // their ROH inputs come from that run rather than disk. blocked: families with no ancestry
        // outputs on disk and no ancestry step to make them.
        roh: [needed: planBucket(), existing: planBucket(), from_ancestry: planBucket(), blocked: planBucket(), need_cohort_merge: false],
        pedigree_drift: [:]
    ]

    def members_by_family = membersByFamily(family_members)

    // The drift check, the deepvariant_sample check and the CRAM-less check below each ask
    // whether a sample's gVCF exists. On network storage every one of those is a round trip,
    // and at a cohort of 10^5 samples asking three times is two passes of pure waste - so the
    // first answer is kept for the rest of this plan.
    def gvcf_exists = [:]
    def hasGvcf = { barcode ->
        def known = gvcf_exists[barcode]
        if (known == null) {
            known = new File("${Sharding.getSampleDir(params.data, barcode)}/deepvariant/${barcode}.g.vcf.gz").exists()
            gvcf_exists[barcode] = known
        }
        return known
    }

    // Unset means the default, not "off": asBoolean reads a missing value as false
    def require_vaf_bedgraph = params.vaf_bedgraph == null ? true : asBoolean(params.vaf_bedgraph)
    
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
        if (new File(pedigree_path).exists() && !pedigree_ok && !quiet) {
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
            def members = members_by_family[fid] ?: []
            def without_gvcf = members.findAll { barcode -> !hasGvcf(barcode) }.sort()
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

    // Check the existing cohort common BCF
    def cohort_bcf_path = "${params.data}/cohorts/${params.cohort_name}/vcfs/${params.cohort_name}.common_gt.bcf"
    def cohort_csi_path = "${params.data}/cohorts/${params.cohort_name}/vcfs/${params.cohort_name}.common_gt.bcf.csi"
    
    def cohort_bcf_exists = new File(cohort_bcf_path).exists() && new File(cohort_csi_path).exists()
    
    if (cohort_bcf_exists) {
        plan.snvs_cohort.existing.add('cohort')  // Single cohort entry
    } else {
        // Need cohort merge if any families have common filtered BCFs or will create them
        def families_with_common_bcfs = plan.annotation.existing + plan.annotation.needed
        if (families_with_common_bcfs.size() > 0) {
            plan.snvs_cohort.needed.add('cohort')  // Single cohort entry
            plan.snvs_cohort.need_bcf_merge = true
        }
    }
    
    // Check existing ancestry outputs (per-sample panel genotypes, family panel
    // genotypes, family tables, cohort tables - all part of ancestry)
    //
    // The panel label is part of every file name on purpose: the depth/quality
    // thresholds and the reference bundle are not recorded anywhere this check can
    // see, so bumping ancestry_panel_name is what invalidates the old extractions.
    if ('ancestry' in pipeline_steps) {
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
            def members = members_by_family[fid] ?: []
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

    // Check existing ROH outputs. They are computed from the ancestry step's family panel and
    // labels, so a family the ancestry step rebuilds in this run is re-called too: its old ROH
    // tables describe a panel or an ancestry label that is about to change. Named by the
    // ancestry panel label for the same reason, which also pins the bundle and its roh/ tier.
    if ('roh' in pipeline_steps) {
        def panel_name = params.ancestry_panel_name
        def kinds = rohFamilyTableKinds()

        families.each { fid ->
            def fam_dir = Sharding.getFamilyDir(params.data, fid)
            def ancestry_redo = plan.ancestry.need_family_merge[fid] == true ||
                                plan.ancestry.need_family_score[fid] == true
            def outputs_exist = kinds.every { kind ->
                new File("${fam_dir}/roh/${fid}.${panel_name}.${kind}.tsv").exists()
            }
            if (outputs_exist && !ancestry_redo) {
                plan.roh.existing.add(fid)
                return
            }
            plan.roh.needed.add(fid)
            if (ancestry_redo) {
                plan.roh.from_ancestry.add(fid)
                return
            }
            def panel_bcf = "${fam_dir}/ancestry/${fid}.panel_gt.${panel_name}.bcf"
            def inputs_exist = new File(panel_bcf).exists() && new File("${panel_bcf}.csi").exists() &&
                new File("${fam_dir}/ancestry/${fid}.${panel_name}.ancestry.tsv").exists()
            if (!inputs_exist) {
                plan.roh.blocked.add(fid)
            }
        }

        def cohort_tables_exist = (kinds + ['roh_windows']).every { kind ->
            new File("${params.data}/cohorts/${params.cohort_name}/roh/${params.cohort_name}.${panel_name}.${kind}.tsv").exists()
        }
        plan.roh.need_cohort_merge = !cohort_tables_exist || plan.roh.needed.size() > 0
    }

    // Check existing individual gVCF files and VAF bedgraphs (both outputs of deepvariant_sample)
    individuals.each { barcode ->
        def smp_dir = Sharding.getSampleDir(params.data, barcode)
        def gvcf_tbi_path = "${smp_dir}/deepvariant/${barcode}.g.vcf.gz.tbi"
        def vaf_bedgraph_path = "${smp_dir}/sequences/${barcode}.vaf.bedgraph.gz"
        def vaf_bedgraph_tbi_path = "${smp_dir}/sequences/${barcode}.vaf.bedgraph.gz.tbi"
        
        // Without vaf_bedgraph, the gVCF alone is the sample's output: a cohort imported as gVCFs
        // has no VCF to build a bedgraph from, so requiring one would send it back to alignment
        def has_vaf_bedgraph = !require_vaf_bedgraph ||
            (new File(vaf_bedgraph_path).exists() && new File(vaf_bedgraph_tbi_path).exists())

        if (hasGvcf(barcode) && new File(gvcf_tbi_path).exists() && has_vaf_bedgraph) {
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
    //
    // With no alignment input configured at all, a CRAM-less sample with a gVCF is the cohort's
    // normal state rather than something to look into - an imported gVCF cohort is nothing
    // but - so those are counted into one line instead of one per sample.
    def expects_crams = params.fastq_pattern || params.old_cram_37 || params.old_cram_38
    def gvcf_only_count = 0
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
            plan.alignment.no_cram.add(barcode)
            // Hoisted out of the two branches below, which each needed it anyway, so the
            // three CRAM-less buckets can be told apart later without a second stat
            def has_gvcf = hasGvcf(barcode)
            if (has_gvcf) {
                plan.alignment.no_cram_with_gvcf.add(barcode)
            }

            // Only need alignment if the individual needs DeepVariant
            if (barcode in plan.deepvariant_sample.needed) {
                plan.alignment.needed.add(barcode)
            } else if (has_gvcf) {
                // Benign: nothing downstream reads the CRAM. Family calling consumes the
                // gVCF and the ancestry panel extraction reads it directly, so a missing
                // CRAM here costs nothing but the ability to regenerate coverage bedgraphs.
                def bedgraph_note = has_bedgraph ? "" : " (its coverage bedgraph is also absent and cannot be regenerated without the CRAM)"
                if (!expects_crams) {
                    gvcf_only_count++
                } else if (!quiet) {
                    log.info "Individual ${barcode} has no ${params.ref_name} CRAM, but its gVCF is present, so nothing downstream needs it${bedgraph_note}"
                }
            } else {
                // Genuinely stuck: no CRAM, no gVCF, and the family is already called, so
                // the plan will never schedule anything for this individual and the family's
                // outputs cannot contain it. See the STALE FAMILY OUTPUTS summary.
                if (!quiet) log.warn "Individual ${barcode} has neither a ${params.ref_name} CRAM nor a gVCF, and family ${family_members[barcode]} is already called - nothing will be scheduled for it and the family's outputs cannot include it. Provide input data for it, or remove it from the pedigree"
            }
        }
    }
    if (gvcf_only_count && !quiet) {
        log.info "${gvcf_only_count} individuals have a gVCF and no ${params.ref_name} CRAM, and no alignment input is configured - nothing downstream needs their CRAMs, but their coverage bedgraphs cannot be generated"
    }
    
    // Only when the step is requested, like ancestry above. Three checks per individual and as
    // many per family is a large share of a planning pass that nothing reads otherwise: every
    // consumer of plan.wisecondorx - the WISECONDORX call, the summary and the completion block -
    // is itself gated on the step.
    if ('wisecondorx' in pipeline_steps) {
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
                def family_has_aberrations = (members_by_family[fid] ?: []).any { barcode ->
                    barcode in plan.wisecondorx.existing || barcode in plan.wisecondorx.needed
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

    // Check the existing cohort wombat tables, one per config. Like the cohort BCF, a table that
    // exists is not rebuilt when families are added - delete it to have it merged again.
    (params.wombat_config_list ?: []).each { config_file ->
        def config_name = config_file.replaceAll(/\.ya?ml$/, '')
        def cohort_wombat_path = "${params.data}/cohorts/${params.cohort_name}/wombat/${params.cohort_name}.rare.${params.vep_config_name}.annotated.${config_name}.results.tsv"
        plan.wombat_cohort.need_merges[config_name] = !new File(cohort_wombat_path).exists()
    }
    if (plan.wombat_cohort.need_merges) {
        if (!plan.wombat_cohort.need_merges.any { k, v -> v == true }) {
            plan.wombat_cohort.existing.add('cohort')  // Single cohort entry
        } else if (plan.wombat.existing || plan.wombat.needed) {
            plan.wombat_cohort.needed.add('cohort')  // Single cohort entry
        }
    }
    
    return plan
}

def displayAnalysisSummary(analysis_plan, resolved_inputs = null, clean_report = null) {
    // Every CRAM-less individual, split the way an operator has to act on them: the ones with
    // a gVCF need nothing, and the rest are only fine if an input actually resolves for them
    def no_cram = analysis_plan.alignment.no_cram ?: []
    def no_cram_gvcf = (analysis_plan.alignment.no_cram_with_gvcf ?: []).size()
    def no_cram_blocked = resolved_inputs == null ? null : no_cram.count { barcode ->
        !(barcode in analysis_plan.alignment.no_cram_with_gvcf) &&
        !(resolved_inputs.by_barcode?.get(barcode)?.usable)
    }

    log.info """
    ========================================================================================
                                    ANALYSIS SUMMARY
    ========================================================================================
    ALIGNMENT: ${analysis_plan.alignment.existing.size()} individuals done, ${analysis_plan.alignment.needed.size()} to align, ${analysis_plan.alignment.need_bedgraph.size()} needing bedgraph only
    SAMPLES WITHOUT CRAM: ${no_cram.size()}${no_cram ? " (${no_cram_gvcf} with a gVCF so nothing needs one, ${no_cram_blocked == null ? 'inputs not scanned' : "${no_cram_blocked} with no usable input"})" : ''}
    == SNVs/INDELs Calling ==
    DEEPVARIANT_SAMPLE: ${analysis_plan.deepvariant_sample.existing.size()} individuals done and ${analysis_plan.deepvariant_sample.needed.size()} to do
    DEEPVARIANT_FAMILY: ${analysis_plan.deepvariant_family.existing.size()} families done and ${analysis_plan.deepvariant_family.needed.size()} to do
    ANNOTATION: ${analysis_plan.annotation.existing.size()} families done and ${analysis_plan.annotation.needed.size()} to do
    WOMBAT: ${analysis_plan.wombat.existing.size()} families done and ${analysis_plan.wombat.needed.size()} to do
    WOMBAT_COHORT: ${'wombat_cohort' in pipeline_steps ? "cohort wombat merges due: ${analysis_plan.wombat_cohort.need_merges.count { _k, v -> v == true }}" : 'Skipped (step not requested)'}
    == Common Variants ==
    SNVS_COHORT: ${'snvs_cohort' in pipeline_steps ? "common variants cohort bcf merge: ${analysis_plan.snvs_cohort.need_bcf_merge ? 'Yes' : 'No'}" : 'Skipped (step not requested)'}
    == SVs Calling ==
    WISECONDORX PREDICT: ${'wisecondorx' in pipeline_steps ? "${analysis_plan.wisecondorx.existing.size()} individuals done and ${analysis_plan.wisecondorx.needed.size()} to do" : 'Skipped (step not requested)'}
    == Ancestry / PGS ==
    ANCESTRY: ${'ancestry' in pipeline_steps ? "${analysis_plan.ancestry.existing.size()} families done and ${analysis_plan.ancestry.needed.size()} to do (${analysis_plan.ancestry.need_extract.size()} samples needing panel extraction)" : 'Skipped (step not requested)'}
    ROH: ${'roh' in pipeline_steps ? "${analysis_plan.roh.existing.size()} families done and ${analysis_plan.roh.needed.size()} to do (${analysis_plan.roh.from_ancestry.size()} waiting on this run's ancestry)" : 'Skipped (step not requested)'}
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
    if (analysis_plan.snvs_cohort.need_bcf_merge && 'snvs_cohort' in pipeline_steps) {
        log.info "Cohort needing common variant merge: Yes"
    }
    def due_wombat_merges = analysis_plan.wombat_cohort.need_merges.findAll { _k, v -> v == true }.keySet()
    if (due_wombat_merges && 'wombat_cohort' in pipeline_steps) {
        log.info "Cohort wombat merges needed for: ${due_wombat_merges.join(', ')}"
    }
    if (analysis_plan.ancestry.needed) {
        log.info "Families needing ancestry/PGS: ${analysis_plan.ancestry.needed.join(', ')}"
    }
    if (analysis_plan.roh.needed) {
        log.info "Families needing ROH/inbreeding: ${analysis_plan.roh.needed.join(', ')}"
    }
    if (analysis_plan.pedigree_drift) {
        def n_drift = analysis_plan.pedigree_drift.size()
        def drift_noun = n_drift == 1 ? "FAMILY" : "FAMILIES"
        // After a real clean, every family still listed here is one the clean refused, so
        // say why rather than repeating advice the operator has already acted on. After a
        // rehearsal nothing was deleted, so the cleanable ones are still listed too and have
        // to be told apart from the refusals.
        def skip_reasons = [:]
        (clean_report?.families_skipped ?: []).each { entry -> skip_reasons[entry.family_id] = entry.reason }
        def would_clean = (clean_report?.rehearsal ? clean_report.families_cleaned : []).collect { it.family_id } as Set

        def drift_lines = analysis_plan.pedigree_drift.collect { fid, drift ->
            def line = "${fid}: ${drift.missing.size()} of ${drift.members} members have no gVCF (${drift.missing.join(', ')})"
            if (skip_reasons[fid]) return "${line}\n        not cleaned: ${skip_reasons[fid]}"
            if (fid in would_clean) return "${line}\n        would be cleaned by this flag"
            return line
        }

        def remedy = clean_report?.rehearsal
            ? """This was a rehearsal: nothing was deleted. The families marked above would be cleaned;
    run the same command without --dry-run to do it. The rest were refused for the reason
    given above them - supply the input data, index the CRAM, or add the missing step to
    steps: - or remove the members that have no data from the pedigree."""
            : clean_report != null
            ? """Each of these was left alone for the reason given above it. Fix that - supply the input
    data, index the CRAM, add the missing step to steps: - and run again with the same flag;
    or remove the members that have no data from the pedigree, if they were never sequenced."""
            : """Then either
      - run again with --clean-stale-families, which deletes each family's norm.bcf and the
        annotation and wombat outputs built from it so the family is re-called in full. It
        only touches families whose missing members can actually be re-called from data on
        disk and from the steps you asked for, and reports the rest. Add --dry-run first to
        see exactly what it would remove; or
      - remove the members that have no data from the pedigree, if they were never sequenced.

    Ancestry outputs need no action: a member with no gVCF has no panel genotypes either, so
    the family's ancestry merge and scores are already re-scheduled on their own."""

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

    ${remedy}

    Set pedigree_strict: true to stop the run on this instead of warning.
    ========================================================================================
    """
    }
}

/*
========================================================================================
    STALE FAMILY CLEAN
========================================================================================
*/

// Delete one file, but only from inside the data tree.
//
// The containment check is not defensive programming for its own sake: every path handed here
// is built by string interpolation from params.data, and this is the only code in the pipeline
// that unlinks anything. A malformed data value must fail loudly rather than quietly become an
// unlink somewhere else on shared project storage.
//
// Records into `deleted` / `failed` rather than throwing, so one bad path cannot abandon a
// family half-cleaned.
def deleteWithinData(String path, String data_root, boolean rehearsal, List deleted, Map failed) {
    try {
        def target = new File(path)
        def canonical = target.canonicalPath
        if (canonical != data_root && !canonical.startsWith(data_root + File.separator)) {
            failed[path] = 'refusing to delete a path outside the data directory'
            return false
        }
        if (rehearsal) {
            if (target.exists()) {
                log.info "  would delete: ${path}"
                deleted.add(path)
            }
            return true
        }
        if (java.nio.file.Files.deleteIfExists(java.nio.file.Paths.get(path))) {
            log.info "  deleted: ${path}"
            deleted.add(path)
        }
        return true
    }
    catch (Exception e) {
        failed[path] = e.toString()
        return false
    }
}

// Delete the derived outputs of families whose pedigree has drifted, so the next plan
// re-schedules them. This is the only code in the pipeline that destroys data, so the whole
// function is built around one rule: never delete what this run cannot rebuild.
//
// A family that loses its norm.bcf and then cannot be re-called is strictly worse off than
// the warning this replaces - the call set is gone, nothing schedules the missing member, and
// the run aborts in validateStepsAvailability having already destroyed the evidence. So every
// family is gated on both halves of "can we rebuild it": the data on disk, and the steps this
// run was actually asked to perform. Families that fail are reported, with the reason, and
// left completely untouched.
//
// Returns a JSON-safe report for the run-state file. Never throws.
def cleanStaleFamilies(analysis_plan, resolved_inputs, family_members) {
    def rehearsal = isRehearsalRun()
    def report = [
        // On a rehearsal these say what *would* happen. The flag is what stops the caller
        // re-planning as though the files were gone, and what lets the summary say so.
        rehearsal             : rehearsal,
        families_cleaned      : [],
        families_skipped      : [],
        cohort_outputs_deleted: [],
        cohort_outputs_stale  : [],
        errors                : []
    ]

    def verb = rehearsal ? 'would delete' : 'deleted'

    def wombat_configs = (params.wombat_config_list ?: []) as List
    def with_wisecondorx = 'wisecondorx' in pipeline_steps

    if (!wombat_configs) {
        log.info "No wombat_config_list is set, so no wombat result TSVs will be removed"
    }
    if (!with_wisecondorx) {
        log.info "The wisecondorx step was not requested, so its family aberration BEDs are left in place"
    }

    // Every path below is built by interpolation from params.data, which is user-supplied
    def data_root = new File(params.data as String).canonicalPath

    analysis_plan.pedigree_drift.each { fid, drift ->
        def verdict = StaleFamily.recovery(drift.missing,
                                           analysis_plan.alignment.existing as Set,
                                           resolved_inputs.resolvable as Set,
                                           resolved_inputs.index_missing ?: [:],
                                           pipeline_steps)
        if (!verdict.recoverable) {
            report.families_skipped.add([family_id: fid, members: drift.members,
                                         missing: drift.missing,
                                         reason: verdict.reasons.join('; ')])
            return
        }

        def paths = StaleFamily.familyOutputs(params.data as String, fid,
                                              params.vep_config_name as String,
                                              wombat_configs, with_wisecondorx)
        def deleted = []
        def failed = [:]
        log.info "Family ${fid}: ${verb} the outputs of a joint call that is missing ${drift.missing.join(', ')}"
        paths.each { path -> deleteWithinData(path, data_root, rehearsal, deleted, failed) }

        // The norm.bcf is what makes the family look done. If it survived, the family was not
        // cleaned however many other files went, and calling it a success would be a lie the
        // next run pays for - it would find the family "existing" with half its outputs gone.
        def norm_bcf = "${Sharding.getFamilyDir(params.data as String, fid)}/vcfs/${fid}.norm.bcf".toString()
        if (failed[norm_bcf] || (!rehearsal && new File(norm_bcf).exists())) {
            def why = failed[norm_bcf] ?: 'it is still present after the delete'
            log.warn "Family ${fid} was NOT cleaned: its norm.bcf could not be removed (${why}) - the family still looks complete and the other files removed from it are now missing"
            report.errors.add("${fid}: norm.bcf could not be removed (${why})".toString())
            report.families_skipped.add([family_id: fid, members: drift.members,
                                         missing: drift.missing,
                                         reason: "norm.bcf could not be removed: ${why}".toString()])
            return
        }

        failed.each { path, why ->
            log.warn "Could not delete ${path}: ${why}"
            report.errors.add("${fid}: ${path}: ${why}".toString())
        }
        report.families_cleaned.add([family_id: fid, members: drift.members,
                                     missing: drift.missing, files_deleted: deleted.size()])
    }

    // The cohort merges keep data from call sets that no longer exist, and the plan only
    // checks that they exist - so nothing would ever rebuild them. But deleting one when this
    // run cannot re-merge it is worse than leaving it: the merge would rebuild from whatever
    // families happen to be done right now and write a *new* wrong file with a fresh
    // timestamp, which looks current and is not. Each merge is rebuilt by its own steps.
    if (report.families_cleaned) {
        def deleted = []
        def stale = []
        StaleFamily.cohortOutputs(params.data as String, params.cohort_name as String,
                                  params.vep_config_name as String,
                                  wombat_configs, with_wisecondorx).each { group ->
            if (group.steps.every { step -> step in pipeline_steps }) {
                def failed = [:]
                log.info "Cohort ${params.cohort_name}: ${verb} the cohort merges rebuilt by ${group.steps.join(', ')}, which contain data from the call sets just removed"
                group.paths.each { path -> deleteWithinData(path, data_root, rehearsal, deleted, failed) }
                failed.each { path, why ->
                    log.warn "Could not delete ${path}: ${why}"
                    report.errors.add("cohort: ${path}: ${why}".toString())
                }
            } else {
                def present = group.paths.findAll { new File(it).exists() }
                stale.addAll(present)
                if (present) {
                    log.warn """The cohort merges below still contain data from the call sets just removed, and were
    left in place because this run cannot rebuild them - ${group.steps.collect { "'${it}'" }.join(', ')} must
    all be in steps. Add them and run again, or delete these by hand once the families are
    re-annotated:
      ${present.join('\n      ')}"""
                }
            }
        }
        report.cohort_outputs_deleted = deleted.collect { new File(it).name }
        report.cohort_outputs_stale = stale.collect { new File(it).name }
    }

    def n_clean = report.families_cleaned.size()
    def n_skip = report.families_skipped.size()
    def n_files = report.families_cleaned.sum { it.files_deleted } ?: 0
    def to_align = report.families_cleaned.collectMany { entry ->
        entry.missing.findAll { barcode -> !(barcode in analysis_plan.alignment.existing) }
    }.unique()

    log.info """
    ========================================================================================
                    ${rehearsal ? 'STALE FAMILY CLEAN (rehearsal - nothing was deleted)' : 'CLEANED STALE FAMILY OUTPUTS'}
    ========================================================================================
    ${n_clean} ${n_clean == 1 ? 'family' : 'families'} ${rehearsal ? 'would be cleaned' : 'cleaned'} (${n_files} ${rehearsal ? 'files would be removed' : 'files removed'}), ${n_skip} left alone.
    ${report.cohort_outputs_deleted ? "Cohort merges ${rehearsal ? 'that would be removed' : 'removed'}: ${report.cohort_outputs_deleted.join(', ')}" : "No cohort merges ${rehearsal ? 'would be removed' : 'removed'}."}
    ${!n_clean ? 'Nothing was cleaned, so nothing new needs aligning.' : to_align ? "${to_align.size()} sample(s) must be aligned from scratch before their families can be re-called: ${to_align.join(', ')}" : 'Every missing member already has a CRAM; only variant calling is needed.'}
    ${report.families_skipped ? "\n    Left alone:\n      " + report.families_skipped.collect { "${it.family_id}: ${it.reason}" }.join('\n      ') : ''}
    ${rehearsal ? '\n    This was a rehearsal run (-preview or -stub-run). Nothing on disk was changed.' : ''}
    ========================================================================================
    """

    return report
}

def validateStepsAvailability(analysis_plan, resolved_inputs, family_members) {
    def errors = []
    
    // Check if alignment is needed but not available
    if (analysis_plan.alignment.needed.size() > 0 && !('alignment' in pipeline_steps)) {
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
    if (analysis_plan.deepvariant_sample.needed.size() > 0 && !('deepvariant_sample' in pipeline_steps)) {
        errors.add("DeepVariant sample step is required for ${analysis_plan.deepvariant_sample.needed.size()} individuals but not included in steps parameter")
    }
    
    // A missing prerequisite is only an error when something in this run would
    // consume it. The ancestry step reads gVCFs straight from deepvariant_sample and
    // never touches the normalized or annotated call sets, so `steps: ["ancestry"]`
    // must not be blocked by families that have no norm.bcf yet. Every steps list
    // that includes one of the consumers below behaves exactly as before.
    def family_consumers = ['deepvariant_family', 'annotation', 'wombat', 'wombat_cohort', 'snvs_cohort', 'extractor']
    def annotation_consumers = ['annotation', 'wombat', 'wombat_cohort', 'snvs_cohort']

    // Check if deepvariant_family is needed but not available
    if (analysis_plan.deepvariant_family.needed.size() > 0 && !('deepvariant_family' in pipeline_steps) &&
        family_consumers.any { step -> step in pipeline_steps }) {
        errors.add("DeepVariant family step is required for ${analysis_plan.deepvariant_family.needed.size()} families but not included in steps parameter")
    }
    
    // Check if annotation is needed but not available
    if (analysis_plan.annotation.needed.size() > 0 && !('annotation' in pipeline_steps) &&
        annotation_consumers.any { step -> step in pipeline_steps }) {
        errors.add("Annotation step is required for ${analysis_plan.annotation.needed.size()} families but not included in steps parameter")
    }

    // roh reads the ancestry step's outputs and nothing else, so it is blocked exactly when they
    // are absent and this run will not make them
    if (analysis_plan.roh.blocked.size() > 0 && 'roh' in pipeline_steps) {
        errors.add("Ancestry step is required for ${analysis_plan.roh.blocked.size()} families before ROH can be called (no family panel or ancestry table under the label ${params.ancestry_panel_name}) but not included in steps parameter")
    }

    // A cohort wombat table is only merged from every family's table, so without the wombat
    // step a family that has none would only make the merge skip itself with a warning
    if (analysis_plan.wombat.needed.size() > 0 && !('wombat' in pipeline_steps) && 'wombat_cohort' in pipeline_steps &&
        analysis_plan.wombat_cohort.needed) {
        errors.add("Wombat step is required for ${analysis_plan.wombat.needed.size()} families but not included in steps parameter")
    }
    
    if (errors) {
        log.error """
        ========================================================================================
                                        ERRORS DETECTED
        ========================================================================================
        ${errors.join('\n        ')}
        
        Please add the required steps to your parameters or ensure all required files exist.
        Available steps: alignment, deepvariant_sample, deepvariant_family, annotation, wombat, wombat_cohort, snvs_cohort, wisecondorx, extractor, ancestry, roh
        ========================================================================================
        """
        recordFailedRun("${errors.size()} unmet step requirement${errors.size() == 1 ? '' : 's'}")
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

// Record which source a barcode's alignment input came from, for the run-state report.
//
// A barcode can sit in more than one source, and the sources are scanned in a fixed order
// that has nothing to do with which one is usable. An unindexed CRAM found first must not be
// the answer for a barcode that also has a good FASTQ pair, or the report would call a
// perfectly alignable sample blocked - so a usable source always wins, whenever it turns up.
def noteAlignmentSource(resolved, String barcode, String source, String path, boolean usable) {
    def existing = resolved.by_barcode[barcode]
    if (existing && (existing.usable || !usable)) return
    resolved.by_barcode[barcode] = [source: source, path: path, usable: usable]
}

// Resolve the actual input files behind every individual that has no CRAM.
//
// alignment.needed means "the output CRAM is missing and something downstream wants it" -
// it says nothing about an input existing. Resolving inputs here, eagerly and once, gives
// validateStepsAvailability something real to check and gives createChannels its rows, so
// the plan and the channels describe the same set of work.
//
// The probe set is deliberately wider than alignment.needed: it is every CRAM-less individual.
// Two things need that. The run-state report answers "can this sample be aligned?" for samples
// nothing is scheduled for, and --clean_stale_families cannot decide whether a drifted family
// is recoverable without it - the missing member of an already-called family is precisely the
// case that never reaches alignment.needed, which is usually empty exactly when it matters.
// The globs are whole-directory scans either way, so the wider set costs only in-memory work.
//
// Rows in cram_37/cram_38 therefore cover barcodes that must NOT be realigned; createChannels
// filters them down to alignment.needed.
def resolveAlignmentInputs(analysis_plan) {
    def resolved = [
        cram_37       : [],
        cram_38       : [],
        resolvable    : [] as Set,
        fastq_barcodes: [] as Set,
        index_missing : [:],
        unpaired_fastq: [:],
        searched      : [],
        // barcode -> [source: String, path: String|null, usable: boolean], the JSON-safe
        // verdict the state file records. One answer per barcode, even when several sources
        // hold it.
        by_barcode    : [:]
    ]

    if (analysis_plan.alignment.no_cram.size() == 0) {
        return resolved
    }

    def needed = analysis_plan.alignment.no_cram as Set

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
            def source = key == 'cram_37' ? 'old_cram_37' : 'old_cram_38'
            if (new File(crai_path).exists()) {
                resolved[key].add([barcode, cram, file(crai_path)])
                resolved.resolvable.add(barcode)
                noteAlignmentSource(resolved, barcode, source, cram.toString(), true)
            } else {
                // Recorded rather than dropped: the data is there, only the index is not,
                // and that needs a different fix than missing input
                resolved.index_missing[barcode] = [source: source, path: cram.toString()]
                noteAlignmentSource(resolved, barcode, "${source}_unindexed".toString(), cram.toString(), false)
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
            def barcode = unit_barcode[unit]
            if (n == 2) {
                resolved.resolvable.add(barcode)
                resolved.fastq_barcodes.add(barcode)
                noteAlignmentSource(resolved, barcode, 'fastq', null, true)
            } else {
                resolved.unpaired_fastq[unit] = n
                noteAlignmentSource(resolved, barcode, 'fastq_unpaired', null, false)
            }
        }
        // Recorded for every probed barcode above, but only worth saying out loud for the ones
        // this run would actually align - a migrated cohort would otherwise bury the log
        resolved.unpaired_fastq.each { unit, n ->
            if (unit_barcode[unit] in analysis_plan.alignment.needed) {
                log.warn "FASTQ unit ${unit} has ${n} read file(s) instead of 2 and will be ignored by fromFilePairs"
            }
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
                   index_missing.collect { barcode, info -> "${barcode}: ${info.path}" }.join("\n        "))
    }

    // Name the downstream work that cannot happen either, so this report covers everything
    // the ANALYSIS SUMMARY above promised
    def blocked = unresolved as Set

    def blocked_samples = analysis_plan.deepvariant_sample.needed.findAll { barcode -> barcode in blocked }
    if (blocked_samples) {
        errors.add("Blocked by the above - variant calling for ${blocked_samples.size()} individual(s): ${blocked_samples.join(', ')}")
    }

    def members_by_family = membersByFamily(family_members)

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
        // Narrowed to what this run will align: the scan probes every CRAM-less individual,
        // and the channel below filters to alignment.needed, so comparing against the raw
        // scan would report every unscheduled sample that happens to have FASTQ as "dropped"
        def fastq_expected = resolved_inputs.fastq_barcodes.findAll { barcode ->
            barcode in analysis_plan.alignment.needed
        }
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
    // second glob - re-globbing here is what let the channels and the plan drift apart.
    //
    // The scan resolves inputs for every CRAM-less individual, not just the ones being
    // aligned, so that the run can report on samples nothing is scheduled for. Those extra
    // rows must not reach ALIGNMENT: realigning a sample that already has what it needs is
    // days of compute for no change.
    channels.cram_37_files = Channel.fromList(resolved_inputs.cram_37)
        .filter { barcode, _cram, _crai -> barcode in analysis_plan.alignment.needed }
    channels.cram_38_files = Channel.fromList(resolved_inputs.cram_38)
        .filter { barcode, _cram, _crai -> barcode in analysis_plan.alignment.needed }
    
    // Every "existing" channel below is built from the plan's own ID buckets and the sharded
    // path each ID implies, never by globbing samples/ or families/. The planner has already
    // proved each of these files exists - that is what put the ID in the bucket - so a glob
    // could only re-discover the same set the slow way: each one was a walk of the whole tree,
    // a dozen of them a run, over every cohort sharing the data directory, before the first
    // task could start. Built this way a path is also the one the planner checked, rather
    // than any file under any directory whose name happens to match.
    def sampleFile = { String barcode, String rel -> file("${Sharding.getSampleDir(params.data, barcode)}/${rel}") }
    def familyFile = { String fid, String rel -> file("${Sharding.getFamilyDir(params.data, fid)}/${rel}") }
    def cramRow = { String barcode ->
        def cram = sampleFile(barcode, "sequences/${barcode}.${params.ref_name}.cram")
        [barcode, cram, file("${cram}.crai")]
    }

    // CRAMs on disk, and the subset that only needs a coverage bedgraph (no realignment).
    // alignment.existing requires the .crai, and need_bedgraph is a subset of it.
    channels.existing_crams = Channel.fromList(analysis_plan.alignment.existing.collect { cramRow(it) })
    channels.bedgraph_only_crams = Channel.fromList(analysis_plan.alignment.need_bedgraph.collect { cramRow(it) })

    // gVCFs: deepvariant_sample.existing requires both the gVCF and its .tbi
    channels.existing_gvcfs = Channel.fromList(analysis_plan.deepvariant_sample.existing.collect { barcode ->
        def gvcf = sampleFile(barcode, "deepvariant/${barcode}.g.vcf.gz")
        [barcode, gvcf, file("${gvcf}.tbi")]
    })

    // Normalized BCFs and family pedigrees: deepvariant_family.existing requires the norm.bcf,
    // its .csi and a non-empty {FID}.pedigree.tsv
    channels.existing_normalized_bcfs = Channel.fromList(analysis_plan.deepvariant_family.existing.collect { fid ->
        def bcf = familyFile(fid, "vcfs/${fid}.norm.bcf")
        [fid, bcf, file("${bcf}.csi")]
    })
    channels.existing_family_pedigrees = Channel.fromList(analysis_plan.deepvariant_family.existing.collect { fid ->
        [fid, familyFile(fid, "${fid}.pedigree.tsv")]
    })

    // Common filtered BCFs (used by snvs_cohort) and annotated BCFs (used by wombat):
    // annotation.existing requires every annotation output and its index
    channels.existing_common_filtered_bcfs = Channel.fromList(analysis_plan.annotation.existing.collect { fid ->
        def bcf = familyFile(fid, "vcfs/${fid}.common_gt.bcf")
        [fid, bcf, file("${bcf}.csi")]
    })
    channels.existing_annotated_bcfs = Channel.fromList(analysis_plan.annotation.existing.collect { fid ->
        def bcf = familyFile(fid, "vcfs/${fid}.rare.${params.vep_config_name}.annotated.bcf")
        [fid, bcf, file("${bcf}.csi")]
    })

    // Create channel for existing Wombat files (output from wombat, used by wombat_cohort)
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
    
    // NPZ files that exist but still need predict: need_predict minus need_npz is exactly the
    // samples whose NPZ the planner found
    channels.existing_npz_files = Channel.fromList(analysis_plan.wisecondorx.need_predict
        .findAll { barcode -> !(barcode in analysis_plan.wisecondorx.need_npz) }
        .collect { barcode -> [barcode, sampleFile(barcode, "svs/wisecondorx/${barcode}.${params.wisecondorx_binsize}.npz")] })
    
    return channels
}

/*
========================================================================================
    COMPLETION HANDLER
========================================================================================
*/

// Registered at script level on purpose. An onComplete closure is delegated to the *script*
// binding, so one registered inside `workflow {}` would see a per-run binding in which even
// `workflow` and `params` read back as null - and, because the delegate is a Map, those
// misses return null silently instead of raising.
//
// This fires on success and on task failure. It does NOT fire for `exit 1`, Ctrl-C or a
// SLURM kill; those are covered by recordFailedRun and by the 'running' marker.
workflow.onComplete {
    def pedigree_file = ghfc_run_state.pedigree_file
    def pedigree_data = ghfc_run_state.pedigree_data
    def plan = ghfc_run_state.plan
    def scan = ghfc_run_state.alignment_scan
    def measured = 'before'

    try {
        if (pedigree_file) {
            // Re-scan now that publishing has finished. This is the only honest answer to
            // "how complete is this cohort?" - the plan built at the start of the run
            // describes the tree as it was before this run published anything.
            pedigree_data = parsePedigreeFile(pedigree_file, true)
            plan = createAnalysisPlan(pedigree_data.families, pedigree_data.individuals,
                                      pedigree_data.family_members, true)
            // And re-resolve inputs against it, so a sample this run aligned drops out of the
            // CRAM-less report instead of being carried over from the start of the run
            scan = resolveAlignmentInputs(plan)
            measured = 'after'
        }
    }
    catch (Throwable t) {
        log.warn "Could not re-scan outputs for ${CohortState.FILE_NAME} (${t}) - recording the counts from the start of the run instead"
        pedigree_data = ghfc_run_state.pedigree_data
        plan = ghfc_run_state.plan
        scan = ghfc_run_state.alignment_scan
        measured = 'before'
    }

    try {
        // Last word of the live file, with the same re-scan the state file gets below
        ghfc_run_state.progress?.stop(workflow.success ? 'success' : 'failed',
            measured == 'after'
                ? buildCompletionBlock(plan, pedigree_data?.families?.size(), pedigree_data?.individuals?.size())
                : null)
    }
    catch (Throwable t) {
        log.warn "Could not finish ${RunProgress.FILE_NAME}: ${t}"
    }

    recordRunState(
        status: workflow.success ? 'success' : 'failed',
        plan: plan,
        pedigree_data: pedigree_data,
        pedigree_file: pedigree_file,
        alignment_scan: scan,
        measured: measured,
        // On a failed run Nextflow cancels in-flight publishDir copies, so outputs this run
        // produced may not have landed yet and the counts above can under-report.
        incomplete: !workflow.success
    )
}
