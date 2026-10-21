/**
 * What main.nf's functions share over one run: the steps it was asked for, and what it has
 * learnt so far that its state file reports (see the COHORT RUN STATE section of main.nf).
 *
 * Static fields, and this is not a style choice. The strict syntax parser, the default from
 * Nextflow 26.04, refuses script-level variables, and a function cannot see a variable that
 * `workflow {}` assigns. A class field is the same one from the workflow, from every function,
 * from the completion handler and from RunProgress's rescan thread. Nextflow runs one pipeline
 * per JVM, so a run never sees another's values.
 */
class RunContext {

    /** params.steps as a list, whichever form it came in; main.nf sets it before any check. */
    static List<String> steps = []

    /**
     * Mutated, never reassigned: the first terminal record and the no-cohort warning are each
     * written once; the rest is what the run stashed for its state and progress files.
     */
    static final Map state = [terminal_written: false, warned_no_cohort: false,
                              pedigree_file: null, pedigree_data: null, plan: null,
                              alignment_scan: null, stale_family_clean: null, progress: null]
}
