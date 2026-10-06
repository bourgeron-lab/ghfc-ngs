import java.nio.file.Path

/**
 * Writes .ghfc-ngs.progress.json, a live view of a run in progress, next to the cohort's state
 * file. Two things go in it, refreshed on two clocks:
 *
 *  - per process, the counters behind Nextflow's own console line ("[a1/b2c3] ALIGN | 12 of 40"):
 *    pending, submitted, running, succeeded, cached, failed. Read from the session's stats
 *    observer, which is what the console renders, every `writeEveryMs`.
 *  - per step, done/total against the pedigree, from the same existence checks as the state
 *    file's `completion` block. That re-scan touches the filesystem for every individual and
 *    family, so it runs every `rescanEveryMs` only, on a thread of its own: on SPARK-GRCh38
 *    (142,357 individuals) one pass outlasted 28 minutes, and the counters must keep coming
 *    meanwhile. A pass also waits RESCAN_SPACING times its own duration before the next one,
 *    so a slow tree is not re-scanned back to back.
 *
 * The first answers "what is Nextflow doing", the second "how far along is the cohort". The
 * process totals cannot answer the second: Nextflow only counts the tasks it has created so far,
 * so a process's total grows as upstream ones emit.
 *
 * The Nextflow dependency is confined to readProcesses(), the default `processes` source, so the
 * rest can be exercised with a plain `groovy -cp lib` script given a fake one. See
 * COHORT_STATE.md, "Live progress", for the schema.
 */
class RunProgress {

    static final String FILE_NAME = '.ghfc-ngs.progress.json'
    static final String TMP_PREFIX = '.ghfc-ngs.progress.'
    static final int SCHEMA_VERSION = 1
    static final int RESCAN_SPACING = 3

    private final Path dir
    private final Map header
    private final Closure processes
    private final Closure rescan
    private final Closure warn
    private final long writeEveryMs
    private final long rescanEveryMs

    private volatile Map completion = null
    private volatile String completionMeasuredAt = null
    private volatile Double completionSeconds = null
    private volatile String scanStartedAt = null
    private volatile boolean stopped = false
    private Thread thread
    private Thread scanner

    /**
     * @param opts.dir           the cohort directory
     * @param opts.header        run_id, cohort_name, started_at, slurm_job_id, launch_dir, host
     * @param opts.rescan        returns a completion block, as CohortState records hold
     * @param opts.processes     returns the per-process records; defaults to Nextflow's
     * @param opts.warn          called with a message when a write or a re-scan fails
     * @param opts.writeEveryMs  default 30 s
     * @param opts.rescanEveryMs default 10 min
     */
    RunProgress(Map opts) {
        this.dir = opts.dir as Path
        this.header = (opts.header ?: [:]) as Map
        this.rescan = opts.rescan as Closure
        this.processes = (opts.processes ?: { boolean exact -> readProcesses(exact) }) as Closure
        this.warn = (opts.warn ?: { String msg -> }) as Closure
        this.writeEveryMs = (opts.writeEveryMs ?: 30_000L) as long
        this.rescanEveryMs = (opts.rescanEveryMs ?: 600_000L) as long
    }

    /** Seed the step counts with the plan the run starts from, then refresh in the background. */
    RunProgress start(Map initialCompletion) {
        setCompletion(initialCompletion, null)
        writeNow('running', false)
        // Daemons, so a run that ends without reaching stop() - an `exit`, a kill - is never
        // held open by them
        thread = Thread.startDaemon('ghfc-ngs-progress') {
            while (!stopped) {
                try {
                    Thread.sleep(writeEveryMs)
                }
                catch (InterruptedException ignored) {
                    break
                }
                if (!stopped) writeNow('running', false)
            }
        }
        if (rescan) {
            scanner = Thread.startDaemon('ghfc-ngs-progress-rescan') {
                long wait = rescanEveryMs
                while (!stopped) {
                    try {
                        Thread.sleep(wait)
                    }
                    catch (InterruptedException ignored) {
                        break
                    }
                    if (stopped) break
                    long took = rescanNow()
                    wait = Math.max(rescanEveryMs, RESCAN_SPACING * took)
                }
            }
        }
        return this
    }

    /**
     * The last write, from the completion handler. `finalCompletion` is the handler's own
     * re-scan, so the file's last word agrees with the state file's.
     */
    void stop(String status, Map finalCompletion) {
        stopped = true
        thread?.interrupt()
        // Not waited for: a pass over a large tree may have minutes to go, and its result
        // would only be overwritten by the completion handler's own re-scan below
        scanner?.interrupt()
        scanStartedAt = null
        if (finalCompletion != null) setCompletion(finalCompletion, null)
        writeNow(status, true)
    }

    /** One re-scan; returns how long it took, in ms. */
    long rescanNow() {
        long t0 = System.currentTimeMillis()
        scanStartedAt = now()
        try {
            def block = rescan.call()
            // A pass that ends after stop() describes the tree before the final one does
            if (!stopped) setCompletion(block, (System.currentTimeMillis() - t0) / 1000.0d)
        }
        catch (Throwable t) {
            if (!stopped) warn.call("Could not re-scan outputs for ${FILE_NAME}: ${t}")
        }
        finally {
            scanStartedAt = null
        }
        return System.currentTimeMillis() - t0
    }

    private void setCompletion(Map block, Double seconds) {
        completion = block
        completionMeasuredAt = now()
        completionSeconds = seconds
    }

    synchronized void writeNow(String status, boolean exact) {
        try {
            CohortState.writeJsonAtomically(dir, FILE_NAME, TMP_PREFIX, snapshot(status, exact))
        }
        catch (Throwable t) {
            // A progress file is a convenience. It must never take down the run it describes.
            warn.call("Could not write ${FILE_NAME}: ${t}")
        }
    }

    Map snapshot(String status, boolean exact) {
        def procs = (processes.call(exact) ?: []) as List<Map>
        def totals = [:]
        ['pending', 'submitted', 'running', 'succeeded', 'cached', 'failed', 'ignored',
         'retries', 'total', 'completed'].each { key ->
            totals[key] = procs.sum { (it[key] ?: 0) as int } ?: 0
        }
        def out = [schema_version: SCHEMA_VERSION]
        out.putAll(header)
        out.status = status
        out.updated_at = now()
        out.totals = totals
        out.processes = procs
        out.completion = completion
        out.completion_measured_at = completionMeasuredAt
        out.completion_scan_seconds = completionSeconds
        // Set while a re-scan is under way: `completion` is then older than this
        out.completion_scan_started_at = scanStartedAt
        return out
    }

    /**
     * Nextflow's per-process records, as its console shows them. `exact` waits for the stats
     * agent to drain; otherwise the last published value is taken without blocking.
     *
     * total = pending + submitted + running + succeeded + failed - retries + cached + stored +
     * aborted, and completed = succeeded + ignored + cached + stored: the console's "N of M".
     */
    static List<Map> readProcesses(boolean exact) {
        def session = nextflow.Global.session
        def observer = session?.statsObserver
        if (observer == null) return []
        def stats = exact ? observer.stats : observer.quickStats
        def records = (stats?.processes ?: []).sort { it.index }
        return records.collect { r ->
            [
                name      : r.name,
                pending   : r.pending,
                submitted : r.submitted,
                running   : r.running,
                succeeded : r.succeeded,
                cached    : r.cached,
                failed    : r.failed,
                ignored   : r.ignored,
                retries   : r.retries,
                total     : r.totalCount,
                completed : r.completedCount,
                terminated: r.terminated,
                errored   : r.errored
            ]
        }
    }

    private static String now() {
        // An explicit pattern: OffsetDateTime.toString() drops ':00' seconds
        return java.time.OffsetDateTime.now().format(
            java.time.format.DateTimeFormatter.ofPattern("yyyy-MM-dd'T'HH:mm:ssXXX"))
    }
}
