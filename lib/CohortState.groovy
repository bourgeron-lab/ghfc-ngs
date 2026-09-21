import groovy.json.JsonOutput
import groovy.json.JsonSlurper

import java.nio.charset.StandardCharsets
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardCopyOption
import java.nio.file.attribute.PosixFilePermissions
import java.security.MessageDigest

/**
 * Reads and writes .ghfc-ngs.state.json, the per-cohort record of what was last run against a
 * cohort, with which pedigree and parameters, and how far each step has got.
 *
 * Deliberately free of any Nextflow dependency - no `log`, no `params`, no `file()` - so that
 * every rule in here can be exercised with a plain `groovy -cp lib` script. The caller in
 * main.nf owns all the logging and all the knowledge of what a "step" is.
 *
 * See COHORT_STATE.md for the schema and the guarantees this class is expected to provide.
 */
class CohortState {

    static final String FILE_NAME = '.ghfc-ngs.state.json'
    static final String TMP_PREFIX = '.ghfc-ngs.state.'
    static final int SCHEMA_VERSION = 1
    static final int HISTORY_LIMIT = 10

    /** Key order in the written file, which is also the order a human wants to read it in. */
    private static final List TOP_LEVEL_KEYS =
        ['schema_version', 'cohort_name', 'last_run', 'last_successful_run', 'history']

    // ------------------------------------------------------------------------------------
    // Checksums
    // ------------------------------------------------------------------------------------

    /**
     * SHA-256 of a file, streamed so a large pedigree costs no memory. Returns null when the
     * path is missing or is not a regular file, because "we could not hash it" and "it hashed
     * to something" must stay distinguishable in the state file.
     */
    static String sha256(String path) {
        if (!path) return null
        def source = new File(path)
        if (!source.isFile()) return null

        def digest = MessageDigest.getInstance('SHA-256')
        source.withInputStream { stream ->
            byte[] buffer = new byte[65536]
            int read
            while ((read = stream.read(buffer)) > 0) {
                digest.update(buffer, 0, read)
            }
        }
        return digest.digest().encodeHex().toString()
    }

    /** SHA-256 of a string, UTF-8 encoded. Used for the serialized effective params map. */
    static String sha256OfString(String text) {
        if (text == null) return null
        def digest = MessageDigest.getInstance('SHA-256')
        digest.update(text.getBytes(StandardCharsets.UTF_8))
        return digest.digest().encodeHex().toString()
    }

    // ------------------------------------------------------------------------------------
    // Completion arithmetic
    // ------------------------------------------------------------------------------------

    /**
     * One step's completion block. A null total means the step was not measured at all, which
     * is not the same as 0% done, so it yields null rather than a zero block. A total of 0
     * (an empty pedigree) yields a null pct rather than dividing by zero.
     */
    static Map progress(Integer done, Integer total) {
        if (done == null || total == null) return null
        if (total <= 0) return [done: done, total: total, pct: null]
        return [done: done, total: total, pct: Math.round(done * 10000.0d / total) / 100.0d]
    }

    // ------------------------------------------------------------------------------------
    // Read / merge / write
    // ------------------------------------------------------------------------------------

    /**
     * Parse the existing state file. A missing file, an unreadable one and a corrupt one all
     * return an empty map: losing the record of the run we are in the middle of writing would
     * be a worse outcome than losing the history of older ones.
     */
    static Map read(Path dir) {
        if (dir == null) return [:]
        def target = dir.resolve(FILE_NAME)
        if (!Files.isRegularFile(target)) return [:]
        try {
            def parsed = new JsonSlurper().parse(target.toFile(), 'UTF-8')
            return (parsed instanceof Map) ? (Map) parsed : [:]
        }
        catch (Exception ignored) {
            return [:]
        }
    }

    /**
     * Fold a new record into the existing state. Pure: no I/O, no clock, no globals, so the
     * whole lifecycle can be tested by calling this with hand-built maps.
     */
    static Map merge(Map existing, Map record, int historyLimit = HISTORY_LIMIT) {
        def state = existing ? new LinkedHashMap(existing) : [:]
        def history = (state.history instanceof List) ? new ArrayList((List) state.history) : []

        // A 'running' record left behind by a *different* run is a run that died without ever
        // reaching the completion handler - a SLURM walltime kill, a Ctrl-C, a lost launch
        // node. Nothing in that run will ever close it out, so the next run to touch the
        // cohort is the only thing that can turn it into a permanent record.
        def previous = state.last_run
        if (previous instanceof Map && previous.status == 'running' && previous.run_id != record?.run_id) {
            def orphan = new LinkedHashMap((Map) previous)
            orphan.status = 'interrupted'
            history.add(0, orphan)
        }

        state.schema_version = SCHEMA_VERSION
        if (record?.cohort_name) state.cohort_name = record.cohort_name
        state.last_run = record

        if (record?.status == 'success') {
            state.last_successful_run = record
        }

        if (record?.status && record.status != 'running') {
            // This run may already be in history if a terminal record was written twice;
            // replace it rather than listing the same run_id under two entries.
            if (record.run_id != null) {
                history.removeAll { it instanceof Map && it.run_id == record.run_id }
            }
            history.add(0, record)
        }

        state.history = (historyLimit > 0 && history.size() > historyLimit)
            ? history[0..<historyLimit]
            : history

        return reorder(state)
    }

    /**
     * Replace the state file atomically. Readers either see the whole previous version or the
     * whole new one, never a half-written file - which matters because this is a shared path
     * that other people's tooling reads while runs are in flight.
     */
    static void write(Path dir, Map state) {
        Files.createDirectories(dir)
        def target = dir.resolve(FILE_NAME)

        // The temp file has to live in the target directory. Files.move can only be atomic
        // within a single filesystem, and $TMPDIR is never the same mount as the project
        // storage, so a temp file anywhere else turns this into a cross-device copy.
        def tmp = Files.createTempFile(dir, TMP_PREFIX, '.tmp')
        try {
            tmp.toFile().withWriter('UTF-8') { writer ->
                writer.write(JsonOutput.prettyPrint(JsonOutput.toJson(state)))
                writer.write('\n')
            }
            // Shared project storage: under a restrictive umask the file would land 0600 and
            // be both unreadable and unreplaceable for everyone else working on the cohort.
            try {
                Files.setPosixFilePermissions(tmp, PosixFilePermissions.fromString('rw-rw-r--'))
            }
            catch (Exception ignored) {
                // Non-POSIX filesystem, or no permission to chmod. The file is still written.
            }
            Files.move(tmp, target, StandardCopyOption.REPLACE_EXISTING, StandardCopyOption.ATOMIC_MOVE)
        }
        finally {
            // A no-op after a successful move; on failure it stops a stray .tmp being left behind
            Files.deleteIfExists(tmp)
        }
    }

    private static Map reorder(Map state) {
        def ordered = [:]
        TOP_LEVEL_KEYS.each { key -> if (state.containsKey(key)) ordered[key] = state[key] }
        state.each { key, value -> if (!ordered.containsKey(key)) ordered[key] = value }
        return ordered
    }
}
