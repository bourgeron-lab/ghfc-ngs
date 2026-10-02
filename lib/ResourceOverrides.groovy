import nextflow.util.Duration
import nextflow.util.MemoryUnit

/**
 * Checks the `resources:` block of a params file when the run starts.
 *
 * The block is applied by sized() in nextflow.config, one task at a time, and nothing there can
 * report a mistake before the first task it concerns is due - possibly hours into a run, and
 * a key that matches nothing is never reported at all. So every rule sized() relies on is
 * checked here first. The two must agree: the schema is in documentation/params.md.
 *
 * nextflow.config cannot call this class (the config is compiled before lib/ is on the
 * classpath), which is why the lookup itself lives there and only the checking lives here.
 */
class ResourceOverrides {

    static final List<String> FIELDS = ['cpus', 'memory', 'memory_step', 'time']
    static final List<String> SECTIONS = ['families', 'samples']

    /**
     * The names a `resources:` key can address, mapped to the process each one runs: every
     * process of modules/, plus every `include { X as Y }` alias, kept in `aliases` too. A task
     * of an aliased process reports the alias as its name, and that is all sized() can see.
     */
    static Map processNames(File projectDir) {
        def declared = [] as Set
        new File(projectDir, 'modules').eachFileRecurse { f ->
            if (f.name.endsWith('.nf')) {
                (f.text =~ /(?m)^\s*process\s+(\w+)\s*\{/).each { m -> declared << m[1] }
            }
        }
        def aliases = [:]
        def scripts = [new File(projectDir, 'main.nf'), new File(projectDir, 'migrate.nf')]
        new File(projectDir, 'workflows').eachFile { f -> if (f.name.endsWith('.nf')) scripts << f }
        scripts.findAll { it.exists() }.each { f ->
            (f.text =~ /(?m)^\s*include\s*\{([^}]*)\}/).each { m ->
                m[1].split(';').each { item ->
                    def parts = item.trim().split(/\s+as\s+/)
                    if (parts.size() == 2) aliases[parts[1].trim()] = parts[0].trim()
                }
            }
        }
        def names = [:]
        declared.each { names[it] = it }
        aliases.each { alias, process -> if (process in declared) names[alias] = process }
        return [names: names, aliases: aliases]
    }

    /**
     * Problems with `resources`, as [errors: [...], warnings: [...]]. An error stops the run; a
     * warning is a value that will be clamped to a ceiling (`max_memory`, `max_cpus`,
     * `max_time`), which sized() does anyway but which is almost certainly not what was meant.
     */
    static Map check(Object resources, Map known, Map ceilings) {
        def errors = []
        def warnings = []
        if (resources == null || resources == '' || (resources instanceof Map && resources.isEmpty())) {
            return [errors: errors, warnings: warnings]
        }
        if (!(resources instanceof Map)) {
            errors << "resources must be a map of process names, not ${describe(resources)}"
            return [errors: errors, warnings: warnings]
        }
        Map names = known.names
        Map aliases = known.aliases
        resources.each { key, entry ->
            def name = key.toString()
            def where = "resources.${name}"
            // BAZAM_BWA_MEM2_REALIGN runs only as its _37 and _38 aliases, so an entry under its
            // own name would silently never apply. withName in a config file matches both, which
            // makes this the one place the two kinds of key differ.
            def runs_as = aliases.findAll { alias, process -> process == name }.keySet().sort()
            if (runs_as) {
                errors << "${where}: ${name} runs as ${runs_as.join(' and ')}, and an entry here has to name the process as it runs. Give those names, or the regex '${name}_.*'"
                return
            }
            if (!names.containsKey(name)) {
                def matched = null
                try {
                    matched = names.keySet().findAll { it ==~ name }
                }
                catch (java.util.regex.PatternSyntaxException e) {
                    errors << "${where}: not a process name, and not a valid regex either (${e.description})"
                    return
                }
                if (!matched) {
                    errors << "${where}: no process is called ${name}, and as a regex it matches none"
                    return
                }
            }
            if (!(entry instanceof Map)) {
                errors << "${where} must be a map of ${FIELDS.join(', ')}, families or samples, not ${describe(entry)}"
                return
            }
            entry.each { field, value ->
                def f = field.toString()
                if (f in SECTIONS) {
                    if (!(value instanceof Map)) {
                        errors << "${where}.${f} must map ids or globs to settings, not ${describe(value)}"
                        return
                    }
                    value.each { id, settings ->
                        checkFields("${where}.${f}.${id}", settings, ceilings, errors, warnings)
                    }
                }
                else if (!(f in FIELDS)) {
                    errors << "${where}.${f}: unknown setting. Expected ${FIELDS.join(', ')}, families or samples"
                }
            }
            checkFields(where, entry.findAll { k, v -> k.toString() in FIELDS }, ceilings, errors, warnings)
        }
        return [errors: errors, warnings: warnings]
    }

    private static void checkFields(String where, Object settings, Map ceilings, List errors, List warnings) {
        if (!(settings instanceof Map)) {
            errors << "${where} must be a map of ${FIELDS.join(', ')}, not ${describe(settings)}"
            return
        }
        settings.each { field, value ->
            def f = field.toString()
            def at = "${where}.${f}"
            switch (f) {
                case 'cpus':
                    if (!(value instanceof Integer || value instanceof Long) || value < 1) {
                        errors << "${at}: expected a whole number of CPUs, got ${describe(value)}"
                    }
                    else if (ceilings.cpus != null && value > (ceilings.cpus as int)) {
                        warnings << "${at}: ${value} is above max_cpus (${ceilings.cpus}) and will be clamped to it"
                    }
                    break
                case 'memory':
                case 'memory_step':
                    def memory = parse(value, MemoryUnit)
                    if (memory == null) {
                        errors << "${at}: expected a memory with its unit, such as 8.GB, got ${describe(value)}"
                    }
                    else if (f == 'memory' && ceilings.memory && memory > (ceilings.memory as MemoryUnit)) {
                        warnings << "${at}: ${memory} is above max_memory (${ceilings.memory}) and will be clamped to it"
                    }
                    break
                case 'time':
                    def time = parse(value, Duration)
                    if (time == null) {
                        errors << "${at}: expected a duration with its unit, such as 4.h, got ${describe(value)}"
                    }
                    else if (ceilings.time && time > (ceilings.time as Duration)) {
                        warnings << "${at}: ${time} is above max_time (${ceilings.time}) and will be clamped to it"
                    }
                    break
                default:
                    if (!(f in SECTIONS)) errors << "${at}: unknown setting. Expected ${FIELDS.join(', ')}"
                    else errors << "${at}: ${f} belongs directly under a process, not under an id"
            }
        }
    }

    // A bare number is refused rather than read as bytes or milliseconds: 64 meaning 64 bytes
    // is never what was written
    private static Object parse(Object value, Class type) {
        if (!(value instanceof CharSequence)) return null
        try {
            return type == MemoryUnit ? new MemoryUnit(value.toString()) : new Duration(value.toString())
        }
        catch (IllegalArgumentException e) {
            return null
        }
    }

    private static String describe(Object value) {
        value == null ? 'nothing' : value instanceof CharSequence ? "'${value}'" : "${value} (${value.getClass().simpleName})"
    }
}
