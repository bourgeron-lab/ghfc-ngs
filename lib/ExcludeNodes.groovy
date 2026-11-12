/**
 * params.exclude_nodes, the Slurm nodes a run's tasks must not be placed on, as a list.
 *
 * It arrives as a comma-separated string from the command line (`--exclude_nodes a,b`), as a
 * string or a list from a parameters file, as a number when a node name is all digits, and as
 * the Boolean true when the flag is given with no value. nextflow.config turns the same value
 * into the `--exclude=` of every task's clusterOptions (process.ext.slurm_options), and cannot
 * call this class - a config is compiled before lib/ is on the classpath - so the two must agree
 * on how a value is split: on commas, each item trimmed, empty ones dropped.
 */
class ExcludeNodes {

    /** A Slurm node name; no hostlist expressions (`maestro-[1082,1096]`), no spaces. */
    static final String NODE_PATTERN = '^[A-Za-z0-9._-]+$'

    /** The nodes, in the order given; empty when none. A bare flag (true) is no node at all. */
    static List<String> parse(Object raw) {
        if (raw == null || raw instanceof Boolean) return []
        def items = (raw instanceof Collection) ? (raw as List) : raw.toString().split(',') as List
        return items.collect { it?.toString()?.trim() }.findAll { it } as List<String>
    }

    /**
     * What is wrong with the value, or null when it is usable. Empty, missing and false all mean
     * "exclude nothing", the default.
     */
    static String check(Object raw) {
        if (raw == true) return 'it was given with no value; expected nodes such as maestro-1096,maestro-1082'
        if (raw instanceof Collection && raw.any { it instanceof Collection || it instanceof Map }) {
            return 'it is a nested list; expected a list of node names'
        }
        if (raw instanceof Map) return 'it is a map; expected node names separated by commas'
        def bad = parse(raw).findAll { !(it ==~ NODE_PATTERN) }
        if (bad) {
            return "${bad.collect { "'${it}'" }.join(', ')} ${bad.size() == 1 ? 'is not a node name' : 'are not node names'}" +
                ' (letters, digits, ., _ and - only; separate nodes with commas)'
        }
        return null
    }
}
