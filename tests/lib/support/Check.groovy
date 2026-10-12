/**
 * The smallest harness that will do: named checks, a power-assert message on failure, and an
 * exit status. Each *Test.groovy calls Check.that(...) for every case, then Check.done().
 */
class Check {

    static int passed = 0
    static List<String> failed = []

    static void that(String name, Closure body) {
        try {
            body()
            passed++
            println "  ok    ${name}"
        }
        catch (Throwable e) {
            failed << name
            println "  FAIL  ${name}"
            println e.message.readLines().collect { "        ${it}" }.join('\n')
        }
    }

    static void done() {
        println "${passed} passed, ${failed.size()} failed"
        System.exit(failed ? 1 : 0)
    }
}
