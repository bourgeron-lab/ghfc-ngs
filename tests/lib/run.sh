#!/usr/bin/env bash

# Compiles every lib/*.groovy and runs the tests in tests/lib/*Test.groovy against it, with the
# Groovy runtime and the Nextflow classes (ResourceOverrides needs nextflow.util) taken from a
# Nextflow framework jar. Nextflow itself is not run.
#
#   tests/lib/run.sh                  # the jar of $NXF_VER, default 25.10.4, fetched if missing
#   NXF_JAR=/path/to/nextflow-one.jar tests/lib/run.sh

set -euo pipefail

root=$(cd "$(dirname "$0")/../.." && pwd)
ver=${NXF_VER:-25.10.4}
jar=${NXF_JAR:-${NXF_HOME:-${HOME}/.nextflow}/framework/${ver}/nextflow-${ver}-one.jar}
if [[ ! -f "${jar}" ]]; then
    # The launcher downloads the framework jar of NXF_VER on first use
    launcher=$(mktemp -d)/nextflow
    curl -fsSL https://get.nextflow.io -o "${launcher}" && chmod +x "${launcher}"
    NXF_VER=${ver} "${launcher}" -version > /dev/null
fi
[[ -f "${jar}" ]] || { echo "no Nextflow jar at ${jar}" >&2; exit 1; }

classes=$(mktemp -d)
trap 'rm -rf "${classes}"' EXIT

echo "=== compile lib/ against $(basename "${jar}")"
java -cp "${jar}" org.codehaus.groovy.tools.FileSystemCompiler -d "${classes}" \
    "${root}"/lib/*.groovy "${root}"/tests/lib/support/*.groovy

status=0
for t in "${root}"/tests/lib/*Test.groovy; do
    echo "=== ${t#"${root}"/}"
    java -cp "${jar}:${classes}" groovy.ui.GroovyMain "${t}" || status=1
done
exit "${status}"
