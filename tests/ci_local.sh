#!/usr/bin/env bash

# The CI checks, locally, in the order CI runs them. Stops at the first failing group.
#
#   tests/ci_local.sh                 # Nextflow 25.10.4, the cluster's version
#   NXF_VER=24.10.5 tests/ci_local.sh
#
# Nextflow is run through its official launcher, which honours NXF_VER - a Homebrew nextflow
# does not, and 26.04's strict parser cannot read this repository. The launcher is fetched once
# into ~/.cache/ghfc-ngs. Linting is skipped if shellcheck is not installed, and the Python 3.6
# run of scripts/ if docker is not; CI runs both.

set -euo pipefail

cd "$(dirname "$0")/.."

export NXF_VER=${NXF_VER:-25.10.4}
export NXF_ANSI_LOG=false
if [[ -z "${NEXTFLOW:-}" ]]; then
    NEXTFLOW="${XDG_CACHE_HOME:-${HOME}/.cache}/ghfc-ngs/nextflow"
    if [[ ! -x "${NEXTFLOW}" ]]; then
        mkdir -p "$(dirname "${NEXTFLOW}")"
        curl -fsSL https://get.nextflow.io -o "${NEXTFLOW}"
        chmod +x "${NEXTFLOW}"
    fi
fi
export NEXTFLOW

section() { printf '\n##### %s\n' "$1"; }

section "lint"
if command -v shellcheck > /dev/null; then tests/lint.sh; else echo "skipped: shellcheck not installed"; fi

section "consistency"
tests/checks/params_doc.py
tests/checks/config_refs.py
tests/checks/repo_hygiene.sh

section "script tests"
tests/run_script_tests.sh

section "script tests, scripts/ on Python 3.6"
if command -v docker > /dev/null; then
    docker run --rm -v "${PWD}:/repo" -w /repo python:3.6-slim bash tests/run_script_tests.sh scripts
else
    echo "skipped: docker not installed"
fi

section "lib/ unit tests (Nextflow ${NXF_VER})"
tests/lib/run.sh

section "stub runs (Nextflow ${NXF_VER})"
tests/stub/run.sh

section "all checks passed"
