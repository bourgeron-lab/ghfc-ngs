#!/usr/bin/env bash

# Lints every tracked script, found by its shebang so a new one is covered without being listed:
#   - shell scripts with shellcheck, at warning level;
#   - Python scripts with ruff, syntax errors and pyflakes only (undefined names, unused
#     imports...) - no style rules. Python 3.7 is the oldest target ruff knows; the operator
#     tools under scripts/ are also run on 3.6 by CI.
# Needs shellcheck and ruff on PATH (or uvx, for ruff).
#
#   tests/lint.sh

set -uo pipefail

cd "$(dirname "$0")/.." || exit 1

shell=() python=()
while IFS= read -r path; do
    [[ -f "${path}" && "${path}" != *.nf ]] || continue
    case "$(head -n 1 "${path}")" in
        '#!'*python*) python+=("${path}") ;;
        '#!'*bash*|'#!'*/sh) shell+=("${path}") ;;
    esac
done < <(git ls-files)

ruff=(ruff)
command -v ruff > /dev/null || ruff=(uvx ruff)

status=0
echo "=== shellcheck: ${#shell[@]} scripts"
shellcheck -S warning "${shell[@]}" || status=1
echo "=== ruff: ${#python[@]} scripts"
"${ruff[@]}" check --select E9,F --target-version py37 "${python[@]}" || status=1
exit "${status}"
