#!/usr/bin/env python3
"""Check that documentation/params.md still describes the parameters the code reads.

The "Every key, A to Z" table is the contract: a key marked ✗ is inert, any other key is live.
Fails when
  - the code reads a key the table does not list, or lists as inert;
  - the table lists a live key that nothing reads any more;
  - a key set in nextflow.config's params block, params.yml or params_example/*.yml is not in
    the table, so an operator meets a key nobody explained;
  - the counts in the sentence above the table no longer match it.

    tests/checks/params_doc.py
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOC = ROOT / "documentation" / "params.md"

# Spellings the code accepts that are deliberately not documented as keys of their own
ALIASES = {"cleanStaleFamilies"}  # camelCase form of clean_stale_families, main.nf

# Methods of the params map, not keys
MAP_METHODS = {"containsKey", "get", "keySet", "each", "collect", "findAll", "getOrDefault",
               "size", "isEmpty", "toString", "entrySet", "values"}


def strip_comments(text):
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"(?m)(^|[^:\"'])//.*$", r"\1", text)


def params_block(config):
    """Text of the top-level `params { ... }` block of nextflow.config, by brace matching."""
    start = re.search(r"(?m)^params\s*\{", config)
    if not start:
        return "", config
    depth, i = 0, start.end() - 1
    while i < len(config):
        if config[i] == "{":
            depth += 1
        elif config[i] == "}":
            depth -= 1
            if depth == 0:
                return config[start.end():i], config[:start.start()] + config[i + 1:]
        i += 1
    sys.exit("params_doc: unbalanced params block in nextflow.config")


def keys_read(text):
    keys = set(re.findall(r"\bparams\.([A-Za-z_]\w*)", text))
    keys |= set(re.findall(r"\bparams\[\s*['\"](\w+)['\"]\s*\]", text))
    keys |= set(re.findall(r"\bparams\.(?:get|containsKey)\(\s*['\"](\w+)['\"]", text))
    # Lists of names checked with params[key], e.g. ['ancestry_reference', ...].findAll { key -> !params[key] }
    for names, _var in re.findall(r"\[([^\[\]]*)\]\s*\.findAll\s*\{\s*(\w+)\s*->\s*!?params\[\2\]", text):
        keys |= set(re.findall(r"['\"](\w+)['\"]", names))
    return keys - MAP_METHODS


def yaml_top_keys(path):
    return {m.group(1) for m in re.finditer(r"(?m)^([A-Za-z_]\w*)\s*:", path.read_text())}


def main():
    doc = DOC.read_text()
    section = re.search(r"(?ms)^## Every key, A to Z\n(.*?)^## ", doc)
    if not section:
        sys.exit("params_doc: no '## Every key, A to Z' section in documentation/params.md")
    rows = re.findall(r"(?m)^\| `(\w+)`( [^|]*)?\|", section.group(1))
    live = {key for key, marks in rows if "✗" not in (marks or "")}
    inert = {key for key, marks in rows if "✗" in (marks or "")}

    config = (ROOT / "nextflow.config").read_text()
    block, rest_of_config = params_block(strip_comments(config))
    sources = [rest_of_config]
    for pattern in ("*.nf", "workflows/*.nf", "modules/**/*.nf", "lib/*.groovy"):
        sources += [strip_comments(p.read_text()) for p in ROOT.glob(pattern)]
    read = set().union(*(keys_read(s) for s in sources)) - ALIASES

    declared = {m.group(1) for m in re.finditer(r"(?m)^\s*([A-Za-z_]\w*)\s*=", block)}
    in_files = {"nextflow.config params block": declared}
    for path in [ROOT / "params.yml", *sorted((ROOT / "params_example").glob("*.yml"))]:
        in_files[str(path.relative_to(ROOT))] = yaml_top_keys(path)

    errors = []
    for key in sorted(read - live - inert):
        errors.append(f"`{key}` is read by the code but missing from the A-Z table")
    for key in sorted(read & inert):
        errors.append(f"`{key}` is marked inert (✗) but the code reads it")
    for key in sorted(live - read):
        errors.append(f"`{key}` is documented as live but nothing reads it: mark it ✗ or remove it")
    for where, keys in in_files.items():
        for key in sorted(keys - live - inert):
            errors.append(f"`{key}` is set in {where} but missing from the A-Z table")

    counts = re.search(r"(\d+) keys are read by the code; (\d+) more", section.group(1))
    if not counts:
        errors.append("the 'N keys are read by the code; M more' sentence is gone from the A-Z section")
    elif (int(counts.group(1)), int(counts.group(2))) != (len(live), len(inert)):
        errors.append(f"the A-Z section says {counts.group(1)} live and {counts.group(2)} inert keys; "
                      f"the table has {len(live)} and {len(inert)}")

    if errors:
        print("documentation/params.md is out of step with the code:")
        for e in errors:
            print(f"  - {e}")
        return 1
    print(f"params_doc: ok ({len(live)} live keys, {len(inert)} inert, all read keys documented)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
