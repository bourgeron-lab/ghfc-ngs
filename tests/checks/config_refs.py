#!/usr/bin/env python3
"""Check the references nextflow.config and main.nf make by name, which nothing else verifies.

Fails when
  - a container image is unpinned (no tag, or `latest`) and not in unpinned_containers.txt,
    or that list names an image no longer used (so it only ever shrinks);
  - a `withName:` selector matches no process or include alias - Nextflow ignores it silently,
    so a typo there quietly drops a container, cpus or memory setting;
  - the step list in main.nf's "--steps parameter is required" message differs from
    `valid_steps`, or a valid step has no workflows/<step>.nf.
Warns when a `withLabel:` selector matches no label.

    tests/checks/config_refs.py
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ALLOWLIST = Path(__file__).with_name("unpinned_containers.txt")


def strip_comments(text):
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"(?m)(^|[^:\"'])//.*$", r"\1", text)


def nf_files():
    return [ROOT / "main.nf", ROOT / "migrate.nf", *sorted(ROOT.glob("workflows/*.nf")),
            *sorted(ROOT.glob("modules/**/*.nf"))]


def process_names():
    """Process names and include aliases, as lib/ResourceOverrides.groovy processNames() finds them."""
    names, labels = set(), set()
    for f in nf_files():
        text = strip_comments(f.read_text())
        names |= set(re.findall(r"(?m)^\s*process\s+(\w+)\s*\{", text))
        labels |= set(re.findall(r"(?m)^\s*label\s+['\"]([^'\"]+)['\"]", text))
        for block in re.findall(r"(?m)^\s*include\s*\{([^}]*)\}", text):
            for item in block.split(";"):
                parts = re.split(r"\s+as\s+", item.strip())
                if len(parts) == 2:
                    names.add(parts[1].strip())
    return names, labels


def image_tag(image):
    image = re.sub(r"^\w+://", "", image)
    if "@sha256:" in image:
        return image, "digest"
    last = image.rsplit("/", 1)[-1]
    return image, last.split(":", 1)[1] if ":" in last else None


def main():
    errors, warnings = [], []
    config = strip_comments((ROOT / "nextflow.config").read_text())
    processes, labels = process_names()

    # Containers, in the config and in the modules
    images = {}
    for path, text in [(ROOT / "nextflow.config", config)] + [(f, strip_comments(f.read_text())) for f in nf_files()]:
        for m in re.finditer(r"(?m)^\s*container\s*=?\s*['\"]([^'\"$]+)['\"]", text):
            image, tag = image_tag(m.group(1))
            line = text[:m.start()].count("\n") + 1
            images.setdefault(image, (tag, f"{path.relative_to(ROOT)}:{line}"))
    allowed = {l.strip() for l in ALLOWLIST.read_text().splitlines() if l.strip() and not l.startswith("#")}
    unpinned = {image: where for image, (tag, where) in images.items() if tag in (None, "latest")}
    for image in sorted(set(unpinned) - allowed):
        errors.append(f"container '{image}' ({unpinned[image]}) is not pinned to a version tag or digest")
    for image in sorted(allowed - set(unpinned)):
        errors.append(f"'{image}' is in {ALLOWLIST.name} but no longer used unpinned: remove it from the list")

    # withName / withLabel selectors
    for m in re.finditer(r"\bwithName\s*:\s*['\"]([^'\"]+)['\"]", config):
        pattern = m.group(1)
        line = config[:m.start()].count("\n") + 1
        try:
            regex = re.compile(pattern)
        except re.error as e:
            errors.append(f"nextflow.config:{line} withName '{pattern}' is not a valid pattern: {e}")
            continue
        # Nextflow matches the simple name and the fully qualified WORKFLOW:NAME
        if not any(regex.fullmatch(p) for p in processes):
            errors.append(f"nextflow.config:{line} withName '{pattern}' matches no process")
    for m in re.finditer(r"\bwithLabel\s*:\s*['\"]?([\w|.*!]+)['\"]?", config):
        line = config[:m.start()].count("\n") + 1
        if not any(re.fullmatch(m.group(1), l) for l in labels):
            warnings.append(f"nextflow.config:{line} withLabel '{m.group(1)}' matches no process label")

    # Steps
    main_nf = (ROOT / "main.nf").read_text()
    valid = re.search(r"def valid_steps = \[([^\]]*)\]", main_nf)
    message = re.search(r"--steps parameter is required\. Available steps: ([\w, ]+)", main_nf)
    if not valid or not message:
        errors.append("main.nf: could not find valid_steps or the 'Available steps' message")
    else:
        valid_steps = re.findall(r"'(\w+)'", valid.group(1))
        listed = [s.strip() for s in message.group(1).split(",")]
        if valid_steps != listed:
            errors.append(f"main.nf: 'Available steps' lists {listed} but valid_steps is {valid_steps}")
        for step in valid_steps:
            if not (ROOT / "workflows" / f"{step}.nf").exists():
                errors.append(f"step '{step}' has no workflows/{step}.nf")

    for w in warnings:
        print(f"warning: {w}")
    if errors:
        print("config references are broken:")
        for e in errors:
            print(f"  - {e}")
        return 1
    print(f"config_refs: ok ({len(images)} images, {len(unpinned)} allow-listed unpinned, "
          f"{len(processes)} processes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
