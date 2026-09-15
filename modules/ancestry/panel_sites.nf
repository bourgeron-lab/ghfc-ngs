#!/usr/bin/env nextflow

/*
 * Ancestry: reference panel site list module
 * Builds the union of the ancestry panel and the PGS catalog site lists
 */

nextflow.enable.dsl=2

process PANEL_SITES {
  /*
  Verify the reference bundle and its fitted model, and build the union site list

  The two site lists are NOT nested: the LD-pruned ancestry panel used by pcs and
  admixture and the PGS catalog list used by pgs-raw overlap only partially, so a
  single per-sample extraction has to cover their union or PGS cannot be computed
  from it.

  Parameters
  ----------
  cohort_name : val
    Name of the cohort
  reference : val
    Path to the mounted ancestry-pgs reference bundle
  catalog : val
    Path to the PGS weight catalog, checked against the fitted model's checksum
  panel_name : val
    Label identifying this panel/threshold combination in output file names

  Returns
  -------
  Tuple of panel name, union site list, bcftools targets file, and bundle report
  */

  tag "${panel_name}"

  publishDir "${params.data}/cohorts/${cohort_name}/ancestry",
    mode: 'copy',
    pattern: "${panel_name}.{sites.tsv.gz,regions.tsv.gz,bundle.json}"

  label 'panel_sites'

  input:
  tuple val(cohort_name), val(reference), val(catalog), val(panel_name)

  output:
  tuple val(panel_name), path("${sites_file}"), path("${regions_file}"), emit: panel_sites
  path "${bundle_report}", emit: bundle_report

  script:
  model_override = params.ancestry_model ?: ""
  sites_file = "${panel_name}.sites.tsv.gz"
  regions_file = "${panel_name}.regions.tsv.gz"
  bundle_report = "${panel_name}.bundle.json"

  """
  set -euo pipefail

  # Fail loudly and early on a missing or truncated bundle rather than part way
  # through a cohort's worth of extractions. --fast skips the checksums, which on a
  # 771 MB bundle would dominate this task.
  ancestry-pgs verify-bundle --reference ${reference} --fast > ${bundle_report}

  python3 <<'PYEOF'
import gzip, hashlib, json, os, sys

reference = "${reference}"
manifest = json.load(open(f"{reference}/manifest.json"))
k = manifest["admixture"]["k"]


def die(*lines):
    for line in lines:
        print(line, file=sys.stderr)
    raise SystemExit(1)


# --- preflight: the fitted z-score model -------------------------------------
#
# models/ is deliberately absent from the bundle manifest, so the verify-bundle
# call above exits 0 even when the directory is missing entirely. Checking it here
# costs one checksum of the catalog; skipping it means the run gets all the way
# through per-sample extraction and family scoring before pgs-adjusted fails on
# every family.
catalog = "${catalog}"
model_override = "${model_override}"

h = hashlib.sha256()
with open(catalog, "rb") as fh:
    for chunk in iter(lambda: fh.read(1 << 20), b""):
        h.update(chunk)
catalog_sha = h.hexdigest()

model_dir = model_override or f"{reference}/models/{catalog_sha}"
absent = [f for f in ("model.npz", "meta.json")
          if not os.path.exists(f"{model_dir}/{f}")]
if absent:
    die(f"ERROR: no fitted z-score model at {model_dir} "
        f"(missing {', '.join(absent)}).",
        f"  catalog: {catalog}",
        f"  sha256:  {catalog_sha}",
        "Fitting is a one-off, not a pipeline step, and the bundle manifest does",
        "not cover models/ - so verify-bundle cannot report this. Either stage the",
        "bundle's models/ directory along with the rest of it, or run:",
        f"  ancestry-pgs fit-zscore-model --reference {reference} --catalog {catalog}",
        "or set ancestry_model to a model directory fitted from this catalog.")

meta = json.load(open(f"{model_dir}/meta.json"))
if meta.get("catalog_sha256") != catalog_sha:
    die(f"ERROR: the model at {model_dir} was fitted from a different catalog.",
        f"  model was fitted from sha256 {meta.get('catalog_sha256')}",
        f"  catalog supplied has  sha256 {catalog_sha}",
        "Refit the model for this catalog, or correct ancestry_catalog.")

print(f"z-score model ok: {model_dir} ({len(meta.get('traits', []))} traits, "
      f"{meta.get('n_pcs')} PCs, {meta.get('n_reference_samples')} reference samples)")


# --- the union site list -----------------------------------------------------
def read(path):
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\\n").split("\\t")
        if header[:5] != ["rsid", "chrom", "pos", "ref", "alt"]:
            die(f"ERROR: unexpected site-list header in {path}: {header[:5]}")
        return [tuple(l.rstrip("\\n").split("\\t")[:5]) for l in fh]

panel = read(f"{reference}/admix/k{k}/snps_in_order.tsv.gz")
catalog_sites = read(f"{reference}/xwalk/catalog_sites.tsv.gz")

# Dedupe on the identity the projection actually matches on - position and allele
# set - not on the rsID, which the two lists do not always agree about.
union = {}
for rsid, chrom, pos, ref, alt in panel + catalog_sites:
    union.setdefault((chrom, int(pos), ref, alt), rsid)

rows = sorted(union.items(), key=lambda kv: (int(kv[0][0]), kv[0][1]))
if not rows:
    die("ERROR: the reference bundle yielded no panel sites")

with gzip.open("${sites_file}", "wt") as sites, gzip.open("${regions_file}", "wt") as regions:
    sites.write("rsid\\tchrom\\tpos\\tref\\talt\\n")
    for (chrom, pos, ref, alt), rsid in rows:
        sites.write(f"{rsid}\\t{chrom}\\t{pos}\\t{ref}\\t{alt}\\n")
        # 1-based CHROM/POS targets, not BED: bcftools infers BED from the file
        # suffix, which a .bed.gz would not reliably trigger.
        #
        # Both contig namings are emitted. A streaming -T never matches a contig
        # the file does not contain, so the unused half is inert - and this avoids
        # assuming whether the gVCFs say "chr1" or "1".
        regions.write(f"{chrom}\\t{pos}\\nchr{chrom}\\t{pos}\\n")

print(f"panel={len(panel)} catalog={len(catalog_sites)} union={len(rows)} "
      f"shared={len(panel) + len(catalog_sites) - len(rows)}")
PYEOF
  """

  stub:
  sites_file = "${panel_name}.sites.tsv.gz"
  regions_file = "${panel_name}.regions.tsv.gz"
  bundle_report = "${panel_name}.bundle.json"
  """
  printf 'rsid\\tchrom\\tpos\\tref\\talt\\nrs1\\t1\\t1000\\tA\\tG\\n' | gzip -c > ${sites_file}
  printf '1\\t1000\\nchr1\\t1000\\n' | gzip -c > ${regions_file}
  echo '{}' > ${bundle_report}
  """
}
