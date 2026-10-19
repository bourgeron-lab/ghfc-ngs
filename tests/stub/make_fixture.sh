#!/usr/bin/env bash

# Builds a fake cohort for the stub runs, all empty files except where the pipeline reads
# content at planning time (the pedigree, the wombat configs). Writes, under <dir>:
#   data/      the cohort tree: a pedigree, FASTQs and old CRAMs
#   refs/      every reference and resource file a step checks for
#   base.yml   the parameters shared by every scenario, with absolute paths
#
# The cohort, STUB:
#   FAM01  BC001 (father, FASTQ, two lanes)  BC002 (mother, FASTQ)  BC003 (child, GRCh38 CRAM)
#   FAM02  BC004 (singleton, GRCh37 CRAM)
# so alignment runs from FASTQ with a multi-unit merge, and realigns from both old references.
#
#   tests/stub/make_fixture.sh <dir>

set -euo pipefail

dir=$(mkdir -p "$1" && cd "$1" && pwd)
data="${dir}/data"
refs="${dir}/refs"

mkdir -p "${data}/cohorts/STUB" "${data}/fastq" "${data}/old_cram_38" "${data}/old_cram_37" \
    "${refs}/bundle/roh" "${refs}/wombat" "${refs}/annotation/gencode" "${dir}/scratch" \
    "${refs}/melt/me_refs/Hg38" "${refs}/melt/add_bed_files/Hg38"

printf 'FID\tbarcode\tfather\tmother\tsex\tphenotype\n' > "${data}/cohorts/STUB/STUB.pedigree.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    FAM01 BC001 0 0 1 1 \
    FAM01 BC002 0 0 2 1 \
    FAM01 BC003 BC001 BC002 1 2 \
    FAM02 BC004 0 0 2 2 \
    >> "${data}/cohorts/STUB/STUB.pedigree.tsv"

# <run>_<x>_<project>_<x>_<barcode>_<lane>_<flowcell>.<dual>.R{1,2}.fastq.gz, as main.nf parses it
for unit in BC001_L001 BC001_L002 BC002_L001; do
    for read in R1 R2; do
        touch "${data}/fastq/RUN_X_PRJ001_X_${unit}_FC01.dual.${read}.fastq.gz"
    done
done
touch "${data}/old_cram_38/BC003.cram" "${data}/old_cram_38/BC003.cram.crai"
touch "${data}/old_cram_37/BC004.cram" "${data}/old_cram_37/BC004.cram.crai"

for f in ref.fa ref.fa.fai old37.fa old38.fa bazam.jar vep.ini gnomad.vcf.gz \
         wisecondorx.100000.npz eh_thresholds.tsv \
         ancestry_catalog.tsv bundle/roh/sites.tsv.gz extract_me.tsv hg19ToHg38.over.chain.gz \
         annotation/gencode/gencode.stub.symbol.genes.bed.gz gnomad_sv.vcf.gz ndd_genes.tsv \
         melt/MELT.jar melt/me_refs/Hg38/ALU_MELT.zip melt/me_refs/Hg38/LINE1_MELT.zip \
         melt/me_refs/Hg38/SVA_MELT.zip melt/add_bed_files/Hg38/Hg38.genes.bed; do
    touch "${refs}/${f}"
done
# A CNV panel for the cohort fit: main.nf reads its panel.json before the run
mkdir -p "${refs}/cnv_panel/v0"
printf '{"panel_version": "0", "tag": "cov_b250q1"}\n' > "${refs}/cnv_panel/v0/panel.json"
echo stub > "${refs}/cnv_panel/v0/panel.cov_b250q1.autosomes.median.bw"
# The expansionhunter workflow reads the catalog itself, for the loci of multi-variant entries
cat > "${refs}/eh_catalog.json" <<'JSON'
[
  {"LocusId": "HTT", "LocusStructure": "(CAG)*", "ReferenceRegion": "chr4:3074876-3074933", "VariantType": "Repeat"},
  {"LocusId": "FXN", "LocusStructure": "(A)*(GAA)*", "ReferenceRegion": ["chr9:69037261-69037286", "chr9:69037286-69037304"],
   "VariantId": ["FXN_A", "FXN"], "VariantType": ["Repeat", "RareRepeat"]}
]
JSON
# An empty config runs the standard wombat path; this one adds the external MNV annotation
: > "${refs}/wombat/standard.yml"
printf 'mnv:\n  annotate: external\n' > "${refs}/wombat/mnv.yml"

cat > "${dir}/base.yml" <<YAML
data: "${data}"
scratch: "${dir}/scratch"
work_dir: "${dir}/work"
ref: "${refs}/ref.fa"
ref_name: "GRCh38_stub"
fastq_pattern: "*.R{1,2}.fastq.gz"
old_cram_38: "${data}/old_cram_38"
old_ref_38: "${refs}/old38.fa"
old_cram_37: "${data}/old_cram_37"
old_ref_37: "${refs}/old37.fa"
bazam: "${refs}/bazam.jar"
deepvariant_threads: 2
vep_config: "${refs}/vep.ini"
vep_config_name: "vep_stub"
gnomad_file: "${refs}/gnomad.vcf.gz"
annotation_annotation_path: "${refs}/annotation"
annotation_gencode: "gencode.stub"
wombat_config_path: "${refs}/wombat"
wombat_config_list: ["standard.yml", "mnv.yml"]
wisecondorx_binsize: 100000
wisecondorx_reference: "${refs}/wisecondorx.{binsize}.npz"
extractor_tsvs_list: ["${refs}/extract_me.tsv"]
liftover_chain: "${refs}/hg19ToHg38.over.chain.gz"
ancestry_catalog: "${refs}/ancestry_catalog.tsv"
ancestry_panel_name: "apgs_stub"
expansionhunter_catalog: "${refs}/eh_catalog.json"
expansionhunter_thresholds: "${refs}/eh_thresholds.tsv"
melt_dir: "${refs}/melt"
melt_gnomad_sv: "${refs}/gnomad_sv.vcf.gz"
melt_ndd_genes: "${refs}/ndd_genes.tsv"
max_cpus: 2
max_memory: "2.GB"
max_time: "1.h"
YAML
