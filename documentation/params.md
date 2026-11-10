# The parameters file (`params.yml`)

A cohort's parameters file is the only thing an operator writes to run this pipeline. It is
a flat YAML mapping — every key is a scalar or a list of scalars, with one exception,
[`resources`](#resources-per-cohort-cpus-memory-and-time) — and it is read by `nextflow run`
through `-params-file`.

This document is meant to be usable **next to a parameters file that has no comments**: it
lists every key the code actually reads, its default, when it is required, and what it
changes. It also lists the keys that appear in existing parameters files and do nothing, so
you can recognise them for what they are.

> **The pipeline does not validate key names or types** — except inside `resources`, which is
> checked in full when the run starts. There is no JSON schema. An unknown
> or misspelled key is accepted in silence, and the code falls back to the default in
> `nextflow.config` or to `null`. Nothing will tell you that `bin_size` is not `bin`. This is
> the single most important thing to know when reading a file someone else wrote.
>
> The corollary is the rule to read this document by: **if a key in your file is not listed
> in the [index](#every-key-a-to-z), it has no effect.** Several such keys are in circulation —
> they are listed in [Keys that look live but are inert](#keys-that-look-live-but-are-inert).

## Reading a parameters file you did not write

In this order:

1. **`steps`** — everything else is conditional on this. Keys belonging to a step that is not
   listed are inert for that run, whatever they say.
2. **`data` and `cohort_name`** — together they fix where every input is looked up and every
   output is written.
3. **`ref_name`, `vep_config_name`, `wisecondorx_binsize`, `bin`, `ancestry_panel_name`** —
   these are *baked into output filenames*. They tell you which existing results the run will
   find and reuse, and which it will ignore and recompute.
4. Then only the groups belonging to the steps in (1).

Anything left over is either a tuning knob with a sane default or, quite possibly, dead. Look
it up in the [index](#every-key-a-to-z): absent from it means it does nothing.

### What an empty string means

The example files are full of `""`, and it means three different things:

| Case | Meaning | Examples |
|---|---|---|
| A default meaning "use the tool's own default" | The flag is omitted from the command line | `ancestry_min_coverage`, `ancestry_model` |
| A default meaning "this input source is not in use" | The source is skipped | `fastq_pattern`, `old_cram_37`, `old_cram_38` |
| Unset, and nothing checks it | The run proceeds and fails later, or produces wrong paths | `ref`, `ref_name`, `data` (this one is checked) |

## Every key, A to Z

74 keys are read by the code; 17 more are in circulation and do nothing. Everything either
program can see is in this table. A key that is not here is not read by anything.

📄 marks a key whose value becomes part of output filenames — see
[Filename keys are the cache](#filename-keys-are-the-cache). ✗ marks an inert key.

| Key | Read when | Documented under |
|---|---|---|
| `ancestry_catalog` | `ancestry` — **required** | [Ancestry and PGS](#ancestry-and-polygenic-scores) |
| `ancestry_min_coverage` | `ancestry` | [Ancestry and PGS](#ancestry-and-polygenic-scores) |
| `ancestry_min_dp` | `ancestry` | [Ancestry and PGS](#ancestry-and-polygenic-scores) |
| `ancestry_min_gq` | `ancestry` | [Ancestry and PGS](#ancestry-and-polygenic-scores) |
| `ancestry_model` | `ancestry` | [Ancestry and PGS](#ancestry-and-polygenic-scores) |
| `ancestry_panel_name` 📄 | `ancestry`, `roh` — **required**; `upd`, to find the roh tables | [Ancestry and PGS](#ancestry-and-polygenic-scores) |
| `ancestry_reference` | `ancestry`, `roh` — **required**; `upd`, for its mask | [Ancestry and PGS](#ancestry-and-polygenic-scores) |
| `annotation_annotation_list` | `annotation` | [Annotation](#annotation) |
| `annotation_annotation_path` | `annotation`, `wisecondorx`, `roh` | [Annotation](#annotation) |
| `annotation_dnm_min_callrate` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `annotation_dnm_min_DP` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `annotation_dnm_min_GQ` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `annotation_dnm_min_VAF` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `annotation_gencode` | `wisecondorx`, `roh` | [Annotation](#annotation) |
| `apptainer_cache` | always, at config-parse time | [Containers and resources](#containers-and-resources) |
| `bazam` | `alignment`, realignment only | [Alignment tools and coverage](#alignment-tools-and-coverage) |
| `bin` 📄 | `alignment` | [Alignment tools and coverage](#alignment-tools-and-coverage) |
| `bin_size` ✗ | never — you want `bin` | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `bwa_mem2` | `alignment` | [Alignment tools and coverage](#alignment-tools-and-coverage) |
| `clean_stale_families` | optional — a repair flag, not a run setting | [Core paths and identity](#core-paths-and-identity) |
| `cnv_cov_bin` 📄 | `cnv_coverage`, `cnv_panel`; `sampleqc`, `upd`, to find the coverage JSONs | [CNV coverage and control panel](#cnv-coverage-and-control-panel) |
| `cnv_cov_mapq` 📄 | `cnv_coverage`, `cnv_panel`; `sampleqc`, `upd`, to find the coverage JSONs | [CNV coverage and control panel](#cnv-coverage-and-control-panel) |
| `cnv_cov_round` | `cnv_coverage`, `cnv_panel` | [CNV coverage and control panel](#cnv-coverage-and-control-panel) |
| `cnv_panel_dir` | `cnv_coverage` | [CNV coverage and control panel](#cnv-coverage-and-control-panel) |
| `cnv_panel_max_fit_mad` | `cnv_panel` | [CNV coverage and control panel](#cnv-coverage-and-control-panel) |
| `cnv_panel_version` 📄 | `cnv_panel` — **required** | [CNV coverage and control panel](#cnv-coverage-and-control-panel) |
| `cohort_name` 📄 | always — **required** | [Core paths and identity](#core-paths-and-identity) |
| `cram4realignment_pattern` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `data` | always — **required** | [Core paths and identity](#core-paths-and-identity) |
| `deepvariant_threads` | `deepvariant_sample` | [DeepVariant and GLnexus](#deepvariant-and-glnexus) |
| `enable_conda` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `expansionhunter_catalog` | `expansionhunter` — **required** | [ExpansionHunter](#expansionhunter) |
| `expansionhunter_catalog_name` 📄 | `expansionhunter` — **required** | [ExpansionHunter](#expansionhunter) |
| `expansionhunter_reviewer` | `expansionhunter` | [ExpansionHunter](#expansionhunter) |
| `expansionhunter_thresholds` | `expansionhunter` | [ExpansionHunter](#expansionhunter) |
| `extractor_tsvs_list` | `extractor` | [Extractor](#extractor) |
| `fastq_pattern` | `alignment` | [Input discovery](#input-discovery) |
| `glnexus_bed` | `deepvariant_family` | [DeepVariant and GLnexus](#deepvariant-and-glnexus) |
| `glnexus_config` | `deepvariant_family` | [DeepVariant and GLnexus](#deepvariant-and-glnexus) |
| `gnomad_file` | `annotation`; `upd` — **required** | [Annotation](#annotation) |
| `gnomad_filter_field` | `annotation` | [Annotation](#annotation) |
| `gnomad_filter_threshold` | `annotation` | [Annotation](#annotation) |
| `liftover_chain` | `extractor`, GRCh37 inputs | [Extractor](#extractor) |
| `max_cpus` | always, at config-parse time | [Containers and resources](#containers-and-resources) |
| `max_memory` | always, at config-parse time | [Containers and resources](#containers-and-resources) |
| `max_time` | always, at config-parse time | [Containers and resources](#containers-and-resources) |
| `melt_dir` | `melt` — **required** | [MELT](#melt) |
| `melt_gnomad_sv` | `melt` — **required** | [MELT](#melt) |
| `melt_label` 📄 | `melt` — **required** | [MELT](#melt) |
| `melt_ndd_columns` | `melt` | [MELT](#melt) |
| `melt_ndd_genes` | `melt` — **required** | [MELT](#melt) |
| `old_cram_37` | `alignment`, realignment only | [Reference genomes](#reference-genomes) |
| `old_cram_38` | `alignment`, realignment only | [Reference genomes](#reference-genomes) |
| `old_ref_37` | `alignment`, realignment only | [Reference genomes](#reference-genomes) |
| `old_ref_38` | `alignment`, realignment only | [Reference genomes](#reference-genomes) |
| `oldref` ✗ | never — you want `old_ref_38` | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `pedigree` | optional — defaults to the cohort's pedigree | [Core paths and identity](#core-paths-and-identity) |
| `pedigree_strict` | always | [Core paths and identity](#core-paths-and-identity) |
| `ref` | most steps | [Reference genomes](#reference-genomes) |
| `ref_name` 📄 | `alignment` and everything downstream | [Reference genomes](#reference-genomes) |
| `ref_par1_end` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `ref_par1_start` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `ref_par2_end` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `ref_par2_start` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `resources` | always, checked at startup | [`resources`](#resources-per-cohort-cpus-memory-and-time) |
| `roh_min_gene_roh_mb` | `roh` | [ROH and inbreeding](#roh-and-inbreeding) |
| `sambamba` | `alignment` | [Alignment tools and coverage](#alignment-tools-and-coverage) |
| `samblaster` | `alignment` | [Alignment tools and coverage](#alignment-tools-and-coverage) |
| `sampleqc_hom_cutoff` | `sampleqc` | [Sample QC](#sample-qc) |
| `sampleqc_label` 📄 | `sampleqc` — **required** | [Sample QC](#sample-qc) |
| `samtools` | `alignment` | [Alignment tools and coverage](#alignment-tools-and-coverage) |
| `scratch` | `alignment` | [Core paths and identity](#core-paths-and-identity) |
| `singularity_cache` | always, at config-parse time | [Containers and resources](#containers-and-resources) |
| `singularity_pull_docker_container` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `skip_alignment` ✗ | never — use `steps` | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `skip_deepvariant` ✗ | never — use `steps` | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `slurm_account` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `slurm_partition` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `somalier_sites` | `sampleqc` — **required** | [Sample QC](#sample-qc) |
| `steps` | always — **required** | [`steps`](#steps-the-key-that-decides-everything-else) |
| `upd_af_field` | `upd` | [UPD](#upd) |
| `upd_label` 📄 | `upd` | [UPD](#upd) |
| `vaf_bedgraph` | always — decides what counts as a called sample | [DeepVariant and GLnexus](#deepvariant-and-glnexus) |
| `vep_config` | `annotation` | [Annotation](#annotation) |
| `vep_config_name` 📄 | `annotation` and everything downstream | [Annotation](#annotation) |
| `wisecondorx_binsize` 📄 | `wisecondorx` | [WisecondorX](#wisecondorx) |
| `wisecondorx_predict_args` | `wisecondorx` | [WisecondorX](#wisecondorx) |
| `wisecondorx_reference` | `wisecondorx` | [WisecondorX](#wisecondorx) |
| `wombat_config_list` 📄 | `wombat` | [Wombat](#wombat) |
| `wombat_config_path` | `wombat` | [Wombat](#wombat) |
| `work_dir` | always, at config-parse time | [Core paths and identity](#core-paths-and-identity) |

## How a parameters file reaches the pipeline

A run is configured in two layers:

| Layer | Where | Notes |
|---|---|---|
| 1. Built-in defaults | the `params { }` block in [`nextflow.config`](../nextflow.config) | Covers 32 of the 50 keys. The other 18 are `null` unless you set them. |
| 2. The parameters file | `-params-file`, via a cohort name or `--params-file` | Where a cohort is configured. |

> **Configure everything in the parameters file.** Nextflow will also accept individual
> parameters on the command line, but do not use them. The parameters file is what gets
> committed, shared and checksummed, so a run configured entirely from it can be reproduced
> from it. A value passed on the command line leaves no trace in the file, and the cohort's
> state record will show a `params_effective_sha256` that does not match the file it names —
> which is exactly the situation that field exists to detect. If a setting needs to change,
> change it in the file.

### Giving the file by cohort name

```bash
ghfc-ngs CANDY_mpx
```

A bare **first** argument is a cohort name, and stands in for its parameters file. Two paths
are tried, in order:

1. `./cohorts/<NAME>/<NAME>.params.yml`, relative to the current directory;
2. `$GHFC_NGS_COHORTS/<NAME>/<NAME>.params.yml`, defaulting to
   `/pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/cohorts`.

A local copy therefore wins over the shared one. If neither exists the runner stops before
launching Nextflow and prints both paths it tried.

### Giving the file by path

```bash
ghfc-ngs --params-file my_params.yml
```

`--params-file` takes precedence over a cohort name given on the same command line, and the
cohort name is ignored with a warning.

### Keys read before the run starts

Most keys are read while the workflow executes. These are read earlier, while Nextflow parses
its configuration, and so cannot be changed by anything that happens during the run:

- `work_dir` — becomes Nextflow's `workDir`, falling back to `work`
- `apptainer_cache`, `singularity_cache` — the container cache directories
- `max_memory`, `max_cpus`, `max_time` — the ceilings applied to every process

## `steps`: the key that decides everything else

```yaml
steps: ["alignment", "deepvariant_sample", "deepvariant_family", "annotation", "wombat", "wombat_cohort", "snvs_cohort"]
```

Required, and must be non-empty. A list is the form to use; a comma-separated string
(`steps: "alignment,wombat"`) is accepted and split, which is also what makes `--steps` work on
the command line. Every entry must be one of these seventeen:

| Step | Produces |
|---|---|
| `alignment` | CRAMs from FASTQ, or realigned from existing GRCh37/GRCh38 CRAMs, plus coverage bedgraphs |
| `deepvariant_sample` | Per-individual gVCF and VCF, plus VAF bedgraphs |
| `deepvariant_family` | Joint family calls via GLnexus, normalised, plus per-family pedigrees |
| `annotation` | gnomAD frequencies, the rare/common split, VEP, and extra bcftools annotations |
| `wombat` | PyWombat filtering of the annotated BCF, one TSV per config |
| `wombat_cohort` | Every family's wombat TSV concatenated into one cohort table per config |
| `snvs_cohort` | Cohort-level merge of every family's common variants into one BCF |
| `wisecondorx` | CNV/SV calling, per sample then merged per family and per cohort |
| `extractor` | Specific variants pulled out of family BCFs, wombat TSVs or gVCFs |
| `ancestry` | Ancestry components, admixture and polygenic scores, per family and per cohort |
| `roh` | Runs of homozygosity and inbreeding (F_ROH) per individual, genes in ROH, per family and per cohort |
| `expansionhunter` | Short tandem repeat genotypes per sample, merged and classed per family, merged and summarised per cohort |
| `melt` | Mobile element insertions (ALU, LINE1, SVA) per sample, merged per family, with their cohort frequency and a ranking of every insertion per cohort |
| `cnv_coverage` | Normalised coverage per sample on a fixed bin grid, as a bigWig for genetics-viz, and with `cnv_panel_dir` the cohort's fit against a control panel |
| `cnv_panel` | A versioned control panel, per-bin median and MAD of the controls' normalised coverage, built once from a dedicated controls cohort |
| `sampleqc` | Relatedness (IBS0, IBS2, KING kinship) against the pedigree within each family and between families, sex and sex-chromosome karyotype, and Mendelian errors, as per-family and cohort tables of flags |
| `upd` | Uniparental disomy segments of every member, with coordinates, type (isodisomy, heterodisomy or mixed), parent of origin and a score, from the family BCF with both parents, one, or none (isodisomy only, from the `roh` tables), per family and per cohort |

### The rule that catches people out

The pipeline looks at what already exists on disk and only does the missing work. But it will
**not silently run a step you did not ask for**: if work is needed for a step that is not in
`steps`, the run stops.

There is one refinement worth knowing, because it is what makes partial runs possible: a
missing prerequisite is only an error when something in *this* run would consume it.

- `alignment` and `deepvariant_sample` are checked unconditionally — if any individual needs
  them and they are not listed, the run fails.
- `deepvariant_family` is only demanded when one of `deepvariant_family`, `annotation`,
  `wombat`, `wombat_cohort`, `snvs_cohort`, `extractor` or `upd` is listed.
- `annotation` is only demanded when one of `annotation`, `wombat`, `wombat_cohort` or
  `snvs_cohort` is listed.
- `wombat` is only demanded when `wombat_cohort` is listed and a cohort table is due, which it
  is whenever a family has no wombat tables yet. A cohort table is only ever merged from every
  family's table, never from part of the cohort.
- `ancestry` is demanded by `roh` for families with no family panel or ancestry table on disk
  under the current `ancestry_panel_name`.
- `expansionhunter` demands nothing. Like `wisecondorx`, it reads the CRAMs, and a sample with
  no CRAM is left out of its family's files, with a warning, rather than sent to alignment.
- `melt` demands nothing either, for the same reason.
- `cnv_coverage` demands nothing either, for the same reason. `cnv_panel` demands
  `cnv_coverage` when the panel is still to be built and a sample with a CRAM has no coverage
  file yet: a panel is only ever built from every sample of the cohort.
- `sampleqc` demands nothing either: it reads the CRAMs, and uses the cnv_coverage JSONs and
  the family `norm.bcf` when they exist or are made in the run, without them otherwise.
- `upd` demands `deepvariant_family`, as above. It uses the cnv_coverage JSONs and the `roh`
  tables when they exist or are made in the run, and does without them otherwise: a member
  with no parent in the family BCF is then reported as skipped, not an error.

`snvs_cohort` used to merge the wombat tables too. A parameters file that lists `snvs_cohort`
and not `wombat_cohort` gets a warning at startup, because its cohort wombat tables are no
longer produced.

Every cohort-level file is merged again, from every family, whenever a family that feeds it is
built or rebuilt in the run: the common-variant BCF when a family is annotated, the cohort
wombat tables when a family's wombat tables are written, the cohort WisecondorX BED when a
family's aberrations are annotated, as the ancestry, roh and expansionhunter cohort files
already were. Adding a family to the pedigree is therefore enough to have it in every cohort
file. On a cohort the size of SPARK that means the common-variant merge is redone for each
batch of families added; leave `snvs_cohort` out of `steps` until the last batch to merge it
once. A merge that is skipped, because its step is not listed or the run stopped before it,
is not caught up by a later run that rebuilds no family: delete the cohort file to have it
merged again.

This is why `steps: ["ancestry"]` runs happily on a cohort whose annotation is incomplete:
ancestry reads gVCFs directly and never touches the normalised or annotated call sets. The same
holds for `steps: ["ancestry", "roh"]`, since `roh` reads only what `ancestry` writes.

### Failures you will see

| Message | Cause |
|---|---|
| `ERROR: --data parameter is required` | `data` unset or empty |
| `ERROR: cohort_name parameter is required` | `cohort_name` unset or empty |
| `ERROR: --steps parameter is required` | `steps` missing or `[]` |
| `ERROR: Invalid steps specified: ...` | a value outside the seventeen above; the message names the offending entries |
| `ERROR: the 'ancestry' step requires ...` | `ancestry` listed without its three required keys |
| `ERROR: the 'roh' step requires ...` | `roh` listed without `ancestry_reference` and `ancestry_panel_name` |
| `ERROR: the 'expansionhunter' step requires ...` | `expansionhunter` listed without `expansionhunter_catalog`, `expansionhunter_catalog_name` or `ref` |
| `ERROR: the 'expansionhunter' step cannot find ...` | the catalog, the thresholds table or the reference's `.fai` does not exist |
| `ERROR: expansionhunter_catalog_name '...' may only hold ...` | a label with a character that cannot go in a file name |
| `ERROR: the 'melt' step requires ...` | `melt` listed without `melt_dir`, `melt_label`, `melt_gnomad_sv`, `melt_ndd_genes`, `ref`, or the gencode settings `annotation_annotation_path` and `annotation_gencode` |
| `ERROR: the 'melt' step cannot find ...` | a file of the MELT distribution, gnomAD-SV, the NDD table, the gencode gene BED or the reference's `.fai` does not exist |
| `ERROR: melt_label '...' may only hold ...` | as for `expansionhunter_catalog_name` |
| `ERROR: cnv_cov_bin '...' must be ...`, and the same for `cnv_cov_mapq`, `cnv_cov_round` | `cnv_coverage` or `cnv_panel` listed with a bin, MAPQ or rounding that is not a number of the right kind |
| `ERROR: the reference ... or its .fai does not exist` | `cnv_coverage` or `cnv_panel` listed, and `ref` or its `.fai` is missing |
| `ERROR: the 'sampleqc' step requires ...` | `sampleqc` listed without `somalier_sites`, `sampleqc_label` or `ref` |
| `ERROR: the 'sampleqc' step cannot find ...` | the sites VCF, the reference or its `.fai` does not exist |
| `ERROR: sampleqc_label '...' may only hold ...` | as for `expansionhunter_catalog_name` |
| `ERROR: the 'upd' step requires ...` | `upd` listed without `upd_label`, `upd_af_field` or `gnomad_file` |
| `ERROR: upd_label '...' may only hold ...` | as for `expansionhunter_catalog_name` |
| `ERROR: upd_af_field '...' is not an INFO field name` | a field name with a character other than letters, digits and `_` |
| `ERROR: the 'upd' step cannot find gnomad_file ...` | the gnomAD file, or its `.csi` or `.tbi` index, does not exist |
| `ERROR: ... has no INFO/<field> (upd_af_field)`, in a `UPD_FAMILY` task | `upd_af_field` names a field the gnomAD file does not have |
| `ERROR: the 'cnv_panel' step requires cnv_panel_version to be set` | `cnv_panel` listed without `cnv_panel_version` |
| `ERROR: cnv_panel_version '...' may only hold ...` / `... is written without its v` | a version that cannot be a directory name, or `v1` where `1` was meant |
| `ERROR: the CNV panel directory ... exists but is incomplete` | a partial `cnv_panel/v<version>/`, from a build that died or a hand edit: it is never overwritten |
| `ERROR: cnv_panel_dir: ...` | `cnv_panel_dir` set to a directory with no readable `panel.json` or no autosomal median bigWig |
| `ERROR: Pedigree file not found: ...` | `pedigree`, or the cohort's default pedigree path, does not exist |
| `... is required for N individuals but not included in steps parameter` | the rule above |

Two of those messages name `--data` and `--steps`. Both do exist as command-line flags, but the
keys they are telling you about are `data:` and `steps:` in the parameters file, which is where
they should be set.

## Key reference

74 keys are read by the code. In the tables below, **no default** means the key is `null` when
omitted — see [the `null` trap](#keys-with-no-default-become-the-string-null). A 📄 marks a key
whose value is **part of output filenames**.

### Core paths and identity

| Key | Type | Default | Needed by | Effect |
|---|---|---|---|---|
| `data` | path | `""` | **always** | Root of the cohort tree. Prefix of nearly every input lookup and every published output. |
| `cohort_name` 📄 | string | *no default* | **always** | Names `<data>/cohorts/<cohort_name>/`, every cohort-level file, and the default pedigree path. |
| `scratch` | path | `""` | `alignment` | Per-unit CRAMs are *published* to `<scratch>/cram`, and sorting uses `<scratch>/tmp`. Despite the name, this holds output. |
| `work_dir` | path | *no default* → `work` | — | Nextflow's work directory. Read at config-parse time. |
| `pedigree` | path | *no default* → the cohort's own pedigree | rarely | Only needed for a pedigree kept outside the cohort directory. See below. |
| `pedigree_strict` | boolean | `false` | — | Stop the run on stale family outputs instead of warning. |
| `clean_stale_families` | boolean | `false` | — | Delete stale family outputs so they are re-called, instead of only reporting them. See below. |
| `steps` | list of strings | `['alignment']` | **always** | See [above](#steps-the-key-that-decides-everything-else). |

`cohort_name` is **required**: the run aborts immediately without it. It names the cohort
directory, is baked into every cohort-level filename, is where the run-state file is written,
and supplies the default pedigree path — so a missing value would otherwise produce
`cohorts/null/null.*` outputs and record nothing about the run.

`pedigree` defaults to the cohort's own pedigree, by convention:

```
<data>/cohorts/<cohort_name>/<cohort_name>.pedigree.tsv
```

This is the same directory that holds the cohort's parameters file and its outputs, so for a
cohort laid out the usual way the key can be left out entirely. Set it only for a pedigree kept
somewhere else. Either way the file must exist, or the run aborts naming the path it tried and
saying which of the two it was.

The list of samples is **not** in the parameters file — it comes from the pedigree.

#### `clean_stale_families`: repairing a drifted cohort

A family whose pedigree gains a member who was never called has *stale* outputs: the joint call
on disk cannot contain that member, and the run reports this as `STALE FAMILY OUTPUTS`. It does
not fix itself — the presence of the family's `norm.bcf` is exactly what stops the missing
member being scheduled for anything.

`clean_stale_families` deletes those families' derived outputs so they are re-called in full:
the `norm.bcf` and its index, the family pedigree, every annotation output, and the wombat
parquet and result TSVs — plus the cohort-level merges built from them. It leaves the
`ancestry/` outputs alone, because those already rebuild themselves, and it never touches
gVCFs, CRAMs or `extractor/` output.

It is a repair flag, not a run setting. Pass it on the command line for the one run that needs
it rather than putting it in a parameters file, where it would silently re-delete on every
future run:

```bash
ghfc-ngs MY_COHORT --clean-stale-families --dry-run   # print every path it would remove
ghfc-ngs MY_COHORT --clean-stale-families             # then actually do it
```

**It refuses families it cannot rebuild.** For each drifted family it checks that every missing
member can reach a gVCF — the member already has a CRAM, or an alignment input resolves for it
in `fastq_pattern` / `old_cram_37` / `old_cram_38` — *and* that `steps` contains the steps that
would do the work. A family failing either check is left completely untouched and reported with
the reason. This is the point of the flag: deleting a call set that cannot be rebuilt would
leave the cohort worse off than the warning it replaces.

Two things to know about how it is spelled:

- **`--dry-run` is a real dry run.** It maps to Nextflow's `-preview`, and the clean prints
  every path it would remove without deleting anything. `-stub-run` is inert in the same way.
- **The hyphenated form only works through the `ghfc-ngs` wrapper.** Nextflow maps a
  command-line `--clean-stale-families` to `params.cleanStaleFamilies`, never to
  `params.clean_stale_families`; the wrapper translates it. Calling `nextflow run` directly,
  use `--clean_stale_families`. Either spelling is read, so neither silently does nothing.
- **Do not write `--clean-stale-families false`.** A value given to a flag on the command line
  arrives as the *string* `"false"`, which is true. The flag reads `false`, `no`, `0` and the
  empty string as off wherever it can, but the reliable way to not clean is to omit the flag.

### Reference genomes

| Key | Type | Default | Needed by | Effect |
|---|---|---|---|---|
| `ref` | path | `""` | alignment, both deepvariant steps, wombat, wisecondorx | The GRCh38 FASTA. |
| `ref_name` 📄 | string | `""` | alignment, anything reading CRAMs | Infix in CRAM names: `<barcode>.<ref_name>.cram`. |
| `old_ref_38` | path | `""` | realignment from legacy GRCh38 CRAMs | Reference the old CRAMs were aligned to, for bazam re-extraction. |
| `old_cram_38` | dir | `""` | as above | Flat directory holding those CRAMs. |
| `old_ref_37` | path | `""` | realignment from GRCh37 CRAMs | As `old_ref_38`, for GRCh37. |
| `old_cram_37` | dir | `""` | as above | Flat directory holding those CRAMs. |

### Input discovery

| Key | Type | Default | Needed by | Effect |
|---|---|---|---|---|
| `fastq_pattern` | glob | `""` | alignment from FASTQ | Glob under `<data>/fastq/`, read as paired-end pairs. |

Alignment needs **at least one** of `fastq_pattern`, `old_cram_37`, `old_cram_38`. With none of
them set and work to do, the run stops rather than quietly aligning nothing. Individuals the
plan schedules but cannot find an input for are named, and the run fails.

Nothing, however, checks that a CRAM directory is paired with its reference: setting
`old_cram_37` without `old_ref_37` realigns against an empty reference path instead of being
rejected. Set the two together.

### Alignment tools and coverage

| Key | Type | Default | Effect |
|---|---|---|---|
| `bwa_mem2` | command | `"bwa-mem2"` | Aligner executable. |
| `samblaster` | command | `"samblaster"` | Duplicate marking. |
| `sambamba` | command | `"sambamba"` | Sorting. |
| `samtools` | command | `"samtools"` | samtools, including indexing. |
| `bazam` | path | `"/path/to/bazam.jar"` | The bazam JAR, for re-extracting FASTQ from CRAM. The default is a placeholder: set it if you realign. |
| `bin` 📄 | integer | `1000` | Coverage bin size, and part of the bedgraph name `<barcode>.by<bin>.bedgraph.gz`. |

### DeepVariant and GLnexus

| Key | Type | Default | Effect |
|---|---|---|---|
| `deepvariant_threads` | integer | `96` | `--num_shards` for the three DeepVariant stages. |
| `glnexus_config` | string | `"DeepVariant_unfiltered"` | GLnexus preset for joint calling. `DeepVariant` applies quality filters; `gatk` is for GATK gVCFs. |
| `glnexus_bed` | path | `""` | BED of the ranges to joint-call, passed to GLnexus as `--bed`. Empty calls the whole genome. For exome cohorts. |
| `vaf_bedgraph` | boolean | `true` | Whether a sample's VAF bedgraph is one of its `deepvariant_sample` outputs. Set `false` for a cohort that arrives as gVCFs only. |

`deepvariant_threads` is **independent of `max_cpus` and of the SLURM allocation**. It exists
because DeepVariant's throughput is hardware-sensitive, and it overrides the process CPU
setting without changing what SLURM actually reserves. Setting it above the cores you were
allocated oversubscribes the node. Keep it, `max_cpus`, and the `withName` blocks in
`nextflow.config` consistent by hand.

**Exome cohorts and `glnexus_bed`.** A DeepVariant gVCF covers the whole genome with reference
blocks, so without a BED every off-target call reaches annotation and wombat. Give the capture
targets padded by the distance you want to keep around them. GLnexus **rejects overlapping
ranges**, so the BED must be sorted and merged after padding, and if the cohort was captured
with several kits, use the union:

```bash
bedtools slop -i targets.bed -g ref.fasta.fai -b 100 | sort -k1,1 -k2,2n | bedtools merge > targets.pad100.bed
```

The BED is read in place rather than staged, so leaving it empty keeps every existing cohort's
GLnexus tasks exactly as they were. Changing it does not re-call families that already have a
`norm.bcf`: like every other setting, it only applies to what is still to be done.

**gVCF-only cohorts and `vaf_bedgraph`.** The VAF bedgraph is built from the per-sample VCF that
DeepVariant writes next to the gVCF. A cohort imported as gVCFs alone has neither that VCF nor a
CRAM, so with the default a sample can never count as called: the plan sends every one of them
back to alignment and the run stops with "Alignment step is required". With `vaf_bedgraph: false`
a gVCF and its `.tbi` are the whole of `deepvariant_sample`'s output, and a pedigree member with
no gVCF still fails the run, by name. See [`params_example/SPARK.params.yml`](../params_example/SPARK.params.yml).

### Annotation

| Key | Type | Default | Effect |
|---|---|---|---|
| `gnomad_file` | path | *no default* | gnomAD BCF used to add population frequencies. |
| `gnomad_filter_field` | string | *no default* | The INFO field to split on, e.g. `AF` or `AF_genomes`. |
| `gnomad_filter_threshold` | string | *no default* | Frequency cutoff for the rare/common split, e.g. `"0.01"`. |
| `vep_config` | path | *no default* | The VEP `.ini` configuration file. |
| `vep_config_name` 📄 | string | *no default* | Infix in every VEP-derived filename. Appears in more than twenty path constructions. |
| `annotation_annotation_path` | dir | *no default* | Directory holding the extra annotation files. **Also used by the `wisecondorx` step.** |
| `annotation_annotation_list` | list of filenames | *no default* | Files inside that directory, applied with `bcftools annotate`. Names only, not paths. |
| `annotation_gencode` | string | *no default* | Gencode basename used to annotate **WisecondorX** aberrations with genes and exons. |

Despite its prefix, `annotation_gencode` belongs to the `wisecondorx` and `roh` steps, not to `annotation`.
`annotation_annotation_path` is shared by both.

### WisecondorX

| Key | Type | Default | Effect |
|---|---|---|---|
| `wisecondorx_reference` | path template | *no default* | Reference NPZ. The literal text `{binsize}` in it is replaced with `wisecondorx_binsize`. |
| `wisecondorx_binsize` 📄 | integer | *no default* | Bin size, passed to the NPZ conversion and substituted into the path above. |
| `wisecondorx_predict_args` | string | *no default* | Extra flags passed verbatim to `wisecondorx predict`, e.g. `"--zscore 5 --minrefbins 50"`. |

`wisecondorx_reference` and `wisecondorx_binsize` are two halves of one setting:

```yaml
wisecondorx_reference: ".../reference_wisecondorX_hg38_{binsize}_refsize100.npz"
wisecondorx_binsize: 1000
```

resolves to `reference_wisecondorX_hg38_1000_refsize100.npz`. Change the bin size and you are
pointing at a different file, which must exist.

### Wombat

| Key | Type | Default | Effect |
|---|---|---|---|
| `wombat_config_path` | dir | *no default* | Directory of PyWombat YAML configs. |
| `wombat_config_list` 📄 | list of filenames | *no default* | One output TSV per entry. Basenames, resolved under `wombat_config_path`. Each **stem becomes part of its output filename**, so renaming a config orphans its previous results. |

### Extractor

| Key | Type | Default | Effect |
|---|---|---|---|
| `extractor_tsvs_list` | list of paths | *no default* | Variant TSVs to extract. Full paths, unlike the two lists above. An empty list leaves the step with nothing to do. |
| `liftover_chain` | path | *no default* | Chain file used when an input TSV is on GRCh37. |

### Ancestry and polygenic scores

| Key | Type | Default | Effect |
|---|---|---|---|
| `ancestry_reference` | dir | `""` | The ancestry-pgs reference bundle. **Required** when `ancestry` is in `steps`. |
| `ancestry_catalog` | path | `""` | PGS weight catalog, SbayesRC layout. **Required** when `ancestry` is in `steps`. |
| `ancestry_panel_name` 📄 | string | `""` | Label baked into every ancestry filename. **Required** when `ancestry` is in `steps`. |
| `ancestry_min_dp` | integer | `10` | Minimum DP (variant records) / MIN_DP (reference blocks) to call a panel site. |
| `ancestry_min_gq` | integer | `20` | Minimum GQ to call a panel site. |
| `ancestry_min_coverage` | string | `""` | Lowers admixture's 0.90 per-sample coverage floor. Empty omits the flag. |
| `ancestry_model` | dir | `""` | A fitted z-score model outside the bundle. Empty uses the bundle's own. |

These three are the only keys in the whole file with a **hard requirement check**: listing
`ancestry` in `steps` without `ancestry_reference`, `ancestry_catalog` and
`ancestry_panel_name` aborts the run before any work starts. Listing `roh` without
`ancestry_reference` and `ancestry_panel_name` does too.

### ROH and inbreeding

| Key | Type | Default | Effect |
|---|---|---|---|
| `roh_min_gene_roh_mb` | number | `1.0` | Only ROH at least this long are intersected with genes. |

`roh` has no keys of its own beyond this one:
- **Reference.** It calls against `ancestry_reference`, which must carry the roh/ tier (ancestry-pgs `bundle/add_roh_tier.py`; bundle ≥ 1.1.0).
- **Inputs.** It reads the ancestry step's family panel and ancestry table, named by `ancestry_panel_name`.
- **Output names.** Its outputs carry `ancestry_panel_name` too. Adding the roh/ tier changes no ancestry result, so it needs no new label. To recompute ROH alone, delete the families' `roh/` directories.
- **Gene table.** It is written only when `annotation_gencode` and `annotation_annotation_path` are both set, from `<annotation_annotation_path>/gencode/<annotation_gencode>.{symbol,ensg}.genes.bed.gz`, the same files `wisecondorx` annotates with.

### ExpansionHunter

| Key | Type | Default | Effect |
|---|---|---|---|
| `expansionhunter_catalog` | file | `""` | ExpansionHunter variant catalog (JSON). **Required** by `expansionhunter`. |
| `expansionhunter_catalog_name` 📄 | string | `""` | Label for that catalog in every output file name. **Required** by `expansionhunter`. Letters, digits, `.`, `_` and `-` only. |
| `expansionhunter_thresholds` | file | `""` | Per-locus thresholds table, used to class each allele. Empty leaves the alleles unclassified. |
| `expansionhunter_reviewer` | boolean | `true` | Draw a REViewer read pileup of every genotype of severity high or medium. Needs `expansionhunter_thresholds`, without which nothing is flagged. |

- **The label is the cache.** Nothing on disk records which catalog a call was made with. A new catalog with an existing label reuses the old calls as if nothing had changed, and a new label recomputes every sample. The catalog of the cohort run, `references/GRCh38/expansionhunter/variant_catalog_custom_JMSA_offtargets.json`, goes with `jmsa50`.
- **Thresholds.** A TSV keyed on the catalog's VariantId, with the columns `locus`, `disease`, `inheritance`, `normal_max`, `intermediate_min` and `pathogenic_min` (`.` for none), and two optional ones: `premutation_min`, a class between intermediate and pathogenic (FMR1's 55-200, which causes FXTAS but not fragile X), and `report`, `flag` (the default) or `risk` for a common risk allele that is never more than `info` (TCF4's for Fuchs dystrophy, NIPA1's for ALS). `references/GRCh38/expansionhunter/str_thresholds.stripy.tsv` holds the STRipy values for the cohort catalog, and its README says how they were converted. It is in no file name, but its checksum is in each family's `inputs.tsv`: after changing it, the next run rebuilds every family and the cohort from the existing calls.
- **Severity.** Each genotype row is `high` (pathogenic status), `medium` (premutation) or `info` (intermediate, carrier, or a `risk` locus), one rank lower when an allele past normal is of low confidence: not spanned by 2 reads though short enough to be, longer with fewer than 3 in-repeat reads, at a degenerate motif (GCN) without 2 spanning reads, or at a locus whose reference allele the sample's reads cannot span and whose pathogenic range starts within two read lengths. Older genomes with 100 bp reads (the IP* samples of GHFC-GRCh38) are where most low-confidence calls come from. `flags.tsv` holds the `high` and `medium` rows; the others stay in `genotypes.tsv`. See the README's ExpansionHunter section for every column.
- **Reference and sex.** Each sample is genotyped against `ref`, which must be the reference its CRAM was aligned to and have a `.fai`. The sex comes from the pedigree; unknown sex is genotyped as female, which only matters at chrX loci.
- **REViewer pileups.** A genotype is flagged when its `severity` is `high` or `medium`: the rows of the cohort's `flags.tsv`. Each flagged locus gets `<ID>.<label>.<LocusId>.svg` in the sample's `str/expansionhunter/reviewer/`, beside `<ID>.<label>.metrics.tsv` (REViewer's allele depths) and `<ID>.<label>.phasing.tsv`. The metrics table is written even when nothing is flagged, and is what marks a sample done: a sample is drawn again when it is called or its family table rebuilt, so after a thresholds change the rebuilt tables redraw their members. SVGs of loci no longer flagged are left in place. Calls made before the step are drawn from the realigned reads beside them, which must be sorted and indexed; one without them is reported and skipped.
- **Calls made before the step.** `scripts/expansionhunter/eh_label_existing` gives unlabelled calls (`<ID>.vcf`, `<ID>.json`, `<ID>_realigned.bam`) the label, after checking that each VCF holds exactly the catalog's variants. Without `--apply` it only reports.

### MELT

| Key | Type | Default | Effect |
|---|---|---|---|
| `melt_dir` | directory | `""` | The MELT v2.2.2 distribution: `MELT.jar`, `me_refs/Hg38/{ALU,LINE1,SVA}_MELT.zip`, `add_bed_files/Hg38/Hg38.genes.bed`. **Required** by `melt`. |
| `melt_label` 📄 | string | `""` | Label for the calls in every output file name. **Required** by `melt`. Letters, digits, `.`, `_` and `-` only. |
| `melt_gnomad_sv` | file | `""` | gnomAD-SV sites VCF, for the population frequency of each insertion (`GRPMAX_AF`, else `AF`, within 50 bp). **Required** by `melt`. |
| `melt_ndd_genes` | file | `""` | NDD gene table: a TSV with a `Symbol` column and columns marking genes `True`. **Required** by `melt`. |
| `melt_ndd_columns` | string | `"High Confidence NDD genes v3"` | The columns of `melt_ndd_genes` that mark a gene NDD, `\|`-separated; a gene marked in any of them is. |

- **MELT is not in any image.** Its licence (University of Maryland, academic use) is non-transferable, so it cannot be redistributed; the image `fcliquet/melt-deps` holds what it runs with, and each task stages `melt_dir`. The cohort's copy is `/pasteur/helix/projects/ghfc_wgs/tools/MELTv2.2.2`.
- **The label is the cache**, as for ExpansionHunter. Give a new label to a new MELT version or a new set of mobile elements. The cohort's first calls, made with MELT 2.2.2 Single on ALU, LINE1 and SVA, go with `melt222`.
- **Reference.** Each sample is called against `ref`, which must be the reference its CRAM was aligned to. A CRAM whose contigs are not all in `ref` (aligned to `1` rather than `chr1`, say) fails at once and is not retried.
- **Genes.** The cohort scoring places insertions in genes with `<annotation_annotation_path>/gencode/<annotation_gencode>.symbol.genes.bed.gz`, the BED the `roh` step uses, so `melt` needs those two keys set too.
- **Scoring.** The cohort's candidates table ranks every PASS insertion by `2(1 - gnomAD AF) + 2(1 - cohort AF) + impact (1 exonic, 0.5 genic) + 2 NDD + 0.5 ASSESS/5`, the weights of the cohort's first MEI analysis. None of these settings is in a file name: the cohort files are rebuilt whenever a sample is called, or after you delete them.
- **Calls made before the step.** `scripts/melt/melt_import_existing` copies the calls of the first analysis (`<ID>/meis/{ALU,LINE1,SVA}.final_comp.vcf`) to `mei/melt/<ID>.<label>.<CLASS>.vcf`, after checking each one. Without `--apply` it only reports.

### CNV coverage and control panel

| Key | Type | Default | Effect |
|---|---|---|---|
| `cnv_cov_bin` 📄 | integer | `250` | Bin size in bp. Part of the tag `cov_b<bin>q<mapq>` in every file name. |
| `cnv_cov_mapq` 📄 | integer | `1` | Minimum mapping quality of the reads counted (mosdepth `-Q`). Part of the tag. |
| `cnv_cov_round` | number | `0.01` | Precision the normalised values are stored at. Coarser makes smaller bigWigs. Not in the tag. |
| `cnv_panel_dir` | directory | `""` | A panel version directory, `.../cnv_panel/v<version>`, to fit the cohort against. Empty: no fit. |
| `cnv_panel_version` 📄 | string | `""` | The version `cnv_panel` builds, as `cohorts/<cohort_name>/cnv_panel/v<version>/`. **Required** by `cnv_panel`. Letters, digits, `.`, `_` and `-`, without a leading `v`. |
| `cnv_panel_max_fit_mad` | number | `0.3` | Controls whose fit MAD is above this are left out of the panel, which is then computed again. |

genetics-viz reads these files directly, so their names, keys and values are a contract with it.

- **Per sample** (`cnv_coverage`, any cohort): `samples/{S1}/{S2}/<ID>/coverage/<ID>.cov_b<bin>q<mapq>.bw` and `.json`. mosdepth bins the CRAM (`-n --fast-mode --by <bin> -Q <mapq>`, its default flag filter, so duplicates, secondary and QC-fail reads are out). Only chr1-22, chrX and chrY are kept, and every bin is written, zeros included, on mosdepth's grid: from 0, the last bin of each contig truncated. A bin holds `2 × depth / autosomal_median`, rounded to `cnv_cov_round`, where the autosomal median is that of the chr1-22 bins with depth above 0. So a normal autosome reads about 2, a heterozygous deletion about 1, a homozygous one 0. The JSON gives the median, `frac_autosomal_bins_zero`, and `chrX_ratio` and `chrY_ratio`, about 1/0 for XX and 0.5/0.5 for XY: the median normalised value over non-PAR chrX bins with coverage, and over the non-PAR chrY bins that are less than half N in the reference, each halved.
- **The panel** (`cnv_panel`, only in the controls cohort's file): `cohorts/<cohort_name>/cnv_panel/v<version>/` holds `panel.cov_b<bin>q<mapq>.<group>.<stat>.bw` for the groups `autosomes`, `XX` and `XY` and the stats `median` and `mad` (unscaled), with `panel.json` and `controls.tsv`. Each control is fitted against the median of the others (the MAD of the log2 ratio over the autosomes). Controls above `cnv_panel_max_fit_mad` are excluded, and the panel is computed once more without them. A control enters `XX` or `XY` only when its pedigree sex and its coverage sex agree; a mismatch, an aneuploidy or an unknown sex keeps it in `autosomes` only. `controls.tsv` says what happened to each one, and why.
- **A panel version is built once.** A complete `v<version>/` counts as done and is never rebuilt or overwritten, even if the cohort has gained samples since (the run warns about them): set a new `cnv_panel_version` for that. A partial one stops the run. The directory is self-contained, to be referenced as it is by other cohorts and genetics-viz.
- **The fit** (`cnv_coverage` with `cnv_panel_dir`): `cohorts/<cohort_name>/coverage/<cohort_name>.cov_b<bin>q<mapq>.panel_v<version>.fit.tsv`, one row per family and sample, with its fit MAD against the panel, sex and ratios. It is rebuilt whenever a sample's coverage is. A panel of another tag is skipped, with a warning.
- **The tag is the cache.** Changing `cnv_cov_bin` or `cnv_cov_mapq` makes new files beside the old ones. `cnv_cov_round` is not in the tag, so changing it leaves the existing samples as they are: delete their `coverage/` files to have them written again.
- **Reference.** The coverage is measured against `ref`, which must be the reference the CRAMs were aligned to and have a `.fai`. The PAR coordinates are those of GRCh38.

### Sample QC

| Key | Type | Default | Effect |
|---|---|---|---|
| `somalier_sites` | file | `""` | somalier's sites VCF for `ref`: for GRCh38, `sites.hg38.vcf.gz` from somalier's README (17,384 autosomal SNPs, 365 on chrX, 17 on chrY). **Required** by `sampleqc`. |
| `sampleqc_label` 📄 | string | `""` | Label for the extracts and tables in every output file name. **Required** by `sampleqc`. Letters, digits, `.`, `_` and `-` only. One label per sites file, e.g. `somalier_hg38`. |
| `sampleqc_hom_cutoff` | number | `0.1` | `somalier relate` calls a site homozygous when its alternate allele fraction is below this, or above 1 minus it (somalier's `SOMALIER_AB_HOM_CUTOFF`, 0.01 by default). Above 0 and below 0.3, where hets start. Changing it redoes every family and the cohort. |

The step is report-only: nothing downstream reads its flags. Read the `flags.tsv` tables, fix the
pedigree, and the next run redoes the families that changed.

- **Per sample:** `samples/{S1}/{S2}/<ID>/qc/somalier/<ID>.<label>.somalier`, somalier's counts at its sites, from the CRAM, independent of the variant calls. The sample name inside is the barcode, whatever the CRAM's read group says.
- **Per family:** `families/{S1}/{S2}/<FID>/qc/<FID>.<label>.{samples,pairs,flags,inputs}.tsv`.
  - `pairs.tsv` has one row per pair of members, with IBS0, IBS2, their ratio, the KING-robust kinship `(shared_hets - 2 ibs0) / hets_ab`, both samples' heterozygous sites at the sites both are called, somalier's relatedness, the expected kinship and relationship from the pedigree, the observed relationship, and, for a parent and child, the opposite-homozygote rate in the family BCF.
  - `samples.tsv` has one row per member, with the pedigree sex, the inferred sex and karyotype, the X and Y copy numbers (from the cnv_coverage JSON when there is one, from somalier's depth otherwise), the X heterozygosity, the depth, the het/hom-alt ratio, the CHARR contamination estimate with the sites it used and somalier's own, the trio's Mendelian error rate, and the flags it carries.
  - `flags.tsv` has one row per flag: `level, family_id, sample_a, sample_b, flag, severity, evidence`.
- **Per cohort:** `cohorts/<cohort_name>/qc/<cohort_name>.<label>.{samples,pairs,flags,inputs}.tsv` and `.html`. These are the family tables concatenated, plus the pairs of different families related at the third degree or closer (but those of a sample with CHARR of 0.05 or more) and the cohort-level flags. The HTML is somalier's interactive plot of every pair.
- **Genotypes.** somalier relate calls a site at depth 7 or more: het when the alternate allele fraction is within 0.3-0.7, homozygous when it is below `sampleqc_hom_cutoff` or above 1 minus it, unknown otherwise. With somalier's own 0.01, a single stray read at 40x loses a homozygous site; the fewer a sample keeps, the fewer opposite homozygotes it shows and the higher its kinship with everyone, so noise and contamination read as relatedness. On a 2,672-sample cohort it made 142,000 cross-family pairs third-degree relatives; at 0.1, the unrelated pairs centre on 0, and the third degree is reached only by contaminated samples and by about 120 pairs of admixed ones.
- **Relationship classes** (KING): duplicate/MZ above 0.354; first degree 0.177-0.354, parent-offspring when IBS0/IBS2 is below 0.01 and full siblings otherwise; second degree 0.0884-0.177; third degree 0.0442-0.0884; unrelated below. A pair with fewer than 2,000 sites called in both is `low_confidence` and is never flagged.
- **Missing parents.** A parent written `0` is one not in the pedigree, not someone else: the expected kinship is the least the pair can have. A pair observed more distant than expected contradicts the pedigree; one observed closer is related beyond what it records (`RELATED_BEYOND_PEDIGREE`, `CONSANGUINITY`), unless they are first-degree relatives who both have both parents recorded. Siblings who share one recorded parent, the other missing for either, may be full or half siblings.
- **Dummy members.** A `DUMMY_` member of the pedigree (see the README's pedigree section) has no extract and no row in any table, but the expected kinships go through it: MZ twins given the same dummy parents are siblings, so `MZ_TWIN_CANDIDATE` rather than `DUPLICATE_OR_MZ`, and an aunt joined to a parent by dummy grandparents is expected at the second degree. A dummy parent is written `0` in somalier's own pedigree, and has no genotypes for the Mendelian counts.
- **Across families.** KING's estimator assumes one homogeneous population: two admixed samples of similar ancestry reach 0.08 unrelated, so a cross-family third degree is information. A pair at the second degree or closer is `high` only when its IBS0/IBS2 is at most 0.6 times the median of the cohort's clearly unrelated pairs (kinship below 0.0221): second-degree relatives share no allele IBD at about half of their sites, so they have about half the opposite homozygotes of unrelated pairs. With fewer than 100 such pairs there is no yardstick, and the check is skipped.
- **Contamination.** CHARR as Lu et al. (AJHG 2023) define it, computed from the extract: at autosomal sites of depth 20-100 that look homozygous (minor reads below 20%, other alleles at most 4%), whose absent allele has a frequency of 10-90% in the sites VCF, the minor reads over the depth divided by that frequency, averaged over at least 500 sites. It estimates the fraction of reads from another person: below 1% in clean WGS; the usual exclusion threshold is 5%. somalier's own `contamination_charr`, kept as `somalier_charr`, calls homozygous sites more loosely: it reads several times higher and mostly measures low depth. Contamination turns homozygous sites into hets and so drags CHARR down at high levels: read it with `het_hom_alt_ratio` and `p_middling_ab`.
- **Flags.**

  | Flag | Severity | Raised when |
  |---|---|---|
  | `NOT_PARENT` | high | a recorded parent and child are not parent-offspring; the evidence names which parent |
  | `HALF_SIB_NOT_FULL`, `FULL_SIB_NOT_HALF` | high | siblings are of the other kind |
  | `EXPECTED_RELATED_UNRELATED` | high, medium past the second degree | recorded relatives are unrelated |
  | `UNEXPECTED_RELATED_IN_FAMILY` | high | first-degree relatives by their genotypes, though both have both parents recorded and the pedigree does not make them so |
  | `RELATED_BEYOND_PEDIGREE` | info | members closer than the pedigree says, which its missing parents allow: e.g. siblings whose parents are not in it |
  | `CONSANGUINITY` | info | the two parents of a member are related |
  | `RELATIONSHIP_MISMATCH` | high or medium | any other disagreement of degree |
  | `UPD_CANDIDATE` | high | a recorded parent and child whose opposite homozygotes, in the extracts, lie on one chromosome (at least 20, 80% of them): the child has no allele from that parent over their span, a whole chromosome or a segment. No het at the child's normal depth is uniparental isodisomy, three quarters of its depth or less a deletion of the parent's copy; the evidence gives the span. Replaces `NOT_PARENT` |
  | `MZ_TWIN_CANDIDATE` | info | siblings of the same sex are one genome: MZ twins, or one sample twice |
  | `DUPLICATE_OR_MZ` | high | any other two samples that are the same person or MZ twins, in a family or across families |
  | `SAMPLE_UNRELATED_TO_FAMILY` | high | a member is related to none of its relatives within the second degree |
  | `SAMPLE_SWAP_CANDIDATE` | high | such a member is related to someone of another family |
  | `CRYPTIC_RELATED` | high at the second degree or closer with the IBS0 of relatives; info at the third degree, or when the IBS0 is that of unrelated pairs | members of different families are related |
  | `RELATEDNESS_UNRELIABLE` | info | a sample with CHARR of 0.05 or more: its pairs with other families, but duplicates, are left out |
  | `MENDELIAN_INCONSISTENT` | medium | the CRAMs say parent and child, but the family BCF has too many opposite homozygotes: the CRAM and the gVCF may not be the same person |
  | `SEX_MISMATCH` | high | the inferred sex is not the pedigree's |
  | `KARYOTYPE_XXY`, `KARYOTYPE_X0`, `KARYOTYPE_XYY`, `KARYOTYPE_XXX`, ... | high | the copy numbers make a karyotype other than XX or XY |
  | `SEX_AMBIGUOUS` | medium | the copy numbers are not near whole numbers: mosaicism or a coverage artefact |
  | `X_HET_INCONSISTENT` | medium | the X heterozygosity disagrees with the number of X: two X and almost no heterozygous sites (X isodisomy, LOH), or one X and many |
  | `SEX_UNKNOWN_IN_PED` | info | the pedigree has no sex; the evidence gives the inferred one |
  | `LOW_INFORMATIVE_SITES` | medium | fewer than 5,000 somalier sites are called |
  | `POSSIBLE_CONTAMINATION` | medium from 0.02, high from 0.05 | CHARR, the estimated fraction of reads from another person |

- **Sex.** The copy numbers are `2 × chrX_ratio` and `2 × chrY_ratio` of the cnv_coverage JSON of the current `cnv_cov_bin` and `cnv_cov_mapq` when the sample has one, and somalier's X and Y depth over its autosomal depth otherwise. The karyotype is cnv_coverage's rule: each copy number rounded, `ambiguous` when either is more than 0.25 from a whole number. A karyotype with a Y is male.
- **Mendelian errors** are counted from the family's `norm.bcf` when it has a parent and child: biallelic autosomal SNVs, genotypes of GQ ≥ 20 and DP ≥ 10.
- **What reruns.** The label is the cache for the extracts. The family and cohort tables are also checked against `inputs.tsv`, the snapshot of what they were built from: the family's pedigree rows (without the phenotype), which members have an extract and a coverage JSON, and whether the BCF was read. A family is redone when its snapshot changes, so a corrected parent or sex reruns it even though every file is in place. The cohort is redone when any family is, or when a family is added or removed. A member with no CRAM is left out, with a warning.

### UPD

| Key | Type | Default | Effect |
|---|---|---|---|
| `upd_label` 📄 | string | `"upd_v4"` | Label for the tables in every output file name. Letters, digits, `.`, `_` and `-` only. A new label reruns every family, and is the only way a new version of the analysis reaches families already done: their snapshot does not change with the code. `upd_v1` was the first version (#50), `upd_v2` the second (#51), `upd_v3` the third (#53), `upd_v4` the fourth (#54), which also writes the windows table (#55): no cohort had run it before. |
| `upd_af_field` | string | `"AF_genomes"` | The INFO field of `gnomad_file` giving the frequency of the missing parent's allele, in a duo, and of the homozygosity score, in solo mode. The lab's `gnomad_v4.1_allChroms.bcf` has `AF_genomes` and no plain `AF`. A field the file lacks fails the family's task. Changing it redoes every family and the cohort. |

Also read: `gnomad_file` (**required**), `ancestry_reference` (the mask), `ancestry_panel_name`
(the roh tables), `cnv_cov_bin` and `cnv_cov_mapq` (the coverage JSONs). The step is report-only.
It runs in the sample QC image.

- **Per family:** `families/{S1}/{S2}/<FID>/upd/<FID>.<label>.{segments,events,members,windows,inputs}.tsv`.
  - `segments.tsv`: one row per segment of an event, a run of one UPD state, with the member, the mode and the parents used, `chrom`, `start`, `end`, `length_mb`, `state` (`iso` or `het`), `origin` (`maternal`, `paternal` or `unknown`), `source` (`hmm` or `roh`), the score columns below, `depth_ratio`, `roh_fraction` (the share of the segment in the member's ROH, when the roh tables are read), `interpretation` (`isodisomy` or `heterodisomy`), `confidence` and `event_id`.
  - `events.tsv`: one row per event, its segments collapsed: `composition` (`iso`, `het` or `mixed`), `extent` (`whole_chromosome` or `segmental`), `terminal`, the summed score and support fraction, `interpretation` (`UPD`), `imprinting` (the disorder the parent of origin implies at each imprinted region the event spans; for a solo event, at every imprinted region of its chromosome), `confidence`, `shared_with` (the siblings with an overlapping event) and `notes`.
  - `members.tsv`: one row per member in the BCF: its mode, the parents used, whether chrX was read, the sites, the error rate (sites normal inheritance does not explain, outside the segments), its `status` and its number of events.
  - `windows.tsv`: the evidence behind the calls, for genetics-viz to draw: one row per member in the BCF with sites (whatever its `status`), chromosome and 500 kb window, fixed from position 1 (`start` 1, 500001, …, `end` the window's last base), a window with no site left out. Its sites are those of the member's model, with the parents it used (`father_used`, `mother_used`): in a trio or duo the sites its HMM reads, in solo mode its called sites under its depth cap; the autosomes, and chrX outside the PARs for a member with two X. Columns `family_id`, `sample`, `mode`, `chrom`, `start`, `end`, `n_sites`, `n_het` and `n_hom_alt` (the member's calls), `n_mat` (father homozygous and the member het: its other allele is the mother's), `n_pat` (the same with the mother homozygous), `n_lacks_father` and `n_lacks_mother` (the member and that parent opposite homozygotes), and `depth_ratio` (the member's mean DP over the sites, over its genome-wide mean). A parent the model does not use is `.` in its two columns: a duo fills only the present parent's, solo mode none. Normal inheritance gives `n_mat` about `n_pat` and almost no opposite homozygote; a maternal isodisomy no het, no `n_mat` or `n_pat`, and many `n_lacks_father`; a maternal heterodisomy no `n_pat`; a deletion of the father's copy looks like a maternal isodisomy at a `depth_ratio` of about 0.5. On chrX the father is hemizygous, so every het of a daughter counts in `n_mat`.
- **Per cohort:** `cohorts/<cohort_name>/upd/<cohort_name>.<label>.{segments,events,members,inputs}.tsv`: the family tables concatenated; the windows stay per family. The events table adds `recurrent_families` and `recurrent`: an event overlapped, by at least half of both, by events of three other families or more is `recurrent`, most likely an artefact region.
- **Modes.** A parent counts only when it is in the family BCF: a parent written `0`, a `DUMMY_` parent or one without a gVCF is missing.
  - **trio**, both parents in the BCF: a hidden Markov model of the child's genotype given its parents', with the states N (one allele from each parent), `iso_mat`, `iso_pat` (two copies of one homolog of that parent), `het_mat`, `het_pat` (both homologs of that parent), and three of normal inheritance with a haplotype the parents share: `auto` (the child got the same haplotype from both: autozygosity), `ibd_mat` and `ibd_pat` (the haplotype the child got from one parent is the other parent's untransmitted one, so the child's genotype is that parent's). An isodisomy or heterodisomy is therefore called on the Mendelian errors it explains, not on a ROH or a run of genotypes identical to a parent's, which short haplotypes shared by unrelated parents give all over the genome.
  - **duo**, one parent in the BCF: the five UPDhmm states, the missing parent transmitting the ALT at its gnomAD frequency. A UPD of the missing parent rests on opposite homozygotes with the present one. A UPD of the present parent looks like autozygosity, related parents or a shared haplotype, so it is held to the extents of solo mode: `high` over a whole chromosome, `medium` for a terminal one of 5 Mb or more, `low` for an interstitial one of 15 Mb or more (10 Mb on a chromosome with an imprinted region), lengths of covered genome, and not written below them; it is `low` unless the member is isolated (below).
  - **solo**, no parent in the BCF, founders included: runs of the `roh` step's segments of 1.5 Mb or more, at most 5 Mb of covered genome apart. A run at 0.75 of the depth or less is a deletion, a lost copy reading as a ROH (a 15q11-q13 deletion as one of 4.8 Mb), and is set aside. The other runs are isodisomies of unknown origin, of a member whose rest of the genome is outbred (isolated): `whole_chromosome` when a run spans 90% of the chromosome (`high`); two runs or more of 5 Mb on one chromosome, made one event (`medium`): the isodisomic parts of a mixed UPD, which leaves heterodisomy between them; `terminal` from 5 Mb (`medium`); `interstitial` from 15 Mb (`low`), or from 10 Mb and `medium` on a chromosome with an imprinted region. Lengths, the run's and the chromosome's, are of covered genome: a stretch the mask left without sites, a centromere say, does not count. `medium` or `high` also needs a support fraction of 0.9: a real ROH holds almost no het call. On those chromosomes any of them marks a possible UPD of the whole chromosome, whose heterodisomic parts a ROH cannot show (a maternal UPD15 is mostly heterodisomy), so the event carries every imprinted region of its chromosome: an interstitial ROH on distal 15q reads `Prader-Willi if maternal or Angelman if paternal (15q11-q13, if the UPD reaches it)`. The roh step's `upd_candidate_chrom`, the chromosome holding most of the long ROH, is noted, not taken for a whole-chromosome isodisomy: one long interstitial ROH makes it. Heterodisomy cannot be seen without a parent. Without roh tables the member's status is `no_roh_tables`.
- **Sites.** Biallelic SNVs, genotypes of GQ ≥ 20 and DP ≥ 10 and at most 2.5 times the member's median DP; a position held by several records (a split multiallelic site) is dropped. When `ancestry_reference` has its roh/ tier, a site outside the roh panel's covered segments (runs of panel sites at most 300 kb apart) is dropped: the panel leaves out the centromeres, heterochromatin, the MHC, the acrocentric arms and segmental duplications, where mismapped reads give clusters of Mendelian errors that read as short UPDs. Without it the run warns. chrX, which the panel does not hold, is not masked. Autosomes, and chrX outside the PARs for a member with two X: karyotype XX from its coverage JSON, or with none, a female pedigree sex and an X heterozygosity of 0.1 or more. The father is hemizygous on X.
- **Score.** For each site, LR is the likelihood of the segment's state over normal inheritance. `lod` is the sum of log10 LR. `n_supporting` and `n_contradicting` count the sites favouring the state or normal inheritance at least 2-fold, `n_strong` those favouring the state at least 10-fold, which normal inheritance explains only as an error: in a trio, the Mendelian errors the state explains. `support_fraction` is reported, not used.
- **Segments** are runs of one UPD state, cut wherever two of their sites are more than 300 kb apart: a segment never spans a region without sites, such as one the mask removed, so its length is covered genome.
- **Confidence.** `high`: `lod` ≥ 20, at least 2 Mb, a depth ratio of 0.8-1.25, and `n_strong` ≥ 20 (`n_supporting` ≥ 50 for a UPD of the present parent of a duo, whose sites are each weak). `medium`: `lod` ≥ 10, at least 0.5 Mb, and `n_strong` ≥ 10 (or `n_supporting` ≥ 20). Both also need a `support_fraction` of 0.9: a trio state, a duo's present-parent state, or an isodisomy of its missing parent, predicts the child's genotypes closely, so a real UPD is all but never contradicted (1.00 on the real events seen) where mapping artefacts are (0.11-0.85). A heterodisomy of a duo's missing parent, whose alleles come from the population, needs 0.5. `low` otherwise. A depth ratio of 0.75 or less is a deletion, one copy looking like an isodisomy of the other parent: it is set aside, never an event. This step calls UPD; CNVs, the 15q11-q13 deletions of Prader-Willi and Angelman among them, are left to WisecondorX.
- **Events.** Segments of one member, chromosome and parent of origin, each `medium` or `high` and at most 5 Mb of covered genome apart (a stretch without sites does not count, so the arms of a chromosome join across its centromere), make an event, whose confidence is its best segment's (capped for a duo as above). A `low` segment makes no event and is not written.
- **Siblings.** Two children sharing a parent, with events on one chromosome overlapping by half of the shorter, share an artefact (a parental CNV, or reads mismapped alike), not a UPD, which is sporadic: both events are `low`, with `shared_with` and a note. MZ twins (genotypes discordant at 2% or less of the sites both are called at) keep the UPDs they share.
- **Imprinting** is by region: 6q24 (PLAGL1), 7p12 (GRB10), 7q21 (PEG10), 7q32 (MEST), 11p15.5 (H19/IGF2 to CDKN1C), 14q32 (DLK1-MEG3), 15q11-q13 (MKRN3 to UBE3A) and 20q13 (GNAS), GRCh38 windows around their genes. An event that spans one gets the disorder of its parent of origin there, e.g. `Prader-Willi (15q11-q13)`; a segmental event elsewhere on the chromosome gets none.
- **Isolated** means F_ROH, at segments of 1.5 Mb or more on the other autosomes, below 1%, and parents not related beyond the 4th degree (kinship ≤ 0.0221) when the roh step measured them: the roh step's own rule for `upd_candidate_chrom`.
- **Parentage.** A child and a parent with opposite homozygotes at more than 1% of the sites outside the child's segments are not parent and child (`sampleqc`'s `NOT_PARENT`): that parent is dropped and the child redone, with the reason in its `status`.
- **Gaps.** Mosaic UPD, which only shows as skewed allele fractions, is not seen. Small events and repeat-rich regions give false calls, which `recurrent`, the sibling rule and the thresholds are for; the immunoglobulin loci (IGH, IGK, IGL), where blood DNA rearranges somatically, are among them. Exome data are sparse for heterodisomies under 10 Mb.
- **What reruns.** The label is the cache. The family and cohort tables are also checked against `inputs.tsv`: the family's pedigree rows, which members are in the BCF and have a coverage JSON, which roh tables are read, and the gnomAD file and field. A family is also redone when its BCF, its roh tables or a member's coverage JSON is made in the run. The cohort is redone when any family is, or when a family is added or removed.

### Containers and resources

| Key | Type | Default | Effect |
|---|---|---|---|
| `apptainer_cache` | dir | `""` | Apptainer image cache. |
| `singularity_cache` | dir | `""` | Singularity image cache. |
| `max_memory` | memory string | `'460.GB'` | Per-process memory ceiling. |
| `max_cpus` | integer | `95` | Per-process CPU ceiling. |
| `max_time` | duration string | `'240.h'` | Per-process wall-clock ceiling. |

A task asking for more than a ceiling — by default, through `resources`, or after a retry has
added memory — is clamped down to it, not failed: the ceilings are every process's
`resourceLimits`. A `resources` value above a ceiling is also warned about when the run starts,
and a ceiling that does not parse stops it there.

### `resources`: per-cohort cpus, memory and time

The cpus, memory and time of every process have defaults in `nextflow.config`, sized for a WGS
family. A cohort that needs something else — an exome cohort needs far less, a 40-member
family far more — says so here rather than in a config file of its own, so that it is part of
the file every run of the cohort is recorded against.

```yaml
resources:
  PYWOMBAT:                    # a process, by name
    memory: 4.GB               # memory for the first attempt
    memory_step: 32.GB         # added for each retry
    time: 4.h
    families:                  # by fid: an exact id, or a glob with *
      FAM042: { memory: 200.GB }
  GLNEXUS_FAMILY:
    families:
      FAM042: { memory: 300.GB, cpus: 8 }
  MAKE_EXAMPLES:
    samples:                   # by barcode
      'IP*': { memory: 460.GB }
  'VEP_ANNOTATION|VEP_MNV_ANNOTATION':   # a regex, like withName
    cpus: 4
```

**Settings.** `cpus` (a whole number), `memory` and `memory_step` (with a unit: `8.GB`,
`500.MB`), `time` (with a unit: `4.h`, `30.m`). Any of them can be given for the process, for
a glob, or for an id. `families:` applies to the processes that run once per family, keyed by
`fid`, and `samples:` to those that run once per sample, keyed by `barcode`. A process that
runs once for the whole cohort (`SNVS_COHORT_MERGE`, `MERGE_WOMBAT`...) only takes the
process-level settings.

**What wins.** For each setting, the most specific value that is given: an exact id, then the
first glob that matches it (in the order of the file), then the process, then the default in
`nextflow.config`. A level that leaves a setting out keeps the one from the level below, so
`FAM042: { memory: 200.GB }` keeps the process's `cpus`, `time` and `memory_step`. Ids are
compared as text, so a numeric fid matches whether or not it is quoted.

**Retries.** Memory on attempt *n* is `memory + (n - 1) × memory_step`. Every default keeps
its retry step — `{ 8.GB * task.attempt }` became `memory 8.GB, memory_step 8.GB` — so a
cohort that overrides only `memory` still gets more on each retry. `cpus` and `time` do not
grow. A task killed by its time limit therefore fails again on every retry, so leave a
generous margin on `time`.

**Process names.** A process key is a name as the process runs, matched in full, or a regex
matched against the whole name. A process included under another name runs under that name:
`BAZAM_BWA_MEM2_REALIGN` runs only as `BAZAM_BWA_MEM2_REALIGN_37` and `_38`, so those (or
`'BAZAM_BWA_MEM2_REALIGN_.*'`) are what an entry has to say. This is the one place where an
entry differs from a `withName` block in a config file, which matches both.

**Checked before anything runs.** An unknown setting, a process name that matches nothing, a
setting under the wrong parent, or a value that does not parse stops the run at startup with
the full list. Nothing else in the parameters file is checked this way.

**Safe to change mid-cohort.** cpus, memory and time are not part of a task's cache key: a
`-resume` after changing `resources` keeps every task that already succeeded, and the new
values apply to the rest.

**DeepVariant.** `MAKE_EXAMPLES` shards its work over `deepvariant_threads`, not over its
cpus. Change one and change the other; see
[DeepVariant and GLnexus](#deepvariant-and-glnexus).

**Sizing it from the cohort's own runs.** Every run leaves a trace in
`runs/<COHORT>/reports/`. `scripts/resources_report` reads them all and shows, per process,
what was asked for against what was used, the tasks killed by a limit, and the families or
samples that used the most:

```bash
scripts/resources_report /pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/runs/CANDY_mpx --suggest
```

`--suggest` prints a `resources:` block to review: each process sized to cover 99% of its
tasks with a 25% margin, and the families and samples above that listed by id. The block is
printed, never written anywhere.

**Sizing from the data.** A default in `nextflow.config`, a process's `ext`, can follow the
size of the task's inputs: `ext = [memory: { 8.GB + 2.GB * barcodes.size() }, memory_step: 8.GB]`
grows with the size of the family, and `memory: { nextflow.util.MemoryUnit.of(4 * annotated_parquet.size()) }`
with the size of a file. The settings in `resources` still override it.

**Beyond cpus, memory and time.** For anything else — a container, an `errorStrategy` — give
a config file with `--config FILE`. It is loaded on top of the pipeline's own configuration.
A `cpus`, `memory` or `time` set in such a file replaces the pipeline's, and `resources` then
no longer applies to that process — or to any process, if it is set for the whole `process`
scope rather than in a `withName` block. To change a default and keep `resources`, set the
process's `ext` instead: `withName: 'PYWOMBAT' { ext = [memory: 16.GB] }` changes its memory and
keeps the pipeline's other defaults for it.

## Gotchas

### Keys with no default become the string `null`

Eighteen keys have no default anywhere:

`cohort_name`, `pedigree`, `work_dir`, `gnomad_file`, `gnomad_filter_field`,
`gnomad_filter_threshold`, `vep_config`, `vep_config_name`, `annotation_annotation_path`,
`annotation_annotation_list`, `annotation_gencode`, `wisecondorx_reference`,
`wisecondorx_binsize`, `wisecondorx_predict_args`, `wombat_config_path`, `wombat_config_list`,
`extractor_tsvs_list`, `liftover_chain`.

`pedigree` and `work_dir` have sensible fallbacks in the code, and `cohort_name` is checked at
startup. The rest have neither: several are interpolated straight into paths, so omitting one
produces a file with `null` in its name rather than an error.

### Filename keys are the cache

The pipeline decides what to recompute by looking for files on disk. Any key marked 📄 above is
part of a filename, so changing it makes the pipeline stop finding the old results and produce
new ones alongside them. Changing anything *not* in a filename does the opposite: the old
results stay, and are reused, even though they were produced with different settings.

This matters most for **`ancestry_panel_name`**, which is the only thing that invalidates
ancestry results. Raising `ancestry_min_dp` or `ancestry_min_gq`, or swapping the reference
bundle, without bumping the label leaves every existing extraction in place — so families
processed before the change and families processed after it are scored differently, in the same
cohort, with nothing to indicate it. Encode the bundle version and the thresholds in the label:

```yaml
ancestry_panel_name: "apgs_b1.0.0_dp10gq20"
```

**`expansionhunter_catalog_name`** works the same way for the ExpansionHunter calls. Pointing
`expansionhunter_catalog` at another catalog without a new label reuses every existing call,
made with the old catalog, and mixes in new ones made with the new catalog. **`melt_label`**
does the same for the MELT calls. **`cnv_cov_bin`** and **`cnv_cov_mapq`** make the tag
`cov_b<bin>q<mapq>` of every CNV coverage file, and **`cnv_panel_version`** the directory of a
panel, which is built once and never rebuilt under the same version. **`sampleqc_label`** does
the same for the somalier extracts, and **`upd_label`** for the UPD tables. The sample QC and
UPD tables are the exception to the rule: they are also rebuilt when the pedigree or the
settings they were built from change, see [Sample QC](#sample-qc) and [UPD](#upd).

### Three different directories

| Key | Holds |
|---|---|
| `data` | The published cohort tree: samples, families, cohorts. The results. |
| `scratch` | Intermediate CRAMs during alignment. |
| `work_dir` | Nextflow's own work directory, needed for `-resume`. |

## Keys that look live but are inert

These **have no effect whatsoever** in the current code. Most were removed from the repository's
own parameters files, but they carried confident comments there, and cohort files copied from
those templates still have them. They are listed here so you can recognise them rather than tune
them.

| Key(s) | Why it does nothing |
|---|---|
| `annotation_dnm_min_callrate`, `annotation_dnm_min_DP`, `annotation_dnm_min_GQ`, `annotation_dnm_min_VAF` | Were read only by the `DNM_EXTRACTION` process, which no workflow ever invoked and which has been removed. The de novo thresholds in effect are whatever a wombat config specifies. |
| `ref_par1_start`, `ref_par1_end`, `ref_par2_start`, `ref_par2_end` | Same: the pseudo-autosomal coordinates were only consumed by the removed `DNM_EXTRACTION`. |
| `slurm_account`, `slurm_partition` | SLURM options are hardcoded in the `clusterOptions` strings in `nextflow.config`. Change them there. |
| `bin_size` | Was in two of the files in `params_example/`. The code reads `bin`, so a cohort file that still sets `bin_size` silently runs with `bin = 1000`. |
| `oldref`, `cram4realignment_pattern` | Superseded by `old_ref_37` / `old_ref_38`. `cram4realignment_pattern` is still defined in `nextflow.config`. |
| `enable_conda`, `singularity_pull_docker_container`, `skip_alignment`, `skip_deepvariant` | Defined in `nextflow.config` and never read. |

## Examples

### Minimal: align a cohort and stop

```yaml
steps: ["alignment"]

data: "/pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/"
scratch: "/pasteur/appa/scratch/ghfc"
# The pedigree is read from
# <data>/cohorts/MY_COHORT/MY_COHORT.pedigree.tsv
cohort_name: "MY_COHORT"

ref: "/pasteur/helix/projects/ghfc_wgs/references/GRCh38/chr/GRCh38.fasta"
ref_name: "GRCh38_GIABv3"

fastq_pattern: "*_R{1,2}.fastq.gz"
```

Everything else falls back to its default. No variant calling happens, and none of the
annotation keys are read.

### Typical: a full SNV cohort

This is the shape of the two files in [`params_example/`](../params_example/).

```yaml
steps: ["alignment", "deepvariant_sample", "deepvariant_family", "annotation", "wombat", "wombat_cohort", "snvs_cohort", "wisecondorx"]

data: "/pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/"
scratch: "/pasteur/appa/scratch/ghfc"
work_dir: "/pasteur/appa/scratch/ghfc/work-ghfc-grch38"
cohort_name: "MY_COHORT"

ref: "/pasteur/helix/projects/ghfc_wgs/references/GRCh38/chr/GRCh38.fasta"
ref_name: "GRCh38_GIABv3"

# Realignment sources, in place of fastq_pattern
old_ref_38: "/pasteur/helix/projects/ghfc_wgs/references/GRCh38/chr/hs38DH.fa"
old_cram_38: "/pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/old_cram_38/"
old_ref_37: "/pasteur/helix/projects/ghfc_wgs/references/GRCh37/chr/hs37d5.fa"
old_cram_37: "/pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/old_cram_37/"
fastq_pattern: ""

deepvariant_threads: 96
glnexus_config: "DeepVariant_unfiltered"

gnomad_file: "/pasteur/appa/scratch/ghfc/scratch-cache/gnomad/gnomad_v4.1_allChroms.bcf"
gnomad_filter_field: "AF_genomes"
gnomad_filter_threshold: "0.01"

vep_config: "/pasteur/appa/scratch/ghfc/scratch-cache/GRCh38_112_standard.ini"
vep_config_name: "VEP_GRCh38_112_standard"

annotation_annotation_path: "/pasteur/appa/scratch/ghfc/scratch-cache/annotations"
annotation_annotation_list: ["AlphaMissense_hg38.bcf", "revel.grch38.bcf", "LCR.grch38.bed.gz"]
annotation_gencode: "gencode.v47.basic"

wombat_config_path: "/pasteur/appa/scratch/ghfc/scratch-cache/pywombat"
wombat_config_list: ["rare_high_impact.yml", "de_novo_mutations.yml"]

wisecondorx_reference: "/pasteur/helix/projects/ghfc_wgs/references/GRCh38/WisecondorX/reference_wisecondorX_hg38_{binsize}_refsize100.npz"
wisecondorx_binsize: 1000
wisecondorx_predict_args: "--zscore 5 --minrefbins 50 --resegment-all --fix-zero-bins"

apptainer_cache: "/pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/container-cache/"
max_memory: "460.GB"
max_cpus: 95
max_time: "240.h"
```

### Ancestry on an existing cohort

The ancestry step reads gVCFs only, so it can be run on its own, even where annotation is
incomplete. It needs `data`, `cohort_name` and its own three required keys:

```yaml
steps: ["ancestry"]

data: "/pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/"
cohort_name: "MY_COHORT"
ref: "/pasteur/helix/projects/ghfc_wgs/references/GRCh38/chr/GRCh38.fasta"
ref_name: "GRCh38_GIABv3"

ancestry_reference: "/pasteur/helix/projects/ghfc_wgs/references/ancestry-pgs/b1.0.0"
ancestry_catalog: "/pasteur/helix/projects/ghfc_wgs/references/ancestry-pgs/catalog.tsv"
ancestry_panel_name: "apgs_b1.0.0_dp10gq20"
ancestry_min_dp: 10
ancestry_min_gq: 20
ancestry_min_coverage: ""
ancestry_model: ""
```

Add `deepvariant_sample` to `steps` if some individuals have no gVCF yet; otherwise the run
stops and names them.

## Where the results land

Output paths are built from `data`, `cohort_name`, and the filename keys marked 📄 above, under
a two-level shard derived from each sample or family ID. The full layout, and the exact list of
files the pipeline looks for when deciding what to skip, are in the README:
[Input Data Structure](../README.md#input-data-structure) and
[Smart File Detection](../README.md#smart-file-detection).

One correction to that list: the WisecondorX NPZ is
`<barcode>.<wisecondorx_binsize>.npz`, not `<barcode>.npz` — the bin size is part of the name,
which is why changing it recomputes every sample.

## See also

- [HOWTO.md](../HOWTO.md) — running the pipeline on the Institut Pasteur cluster
- [COHORT_STATE.md](../COHORT_STATE.md) — what each run records about the parameters it used,
  including the checksum of the resolved parameter map
- [README.md](../README.md) — pipeline stages, input layout and outputs
- [`params_example/`](../params_example/) — real cohort files, subject to the caveats in
  [Keys that look live but are inert](#keys-that-look-live-but-are-inert)
