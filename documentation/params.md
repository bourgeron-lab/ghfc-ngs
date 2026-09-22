# The parameters file (`params.yml`)

A cohort's parameters file is the only thing an operator writes to run this pipeline. It is
a flat YAML mapping — **no nesting anywhere**, every key is a scalar or a list of scalars —
and it is read by `nextflow run` through `-params-file`.

This document is meant to be usable **next to a parameters file that has no comments**: it
lists every key the code actually reads, its default, when it is required, and what it
changes. It also lists the keys that appear in existing parameters files and do nothing, so
you can recognise them for what they are.

> **The pipeline does not validate key names or types.** There is no JSON schema. An unknown
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

49 keys are read by the code; 17 more are in circulation and do nothing. Everything either
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
| `ancestry_panel_name` 📄 | `ancestry` — **required** | [Ancestry and PGS](#ancestry-and-polygenic-scores) |
| `ancestry_reference` | `ancestry` — **required** | [Ancestry and PGS](#ancestry-and-polygenic-scores) |
| `annotation_annotation_list` | `annotation` | [Annotation](#annotation) |
| `annotation_annotation_path` | `annotation`, `wisecondorx` | [Annotation](#annotation) |
| `annotation_dnm_min_callrate` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `annotation_dnm_min_DP` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `annotation_dnm_min_GQ` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `annotation_dnm_min_VAF` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `annotation_gencode` | `wisecondorx` | [Annotation](#annotation) |
| `apptainer_cache` | always, at config-parse time | [Containers and resources](#containers-and-resources) |
| `bazam` | `alignment`, realignment only | [Alignment tools and coverage](#alignment-tools-and-coverage) |
| `bin` 📄 | `alignment` | [Alignment tools and coverage](#alignment-tools-and-coverage) |
| `bin_size` ✗ | never — you want `bin` | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `bwa_mem2` | `alignment` | [Alignment tools and coverage](#alignment-tools-and-coverage) |
| `cohort_name` 📄 | always — **required** | [Core paths and identity](#core-paths-and-identity) |
| `cram4realignment_pattern` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `data` | always — **required** | [Core paths and identity](#core-paths-and-identity) |
| `deepvariant_threads` | `deepvariant_sample` | [DeepVariant and GLnexus](#deepvariant-and-glnexus) |
| `enable_conda` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `extractor_tsvs_list` | `extractor` | [Extractor](#extractor) |
| `fastq_pattern` | `alignment` | [Input discovery](#input-discovery) |
| `glnexus_config` | `deepvariant_family` | [DeepVariant and GLnexus](#deepvariant-and-glnexus) |
| `gnomad_file` | `annotation` | [Annotation](#annotation) |
| `gnomad_filter_field` | `annotation` | [Annotation](#annotation) |
| `gnomad_filter_threshold` | `annotation` | [Annotation](#annotation) |
| `liftover_chain` | `extractor`, GRCh37 inputs | [Extractor](#extractor) |
| `max_cpus` | always, at config-parse time | [Containers and resources](#containers-and-resources) |
| `max_memory` | always, at config-parse time | [Containers and resources](#containers-and-resources) |
| `max_time` | always, at config-parse time | [Containers and resources](#containers-and-resources) |
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
| `sambamba` | `alignment` | [Alignment tools and coverage](#alignment-tools-and-coverage) |
| `samblaster` | `alignment` | [Alignment tools and coverage](#alignment-tools-and-coverage) |
| `samtools` | `alignment` | [Alignment tools and coverage](#alignment-tools-and-coverage) |
| `scratch` | `alignment` | [Core paths and identity](#core-paths-and-identity) |
| `singularity_cache` | always, at config-parse time | [Containers and resources](#containers-and-resources) |
| `singularity_pull_docker_container` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `skip_alignment` ✗ | never — use `steps` | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `skip_deepvariant` ✗ | never — use `steps` | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `slurm_account` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `slurm_partition` ✗ | never | [**Inert — does nothing**](#keys-that-look-live-but-are-inert) |
| `steps` | always — **required** | [`steps`](#steps-the-key-that-decides-everything-else) |
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
| 1. Built-in defaults | the `params { }` block in [`nextflow.config`](../nextflow.config) | Covers 31 of the 49 keys. The other 18 are `null` unless you set them. |
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
steps: ["alignment", "deepvariant_sample", "deepvariant_family", "annotation", "wombat", "snvs_cohort"]
```

Required, and must be non-empty. Every entry must be one of these nine:

| Step | Produces |
|---|---|
| `alignment` | CRAMs from FASTQ, or realigned from existing GRCh37/GRCh38 CRAMs, plus coverage bedgraphs |
| `deepvariant_sample` | Per-individual gVCF and VCF, plus VAF bedgraphs |
| `deepvariant_family` | Joint family calls via GLnexus, normalised, plus per-family pedigrees |
| `annotation` | gnomAD frequencies, the rare/common split, VEP, and extra bcftools annotations |
| `wombat` | PyWombat filtering of the annotated BCF, one TSV per config |
| `snvs_cohort` | Cohort-level merge of common variants and of the wombat results |
| `wisecondorx` | CNV/SV calling, per sample then merged per family and per cohort |
| `extractor` | Specific variants pulled out of family BCFs, wombat TSVs or gVCFs |
| `ancestry` | Ancestry components, admixture and polygenic scores, per family and per cohort |

### The rule that catches people out

The pipeline looks at what already exists on disk and only does the missing work. But it will
**not silently run a step you did not ask for**: if work is needed for a step that is not in
`steps`, the run stops.

There is one refinement worth knowing, because it is what makes partial runs possible: a
missing prerequisite is only an error when something in *this* run would consume it.

- `alignment` and `deepvariant_sample` are checked unconditionally — if any individual needs
  them and they are not listed, the run fails.
- `deepvariant_family` is only demanded when one of `deepvariant_family`, `annotation`,
  `wombat`, `snvs_cohort` or `extractor` is listed.
- `annotation` is only demanded when one of `annotation`, `wombat` or `snvs_cohort` is listed.

This is why `steps: ["ancestry"]` runs happily on a cohort whose annotation is incomplete:
ancestry reads gVCFs directly and never touches the normalised or annotated call sets.

### Failures you will see

| Message | Cause |
|---|---|
| `ERROR: --data parameter is required` | `data` unset or empty |
| `ERROR: cohort_name parameter is required` | `cohort_name` unset or empty |
| `ERROR: --steps parameter is required` | `steps` missing or `[]` |
| `ERROR: Invalid steps specified: ...` | a value outside the nine above — **or** `steps` given as a string (`steps: "alignment"`) instead of a list |
| `ERROR: the 'ancestry' step requires ...` | `ancestry` listed without its three required keys |
| `ERROR: Pedigree file not found: ...` | `pedigree`, or the cohort's default pedigree path, does not exist |
| `... is required for N individuals but not included in steps parameter` | the rule above |

Two of those messages name `--data` and `--steps` as though they were command-line flags. They
are not: the keys they refer to are `data:` and `steps:` in the parameters file.

## Key reference

49 keys are read by the code. In the tables below, **no default** means the key is `null` when
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

`deepvariant_threads` is **independent of `max_cpus` and of the SLURM allocation**. It exists
because DeepVariant's throughput is hardware-sensitive, and it overrides the process CPU
setting without changing what SLURM actually reserves. Setting it above the cores you were
allocated oversubscribes the node. Keep it, `max_cpus`, and the `withName` blocks in
`nextflow.config` consistent by hand.

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

Despite its prefix, `annotation_gencode` belongs to the `wisecondorx` step, not to `annotation`.
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
`ancestry_panel_name` aborts the run before any work starts.

### Containers and resources

| Key | Type | Default | Effect |
|---|---|---|---|
| `apptainer_cache` | dir | `""` | Apptainer image cache. |
| `singularity_cache` | dir | `""` | Singularity image cache. |
| `max_memory` | memory string | `'460.GB'` | Per-process memory ceiling. |
| `max_cpus` | integer | `95` | Per-process CPU ceiling. Also caps GLnexus at `min(4, max_cpus)`. |
| `max_time` | duration string | `'240.h'` | Per-process wall-clock ceiling. |

A process asking for more than a ceiling is clamped down to it with a warning, not failed.

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

### Three different directories

| Key | Holds |
|---|---|
| `data` | The published cohort tree: samples, families, cohorts. The results. |
| `scratch` | Intermediate CRAMs during alignment. |
| `work_dir` | Nextflow's own work directory, needed for `-resume`. |

## Keys that look live but are inert

These appear in the repository's own parameters files, carry confident comments, and **have no
effect whatsoever** in the current code. They are listed here so you can recognise them rather
than tune them.

| Key(s) | Why it does nothing |
|---|---|
| `annotation_dnm_min_callrate`, `annotation_dnm_min_DP`, `annotation_dnm_min_GQ`, `annotation_dnm_min_VAF` | Read only by the `DNM_EXTRACTION` process, which is **defined but never invoked** by any workflow. The de novo thresholds in effect are whatever a wombat config specifies. |
| `ref_par1_start`, `ref_par1_end`, `ref_par2_start`, `ref_par2_end` | Same: the pseudo-autosomal coordinates are only consumed by `DNM_EXTRACTION`. |
| `slurm_account`, `slurm_partition` | SLURM options are hardcoded in the `clusterOptions` strings in `nextflow.config`. Change them there. |
| `bin_size` | Present in both files in `params_example/`. The code reads `bin`. Those cohorts silently run with `bin = 1000`. |
| `oldref`, `cram4realignment_pattern` | Superseded by `old_ref_37` / `old_ref_38`. |
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
steps: ["alignment", "deepvariant_sample", "deepvariant_family", "annotation", "wombat", "snvs_cohort", "wisecondorx"]

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
