# GHFC WGS Family-based Variant Calling Pipeline (Nextflow)

This is a Nextflow implementation of the GHFC WGS family-based variant calling pipeline. The pipeline supports alignment from FASTQ files, realignment from existing CRAM files, individual variant calling with DeepVariant, and family-based joint calling with GLnexus.

## Features

- **Family-aware pipeline** processing based on pedigree structure
- **Smart dependency resolution** - automatically determines what needs to be run
- **BWA-MEM2 alignment** from paired-end FASTQ files
- **Realignment** from existing CRAM files using Bazam
- **DeepVariant variant calling** with 3-stage processing pipeline (individual and family-based)
- **GLnexus family calling** for joint variant calling per family
- **gnomAD frequency annotation** and rare/common variant filtering
- **VEP annotation** with custom configuration support
- **WisecondorX CNV/SV calling** with family and cohort-level analysis
- **PyWombat variant filtering** and prioritization with custom configurations
- **Cohort-level merging** for common variants and Wombat results
- **Variant extraction** from custom TSV lists with liftover support
- **Genetic ancestry and polygenic scores** per family and per cohort, projected onto a global reference panel
- **Unit merging** for samples with multiple sequencing runs
- **SLURM integration** with configurable resource allocation
- **Container support** via Apptainer/Singularity
- **Automatic file detection** - skips steps when outputs already exist
- **Comprehensive error checking** and validation

## Requirements

- Nextflow (≥23.04.0)
- SLURM workload manager (for cluster execution)
- Apptainer or Singularity (recommended for containerized execution)
- A pedigree file in TSV format

### Tool Dependencies (when not using containers)

- BWA-MEM2
- SAMtools
- Sambamba
- SAMblaster
- Bazam (Java jar file)
- Java (for Bazam)
- DeepVariant
- GLnexus
- BCFtools
- tabix (for VCF indexing)

## Pipeline Overview

This pipeline provides a comprehensive analysis framework for whole-genome sequencing data:

### Typical Usage Scenarios

**Complete SNV/INDEL Analysis:**
```bash
steps: ["alignment", "deepvariant_sample", "deepvariant_family", "annotation", "wombat", "snvs_cohort"]
```

**CNV/SV Analysis Only:**
```bash
steps: ["alignment", "wisecondorx"]
```

**Full Analysis (SNVs + SVs):**
```bash
steps: ["alignment", "deepvariant_sample", "deepvariant_family", "annotation", "wombat", "snvs_cohort", "wisecondorx"]
```

**Variant Extraction from Lists:**
```bash
steps: ["extractor"]  # Requires existing data and extractor_tsvs_list
```

**Incremental Analysis:**
- Start with `alignment` only, then add `deepvariant_sample` when ready
- Pipeline automatically detects existing files and resumes from where it left off
- Can run different workflow branches (SNVs vs SVs) independently

## Quick Start

### 1. Clone/Download the Pipeline

```bash
cd /path/to/your/workspace
# Pipeline files should be in the current directory
```

### 2. Migrating from Legacy Formats (Optional)

If you have existing data from an older version of the pipeline with VCF.gz files, you'll need to run the migration workflow once before running the main pipeline:

```bash
nextflow run migrate.nf -params-file params.yml
```

See [MIGRATION.md](MIGRATION.md) for detailed migration instructions.

### 3. Prepare Your Pedigree File

Create a TSV file with 6 columns (FID, barcode, father, mother, sex, phenotype), and save it in
the cohort's own directory, where the pipeline looks for it by default:

```
${data}/cohorts/<COHORT_NAME>/<COHORT_NAME>.pedigree.tsv
```

A pedigree kept anywhere else has to be named with the `pedigree` parameter.

```tsv
FID barcode father mother sex phenotype
FAM001 C000F2W 0 0 1 2
FAM001 C000F2X 0 0 2 1
FAM001 C000F2Y C000F2W C000F2X 1 2
FAM002 C000F3A 0 0 1 1
FAM002 C000F3B 0 0 2 1
FAM002 C000F3C C000F3A C000F3B 2 2
```

- **FID**: Family ID
- **barcode**: Individual sample ID
- **father**: Father's barcode (0 if founder)
- **mother**: Mother's barcode (0 if founder)
- **sex**: 1=male, 2=female
- **phenotype**: 1=unaffected, 2=affected

### 3. Configure Parameters

Copy and modify the example parameters file:

```bash
cp params.yml my_params.yml
# Edit my_params.yml with your specific paths and settings
```

### 4. Run the Pipeline

#### Cohort Shorthand

A cohort name given as the **first** argument stands in for its parameters file:

```bash
# Equivalent to --params-file cohorts/CANDY_mpx/CANDY_mpx.params.yml
./run_pipeline.sh CANDY_mpx

# Composes with any other option
./run_pipeline.sh CANDY_mpx --resume
```

The file is looked up under `./cohorts/` in the current directory first, then under
`$GHFC_NGS_COHORTS` (default `/pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/cohorts`, the same
per-cohort directory described in [Input Data Structure](#input-data-structure)). If neither holds
`<NAME>.params.yml`, the runner reports both paths it tried and exits without launching Nextflow.

`--params-file` remains available for parameters files outside that layout, and takes precedence
over a cohort name on the same command line.

#### Full Pipeline (Default)

```bash
# With Apptainer (runs all steps)
./run_pipeline.sh --profile slurm,apptainer --params-file my_params.yml

# With Singularity
./run_pipeline.sh --profile slurm,singularity --params-file my_params.yml
```

#### Step-by-Step Execution

```bash
# Run only alignment step
./run_pipeline.sh --profile slurm,apptainer --steps "alignment" --params-file my_params.yml

# Run alignment and individual variant calling
./run_pipeline.sh --profile slurm,apptainer --steps "alignment,deepvariant_sample" --params-file my_params.yml

# Run only individual variant calling (assumes CRAM files exist)
./run_pipeline.sh --profile slurm,apptainer --steps "deepvariant_sample" --params-file my_params.yml

# Run only family calling (assumes gVCF files exist)
./run_pipeline.sh --profile slurm,apptainer --steps "deepvariant_family" --params-file my_params.yml

# Run full SNV/INDEL pipeline with annotation
./run_pipeline.sh --profile slurm,apptainer --steps "alignment,deepvariant_sample,deepvariant_family,annotation" --params-file my_params.yml

# Run CNV/SV calling with WisecondorX
./run_pipeline.sh --profile slurm,apptainer --steps "alignment,wisecondorx" --params-file my_params.yml

# Full pipeline (alignment through annotation and Wombat)
./run_pipeline.sh --profile slurm,apptainer --steps "alignment,deepvariant_sample,deepvariant_family,annotation,wombat,snvs_cohort,wisecondorx" --params-file my_params.yml
```

#### Testing and Debugging

```bash
# Dry run (check configuration without execution)
./run_pipeline.sh --profile slurm,apptainer --dry-run --params-file my_params.yml

# Stub run (fast execution for testing workflow logic)
./run_pipeline.sh --profile test --stub-run --params-file my_params.yml
```

#### Resume Failed Runs

```bash
./run_pipeline.sh --profile slurm,apptainer --resume --params-file my_params.yml
```

## Configuration

### Pipeline Architecture

The pipeline operates at the **family level** and automatically determines what needs to be done:

1. **Reads pedigree file** to identify families and individuals
2. **Checks existing outputs** to determine what steps are needed
3. **Validates dependencies** - ensures required steps are available
4. **Runs only necessary processes** based on missing files

### Pipeline Steps

The pipeline supports nine main steps that must be explicitly listed in the `steps` parameter:

**Available Steps:**

- **`alignment`**: BWA-MEM2 alignment from FASTQ files and realignment from CRAM files. Produces CRAM files and coverage bedgraph files.
- **`deepvariant_sample`**: Individual variant calling using Google DeepVariant (3-stage process). Produces gVCF and VCF files, plus VAF bedgraph files.
- **`deepvariant_family`**: Family-based joint calling using GLnexus, followed by normalization and pedigree extraction. Produces family VCF/BCF files and family-specific pedigree files.
- **`annotation`**: Multi-step annotation workflow including gnomAD frequency annotation, rare/common variant filtering, VEP annotation, and additional bcftools annotations. Produces separate rare and common variant files.
- **`wombat`**: PyWombat variant filtering and prioritization using custom YAML configurations. Converts BCF to TSV and runs user-defined filtering rules.
- **`snvs_cohort`**: Cohort-level merging of common variants and Wombat results across all families.
- **`wisecondorx`**: CNV/SV calling using WisecondorX. Includes NPZ conversion, prediction, family/cohort merging, and gene annotation.
- **`extractor`**: Extract specific variants from TSV lists across family BCFs, Wombat outputs, or individual gVCFs. Supports GRCh37→GRCh38 liftover.
- **`ancestry`**: Genetic ancestry and polygenic scores using [ancestry-pgs](https://github.com/bourgeron-lab/ancestry-pgs). Genotypes the reference panel sites directly from each sample's gVCF, merges them per family, then projects each family onto the panel and scores the PGS catalog. Cohort tables are concatenations of the family tables.

**Step Dependencies:**

- `annotation` requires normalized family BCF files (triggers `deepvariant_family` if missing)
- `deepvariant_family` requires gVCF files (triggers `deepvariant_sample` if missing)
- `deepvariant_sample` requires CRAM files (triggers `alignment` if missing)
- `wombat` requires annotated BCF files (triggers `annotation` if missing)
- `snvs_cohort` requires common filtered BCFs and/or Wombat outputs (triggers upstream steps if missing)
- `wisecondorx` requires CRAM files (triggers `alignment` if missing)
- `extractor` requires normalized BCFs, Wombat outputs, or gVCFs depending on extraction mode
- `ancestry` requires gVCF files only (triggers `deepvariant_sample` if missing). It does **not** use the normalized, annotated or cohort-merged call sets, so it can be run on its own with `steps: ["ancestry"]` even on a cohort whose annotation is incomplete
- The pipeline will error if required steps are not listed in parameters

### Workflow Details

#### Alignment Workflow

Handles BWA-MEM2 alignment from FASTQ files and realignment from existing CRAM files. Supports:
- Alignment from paired-end FASTQ files with multiple sequencing units per sample
- Realignment from GRCh37 or GRCh38 CRAM files using Bazam
- Automatic unit merging for samples with multiple sequencing runs
- Generation of coverage bedgraph files (binned by `bin` parameter)

**Outputs:** CRAM files, CRAM indices, coverage bedgraph files

#### DeepVariant Sample Workflow

Runs individual variant calling using Google DeepVariant's 3-stage pipeline:
1. `make_examples` - Extract candidate variant regions
2. `call_variants` - Call variants using deep learning model
3. `postprocess_variants` - Generate final VCF and gVCF files

Also generates VAF (Variant Allele Frequency) bedgraph files for visualization.

**Outputs:** VCF files, gVCF files, VAF bedgraph files

#### DeepVariant Family Workflow

Performs family-based joint calling and post-processing:
1. Run GLnexus to jointly call variants across family members
2. Normalize variants with bcftools (left-align and normalize indels)
3. Extract family-specific pedigree from main pedigree file

**Outputs:** Normalized BCF files, family pedigree files

#### Annotation Workflow

Multi-step annotation and filtering pipeline:
1. **gnomAD annotation** - Add population frequency information
2. **Rare/Common filtering** - Split variants based on frequency threshold
3. **Common variant processing** - Keep only GT field for common variants
4. **VEP annotation** - Annotate rare variants with Ensembl VEP
5. **Additional annotations** - Add custom annotations with bcftools

**Outputs:** Rare VCF (pre-VEP), rare VCF (VEP annotated), rare BCF (fully annotated), common BCF, common BCF (GT only)

#### Wombat Workflow

PyWombat-based variant filtering and prioritization:
1. **BCF to TSV conversion** - Convert annotated BCF to TSV format
2. **PyWombat filtering** - Apply user-defined filtering rules from YAML configurations

Supports multiple configuration files for different filtering strategies (e.g., rare high-impact variants, de novo mutations, loss-of-function variants).

**Outputs:** TSV.gz file (BCF converted), filtered TSV files (one per config)

#### SNVs Cohort Workflow

Cohort-level merging of results:
1. **Common variants merge** - Merge common variant BCFs across all families
2. **Wombat results merge** - Concatenate Wombat filtered results across families

**Outputs:** Cohort-level common variant BCF, cohort-level Wombat TSV files

#### WisecondorX Workflow

CNV/SV calling using WisecondorX:
1. **NPZ conversion** - Convert CRAM to NPZ format for WisecondorX
2. **Predict** - Call CNVs using WisecondorX predict
3. **Chr reformat** - Reformat chromosome names
4. **Family merge** - Merge aberrations within families
5. **Annotation** - Annotate aberrations with gene and exon information
6. **Cohort merge** - Merge aberrations across all families

**Outputs:** NPZ files, individual aberrations BED, family aberrations BED (annotated), cohort aberrations BED

#### Extractor Workflow

Extract specific variants from TSV lists across the cohort:
1. **Refactor TSV** - Parse and validate input TSV (supports GRCh37/38 with liftover)
2. **Extract from sources**:
   - Family BCFs (normalized joint calls)
   - Family TSVs (Wombat annotated variants)
   - Individual gVCFs (sample-specific calls)
3. **Aggregate** - Combine results per family and across cohort

Useful for validating specific variants, extracting variants of interest, or comparing calls across sources.

**Outputs:** Extracted variants TSV, per-family aggregated TSV, cohort-aggregated TSV

#### Ancestry / PGS Workflow

Genetic ancestry and polygenic scores per family, and a cohort-level aggregate:
1. **Panel sites** - Verify the reference bundle and build the union of the two site lists it carries (once per run)
2. **Panel extraction** - For each sample, genotype the panel sites straight from the DeepVariant gVCF
3. **Family merge** - Join the family's per-sample panel genotypes into one family BCF
4. **Scoring** - Per family: projected principal components with an ancestry label, admixture proportions, raw polygenic scores, then the ancestry-adjusted scores and z-scores
5. **Cohort merge** - Concatenate the family tables into cohort tables, with a `family_id` column

**Why the extraction reads gVCFs.** The two site lists are not nested: the LD-pruned
ancestry panel used for the components and admixture, and the PGS catalog list, overlap
only partially, so a single extraction covers their union. More importantly, a gVCF
carries reference blocks, so a panel site with coverage and no variant is reported as
`0/0` while a site with no coverage stays `./.`. A family's own `common_gt.bcf` holds
only the sites where the family carries an alt allele — roughly 57% of the panel for a
trio, less for a duo — which is well below the 90% per-sample coverage that the
admixture projection requires, and enough to distort the projected components. Reading
the gVCF recovers the rest.

Two consequences worth knowing:

- **Results do not depend on cohort composition.** Each sample is genotyped against the
  reference bundle alone, so adding a family never changes another family's components
  or scores.
- **Cohort tables are exact concatenations.** Every value is a per-sample projection
  against the bundle, so a sample's row is identical whether its family was scored alone
  or as part of a cohort-wide run — there is nothing to recompute at cohort level.

The `ancestry_panel_name` label is part of every output file name, and since the
pipeline decides what to recompute from what exists on disk, that label is the only
thing that invalidates earlier results. Bump it whenever the depth/quality thresholds
or the reference bundle change.

The extraction logic ships with a self-contained branch-coverage test that needs only
`bcftools` and `python3` — it builds a synthetic gVCF exercising every decision branch
and asserts the genotype expected at each panel site:

```bash
modules/ancestry/scripts/panel_genotype_test
```

**Outputs:** Per-sample panel genotype BCF and call-rate stats, family panel genotype BCF, per-family PCs / ancestry labels / admixture proportions / raw, adjusted and z-scored PGS with a QC JSON per command, and the cohort-level concatenation of each table

### Parameters File (params.yml)

Every parameter the workflow reads is documented in
**[documentation/params.md](documentation/params.md)** — what each key does, its default,
which step reads it, and which keys found in older parameters files no longer do anything.

The smallest file that runs:

```yaml
steps: ["alignment"]

data: "/path/to/your/data/"
scratch: "/path/to/your/scratch/"
# The pedigree defaults to
# <data>/cohorts/my_cohort/my_cohort.pedigree.tsv
cohort_name: "my_cohort"

ref: "/path/to/reference/genome.fa"
ref_name: "GRCh38_GIABv3"

# One input source is required: fastq_pattern, old_cram_38 or old_cram_37
fastq_pattern: "*_R{1,2}.fastq.gz"
```

Realistic, complete files live in [`params_example/`](params_example/). Configure a run from
this file alone — parameters passed on the command line leave no trace in it, which is what
the cohort's state record is checksumming against. One thing that catches people out:

- Several keys are baked into output filenames (`ref_name`, `bin`, `vep_config_name`,
  `wisecondorx_binsize`, `ancestry_panel_name`, `cohort_name`). Changing one makes the
  pipeline recompute; changing a setting that is *not* in a filename leaves existing results
  in place and reuses them.

### Smart File Detection

The pipeline automatically detects existing files and skips unnecessary work. All sample and family paths are resolved using the two-level shard scheme (`{S1}/{S2}` computed from the entity ID — see [Sharding scheme](#sharding-scheme)).

- **CRAM files**: `${data}/samples/{S1}/{S2}/${barcode}/sequences/${barcode}.${ref_name}.cram` (and `.crai`)
- **Coverage bedgraph files**: `${data}/samples/{S1}/{S2}/${barcode}/sequences/${barcode}.by${bin}.bedgraph.gz` (and `.tbi`)
- **Individual gVCF files**: `${data}/samples/{S1}/{S2}/${barcode}/deepvariant/${barcode}.g.vcf.gz` (and `.tbi`)
- **VAF bedgraph files**: `${data}/samples/{S1}/{S2}/${barcode}/sequences/${barcode}.vaf.bedgraph.gz` (and `.tbi`)
- **Normalized family BCF files**: `${data}/families/{S1}/{S2}/${FID}/vcfs/${FID}.norm.bcf` (and `.csi`)
- **Family pedigree files**: `${data}/families/{S1}/{S2}/${FID}/${FID}.pedigree.tsv`
- **Rare variant VCF files**: `${data}/families/{S1}/{S2}/${FID}/vcfs/${FID}.rare.vcf.gz` (and `.tbi`)
- **Common variant BCF files**: `${data}/families/{S1}/{S2}/${FID}/vcfs/${FID}.common.bcf` and `${FID}.common_gt.bcf` (with `.csi`)
- **VEP annotated VCF files**: `${data}/families/{S1}/{S2}/${FID}/vcfs/${FID}.rare.${vep_config_name}.vcf.gz` (and `.tbi`)
- **Fully annotated BCF files**: `${data}/families/{S1}/{S2}/${FID}/vcfs/${FID}.rare.${vep_config_name}.annotated.bcf` (and `.csi`)
- **Wombat TSV files**: `${data}/families/{S1}/{S2}/${FID}/wombat/${FID}.rare.${vep_config_name}.annotated.${config_name}.tsv`
- **WisecondorX NPZ files**: `${data}/samples/{S1}/{S2}/${barcode}/svs/wisecondorx/${barcode}.${wisecondorx_binsize}.npz`
- **WisecondorX aberrations**: `${data}/samples/{S1}/{S2}/${barcode}/svs/wisecondorx/${barcode}_aberrations.chr.bed`
- **Cohort BCF files**: `${data}/cohorts/${cohort_name}/vcfs/${cohort_name}.common_gt.bcf` (and `.csi`)
- **Panel genotype files (sample)**: `${data}/samples/{S1}/{S2}/${barcode}/ancestry/${barcode}.panel_gt.${ancestry_panel_name}.bcf` (and `.csi`)
- **Panel genotype files (family)**: `${data}/families/{S1}/{S2}/${FID}/ancestry/${FID}.panel_gt.${ancestry_panel_name}.bcf` (and `.csi`)
- **Ancestry/PGS tables (family)**: `${data}/families/{S1}/{S2}/${FID}/ancestry/${FID}.${ancestry_panel_name}.{pcs,ancestry,Q,pgs_raw,pgs_adjusted,pgs_zscore}.tsv`
- **Ancestry/PGS tables (cohort)**: `${data}/cohorts/${cohort_name}/ancestry/${cohort_name}.${ancestry_panel_name}.{pcs,ancestry,Q,pgs_raw,pgs_adjusted,pgs_zscore}.tsv`

If these files exist with their indices (where applicable), the corresponding steps are skipped.

## Input Data Structure

The pipeline expects the following directory structure:

### Sharding scheme

`samples/` and `families/` use a **two-level sharding** scheme to avoid large flat directories. Given an entity ID:

1. Strip all `-`, `.`, and `_` characters.
2. **Shard 1** = last character of the stripped ID (uppercased).
3. **Shard 2** = second-to-last character (uppercased).

Examples: `C000EZJ` → `samples/J/Z/C000EZJ/`; `C0733-011-068` → `families/8/6/C0733-011-068/`

Cohort directories are **not** sharded.

```
data/
├── fastq/                         # FASTQ files for alignment (optional)
│   ├── A001_DA_SAMPLE1_L001_1_001.HG7T2.dual.fastq.gz
│   ├── A001_DA_SAMPLE1_L001_2_001.HG7T2.dual.fastq.gz
│   └── ...
├── cram4realignment/              # Existing CRAM files for realignment (optional)
│   ├── SAMPLE1.cram
│   ├── SAMPLE1.cram.crai
│   └── ...
├── samples/                       # Sample-specific output directories (sharded)
│   └── {S1}/                      # Shard level 1 (single character)
│       └── {S2}/                  # Shard level 2 (single character)
│           └── BC001/             # Sample barcode directory
│               ├── sequences/     # CRAM and bedgraph files
│               │   ├── BC001.GRCh38_GIABv3.cram
│               │   ├── BC001.GRCh38_GIABv3.cram.crai
│               │   ├── BC001.by1000.bedgraph.gz        # Coverage bedgraph
│               │   ├── BC001.by1000.bedgraph.gz.tbi
│               │   ├── BC001.vaf.bedgraph.gz           # VAF bedgraph
│               │   └── BC001.vaf.bedgraph.gz.tbi
│               ├── deepvariant/   # DeepVariant outputs
│               │   ├── BC001.g.vcf.gz
│               │   ├── BC001.g.vcf.gz.tbi
│               │   ├── BC001.vcf.gz
│               │   └── BC001.vcf.gz.tbi
│               ├── svs/           # SV calling outputs
│               │   └── wisecondorx/
│               │       ├── BC001.npz
│               │       ├── BC001_aberrations.bed
│               │       └── BC001_aberrations.chr.bed
│               └── ancestry/      # Panel genotypes for ancestry/PGS
│                   ├── BC001.panel_gt.apgs_b1.0.0_dp10gq20.bcf
│                   ├── BC001.panel_gt.apgs_b1.0.0_dp10gq20.bcf.csi
│                   └── BC001.panel_gt.apgs_b1.0.0_dp10gq20.stats.tsv   # Per-sample call rate
├── families/                      # Family-specific output directories (sharded)
│   └── {S1}/                      # Shard level 1 (single character)
│       └── {S2}/                  # Shard level 2 (single character)
│           └── FID001/            # Family directory
│               ├── FID001.pedigree.tsv    # Family-specific pedigree
│               ├── vcfs/
│               │   ├── FID001.norm.bcf                                      # Normalized family BCF
│               │   ├── FID001.norm.bcf.csi
│               │   ├── FID001.rare.vcf.gz                                   # Rare variants (pre-VEP)
│               │   ├── FID001.rare.vcf.gz.tbi
│               │   ├── FID001.common.bcf                                    # Common variants
│               │   ├── FID001.common.bcf.csi
│               │   ├── FID001.common_gt.bcf                                 # Common variants (GT only)
│               │   ├── FID001.common_gt.bcf.csi
│               │   ├── FID001.rare.ensembl_vep_115.vcf.gz                   # VEP annotated
│               │   ├── FID001.rare.ensembl_vep_115.vcf.gz.tbi
│               │   ├── FID001.rare.ensembl_vep_115.annotated.bcf            # Fully annotated
│               │   └── FID001.rare.ensembl_vep_115.annotated.bcf.csi
│               ├── wombat/
│               │   ├── FID001.rare.ensembl_vep_115.annotated.tsv.gz                        # BCF to TSV
│               │   ├── FID001.rare.ensembl_vep_115.annotated.de_novo_mutations.tsv         # Wombat filtered
│               │   └── FID001.rare.ensembl_vep_115.annotated.rare_variants_high_impact.tsv
│               ├── svs/
│               │   └── wisecondorx/
│               │       ├── FID001_aberrations.bed           # Family merged
│               │       └── FID001_aberrations.annotated.bed # Gene annotated
│               └── ancestry/
│                   ├── FID001.panel_gt.apgs_b1.0.0_dp10gq20.bcf        # Family panel genotypes
│                   ├── FID001.panel_gt.apgs_b1.0.0_dp10gq20.bcf.csi
│                   ├── FID001.apgs_b1.0.0_dp10gq20.pcs.tsv             # Projected components
│                   ├── FID001.apgs_b1.0.0_dp10gq20.ancestry.tsv        # Region/population label
│                   ├── FID001.apgs_b1.0.0_dp10gq20.Q.tsv               # Admixture proportions
│                   ├── FID001.apgs_b1.0.0_dp10gq20.pgs_raw.tsv         # Raw scores
│                   ├── FID001.apgs_b1.0.0_dp10gq20.pgs_adjusted.tsv    # Ancestry-adjusted
│                   ├── FID001.apgs_b1.0.0_dp10gq20.pgs_zscore.tsv      # Z-scored
│                   └── FID001.apgs_b1.0.0_dp10gq20.*.qc.json           # One QC report per command
├── cohorts/                       # Cohort-specific directories (not sharded)
│   └── COHORT_NAME/
│       ├── COHORT_NAME.pedigree.tsv                 # Family structure (required)
│       ├── COHORT_NAME.params.yml                   # This cohort's parameters
│       ├── vcfs/
│       │   ├── COHORT_NAME.common_gt.bcf            # Cohort common variants
│       │   └── COHORT_NAME.common_gt.bcf.csi
│       ├── wombat/
│       │   ├── COHORT_NAME.rare.ensembl_vep_115.annotated.de_novo_mutations.results.tsv
│       │   └── COHORT_NAME.rare.ensembl_vep_115.annotated.rare_variants_high_impact.results.tsv
│       ├── svs/
│       │   └── wisecondorx/
│       │       └── COHORT_NAME_aberrations.bed      # Cohort merged aberrations
│       └── ancestry/
│           ├── apgs_b1.0.0_dp10gq20.sites.tsv.gz    # Union panel site list
│           ├── apgs_b1.0.0_dp10gq20.regions.tsv.gz  # bcftools targets file
│           ├── apgs_b1.0.0_dp10gq20.bundle.json     # Reference bundle check
│           └── COHORT_NAME.apgs_b1.0.0_dp10gq20.{pcs,ancestry,Q,pgs_raw,pgs_adjusted,pgs_zscore}.tsv
└── extractor/                     # Extractor outputs (if TSV lists provided)
    └── VARIANT_LIST/
        ├── VARIANT_LIST.extracted.tsv               # Extracted variants
        ├── VARIANT_LIST.aggregated.tsv              # Fully aggregated
        └── families/
            ├── FID001.extracted.tsv
            └── ...
```

### FASTQ File Naming Convention

The pipeline expects FASTQ files to follow this naming pattern:

```
A{PROJECT}_DA_{BARCODE}_{LANE}_{READ}_{INDEX}.{FLOWCELL}.{DUAL}.fastq.gz
```

Example: `A001_DA_BC001_L001_1_001.HG7T2.dual.fastq.gz`

Where:

- `{BARCODE}`: Sample identifier (must match pedigree file)
- `{LANE}`: Sequencing lane
- `{READ}`: 1 or 2 for paired-end reads
- `{FLOWCELL}`: Flowcell identifier
- `{DUAL}`: Dual index information

## Output

The pipeline generates:

All sample and family output paths include two shard levels (`{S1}/{S2}`) computed from the entity ID. See [Sharding scheme](#sharding-scheme).

### Alignment Outputs

- **Final CRAM files**: `${data}/samples/{S1}/{S2}/${barcode}/sequences/${barcode}.${ref_name}.cram`
- **CRAM indices**: `${data}/samples/{S1}/{S2}/${barcode}/sequences/${barcode}.${ref_name}.cram.crai`
- **Bedgraph files**: `${data}/samples/{S1}/{S2}/${barcode}/sequences/${barcode}.by${bin}.bedgraph.gz`
- **Bedgraph indices**: `${data}/samples/{S1}/{S2}/${barcode}/sequences/${barcode}.by${bin}.bedgraph.gz.tbi`

### DeepVariant Sample Outputs

- **Individual VCF files**: `${data}/samples/{S1}/{S2}/${barcode}/deepvariant/${barcode}.vcf.gz`
- **Individual gVCF files**: `${data}/samples/{S1}/{S2}/${barcode}/deepvariant/${barcode}.g.vcf.gz`
- **VCF indices**: `*.vcf.gz.tbi` and `*.g.vcf.gz.tbi`
- **VAF bedgraph files**: `${data}/samples/{S1}/{S2}/${barcode}/sequences/${barcode}.vaf.bedgraph.gz`
- **VAF bedgraph indices**: `${data}/samples/{S1}/{S2}/${barcode}/sequences/${barcode}.vaf.bedgraph.gz.tbi`

### DeepVariant Family Outputs

- **Normalized BCF files**: `${data}/families/{S1}/{S2}/${FID}/vcfs/${FID}.norm.bcf`
- **Normalized BCF indices**: `${data}/families/{S1}/{S2}/${FID}/vcfs/${FID}.norm.bcf.csi`
- **Family pedigree files**: `${data}/families/{S1}/{S2}/${FID}/${FID}.pedigree.tsv`

### Annotation Outputs

- **Rare variant VCF files** (pre-VEP): `${data}/families/{S1}/{S2}/${FID}/vcfs/${FID}.rare.vcf.gz`
- **Common variant BCF files**: `${data}/families/{S1}/{S2}/${FID}/vcfs/${FID}.common.bcf`
- **Common variant BCF (GT only)**: `${data}/families/{S1}/{S2}/${FID}/vcfs/${FID}.common_gt.bcf`
- **VEP annotated rare VCF files**: `${data}/families/{S1}/{S2}/${FID}/vcfs/${FID}.rare.${vep_config_name}.vcf.gz`
- **Fully annotated BCF files**: `${data}/families/{S1}/{S2}/${FID}/vcfs/${FID}.rare.${vep_config_name}.annotated.bcf`
- **All indices**: `*.tbi` for VCF.gz and `*.csi` for BCF files

### Wombat Outputs

- **BCF to TSV conversion**: `${data}/families/{S1}/{S2}/${FID}/wombat/${FID}.rare.${vep_config_name}.annotated.tsv.gz`
- **PyWombat filtered results**: `${data}/families/{S1}/{S2}/${FID}/wombat/${FID}.rare.${vep_config_name}.annotated.${config_name}.tsv`
  - One file per configuration in `wombat_config_list`

### SNVs Cohort Outputs

- **Cohort common variants BCF**: `${data}/cohorts/${cohort_name}/vcfs/${cohort_name}.common_gt.bcf`
- **Cohort Wombat results**: `${data}/cohorts/${cohort_name}/wombat/${cohort_name}.rare.${vep_config_name}.annotated.${config_name}.results.tsv`
  - One file per configuration in `wombat_config_list`

### WisecondorX Outputs

- **Individual NPZ files**: `${data}/samples/{S1}/{S2}/${barcode}/svs/wisecondorx/${barcode}.npz`
- **Individual aberrations**: `${data}/samples/{S1}/{S2}/${barcode}/svs/wisecondorx/${barcode}_aberrations.bed`
- **Individual aberrations (chr format)**: `${data}/samples/{S1}/{S2}/${barcode}/svs/wisecondorx/${barcode}_aberrations.chr.bed`
- **Family merged aberrations**: `${data}/families/{S1}/{S2}/${FID}/svs/wisecondorx/${FID}_aberrations.bed`
- **Family annotated aberrations**: `${data}/families/{S1}/{S2}/${FID}/svs/wisecondorx/${FID}_aberrations.annotated.bed`
- **Cohort merged aberrations**: `${data}/cohorts/${cohort_name}/svs/wisecondorx/${cohort_name}_aberrations.bed`

### Extractor Outputs

- **Extracted variants**: `${data}/extractor/${original_filename}/${original_filename}.extracted.tsv`
- **Aggregated per family**: `${data}/extractor/${original_filename}/families/${FID}.extracted.tsv`
- **Fully aggregated**: `${data}/extractor/${original_filename}.aggregated.tsv`

### Ancestry / PGS Outputs

Paths use `${P}` as shorthand for `${ancestry_panel_name}`.

- **Panel site list**: `${data}/cohorts/${cohort_name}/ancestry/${P}.sites.tsv.gz` (plus `${P}.regions.tsv.gz` and the `${P}.bundle.json` bundle check)
- **Per-sample panel genotypes**: `${data}/samples/{S1}/{S2}/${barcode}/ancestry/${barcode}.panel_gt.${P}.bcf` (and `.csi`)
- **Per-sample call-rate stats**: `${data}/samples/{S1}/{S2}/${barcode}/ancestry/${barcode}.panel_gt.${P}.stats.tsv`
- **Family panel genotypes**: `${data}/families/{S1}/{S2}/${FID}/ancestry/${FID}.panel_gt.${P}.bcf` (and `.csi`)
- **Family tables**: `${data}/families/{S1}/{S2}/${FID}/ancestry/${FID}.${P}.{pcs,ancestry,Q,pgs_raw,pgs_adjusted,pgs_zscore}.tsv`
- **Family QC reports**: `${data}/families/{S1}/{S2}/${FID}/ancestry/${FID}.${P}.{pcs,admixture,pgs-raw,pgs-adjusted,pgs-zscore}.qc.json`
- **Cohort tables**: `${data}/cohorts/${cohort_name}/ancestry/${cohort_name}.${P}.{pcs,ancestry,Q,pgs_raw,pgs_adjusted,pgs_zscore}.tsv`
  - Concatenations of the family tables with a `family_id` column appended

**What to check in the QC reports.** These are the fields that catch a quietly degraded
result:

- `site_coverage` — the fraction of panel sites present; should be at or near 1.0, since the extraction emits a record for every site
- `min_per_sample_coverage` — the worst per-sample call rate; must clear 0.90 or `admixture` refuses
- `n_outliers` — samples lying outside the reference panel's coverage in component space, labelled `outlier` rather than assigned a population
- `median_knn_distance_to_reference` — how much reference data supports these z-scores
- `distinct_allele_ct` — **expect this to exceed 1, with its warning.** Per-sample missingness genuinely differs between samples, so the raw score sums are not directly comparable between individuals; the adjusted and z-scored tables are what you compare. This is a consequence of recording real coverage instead of assuming every uncalled site is homozygous reference.
- `n_traits_constant` — catalog columns carrying no non-zero weight at any scored site, emitted as `nan` rather than a number

### Pipeline Reports

- **Execution timeline**: `reports/timeline.html`
- **Execution report**: `reports/report.html`
- **Process trace**: `reports/trace.txt`
- **Workflow diagram**: `reports/dag.svg`

## Advanced Usage

### Pipeline Analysis Summary

Before running any processes, the pipeline displays:

```
========================================================================================
                                ANALYSIS SUMMARY
========================================================================================
ALIGNMENT: 2 individuals done, 1 to align, 0 needing bedgraph only
== SNVs/INDELs Calling ==
DEEPVARIANT_SAMPLE: 1 individuals done and 2 to do
DEEPVARIANT_FAMILY: 0 families done and 2 to do
ANNOTATION: 0 families done and 2 to do
WOMBAT: 0 families done and 2 to do
== Common Variants ==
SNVS_COHORT: common variants cohort bcf merge: Yes - wombat cohort merges due: 1
== SVs Calling ==
WISECONDORX PREDICT: 0 individuals done and 3 to do
== Ancestry / PGS ==
ANCESTRY: Skipped (step not requested)
== Other ==
EXTRACTOR: Skipped (no TSV files provided)
========================================================================================
```

This helps you understand what work will be performed.

### Cohort Run State

After every run, the pipeline writes a hidden record of what it did into the cohort directory:

```
${data}/cohorts/<COHORT_NAME>/.ghfc-ngs.state.json
```

It holds the last run and the last successful run, each with a timestamp, a completion status,
SHA-256 checksums of the pedigree and parameters file that produced the outputs, the pipeline
version and commit, and per-step completion percentages - plus a short history of earlier runs.

Two further blocks describe the cohort's inputs rather than its outputs.
`samples_without_cram` lists every individual with no CRAM, with its family, whether it already
has a gVCF, and which input source - if any - could align it; its `blocked` count is the number
that can neither progress nor be fixed from what is on disk. `stale_family_clean` records what
a `--clean-stale-families` run removed, and which families it refused and why.

```bash
# Did the last run finish, and how complete is the cohort?
jq -r '.last_run.status' "$data/cohorts/EAGER/.ghfc-ngs.state.json"
jq -r '.last_run.completion' "$data/cohorts/EAGER/.ghfc-ngs.state.json"

# Which samples are stuck, and what would unstick them?
jq -r '.last_run.samples_without_cram.samples[] | select(.blocked)
       | "\(.barcode)\t\(.family_id)\t\(.input_source)"' \
   "$data/cohorts/EAGER/.ghfc-ngs.state.json" | column -t -s $'\t'
```

Stub and preview runs deliberately write nothing, and the file is skipped when `cohort_name` is
not set. See [COHORT_STATE.md](COHORT_STATE.md) for the schema, the status lifecycle and more
recipes.

### Error Checking

The pipeline validates dependencies and will stop with clear error messages:

```bash
ERROR: DeepVariant sample step is required for 3 individuals but not included in steps parameter

Please add the required steps to your parameters or ensure all required files exist.
Available steps: alignment, deepvariant_sample, deepvariant_family, annotation, snvs_cohort, wisecondorx, wombat, extractor
```

### Custom Configuration

You can override any configuration parameter:

```bash
nextflow run main.nf \
    -profile slurm,apptainer \
    --data /my/data \
    --scratch /my/scratch \
    --ref /my/reference.fa \
    --ref_name "MyRef" \
    --cohort_name "my_cohort" \
    --steps "alignment,deepvariant_sample,deepvariant_family,annotation"
```

### Resume Functionality

Nextflow's resume feature works seamlessly with the family-aware logic:

```bash
./run_pipeline.sh --profile slurm,apptainer --resume --params-file my_params.yml
```

Only processes that haven't completed successfully will be re-run.

## Container Configuration

### Default Containers

The pipeline uses the following containers by default:

- **Alignment**: `fcliquet/bioinfo-swissknife:latest` (BWA-MEM2, SAMtools, etc.)
- **DeepVariant**: `google/deepvariant:1.9.0`
- **GLnexus**: `cgrlab/glnexus:v1.4.1`
- **SAMtools**: `biocontainers/samtools:1.19.2--h50ea8bc_1`

### Container Profiles

- **`apptainer`**: Use Apptainer containers
- **`singularity`**: Use Singularity containers (deprecated, use apptainer)
- **Native execution**: Use local tools (not recommended)

## Resource Configuration

### Default Resource Allocations

- **BWA_MEM2_ALIGN**: 95 CPUs, 460GB RAM, 240h
- **BAZAM_BWA_MEM2_REALIGN**: 47 CPUs, 230GB RAM, 240h
- **DV_MAKE_EXAMPLES**: 95 CPUs, 460GB RAM, 240h
- **DV_CALL_VARIANTS**: 95 CPUs, 460GB RAM, 240h
- **DV_POSTPROCESS_VARIANTS**: 95 CPUs, 460GB RAM, 240h
- **GLNEXUS_FAMILY**: 4 CPUs, 50GB RAM, 240h
- **MERGE_UNITS**: 1 CPU, 10GB RAM, 240h
- **INDEX_CRAM**: 1 CPU, 5GB RAM, 240h

### Resource Tuning

Adjust resources in `nextflow.config`:

```groovy
process {
    withName: 'BWA_MEM2_ALIGN' {
        cpus = 64        // Reduce if nodes have fewer cores
        memory = '250.GB' // Adjust based on available memory
    }
    
    withName: 'GLNEXUS_FAMILY' {
        cpus = 8         // Scale based on family size
        memory = '100.GB' // GLnexus can be memory intensive
    }
}
```

## SLURM Configuration

### Account and Partition Settings

SLURM account and partition are set in `nextflow.config`, not in the parameters file:

```groovy
process {
    clusterOptions = '-p ghfc --qos=ghfc --account=your_account'
}
```

> `slurm_account` and `slurm_partition` appear in some parameters files but are **not read by
> anything** — the `clusterOptions` strings above are the only place these are configured.
> See [documentation/params.md](documentation/params.md#keys-that-look-live-but-are-inert).

### Queue Management

```groovy
executor {
    queueSize = 50          // Maximum jobs in queue
    submitRateLimit = '10 sec' // Job submission rate
}
```

## GLnexus Configuration

### Configuration Presets

GLnexus supports different configuration presets:

```yaml
glnexus_config: "DeepVariant_unfiltered"  # Default
# glnexus_config: "DeepVariant"            # With quality filters
# glnexus_config: "gatk"                   # For GATK gVCFs
```

### Family Size Considerations

- **Small families (2-4 individuals)**: Default resources usually sufficient
- **Large families (>10 individuals)**: May need increased memory and time
- **Very large families (>50 individuals)**: Consider splitting or specialized configuration

## Testing

### Test Profile

```bash
# Quick test with reduced resources
./run_pipeline.sh --profile test --stub-run --data test_data/
```

### Stub Run Features

- **Fast execution**: Creates empty output files instead of running tools
- **Workflow validation**: Ensures all file paths and dependencies are correct
- **Resource efficient**: Minimal CPU/memory usage
- **Error detection**: Catches configuration issues quickly

## Troubleshooting

### Common Issues

1. **Missing pedigree file**: Ensure pedigree.tsv exists and is properly formatted
2. **Step dependency errors**: Add required steps to the `steps` parameter using correct step names (`deepvariant_sample`, `deepvariant_family`, `annotation`, etc.)
3. **File not found errors**: Check file paths and ensure barcode names match between pedigree and data files
4. **Out of memory errors**: Increase memory allocation for memory-intensive processes (especially GLnexus, PyWombat, WisecondorX)
5. **GLnexus failures**: Often memory-related, try increasing memory allocation
6. **Wombat config not found**: Ensure `wombat_config_path` points to directory containing YAML files and filenames in `wombat_config_list` are correct
7. **WisecondorX reference missing**: Ensure `wisecondorx_reference` points to valid NPZ reference file
8. **gnomAD annotation failures**: Verify `gnomad_file` path exists and is properly indexed (.csi file)
9. **VEP failures**: Check `vep_config` INI file exists and VEP cache is properly configured
10. **Cohort merge failures**: Ensure `cohort_name` is set when running `snvs_cohort` or `wisecondorx` workflows

### Debug Information

```bash
# Run with detailed logging
nextflow run main.nf -profile slurm,apptainer --params-file params.yml -with-trace -with-timeline -with-dag dag.png
```

### Log Files

- `.nextflow.log`: Main Nextflow execution log
- `work/`: Individual task logs and intermediate files
- `reports/`: Execution reports and timeline

### File Naming Issues

Ensure that:

- Barcode names in pedigree file match those in FASTQ filenames
- Reference name (`ref_name`) is consistent across runs
- File patterns in parameters match your actual file structure

## Migration from Previous Versions

### From Individual-based Pipeline

The new family-based pipeline:

- **Requires a pedigree file** (wasn't needed before)
- **Processes families together** (not individual samples)
- **Automatically manages dependencies** (no manual step skipping)
- **Produces family VCF files** (new output type)

### Key Changes

1. **Add pedigree file**: Create TSV file with family structure
2. **Update step names**: Use new step names (`deepvariant_sample`, `deepvariant_family`, `annotation` instead of `deepvariant`, `family_calling`, `vep_annotation`)
3. **Update parameters**: Remove `skip_*` options, add all required steps to `steps` list. Add new required parameters like `cohort_name`, `gnomad_file`, `vep_config_name`, etc.
4. **Update paths**: Family VCFs now use BCF format in many cases, check output structure
5. **Review naming**: CRAM files now include `ref_name` in filename
6. **File format changes**: Many intermediate files now use BCF format instead of VCF.gz

## Performance Optimization

### For Large Cohorts

1. **Increase queue size**: Allow more parallel jobs
2. **Optimize GLnexus memory**: Scale based on largest family size
3. **Use fast storage**: Place scratch directory on fast filesystem
4. **Monitor resource usage**: Use execution reports to optimize allocations

### For Small Cohorts

1. **Reduce default resources**: Lower CPU/memory allocations
2. **Use test profile**: For development and small datasets
3. **Enable stub runs**: For rapid pipeline testing

## Support

For issues and questions:

1. Check the troubleshooting section above
2. Review execution reports in `reports/` directory
3. Check individual task logs in `work/` directory
4. Verify pedigree file format and content
5. Ensure all required steps are listed in parameters
6. Review Nextflow documentation: <https://www.nextflow.io/docs/latest/>

## Notes and Limitations

### De Novo Mutation Analysis

De novo mutations are now handled through the Wombat workflow using YAML configuration files (e.g., `de_novo_mutations.yml`). The legacy `dnm_report` script has been replaced by this more flexible configuration-based approach that allows custom filtering criteria.

### File Format Considerations

The pipeline uses BCF format for many intermediate and final files to improve performance and reduce disk space usage. BCF files can be converted to VCF with `bcftools view`:

```bash
bcftools view -O z -o output.vcf.gz input.bcf
```

### Cohort Name Requirement

When running `snvs_cohort` or generating cohort-level outputs, the `cohort_name` parameter is required to name the output directory and files.

## Citation

If you use this pipeline in your research, please cite:

- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319.
- **BWA-MEM2**: Vasimuddin, M., et al. (2019). Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. IEEE IPDPS.
- **DeepVariant**: Poplin, R., et al. (2018). A universal SNP and small-indel variant caller using deep neural networks. Nature Biotechnology, 36(10), 983-987.
- **GLnexus**: Yun, T., et al. (2021). Accurate, scalable cohort variant calls using DeepVariant and GLnexus. Bioinformatics, 37(5), 682-685.
- **VEP**: McLaren, W., et al. (2016). The Ensembl Variant Effect Predictor. Genome Biology, 17(1), 122.
- **WisecondorX**: Raman, L., et al. (2019). WisecondorX: improved copy number detection for routine shallow whole-genome sequencing. Nucleic Acids Research, 47(4), 1605-1614.
- **gnomAD**: Chen, S., et al. (2024). A genomic mutational constraint map using variation in 76,156 human genomes. Nature, 625, 92-100.
