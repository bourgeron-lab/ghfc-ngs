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

Create a TSV file with 6 columns (FID, barcode, father, mother, sex, phenotype):

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

The pipeline supports eight main steps that must be explicitly listed in the `steps` parameter:

**Available Steps:**

- **`alignment`**: BWA-MEM2 alignment from FASTQ files and realignment from CRAM files. Produces CRAM files and coverage bedgraph files.
- **`deepvariant_sample`**: Individual variant calling using Google DeepVariant (3-stage process). Produces gVCF and VCF files, plus VAF bedgraph files.
- **`deepvariant_family`**: Family-based joint calling using GLnexus, followed by normalization and pedigree extraction. Produces family VCF/BCF files and family-specific pedigree files.
- **`annotation`**: Multi-step annotation workflow including gnomAD frequency annotation, rare/common variant filtering, VEP annotation, and additional bcftools annotations. Produces separate rare and common variant files.
- **`wombat`**: PyWombat variant filtering and prioritization using custom YAML configurations. Converts BCF to TSV and runs user-defined filtering rules.
- **`snvs_cohort`**: Cohort-level merging of common variants and Wombat results across all families.
- **`wisecondorx`**: CNV/SV calling using WisecondorX. Includes NPZ conversion, prediction, family/cohort merging, and gene annotation.
- **`extractor`**: Extract specific variants from TSV lists across family BCFs, Wombat outputs, or individual gVCFs. Supports GRCh37→GRCh38 liftover.

**Step Dependencies:**

- `annotation` requires normalized family BCF files (triggers `deepvariant_family` if missing)
- `deepvariant_family` requires gVCF files (triggers `deepvariant_sample` if missing)
- `deepvariant_sample` requires CRAM files (triggers `alignment` if missing)
- `wombat` requires annotated BCF files (triggers `annotation` if missing)
- `snvs_cohort` requires common filtered BCFs and/or Wombat outputs (triggers upstream steps if missing)
- `wisecondorx` requires CRAM files (triggers `alignment` if missing)
- `extractor` requires normalized BCFs, Wombat outputs, or gVCFs depending on extraction mode
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

### Parameters File (params.yml)

Key parameters to configure:

```yaml
# Pipeline steps - ALL REQUIRED STEPS MUST BE LISTED
steps: [
    "alignment",
    "deepvariant_sample",
    "deepvariant_family",
    "annotation",
    "snvs_cohort",
    "wisecondorx",
    "wombat"
]

# Data directories
data: "/path/to/your/data/"
scratch: "/path/to/your/scratch/"

# Cohort name for merged files
cohort_name: "my_cohort"

# Pedigree file (optional - defaults to ${data}/pedigree.tsv)
pedigree: "/path/to/your/pedigree.tsv"

# Reference genomes
ref: "/path/to/reference/genome.fa"
ref_name: "GRCh38_GIABv3"  # Used in output file naming
oldref: "/path/to/old/reference/genome.fa"  # For realignment (optional)
wisecondorx_reference: "/path/to/wisecondorx/reference.npz"  # WisecondorX reference

# Bedgraph parameters
bin: 1000  # Bin size for coverage bedgraph generation

# Tool paths (when not using containers)
bwa_mem2: "bwa-mem2"
samtools: "samtools"
sambamba: "sambamba"
samblaster: "samblaster"
bazam: "/path/to/bazam.jar"

# GLnexus configuration
glnexus_config: "DeepVariant_unfiltered"

# DeepVariant configuration
deepvariant_threads: 96  # Number of parallel tasks for DeepVariant make_examples

# gnomAD frequency annotation configuration
gnomad_file: "/path/to/gnomad_v4.1_allChroms.bcf"  # Path to gnomAD annotation file
gnomad_filter_field: "AF"                           # gnomAD field to filter on (e.g., AF, AF_popmax)
gnomad_filter_threshold: "0.01"                     # Frequency threshold for filtering (e.g., 0.01 for 1%)

# VEP annotation configuration
vep_config: "/path/to/vep/config.ini"      # Path to VEP config INI file
vep_config_name: "ensembl_vep_115"         # Name suffix for output files

# Additional annotation configuration
annotation_annotation_path: "/path/to/annotations/"                     # Directory containing annotation BCF files
annotation_annotation_list: ["gnomad_v4.1_allChroms.bcf", "LCR.bed.gz"] # List of BCF/BED files for annotation
annotation_gencode: "gencode.v47.basic"                                  # Gencode version for WisecondorX aberrations annotation

# De novo mutation extraction configuration
annotation_dnm_min_callrate: "0.9"   # Minimum call rate for de novo variants
annotation_dnm_min_DP: "10"          # Minimum depth (DP) for de novo variants in child
annotation_dnm_min_GQ: "19"          # Minimum genotype quality (GQ) for de novo variants in child
annotation_dnm_min_VAF: "0.25"       # Minimum variant allele frequency (VAF) for de novo variants in child

# Wombat configuration
wombat_config_path: "/path/to/wombat/configs"                                      # Directory containing Wombat YAML configuration files
wombat_config_list: ["rare_variants_high_impact.yml", "de_novo_mutations.yml"]   # List of Wombat config files

# Extractor configuration (optional - for variant extraction from TSV lists)
extractor_tsvs_list: []                              # List of TSV files with variants to extract
liftover_chain: "/path/to/hg19ToHg38.over.chain.gz"  # Chain file for GRCh37 to GRCh38 liftover

# Pseudo-autosomal regions (PAR) coordinates for GRCh38
ref_par1_start: "10001"        # PAR1 start position on chrX
ref_par1_end: "2781479"        # PAR1 end position on chrX
ref_par2_start: "155701383"    # PAR2 start position on chrX
ref_par2_end: "156030895"      # PAR2 end position on chrX

# SLURM account and partition settings
slurm_account: "your_account"
slurm_partition: "ghfc"

# Container cache directory
apptainer_cache: "/path/to/apptainer/cache/"
singularity_cache: "/path/to/singularity/cache/"

# Resource limits
max_memory: "460.GB"
max_cpus: 95
max_time: "240.h"
```

### Smart File Detection

The pipeline automatically detects existing files and skips unnecessary work:

- **CRAM files**: `${data}/samples/${barcode}/sequences/${barcode}.${ref_name}.cram` (and `.crai`)
- **Coverage bedgraph files**: `${data}/samples/${barcode}/sequences/${barcode}.by${bin}.bedgraph.gz` (and `.tbi`)
- **Individual gVCF files**: `${data}/samples/${barcode}/deepvariant/${barcode}.g.vcf.gz` (and `.tbi`)
- **VAF bedgraph files**: `${data}/samples/${barcode}/sequences/${barcode}.vaf.bedgraph.gz` (and `.tbi`)
- **Normalized family BCF files**: `${data}/families/${FID}/vcfs/${FID}.norm.bcf` (and `.csi`)
- **Family pedigree files**: `${data}/families/${FID}/${FID}.pedigree.tsv`
- **Rare variant VCF files**: `${data}/families/${FID}/vcfs/${FID}.rare.vcf.gz` (and `.tbi`)
- **Common variant BCF files**: `${data}/families/${FID}/vcfs/${FID}.common.bcf` and `${FID}.common_gt.bcf` (with `.csi`)
- **VEP annotated VCF files**: `${data}/families/${FID}/vcfs/${FID}.rare.${vep_config_name}.vcf.gz` (and `.tbi`)
- **Fully annotated BCF files**: `${data}/families/${FID}/vcfs/${FID}.rare.${vep_config_name}.annotated.bcf` (and `.csi`)
- **Wombat TSV files**: `${data}/families/${FID}/wombat/${FID}.rare.${vep_config_name}.annotated.${config_name}.tsv`
- **WisecondorX NPZ files**: `${data}/samples/${barcode}/svs/wisecondorx/${barcode}.npz`
- **WisecondorX aberrations**: `${data}/samples/${barcode}/svs/wisecondorx/${barcode}_aberrations.chr.bed`
- **Cohort BCF files**: `${data}/cohorts/${cohort_name}/vcfs/${cohort_name}.common_gt.bcf` (and `.csi`)

If these files exist with their indices (where applicable), the corresponding steps are skipped.

## Input Data Structure

The pipeline expects the following directory structure:

```
data/
├── pedigree.tsv                   # Family structure (required)
├── fastq/                         # FASTQ files for alignment (optional)
│   ├── A001_DA_SAMPLE1_L001_1_001.HG7T2.dual.fastq.gz
│   ├── A001_DA_SAMPLE1_L001_2_001.HG7T2.dual.fastq.gz
│   └── ...
├── cram4realignment/              # Existing CRAM files for realignment (optional)
│   ├── SAMPLE1.cram
│   ├── SAMPLE1.cram.crai
│   └── ...
├── samples/                       # Sample-specific output directories
│   ├── BC001/                     # Sample barcode directory
│   │   ├── sequences/             # CRAM and bedgraph files
│   │   │   ├── BC001.GRCh38_GIABv3.cram
│   │   │   ├── BC001.GRCh38_GIABv3.cram.crai
│   │   │   ├── BC001.by1000.bedgraph.gz        # Coverage bedgraph
│   │   │   ├── BC001.by1000.bedgraph.gz.tbi
│   │   │   ├── BC001.vaf.bedgraph.gz           # VAF bedgraph
│   │   │   └── BC001.vaf.bedgraph.gz.tbi
│   │   ├── deepvariant/           # DeepVariant outputs
│   │   │   ├── BC001.g.vcf.gz
│   │   │   ├── BC001.g.vcf.gz.tbi
│   │   │   ├── BC001.vcf.gz
│   │   │   └── BC001.vcf.gz.tbi
│   │   └── svs/                   # SV calling outputs
│   │       └── wisecondorx/
│   │           ├── BC001.npz
│   │           ├── BC001_aberrations.bed
│   │           └── BC001_aberrations.chr.bed
│   └── ...
├── families/                      # Family-specific output directories
│   ├── FID001/                    # Family directory
│   │   ├── FID001.pedigree.tsv    # Family-specific pedigree
│   │   ├── vcfs/
│   │   │   ├── FID001.norm.bcf                                      # Normalized family BCF
│   │   │   ├── FID001.norm.bcf.csi
│   │   │   ├── FID001.rare.vcf.gz                                   # Rare variants (pre-VEP)
│   │   │   ├── FID001.rare.vcf.gz.tbi
│   │   │   ├── FID001.common.bcf                                    # Common variants
│   │   │   ├── FID001.common.bcf.csi
│   │   │   ├── FID001.common_gt.bcf                                 # Common variants (GT only)
│   │   │   ├── FID001.common_gt.bcf.csi
│   │   │   ├── FID001.rare.ensembl_vep_115.vcf.gz                   # VEP annotated
│   │   │   ├── FID001.rare.ensembl_vep_115.vcf.gz.tbi
│   │   │   ├── FID001.rare.ensembl_vep_115.annotated.bcf            # Fully annotated
│   │   │   └── FID001.rare.ensembl_vep_115.annotated.bcf.csi
│   │   ├── wombat/
│   │   │   ├── FID001.rare.ensembl_vep_115.annotated.tsv.gz                        # BCF to TSV
│   │   │   ├── FID001.rare.ensembl_vep_115.annotated.de_novo_mutations.tsv         # Wombat filtered
│   │   │   └── FID001.rare.ensembl_vep_115.annotated.rare_variants_high_impact.tsv
│   │   └── svs/
│   │       └── wisecondorx/
│   │           ├── FID001_aberrations.bed           # Family merged
│   │           └── FID001_aberrations.annotated.bed # Gene annotated
│   └── ...
├── cohorts/                       # Cohort-specific output directories
│   └── COHORT_NAME/
│       ├── vcfs/
│       │   ├── COHORT_NAME.common_gt.bcf            # Cohort common variants
│       │   └── COHORT_NAME.common_gt.bcf.csi
│       ├── wombat/
│       │   ├── COHORT_NAME.rare.ensembl_vep_115.annotated.de_novo_mutations.results.tsv
│       │   └── COHORT_NAME.rare.ensembl_vep_115.annotated.rare_variants_high_impact.results.tsv
│       └── svs/
│           └── wisecondorx/
│               └── COHORT_NAME_aberrations.bed      # Cohort merged aberrations
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

### Alignment Outputs

- **Final CRAM files**: `${data}/samples/${barcode}/sequences/${barcode}.${ref_name}.cram`
- **CRAM indices**: `${data}/samples/${barcode}/sequences/${barcode}.${ref_name}.cram.crai`
- **Bedgraph files**: `${data}/samples/${barcode}/sequences/${barcode}.${ref_name}.bedgraph.gz`
- **Bedgraph indices**: `${data}/samples/${barcode}/sequences/${barcode}.${ref_name}.bedgraph.gz.tbi`

### DeepVariant Sample Outputs

- **Individual VCF files**: `${data}/samples/${barcode}/deepvariant/${barcode}.vcf.gz`
- **Individual gVCF files**: `${data}/samples/${barcode}/deepvariant/${barcode}.g.vcf.gz`
- **VCF indices**: `*.vcf.gz.tbi` and `*.g.vcf.gz.tbi`
- **VAF bedgraph files**: `${data}/samples/${barcode}/sequences/${barcode}.vaf.bedgraph.gz`
- **VAF bedgraph indices**: `${data}/samples/${barcode}/sequences/${barcode}.vaf.bedgraph.gz.tbi`

### DeepVariant Family Outputs

- **Normalized BCF files**: `${data}/families/${FID}/vcfs/${FID}.norm.bcf`
- **Normalized BCF indices**: `${data}/families/${FID}/vcfs/${FID}.norm.bcf.csi`
- **Family pedigree files**: `${data}/families/${FID}/${FID}.pedigree.tsv`

### Annotation Outputs

- **Rare variant VCF files** (pre-VEP): `${data}/families/${FID}/vcfs/${FID}.rare.vcf.gz`
- **Common variant BCF files**: `${data}/families/${FID}/vcfs/${FID}.common.bcf`
- **Common variant BCF (GT only)**: `${data}/families/${FID}/vcfs/${FID}.common_gt.bcf`
- **VEP annotated rare VCF files**: `${data}/families/${FID}/vcfs/${FID}.rare.${vep_config_name}.vcf.gz`
- **Fully annotated BCF files**: `${data}/families/${FID}/vcfs/${FID}.rare.${vep_config_name}.annotated.bcf`
- **All indices**: `*.tbi` for VCF.gz and `*.csi` for BCF files

### Wombat Outputs

- **BCF to TSV conversion**: `${data}/families/${FID}/wombat/${FID}.rare.${vep_config_name}.annotated.tsv.gz`
- **PyWombat filtered results**: `${data}/families/${FID}/wombat/${FID}.rare.${vep_config_name}.annotated.${config_name}.tsv`
  - One file per configuration in `wombat_config_list`

### SNVs Cohort Outputs

- **Cohort common variants BCF**: `${data}/cohorts/${cohort_name}/vcfs/${cohort_name}.common_gt.bcf`
- **Cohort Wombat results**: `${data}/cohorts/${cohort_name}/wombat/${cohort_name}.rare.${vep_config_name}.annotated.${config_name}.results.tsv`
  - One file per configuration in `wombat_config_list`

### WisecondorX Outputs

- **Individual NPZ files**: `${data}/samples/${barcode}/svs/wisecondorx/${barcode}.npz`
- **Individual aberrations**: `${data}/samples/${barcode}/svs/wisecondorx/${barcode}_aberrations.bed`
- **Individual aberrations (chr format)**: `${data}/samples/${barcode}/svs/wisecondorx/${barcode}_aberrations.chr.bed`
- **Family merged aberrations**: `${data}/families/${FID}/svs/wisecondorx/${FID}_aberrations.bed`
- **Family annotated aberrations**: `${data}/families/${FID}/svs/wisecondorx/${FID}_aberrations.annotated.bed`
- **Cohort merged aberrations**: `${data}/cohorts/${cohort_name}/svs/wisecondorx/${cohort_name}_aberrations.bed`

### Extractor Outputs

- **Extracted variants**: `${data}/extractor/${original_filename}/${original_filename}.extracted.tsv`
- **Aggregated per family**: `${data}/extractor/${original_filename}/families/${FID}.extracted.tsv`
- **Fully aggregated**: `${data}/extractor/${original_filename}.aggregated.tsv`

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
ALIGNMENT: 2 individuals done and 1 to do
== SNVs/INDELs Calling ==
DEEPVARIANT_SAMPLE: 1 individuals done and 2 to do
DEEPVARIANT_FAMILY: 0 families done and 2 to do
ANNOTATION: 0 families done and 2 to do
WOMBAT: 0 families done and 2 to do
== Common Variants ==
SNVS_COHORT: common variants cohort bcf is needed: Yes
== SVs Calling ==
WISECONDORX PREDICT: 0 individuals done and 3 to do
== Other ==
EXTRACTOR: Skipped (no TSV files provided)
========================================================================================
```

This helps you understand what work will be performed.

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

Modify the SLURM settings in your parameters:

```yaml
slurm_account: "your_account"
slurm_partition: "ghfc"
```

Or in `nextflow.config`:

```groovy
process {
    clusterOptions = '-p ghfc --qos=ghfc --account=your_account'
}
```

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
