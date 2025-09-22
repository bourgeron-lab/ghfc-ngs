# GHFC WGS Family-based Variant Calling Pipeline (Nextflow)

This is a Nextflow implementation of the GHFC WGS family-based variant calling pipeline. The pipeline supports alignment from FASTQ files, realignment from existing CRAM files, individual variant calling with DeepVariant, and family-based joint calling with GLnexus.

## Features

- **Family-aware pipeline** processing based on pedigree structure
- **Smart dependency resolution** - automatically determines what needs to be run
- **BWA-MEM2 alignment** from paired-end FASTQ files
- **Realignment** from existing CRAM files using Bazam
- **DeepVariant variant calling** with 3-stage processing pipeline
- **GLnexus family calling** for joint variant calling per family
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

## Quick Start

### 1. Clone/Download the Pipeline

```bash
cd /path/to/your/workspace
# Pipeline files should be in the current directory
```

### 2. Prepare Your Pedigree File

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

# Run alignment and variant calling
./run_pipeline.sh --profile slurm,apptainer --steps "alignment,deepvariant" --params-file my_params.yml

# Run only variant calling (assumes CRAM files exist)
./run_pipeline.sh --profile slurm,apptainer --steps "deepvariant" --params-file my_params.yml

# Run only family calling (assumes gVCF files exist)
./run_pipeline.sh --profile slurm,apptainer --steps "family_calling" --params-file my_params.yml

# Full pipeline
./run_pipeline.sh --profile slurm,apptainer --steps "alignment,deepvariant,family_calling" --params-file my_params.yml
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

The pipeline supports three main steps that must be explicitly listed in the `steps` parameter:

**Available Steps:**

- **`alignment`**: BWA-MEM2 alignment from FASTQ files and realignment from CRAM files
- **`deepvariant`**: Individual variant calling using Google DeepVariant (3-stage process)
- **`family_calling`**: Family-based joint calling using GLnexus

**Step Dependencies:**

- `family_calling` requires gVCF files (triggers `deepvariant` if missing)
- `deepvariant` requires CRAM files (triggers `alignment` if missing)
- The pipeline will error if required steps are not listed in parameters

### Parameters File (params.yml)

Key parameters to configure:

```yaml
# Pipeline steps - ALL REQUIRED STEPS MUST BE LISTED
steps: ["alignment", "deepvariant", "family_calling", "vep_annotation"]

# Data directories
data: "/path/to/your/data/"
scratch: "/path/to/your/scratch/"

# Pedigree file (optional - defaults to ${data}/pedigree.tsv)
pedigree: "/path/to/your/pedigree.tsv"

# Reference genomes
ref: "/path/to/reference/genome.fa"
ref_name: "GRCh38_GIABv3"  # Used in output file naming
oldref: "/path/to/old/reference/genome.fa"  # For realignment

# Tool paths (when not using containers)
bwa_mem2: "bwa-mem2"
samtools: "samtools"
sambamba: "sambamba"
samblaster: "samblaster"
bazam: "/path/to/bazam.jar"

# GLnexus configuration
glnexus_config: "DeepVariant_unfiltered"

# VEP annotation configuration
vep_config: "/path/to/vep/config.ini"      # Path to VEP config INI file
vep_config_name: "ensembl_vep_115"         # Name suffix for output files
```

### Smart File Detection

The pipeline automatically detects existing files and skips unnecessary work:

- **Family VCF files**: `${data}/families/${FID}/vcfs/${FID}.vcf.gz`
- **VEP annotated VCF files**: `${data}/families/${FID}/vcfs/${FID}.${vep_config_name}.vcf.gz`
- **Individual gVCF files**: `${data}/samples/${barcode}/deepvariant/${barcode}.g.vcf.gz`
- **CRAM files**: `${data}/samples/${barcode}/sequences/${barcode}.${ref_name}.cram`

If these files exist with their indices, the corresponding steps are skipped.

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
└── samples/                       # Sample-specific output directories
    ├── BC001/                     # Sample barcode directory
    │   ├── sequences/             # CRAM and bedgraph files
    │   │   ├── BC001.hg38.cram
    │   │   ├── BC001.hg38.cram.crai
    │   │   ├── BC001.hg38.bedgraph.gz
    │   │   └── BC001.hg38.bedgraph.gz.tbi
    │   └── deepvariant/           # DeepVariant outputs
    │       ├── BC001.g.vcf.gz
    │       ├── BC001.g.vcf.gz.tbi
    │       ├── BC001.vcf.gz
    │       └── BC001.vcf.gz.tbi
    └── ...
└── families/                      # Family-specific output directories
    └── FID001/                    # Family directory
        └── vcfs/
            ├── FID001.vcf.gz
            ├── FID001.vcf.gz.tbi
            ├── FID001.ensembl_vep_115.vcf.gz
            └── FID001.ensembl_vep_115.vcf.gz.tbi
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

### DeepVariant Outputs

- **Individual VCF files**: `${data}/samples/${barcode}/deepvariant/${barcode}.vcf.gz`
- **Individual gVCF files**: `${data}/samples/${barcode}/deepvariant/${barcode}.g.vcf.gz`
- **VCF indices**: `*.vcf.gz.tbi` and `*.g.vcf.gz.tbi`

### Family Calling Outputs

- **Family VCF files**: `${data}/families/${FID}/vcfs/${FID}.vcf.gz`
- **Family VCF indices**: `${data}/families/${FID}/vcfs/${FID}.vcf.gz.tbi`

### VEP Annotation Outputs

- **VEP annotated VCF files**: `${data}/families/${FID}/vcfs/${FID}.${vep_config_name}.vcf.gz`
- **VEP annotated VCF indices**: `${data}/families/${FID}/vcfs/${FID}.${vep_config_name}.vcf.gz.tbi`

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
ALIGNMENT:
    - Existing CRAM files: 2 individuals
    - Need alignment: 1 individuals

DEEPVARIANT:
    - Existing gVCF files: 1 individuals  
    - Need variant calling: 2 individuals

FAMILY_CALLING:
    - Existing family VCFs: 0 families
    - Need family calling: 2 families
========================================================================================
```

This helps you understand what work will be performed.

### Error Checking

The pipeline validates dependencies and will stop with clear error messages:

```bash
ERROR: DeepVariant step is required for 3 individuals but not included in steps parameter

Please add the required steps to your parameters or ensure all required files exist.
Available steps: alignment, deepvariant, family_calling
```

### Custom Configuration

You can override any configuration parameter:

```bash
nextflow run main.nf \
    -profile slurm,apptainer \
    --data /my/data \
    --scratch /my/scratch \
    --ref /my/reference.fa \
    --ref-name "MyRef" \
    --steps "alignment,deepvariant,family_calling"
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
2. **Step dependency errors**: Add required steps to the `steps` parameter
3. **File not found errors**: Check file paths and ensure barcode names match between pedigree and data files
4. **Out of memory errors**: Increase memory allocation for memory-intensive processes
5. **GLnexus failures**: Often memory-related, try increasing memory allocation

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
2. **Update parameters**: Remove `skip_*` options, add all required steps to `steps` list
3. **Update paths**: Family VCFs now in `vcf/families/` directory
4. **Review naming**: CRAM files now include `ref_name` in filename

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

## Citation

If you use this pipeline in your research, please cite:

- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319.
- **BWA-MEM2**: Vasimuddin, M., et al. (2019). Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. IEEE IPDPS.
- **DeepVariant**: Poplin, R., et al. (2018). A universal SNP and small-indel variant caller using deep neural networks. Nature Biotechnology, 36(10), 983-987.
- **GLnexus**: Yun, T., et al. (2021). Accurate, scalable cohort variant calls using DeepVariant and GLnexus. Bioinformatics, 37(5), 682-685.
