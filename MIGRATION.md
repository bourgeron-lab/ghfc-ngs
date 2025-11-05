# Legacy Format Migration Guide

## Overview

This guide explains how to migrate legacy VCF.gz files to the new BCF format used by the updated GHFC WGS pipeline.

## What This Migration Does

The migration workflow performs the following operations:

### 1. **File Conversions**

Converts legacy VCF.gz files to BCF format:

- `families/{FID}/vcfs/{FID}.norm.vcf.gz` → `families/{FID}/vcfs/{FID}.norm.bcf`
- `families/{FID}/vcfs/{FID}.common.vcf.gz` → `families/{FID}/vcfs/{FID}.common.bcf`

For each conversion:

- Creates new BCF files with proper `.csi` indexes
- Removes old VCF.gz files and their `.tbi` indexes

### 2. **Cleanup of Obsolete Intermediate Files**

Removes files that are no longer published by the pipeline:

- `families/{FID}/vcfs/{FID}.vcf.gz` (GLnexus intermediate output)
- `families/{FID}/vcfs/{FID}.gnomad.vcf.gz` (gnomAD annotation intermediate)
- `families/{FID}/vcfs/{FID}.gnomad.bcf` (gnomAD annotation intermediate)

These files are now kept only in the Nextflow work directory and not published to output folders.

## When to Run Migration

Run this migration workflow **once, before running the main pipeline** if:

1. You have existing data from an older version of the pipeline
2. You have `{FID}.norm.vcf.gz` or `{FID}.common.vcf.gz` files
3. You want to clean up obsolete intermediate files

## How to Run Migration

### Basic Usage

The migration workflow automatically uses SLURM for job execution and Apptainer for containers (profiles: `slurm,apptainer`).

```bash
# Using your existing params.yml file (default: uses SLURM + Apptainer)
nextflow run migrate.nf -params-file params.yml

# Or specify data directory directly
nextflow run migrate.nf --data /path/to/your/data

# Override default profiles (e.g., for local execution)
nextflow run migrate.nf -params-file params.yml -profile local

# Using the run_pipeline.sh script
./run_pipeline.sh --migrate --params-file params.yml

# With custom Nextflow configuration
nextflow run migrate.nf -params-file params.yml -c your_config.config
```

### Resource Requirements

The migration processes use the following resources on SLURM:

- **VCF to BCF Conversion**: 2 CPUs, 8 GB memory, 4 hours
- **Cleanup**: 1 CPU, 2 GB memory, 1 hour
- **Container**: `staphb/bcftools:latest` (via Apptainer)

These settings are defined in `nextflow.config` and can be overridden if needed.

### What Happens During Migration

1. **Discovery Phase**: The workflow scans your data directory for:
   - Legacy `{FID}.norm.vcf.gz` files that need conversion
   - Legacy `{FID}.common.vcf.gz` files that need conversion
   - All family directories for cleanup

2. **Conversion Phase**:
   - Converts each VCF.gz file to BCF format using `bcftools view -O b`
   - Indexes each new BCF file with `bcftools index`

3. **Cleanup Phase**:
   - Removes original VCF.gz files and their indexes after successful conversion
   - Removes obsolete intermediate files

4. **Completion**: You'll see a summary of:
   - Number of files converted
   - Number of families cleaned up
   - Total duration

### Example Output

```text
========================================================================================
                    GHFC WGS Legacy Format Migration
========================================================================================
Data directory   : /pasteur/helix/projects/ghfc_wgs/WGS/test-GRCh38/
========================================================================================
This workflow will:
1. Convert {FID}.norm.vcf.gz → {FID}.norm.bcf
2. Convert {FID}.common.vcf.gz → {FID}.common.bcf
3. Remove obsolete intermediate files
========================================================================================

Found 15 families to check for migration

  Found legacy norm.vcf.gz for family FID001
  Found legacy common.vcf.gz for family FID001
  Found legacy norm.vcf.gz for family FID002
  ...

Migration plan:
  - 10 norm.vcf.gz files to convert
  - 8 common.vcf.gz files to convert
  - 15 families to clean up

[work directory] Converting FID001.norm.vcf.gz to BCF format
[work directory] Converting FID001.common.vcf.gz to BCF format
...

========================================================================================
Migration workflow completed!
========================================================================================
Status: SUCCESS
Duration: 5m 23s

✓ Legacy VCF.gz files have been converted to BCF format
✓ Obsolete intermediate files have been removed

You can now run the main pipeline safely.
========================================================================================
```

## Safety Features

- **No Overwriting**: The migration will not overwrite existing BCF files. If `{FID}.norm.bcf` already exists, it will skip that family.
- **Separate Process**: Migration runs independently of the main pipeline
- **Logged Actions**: All conversions and deletions are logged

## After Migration

Once migration is complete:

1. **Verify Converted Files**: Check that BCF files were created correctly

   ```bash
   ls -lh ${data}/families/*/vcfs/*.bcf
   ```

2. **Run Main Pipeline**: You can now run the main pipeline normally

   ```bash
   nextflow run main.nf -params-file params.yml
   ```

3. **Optional - Clean Work Directory**: Remove migration work files

   ```bash
   rm -rf work/
   ```

## Troubleshooting

### Migration Reports "Nothing to migrate"

- Check that your `--data` parameter points to the correct directory
- Verify that family directories exist in `${data}/families/`
- Check if BCF files already exist (migration skips existing files)

### Conversion Fails

- Ensure bcftools container is accessible
- Check that VCF.gz files are valid and indexed
- Verify sufficient disk space for new BCF files

### Need to Re-run Migration

If you need to re-run migration:

1. Remove the BCF files that were created: `rm ${data}/families/*/vcfs/*.norm.bcf*`
2. Run the migration workflow again

## Technical Details

### File Formats

**VCF.gz Format (Legacy)**:

- Compressed VCF file (gzip)
- Index: `.tbi` (tabix)
- Larger file size
- Slower access

**BCF Format (Current)**:

- Binary compressed VCF
- Index: `.csi` (CSI)
- Smaller file size (~30% smaller)
- Faster access and processing
- Native bcftools format

### Module Structure

The migration workflow consists of:

- `migrate.nf`: Entry point
- `workflows/migrate_legacy_formats.nf`: Main workflow logic
- `modules/migration/convert_norm_vcf_to_bcf.nf`: Normalize file conversion
- `modules/migration/convert_common_vcf_to_bcf.nf`: Common file conversion
- `modules/migration/cleanup_legacy_files.nf`: Legacy file cleanup

## Questions?

If you encounter issues or have questions about migration, please:

1. Check the Nextflow log files
2. Review this documentation
3. Contact the pipeline maintainers
