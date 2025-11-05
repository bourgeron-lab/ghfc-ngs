#!/usr/bin/env nextflow

/*
========================================================================================
    GHFC WGS Legacy Format Migration - Entry Point
========================================================================================
    Run this workflow to migrate legacy VCF.gz files to BCF format before running
    the main pipeline.
    
    Usage:
        nextflow run migrate.nf -params-file params.yml
        
    Or with specific data directory:
        nextflow run migrate.nf --data /path/to/data
    
    The migration workflow automatically uses the slurm,apptainer profiles.
    To override:
        nextflow run migrate.nf -params-file params.yml -profile local
========================================================================================
*/

nextflow.enable.dsl = 2

// Include the migration workflow
include { MIGRATE_LEGACY_FORMATS } from './workflows/migrate_legacy_formats'

// Default parameters
params.data = ""

// Run the migration workflow
workflow {
    MIGRATE_LEGACY_FORMATS()
}
