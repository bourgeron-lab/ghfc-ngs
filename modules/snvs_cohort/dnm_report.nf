process DNM_REPORT {
    tag "${cohort_name}"
    
    publishDir "${params.data}/cohorts/${cohort_name}/reports", mode: 'copy'
    
    input:
    tuple val(cohort_name), val(vep_config_name), path(dnm_tsv)
    
    output:
    path "${cohort_name}.${vep_config_name}.dnm.report.pdf"
    
    script:
    def report_script = "${projectDir}/modules/snvs_cohort/scripts/dnm_report"
    def output_pdf = "${cohort_name}.${vep_config_name}.dnm.report.pdf"
    """
    ${report_script} ${dnm_tsv} ${output_pdf}
    """
}
