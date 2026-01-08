process AGGREGATE {
    tag "${original_filename}"
    
    publishDir "${params.data}/cohorts/${params.cohort_name}/extractor", mode: 'copy'
    
    input:
    tuple val(original_filename), path(summary_files)
    
    output:
    tuple val(original_filename), path("${original_filename}.summary.tsv")
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import os
    import glob
    
    original_filename = '${original_filename}'
    
    # Read all summary files
    all_summaries = []
    for summary_file in glob.glob('*.summary.tsv'):
        df = pd.read_csv(summary_file, sep='\\t')
        # Extract FID from filename
        fid = summary_file.split('.')[0]
        df['family_id'] = fid
        all_summaries.append(df)
    
    if all_summaries:
        # Concatenate all summaries
        combined_df = pd.concat(all_summaries, ignore_index=True)
        
        # Reorder columns to put family_id after sample_id
        base_cols = ['chr', 'position', 'ref', 'alt', 'sample_id', 'family_id']
        extra_cols = ['extractor_gvcf_filter', 'extractor_fam_bcf_found', 'extractor_fam_wombat_found']
        other_cols = [c for c in combined_df.columns if c not in base_cols + extra_cols]
        
        available_cols = [c for c in base_cols + extra_cols + other_cols if c in combined_df.columns]
        combined_df = combined_df[available_cols]
    else:
        combined_df = pd.DataFrame()
    
    # Save output
    output_file = f"{original_filename}.summary.tsv"
    combined_df.to_csv(output_file, sep='\\t', index=False)
    
    print(f"Aggregated {len(all_summaries)} family summaries into cohort summary with {len(combined_df)} total variants")
    """
}
