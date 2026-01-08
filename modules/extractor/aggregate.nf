process AGGREGATE {
    tag "${original_filename}"
    
    publishDir "${params.data}/cohorts/${params.cohort_name}/extractor", mode: 'copy'
    
    input:
    tuple val(original_filename), path(summary_files), path(fam_tsv_found_files), path(fam_tsv_notfound_files)
    
    output:
    tuple val(original_filename), path("${original_filename}.summary.tsv"), path("${original_filename}.fam_tsv.found.tsv"), path("${original_filename}.fam_tsv.notfound.tsv")
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import os
    import glob
    
    original_filename = '${original_filename}'
    
    def aggregate_files(pattern, output_name, file_type='summary'):
        all_dfs = []
        for f in glob.glob(pattern):
            try:
                df = pd.read_csv(f, sep='\\t')
                if len(df) > 0:
                    # Extract FID from filename (e.g., FID.original_filename.fam_tsv.found.tsv)
                    fid = os.path.basename(f).split('.')[0]
                    df['family_id'] = fid
                    all_dfs.append(df)
            except Exception as e:
                print(f"Warning: Could not read {f}: {e}")
        
        if all_dfs:
            combined_df = pd.concat(all_dfs, ignore_index=True)
            # Reorder columns to put family_id early
            if 'sample_id' in combined_df.columns:
                base_cols = ['chr', 'position', 'ref', 'alt', 'sample_id', 'family_id']
            else:
                base_cols = ['chr', 'position', 'ref', 'alt', 'family_id']
            other_cols = [c for c in combined_df.columns if c not in base_cols]
            available_cols = [c for c in base_cols if c in combined_df.columns] + other_cols
            combined_df = combined_df[available_cols]
        else:
            combined_df = pd.DataFrame()
        
        combined_df.to_csv(output_name, sep='\\t', index=False)
        return len(all_dfs), len(combined_df)
    
    # Aggregate summary files
    n_summaries, n_summary_rows = aggregate_files('*.summary.tsv', f"{original_filename}.summary.tsv", 'summary')
    print(f"Aggregated {n_summaries} family summaries into cohort summary with {n_summary_rows} total variants")
    
    # Aggregate fam_tsv found files
    n_found, n_found_rows = aggregate_files('*.fam_tsv.found.tsv', f"{original_filename}.fam_tsv.found.tsv", 'found')
    print(f"Aggregated {n_found} fam_tsv.found files with {n_found_rows} total rows")
    
    # Aggregate fam_tsv notfound files
    n_notfound, n_notfound_rows = aggregate_files('*.fam_tsv.notfound.tsv', f"{original_filename}.fam_tsv.notfound.tsv", 'notfound')
    print(f"Aggregated {n_notfound} fam_tsv.notfound files with {n_notfound_rows} total rows")
    """
}
