process FAM_AGGREGATE {
    tag "${fid}_${original_filename}"
    
    publishDir "${params.data}/families/${fid}/extractor", mode: 'copy'
    
    input:
    tuple val(fid), val(original_filename), path(extractor_tsv), path(fam_tsv_found), path(fam_tsv_notfound), path(fam_bcf_found), path(fam_bcf_notfound), path(ind_gvcf_files)
    
    output:
    tuple val(fid), val(original_filename), path("${fid}.${original_filename}.summary.tsv")
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import os
    import glob
    
    fid = '${fid}'
    original_filename = '${original_filename}'
    
    # Read extractor TSV (original family variants)
    extractor_df = pd.read_csv('${extractor_tsv}', sep='\\t')
    
    # Create unique key for matching
    extractor_df['_match_key'] = (
        extractor_df['chr'].astype(str) + '_' + 
        extractor_df['position'].astype(str) + '_' + 
        extractor_df['ref'].astype(str) + '_' + 
        extractor_df['alt'].astype(str) + '_' +
        extractor_df['sample_id'].astype(str)
    )
    
    # Read fam_tsv results
    fam_tsv_found_df = pd.read_csv('${fam_tsv_found}', sep='\\t')
    fam_tsv_notfound_df = pd.read_csv('${fam_tsv_notfound}', sep='\\t')
    
    # Create keys for fam_tsv found
    fam_tsv_found_keys = set()
    if len(fam_tsv_found_df) > 0 and 'chr' in fam_tsv_found_df.columns:
        fam_tsv_found_df['_match_key'] = (
            fam_tsv_found_df['chr'].astype(str) + '_' + 
            fam_tsv_found_df['position'].astype(str) + '_' + 
            fam_tsv_found_df['ref'].astype(str) + '_' + 
            fam_tsv_found_df['alt'].astype(str) + '_' +
            fam_tsv_found_df['sample_id'].astype(str)
        )
        fam_tsv_found_keys = set(fam_tsv_found_df['_match_key'])
    
    # Read fam_bcf results
    fam_bcf_found_df = pd.read_csv('${fam_bcf_found}', sep='\\t')
    fam_bcf_notfound_df = pd.read_csv('${fam_bcf_notfound}', sep='\\t')
    
    # Create keys for fam_bcf found
    fam_bcf_found_keys = set()
    if len(fam_bcf_found_df) > 0 and 'chr' in fam_bcf_found_df.columns:
        fam_bcf_found_df['_match_key'] = (
            fam_bcf_found_df['chr'].astype(str) + '_' + 
            fam_bcf_found_df['position'].astype(str) + '_' + 
            fam_bcf_found_df['ref'].astype(str) + '_' + 
            fam_bcf_found_df['alt'].astype(str) + '_' +
            fam_bcf_found_df['sample_id'].astype(str)
        )
        fam_bcf_found_keys = set(fam_bcf_found_df['_match_key'])
    
    # Read and combine ind_gvcf results
    gvcf_data = {}
    for gvcf_file in glob.glob('*.ind_gvcf.tsv'):
        gvcf_df = pd.read_csv(gvcf_file, sep='\\t')
        for _, row in gvcf_df.iterrows():
            key = f"{row['chr']}_{row['position']}_{row['ref']}_{row['alt']}_{row['sample_id']}"
            gvcf_data[key] = row.get('gvcf_FILTER', None)
    
    # Build summary
    results = []
    for _, row in extractor_df.iterrows():
        key = row['_match_key']
        
        result_row = row.drop('_match_key').to_dict()
        result_row['extractor_gvcf_filter'] = gvcf_data.get(key, None)
        result_row['extractor_fam_bcf_found'] = key in fam_bcf_found_keys
        result_row['extractor_fam_wombat_found'] = key in fam_tsv_found_keys
        
        results.append(result_row)
    
    # Create output dataframe
    summary_df = pd.DataFrame(results)
    
    # Ensure column order
    base_cols = ['chr', 'position', 'ref', 'alt', 'sample_id']
    extra_cols = ['extractor_gvcf_filter', 'extractor_fam_bcf_found', 'extractor_fam_wombat_found']
    other_cols = [c for c in summary_df.columns if c not in base_cols + extra_cols]
    summary_df = summary_df[base_cols + extra_cols + other_cols]
    
    # Save output
    output_file = f"{fid}.{original_filename}.summary.tsv"
    summary_df.to_csv(output_file, sep='\\t', index=False)
    
    found_bcf = summary_df['extractor_fam_bcf_found'].sum()
    found_wombat = summary_df['extractor_fam_wombat_found'].sum()
    found_gvcf = summary_df['extractor_gvcf_filter'].notna().sum()
    
    print(f"Summary for {fid}: {len(summary_df)} variants")
    print(f"  - Found in BCF: {found_bcf}")
    print(f"  - Found in Wombat: {found_wombat}")
    print(f"  - Found in gVCF: {found_gvcf}")
    """
}
