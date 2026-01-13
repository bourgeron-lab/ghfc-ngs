process PROCESS_FAM_TSV {
    tag "${fid}_${original_filename}"
    
    publishDir "${params.data}/families/${fid}/extractor", mode: 'copy'
    
    input:
    tuple val(fid), val(original_filename), path(extractor_tsv), path(wombat_tsv)
    
    output:
    tuple val(fid), val(original_filename), path("${fid}.${original_filename}.fam_tsv.found.tsv"), path("${fid}.${original_filename}.fam_tsv.notfound.tsv")
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import gzip
    
    fid = '${fid}'
    original_filename = '${original_filename}'
    
    # Read extractor TSV
    extractor_df = pd.read_csv('${extractor_tsv}', sep='\\t')
    
    # Read wombat TSV (gzipped bcf2tsv output)
    with gzip.open('${wombat_tsv}', 'rt') as f:
        wombat_df = pd.read_csv(f, sep='\\t')
    
    # Standardize column names for matching
    # Extractor columns: chr, position, ref, alt, sample_id
    # Wombat columns: CHROM or #CHROM, POS, REF, ALT, ...
    
    # Handle both CHROM and #CHROM column names
    if '#CHROM' in wombat_df.columns:
        wombat_df = wombat_df.rename(columns={'#CHROM': 'CHROM'})
    
    # Create match keys
    extractor_df['_match_key'] = (
        extractor_df['chr'].astype(str) + '_' + 
        extractor_df['position'].astype(str) + '_' + 
        extractor_df['ref'].astype(str) + '_' + 
        extractor_df['alt'].astype(str)
    )
    
    wombat_df['_match_key'] = (
        wombat_df['CHROM'].astype(str) + '_' + 
        wombat_df['POS'].astype(str) + '_' + 
        wombat_df['REF'].astype(str) + '_' + 
        wombat_df['ALT'].astype(str)
    )
    
    # Find matches
    matched_keys = set(extractor_df['_match_key']) & set(wombat_df['_match_key'])
    
    # Identify non-sample columns in wombat (those that don't contain ':')
    wombat_non_sample_cols = [col for col in wombat_df.columns if ':' not in col and col != '_match_key']
    
    # Create found dataframe (merge on match key - may produce multiple rows)
    found_extractor = extractor_df[extractor_df['_match_key'].isin(matched_keys)].copy()
    found_wombat = wombat_df[wombat_df['_match_key'].isin(matched_keys)].copy()
    
    if len(found_extractor) > 0:
        # For each row, extract the sample-specific columns from wombat
        found_rows = []
        for _, ext_row in found_extractor.iterrows():
            sample_id = ext_row['sample_id']
            match_key = ext_row['_match_key']
            
            # Find the matching wombat row
            wombat_row = found_wombat[found_wombat['_match_key'] == match_key]
            
            if len(wombat_row) > 0:
                wombat_row = wombat_row.iloc[0]
                
                # Find the sample column (format: SAMPLE:GT:SAMPLE:DP:SAMPLE:GQ:SAMPLE:AD)
                sample_col = None
                for col in found_wombat.columns:
                    if col.startswith(f"{sample_id}:"):
                        sample_col = col
                        break
                
                # Start with extractor columns
                row_dict = ext_row.drop('_match_key').to_dict()
                
                # Add non-sample wombat columns
                for col in wombat_non_sample_cols:
                    row_dict[f"wombat_{col}"] = wombat_row[col]
                
                # Parse sample-specific columns
                if sample_col and pd.notna(wombat_row[sample_col]):
                    sample_values = str(wombat_row[sample_col]).split(':')
                    # Expected format: GT:DP:GQ:AD (4 values)
                    row_dict['GT'] = sample_values[0] if len(sample_values) > 0 else '.'
                    row_dict['DP'] = sample_values[1] if len(sample_values) > 1 else '.'
                    row_dict['GQ'] = sample_values[2] if len(sample_values) > 2 else '.'
                    row_dict['AD'] = sample_values[3] if len(sample_values) > 3 else '.'
                else:
                    # No sample data found
                    row_dict['GT'] = '.'
                    row_dict['DP'] = '.'
                    row_dict['GQ'] = '.'
                    row_dict['AD'] = '.'
                
                found_rows.append(row_dict)
        
        found_df = pd.DataFrame(found_rows) if found_rows else pd.DataFrame()
    else:
        found_df = pd.DataFrame()
    
    # Create not found dataframe
    notfound_df = extractor_df[~extractor_df['_match_key'].isin(matched_keys)].copy()
    notfound_df = notfound_df.drop(columns=['_match_key'])
    
    # Save outputs
    found_output = f"{fid}.{original_filename}.fam_tsv.found.tsv"
    notfound_output = f"{fid}.{original_filename}.fam_tsv.notfound.tsv"
    
    if len(found_df) > 0:
        found_df.to_csv(found_output, sep='\\t', index=False)
    else:
        # Create empty file with header
        header_cols = ['chr', 'position', 'ref', 'alt', 'sample_id'] + \
                      [f"wombat_{col}" for col in wombat_non_sample_cols] + \
                      ['GT', 'DP', 'GQ', 'AD']
        pd.DataFrame(columns=header_cols).to_csv(found_output, sep='\\t', index=False)
    
    notfound_df.to_csv(notfound_output, sep='\\t', index=False)
    
    print(f"Found {len(found_df)} variant matches, {len(notfound_df)} not found")
    """
}
