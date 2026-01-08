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
    # Wombat columns: CHROM, POS, REF, ALT, ...
    
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
    
    # Create found dataframe (merge on match key - may produce multiple rows)
    found_extractor = extractor_df[extractor_df['_match_key'].isin(matched_keys)].copy()
    found_wombat = wombat_df[wombat_df['_match_key'].isin(matched_keys)].copy()
    
    if len(found_extractor) > 0:
        # Merge to get wombat annotations
        found_df = pd.merge(
            found_extractor,
            found_wombat,
            on='_match_key',
            how='left',
            suffixes=('', '_wombat')
        )
        # Remove temporary key column
        found_df = found_df.drop(columns=['_match_key'])
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
        pd.DataFrame(columns=extractor_df.drop(columns=['_match_key']).columns.tolist() + 
                     wombat_df.drop(columns=['_match_key']).columns.tolist()).to_csv(found_output, sep='\\t', index=False)
    
    notfound_df.to_csv(notfound_output, sep='\\t', index=False)
    
    print(f"Found {len(found_df)} variant matches, {len(notfound_df)} not found")
    """
}
