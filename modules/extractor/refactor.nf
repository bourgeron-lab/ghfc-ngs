process REFACTOR {
    tag "${original_filename}"
    
    publishDir "${params.data}", mode: 'copy', pattern: "families/*/*/*.tsv"
    publishDir "${params.data}", mode: 'copy', pattern: "samples/*/*/*.tsv"
    
    input:
    tuple val(original_filename), path(input_tsv), path(pedigree), path(chain_file)
    
    output:
    tuple val(original_filename), path("families/*/*/*.tsv"), path("samples/*/*/*.tsv"), emit: refactored_files
    tuple val(original_filename), path("family_list.txt"), emit: family_list
    tuple val(original_filename), path("sample_list.txt"), emit: sample_list
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import os
    from pyliftover import LiftOver
    
    # Column name mappings (priority order)
    CHR_COLS = ['chr', 'chrom', '#chrom', 'chromosome']
    POS_COLS_GRCH38 = ['pos_grch38', 'position_grch38', 'pos_hg38', 'position_hg38']
    POS_COLS_GRCH37 = ['pos_grch37', 'position_grch37', 'pos_hg19', 'position_hg19']
    POS_COLS_GENERIC = ['pos', 'position']
    REF_COLS = ['ref', 'reference']
    ALT_COLS = ['alt', 'alternative']
    SAMPLE_COLS = ['sample', 'sample_id', 'barcode']
    PATERNAL_COLS = ['paternal_id', 'father_id', 'father', 'paternal', 'pat']
    MATERNAL_COLS = ['maternal_id', 'mother_id', 'mother', 'maternal', 'mat']
    
    # Find first matching column from priority list (case-insensitive)
    def find_column(df, col_list):
        df_cols_lower = {c.lower(): c for c in df.columns}
        for col in col_list:
            if col.lower() in df_cols_lower:
                return df_cols_lower[col.lower()]
        return None
    
    # Detect position column and whether liftover is needed
    def detect_position_column(df):
        # Check GRCh38 columns first
        col = find_column(df, POS_COLS_GRCH38)
        if col:
            return col, False
        
        # Check GRCh37 columns
        col = find_column(df, POS_COLS_GRCH37)
        if col:
            return col, True
        
        # Check generic columns (assume GRCh38)
        col = find_column(df, POS_COLS_GENERIC)
        if col:
            return col, False
        
        return None, False
    
    # Lift position from GRCh37 to GRCh38
    def liftover_position(lo, chrom, pos):
        # Ensure chromosome has 'chr' prefix
        if not str(chrom).startswith('chr'):
            chrom = f'chr{chrom}'
        
        result = lo.convert_coordinate(chrom, int(pos) - 1)  # 0-based for liftover
        if result and len(result) > 0:
            return result[0][1] + 1  # Convert back to 1-based
        return None
    
    # Read pedigree to get sample-to-family mapping
    # Pedigree format: FID, IID, PAT, MAT, SEX, PHENO (no header or header starting with FID)
    pedigree_df = pd.read_csv('${pedigree}', sep='\\t', header=None, 
                               names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'])
    
    # Skip header if first row starts with 'FID'
    if pedigree_df.iloc[0]['FID'] == 'FID':
        pedigree_df = pedigree_df.iloc[1:].reset_index(drop=True)
    
    # Create sample-to-family mapping using positional columns
    sample_to_family = dict(zip(pedigree_df['IID'].astype(str), pedigree_df['FID'].astype(str)))
    
    # Read input TSV
    df = pd.read_csv('${input_tsv}', sep='\\t')
    
    # Find required columns
    chr_col = find_column(df, CHR_COLS)
    pos_col, needs_liftover = detect_position_column(df)
    ref_col = find_column(df, REF_COLS)
    alt_col = find_column(df, ALT_COLS)
    sample_col = find_column(df, SAMPLE_COLS)
    paternal_col = find_column(df, PATERNAL_COLS)
    maternal_col = find_column(df, MATERNAL_COLS)
    
    missing_cols = []
    if not chr_col: missing_cols.append('chromosome')
    if not pos_col: missing_cols.append('position')
    if not ref_col: missing_cols.append('ref')
    if not alt_col: missing_cols.append('alt')
    if not sample_col: missing_cols.append('sample')
    
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}. Available columns: {list(df.columns)}")
    
    print(f"Found paternal column: {paternal_col if paternal_col else 'None'}")
    print(f"Found maternal column: {maternal_col if maternal_col else 'None'}")
    
    # Initialize liftover if needed
    lo = None
    if needs_liftover:
        lo = LiftOver('${chain_file}')
    
    # Create standardized dataframe
    result_df = pd.DataFrame()
    result_df['chr'] = df[chr_col].apply(lambda x: str(x) if str(x).startswith('chr') else f'chr{x}')
    result_df['ref'] = df[ref_col]
    result_df['alt'] = df[alt_col]
    result_df['sample_id'] = df[sample_col].astype(str)
    
    # Add paternal and maternal IDs if available
    if paternal_col:
        result_df['paternal_id'] = df[paternal_col].astype(str).replace({'0': '', 'nan': '', '.': ''})
    else:
        result_df['paternal_id'] = ''
    
    if maternal_col:
        result_df['maternal_id'] = df[maternal_col].astype(str).replace({'0': '', 'nan': '', '.': ''})
    else:
        result_df['maternal_id'] = ''
    
    # Handle position (with liftover if needed)
    if needs_liftover:
        print(f"Lifting positions from GRCh37 to GRCh38 using {pos_col}")
        positions = []
        for _, row in df.iterrows():
            new_pos = liftover_position(lo, row[chr_col], row[pos_col])
            positions.append(new_pos)
        result_df['position'] = positions
        # Remove rows where liftover failed
        failed_count = result_df['position'].isna().sum()
        if failed_count > 0:
            print(f"Warning: {failed_count} positions failed liftover and will be excluded")
        result_df = result_df.dropna(subset=['position'])
        result_df['position'] = result_df['position'].astype(int)
    else:
        result_df['position'] = df[pos_col].astype(int)
    
    # Reorder columns
    result_df = result_df[['chr', 'position', 'ref', 'alt', 'sample_id', 'paternal_id', 'maternal_id']]
    
    # Map samples to families
    result_df['family_id'] = result_df['sample_id'].map(sample_to_family)
    
    # Check for samples not in pedigree
    unknown_samples = result_df[result_df['family_id'].isna()]['sample_id'].unique()
    if len(unknown_samples) > 0:
        print(f"Warning: {len(unknown_samples)} samples not found in pedigree: {list(unknown_samples)[:5]}...")
        result_df = result_df.dropna(subset=['family_id'])
    
    # Create output directories and files
    families_processed = set()
    samples_processed = set()
    
    original_filename = '${original_filename}'
    
    for family_id in result_df['family_id'].unique():
        family_df = result_df[result_df['family_id'] == family_id].copy()
        family_df = family_df[['chr', 'position', 'ref', 'alt', 'sample_id', 'paternal_id', 'maternal_id']]
        
        # Create family output
        family_dir = f"families/{family_id}/extractor"
        os.makedirs(family_dir, exist_ok=True)
        family_output = f"{family_dir}/{family_id}.{original_filename}.tsv"
        family_df.to_csv(family_output, sep='\\t', index=False)
        families_processed.add(family_id)
        
        # Create sample outputs
        for sample_id in family_df['sample_id'].unique():
            sample_df = family_df[family_df['sample_id'] == sample_id].copy()
            sample_dir = f"samples/{sample_id}/extractor"
            os.makedirs(sample_dir, exist_ok=True)
            sample_output = f"{sample_dir}/{sample_id}.{original_filename}.tsv"
            sample_df.to_csv(sample_output, sep='\\t', index=False)
            samples_processed.add(sample_id)
    
    # Write lists for downstream processing
    with open('family_list.txt', 'w') as f:
        for fid in sorted(families_processed):
            f.write(f"{fid}\\n")
    
    with open('sample_list.txt', 'w') as f:
        for sid in sorted(samples_processed):
            f.write(f"{sid}\\n")
    
    print(f"Processed {len(result_df)} variants into {len(families_processed)} families and {len(samples_processed)} samples")
    """
}
