process PROCESS_FAM_BCF {
    tag "${fid}_${original_filename}"
    
    publishDir "${params.data}/families/${fid}/extractor", mode: 'copy'
    
    input:
    tuple val(fid), val(original_filename), path(extractor_tsv), path(norm_bcf), path(norm_bcf_csi)
    
    output:
    tuple val(fid), val(original_filename), path("${fid}.${original_filename}.fam_bcf.found.tsv"), path("${fid}.${original_filename}.fam_bcf.notfound.tsv")
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import subprocess
    import tempfile
    import os
    
    fid = '${fid}'
    original_filename = '${original_filename}'
    
    # Read extractor TSV
    extractor_df = pd.read_csv('${extractor_tsv}', sep='\\t')
    
    print(f"Read {len(extractor_df)} variants from extractor TSV")
    print(f"Sample chromosomes in extractor: {extractor_df['chr'].unique()[:5].tolist()}")
    
    # Create regions file for bcftools
    regions = []
    for _, row in extractor_df.iterrows():
        regions.append(f"{row['chr']}:{row['position']}-{row['position']}")
    
    # Write regions to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        regions_file = f.name
        for r in regions:
            f.write(r + '\\n')
    
    # Extract variants from BCF using bcftools
    print(f"Running bcftools query with {len(regions)} regions")
    try:
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        if result.stdout:
            print(f"bcftools stdout: {result.stdout[:200]}")
        if result.stderr:
            print(f"bcftools stderr: {result.stderr[:200]}"FO\\t%FORMAT[\\t%SAMPLE=%TGT]\\n' '${norm_bcf}' > {bcf_output}"
    
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"bcftools error: {e.stderr}")
        # Create empty output
        pd.DataFrame().to_csv(bcf_output, sep='\\t', index=False)
    
    os.unlink(regions_file)
    lines = f.readlines()
            print(f"BCF output contains {len(lines)} lines")
            if lines:
                print(f"First line: {lines[0][:200]}")
            for line in lines:
                parts = line.strip().split('\\t')
                if len(parts) >= 4:
                    chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
                    info = parts[4] if len(parts) > 4 else ''
                    fmt = parts[5] if len(parts) > 5 else ''
                    samples = parts[6:] if len(parts) > 6 else []
                    
                    # Handle multiallelic: ALT can be comma-separated
                    for single_alt in alt.split(','):
                        key = f"{chrom}_{pos}_{ref}_{single_alt}"
                        bcf_variants[key] = {
                            'INFO': info,
                            'FORMAT': fmt,
                            'SAMPLES': '\\t'.join(samples)
                        }
            print(f"Parsed {len(bcf_variants)} variant keys from BCF")
            if bcf_variants:
                print(f"Sample keys: {list(bcf_variants.keys())[:3]}")
    else:
        print(f"BCF output file is empty or doesn't exist")ey = f"{chrom}_{pos}_{ref}_{alt}"
                    bcf_variants[key] = {
                        'INFO': info,
                        'FORMAT': fmt,
                        'SAMPLES': '\\t'.join(samples)
                    }
    
    # Create match keys for extractor
    extractor_df['_match_key'] = (
        extractor_df['chr'].astype(str) + '_' + 
        extractor_df['position'].astype(str) + '_' + 
        extractor_df['ref'].astype(str) + '_' + 
        extractor_df['alt'].astype(str)
    )
    
    print(f"Sample extractor keys: {extractor_df['_match_key'].head(3).tolist()}")
    
    # Split into found and not found
    found_rows = []
    notfound_rows = []
    
    for _, row in extractor_df.iterrows():
        key = row['_match_key']
        if key in bcf_variants:
            row_dict = row.drop('_match_key').to_dict()
            row_dict['INFO'] = bcf_variants[key]['INFO']
            row_dict['FORMAT'] = bcf_variants[key]['FORMAT']
            # Add sample columns
            for sample_info in bcf_variants[key]['SAMPLES'].split('\\t'):
                if '=' in sample_info:
                    sample_name, sample_gt = sample_info.split('=', 1)
                    row_dict[sample_name] = sample_gt
            found_rows.append(row_dict)
        else:
            notfound_rows.append(row.drop('_match_key').to_dict())
    
    # Create dataframes
    found_df = pd.DataFrame(found_rows) if found_rows else pd.DataFrame()
    notfound_df = pd.DataFrame(notfound_rows) if notfound_rows else pd.DataFrame()
    
    # Save outputs
    found_output = f"{fid}.{original_filename}.fam_bcf.found.tsv"
    notfound_output = f"{fid}.{original_filename}.fam_bcf.notfound.tsv"
    
    if len(found_df) > 0:
        found_df.to_csv(found_output, sep='\\t', index=False)
    else:
        # Create empty file with header
        pd.DataFrame(columns=['chr', 'position', 'ref', 'alt', 'sample_id', 'INFO', 'FORMAT']).to_csv(found_output, sep='\\t', index=False)
    
    if len(notfound_df) > 0:
        notfound_df.to_csv(notfound_output, sep='\\t', index=False)
    else:
        pd.DataFrame(columns=['chr', 'position', 'ref', 'alt', 'sample_id']).to_csv(notfound_output, sep='\\t', index=False)
    
    # Cleanup
    if os.path.exists(bcf_output):
        os.unlink(bcf_output)
    
    print(f"Found {len(found_df)} variant matches in BCF, {len(notfound_df)} not found")
    """
}
