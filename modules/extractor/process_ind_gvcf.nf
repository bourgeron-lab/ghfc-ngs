process PROCESS_IND_GVCF {
    tag "${barcode}_${original_filename}"
    
    publishDir "${params.data}/samples/${barcode}/extractor", mode: 'copy'
    
    input:
    tuple val(barcode), val(original_filename), path(extractor_tsv), path(gvcf), path(gvcf_tbi)
    
    output:
    tuple val(barcode), val(original_filename), path("${barcode}.${original_filename}.ind_gvcf.tsv")
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import pysam
    import os
    
    barcode = '${barcode}'
    original_filename = '${original_filename}'
    
    # Read extractor TSV
    extractor_df = pd.read_csv('${extractor_tsv}', sep='\\t')
    
    # Filter for this sample only
    extractor_df = extractor_df[extractor_df['sample_id'] == barcode].copy()
    
    # Open gVCF
    gvcf_path = '${gvcf}'
    vcf = pysam.VariantFile(gvcf_path)
    
    # Get sample name from gVCF
    samples_in_vcf = list(vcf.header.samples)
    if barcode not in samples_in_vcf:
        # Try to find a matching sample
        sample_name = samples_in_vcf[0] if samples_in_vcf else None
        print(f"Warning: {barcode} not found in gVCF samples {samples_in_vcf}, using {sample_name}")
    else:
        sample_name = barcode
    
    results = []
    
    for _, row in extractor_df.iterrows():
        chrom = str(row['chr'])
        pos = int(row['position'])
        ref = str(row['ref'])
        alt = str(row['alt'])
        
        result_row = {
            'chr': chrom,
            'position': pos,
            'ref': ref,
            'alt': alt,
            'sample_id': barcode,
            'gvcf_FILTER': None,
            'gvcf_FORMAT': None,
            'gvcf_SAMPLE': None
        }
        
        # Query region
        try:
            # Fetch records at this position
            for record in vcf.fetch(chrom, pos - 1, pos):  # pysam uses 0-based
                # Check if this is the right variant
                if record.pos == pos:
                    # Check if ref/alt match
                    record_alts = [str(a) for a in record.alts] if record.alts else []
                    if str(record.ref) == ref and alt in record_alts:
                        # Found exact match
                        result_row['gvcf_FILTER'] = ';'.join(record.filter.keys()) if record.filter else 'PASS'
                        result_row['gvcf_FORMAT'] = ':'.join(record.format.keys())
                        
                        # Get sample data
                        if sample_name and sample_name in record.samples:
                            sample_data = record.samples[sample_name]
                            sample_values = []
                            for key in record.format.keys():
                                val = sample_data.get(key, '.')
                                if isinstance(val, tuple):
                                    val = ','.join(str(v) if v is not None else '.' for v in val)
                                elif val is None:
                                    val = '.'
                                sample_values.append(str(val))
                            result_row['gvcf_SAMPLE'] = ':'.join(sample_values)
                        break
        except Exception as e:
            print(f"Error fetching {chrom}:{pos}: {e}")
        
        results.append(result_row)
    
    vcf.close()
    
    # Create output dataframe
    result_df = pd.DataFrame(results)
    
    # Save output
    output_file = f"{barcode}.{original_filename}.ind_gvcf.tsv"
    result_df.to_csv(output_file, sep='\\t', index=False)
    
    found_count = result_df['gvcf_FILTER'].notna().sum()
    print(f"Found {found_count}/{len(result_df)} variants in gVCF for {barcode}")
    """
}
