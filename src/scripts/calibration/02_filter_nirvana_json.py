"""
Script to filter Nirvana-annotated ClinVar JSON down to only rare variants (gnomAD allAf <= 0.01) that have AlphaMissense scores.
Intended to be run after 01_filter_clinvar_vcf.py.

Returns final JSON file with only rare variants.

Check out the Nirvana 3.18 Documentation!
https://illumina.github.io/NirvanaDocumentation/3.18/
"""

import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="Path to ClinVar Nirvana-annotated file")
parser.add_argument("--output", required=True, help="Path to filtered ClinVar Nirvana-annotated file")
args = parser.parse_args()

passed_variant_count = 0

with open(args.input, 'r') as infile, open(args.output, 'w') as outfile:
    # Skip the opening header line of the Nirvana JSON
    next(infile)
    for line in infile:
        line = line.strip().rstrip(',')

        # Skip empty lines and Nirvana JSON footer artifacts
        if not line or line.startswith('],"genes":[') or line == ']}':
            continue
        
        clinvar_position = json.loads(line)
        if 'variants' not in clinvar_position:
            continue

        variant = clinvar_position['variants'][0]

        # --- DEFAULT: Use this code block if you would use all allele frequency as the gnomAD filter ---
        if ('gnomad' not in variant or variant['gnomad']['allAf'] <= 0.01)  and 'AlphaMissense' in variant:
            outfile.write(line + '\n')
            passed_variant_count += 1
        
        # --- ALTERNATIVE: Use this code block if you would use the population maximum allele frequency as the gnomAD filter ---
        # max_af = 0
        # if 'gnomad' in variant:
        #     for key, value in variant['gnomad'].items():
        #         if 'Af' in key and key!= 'allAf' and key != 'maleAf' and key != 'femaleAf' and key != 'controlsAllAf':
        #             if value > max_af:
        #                 max_af = value
        
        # if ('gnomad' not in variant or max_af <= 0.01)  and 'AlphaMissense' in variant:
        #     outfile.write(line + '\n')
        #     passed_variant_count += 1
    

print(f'Total size: {passed_variant_count}')
