"""
Script to filter down ClinVar VCF file, loosely modeling the filtering steps in Pejaver et al. 2022. 

Filter 1:
    - Molecular Consequence: Only contain missense variants (SO:0001583)
    - Clinical Significance: Exclude VUS and conflicting interpretations. Only include Benign, Likely Benign, Pathogenic, Likely Pathogenic
    - Review Status: 1+ stars for the calibration set. (remove variants with less than 1 stars)

Filter 2: 
    - Restrict to variants in genes that have at least one P/LP variant

Returns an output VCF ready for annotation.

To run this script, machine must have the ClinVar VCF file downloaded.

Download the ClinVar VCF (GRCh37) from:
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20260530.vcf.gz
"""

import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="Path to ClinVar VCF file")
parser.add_argument("--output", required=True, help="Path to filtered ClinVar VCF file")
args = parser.parse_args()

clnsig_dict = {'Benign': 0, 'Likely_benign': 0, 'Pathogenic': 0, 'Likely_pathogenic': 0}
review_status_allowed = {"criteria_provided,_single_submitter", "criteria_provided,_multiple_submitters,_no_conflicts", "reviewed_by_expert_panel", "practice_guideline"}
path_values = ['Pathogenic', 'Likely_pathogenic']
path_genes = set()
survivor_genes = set()
output_genes = set()
survivors = []
header_lines = []

passed_filter_2_count = 0

# Pass 0: collect all genes that have at least one P/LP variant
# Used in filter 2 to restrict output to genes with known pathogenic variant
with open(args.input, "r") as infile:
    for line in infile:
        if line.startswith("#"):
            continue

        vcf_fields = line.strip().split("\t")
        info_field = vcf_fields[7]
        info_dict = {}

        for item in info_field.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                info_dict[key] = value
                
        clnsig = re.split(r'[/|,]', info_dict.get("CLNSIG", ""))[0]
        gene_ids = [pair.split(':')[0] for pair in info_dict.get("GENEINFO", "").split('|') if ':' in pair]

        if clnsig in path_values:
            path_genes.update(gene_ids)

# Filter 1: keep missense-only, 1+ star, non-VUS variants
with open(args.input, "r") as infile:
    for line in infile:
        if line.startswith("#"):
            header_lines.append(line)
            continue

        vcf_fields = line.strip().split("\t")
        info_field = vcf_fields[7]
        info_dict = {}

        for item in info_field.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                info_dict[key] = value

        # Take only the primary classification, compound entries like "Pathogenic|other" use |, /, or , as delimiters 
        clnsig = re.split(r'[/|,]', info_dict.get("CLNSIG", ""))[0]
        review_status = info_dict.get("CLNREVSTAT", "")
        molecular_consequence = info_dict.get("MC", "")
        gene_ids = [pair.split(':')[0] for pair in info_dict.get("GENEINFO", "").split('|') if ':' in pair]

        if "missense_variant" not in molecular_consequence:
            continue        
        if review_status not in review_status_allowed:
            continue       
        if clnsig not in clnsig_dict:
            continue
        
        survivors.append((line, clnsig, gene_ids))
        survivor_genes.update(gene_ids)

print(f'Passed Filter 1: {len(survivors)}')
print(f'# of Genes: {len(survivor_genes)}')
print('-----------------------------------')

# Filter 2: restrict to variants in genes that have at least one P/LP variant
with open(args.output, "w") as outfile:      
    for header_line in header_lines:
        outfile.write(header_line)  
    for line, clnsig, gene_ids in survivors:
        if any(g in path_genes for g in gene_ids):
            outfile.write(line)
            clnsig_dict[clnsig] += 1
            passed_filter_2_count+=1
            output_genes.update(gene_ids)

        
print(f"Passed Filter 2: {passed_filter_2_count}")
print(f"# of Genes: {len(output_genes)}")
print('-----------------------------------')

for key, value in clnsig_dict.items():
    print(f'{key}: {value}')