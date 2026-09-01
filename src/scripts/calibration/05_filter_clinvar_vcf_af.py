"""
Filters ClinVar VCF to missense variants with B/LB/P/LP significance and 1+
review stars. Outputs a TSV with variant coordinates and ClinVar fields.

Download the ClinVar VCF (GRCh38) from:
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
"""

import argparse
import csv
import re

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="Path to ClinVar VCF file")
parser.add_argument("--output", required=True, help="Path to filtered TSV")
args = parser.parse_args()

clnsig_allowed = {'Benign', 'Likely_benign', 'Pathogenic', 'Likely_pathogenic'}
review_status_stars = {
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "reviewed_by_expert_panel": 3,
    "practice_guideline": 4
}

# Pass 1: parse all variants from VCF.
all_variants = []
with open(args.input, "r") as infile:
    for line in infile:
        if line.startswith("#"):
            continue
        cols = line.strip().split("\t")
        info = cols[7]
        info_dict = {}
        for item in info.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                info_dict[key] = value

        mc = info_dict.get("MC", "")
        clnsig = re.split(r'[/|,]', info_dict.get("CLNSIG", ""))[0]
        review_status = info_dict.get("CLNREVSTAT", "")
        gene_info = info_dict.get("GENEINFO", "")
        gene_symbol = gene_info.split(":")[0] if gene_info else ""

        all_variants.append({
            'chrom': cols[0], 'pos': cols[1], 'ref': cols[3], 'alt': cols[4],
            'mc': mc, 'clnsig': clnsig, 'review_status': review_status,
            'gene_symbol': gene_symbol,
        })

def variant_stats(variants):
    genes = set(v['gene_symbol'] for v in variants if v['gene_symbol'])
    benign = sum(1 for v in variants if v.get('significance') == 'benign')
    path = sum(1 for v in variants if v.get('significance') == 'pathogenic')
    return len(genes), benign, path

# Pass 2: sequential filter cascade.
# Step 1: Missense filter
total = len(all_variants)
after_missense = [v for v in all_variants if 'missense_variant' in v['mc']]
removed_missense = total - len(after_missense)

# Step 2: Clinical significance (B/LB/P/LP)
after_clnsig = []
rejected_clnsig_to_count = {}
for v in after_missense:
    if v['clnsig'] in clnsig_allowed:
        # Collapse Benign/Likely_benign -> benign, Pathogenic/Likely_pathogenic -> pathogenic
        if v['clnsig'] in ('Likely_benign', 'Benign'):
            v['significance'] = 'benign'
        else:
            v['significance'] = 'pathogenic'
        after_clnsig.append(v)
    else:
        rejected_clnsig_to_count[v['clnsig']] = rejected_clnsig_to_count.get(v['clnsig'], 0) + 1
removed_clnsig = len(after_missense) - len(after_clnsig)

# Step 3: Review status (>= 1 star)
# review_status_stars only contains statuses with 1+ stars, so presence in the
# dict is sufficient to pass the filter.
after_review = []
for v in after_clnsig:
    if v['review_status'] in review_status_stars:
        v['review_stars'] = review_status_stars[v['review_status']]
        after_review.append(v)
removed_review = len(after_clnsig) - len(after_review)

# Print sequential cascade.
print(f"{'Step':<35s} {'Removed':>10s} {'Remaining':>10s} {'Genes':>7s} {'Benign':>8s} {'Path':>8s}")
print("-" * 82)
print(f"{'Total VCF variants':<35s} {'':>10s} {total:>10,d}")

g, b, p = variant_stats(after_missense)
print(f"{'1. Missense filter':<35s} {removed_missense:>10,d} {len(after_missense):>10,d} {g:>7,d}")

g, b, p = variant_stats(after_clnsig)
print(f"{'2. Clinical significance (B/LB/P/LP)':<35s} {removed_clnsig:>10,d} {len(after_clnsig):>10,d} {g:>7,d} {b:>8,d} {p:>8,d}")

g, b, p = variant_stats(after_review)
print(f"{'3. Review status (>= 1 star)':<35s} {removed_review:>10,d} {len(after_review):>10,d} {g:>7,d} {b:>8,d} {p:>8,d}")

# Significance breakdown of final set.
clnsig_to_count = {}
for v in after_review:
    clnsig_to_count[v['clnsig']] = clnsig_to_count.get(v['clnsig'], 0) + 1
print(f"\nFinal breakdown:")
for key, value in sorted(clnsig_to_count.items(), key=lambda x: -x[1]):
    print(f"  {key}: {value:,d}")

if rejected_clnsig_to_count:
    print(f"\nRejected clinical classifications (step 2):")
    for key, value in sorted(rejected_clnsig_to_count.items(), key=lambda x: -x[1]):
        print(f"  {key}: {value:,d}")

# Build output rows.
output_rows = []
for v in after_review:
    output_rows.append([v['chrom'], v['pos'], v['ref'], v['alt'],
                        v['gene_symbol'], v['significance'], v['review_stars']])

headers = ["chromosome", "position", "ref_allele", "alt_allele", "gene_symbol",
           "clinical_significance", "review_stars"]

with open(args.output, 'w', newline="", encoding="utf-8") as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    writer.writerow(headers)
    writer.writerows(output_rows)
