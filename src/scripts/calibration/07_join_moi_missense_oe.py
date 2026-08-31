"""
Joins per-variant TSV from script 06 to inheritance mode (AD/AR/XL/mixed/unknown_moi)
using HPO, GenCC, and ClinGen as ordered fallback sources. Also joins missense
O/E and missense O/E upper CI from gnomAD v4.1.1 constraint metrics, prioritizing
MANE select transcript over canonical.

Outputs calibration_dataset.tsv with additional fields: inheritance_mode,
mis_oe, mis_oe_ci_upper.
"""

import argparse
import csv
import os

parser = argparse.ArgumentParser()
parser.add_argument("--input_tsv", required=True, help="Path to per-variant TSV")
parser.add_argument("--input_sources", required=True, help="Path to disease_gene_to_sources")
parser.add_argument("--constraints", required=True, help="Path to gnomad.v4.1.1.constraint_metrics.tsv")
parser.add_argument("--output", required=True, help="Path to output calibration_dataset.tsv")
args = parser.parse_args()

# Maps HPO term IDs to AD/AR/XL
hpo_term_to_moi = {
    "HP:0000006": "AD",
    "HP:0000007": "AR",
    "HP:0001417": "XL",
    "HP:0001419": "XL",
    "HP:0001423": "XL",
    "HP:0012275": "AD",
    "HP:0034341": "AR",
    "HP:0032113": "AD",  # Semidominant
}

# Maps ClinGen two-letter MOI codes to AD/AR/XL
clingen_code_to_moi = {
    "AD": "AD",
    "AR": "AR",
    "XL": "XL",
    "SD": "AD",  # semidominant
}

gene_to_moi_hpo = {}

hpo_path = os.path.join(args.input_sources, "hpo_genes_to_phenotype.txt")
with open(hpo_path, 'r') as f:
    next(f)
    for line in f:
        cols = line.strip().split('\t')
        gene_symbol = cols[1]
        hpo_id = cols[2]

        if hpo_id not in hpo_term_to_moi:
            continue

        moi = hpo_term_to_moi[hpo_id]
        if gene_symbol not in gene_to_moi_hpo:
            gene_to_moi_hpo[gene_symbol] = set()
        gene_to_moi_hpo[gene_symbol].add(moi)

# Genes appearing with multiple MOI labels are marked mixed
for gene, moi_set in gene_to_moi_hpo.items():
    gene_to_moi_hpo[gene] = moi_set.pop() if len(moi_set) == 1 else 'mixed'

print(f"HPO genes loaded: {len(gene_to_moi_hpo)}")

gene_to_moi_gencc = {}

gencc_path = os.path.join(args.input_sources, "gencc_submissions.tsv")
with open(gencc_path, 'r') as f:
    next(f)
    for line in f:
        cols = line.strip().split('\t')
        if len(cols) < 11:
            continue
        gene_symbol = cols[2]
        moi_curie = cols[9]

        if moi_curie not in hpo_term_to_moi:
            continue

        moi = hpo_term_to_moi[moi_curie]
        if gene_symbol not in gene_to_moi_gencc:
            gene_to_moi_gencc[gene_symbol] = set()
        gene_to_moi_gencc[gene_symbol].add(moi)

for gene, moi_set in gene_to_moi_gencc.items():
    gene_to_moi_gencc[gene] = moi_set.pop() if len(moi_set) == 1 else 'mixed'

print(f"GenCC genes loaded: {len(gene_to_moi_gencc)}")

gene_to_moi_clingen = {}

clingen_path = os.path.join(args.input_sources, "clingen_gene_disease_summary.csv")
with open(clingen_path, 'r') as f:
    # Skip header
    for _ in range(4):
        next(f)
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        if len(row) < 5:
            continue
        gene_symbol = row[0].strip()
        moi = row[4].strip()

        if moi not in clingen_code_to_moi:
            continue

        moi = clingen_code_to_moi[moi]
        if gene_symbol not in gene_to_moi_clingen:
            gene_to_moi_clingen[gene_symbol] = set()
        gene_to_moi_clingen[gene_symbol].add(moi)

for gene, moi_set in gene_to_moi_clingen.items():
    gene_to_moi_clingen[gene] = moi_set.pop() if len(moi_set) == 1 else 'mixed'

print(f"ClinGen genes loaded: {len(gene_to_moi_clingen)}")

gene_to_moi = {}

all_disease_genes = set(gene_to_moi_hpo) | set(gene_to_moi_gencc) | set(gene_to_moi_clingen)

# Aggregate three sources into final lookup
# Priority order: HPO > GenCC > ClinGen
# Genes with disease association but no resolvable MOI -> unknown_moi
# Genes absent from all three sources -> unknown_disease (applied at join time)
for gene in all_disease_genes:
    gene_to_moi[gene] = gene_to_moi_hpo.get(gene) or gene_to_moi_gencc.get(gene) or gene_to_moi_clingen.get(gene) or 'unknown_moi'

# Build gene -> (mis.oe, mis.oe_ci.upper) lookup, prioritizing MANE select > canonical
gene_to_mis_oe = {}
gene_to_mis_oe_ci_upper = {}
gene_to_source = {}

with open(args.constraints, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        gene = row['gene']
        mis_oe = row['mis.oe']
        mis_oe_upper = row['mis.oe_ci.upper']
        is_mane = row['mane_select'] == 'true'
        is_canonical = row['canonical'] == 'true'

        if mis_oe == 'NA' or mis_oe == '' or mis_oe_upper == 'NA' or mis_oe_upper == '':
            continue

        if is_mane:
            gene_to_mis_oe[gene] = float(mis_oe)
            gene_to_mis_oe_ci_upper[gene] = float(mis_oe_upper)
            gene_to_source[gene] = 'mane_select'
        elif is_canonical and gene_to_source.get(gene) != 'mane_select':
            gene_to_mis_oe[gene] = float(mis_oe)
            gene_to_mis_oe_ci_upper[gene] = float(mis_oe_upper)
            gene_to_source[gene] = 'canonical'

print(f"\nBuilt mis.oe lookup: {len(gene_to_mis_oe)} genes")
print(f"  MANE select: {sum(1 for v in gene_to_source.values() if v == 'mane_select')}")
print(f"  Canonical (no MANE): {sum(1 for v in gene_to_source.values() if v == 'canonical')}")

# Join MOI and constraint metrics onto calibration dataset
n_joined = 0
n_missing_constraint = 0
missing_genes = set()

with open(args.input_tsv, 'r') as infile, open(args.output, 'w', newline='') as outfile:
    reader = csv.reader(infile, delimiter='\t')
    writer = csv.writer(outfile, delimiter='\t')

    header = next(reader)
    writer.writerow(header + ['inheritance_mode', 'mis_oe', 'mis_oe_ci_upper'])

    for row in reader:
        gene_symbol = row[4]
        inheritance_mode = gene_to_moi.get(gene_symbol, 'unknown_disease')

        mis_oe = gene_to_mis_oe.get(gene_symbol)
        mis_oe_upper = gene_to_mis_oe_ci_upper.get(gene_symbol)

        if mis_oe is not None and mis_oe_upper is not None:
            writer.writerow(row + [inheritance_mode, f"{mis_oe:.6e}", f"{mis_oe_upper:.6e}"])
            n_joined += 1
        else:
            missing_genes.add(gene_symbol)
            n_missing_constraint += 1

print(f"\nJoined: {n_joined} variants")
print(f"Dropped (no mis.oe): {n_missing_constraint} variants from {len(missing_genes)} genes")
if missing_genes and len(missing_genes) <= 20:
    print(f"  Missing genes: {sorted(missing_genes)}")
