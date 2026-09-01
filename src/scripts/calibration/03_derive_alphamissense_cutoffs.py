"""
Script to derive AlphaMissense cutoffs via the local-posterior sliding-window method (Pejaver et al. 2022), reproducing the cutoffs in Bergquist et al. 2025.
This script is a variation of create_alphamissense_bins_json.py, and falls in the same place as this file in our developing pipeline.

Input: data_2024.zip @ https://zenodo.org/records/13766399

Output: Printed PP3/BP4 score cutoffs per evidence strength, adjustable windows + bootstrapping.
A TSV with columns 'score', 'win_lo', 'win_hi', 'n_path', 'n_benign', 'percent_path', 'percent_benign', 'LR' on original variant data.
"""

import argparse
import json
import csv
import random
import numpy as np
import bisect

parser = argparse.ArgumentParser()
parser.add_argument("--clinvar_input", required=True, help="Path to ClinVar variants file")
parser.add_argument("--gnomad_input", required=True, help="Path to gnomAD variants file")
parser.add_argument("--output", required=True, help="Path to output TSV of AlphaMissense-separated variants")
args = parser.parse_args()

# Load labeled ClinVar calibration variants with AlphaMissense scores.
path_labels = {'Pathogenic', 'Likely_pathogenic'}
benign_labels = {'Benign', 'Likely_benign'}
clinvar_variants = []
with open(args.clinvar_input, 'r') as infile:
    reader = csv.DictReader(infile, delimiter = '\t')
    for row in reader:
        if row['clnsig'] not in path_labels and row['clnsig'] not in benign_labels:
            continue
        if row['AlphaMissense_score'] == '.' or row['AlphaMissense_score'] == '':
            continue      
        clinvar_variants.append(row)


# --- ALTERNATIVE: Load labeled ClinVar calibration variants with AlphaMissense scores IN JSON FORMAT. (unused) ---
# path_labels = {'pathogenic', 'likely pathogenic'}
# benign_labels = {'benign', 'likely benign'}
# clinvar_variants = []
# with open(args.clinvar_input, 'r') as infile:
#     for line in infile:
#         line = line.strip().rstrip(',')
#         if not line or line.startswith('],"genes":[') or line == ']}':
#             continue       
#         position = json.loads(line)
#         variant = position['variants'][0]
#         if 'clinvar' not in variant:
#             continue
#         final_call = None
#         for call in variant['clinvar']:
#             if call['id'].startswith('VCV'):
#                 final_call = call
#         if final_call == None:
#             continue
#         significance = final_call['significance'][0]
#         clinvar_variants.append({'clnsig': significance, 'AlphaMissense_score': variant['AlphaMissense']['AM_score']})


# Load rare gnomAD AlphaMissense scores.
# These are used to determine adaptive window width (not to compute LR).
gnomad_scores = []
with open(args.gnomad_input, 'r') as infile:
    next(infile)
    for line in infile:
        am = line.rstrip('\n').split('\t')[27]
        if am == '.' or am == '':
            continue
        gnomad_scores.append(float(am))

gnomad_scores = sorted(gnomad_scores)

def count_between(sorted_scores, low, high):
    """Count values in a sorted score list that fall within [low, high] using bisect library."""
    return bisect.bisect_right(sorted_scores, high) - bisect.bisect_left(sorted_scores, low)

def compute_adaptive_lr_curve(variants):
    """Compute the local LR at each unique AlphaMissense score using an adpative sliding
    window over the labeled ClinVar data. Returns rows of 
    [score, win_lo, win_hi, n_path, n_benign, percent_path, percent_benign, LR]."""
    pathogenic_scores = sorted([float(v['AlphaMissense_score']) for v in variants if v['clnsig'] in path_labels])
    benign_scores = sorted([float(v['AlphaMissense_score']) for v in variants if v['clnsig'] in benign_labels])

    lr_curve_rows = []

    # Center one adaptive sliding window on each unique ClinVar AlphaMissense score.
    for unique_score in sorted(set(pathogenic_scores) | set(benign_scores)):
        window_radius = 0.0
        # Expand the window until it contains enough ClinVar labels and gnomAD score density.
        while True:
            window_low, window_high = unique_score - window_radius, unique_score + window_radius
            full_width = window_high - window_low
            # Near score boundaries, scale down the required counts because the window is clipped at 0 or 1.
            edge_scale = (min(window_high, 1) - max(window_low, 0)) / full_width if full_width > 0 else 1.0
            if (
                count_between(pathogenic_scores, window_low, window_high) + count_between(benign_scores, window_low, window_high) >= 100 * edge_scale and 
                count_between(gnomad_scores, window_low, window_high) >= 0.03 * len(gnomad_scores) * edge_scale
            ):
                break
            window_radius += 0.0001

        pathogenic_count = count_between(pathogenic_scores, window_low, window_high)
        benign_count = count_between(benign_scores, window_low, window_high)
        pathogenic_window_fraction = pathogenic_count / len(pathogenic_scores)
        benign_window_fraction = benign_count / len(benign_scores)
        # Estimate local LR as the fraction of pathogenic variants in the window divided by the fraction of benign variants in the window.
        lr = pathogenic_window_fraction / benign_window_fraction if benign_count > 0 else float("inf")

        lr_curve_rows.append([unique_score, window_low, window_high, pathogenic_count, benign_count, pathogenic_window_fraction * 100, benign_window_fraction * 100, lr])
    
    return lr_curve_rows

n_clinvar_variants = len(clinvar_variants)
thresholds = {'Benign_Strong': 0.0298, 'Benign_Moderate': 0.1727, 'Benign_Supporting': 0.4156, 'Path_Supporting': 2.406, 'Path_Moderate': 5.790, 'Path_Strong': 33.53}
bootstrap_threshold_scores = {'Benign_Strong': [], 'Benign_Moderate': [], 'Benign_Supporting': [], 'Path_Supporting': [], 'Path_Moderate': [], 'Path_Strong': []}

# Bootstrap ClinVar variants to estimate uncertainty in each score cutoff and smooth threshold estimates. 
for i in range(500):
    boot_sample = random.choices(clinvar_variants, k=n_clinvar_variants)
    bootstrap_lr_curve = compute_adaptive_lr_curve(boot_sample)
    for key, value in thresholds.items():
        lr_values = [row[-1] for row in bootstrap_lr_curve]
        scores = [row[0] for row in bootstrap_lr_curve]
        
        if key.startswith('Path'):
            last_failure = -1
            for j in range(len(lr_values)):
                if lr_values[j] < value:
                    last_failure = j
            if last_failure < len(lr_values) - 1:
                bootstrap_threshold_scores[key].append(scores[last_failure + 1])
        else:
            last_failure = len(lr_values)
            for j in range(len(lr_values) - 1, -1, -1):
                if lr_values[j] > value:
                    last_failure = j
            if last_failure > 0:
                bootstrap_threshold_scores[key].append(scores[last_failure - 1])

for key, value in bootstrap_threshold_scores.items():
    pct = 95 if key.startswith('Path') else 5
    # Use conservative bootstrap percentiles for final threshold estimates (95 for Pathogenic and 5 for Benign).
    print(key, np.percentile(value, pct) if value else None)

lr_curve_original_variants = compute_adaptive_lr_curve(clinvar_variants)
headers = ['score', 'win_lo', 'win_hi', 'n_path', 'n_benign', 'percent_path', 'percent_benign', 'LR']

# Write the adaptive LR curve from the original ClinVar data.
with open(args.output, 'w', newline = "", encoding = "utf-8") as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    writer.writerow(headers)
    writer.writerows(lr_curve_original_variants)