"""
Script to plot empirical likelhood ratio (LR) curve from AlphaMissense bins.

Input: TSV file with likelihood ratios per AlphaMissense score (range 0-1)

Output: Two plots - pathogenic LR curve and benign LR curve.
"""

import argparse
import csv
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="Path to AlphaMissense LR TSV")
parser.add_argument("--output_path", required=True, help="Path to Pathogenic LR curve image")
parser.add_argument("--output_benign", required=True, help="Path to Benign LR curve image")
args = parser.parse_args()

am_scores = []
likelihood_ratios = []

with open(args.input, 'r') as infile:
    reader = csv.DictReader(infile, delimiter = '\t')
    for row in reader:
        am_scores.append(float(row['score']))
        likelihood_ratios.append(float(row['LR']) if row['LR'] != '' else None)

plt.figure()
plt.plot(am_scores, likelihood_ratios)
plt.xlabel('AlphaMissense Score')
plt.ylabel('Likelihood Ratio')
plt.title('Pathogenic LR Curve')
plt.axhline(y = 2.406, color = 'green', linestyle = '--', label = 'Supporting (2.406)')
plt.axhline(y = 5.790, color = 'orange', linestyle = '--', label = 'Moderate (5.790)')
plt.axhline(y = 33.53, color = 'red', linestyle = '--', label = 'Strong (33.53)')
plt.legend()
plt.savefig(args.output_path)
plt.clf()

benign_am_scores = []
benign_likelihood_ratios = []
# Subset to LR < 1 scores only for the benign plot
for i in range(1, len(am_scores)):
    if likelihood_ratios[i] == None or likelihood_ratios[i] < 1:
        benign_am_scores.append(am_scores[i])
        benign_likelihood_ratios.append(likelihood_ratios[i])

plt.figure()
plt.plot(benign_am_scores, benign_likelihood_ratios)
plt.xlabel('AlphaMissense Score')
plt.ylabel('Likelihood Ratio')
plt.title('Benign LR Curve')
plt.axhline(y = 0.4156, color = 'green', linestyle = '--', label = 'Supporting (0.416)')
plt.axhline(y = 0.1727, color = 'orange', linestyle = '--', label = 'Moderate (0.173)')
plt.axhline(y = 0.0298, color = 'red', linestyle = '--', label = 'Strong (0.03)')
plt.legend()
plt.savefig(args.output_benign)

# Thresholds derived from P(path) = 0.0441
thresholds = {'Benign_Strong': 0.0298, 'Benign_Moderate': 0.1727, 'Benign_Supporting': 0.4156, 'Path_Supporting': 2.406, 'Path_Moderate': 5.790, 'Path_Strong': 33.53, 'Path_Very_Strong': 1124}

# Identify AM score cutoffs by finding where the empirical LR curve crosses each threshold
for threshold_label, lr_threshold in thresholds.items():
    for i in range(len(likelihood_ratios)):
        if likelihood_ratios[i-1] == None or likelihood_ratios[i] == None:
            continue
        if likelihood_ratios[i] > lr_threshold and likelihood_ratios[i-1] < lr_threshold:
            print(f'LR crosses "{threshold_label}" ({lr_threshold}) threshold at AlphaMissense Score = {am_scores[i]}.')
            break
        