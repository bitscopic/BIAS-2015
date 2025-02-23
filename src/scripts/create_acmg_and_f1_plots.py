"""
Takes in files produced by evaluate_bias_erepo_validation and generates the figures of the style
used in the BIAS v2.0.0 manuscript.
"""
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams["font.family"] = "serif"

def plot_classifier_assignments(bias_filtered, intervar_filtered, output_file_1):
    """ Generate bar chart comparing classifier assignment counts (Total Calls) """
    classifiers = bias_filtered["code"].tolist()
    bias_counts = (bias_filtered["alg_correct"] + bias_filtered["alg_fp"]).tolist()
    intervar_counts = (intervar_filtered["alg_correct"] + intervar_filtered["alg_fp"]).tolist()
    erepo_counts = [bias + fn for bias, fn in zip(bias_filtered["alg_correct"], bias_filtered["alg_fn"])]

    x = np.arange(len(classifiers))
    width = 0.25  # Bar width

    _, ax = plt.subplots(figsize=(12, 6))

    # Bit thicker lines to make it punchier on the page
    ax.spines["top"].set_linewidth(1.5)
    ax.spines["right"].set_linewidth(1.5)
    ax.spines["bottom"].set_linewidth(1.5)
    ax.spines["left"].set_linewidth(1.5)


    
    ax.bar(x - width, bias_counts, width, label="BIAS-2015", color="#2a4db8", alpha=0.7)
    ax.bar(x, intervar_counts, width, label="InterVar", color="#d4af37", alpha=0.7)
    ax.bar(x + width, erepo_counts, width, label="eRepo (Ground Truth)", color="#4b8251", alpha=0.7)

    ax.set_xlabel("ACMG Classifier Codes", fontsize=14)
    ax.set_ylabel("Number of Assignments", fontsize=14)
    ax.set_title("ACMG Classifier Assignments: BIAS-2015 vs. InterVar vs. eRepo", fontsize=16)
    ax.set_xticks(x)
    ax.set_xticklabels(classifiers, rotation=45, ha="right", fontsize=12)
    ax.legend(fontsize=12)

    # Add horizontal grid lines
    ax.yaxis.grid(True, linestyle="--", linewidth=0.7, alpha=0.4)
    
    # Dynamically adjust minor tick spacing
    ax.yaxis.set_minor_locator(plt.MultipleLocator(250))  # Minor ticks every 250 units
    ax.yaxis.grid(True, which="minor", linestyle=":", linewidth=0.5, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file_1, dpi=300)
    print(f"Classifier assignments chart saved to {output_file_1}")

def plot_f1_concordance(bias_filtered, intervar_filtered, output_file_2):
    """ Generate bar chart comparing F1 scores and Concordance values """
    classifiers = bias_filtered["code"].tolist()
    bias_f1 = bias_filtered["f1_score"].tolist()
    intervar_f1 = intervar_filtered["f1_score"].tolist()

    x = np.arange(len(classifiers))
    width = 0.33  # Bar width

    _, ax = plt.subplots(figsize=(12, 6))
    
    # Bit thicker lines to make it punchier on the page
    ax.spines["top"].set_linewidth(1.5)
    ax.spines["right"].set_linewidth(1.5)
    ax.spines["bottom"].set_linewidth(1.5)
    ax.spines["left"].set_linewidth(1.5)
    
    ax.bar(x - width*1.5, bias_f1, width, label="BIAS-2015 F1", color="#2a4db8", alpha=0.7)
    ax.bar(x - width/2, intervar_f1, width, label="InterVar F1", color="#d4af37", alpha=0.7)

    ax.set_xlabel("ACMG Classifier Codes", fontsize=14)
    ax.set_ylabel("Score (0-1)", fontsize=14)
    ax.set_title("ACMG Classifier Performance: F1 Score", fontsize=16)
    ax.set_xticks(x)
    ax.set_xticklabels(classifiers, rotation=45, ha="right", fontsize=12)
    ax.legend(fontsize=12)

    # Add horizontal grid lines
    ax.yaxis.grid(True, linestyle="--", linewidth=0.7, alpha=0.4)
    
    # Dynamically adjust minor tick spacing
    ax.yaxis.set_minor_locator(plt.MultipleLocator(.05))  # Minor ticks every 250 units
    ax.yaxis.grid(True, which="minor", linestyle=":", linewidth=0.5, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file_2, dpi=300)
    print(f"F1 & Concordance chart saved to {output_file_2}")

def main(bias_file, intervar_file, output_file_1, output_file_2):
    """
    Load data then generate graphs
    """
    # Load the data
    bias_df = pd.read_csv(bias_file, sep="\t")
    intervar_df = pd.read_csv(intervar_file, sep="\t")

    # List of ACMG codes to exclude
    excluded_codes = {"PS2", "PM3", "PM6", "PP1", "PP4", "BS3", "BS4", "BP2", "BP5"}

    # Apply filtering to remove unwanted classifiers
    bias_filtered = bias_df[~bias_df["code"].isin(excluded_codes)]
    intervar_filtered = intervar_df[~intervar_df["code"].isin(excluded_codes)]

    # Generate both bar charts
    plot_classifier_assignments(bias_filtered, intervar_filtered, output_file_1)
    
    excluded_codes = {"PS2", "PM3", "PM6", "PP1", "PP4", "PP5", "BS3", "BS4", "BP2", "BP5", "BP6"}
    # Apply filtering to remove unwanted classifiers
    bias_filtered = bias_df[~bias_df["code"].isin(excluded_codes)]
    intervar_filtered = intervar_df[~intervar_df["code"].isin(excluded_codes)]
    plot_f1_concordance(bias_filtered, intervar_filtered, output_file_2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate bar charts comparing ACMG classifier assignments, F1 scores, and concordance between BIAS-2015 and InterVar.")
    parser.add_argument("bias_file", help="Path to bias_conc.tsv")
    parser.add_argument("intervar_file", help="Path to intervar_conc.tsv")
    parser.add_argument("output_file_1", help="Path to save the classifier assignment chart (e.g., assignments.png)")
    parser.add_argument("output_file_2", help="Path to save the F1 & concordance chart (e.g., f1_concordance.png)")

    args = parser.parse_args()
    main(args.bias_file, args.intervar_file, args.output_file_1, args.output_file_2)
