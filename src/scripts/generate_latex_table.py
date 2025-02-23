#!/usr/bin/env python3

"""
Generate a LaTeX table comparing sensitivity, specificity, precision, and F1-score 
between BIAS-2015 and InterVar, highlighting the higher value in bold.
"""

import argparse
import pandas as pd

def load_data(file_path):
    """
    Load sensitivity and specificity data from a tab-separated file.

    Args:
        file_path (str): Path to the TSV file.

    Returns:
        pd.DataFrame: DataFrame containing the parsed data.
    """
    return pd.read_csv(file_path, sep="\t")

def bold_higher_value(bias_value, intervar_value):
    """
    Compare two values and return a LaTeX-formatted string with the higher value bolded.

    Args:
        bias_value (float): Value from the BIAS-2015 dataset.
        intervar_value (float): Value from the InterVar dataset.

    Returns:
        tuple: Formatted strings for LaTeX table.
    """
    bias_formatted = f"\\textbf{{{bias_value:.2%}}}" if bias_value > intervar_value else f"{bias_value:.2%}"
    intervar_formatted = f"\\textbf{{{intervar_value:.2%}}}" if intervar_value > bias_value else f"{intervar_value:.2%}"
    return bias_formatted, intervar_formatted

def generate_latex_table(bias_df, intervar_df, output_file):
    """
    Generate a LaTeX table comparing BIAS-2015 and InterVar classification performance.

    Args:
        bias_df (pd.DataFrame): DataFrame with BIAS-2015 data.
        intervar_df (pd.DataFrame): DataFrame with InterVar data.
        output_file (str): Path to the output LaTeX file.
    """
    category_mapping = {
        "pathogenic": "P/LP",
        "benign": "B/LB",
        "uncertain": "VUS"
    }

    with open(output_file, "w") as f:
        f.write("\\label{tab:sensitivity_specificity}\n")
        f.write("\\centering\n")
        f.write("\\begin{minipage}{\\textwidth}\n")
        f.write("\\small\n")
        f.write("\\tabcolsep=5pt\n")
        f.write("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lcccccccc@{}}\n")
        f.write("\\toprule\n")
        f.write("\\multirow{2}{*}{\\textbf{Category}} & \\multicolumn{4}{c}{\\textbf{BIAS-2015 v2.0.0}} & \\multicolumn{4}{c}{\\textbf{InterVar}} \\\\\n")
        f.write("\\cmidrule(lr){2-5} \\cmidrule(lr){6-9}\n")
        f.write(" & \\textbf{Sensitivity} & \\textbf{Specificity} & \\textbf{Precision} & \\textbf{F1 Score} & ")
        f.write("\\textbf{Sensitivity} & \\textbf{Specificity} & \\textbf{Precision} & \\textbf{F1 Score} \\\\\n")
        f.write("\\midrule\n")

        for category in category_mapping.keys():
            bias_row = bias_df[bias_df["category"] == category].iloc[0]
            intervar_row = intervar_df[intervar_df["category"] == category].iloc[0]

            sens_bias, sens_intervar = bold_higher_value(bias_row["sensitivity"], intervar_row["sensitivity"])
            spec_bias, spec_intervar = bold_higher_value(bias_row["specificity"], intervar_row["specificity"])
            prec_bias, prec_intervar = bold_higher_value(bias_row["precision"], intervar_row["precision"])
            f1_bias, f1_intervar = bold_higher_value(bias_row["f1_score"], intervar_row["f1_score"])

            f.write(f"{category_mapping[category]} & {sens_bias} & {spec_bias} & {prec_bias} & {f1_bias} ".replace("%", "\%"))
            f.write(f"& {sens_intervar} & {spec_intervar} & {prec_intervar} & {f1_intervar} \\\\\n".replace("%","\%"))

        f.write("\\bottomrule\n")
        f.write("\\end{tabular*}\n")
        f.write("\\end{minipage}\n")
        f.write("\\end{table*}\n")

def main():
    """
    Main function to parse arguments and generate the LaTeX table.
    """
    parser = argparse.ArgumentParser(description="Generate LaTeX table comparing BIAS-2015 and InterVar performance.")
    parser.add_argument("bias_file", help="Path to the BIAS-2015 sensitivity/specificity TSV file.")
    parser.add_argument("intervar_file", help="Path to the InterVar sensitivity/specificity TSV file.")
    parser.add_argument("output_file", help="Path to the output LaTeX file.")

    args = parser.parse_args()

    bias_df = load_data(args.bias_file)
    intervar_df = load_data(args.intervar_file)

    generate_latex_table(bias_df, intervar_df, args.output_file)

if __name__ == "__main__":
    main()
