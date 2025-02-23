"""
Generate combined metrics tables for supplementary documents
"""
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def load_and_merge_data(bias_path: str, intervar_path: str) -> pd.DataFrame:
    """Load BIAS-2015 and InterVar comparison files, merge on code, and return a DataFrame."""
    bias_df = pd.read_csv(bias_path, sep="\t")
    intervar_df = pd.read_csv(intervar_path, sep="\t")

    # Rename 'code' to ensure it remains the merge key
    bias_df = bias_df.rename(columns={"code": "code"}).add_prefix("BIAS_")
    intervar_df = intervar_df.rename(columns={"code": "code"}).add_prefix("INTERVAR_")

    # Fix column names for proper merge
    bias_df = bias_df.rename(columns={"BIAS_code": "code"})
    intervar_df = intervar_df.rename(columns={"INTERVAR_code": "code"})

    # Merge on 'code'
    merged_df = pd.merge(bias_df, intervar_df, on="code")
    merged_df = pd.merge(bias_df, intervar_df, on="code").drop(columns=["BIAS_algorithm", "INTERVAR_algorithm"])

    return merged_df


def save_tsv(df: pd.DataFrame, output_path: str) -> None:
    """Save the merged DataFrame as a TSV file."""
    df.to_csv(output_path, sep="\t", index=False)


def save_table_figure(df: pd.DataFrame, output_path: str) -> None:
    """Create and save a table figure from the DataFrame using Matplotlib."""
    _, ax = plt.subplots(figsize=(12, len(df) * 0.4))
    ax.axis("tight")
    ax.axis("off")

    # Define table column names
    columns = df.columns.tolist()

    # Convert to list of lists for matplotlib table
    table_data = [df.columns.tolist()] + df.values.tolist()

    # Create table in figure
    table = ax.table(cellText=table_data, colLabels=None, loc="center", cellLoc="center")

    # Style table
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.auto_set_column_width(list(range(len(columns))))

    plt.savefig(output_path, dpi=300, bbox_inches="tight")


def main():
    """Main function to parse arguments, process files, and save outputs."""
    parser = argparse.ArgumentParser(description="Merge BIAS-2015 and InterVar comparison data, save as TSV and image.")
    parser.add_argument("bias_file", type=str, help="Path to BIAS-2015 comparison TSV file.")
    parser.add_argument("intervar_file", type=str, help="Path to InterVar comparison TSV file.")
    parser.add_argument("output_tsv", type=str, help="Path to save the merged TSV file.")
    parser.add_argument("output_png", type=str, help="Path to save the table figure PNG file.")

    args = parser.parse_args()

    merged_df = load_and_merge_data(args.bias_file, args.intervar_file)
    save_tsv(merged_df, args.output_tsv)
    save_table_figure(merged_df, args.output_png)

    print(f"Supplementary table saved to {args.output_tsv} and {args.output_png}")


if __name__ == "__main__":
    main()
