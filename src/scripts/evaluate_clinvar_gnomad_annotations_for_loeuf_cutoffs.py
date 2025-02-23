"""
Take a clinvar vcf nirvana json file and identify for each gene, the LOEUF and the number of pathogenic/benign variants
found in the gene, along with the gnomad frequencies for each variant.

This will allow us to correlate LOEUF with the observed gnomad frequencies for pathogenic/benign variants. This lets us
assign an average gnomad frequency for pathogenic variants, and an average gnomad frequency for benign variants.

Once per-gene LOEUF and allele frequency data is extracted, we apply clustering methods to identify natural groupings
of allele frequencies. Using log-transformed allele frequencies, we cluster variants into frequency bins using 
K-means. This allows us to define empirical LOEUF thresholds by observing where 
benign and pathogenic variants naturally separate across gene constraint levels.

By analyzing the distribution of LOEUF values across clustered allele frequencies, we derive breakpoints that 
differentiate highly constrained genes from tolerant ones. This removes arbitrary LOEUF cutoffs and instead 
establishes data-driven thresholds.
"""
import json
import sys
import os
import io
import gzip
import argparse
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr
from sklearn.mixture import GaussianMixture


# Add the parent directory of "scripts" to the Python path to make "pipeline_runner_core" importable
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from src.bias_2015 import extract_from_nirvana_json
from src.bias_2015.constants import clinvar_review_status_to_level 


def parseArgs(): 
    """
    Parse the command line arguments into useful python objects.  '--' variables are optional
    set_defaults only applies when the argument is not provided (it won't override)
    """
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument("clinvar_nirvana_json",
                        help = " Input file 1",
                        action = "store")
    parser.add_argument("in_vcf",
                        help = " Input file 1",
                        action = "store")
    parser.add_argument("output_file",
                        help = " Input file 1",
                        action = "store")
    parser.add_argument("--verbose",
                        help = " The verbosity level for stdout messages (default INFO)",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        action = "store")

    parser.set_defaults(verbose = "WARNING")

    options = parser.parse_args()
    logging.basicConfig(level=getattr(logging, options.verbose), format='%(message)s')
    return options


def open_file(file_path, mode):
    """
    Open either a normal or a .gz file
    """
    _, file_extension = os.path.splitext(file_path)
    if file_extension == ".gz":
        return gzip.open(file_path, mode)
    return io.open(file_path, mode, encoding="utf-8")


def get_clinvar_var_to_data(in_vcf):
    """
    Expects a clinvar formatted vcf file with 'chr' prepended to the chromosome name

    Creates a dictionary mapping the variant (chrom, pos, ref, alt) to other additional information
    from ClinVar annotating the variant
    """
    clinvar_var_to_data = {}
    with open(in_vcf) as in_file:
        for line in in_file:
            if not line: continue
            if line.startswith("#"): continue
            split_line = line.strip().split()
            chrom = split_line[0]
            pos = split_line[1]
            clinvar_id = split_line[2]
            ref = split_line[3]
            alt = split_line[4]
            variant = (chrom, pos, ref, alt)
            term_to_value = {}
            elements = split_line[7].split(";")
            for pair in elements:
                term, value = pair.split("=")
                term_to_value[term] = value
            rs_id = term_to_value.get('RS', '')
            criteria = term_to_value.get('CLNREVSTAT', '')
            signif = term_to_value.get('CLNSIG', '')
            if not criteria: continue
            if not signif: continue
            clinvar_var_to_data[variant] = (rs_id, criteria, signif, clinvar_id)
    return clinvar_var_to_data

def get_cons_to_var_list(var_list):
    """
    Reformat into a dictionary mapping consequence to lists of vars
    """
    cons_to_var_list = {}
    for var in var_list:
        consequence = var[6]
        if cons_to_var_list.get(consequence):
            cons_to_var_list[consequence].append(var)
        else:
            cons_to_var_list[consequence] = [var]
    return cons_to_var_list


def identify_loeuf_cutoff_values(var_data_list, plot_filename):
    """
    Clusters pathogenic variant AFs by LOEUF bins to assess correlation.

    Determines LOEUF breakpoints using Gaussian Mixture Model (GMM).
    Computes mean and median AF per LOEUF bin.
    Evaluates distribution of variants across bins to check for imbalances.
    Also computes AF statistics using quantile-based LOEUF binning (qcut).
    """
    # Convert path_vars list to DataFrame
    df = pd.DataFrame(var_data_list, columns=["gene", "af", "loeuf", "clinvar_id", "consequence"])

    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["loeuf", "af"])

    print(f"Total number of variants in df: {len(df)}")

    # Compute correlation between LOEUF and AF
    spearman_corr, spearman_p = spearmanr(df["loeuf"], df["af"])
    pearson_corr, pearson_p = pearsonr(df["loeuf"], df["af"])
    print(f"Spearman correlation: {spearman_corr:.3f}, p-value: {spearman_p:.3g}")
    print(f"Pearson correlation: {pearson_corr:.3f}, p-value: {pearson_p:.3g}")

    # Try different cluster counts (e.g., 2 to 10)
    cluster_range = range(2, 11)
    bic_scores = []
    loeuf_data = df[["loeuf"]].values

    for n in cluster_range:
        gmm = GaussianMixture(n_components=n, random_state=42)
        gmm.fit(loeuf_data)
        bic_scores.append(gmm.bic(loeuf_data))

    # Find the optimal cluster count (lowest BIC)
    optimal_clusters = cluster_range[np.argmin(bic_scores)]
    print(f"Optimal number of clusters based on BIC: {optimal_clusters}")

    # Fit GMM with the optimal cluster count
    gmm = GaussianMixture(n_components=optimal_clusters, random_state=42)
    df["loeuf_cluster"] = gmm.fit_predict(loeuf_data)

    # Extract LOEUF breakpoints
    loeuf_breakpoints = np.sort(gmm.means_.flatten()).tolist()
    print("Dynamically derived LOEUF breakpoints:", loeuf_breakpoints)

    # Assign LOEUF bins using dynamically derived breakpoints
    loeuf_bins = [0] + loeuf_breakpoints + [df["loeuf"].max() + 0.1]
    df["loeuf_bin_gmm"] = pd.cut(df["loeuf"], bins=loeuf_bins, labels=[f"bin_{i}" for i in range(len(loeuf_breakpoints) + 1)])

    # Compute AF stats for GMM-based bins
    af_summary_gmm = df.groupby("loeuf_bin_gmm", observed=True)["af"].agg(["mean", "median", "std", "min", "max"])
    af_summary_gmm["count"] = df["loeuf_bin_gmm"].value_counts().reindex(af_summary_gmm.index, fill_value=0)
    print("GMM-based LOEUF bin AF summary:")
    print(af_summary_gmm)

    # Generate correct number of labels for bins
    df["loeuf_bin_qcut"], bin_edges = pd.qcut(df["loeuf"], q=10, retbins=True, duplicates="drop")

    # Print LOEUF bin edges
    print("Quantile-based LOEUF bin edges:", bin_edges)

    # Compute AF stats for quantile-based bins
    af_summary_qcut = df.groupby("loeuf_bin_qcut", observed=True)["af"].agg(["mean", "median", "std", "min", "max"])
    af_summary_qcut["count"] = df["loeuf_bin_qcut"].value_counts().reindex(af_summary_qcut.index, fill_value=0)
    print("Quantile-based LOEUF bin AF summary:")
    print(af_summary_qcut)

        # Scatter plot of LOEUF vs AF
    # Function to generate scatter plot with trend lines
    def plot_af_vs_loeuf(filename, log_scale=True):
        plt.figure(figsize=(8, 6))
        plt.scatter(df["loeuf"], df["af"], alpha=0.5, s=10, label="Variants")

        mean_loeuf = df["loeuf"].mean()
        mean_af = df["af"].mean()
        plt.axhline(mean_af, color="red", linestyle="--", linewidth=1, label=f"Mean AF: {mean_af:.2e}")
        plt.axvline(mean_loeuf, color="blue", linestyle="--", linewidth=1, label=f"Mean LOEUF: {mean_loeuf:.2f}")

        # Function to compute and plot binned averages
        def plot_binned_means(bin_width, color, label):
            bins = np.arange(0, df["loeuf"].max() + bin_width, bin_width)
            df["loeuf_bin"] = pd.cut(df["loeuf"], bins=bins, right=False)
            binned_af = df.groupby("loeuf_bin", observed=True)["af"].mean().dropna()

            bin_midpoints = [(interval.left + interval.right) / 2 for interval in binned_af.index]
            plt.plot(bin_midpoints, binned_af.values, color=color, linewidth=2, label=label)

        # Plot two different bin sizes
        plot_binned_means(0.05, "green", "Mean AF in 0.05 LOEUF bins")
        plot_binned_means(0.15, "orange", "Mean AF in 0.15 LOEUF bins")

        plt.xlabel("LOEUF")
        plt.ylabel("AF")
        plt.title(f"Scatter Plot of AF vs LOEUF ({'Log Scale' if log_scale else 'Linear Scale'})")

        if log_scale:
            plt.yscale("log")  # Log scale for better visualization

        plt.legend()
        plt.grid(True)
        plt.savefig(filename)
        plt.close()

    # Generate both versions
    if plot_filename:
        plot_af_vs_loeuf(plot_filename, log_scale=True)

# Parse the nirvana output into a more machine friendly format
def extract_aa_information(inJson, clinvar_var_to_data, output_file):
    """
    Parses a Nirvana JSON file and extracts header information, position data, and gene data.

    Args:
        inJson (str): Path to the Nirvana JSON file.

    Returns:
        tuple: A tuple containing the extracted data in the following order:
        - header_json (dict): A dictionary representing the header information.
        - position_list (list): A list of dictionaries representing the position information.
        - hgnc_to_gene_data (dict): A dictionary mapping HGNC IDs to gene data.

    Raises:
        FileNotFoundError: If the input Nirvana JSON file does not exist.
        IOError: If there are any issues while reading or parsing the file.
    """
    skipped = 0 
    gene_to_variant_composition = {}
    hgnc_to_gene_data = extract_from_nirvana_json.load_nirvana_gene_information(inJson) 
    with open_file(inJson, "rt") as f:
        # gather the header line
        _ = f.readline()[10:-15]
        # Read the file line by line
        count = 0
        for line in f:
            # Check if we have reached the start of the genes section
            if '"genes":[' in line:
                break

            # Parse each line as JSON
            if line == ']}\n':
                continue
            if line[-2] == ",":
                try:
                    data = json.loads(line[:-2])
                except:
                    print("skipping line")
            else:
                data = json.loads(line)

            # Extract position information
            position = data
            variant_list = position.get("variants")
            if not data.get("altAlleles"): 
                continue # typically minor reference alleles.
            for variant in variant_list:
                chrom = position["chromosome"]
                pos = str(position["position"])
                var = extract_from_nirvana_json.process_variant(variant, hgnc_to_gene_data, "RefSeq", chrom, pos, position['refAllele'], position['altAlleles'][0])
                v_index = (var.chromosome, var.position, var.refAllele, var.altAllele)
                if not clinvar_var_to_data.get(v_index):
                    skipped += 1
                    continue
                rs_id, criteria, signif, clinvar_id = clinvar_var_to_data[v_index]
                # Don't process intergenic pathogenic variants - they won't have any AA change associated with them
                if var.geneName == 'n/a': 
                    continue
                if not('benign' in signif.lower() or 'pathogenic' in signif.lower()):
                    continue

                values = [
                        var.protein_variant,
                        rs_id,
                        var.gnomad,
                        criteria.replace("_", " "),
                        signif,
                        clinvar_id,
                        var.consequence
                       ]
                if clinvar_review_status_to_level[values[3]] < 1:
                    continue
                if gene_to_variant_composition.get(var.geneName):
                    if 'benign' in signif.lower():
                        gene_to_variant_composition[var.geneName]['benign'].append(values)
                    else:
                        gene_to_variant_composition[var.geneName]['pathogenic'].append(values)
                else:
                    if 'benign' in signif.lower():
                        gene_to_variant_composition[var.geneName] = {'benign': [values], 'pathogenic': []}
                    else:
                        gene_to_variant_composition[var.geneName] = {'pathogenic': [values], 'benign': []}

            count += 1
            if count % 10000 == 0:
                logging.info("processed %i variants", count)
    

    path_vars = []
    path_below_1 = []
    path_above_1 = []
    ben_vars = []
    ben_below_1 = []
    ben_above_1 = []
    gene_to_pheno_to_cons_to_var_list = {}
    for gene, path_ben in gene_to_variant_composition.items():
        path_list = path_ben['pathogenic']
        ben_list = path_ben['benign']
        gene_to_pheno_to_cons_to_var_list[gene] = {
                'pathogenic': get_cons_to_var_list(path_list),
                'benign': get_cons_to_var_list(ben_list)
                     }
        hgnc_data = hgnc_to_gene_data.get(gene)
        for path_var in path_list: 
            if hgnc_data:
                if hgnc_data.get('gnomAD'):
                    loeuf = hgnc_data['gnomAD'].get('loeuf')
                else:
                    continue
            else:
                continue
            if not loeuf: continue
            af = 0
            if path_var[2]:
                af = path_var[2].get('controlsAllAf', 0)
            if loeuf <= 1:
                path_below_1.append((gene, af, loeuf, path_var[5], path_var[6]))
            else:
                path_above_1.append((gene, af, loeuf, path_var[5], path_var[6]))
            path_vars.append((gene, af, loeuf, path_var[5], path_var[6]))
        for ben_var in ben_list: 
            if hgnc_data:
                if hgnc_data.get('gnomAD'):
                    loeuf = hgnc_data['gnomAD'].get('loeuf')
                    if not loeuf: continue
                else:
                    continue
            else:
                continue
            if ben_var[2]:
                af = ben_var[2].get('controlsAllAf')
                if not af: continue # Benign variants really ought to be in gnomad. 
            if loeuf <= 1:
                ben_below_1.append((gene, af, loeuf, ben_var[5], ben_var[6]))
            else:
                ben_above_1.append((gene, af, loeuf, ben_var[5], ben_var[6]))
            ben_vars.append((gene, af, loeuf, ben_var[5], ben_var[6]))
    
    with open(output_file, 'w') as o_file:
        for var in path_vars:
            o_file.write(f"{var[0]}\t{var[1]}\t{var[2]}\t{var[3]}\t{var[4]}\tpathogenic\n")
        for var in ben_vars:
            o_file.write(f"{var[0]}\t{var[1]}\t{var[2]}\t{var[3]}\t{var[4]}\tbenign\n")

    print("full paths")
    identify_loeuf_cutoff_values(path_vars, 'path_full.png')
    print("constrained paths")
    identify_loeuf_cutoff_values(path_below_1, 'path_constrained.png')
    print("tolerant paths")
    identify_loeuf_cutoff_values(path_above_1, 'path_tolerant.png')
    print("full bens")
    identify_loeuf_cutoff_values(ben_vars, 'ben_full.png')
    print("constrained bens")
    identify_loeuf_cutoff_values(ben_below_1, 'ben_constrained.png')
    print("tolerant bens")
    identify_loeuf_cutoff_values(ben_above_1, 'ben_tolerant.png')

    if skipped:
        logging.warning("Skipped %i variants! Review script, this should be 0 or very low.", skipped)

def main():
    """
    main function, calls other functions
    """
    options = parseArgs()
    logging.debug("extracting clinvar data into dictionary")
    clinvar_to_data = get_clinvar_var_to_data(options.in_vcf)
    logging.debug("extracting clinvar data into dictionary")
    extract_aa_information(options.clinvar_nirvana_json, clinvar_to_data, options.output_file)

if __name__ == "__main__":
    sys.exit(main())
