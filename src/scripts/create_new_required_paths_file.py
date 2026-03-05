"""
create_new_required_paths_file.py

This script checks for the presence of specific required files in a given directory.
Users must specify the genome reference build (hg19 or hg38) and the annotator (nirvana or vep).
If all required files are found, it generates a JSON file mapping file types to their absolute paths.
If any files are missing, it reports them and does not generate the JSON output.

Usage:
    python3 create_new_required_paths_file.py <input_dir> <ref_b> <output_json> --annotator <nirvana|vep>
"""

import sys
import os
import json
import argparse
import logging


def parse_args():
    """
    Parse the command line arguments into useful Python objects.

    '--' variables are optional.
    'set_defaults' only applies when the argument is not provided (it won't override).
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir",
                        help="Directory where the files are expected to be found.",
                        action="store")
    parser.add_argument("ref_b",
                        choices=["hg19", "hg38"],
                        help="Reference genome build (hg19 or hg38).",
                        action="store")
    parser.add_argument("output_json",
                        help="Path for the output JSON file.",
                        action="store")
    parser.add_argument("--annotator",
                        choices=["nirvana", "vep"],
                        default="nirvana",
                        help="Annotation tool used during preprocessing (default: nirvana).",
                        action="store")
    parser.add_argument("--verbose",
                        help="The verbosity level for stdout messages (default WARNING).",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        action="store")

    parser.set_defaults(verbose="WARNING")
    options = parser.parse_args()
    logging.basicConfig(level=getattr(logging, options.verbose), format='%(message)s')
    return options


def check_required_files(input_dir, ref_b, annotator):
    """
    Check for the presence of required files in the given directory.

    Some files are annotator-specific (e.g. PS1_PM5 clinvar pathogenic aa data
    differs between nirvana and vep). Other files are shared across annotators.

    Returns a tuple of (found_files, missing_files).
    """
    output_dir = os.path.abspath(input_dir)

    # Files shared across both annotators
    files = {
        "PVS1_ncbi_ref_seq_hgmd_fp": f"{ref_b}_PVS1_ncbiRefSeqHgmd.tsv",
        "PVS1_PP3_BP4_BP7_splice_fp": f"{ref_b}_PVS1_PP3_BP4_BP7_splice_data.tsv",
        "PS1_PM5_clinvar_pathogenic_aa_fp": f"{ref_b}_PS1_PM5_clinvar_pathogenic_aa_{annotator}.tsv",
        "PS3_literature_gene_aa_fp": f"{ref_b}_PS3_lit_gene_aa.tsv",
        "PS3_literature_variant_fp": f"{ref_b}_PS3_lit_variant.tsv",
        "PS4_gwas_dbsnp_fp": f"{ref_b}_PS4_gwasCatalog.txt",
        "PM1_chrom_to_pathogenic_domain_list_fp": f"{ref_b}_PM1_chrom_to_pathogenic_domain_list.tsv",
        "PM4_BP3_coding_repeat_region_fp": f"{ref_b}_PM4_BP3_coding_repeat_regions.tsv",
        "PP2_missense_pathogenic_gene_to_region_list_fp": f"{ref_b}_PP2_missense_pathogenic_genes.tsv",
        "BP1_truncating_gene_to_data_fp": f"{ref_b}_BP1_truncating_genes.tsv",
        "BS2_clingen_gene_disease_validity_fp": f"{ref_b}_clingen_gene_disease_validity.csv",
        "PVS1_PM2_BA1_BS1_BS2_gnomad_gene_constraints_fp": f"{ref_b}_gnomad_gene_constraints.txt",
    }

    # VEP has an additional file for PS4 clinvar submitter counts
    if annotator == "vep":
        files["PS4_clinvar_submitter_counts_fp"] = f"{ref_b}_PS4_clinvar_submitter_counts.tsv"

    found_files = {}
    missing_files = []

    for key, filename in files.items():
        full_path = os.path.join(output_dir, filename)
        if os.path.exists(full_path):
            found_files[key] = full_path
        else:
            missing_files.append(full_path)

    return found_files, missing_files


def generate_output_json(file_dict, output_json):
    """
    Write the dictionary of found file paths to a JSON file.
    """
    with open(output_json, "w") as f:
        json.dump(file_dict, f, indent=4)


def main():
    """
    Parse arguments, check for required files, and handle output accordingly.
    If any required files are missing, they are reported, and no output is generated.
    If all required files are present, a JSON file with their paths is created.
    """
    options = parse_args()

    found_files, missing_files = check_required_files(options.input_dir, options.ref_b, options.annotator)

    if missing_files:
        logging.error("Error: The following required files are missing:")
        for missing in missing_files:
            logging.error(" - %s", missing)
        sys.exit(1)

    generate_output_json(found_files, options.output_json)
    logging.info("Output JSON written to: %s", options.output_json)


if __name__ == "__main__":
    sys.exit(main())
