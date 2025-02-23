"""
Take a clinvar vcf nirvana json file and generate a tsv file with the pathogenic aa changes.

In the Richards et al. guidelines for variant classification, PS1 (Pathogenic Strong 1) is a criterion used when a
missense variant involves the same amino acid change as a known pathogenic variant, but with a different nucleotide
change. This suggests that the variant likely has the same functional impact as the previously established pathogenic
variant. 

PS1 is applied only when the original pathogenic variant's evidence is strong and there is no reason to doubt the new
variant's equivalence in effect.
"""
import json
import sys
import os
import io
import gzip
import argparse
import logging

# Add the parent directory of "scripts" to the Python path to make "pipeline_runner_core" importable
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from src.bias_2015 import extract_from_nirvana_json


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

    parser.set_defaults(verbose = "INFO")

    options = parser.parse_args()
    return options


def open_file(file_path, mode):
    """
    Open either a normal or a .gz file
    """
    _, file_extension = os.path.splitext(file_path)
    if file_extension == ".gz":
        return gzip.open(file_path, mode)
    return io.open(file_path, mode, encoding="utf-8")

# Parse the nirvana output into a more machine friendly format
def extract_aa_information(inJson, inVcf, output_file):
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
    clinvar_var_to_data = {}
    with open(inVcf) as in_file:
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
    skipped = 0 
    hgnc_to_gene_data = extract_from_nirvana_json.load_nirvana_gene_information(inJson) 
    with open(output_file, 'w') as o_file:
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
                    rs_id, criteria, signif, clinvar_id = clinvar_var_to_data[v_index]
                    # Don't process intergenic pathogenic variants - they won't have any AA change associated with them
                    if var.geneName == 'n/a': 
                        continue
                    # There is not an aa change associated with this variant
                    if var.protein_variant == "n/a": 
                        continue
                    values = [
                            var.geneName,
                            var.protein_variant,
                            rs_id,
                            criteria.replace("_", " "),
                            signif,
                            clinvar_id,
                            var.consequence
                           ]
                    if var.protein_variant:
                        o_file.write("\t".join(values) + "\n")
                count += 1
                if count % 10000 == 0:
                    logging.info("processed %i variants", count)
    if skipped:
        logging.warning("Skipped %i variants! Review script, this should be 0 or very low.", skipped)

def main():
    """
    main function, calls other functions
    """
    options = parseArgs()
    extract_aa_information(options.clinvar_nirvana_json, options.in_vcf, options.output_file)

if __name__ == "__main__":
    sys.exit(main())
