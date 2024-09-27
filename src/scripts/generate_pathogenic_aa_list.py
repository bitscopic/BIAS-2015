"""
Take a clinvar vcf nirvana json file and generate a tsv file with the pathogenic
aa changes
"""
import json
import sys
import os
import io
import gzip
import argparse
# Add the parent directory of "scripts" to the Python path to make "pipeline_runner_core" importable
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from src.bias_2015 import extract_from_nirvana_json

# Map ClinVar reviewStatus to an integer value
clinvar_review_status_to_level = {
    "practice guideline": 5,
    "reviewed by expert panel": 5,
    "criteria provided, multiple submitters, no conflicts": 4,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 2,
    "no assertion criteria provided": 0,
    "no assertion provided": 0
}

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
def parse_nirvana_json(inJson, inVcf, output_file):
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
            ref = split_line[3]
            alt = split_line[4]
            if ref[0] == alt[0]:
                if len(ref) > 1:
                    ref = ref[1:]
                else:
                    ref = "-"
                if len(alt) > 1:
                    alt = alt[1:]
                else:
                    alt = "-"
            variant = (chrom, pos, ref, alt)
            term_to_value = {}
            elements = split_line[7].split(";")
            for pair in elements:
                term, value = pair.split("=")
                term_to_value[term] = value
            rs_id = term_to_value.get('RS', '')
            criteria = term_to_value.get('CLNREVSTAT', '')
            signif = term_to_value.get('CLNSIG', '')
            if not rs_id: continue
            if not criteria: continue
            if not signif: continue
            clinvar_var_to_data[variant] = (rs_id, criteria, signif)

    with open(output_file, 'w') as o_file:
        with open_file(inJson, "rt") as f:
            # gather the header line
            _ = f.readline()[10:-15]
            # Read the file line by line
            print("reading variant data")
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
                for variant in variant_list:
                    chrom = position["chromosome"]
                    pos = str(position["position"])
                    var = extract_from_nirvana_json.process_variant(variant, {}, "RefSeq", chrom, pos)
                    v_index = (var.chromosome, var.position, var.refAllele, var.altAllele) 
                    rs_id, criteria, signif = clinvar_var_to_data.get(v_index, ("", "", ""))
                    # Don't process intergenic pathogenic variants - they won't have any AA change associated with them
                    if var.geneName == 'n/a': continue
                    # The variant does not have one review status with a score > 1 (see scoring table)
                    if clinvar_review_status_to_level.get(criteria.lower().replace("_", " "), 0) < 2: continue

                    # There is not an aa change associated with this variant
                    if var.protein_variant == "n/a": continue
                    print(var.geneName)
                    values = [
                            var.geneName,
                            var.protein_variant,
                            rs_id,
                            criteria.replace("_", " "),
                            signif
                           ]
                    if var.protein_variant:
                        o_file.write("\t".join(values) + "\n")
                count += 1
                if count % 10000 == 0:
                    print(f"{count} done")

def main():
    """
    main function, calls other functions
    """
    options = parseArgs()
    parse_nirvana_json(options.clinvar_nirvana_json, options.in_vcf, options.output_file)

if __name__ == "__main__":
    sys.exit(main())
