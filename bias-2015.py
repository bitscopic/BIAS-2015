"""
Take in an Illumina Connected Annotations (aka NIRVANA) .json file and apply ACMG standards using
the BIAS-2015 algorithm to classify the variants.
"""
import json
import sys
import os
import io
import gzip
import argparse
import time
from src.bias_2015 import extract_from_nirvana_json, variant_interpretation, variant_interpretation_dataset_loader
    
# Keep track of runtime
start_time = time.time()

def parseArgs(): 
    """
    Parse the command line arguments into useful python objects.  '--' variables are optional
    set_defaults only applies when the argument is not provided (it won't override)
    """
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument("nirvana_json_file",
                        help = " Variant call format (vcf) file",
                        action = "store")
    parser.add_argument("file_paths_json",
                        help = " Variant call format (vcf) file",
                        action = "store")
    parser.add_argument("output_file",
                        help = " Output file",
                        action = "store")
    parser.add_argument("--user_classifiers",
                        help = " User provided classifiers for each variant",
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


def run_bias(options):
    """
    Load user provided classifications if provided. Load in the data files needed for classification

    Stream through the ICA/NIRVANA json file, keeping only one json element in memory at a time. 

    Classify each variant one at a time, and write it to the output file as it is being classified. 
    """
    # Load the core paths
    config = {}
    with open(options.file_paths_json,'r') as config_file:
        config = json.load(config_file)

    # User provided classifications (optional)
    variant_to_user_classification = {}
    if options.user_classifiers:
        print("User classifications detected")
        with open(options.user_classifiers, 'r') as in_file:
            in_file.readline()
            for line in in_file:
                split_line = line.strip().split("\t")
                chrom, pos, ref, alt = split_line[:4]
                justification = json.loads(split_line[-1]) # The justification is the final column of the output
                variant = (chrom, pos, ref, alt)
                variant_to_user_classification[variant] = justification 

    # Load in the datasets needed for classification
    ###### NOTE: This step loads a lot of data into memory, it can take up to ~5 GB of RAM ########
    print("Loading data sets")
    name_to_dataset = variant_interpretation_dataset_loader.get_name_to_dataset(config)
    end_time = time.time()
    load_time = end_time - start_time
    print(f"Time to load datasets: {load_time:.4f} seconds")
    
    # Open the output file and stream through the input file, loading one variant into memory at a time,
    # classifying it, then writing the result to the output file and getting another variant.
    v_count = 0
    with open(options.output_file, 'w') as o_file:
        headers = ["chromosome", "position", "refAllele", "altAllele", "variantType", "consequence", "acmgClassification",
                   "alleleFreq", "hgvsg", "hgvsc", "hgvsp", "aaChange", "geneName", "pubmedIds", "associatedDiseases", "dbSnpids",
                   "transcript"]

        o_file.write("\t".join(headers) + "\n")
        # Open the input .json file
        with open_file(options.nirvana_json_file, "rt") as f:
            # gather the header line
            _ = f.readline()[10:-15]

            print("Classifying variants")
            genes_section = False
            hgnc_to_gene_data = {}
            # Read the file line by line
            for line in f:
                if genes_section:
                    # Check if we have reached the end of the genes section
                    if line.strip() == "]}," or line.strip() == "]}":
                        break

                    # Parse each line as JSON
                    if line[-2] == ",":
                        data = json.loads(line[:-2])
                    else:
                        data = json.loads(line)

                    # Extract gene information
                    hgnc_id = data["name"]
                    hgnc_to_gene_data[hgnc_id] = data
                    continue
                # Check if we have reached the start of the genes section
                if '"genes":[' in line:
                    genes_section = True
                    continue
                # Parse each line as JSON
                if line == ']}\n':
                    continue
                if line[-2] == ",":
                    try:
                        data = json.loads(line[:-2])
                    except:
                        data = json.loads(line[:-2].replace('""', '"'))
                else:
                    try:
                        data = json.loads(line)
                    except:
                        data = json.loads(line[:-2].replace('""', '"'))
                
                v_count += 1
                if v_count % 1000 == 0:
                    print(f"Processed {v_count} variants")

                # Extract position information
                position = data
                transcript_database = "RefSeq"
                variant_list = data.get("variants")
                for variant in variant_list:
                    chrom = position["chromosome"]
                    pos = str(position["position"])
                    single_variant = extract_from_nirvana_json.process_variant(variant, hgnc_to_gene_data, transcript_database, chrom, pos)
                    variant_key = (chrom, pos, single_variant.refAllele, single_variant.altAllele)
                    supplemental_codes = {}
                    if variant_to_user_classification.get(variant_key):
                        supplemental_codes = variant_to_user_classification[variant_key]
                    single_variant.significance, single_variant.justification = \
                            variant_interpretation.get_variant_classification(single_variant, name_to_dataset, supplemental_codes)
                    # Convert the justification to a JSON string
                    o_file.write(single_variant.to_tsv() + "\n")
    print(f"Processed {v_count} total variants")

def main():
    """
    Go through Nirvana json and perform calculations
    """
    # Parse arguments 
    options = parseArgs()
    
    # Run the BIAS algorithm
    run_bias(options)
    
    # Print out wall runtime
    end_time = time.time()
    total_runtime = end_time - start_time
    print(f"Total runtime: {total_runtime:.4f} seconds")

if __name__ == "__main__":
    sys.exit(main())
