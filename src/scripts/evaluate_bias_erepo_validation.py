"""
Compares algorithm classification results to the expected values. Currently supports comparing
BIAS-2015 and Intervar results to the evRepo dataset.

The evRepo data was downloaded and converted to a .vcf file, which was then processed using either
algorithm. The resulting output files can be used with this script to see how each algorithm performed.

ClinGen evRepo is a curated repository of evidence-level genomic data, where each entry is manually reviewed by experts and
classified according to ACMG guidelines, ensuring a consistent and accurate interpretation of variants. By integrating 
raw evidence, annotations, and metadata, evRepo serves as a source of truth for variant classification. Its rigorous
review process and adherence to ACMG standards make it an invaluable resource for bioinformatics workflows, variant
interpretation, and downstream analysis.

https://erepo.clinicalgenome.org/evrepo/

"""
import argparse
import json

def parseArgs():
    """
    Parse the command line arguments into useful python objects.  '--' variables are optional
    set_defaults only applies when the argument is not provided (it won't override)
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("erepo_file",
                        help="eRepo file containing variant information",
                        action="store")
    parser.add_argument("variant_file",
                        help="File containing Bitscopic/Intervar variant classifications",
                        action="store")
    parser.add_argument("vcf_file",
                        help="Variant call format (VCF) file",
                        action="store")
    parser.add_argument("algorithm",
                        help="Variant call format (VCF) file",
                        action="store",
                        choices=["bias", "intervar"])
    parser.add_argument("--sens_spec_output",
                        help="Output holding sensitivity and specificity values",
                        action="store")
    parser.add_argument("--concordance_output",
                        help="Output holding concordance values",
                        action="store")
    parser.add_argument("--full_code_table_output",
                        help="Output holding the full analytical results for all 28 codes",
                        action="store")
    parser.add_argument("--verbose",
                        help="The verbosity level for stdout messages (default INFO)",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        action="store")

    parser.set_defaults(sens_spec_output="sens_spec.tsv")
    parser.set_defaults(concordance_output="concordance.tsv")
    parser.set_defaults(verbose="INFO")
    options = parser.parse_args()

    return options

def get_variant_to_user_classification(user_file):
    """
    Parse the user file and map each variant to its classification and details.
    
    Arguments:
    user_file -- Path to the user file.
    
    Returns:
    variant_to_bit_class -- Dictionary mapping variants to their classification and details.
    """
    variant_to_bit_class = {}
    with open(user_file, 'r') as user_data:
        user_data.readline()
        for line in user_data:
            split_line = line.strip().split("\t")
            chromosome, position, reference, alternate = split_line[:4]
            classification = split_line[6]
            details = split_line[-1]
            if len(reference) > 1 and len(alternate) == 1:
                reference = reference[1:]
                alternate = "-"
            if len(alternate) > 1 and len(reference) == 1:
                alternate = alternate[1:]
                reference = "-"
            if len(alternate) > 1 and len(reference) > 1:
                if alternate[0] == reference[0]:
                    alternate = alternate[1:]
                    reference = reference[1:]
            variant = (chromosome, position, reference, alternate)
            variant_to_bit_class[variant] = {
                "classification": classification,
                "details": json.loads(details.replace("'", '"'))  # Convert string dict to actual dict
            }
    return variant_to_bit_class

def parse_intervar_data(data):
    """
    Take something like 

    InterVar: Uncertain significance PVS1=0 PS=[0, 0, 0, 0, 0] PM=[0, 0, 0, 0, 0, 0, 0] PP=[0, 0, 0, 0, 0, 0] BA1=0 
                BS=[0, 0, 0, 0, 0] BP=[0, 0, 0, 0, 0, 0, 0, 0] 

    and make something like

    {'pvs': {'pvs1': [0, '']}, 'ps': {'ps1': [0, ''], 'ps2': [0, ''], 'ps3': [0, ''], 'ps4': [0, ''], 'ps': [0, '']}, 
    'pm': {'pm1': [0, ''], 'pm2': [0, ''], 'pm3': [0, ''], 'pm4': [0, ''], 'pm5': [0, ''], 'pm6': [0, ''], 'pm': [0, '']},
    'pp': {'pp1': [0, ''], 'pp2': [0, ''], 'pp3': [0, ''], 'pp4': [0, ''], 'pp5': [0, ''], 'pp': [0, '']}, 'ba': 
    {'ba1': [0, '']}, 'bs': {'bs1': [0, ''], 'bs2': [0, ''], 'bs3': [0, ''], 'bs4': [0, ''], 'bs': [0, '']}, 'bp': 
    {'bp1': [0, ''], 'bp2': [0, ''], 'bp3': [0, ''], 'bp4': [0, ''], 'bp5': [0, ''], 'bp6': [0, ''], 'bp7': [0, ''], 
    'bp': [0, '']}}

    """
    # Split the input string into parts by whitespace
    parts = data.split()

    # Loop through the list and find the index
    index_of_pvs1 = -1
    for i, entry in enumerate(parts):
        if "PVS1" in entry:
            index_of_pvs1 = i
            break

    # Extract the classification string
    classification = " ".join(parts[1:index_of_pvs1])

    # Initialize the JSON-like structure
    parsed_data = {
        "pvs": {"pvs1": [0, ""]},
        "ps": {f"ps{i+1}": [0, ""] for i in range(5)},
        "pm": {f"pm{i+1}": [0, ""] for i in range(7)},
        "pp": {f"pp{i+1}": [0, ""] for i in range(6)},
        "ba": {"ba1": [0, ""]},
        "bs": {f"bs{i+1}": [0, ""] for i in range(5)},
        "bp": {f"bp{i+1}": [0, ""] for i in range(8)}
    }

    # Split the remaining data into key-value pairs
    cur_section = 'ps'
    count = 0
    for entry in parts[index_of_pvs1:]:
        if 'PVS1' in entry:
            key, value = entry.split("=")
            parsed_data['pvs']['pvs1'] = [int(value), ""]
            continue
        if 'BA1' in entry:
            key, value = entry.split("=")
            parsed_data['ba']['ba1'] = [int(value), ""]
            continue
        if "=" in entry:
            key, value = entry.replace("[", "").replace(",","").split("=")
            if key.lower() != cur_section:
                cur_section = key.lower()
                count = 1
            if cur_section == "ps":
                count += 1
            parsed_data[cur_section][f"{cur_section}{count}"] = [int(value), ""]
        else:
            count += 1
            value = int(entry.replace("]", "").replace(",",""))
            parsed_data[cur_section][f"{cur_section}{count}"] = [value, ""]

    if classification == "Uncertain significance":
        classification = "uncertain"
    return classification.lower(), parsed_data

def get_variant_to_intervar_classification(intervar_file):
    """
    Parse the InterVar file and map each variant to its classification and structured details.

    Arguments:
    intervar_file -- Path to the InterVar file.

    Returns:
    variant_to_intrvr_class -- Dictionary mapping variants to their classification and structured details.
    """
    variant_to_intrvr_class = {}
    with open(intervar_file, 'r') as intervar_data:
        header = intervar_data.readline().strip().split("\t")

        # Column indices for required fields
        chr_col = header.index("#Chr")
        start_col = header.index("Start")
        ref_col = header.index("Ref")
        alt_col = header.index("Alt")
        intervar_col = 13
        for line in intervar_data:
            split_line = line.strip().split("\t")

            # Extract variant fields
            chromosome = split_line[chr_col]
            position = split_line[start_col]
            reference = split_line[ref_col]
            alternate = split_line[alt_col]

            classification, details = parse_intervar_data(split_line[intervar_col])

            # Construct variant key
            variant = (f"chr{chromosome}", position, reference, alternate)
            # Add to result dictionary
            variant_to_intrvr_class[variant] = {
                "classification": classification,
                "details": details,
            }
    return variant_to_intrvr_class


def get_variant_to_erepo_line(vcf_file):
    """
    Parse the VCF file and map each variant to its corresponding eRepo line.
    
    Arguments:
    vcf_file -- Path to the VCF file.
    
    Returns:
    variant_to_erepo_line -- Dictionary mapping variants to eRepo line numbers.
    """
    variant_to_erepo_line = {}
    with open(vcf_file) as vcf_data:
        for line in vcf_data:
            if line.startswith('#'):
                continue
            chromosome, position, _, reference, alternate, _, _, info, _, _ = line.split('\t')
            info_dict = dict((item.split('=') for item in info.split(';') if '=' in item))
            if reference[0] == alternate[0]:
                if len(reference) > 1:
                    reference = reference[1:]
                else:
                    reference = '-'
                if len(alternate) > 1:
                    alternate = alternate[1:]
                else:
                    alternate = "-"

            variant = (chromosome, position, reference, alternate)
            erepo_line = int(info_dict['EreppoLine'])
            variant_to_erepo_line[variant] = erepo_line
            if len(reference) > 1 or len(alternate) > 1:
                index = 0
                # Ensure bounds check happens before character comparison
                while index < len(reference) and index < len(alternate) and reference[index] == alternate[index]:
                    index += 1
                if index > 0:
                    # Adjust reference and alternate based on the common prefix
                    new_reference = reference[index:] if reference[index:] else "-"
                    new_alternate = alternate[index:] if alternate[index:] else "-"
                    
                    # Create a new variant and add it to the dictionary
                    new_variant = (chromosome, str(int(position) + index), new_reference, new_alternate)
                    variant_to_erepo_line[new_variant] = erepo_line

    return variant_to_erepo_line

def get_erepo_line_to_data(erepo_file):
    """
    Parse the eRepo file and map each line number to the data for that line.
    
    Arguments:
    erepo_file -- Path to the eRepo file.
    
    Returns:
    erepo_line_to_data -- Dictionary mapping line numbers to eRepo data.
    """
    erepo_line_to_data = {}
    with open(erepo_file) as erepo_data:
        headers = next(erepo_data).split('\t')
        for line_number, line in enumerate(erepo_data, start=1):
            values = line.split('\t')
            erepo_line_to_data[line_number] = dict(zip(headers, values))
    return erepo_line_to_data

class_to_score = {
        'benign': 2,
        'likely benign': 1,
        'uncertain': 0,
        'uncertain significance': 0,
        'likely pathogenic': -1,
        'pathogenic': -2
        }

# Calculate Sensitivity, Specificity, Precision, and F1 Score
def calculate_metrics(confusion_matrix):
    """
    Calculate sensitivity, specificity, precision, and F1 score.
    """
    TP = confusion_matrix['TP']
    FP = confusion_matrix['FP']
    TN = confusion_matrix['TN']
    FN = confusion_matrix['FN']
    
    sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0
    specificity = TN / (TN + FP) if (TN + FP) > 0 else 0
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    f1_score = (2 * precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) > 0 else 0
    
    return sensitivity, specificity, precision, f1_score


def compare_variant_set_to_erepo(variant_to_bit_class, variant_to_erepo_line, erepo_line_to_data, sens_spec_out, conc_out, full_table_out, algorithm):
    """
    Takes a set of properly formatted variants and compares it to the eRepo data.
    Adds concordance calculations for ACMG codes (e.g., PVS1, PP2, etc.).
    """
    # Initialize confusion matrices
    confusion_matrices = {
        'pathogenic': {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0},
        'benign': {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0},
        'uncertain': {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    }
    user_code_to_count = {}
    acmg_code_to_count = {}
    concordance_counts = {code: {'agree': 0, 'alg': 0, 'evrepo':0} for code in [
        'PVS1', 'PS1', 'PS2', 'PS3', 'PS4', 
        'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
        'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
        'BA1', 'BS1', 'BS2', 'BS3', 'BS4', 'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'
    ]}
    full_conc_counts = {}
    compared_variant_count = 0

    my_set = set()

    for variant, bit_class in variant_to_bit_class.items():
        erepo_line = variant_to_erepo_line.get(variant)
        if not erepo_line:  # Handle potential variant representation mismatches
            new_pos = int(variant[1]) - 1
            new_var = (variant[0], str(new_pos), variant[2], variant[3])
            erepo_line = variant_to_erepo_line.get(new_var)
            if not erepo_line:
                new_pos = int(variant[1]) + 1
                new_var = (variant[0], str(new_pos), variant[2], variant[3])
                erepo_line = variant_to_erepo_line.get(new_var)

        if erepo_line:
            erepo_data = erepo_line_to_data.get(erepo_line)
            if erepo_data:
                compared_variant_count += 1
                bit_score = class_to_score[bit_class['classification']]
                erepo_score = class_to_score[erepo_data['Assertion'].lower()]

                # Track user codes
                alg_acmg_codes = []
                full_alg_acmg_codes = []
                for _, code_to_evidence in bit_class['details'].items():
                    for code, evidence in code_to_evidence.items():
                        code_str = code.upper()
                        code_score, code_explain = evidence
                        if code_score > 0:
                            alg_acmg_codes.append(code_str)
                            user_code_to_count[code_str] = user_code_to_count.get(code_str, 0) + 1
                            full_code = code_explain.split(":")[0] # Full acmg code, ex BP3_Strong
                            full_alg_acmg_codes.append(full_code.lower())

                # Track ACMG codes from eRepo
                raw_acmg_codes = erepo_data['Applied Evidence Codes (Met)'].replace(" ", "").split(",")
                true_acmg_codes = []
                full_acmg_codes = []
                for acmg_code in raw_acmg_codes:
                    my_set.add(acmg_code)
                    true_code = acmg_code.split("_")[0] if "_" in acmg_code else acmg_code
                    if not true_code: 
                        continue
                    full_acmg_codes.append(acmg_code.lower().replace("(", "").replace(")", ""))
                    true_acmg_codes.append(true_code)
                    acmg_code_to_count[true_code] = acmg_code_to_count.get(true_code, 0) + 1

                # Calculate concordance for the codes
                best_agree = 0
                for code in full_alg_acmg_codes:
                    if code in full_acmg_codes:
                        best_agree += 1
                        if not full_conc_counts.get(code):
                            full_conc_counts[code] = {'agree': 1, 'alg': 0, 'evrepo': 0}
                        else:
                            full_conc_counts[code]['agree'] += 1
                    else:
                        if not full_conc_counts.get(code):
                            full_conc_counts[code] = {'agree': 0, 'alg': 1, 'evrepo': 0}
                        else:
                            full_conc_counts[code]['alg'] += 1
                for code in full_acmg_codes:
                    if code in full_alg_acmg_codes: # already evaluated
                        continue
                    if not full_conc_counts.get(code):
                        full_conc_counts[code] = {'agree': 0, 'alg': 0, 'evrepo': 1}
                    else:
                        full_conc_counts[code]['evrepo'] += 1

                # Calculate concordance for the codes
                best_agree = 0
                for code in alg_acmg_codes:
                    if code in true_acmg_codes:
                        best_agree += 1
                        #if code == "PS1": print(variant, bit_class)
                        concordance_counts[code]['agree'] += 1
                    else:
                        concordance_counts[code]['alg'] += 1
                for code in true_acmg_codes:
                    #if code == "PS1": print(variant, erepo_data, bit_class)
                    if code in alg_acmg_codes: # already evaluated
                        continue
                    concordance_counts[code]['evrepo'] += 1

                # Update confusion matrices
                if erepo_score < 0:  # True Pathogenic
                    if bit_score < 0:
                        confusion_matrices['pathogenic']['TP'] += 1
                    else:
                        confusion_matrices['pathogenic']['FN'] += 1
                else:  # Not Pathogenic
                    if bit_score < 0:
                        confusion_matrices['pathogenic']['FP'] += 1
                    else:
                        confusion_matrices['pathogenic']['TN'] += 1

                if erepo_score > 0:  # True Benign
                    if bit_score > 0:
                        confusion_matrices['benign']['TP'] += 1
                    else:
                        confusion_matrices['benign']['FN'] += 1
                else:  # Not Benign
                    if bit_score > 0:
                        confusion_matrices['benign']['FP'] += 1
                    else:
                        confusion_matrices['benign']['TN'] += 1

                if erepo_score == 0:  # True Uncertain
                    if bit_score == 0:
                        confusion_matrices['uncertain']['TP'] += 1
                    else:
                        confusion_matrices['uncertain']['FN'] += 1
                else:  # Not Uncertain
                    if bit_score == 0:
                        confusion_matrices['uncertain']['FP'] += 1
                    else:
                        confusion_matrices['uncertain']['TN'] += 1
        else:
            print(variant)

    # Collect concordance values for codes with F1 scores
    concordance_list = []

    for code, counts in full_conc_counts.items():
        TP = counts['agree']
        FP = counts['alg']
        FN = counts['evrepo']
        
        # Compute metrics
        precision = TP / (TP + FP) if (TP + FP) > 0 else 0
        recall = TP / (TP + FN) if (TP + FN) > 0 else 0
        f1_score = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        concordance = TP / (TP + FP + FN) if (TP + FP + FN) > 0 else 0

        # Store results
        concordance_list.append((code, TP, FP, FN, concordance, precision, recall, f1_score))

    # Sort by code
    concordance_list.sort(key=lambda x: x[0])

    # Print results
    with open(full_table_out, 'w') as o_file:
        o_file.write("ACMG_code\ttrue_positive\tfalse_positive\tfalse_negative\tconcordance\tprecision\trecall\tf1\n")
        for code, TP, FP, FN, concordance, precision, recall, f1_score in concordance_list:
            print(f"{code}: TP={TP}, FP={FP}, FN={FN}, Concordance={concordance:.4f}, Precision={precision:.4f}, Recall={recall:.4f}, F1={f1_score:.4f}")
            o_file.write(f"{code}\t{TP}\t{FP}\t{FN}\t{concordance:.4f}\t{precision:.4f}\t{recall:.4f}\t{f1_score:.4f}\n")

    total_correct, total_missed, total_wrong, total_f1, codes_evaluated = 0, 0, 0, 0, 0
    print("Concordance and F1 scores for overlapping codes:")
    total_f1 = 0
    with open(conc_out, 'w') as out_file:
        out_file.write("algorithm\tcode\talg_fp\talg_fn\talg_correct\tconcordance\tprecision\trecall\tf1_score\n")
        for code, counts in concordance_counts.items():
            TP = counts['agree']
            FP = counts['alg']
            FN = counts['evrepo']
            
            # Compute metrics
            precision = TP / (TP + FP) if (TP + FP) > 0 else 0
            recall = TP / (TP + FN) if (TP + FN) > 0 else 0
            f1_score = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
            total_f1 += f1_score
            concordance = TP / (TP + FP + FN) if (TP + FP + FN) > 0 else 0

            total_correct += TP
            total_missed += FN
            total_wrong += FP
            codes_evaluated += 1

            # Write to output file
            out_file.write(f"{algorithm}\t{code}\t{FP}\t{FN}\t{TP}\t{concordance:.4f}\t{precision:.4f}\t{recall:.4f}\t{f1_score:.4f}\n")

            # Print results
            print(f"{code}: TP={TP}, FP={FP}, FN={FN}, Concordance={concordance:.4f}, Precision={precision:.4f}, Recall={recall:.4f}, F1={f1_score:.4f}")

    # Compute total concordance
    tot_conc = total_correct / (total_correct + total_missed + total_wrong)
    print(f"When evaluating {codes_evaluated} codes: Called {total_correct} codes correctly, missed {total_missed}, " +\
            f"called {total_wrong} incorrectly, total concordance {tot_conc:.4f}, total F1 {total_f1:.4f}")

    print(f"Compared {compared_variant_count} variants")
    print("All user codes")
    print(dict(sorted(user_code_to_count.items())))
    print("All ACMG codes")
    print(dict(sorted(acmg_code_to_count.items())))

    # Print sensitivity and specificity
    with open(sens_spec_out, 'w') as out_file:
        out_file.write("algorithm\tcategory\tsensitivity\tspecificity\tprecision\tf1_score\n")
        for category, matrix in confusion_matrices.items():
            sensitivity, specificity, precision, f1_score = calculate_metrics(matrix)
            print(f"{category.capitalize()} Sensitivity: {sensitivity:.4f}")
            print(f"{category.capitalize()} Specificity: {specificity:.4f}")
            out_file.write(f"{algorithm}\t{category}\t{sensitivity}\t{specificity}\t{precision}\t{f1_score}\n")


def main():
    """
    Main function to process variant data and map it to eRepo data using provided files.
    
    Workflow:
    - Parse command line arguments and configuration.
    - Load variant to user classification mappings.
    - Load variant to eRepo line mappings.
    - Load eRepo line to data mappings.
    - Combine data and perform necessary operations.
    """
    options = parseArgs()

    variant_file = options.variant_file
    vcf_file = options.vcf_file
    erepo_file = options.erepo_file

    variant_to_erepo_line = get_variant_to_erepo_line(vcf_file)
    erepo_line_to_data = get_erepo_line_to_data(erepo_file)

    variant_to_acmg_class = {}
    if options.algorithm == "bias":
        variant_to_acmg_class = get_variant_to_user_classification(variant_file)
    elif options.algorithm == "intervar":
        variant_to_acmg_class = get_variant_to_intervar_classification(variant_file)
   
    compare_variant_set_to_erepo(variant_to_acmg_class, variant_to_erepo_line, erepo_line_to_data, options.sens_spec_output,
                                 options.concordance_output, options.full_code_table_output, options.algorithm)

if __name__ == "__main__":
    main()
