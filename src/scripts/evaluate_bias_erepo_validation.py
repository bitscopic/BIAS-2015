"""
Runs some evaluation code
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
    parser.add_argument("bitscopic_file",
                        help="File containing Bitscopic variant classifications",
                        action="store")
    parser.add_argument("vcf_file",
                        help="Variant call format (VCF) file",
                        action="store")
    parser.add_argument("--verbose",
                        help="The verbosity level for stdout messages (default INFO)",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        action="store")

    parser.set_defaults(verbose="INFO")
    options = parser.parse_args()

    return options

def get_variant_to_bitscopic_classification(bitscopic_file):
    """
    Parse the bitscopic file and map each variant to its classification and details.
    
    Arguments:
    bitscopic_file -- Path to the bitscopic file.
    
    Returns:
    variant_to_bit_class -- Dictionary mapping variants to their classification and details.
    """
    variant_to_bit_class = {}
    with open(bitscopic_file, 'r') as bitscopic_data:
        bitscopic_data.readline()
        for line in bitscopic_data:
            split_line = line.strip().split("\t")
            chromosome, position, reference, alternate = split_line[:4]
            classification = split_line[6]
            details = split_line[-1]
            variant = (chromosome, position, reference, alternate)
            variant_to_bit_class[variant] = {
                "classification": classification,
                "details": json.loads(details.replace("'", '"'))  # Convert string dict to actual dict
            }
    return variant_to_bit_class

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


# Calculate Sensitivity and Specificity
def calculate_sensitivity_specificity(confusion_matrix):
    """
    Calculate sensitivity and specificity
    """
    sensitivity = confusion_matrix['TP'] / (confusion_matrix['TP'] + confusion_matrix['FN']) if (confusion_matrix['TP'] + confusion_matrix['FN']) > 0 else 0
    specificity = confusion_matrix['TN'] / (confusion_matrix['TN'] + confusion_matrix['FP']) if (confusion_matrix['TN'] + confusion_matrix['FP']) > 0 else 0
    return sensitivity, specificity


def main():
    """
    Main function to process variant data and map it to eRepo data using provided files.
    
    Workflow:
    - Parse command line arguments and configuration.
    - Load variant to bitscopic classification mappings.
    - Load variant to eRepo line mappings.
    - Load eRepo line to data mappings.
    - Combine data and perform necessary operations.
    """
    options = parseArgs()

    bitscopic_file = options.bitscopic_file
    vcf_file = options.vcf_file
    erepo_file = options.erepo_file

    variant_to_bit_class = get_variant_to_bitscopic_classification(bitscopic_file)
    variant_to_erepo_line = get_variant_to_erepo_line(vcf_file)
    erepo_line_to_data = get_erepo_line_to_data(erepo_file)

    # Initialize confusion matrices
    confusion_matrices = {
        'pathogenic': {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0},
        'benign': {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0},
        'uncertain': {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    }

    bitscopic_code_to_count = {}
    acmg_code_to_count = {}
    acmg_code_set = set()
    compared_variant_count = 0
    for variant, bit_class in variant_to_bit_class.items():
        erepo_line = variant_to_erepo_line.get(variant)
        if erepo_line:
            erepo_data = erepo_line_to_data.get(erepo_line)
            if erepo_data:
                compared_variant_count += 1
                bit_score = class_to_score[bit_class['classification']]
                erepo_score = class_to_score[erepo_data['Assertion'].lower()]
                for code_class, code_list in bit_class['details'].items():
                    count = 0
                    for code_tup in code_list:
                        count += 1
                        score, _ = code_tup
                        if score > 0:
                            code_str = f"{code_class.upper()}{count}"
                            if bitscopic_code_to_count.get(code_str):
                                bitscopic_code_to_count[code_str] += 1
                            else:
                                bitscopic_code_to_count[code_str] = 1
                true_acmg_codes = erepo_data['Applied Evidence Codes (Met)'].replace(" ", "").split(",")
                for acmg_code in true_acmg_codes:
                    acmg_code_set.add(acmg_code)
                    true_code = ""
                    if "_" in acmg_code:
                        true_code = acmg_code.split("_")[0]
                    else:
                        true_code = acmg_code
                    if acmg_code_to_count.get(true_code):
                        acmg_code_to_count[true_code] += 1
                    else:
                        acmg_code_to_count[true_code] = 1
                # Pathogenic confusion matrix
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

                # Benign confusion matrix
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

                # Uncertain confusion matrix
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
    
    print(f"Compared {compared_variant_count} variants")
    print("All bitscopic codes")
    sorted_data = dict(sorted(bitscopic_code_to_count.items()))
    print(sorted_data)
    print(len(sorted_data))
    print("All acmg codes")
    del acmg_code_to_count[''] # This clutters output
    sorted_data = dict(sorted(acmg_code_to_count.items()))
    print(sorted_data)
    print(len(sorted_data))

    # Calculate sensitivity and specificity for each category
    sens_spec = {}
    for category, matrix in confusion_matrices.items():
        sens_spec[category] = calculate_sensitivity_specificity(matrix)

    # Print results
    for category, (sensitivity, specificity) in sens_spec.items():
        print(f"{category.capitalize()} Sensitivity: {sensitivity:.4f}")
        print(f"{category.capitalize()} Specificity: {specificity:.4f}")

if __name__ == "__main__":
    main()
