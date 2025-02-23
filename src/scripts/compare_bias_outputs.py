'''
This script compares two TSV files with the same headers and outputs the differences in the ACMG classification and rationale columns.
Usage: python compare_bias.py <file1> <file2> > output.txt
This is to be used for validation only.

Steps to use:
1) Pull two different versions of BIAS.
2) Run the BIAS pipeline on the same VCF file using both versions.
3) Run the compare_bias.py script on the two output files.

Recommended to use with run_full_bias_eval.sh to get the outputs. 
Change the bias script locations in that .sh script to the correct locations of the two versions of BIAS, one for each run.
'''


import csv
import sys
import json

def read_tsv(file_path):
    data = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        headers = next(reader)
        for row in reader:
            data.append(row)
    return headers, data

def compare_files(file1, file2):
    csv.field_size_limit(sys.maxsize)
    headers1, data1 = read_tsv(file1)    
    headers2, data2 = read_tsv(file2)

    if headers1 != headers2:
        print("Files have different headers.")
        return

    acmg_index = headers1.index('acmgClassification')
    chr_index = headers1.index('chromosome')
    pos_index = headers1.index('position')
    ref_index = headers1.index('refAllele')
    alt_index = headers1.index('altAllele')
    rationale_index = len(headers1) - 1

    differences = []
    data2_dict = {(row[chr_index], row[pos_index], row[ref_index], row[alt_index]): row for row in data2}

    for row1 in data1:
        key = (row1[chr_index], row1[pos_index], row1[ref_index], row1[alt_index])
        if key in data2_dict:
            row2 = data2_dict[key]
            if row1[acmg_index] != row2[acmg_index]:                
                rationale1 = row1[rationale_index] if rationale_index < len(row1) else ""
                rationale2 = row2[rationale_index] if rationale_index < len(row2) else ""
                rationale_diff = extract_rationale_diff(rationale1, rationale2)
                differences.append((row1[chr_index], row1[pos_index], row1[ref_index], row1[alt_index], row1[acmg_index], row2[acmg_index], rationale_diff))

    return differences

def extract_rationale_diff(rationale1, rationale2):
    try:
        rationale_dict1 = json.loads(rationale1)
        rationale_dict2 = json.loads(rationale2)
    except json.JSONDecodeError:
        return f"Invalid JSON in rationale: {rationale1} vs {rationale2}"

    def compare_dicts(dict1, dict2):
        diff_str = ""
        for key in set(dict1.keys()).union(dict2.keys()):
            value1 = dict1.get(key, "")
            value2 = dict2.get(key, "")
            if isinstance(value1, dict) and isinstance(value2, dict):
                nested_diff_str = compare_dicts(value1, value2)
                if nested_diff_str:
                    diff_str += f"{key}: {nested_diff_str},"
            elif value1 != value2:
                diff_str += f"{key}: {value1} vs {value2};"
        return diff_str


    differences = compare_dicts(rationale_dict1, rationale_dict2)
    return differences

def main(file1, file2):
    differences = compare_files(file1, file2)
    if differences:
        print(f"Chromosome\tPosition\tACMGclassification ({file1})\tACMGclassification ({file2})\tRationale Differences")
        for diff in differences:
            print("\t".join(map(str, diff)))
    else:
        print("No differences found in ACMGclassification.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare_bias.py <file1> <file2>")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    main(file1, file2)