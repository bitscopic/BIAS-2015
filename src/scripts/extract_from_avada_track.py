import argparse
import csv

def clean_and_extract_columns(input_file, output_file):
    # Read the file, clean corrupted characters, and extract specified columns
    with open(input_file, 'rb') as f:
        content = f.read().decode('utf-8', errors='ignore')

    lines = content.splitlines()
    cleaned_data = [line.split('\t') for line in lines]

    # Extract the desired columns (1-3, 15)
    extracted_data = [[row[0], row[1], row[2], row[14]] for row in cleaned_data]

    # Write the cleaned and extracted data to the output file
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(extracted_data)

def main():
    parser = argparse.ArgumentParser(description="Clean a TSV file and extract specified columns.")
    parser.add_argument('input_file', type=str, help='Input TSV file with potential corrupted characters.')
    parser.add_argument('output_file', type=str, help='Output TSV file with cleaned and extracted data.')

    args = parser.parse_args()
    
    clean_and_extract_columns(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
