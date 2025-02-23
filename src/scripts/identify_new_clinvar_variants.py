#!/usr/bin/env python3

"""
Helper script that identifies new variants when comparing two ClinVar datasets.

Usage: python3 identify_new_clinvar_variants.py old_clinvar.vcf new_clinvar.vcf output.vcf
"""

import argparse

def parse_arguments():
    """
    Parse command-line arguments for input and output VCF files.
    """
    parser = argparse.ArgumentParser(description='Identify new variants between two ClinVar VCF files.')
    parser.add_argument('old_vcf', help='Path to the old ClinVar VCF file.')
    parser.add_argument('new_vcf', help='Path to the new ClinVar VCF file.')
    parser.add_argument('output_vcf', help='Path to the output VCF file for new variants.')
    return parser.parse_args()

def read_vcf(file_path):
    """
    Read a VCF file and return its header lines and variant entries.
    """
    header_lines = []
    variants = []
    with open(file_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                variants.append(line.strip())
    return header_lines, variants

def extract_variant_info(variant_line):
    """
    Extract (chrom, pos, rs_id, ref, alt) from a VCF variant line.
    """
    fields = variant_line.split('\t')
    chrom = fields[0]
    pos = fields[1]
    rs_id = fields[2]
    ref = fields[3]
    alt = fields[4]
    return chrom, pos, rs_id, ref, alt

def build_variant_sets(variants):
    """
    Build sets of (chrom, pos, ref, alt) tuples and rsIDs from variant entries.
    """
    variant_set = set()
    rsid_set = set()
    for variant in variants:
        chrom, pos, rs_id, ref, alt = extract_variant_info(variant)
        variant_set.add((chrom, pos, ref, alt))
        if rs_id != '.':
            rsid_set.add(rs_id)
    return variant_set, rsid_set

def identify_new_variants(old_variants, new_variants):
    """
    Identify variants present in new_variants but not in old_variants.
    """
    old_variant_set, old_rsid_set = build_variant_sets(old_variants)
    new_variant_set, new_rsid_set = build_variant_sets(new_variants)

    new_unique_variants = []
    for variant in new_variants:
        chrom, pos, rs_id, ref, alt = extract_variant_info(variant)
        if ((chrom, pos, ref, alt) not in old_variant_set and
            (rs_id == '.' or rs_id not in old_rsid_set)):
            new_unique_variants.append(variant)
    return new_unique_variants

def write_vcf(header_lines, variants, output_path):
    """
    Write header lines and variant entries to an output VCF file.
    """
    with open(output_path, 'w') as output_file:
        for header in header_lines:
            output_file.write(header)
        for variant in variants:
            output_file.write(variant + '\n')

def main():
    """
    Main function to identify and write new variants between two VCF files.
    """
    args = parse_arguments()
    old_header, old_variants = read_vcf(args.old_vcf)
    new_header, new_variants = read_vcf(args.new_vcf)

    # Ensure the new VCF header is used for the output
    output_header = new_header

    new_unique_variants = identify_new_variants(old_variants, new_variants)
    write_vcf(output_header, new_unique_variants, args.output_vcf)

if __name__ == '__main__':
    main()
