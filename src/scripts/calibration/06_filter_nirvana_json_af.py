"""
Joins gnomAD allele frequency data from a Nirvana-annotated JSON onto the
filtered TSV from script 05. No re-filtering - script 05's decisions are final.

Outputs a TSV with gnomAD popmax AF, sex-specific AFs, AC, and AN added.
"""

import argparse
import json
import csv

parser = argparse.ArgumentParser()
parser.add_argument("--input_tsv", required=True, help="Path to filtered TSV from script 05")
parser.add_argument("--input_json", required=True, help="Path to Nirvana-annotated JSON")
parser.add_argument("--output", required=True, help="Path to output TSV with gnomAD AF")
args = parser.parse_args()

# Build gnomAD AF dict from Nirvana JSON: (chrom, pos, ref, alt) -> gnomAD fields
variant_to_gnomad_afs = {}

print("Building gnomAD lookup from Nirvana JSON...")

with open(args.input_json, 'r') as infile:
    next(infile)  # Skip Nirvana JSON header line
    for line in infile:
        line = line.strip().rstrip(',')

        if not line or line.startswith('],"genes":[') or line == ']}' or line.startswith('{"name":'):
            continue

        parsed_nirvana_variant_information = json.loads(line)
        if 'variants' not in parsed_nirvana_variant_information:
            continue

        chrom = parsed_nirvana_variant_information['chromosome']
        pos = str(parsed_nirvana_variant_information['position'])
        ref = parsed_nirvana_variant_information['refAllele']

        for variant in parsed_nirvana_variant_information['variants']:
            alt = variant['altAllele']

            popmax_af = 0
            male_popmax_af = None
            female_popmax_af = None
            max_an = 0
            all_ac = 0
            all_an = 0

            for gnomad_key in ('gnomad', 'gnomad-exome'):
                if gnomad_key not in variant:
                    continue
                if variant[gnomad_key].get('failedFilter', False):
                    continue
                an = variant[gnomad_key].get('allAn', 0)
                if an > max_an:
                    max_an = an
                    all_ac = variant[gnomad_key].get('allAc', 0)
                    all_an = an
                for key, value in variant[gnomad_key].items():
                    if 'Af' in key and key not in ('allAf', 'maleAf', 'femaleAf', 'controlsAllAf'):
                        if value > popmax_af:
                            popmax_af = value
                m = variant[gnomad_key].get('maleAf')
                f = variant[gnomad_key].get('femaleAf')
                if m is not None and (male_popmax_af is None or m > male_popmax_af):
                    male_popmax_af = m
                if f is not None and (female_popmax_af is None or f > female_popmax_af):
                    female_popmax_af = f

            variant_to_gnomad_afs[(chrom, pos, ref, alt)] = (popmax_af, male_popmax_af, female_popmax_af, all_ac, all_an)

print(f"gnomAD lookup built: {len(variant_to_gnomad_afs)} variants")

# Read script 05 TSV and join gnomAD AF
matched = 0
unmatched = 0
output_rows = []

with open(args.input_tsv, 'r') as infile:
    reader = csv.DictReader(infile, delimiter='\t')
    for row in reader:
        key = (row['chromosome'], row['position'], row['ref_allele'], row['alt_allele'])
        gnomad = variant_to_gnomad_afs.get(key)

        if gnomad:
            popmax_af, male_popmax_af, female_popmax_af, all_ac, all_an = gnomad
            matched += 1
        else:
            popmax_af, male_popmax_af, female_popmax_af, all_ac, all_an = 0, None, None, 0, 0
            unmatched += 1

        output_rows.append([
            row['chromosome'], row['position'], row['ref_allele'], row['alt_allele'],
            row['gene_symbol'], row['clinical_significance'], row['review_stars'],
            popmax_af, male_popmax_af, female_popmax_af, all_ac, all_an
        ])

print(f"\nMatched to gnomAD: {matched}")
print(f"No gnomAD match: {unmatched}")
print(f"Total output: {len(output_rows)}")

headers = ["chromosome", "position", "ref_allele", "alt_allele", "gene_symbol",
           "clinical_significance", "review_stars", "gnomad_popmax_af",
           "gnomad_male_popmax_af", "gnomad_female_popmax_af", "gnomad_all_ac", "gnomad_all_an"]

with open(args.output, 'w', newline="", encoding="utf-8") as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    writer.writerow(headers)
    writer.writerows(output_rows)
