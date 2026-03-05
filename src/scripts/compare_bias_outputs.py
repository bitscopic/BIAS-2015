#!/usr/bin/env python3
"""
Compare two BIAS output files and report differences.

This script compares variant classifications and ACMG codes between two BIAS output files
(e.g., from different annotators like Nirvana vs VEP, or different BIAS versions).

Outputs:
  - High-level summary to stdout: classification counts, concordance matrix
  - Detailed per-variant comparison to output file

Examples:
    python3 src/scripts/compare_bias_outputs.py nirvana_output.tsv vep_output.tsv comparison_results.tsv
    python3 src/scripts/compare_bias_outputs.py file1.tsv file2.tsv results.tsv --label1 Nirvana --label2 VEP
"""
import argparse
from collections import defaultdict

from bias_output_utils import (
    CLASSIFICATIONS,
    CLASS_ABBREV,
    parse_bias_output,
    parse_rationale,
    get_active_codes,
    compare_codes,
    match_variants,
)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("file1", help="First BIAS output file")
    parser.add_argument("file2", help="Second BIAS output file")
    parser.add_argument("output_file", help="Output file for detailed variant comparison")
    parser.add_argument("--label1", default="File1", help="Label for first file (default: File1)")
    parser.add_argument("--label2", default="File2", help="Label for second file (default: File2)")
    return parser.parse_args()


def print_summary(variants1, variants2, matched, only1, only2, label1, label2):
    """Print high-level summary to stdout."""
    print("=" * 80)
    print(f"BIAS Output Comparison: {label1} vs {label2}")
    print("=" * 80)

    # Count classifications
    counts1 = defaultdict(int)
    counts2 = defaultdict(int)

    for v in variants1.values():
        counts1[v['classification']] += 1
    for v in variants2.values():
        counts2[v['classification']] += 1

    print(f"\n{'Classification Counts':^40}")
    print("-" * 50)
    print(f"{'Classification':<20} {label1:>12} {label2:>12}")
    print("-" * 50)
    for cls in CLASSIFICATIONS:
        print(f"{cls:<20} {counts1[cls]:>12} {counts2[cls]:>12}")
    print("-" * 50)
    print(f"{'TOTAL':<20} {sum(counts1.values()):>12} {sum(counts2.values()):>12}")

    # Variant matching summary
    print(f"\n{'Variant Matching':^40}")
    print("-" * 50)
    print(f"Matched variants:      {len(matched):>6}")
    print(f"Only in {label1}:        {len(only1):>6}")
    print(f"Only in {label2}:        {len(only2):>6}")

    # Classification concordance for matched variants
    print(f"\n{'Classification Concordance (Matched Variants)':^50}")
    print("-" * 60)

    # Build concordance matrix
    concordance = defaultdict(lambda: defaultdict(int))
    exact_match = 0

    for key1, key2 in matched:
        cls1 = variants1[key1]['classification']
        cls2 = variants2[key2]['classification']
        concordance[cls1][cls2] += 1
        if cls1 == cls2:
            exact_match += 1

    # Print matrix header
    header = f"{'':>15}"
    for cls in CLASSIFICATIONS:
        header += f" {CLASS_ABBREV[cls]:>6}"
    print(f"\n{label2} ->")
    print(header)
    print(f"{label1} v")

    # Print matrix rows
    for cls1 in CLASSIFICATIONS:
        row = f"{CLASS_ABBREV[cls1]:>15}"
        for cls2 in CLASSIFICATIONS:
            count = concordance[cls1][cls2]
            row += f" {count:>6}"
        print(row)

    print("-" * 60)
    concordance_rate = (exact_match / len(matched) * 100) if matched else 0
    print(f"Exact classification match: {exact_match}/{len(matched)} ({concordance_rate:.1f}%)")

    # Pathogenic agreement (P or LP)
    path_agree = sum(concordance[c1][c2] for c1 in ['pathogenic', 'likely pathogenic']
                     for c2 in ['pathogenic', 'likely pathogenic'])
    benign_agree = sum(concordance[c1][c2] for c1 in ['benign', 'likely benign']
                       for c2 in ['benign', 'likely benign'])

    print(f"Pathogenic direction agreement (P/LP in both): {path_agree}")
    print(f"Benign direction agreement (B/LB in both): {benign_agree}")


def print_code_summary(variants1, variants2, matched, label1, label2):
    """Print summary of ACMG code concordance."""
    print(f"\n{'ACMG Code Concordance (Matched Variants)':^50}")
    print("-" * 70)

    # Aggregate code statistics
    code_stats = defaultdict(lambda: {'both': 0, 'only1': 0, 'only2': 0})

    for key1, key2 in matched:
        v1 = variants1[key1]
        v2 = variants2[key2]

        codes1 = parse_rationale(v1['rationale'])
        codes2 = parse_rationale(v2['rationale'])
        comparison = compare_codes(codes1, codes2)

        for code in comparison['both']:
            code_stats[code]['both'] += 1
        for code in comparison['only1']:
            code_stats[code]['only1'] += 1
        for code in comparison['only2']:
            code_stats[code]['only2'] += 1

    # Print code concordance table
    print(f"{'Code':<8} {'Both':>8} {label1+' only':>12} {label2+' only':>12} {'Concordance':>12}")
    print("-" * 70)

    for code in sorted(code_stats.keys()):
        stats = code_stats[code]
        total = stats['both'] + stats['only1'] + stats['only2']
        concordance = (stats['both'] / total * 100) if total > 0 else 0
        print(f"{code:<8} {stats['both']:>8} {stats['only1']:>12} {stats['only2']:>12} {concordance:>11.1f}%")


def write_detailed_output(output_file, variants1, variants2, matched, only1, only2, label1, label2):
    """Write detailed per-variant comparison to output file."""

    with open(output_file, 'w') as f:
        # Header
        headers = [
            'chromosome', 'position', 'refAllele', 'altAllele',
            'gene', 'consequence', 'variantType',
            f'classification_{label1}', f'classification_{label2}', 'classification_match',
            f'transcript_{label1}', f'transcript_{label2}',
            'codes_both', f'codes_only_{label1}', f'codes_only_{label2}', 'code_score_diffs',
            f'rationale_{label1}', f'rationale_{label2}',
            'match_status'
        ]
        f.write('\t'.join(headers) + '\n')

        # Write matched variants
        for key1, key2 in matched:
            v1 = variants1[key1]
            v2 = variants2[key2]

            codes1 = parse_rationale(v1['rationale'])
            codes2 = parse_rationale(v2['rationale'])
            code_comparison = compare_codes(codes1, codes2)

            cls_match = "yes" if v1['classification'] == v2['classification'] else "no"

            row = [
                v1['chromosome'],
                v1['position'],
                v1['refAllele'],
                v1['altAllele'],
                v1['geneName'],
                v1['consequence'],
                v1['variantType'],
                v1['classification'],
                v2['classification'],
                cls_match,
                v1['transcript'],
                v2['transcript'],
                ','.join(sorted(code_comparison['both'])) or '-',
                ','.join(sorted(code_comparison['only1'])) or '-',
                ','.join(sorted(code_comparison['only2'])) or '-',
                ';'.join(f"{c}:{s1}->{s2}" for c, (s1, s2) in code_comparison['score_diffs'].items()) or '-',
                v1['rationale'],
                v2['rationale'],
                'matched'
            ]
            f.write('\t'.join(row) + '\n')

        # Write variants only in file1
        for key1 in sorted(only1):
            v1 = variants1[key1]
            codes1 = parse_rationale(v1['rationale'])
            active_codes = ','.join(sorted(get_active_codes(codes1))) or '-'

            row = [
                v1['chromosome'],
                v1['position'],
                v1['refAllele'],
                v1['altAllele'],
                v1['geneName'],
                v1['consequence'],
                v1['variantType'],
                v1['classification'],
                '-',
                'n/a',
                v1['transcript'],
                '-',
                '-',
                active_codes,
                '-',
                '-',
                v1['rationale'],
                '-',
                f'only_{label1}'
            ]
            f.write('\t'.join(row) + '\n')

        # Write variants only in file2
        for key2 in sorted(only2):
            v2 = variants2[key2]
            codes2 = parse_rationale(v2['rationale'])
            active_codes = ','.join(sorted(get_active_codes(codes2))) or '-'

            row = [
                v2['chromosome'],
                v2['position'],
                v2['refAllele'],
                v2['altAllele'],
                v2['geneName'],
                v2['consequence'],
                v2['variantType'],
                '-',
                v2['classification'],
                'n/a',
                '-',
                v2['transcript'],
                '-',
                '-',
                active_codes,
                '-',
                '-',
                v2['rationale'],
                f'only_{label2}'
            ]
            f.write('\t'.join(row) + '\n')

    print(f"\nDetailed comparison written to: {output_file}")


def main():
    options = parse_args()

    print(f"Loading {options.label1}: {options.file1}")
    variants1, keys_map1 = parse_bias_output(options.file1)
    print(f"  Loaded {len(variants1)} variants")

    print(f"Loading {options.label2}: {options.file2}")
    variants2, keys_map2 = parse_bias_output(options.file2)
    print(f"  Loaded {len(variants2)} variants")

    # Match variants
    matched, only1, only2 = match_variants(variants1, keys_map1, variants2, keys_map2)

    # Print summaries
    print_summary(variants1, variants2, matched, only1, only2, options.label1, options.label2)
    print_code_summary(variants1, variants2, matched, options.label1, options.label2)

    # Write detailed output
    write_detailed_output(options.output_file, variants1, variants2, matched, only1, only2,
                          options.label1, options.label2)


if __name__ == "__main__":
    main()
