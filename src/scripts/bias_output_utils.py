"""
Shared utilities for parsing and comparing BIAS output files.

This module provides common functions used by:
  - compare_bias_outputs.py (comparing two BIAS outputs)
  - evaluate_bias_erepo_validation.py (comparing BIAS to eRepo ground truth)

Functions include:
  - BIAS output file parsing
  - Variant key normalization for indel matching
  - Rationale/code parsing
  - Classification utilities
  - Concordance and metrics calculations
"""
import json
from collections import defaultdict


# Classification categories
CLASSIFICATIONS = ['pathogenic', 'likely pathogenic', 'uncertain', 'likely benign', 'benign']

CLASS_ABBREV = {
    'pathogenic': 'P',
    'likely pathogenic': 'LP',
    'uncertain': 'VUS',
    'likely benign': 'LB',
    'benign': 'B'
}

# Numeric scores for classifications (negative = pathogenic direction)
CLASS_TO_SCORE = {
    'benign': 2,
    'likely benign': 1,
    'uncertain': 0,
    'uncertain significance': 0,
    'likely pathogenic': -1,
    'pathogenic': -2
}


def normalize_variant_key(chrom, pos, ref, alt):
    """
    Normalize variant representation to handle different indel formats.

    Different annotators (VEP, Nirvana, etc.) represent indels differently:
    - Nirvana/VCF-style: chr1, 171605496, CAG, C (with padding base)
    - VEP simplified: chr1, 171605497, AG, - (no padding base)

    Args:
        chrom: Chromosome
        pos: Position (string)
        ref: Reference allele
        alt: Alternate allele

    Returns:
        List of possible variant key tuples to enable matching
    """
    keys = [(chrom, pos, ref, alt)]

    # Generate alternative representations for indels
    if len(ref) > 1 and len(alt) == 1 and ref[0] == alt[0]:
        # Deletion: CAG>C -> AG>- at pos+1
        new_pos = str(int(pos) + 1)
        new_ref = ref[1:]
        new_alt = "-"
        keys.append((chrom, new_pos, new_ref, new_alt))

    if len(alt) > 1 and len(ref) == 1 and alt[0] == ref[0]:
        # Insertion: C>CAG -> ->AG at pos+1
        new_pos = str(int(pos) + 1)
        new_ref = "-"
        new_alt = alt[1:]
        keys.append((chrom, new_pos, new_ref, new_alt))

    # Handle reverse direction (VEP -> VCF style)
    if alt == "-" and len(ref) > 0:
        new_pos = str(int(pos) - 1)
        keys.append((chrom, new_pos, "X" + ref, "X"))

    if ref == "-" and len(alt) > 0:
        new_pos = str(int(pos) - 1)
        keys.append((chrom, new_pos, "X", "X" + alt))

    return keys


def get_variant_type(ref, alt):
    """
    Classify a variant based on its REF and ALT alleles.

    Args:
        ref: Reference allele
        alt: Alternate allele

    Returns:
        One of: 'snv', 'mnv', 'insertion', 'deletion', 'indel'
    """
    if ref == "-":
        return "insertion"
    if alt == "-":
        return "deletion"

    ref_len = len(ref)
    alt_len = len(alt)

    if ref_len == 1 and alt_len == 1:
        return "snv"
    elif ref_len == alt_len:
        return "mnv"
    elif ref_len < alt_len:
        return "insertion"
    elif ref_len > alt_len:
        return "deletion"
    else:
        return "indel"


def parse_bias_output(filepath):
    """
    Parse a BIAS output file and return variant data.

    Args:
        filepath: Path to BIAS output TSV file

    Returns:
        Tuple of (variants_dict, keys_map)
        - variants_dict: variant_key -> {classification, rationale, gene, ...}
        - keys_map: normalized_key -> original_key (for indel matching)
    """
    variants = {}
    variant_keys_map = {}

    with open(filepath, 'r') as f:
        header = f.readline().strip().split('\t')
        col_idx = {name: i for i, name in enumerate(header)}

        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < len(header):
                continue

            chrom = fields[col_idx['chromosome']]
            pos = fields[col_idx['position']]
            ref = fields[col_idx['refAllele']]
            alt = fields[col_idx['altAllele']]

            primary_key = (chrom, pos, ref, alt)

            data = {
                'chromosome': chrom,
                'position': pos,
                'refAllele': ref,
                'altAllele': alt,
                'variantType': fields[col_idx.get('variantType', 4)] if 'variantType' in col_idx else '',
                'consequence': fields[col_idx.get('consequence', 5)] if 'consequence' in col_idx else '',
                'classification': fields[col_idx.get('acmgClassification', 6)] if 'acmgClassification' in col_idx else '',
                'geneName': fields[col_idx.get('geneName', 12)] if 'geneName' in col_idx else '',
                'transcript': fields[col_idx.get('transcript', 16)] if 'transcript' in col_idx else '',
                'rationale': fields[col_idx.get('rationale', -1)] if 'rationale' in col_idx else '{}'
            }

            variants[primary_key] = data

            # Store normalized keys for matching
            for norm_key in normalize_variant_key(chrom, pos, ref, alt):
                variant_keys_map[norm_key] = primary_key

    return variants, variant_keys_map


def parse_rationale(rationale_str):
    """
    Parse the rationale JSON and extract code scores.

    Args:
        rationale_str: JSON string containing ACMG code evidence

    Returns:
        dict: code -> (score, explanation)
    """
    codes = {}
    try:
        data = json.loads(rationale_str)
        for category, cat_codes in data.items():
            for code, value in cat_codes.items():
                if isinstance(value, list) and len(value) >= 2:
                    codes[code] = (value[0], value[1])
                elif isinstance(value, list) and len(value) == 1:
                    codes[code] = (value[0], "")
    except (json.JSONDecodeError, TypeError):
        pass
    return codes


def get_active_codes(codes_dict):
    """
    Return set of codes with score > 0.

    Args:
        codes_dict: dict from parse_rationale()

    Returns:
        Set of active code names
    """
    return {code for code, (score, _) in codes_dict.items() if score > 0}


def compare_codes(codes1, codes2):
    """
    Compare two sets of ACMG codes.

    Args:
        codes1: dict from parse_rationale() for first variant
        codes2: dict from parse_rationale() for second variant

    Returns:
        dict with keys: 'both', 'only1', 'only2', 'score_diffs'
    """
    active1 = get_active_codes(codes1)
    active2 = get_active_codes(codes2)

    both = active1 & active2
    only1 = active1 - active2
    only2 = active2 - active1

    # Check for score differences in codes present in both
    score_diffs = {}
    for code in both:
        score1 = codes1.get(code, (0, ""))[0]
        score2 = codes2.get(code, (0, ""))[0]
        if score1 != score2:
            score_diffs[code] = (score1, score2)

    return {
        'both': both,
        'only1': only1,
        'only2': only2,
        'score_diffs': score_diffs
    }


def match_variants(variants1, keys_map1, variants2, keys_map2):
    """
    Match variants between two datasets, handling different indel representations.

    Args:
        variants1: dict from parse_bias_output() for first file
        keys_map1: keys map from parse_bias_output() for first file
        variants2: dict from parse_bias_output() for second file
        keys_map2: keys map from parse_bias_output() for second file

    Returns:
        Tuple of (matched, only1, only2)
        - matched: list of (key1, key2) tuples for matched variants
        - only1: set of keys only in file1
        - only2: set of keys only in file2
    """
    matched = []
    matched_keys1 = set()
    matched_keys2 = set()

    for key1 in variants1.keys():
        for norm_key in normalize_variant_key(*key1):
            if norm_key in keys_map2:
                key2 = keys_map2[norm_key]
                if key2 not in matched_keys2:
                    matched.append((key1, key2))
                    matched_keys1.add(key1)
                    matched_keys2.add(key2)
                    break
            elif norm_key in variants2:
                if norm_key not in matched_keys2:
                    matched.append((key1, norm_key))
                    matched_keys1.add(key1)
                    matched_keys2.add(norm_key)
                    break

    only1 = set(variants1.keys()) - matched_keys1
    only2 = set(variants2.keys()) - matched_keys2

    return matched, only1, only2


def calculate_metrics(confusion_matrix):
    """
    Calculate sensitivity, specificity, precision, and F1 score.

    Args:
        confusion_matrix: dict with keys 'TP', 'FP', 'TN', 'FN'

    Returns:
        Tuple of (sensitivity, specificity, precision, f1_score)
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


def calculate_concordance_metrics(agree, fp, fn):
    """
    Calculate concordance metrics for code-level comparison.

    Args:
        agree: Number of agreements (true positives)
        fp: Number of false positives (called by algorithm but not in reference)
        fn: Number of false negatives (in reference but not called)

    Returns:
        Tuple of (concordance, precision, recall, f1_score)
    """
    total = agree + fp + fn
    concordance = agree / total if total > 0 else 0
    precision = agree / (agree + fp) if (agree + fp) > 0 else 0
    recall = agree / (agree + fn) if (agree + fn) > 0 else 0
    f1_score = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    return concordance, precision, recall, f1_score
