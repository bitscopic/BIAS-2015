"""
Unit testing for the Bitscopic Interpreting ACMG Standards (BIAS) benign classifiers
"""
import os
import sys
# NOTE: This should be handled more elegantly for public facing code
# Some hacky path adjustments so users dont have to fiddle with env variables
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from src.bias_2015 import benign_classifiers, extract_from_nirvana_json

clinvar_review_status_to_level = {
    "practice guideline": 5,
    "reviewed by expert panel": 5,
    "criteria provided, multiple submitters, no conflicts": 3,
    "criteria provided, single submitter": 2,
    "criteria provided, conflicting interpretations": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0
}

def test_get_ba1():
    """
    Test the get_ba1 function
    """
    # Variant with One Thousand Genome AF greater than 5%
    variant = type('Variant', (object,), {})()
    variant.oneKg = {'allAf': 0.06}
    variant.gnomad = None
    score, ba1 = benign_classifiers.get_ba1(variant)
    expected_score = 1
    expected_ba1 = "BA1 (1): One Thousand Genome general population AF=6.00% which is greater than 5%."
    assert (score, ba1) == (expected_score, expected_ba1)

    # Variant with Gnomad AF greater than 5%
    variant = type('Variant', (object,), {})()
    variant.oneKg = None
    variant.gnomad = {'allAf': 0.07}
    score, ba1 = benign_classifiers.get_ba1(variant)
    expected_score = 1
    expected_ba1 = "BA1 (1): Gnomad general population AF=7.00% which is greater than 5%."
    assert (score, ba1) == (expected_score, expected_ba1)

    # Variant with One Thousand Genome AF less than or equal to 5%
    variant = type('Variant', (object,), {})()
    variant.oneKg = {'allAf': 0.05}
    variant.gnomad = None
    score, ba1 = benign_classifiers.get_ba1(variant)
    expected_score = 0
    expected_ba1 = ""
    assert (score, ba1) == (expected_score, expected_ba1)

    # Variant with Gnomad AF less than or equal to 5%
    variant = type('Variant', (object,), {})()
    variant.oneKg = None
    variant.gnomad = {'allAf': 0.04}
    score, ba1 = benign_classifiers.get_ba1(variant)
    expected_score = 0
    expected_ba1 = ""
    assert (score, ba1) == (expected_score, expected_ba1)

    # Variant with neither One Thousand Genome nor Gnomad AF data
    variant = type('Variant', (object,), {})()
    variant.oneKg = None
    variant.gnomad = None
    score, ba1 = benign_classifiers.get_ba1(variant)
    expected_score = 0
    expected_ba1 = ""
    assert (score, ba1) == (expected_score, expected_ba1)

def test_get_ba1_with_high_one_thousand_genome_af():
    """
    Test get_ba1 function with a variant having One Thousand Genome AF greater than 5%.
    """
    variant = type('Variant', (object,), {})()
    variant.oneKg = {'allAf': 0.06}
    variant.gnomad = None
    score, ba1 = benign_classifiers.get_ba1(variant)
    expected_score = 1
    expected_ba1 = "BA1 (1): One Thousand Genome general population AF=6.00% which is greater than 5%."
    assert (score, ba1) == (expected_score, expected_ba1)

def test_get_ba1_with_high_gnomad_af():
    """
    Test get_ba1 function with a variant having Gnomad AF greater than 5%.
    """
    variant = type('Variant', (object,), {})()
    variant.oneKg = None
    variant.gnomad = {'allAf': 0.07}
    score, ba1 = benign_classifiers.get_ba1(variant)
    expected_score = 1
    expected_ba1 = "BA1 (1): Gnomad general population AF=7.00% which is greater than 5%."
    assert (score, ba1) == (expected_score, expected_ba1)

def test_get_ba1_with_low_one_thousand_genome_af():
    """
    Test get_ba1 function with a variant having One Thousand Genome AF less than or equal to 5%.
    """
    variant = type('Variant', (object,), {})()
    variant.oneKg = {'allAf': 0.05}
    variant.gnomad = None
    score, ba1 = benign_classifiers.get_ba1(variant)
    expected_score = 0
    expected_ba1 = ""
    assert (score, ba1) == (expected_score, expected_ba1)

def test_get_ba1_with_low_gnomad_af():
    """
    Test get_ba1 function with a variant having Gnomad AF less than or equal to 5%.
    """
    variant = type('Variant', (object,), {})()
    variant.oneKg = None
    variant.gnomad = {'allAf': 0.04}
    score, ba1 = benign_classifiers.get_ba1(variant)
    expected_score = 0
    expected_ba1 = ""
    assert (score, ba1) == (expected_score, expected_ba1)

def test_get_ba1_with_no_af_data():
    """
    Test get_ba1 function with a variant having neither One Thousand Genome nor Gnomad AF data.
    """
    variant = type('Variant', (object,), {})()
    variant.oneKg = None
    variant.gnomad = None
    score, ba1 = benign_classifiers.get_ba1(variant)
    expected_score = 0
    expected_ba1 = ""
    assert (score, ba1) == (expected_score, expected_ba1)

def test_get_bp3_inframe_insertion_in_repeat_region():
    """
    Test get_bp3 function with an in-frame insertion in a repeat region.
    """
    chrom_to_repeat_regions = {
        "chr1": [(1000, 2000, "repeat_region")]
    }
    variant = type('Variant', (object,), {})()
    variant.chromosome = "chr1"
    variant.position = 1500
    variant.refAllele = "A"
    variant.altAllele = "AAT"
    variant.consequence = "inframe_insertion"
    score, bp3 = benign_classifiers.get_bp3(variant, chrom_to_repeat_regions)
    expected_score = 1
    expected_bp3 = "BP3 (1): In-frame INDEL of length 2 in repeat region chr1 1000-2000"
    assert (score, bp3) == (expected_score, expected_bp3)

def test_get_bp3_inframe_deletion_in_repeat_region_with_additional_score(): #fix
    """
    Test get_bp3 function with an in-frame deletion in a repeat region triggering additional score.
    """
    chrom_to_repeat_regions = {
        "chr1": [(1000, 2000, "repeat_region")]
    }
    variant = type('Variant', (object,), {})()
    variant.chromosome = "chr1"
    variant.position = 1500
    variant.refAllele = "AAT"
    variant.altAllele = "A"
    variant.consequence = "inframe_deletion"
    score, bp3 = benign_classifiers.get_bp3(variant, chrom_to_repeat_regions)
    expected_score = 3
    expected_bp3 = "BP3 (3): In-frame INDEL of length 2 in repeat region chr1 1000-2000"
    assert (score, bp3) == (expected_score, expected_bp3)

def test_get_bp3_out_of_frame_insertion():
    """
    Test get_bp3 function with an out-of-frame insertion (length not divisible by 3).
    """
    chrom_to_repeat_regions = {
        "chr1": [(1000, 2000, "repeat_region")]
    }
    variant = type('Variant', (object,), {})()
    variant.chromosome = "chr1"
    variant.position = 1500
    variant.refAllele = "A"
    variant.altAllele = "AA"
    variant.consequence = "inframe_insertion"
    score, bp3 = benign_classifiers.get_bp3(variant, chrom_to_repeat_regions)
    expected_score = 0
    expected_bp3 = ""
    assert (score, bp3) == (expected_score, expected_bp3)

def test_get_bp3_inframe_insertion_not_in_repeat_region():
    """
    Test get_bp3 function with an in-frame insertion that is not in a repeat region.
    """
    chrom_to_repeat_regions = {
        "chr1": [(1000, 2000, "repeat_region")]
    }
    variant = type('Variant', (object,), {})()
    variant.chromosome = "chr1"
    variant.position = 2500
    variant.refAllele = "A"
    variant.altAllele = "AAT"
    variant.consequence = "inframe_insertion"
    score, bp3 = benign_classifiers.get_bp3(variant, chrom_to_repeat_regions)
    expected_score = 0
    expected_bp3 = ""
    assert (score, bp3) == (expected_score, expected_bp3)

def test_get_bp3_inframe_deletion_small_length_in_repeat_region(): #fix
    """
    Test get_bp3 function with an in-frame deletion of a very small length in a repeat region.
    """
    chrom_to_repeat_regions = {
        "chr1": [(1000, 2000, "repeat_region")]
    }
    variant = type('Variant', (object,), {})()
    variant.chromosome = "chr1"
    variant.position = 1500
    variant.refAllele = "ATG"
    variant.altAllele = "A"
    variant.consequence = "inframe_deletion"
    score, bp3 = benign_classifiers.get_bp3(variant, chrom_to_repeat_regions)
    expected_score = 4
    expected_bp3 = "BP3 (4): In-frame INDEL of length 2 in repeat region chr1 1000-2000"
    assert (score, bp3) == (expected_score, expected_bp3)

def test_get_bp3_variant_not_in_repeat_region():
    """
    Test get_bp3 function with a variant that is not in any repeat region.
    """
    chrom_to_repeat_regions = {
        "chr1": [(1000, 2000, "repeat_region")]
    }
    variant = type('Variant', (object,), {})()
    variant.chromosome = "chr1"
    variant.position = 500
    variant.refAllele = "A"
    variant.altAllele = "AA"
    variant.consequence = "inframe_insertion"
    score, bp3 = benign_classifiers.get_bp3(variant, chrom_to_repeat_regions)
    expected_score = 0
    expected_bp3 = ""
    assert (score, bp3) == (expected_score, expected_bp3)

def test_get_bp4_all_computational_evidence_benign():
    """
    Test the get_bp4 function when all computational evidence suggests a benign impact.
    """
    variant = type('Variant', (object,), {})()
    variant.revel = 0.20
    variant.dann = 0.20
    variant.gerp = 1.50
    score, bp4 = benign_classifiers.get_bp4(variant)
    expected_score = 3
    expected_bp4 = "BP4 (3): Computational evidence support a benign effect; gerp 1.5 | dann 0.2 | revel 0.2"

    assert (score, bp4) == (expected_score, expected_bp4)

def test_get_bp4_only_revel_suggests(): 
    """
    Test the get_bp4 function when only revel suggests a benign impact.
    """
    variant = type('Variant', (object,), {})()
    variant.revel = 0.20
    variant.dann = 0.30
    variant.gerp = 2.50
    score, bp4 = benign_classifiers.get_bp4(variant)
    expected_score = 1
    expected_bp4 = "BP4 (1): Computational evidence support a benign effect; gerp 2.5 | dann 0.3 | revel 0.2"
    assert (score, bp4) == (expected_score, expected_bp4)

def test_get_bp4_only_dann_suggests(): 
    """
    Test the get_bp4 function when only dann suggests a benign impact.
    """
    variant = type('Variant', (object,), {})()
    variant.revel = 0.30
    variant.dann = 0.20
    variant.gerp = 2.50
    score, bp4 = benign_classifiers.get_bp4(variant)
    expected_score = 1
    expected_bp4 = "BP4 (1): Computational evidence support a benign effect; gerp 2.5 | dann 0.2 | revel 0.3"
    assert (score, bp4) == (expected_score, expected_bp4)

def test_get_bp4_only_gerp_suggests(): 
    """
    Test the get_bp4 function when only gerp suggests a benign impact.
    """
    variant = type('Variant', (object,), {})()
    variant.revel = 0.30
    variant.dann = 0.30
    variant.gerp = 1.50
    score, bp4 = benign_classifiers.get_bp4(variant)
    expected_score = 1
    expected_bp4 = "BP4 (1): Computational evidence support a benign effect; gerp 1.5 | dann 0.3 | revel 0.3"
    assert (score, bp4) == (expected_score, expected_bp4)

def test_get_bp4_no_evidence_suggests():
    """
    Test the get_bp4 function when no evidence suggests a benign impact.
    """
    variant = type('Variant', (object,), {})()
    variant.revel = 0.30
    variant.dann = 0.30
    variant.gerp = 2.50
    score, bp4 = benign_classifiers.get_bp4(variant)
    expected_score = 0
    expected_bp4 = ""
    assert (score, bp4) == (expected_score, expected_bp4)

def test_get_bp6_benign_high_review_status():
    """
    Test get_bp6 with a variant that is benign with high review status.
    """
    variant = type('Variant', (object,), {})()
    variant.clinvar_review_status = "criteria provided, multiple submitters, no conflicts"
    variant.clinvar_significance = "benign"
    score, bp6 = benign_classifiers.get_bp6(variant)
    expected_score = 3
    expected_bp6 = ("BP6 (3): Variant was found in ClinVar as benign with review status "
                    "of criteria provided, multiple submitters, no conflicts and given a "
                    "weighted PP5 value of 3")
    assert (score, bp6) == (expected_score, expected_bp6)

def test_get_bp6_likely_benign_lower_review_status():
    """
    Test get_bp6 with a variant that is likely benign with lower review status.
    """
    variant = type('Variant', (object,), {})()
    variant.clinvar_review_status = "criteria provided, single submitter"
    variant.clinvar_significance = "likely benign"
    score, bp6 = benign_classifiers.get_bp6(variant)
    expected_score = 2
    expected_bp6 = ("BP6 (2): Variant was found in ClinVar as likely benign with review status "
                    "of criteria provided, single submitter and given a weighted PP5 value of 2")
    assert (score, bp6) == (expected_score, expected_bp6)

def test_get_bp6_uncertain_significance_high_review_status():
    """
    Test get_bp6 with a variant that has uncertain significance with high review status.
    """
    variant = type('Variant', (object,), {})()
    variant.clinvar_review_status = "criteria provided, multiple submitters, no conflicts"
    variant.clinvar_significance = "uncertain significance"
    score, bp6 = benign_classifiers.get_bp6(variant)
    expected_score = 0
    expected_bp6 = ""
    assert (score, bp6) == (expected_score, expected_bp6)

def test_get_bp6_benign_low_review_status():
    """
    Test get_bp6 with a variant that is benign but with a low review status.
    """
    variant = type('Variant', (object,), {})()
    variant.clinvar_review_status = "no assertion provided"
    variant.clinvar_significance = "benign"
    score, bp6 = benign_classifiers.get_bp6(variant)
    expected_score = 0
    expected_bp6 = ""
    assert (score, bp6) == (expected_score, expected_bp6)

def test_get_bp6_benign_unknown_review_status():
    """
    Test get_bp6 with a variant that is benign but with an unknown review status.
    """
    variant = type('Variant', (object,), {})()
    variant.clinvar_review_status = "unknown status"
    variant.clinvar_significance = "benign"
    score, bp6 = benign_classifiers.get_bp6(variant)
    expected_score = 0
    expected_bp6 = ""
    assert (score, bp6) == (expected_score, expected_bp6)

import pytest

def test_get_bp7_synonymous_consequence():
    """
    Test get_bp7 with a synonymous consequence without splice involvement.
    """
    variant = type('Variant', (object,), {})()
    variant.consequence = "synonymous_variant"
    variant.geneName = "gene1"
    score, bp7 = benign_classifiers.get_bp7(variant)
    expected_score = 1
    expected_bp7 = "BP7 (1): Variant has synonymous associated consequence synonymous_variant"
    assert (score, bp7) == (expected_score, expected_bp7)

def test_get_bp7_intronic_consequence():
    """
    Test get_bp7 with an intronic consequence without splice involvement.
    """
    variant = type('Variant', (object,), {})()
    variant.consequence = "intron_variant"
    variant.geneName = "gene2"
    score, bp7 = benign_classifiers.get_bp7(variant)
    expected_score = 1
    expected_bp7 = "BP7 (1): Variant has intronic associated consequence intron_variant"
    assert (score, bp7) == (expected_score, expected_bp7)

def test_get_bp7_intergenic_consequence():
    """
    Test get_bp7 with an intergenic consequence.
    """
    variant = type('Variant', (object,), {})()
    variant.consequence = "intergenic_variant"
    variant.geneName = "intergenic"
    score, bp7 = benign_classifiers.get_bp7(variant)
    expected_score = 1  
    expected_bp7 = "BP7 (1): Intergenic variant" 
    assert (score, bp7) == (expected_score, expected_bp7)

def test_get_bp7_synonymous_with_splice():
    """
    Test get_bp7 with a synonymous consequence that involves splice.
    """
    variant = type('Variant', (object,), {})()
    variant.consequence = "synonymous_variant&splice_region_variant"
    variant.geneName = "gene1"
    score, bp7 = benign_classifiers.get_bp7(variant)
    expected_score = 0
    expected_bp7 = ""
    assert (score, bp7) == (expected_score, expected_bp7)

def test_get_bp7_intronic_with_splice():
    """
    Test get_bp7 with an intronic consequence that involves splice.
    """
    variant = type('Variant', (object,), {})()
    variant.consequence = "intron_variant&splice_region_variant"
    variant.geneName = "gene2"
    score, bp7 = benign_classifiers.get_bp7(variant)
    expected_score = 0
    expected_bp7 = ""
    assert (score, bp7) == (expected_score, expected_bp7)

def test_get_bp7_non_synonymous_non_intronic_non_intergenic():
    """
    Test get_bp7 with a non-synonymous, non-intronic, non-intergenic consequence.
    """
    variant = type('Variant', (object,), {})()
    variant.consequence = "missense_variant"
    variant.geneName = "gene3"
    score, bp7 = benign_classifiers.get_bp7(variant)
    expected_score = 0
    expected_bp7 = ""
    assert (score, bp7) == (expected_score, expected_bp7)

def test_get_bp7_intergenic_with_splice():
    """
    Test get_bp7 with an intergenic consequence but with splice involvement.
    """
    variant = type('Variant', (object,), {})()
    variant.consequence = "splice_region_variant"
    variant.geneName = "intergenic"
    score, bp7 = benign_classifiers.get_bp7(variant)
    expected_score = 1
    expected_bp7 = "BP7 (1): Intergenic variant"
    assert (score, bp7) == (expected_score, expected_bp7)

#test_get_ba1 (done)
#test_get_bs1 (done)
#test_get_bp1 (done)
#test_get_bp3 (failing, repeat region) 
#test_get_bp4 (done)
#test_get_bp6 (done)
#test_get_bp7 (done)