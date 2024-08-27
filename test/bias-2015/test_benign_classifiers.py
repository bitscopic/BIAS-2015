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

def test_get_bs1():
    """
    Test the get_bs1 function
    """
    # Variant with One Thousand Genome AF between 1% and 5%
    variant = type('Variant', (object,), {})()
    variant.oneKg = {'allAf': 0.03}
    variant.gnomad = None
    score, bs1 = benign_classifiers.get_bs1(variant)
    expected_score = 1
    expected_bs1 = "BS1 (1): One Thousand Genome general population AF=3.00% which is between 1% and 5%."
    assert (score, bs1) == (expected_score, expected_bs1)

    # Variant with Gnomad AF between 1% and 5%
    variant = type('Variant', (object,), {})()
    variant.oneKg = None
    variant.gnomad = {'allAf': 0.02}
    score, bs1 = benign_classifiers.get_bs1(variant)
    expected_score = 1
    expected_bs1 = "BS1 (1): Gnomad general population AF=2.00% which is between 1% and 5%."
    assert (score, bs1) == (expected_score, expected_bs1)

    # Variant with One Thousand Genome AF greater than 5%
    variant = type('Variant', (object,), {})()
    variant.oneKg = {'allAf': 0.06}
    variant.gnomad = None
    score, bs1 = benign_classifiers.get_bs1(variant)
    expected_score = 0
    expected_bs1 = ""
    assert (score, bs1) == (expected_score, expected_bs1)

    # Variant with Gnomad AF greater than 5%
    variant = type('Variant', (object,), {})()
    variant.oneKg = None
    variant.gnomad = {'allAf': 0.06}
    score, bs1 = benign_classifiers.get_bs1(variant)
    expected_score = 0
    expected_bs1 = ""
    assert (score, bs1) == (expected_score, expected_bs1)

    # Variant with One Thousand Genome AF less than or equal to 1%
    variant = type('Variant', (object,), {})()
    variant.oneKg = {'allAf': 0.005}
    variant.gnomad = None
    score, bs1 = benign_classifiers.get_bs1(variant)
    expected_score = 0
    expected_bs1 = ""
    assert (score, bs1) == (expected_score, expected_bs1)

    # Variant with Gnomad AF less than or equal to 1%
    variant = type('Variant', (object,), {})()
    variant.oneKg = None
    variant.gnomad = {'allAf': 0.008}
    score, bs1 = benign_classifiers.get_bs1(variant)
    expected_score = 0
    expected_bs1 = ""
    assert (score, bs1) == (expected_score, expected_bs1)

    # Variant with neither One Thousand Genome nor Gnomad AF data
    variant = type('Variant', (object,), {})()
    variant.oneKg = None
    variant.gnomad = None
    score, bs1 = benign_classifiers.get_bs1(variant)
    expected_score = 0
    expected_bs1 = ""
    assert (score, bs1) == (expected_score, expected_bs1)

def test_get_bp1():
    """
    Test the get_bp1 function
    """
    # Variant is a missense variant in a gene primarily associated with truncating variants
    truncating_genes = {"GENE1"}
    variant = type('Variant', (object,), {})()
    variant.geneName = "GENE1"
    variant.consequence = "missense_variant"
    score, bp1 = benign_classifiers.get_bp1(variant, truncating_genes)
    expected_score = 1
    expected_bp1 = "BP1: Missense variant type missense_variant in gene GENE1 which has over 80% truncating pathogenic variants"
    assert (score, bp1) == (expected_score, expected_bp1)

    # Variant is a truncating variant in a gene primarily associated with truncating variants
    truncating_genes = {"GENE1"}
    variant = type('Variant', (object,), {})()
    variant.geneName = "GENE1"
    variant.consequence = "stop_gained"
    score, bp1 = benign_classifiers.get_bp1(variant, truncating_genes)
    expected_score = 0
    expected_bp1 = ""
    assert (score, bp1) == (expected_score, expected_bp1)

    # Variant is a missense variant in a gene not associated with truncating variants
    truncating_genes = {"GENE2"}
    variant = type('Variant', (object,), {})()
    variant.geneName = "GENE1"
    variant.consequence = "missense_variant"
    score, bp1 = benign_classifiers.get_bp1(variant, truncating_genes)
    expected_score = 0
    expected_bp1 = ""
    assert (score, bp1) == (expected_score, expected_bp1)

    # Variant is a missense variant, but geneName is None
    truncating_genes = {"GENE1"}
    variant = type('Variant', (object,), {})()
    variant.geneName = None
    variant.consequence = "missense_variant"
    score, bp1 = benign_classifiers.get_bp1(variant, truncating_genes)
    expected_score = 0
    expected_bp1 = ""
    assert (score, bp1) == (expected_score, expected_bp1)

    # Variant is a missense variant, but consequence is None
    truncating_genes = {"GENE1"}
    variant = type('Variant', (object,), {})()
    variant.geneName = "GENE1"
    variant.consequence = None
    score, bp1 = benign_classifiers.get_bp1(variant, truncating_genes)
    expected_score = 0
    expected_bp1 = ""
    assert (score, bp1) == (expected_score, expected_bp1)

def test_get_bp3():
    """
    Test the get_bp3 function.
    """

    # Define the repeat regions
    chrom_to_repeat_regions = {
        "chr1": [(1000, 2000, "repeat_region")]
    }

    # In-frame insertion in a repeat region
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

    # In-frame deletion in a repeat region with a length triggering additional score
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

    # Out-of-frame insertion (length not divisible by 3)
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

    # In-frame insertion, but not in a repeat region
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

    # In-frame deletion of a very small length (e.g., 1 bp) in a repeat region
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

    # Variant not in any repeat region
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

def test_get_bp4():
    """
    Test the get_bp4 function.
    """

    # All computational evidence suggests benign impact
    variant = type('Variant', (object,), {})()
    variant.revel = 0.20
    variant.dann = 0.20
    variant.gerp = 1.50
    score, bp4 = benign_classifiers.get_bp4(variant)
    expected_score = 3
    expected_bp4 = "BP4 (3): Computational evidence support a benign effect; gerp 1.5 | dann 0.2 | revel 0.2"
    assert (score, bp4) == (expected_score, expected_bp4)

    # Only revel suggests benign impact
    variant = type('Variant', (object,), {})()
    variant.revel = 0.20
    variant.dann = 0.30
    variant.gerp = 2.50
    score, bp4 = benign_classifiers.get_bp4(variant)
    expected_score = 1
    expected_bp4 = "BP4 (1): Computational evidence support a benign effect; gerp 2.5 | dann 0.3 | revel 0.2"
    assert (score, bp4) == (expected_score, expected_bp4)

    # Only dann suggests benign impact
    variant = type('Variant', (object,), {})()
    variant.revel = 0.30
    variant.dann = 0.20
    variant.gerp = 2.50
    score, bp4 = benign_classifiers.get_bp4(variant)
    expected_score = 1
    expected_bp4 = "BP4 (1): Computational evidence support a benign effect; gerp 2.5 | dann 0.2 | revel 0.3"
    assert (score, bp4) == (expected_score, expected_bp4)

    # Only gerp suggests benign impact
    variant = type('Variant', (object,), {})()
    variant.revel = 0.30
    variant.dann = 0.30
    variant.gerp = 1.50
    score, bp4 = benign_classifiers.get_bp4(variant)
    expected_score = 1
    expected_bp4 = "BP4 (1): Computational evidence support a benign effect; gerp 1.5 | dann 0.3 | revel 0.3"
    assert (score, bp4) == (expected_score, expected_bp4)

    # No evidence suggests benign impact
    variant = type('Variant', (object,), {})()
    variant.revel = 0.30
    variant.dann = 0.30
    variant.gerp = 2.50
    score, bp4 = benign_classifiers.get_bp4(variant)
    expected_score = 0
    expected_bp4 = ""
    assert (score, bp4) == (expected_score, expected_bp4)

def test_get_bp6():
    """
    Test the get_bp6 function.
    """
    # Variant is benign with high review status
    variant = type('Variant', (object,), {})()
    variant.clinvar_review_status = "criteria provided, multiple submitters, no conflicts"
    variant.clinvar_significance = "benign"
    score, bp6 = benign_classifiers.get_bp6(variant)
    expected_score = 3
    expected_bp6 = ("BP6 (3): Variant was found in ClinVar as benign with review status "
                    "of criteria provided, multiple submitters, no conflicts and given a "
                    "weighted PP5 value of 3")
    assert (score, bp6) == (expected_score, expected_bp6)

    # Variant is likely benign with lower review status
    variant = type('Variant', (object,), {})()
    variant.clinvar_review_status = "criteria provided, single submitter"
    variant.clinvar_significance = "likely benign"
    score, bp6 = benign_classifiers.get_bp6(variant)
    expected_score = 2
    expected_bp6 = ("BP6 (2): Variant was found in ClinVar as likely benign with review status of criteria provided, single submitter and given a weighted PP5 value of 2")
    assert (score, bp6) == (expected_score, expected_bp6)

    # Variant has uncertain significance with high review status
    variant = type('Variant', (object,), {})()
    variant.clinvar_review_status = "criteria provided, multiple submitters, no conflicts"
    variant.clinvar_significance = "uncertain significance"
    score, bp6 = benign_classifiers.get_bp6(variant)
    expected_score = 0
    expected_bp6 = ""
    assert (score, bp6) == (expected_score, expected_bp6)

    # Variant is benign but with a low review status
    variant = type('Variant', (object,), {})()
    variant.clinvar_review_status = "no assertion provided"
    variant.clinvar_significance = "benign"
    score, bp6 = benign_classifiers.get_bp6(variant)
    expected_score = 0
    expected_bp6 = ""
    assert (score, bp6) == (expected_score, expected_bp6)

    # Variant is benign but with an unknown review status
    variant = type('Variant', (object,), {})()
    variant.clinvar_review_status = "unknown status"
    variant.clinvar_significance = "benign"
    score, bp6 = benign_classifiers.get_bp6(variant)
    expected_score = 0
    expected_bp6 = ""
    assert (score, bp6) == (expected_score, expected_bp6)

def test_get_bp7():
    """
    Test the get_bp7 function.
    """
    # Test case 1: Variant has a synonymous consequence without splice involvement
    variant = type('Variant', (object,), {})()
    variant.consequence = "synonymous_variant"
    variant.geneName = "gene1"
    score, bp7 = benign_classifiers.get_bp7(variant)
    expected_score = 1
    expected_bp7 = "BP7 (1): Variant has synonymous associated consequence synonymous_variant"
    assert (score, bp7) == (expected_score, expected_bp7)

    # Test case 2: Variant has an intronic consequence without splice involvement
    variant = type('Variant', (object,), {})()
    variant.consequence = "intron_variant"
    variant.geneName = "gene2"
    score, bp7 = benign_classifiers.get_bp7(variant)
    expected_score = 1
    expected_bp7 = "BP7 (1): Variant has intronic associated consequence intron_variant"
    assert (score, bp7) == (expected_score, expected_bp7)

    # Variant is intergenic 
    variant = type('Variant', (object,), {})()
    variant.consequence = "intergenic_variant"
    variant.geneName = "intergenic"
    score, bp7 = benign_classifiers.get_bp7(variant)
    expected_score = 1  
    expected_bp7 = "BP7 (1): Intergenic variant" 
    assert (score, bp7) == (expected_score, expected_bp7)

    # Test case 4: Variant has a synonymous consequence but involves splice
    variant = type('Variant', (object,), {})()
    variant.consequence = "synonymous_variant&splice_region_variant"
    variant.geneName = "gene1"
    score, bp7 = benign_classifiers.get_bp7(variant)
    expected_score = 0
    expected_bp7 = ""
    assert (score, bp7) == (expected_score, expected_bp7)

    # Test case 5: Variant has an intronic consequence but involves splice
    variant = type('Variant', (object,), {})()
    variant.consequence = "intron_variant&splice_region_variant"
    variant.geneName = "gene2"
    score, bp7 = benign_classifiers.get_bp7(variant)
    expected_score = 0
    expected_bp7 = ""
    assert (score, bp7) == (expected_score, expected_bp7)

    # Test case 6: Variant has a non-synonymous, non-intronic, non-intergenic consequence
    variant = type('Variant', (object,), {})()
    variant.consequence = "missense_variant"
    variant.geneName = "gene3"
    score, bp7 = benign_classifiers.get_bp7(variant)
    expected_score = 0
    expected_bp7 = ""
    assert (score, bp7) == (expected_score, expected_bp7)

    # Test case 7: Variant is intergenic but with splice involvement
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