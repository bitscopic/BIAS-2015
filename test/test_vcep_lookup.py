"""
Smoke tests for src.bias_2015.vcep_lookup and the CSpec preprocessing parser.

Classifiers call `get_vcep_rule` at the top of each AF-based code; when a rule
applies its VCEP threshold/strength replaces the tiered evaluation, otherwise
the classifier falls through to the mis_oe × MOI cascade and the flat ACMG
default. These tests exercise both the resolver and end-to-end classifier
behavior.

Run: python3 test/test_vcep_lookup.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))

from src.bias_2015.vcep_lookup import VcepAfRule, get_vcep_rule
from src.bias_2015.benign_classifiers import get_ba1, get_bs1
from src.bias_2015.pathogenic_classifiers import get_pm2
from src.preprocessing.generate_vcep_af_cutoffs import parse_af_description


class FakeVariant:
    def __init__(self, gene_name, gnomad_af=None, all_an=None, popmax=None,
                 inheritance=None, loeuf=1.0, onekg_af=None):
        self.geneName = gene_name
        self.gnomad = {}
        if gnomad_af is not None:
            self.gnomad["allAf"] = gnomad_af
        if all_an is not None:
            self.gnomad["allAn"] = all_an
        if popmax is not None:
            self.gnomad["popmax"] = popmax
        self.oneKg = {}
        if onekg_af is not None:
            self.oneKg["allAf"] = onekg_af
        self.gene_gnomad = {"loeuf": loeuf}
        self.clingen_gene_validity = []
        if inheritance:
            for inh in ([inheritance] if isinstance(inheritance, str) else inheritance):
                self.clingen_gene_validity.append({
                    "inheritance": inh, "disease": "test",
                    "diseaseId": "MONDO:0000000",
                })


def acadvl_rules():
    return {
        "ACADVL": {
            "BA1": [VcepAfRule("StandAlone", 0.007, ">=", "raw_af", "any_continental",
                               2000, "", "GN021", "2.2")],
            "BS1": [VcepAfRule("Strong", 0.0035, ">=", "raw_af", "any_continental",
                               2000, "", "GN021", "2.2")],
            "PM2": [VcepAfRule("Supporting", 0.001, "<", "raw_af", "any_continental",
                               2000, "", "GN021", "2.2")],
        }
    }


# ---------- get_vcep_rule tests ----------

def test_get_vcep_rule_returns_rule_when_gene_covered():
    v = FakeVariant("ACADVL", gnomad_af=0.01, all_an=3000)
    rule = get_vcep_rule(acadvl_rules(), "ACADVL", "BA1", v)
    assert rule is not None
    assert rule.gn_id == "GN021"
    assert rule.strength == "StandAlone"
    assert rule.threshold == 0.007
    print("PASS: get_vcep_rule returns the VCEP rule when gene is covered")


def test_get_vcep_rule_returns_none_when_gene_not_covered():
    v = FakeVariant("UNKNOWN", gnomad_af=0.01)
    assert get_vcep_rule(acadvl_rules(), "UNKNOWN", "BA1", v) is None
    print("PASS: get_vcep_rule returns None when gene has no VCEP entry")


def test_get_vcep_rule_returns_not_applicable():
    rules = {"SOMEGENE": {
        "BA1": [VcepAfRule("NotApplicable", None, "", "", "", 0, "", "GN999", "1.0")]
    }}
    v = FakeVariant("SOMEGENE", gnomad_af=0.9)
    rule = get_vcep_rule(rules, "SOMEGENE", "BA1", v)
    assert rule is not None
    assert rule.strength == "NotApplicable"
    assert rule.gn_id == "GN999"
    print("PASS: get_vcep_rule returns the NotApplicable rule for suppressed codes")


def test_get_vcep_rule_filters_by_moi_for_recessive_gene():
    rules = {"USH2A": {"PM2": [
        VcepAfRule("Supporting", 0.00007, "<=", "raw_af", "", 0, "AR", "GN005", "2.0"),
        VcepAfRule("Supporting", 0.00002, "<=", "raw_af", "", 0, "AD", "GN005", "2.0"),
    ]}}
    ar_variant = FakeVariant("USH2A", gnomad_af=0.00005, inheritance="autosomal recessive")
    rule = get_vcep_rule(rules, "USH2A", "PM2", ar_variant)
    assert rule.mode_of_inheritance == "AR"
    assert rule.threshold == 0.00007
    print("PASS: get_vcep_rule keeps only the MoI-matching rule for an AR gene")


def test_get_vcep_rule_returns_strongest():
    """Hearing loss GN005 has both Strong and Supporting BS1 rules; returns Strong."""
    rules = {"USH2A": {"BS1": [
        VcepAfRule("Strong", 0.003, ">=", "raw_af", "", 0, "AR", "GN005", "2.0"),
        VcepAfRule("Supporting", 0.0007, ">=", "raw_af", "", 0, "", "GN005", "2.0"),
    ]}}
    v = FakeVariant("USH2A", gnomad_af=0.001, inheritance="autosomal recessive")
    rule = get_vcep_rule(rules, "USH2A", "BS1", v)
    assert rule.strength == "Strong"
    assert rule.threshold == 0.003
    print("PASS: get_vcep_rule returns the strongest matching rule")


# ---------- classifier end-to-end tests ----------

def test_ba1_fires_at_vcep_threshold():
    v = FakeVariant("ACADVL", gnomad_af=0.01, all_an=3000)
    score, rationale = get_ba1(v, acadvl_rules())
    assert score == 5, f"expected BA1 fire (score 5), got {score}"
    assert "GN021" in rationale
    print("PASS: get_ba1 fires when variant exceeds VCEP threshold")


def test_ba1_does_not_fire_below_vcep_threshold():
    v = FakeVariant("ACADVL", gnomad_af=0.001, all_an=3000)
    score, _ = get_ba1(v, acadvl_rules())
    assert score == 0
    print("PASS: get_ba1 stays at 0 when AF below VCEP threshold")


def test_ba1_falls_back_to_acmg_default_for_uncovered_gene():
    # AF above ACMG BA1 default (5%), no VCEP rule, no mis_oe/MOI coverage → ACMG default fires
    v = FakeVariant("UNKNOWN_GENE", gnomad_af=0.10, loeuf=1.5)
    score, rationale = get_ba1(v, acadvl_rules())
    assert score == 5
    assert "ACMG default" in rationale
    print("PASS: get_ba1 uses the ACMG default when the gene has no VCEP or mis_oe rule")


def test_pm2_fires_on_absent_variant_via_vcep():
    """VCEP absent rule (threshold=0, <=) fires on variant with no gnomAD entry."""
    rules = {"KITGENE": {"PM2": [
        VcepAfRule("Supporting", 0.0, "<=", "raw_af", "", 0, "", "GN999", "1.0")
    ]}}
    v = FakeVariant("KITGENE")  # variant.gnomad = {}
    score, rationale = get_pm2(v, rules)
    assert score == 1, f"expected PM2_Supporting, got {score}"
    assert "GN999" in rationale
    print("PASS: get_pm2 fires on absent variant via VCEP absent-rule")


def test_pm2_not_applicable_suppresses_code():
    rules = {"GENEX": {"PM2": [
        VcepAfRule("NotApplicable", None, "", "", "", 0, "", "GN999", "1.0")
    ]}}
    v = FakeVariant("GENEX", gnomad_af=1e-9)
    score, rationale = get_pm2(v, rules)
    assert score == 0
    assert "Not Applicable" in rationale
    print("PASS: NotApplicable suppresses PM2 even for effectively-absent variants")


def test_bs1_supporting_downgrade_via_vcep():
    """VCEP-labeled Supporting strength returns score=1 even when AF is well above threshold."""
    rules = {"GJB2": {"BS1": [
        VcepAfRule("Supporting", 0.0007, ">=", "raw_af", "", 0, "", "GN005", "2.0")
    ]}}
    v = FakeVariant("GJB2", gnomad_af=0.01)  # well above 0.0007 threshold
    score, rationale = get_bs1(v, rules)
    assert score == 1, f"expected BS1_Supporting (score 1), got {score}"
    assert "supporting" in rationale
    print("PASS: VCEP Supporting-labeled BS1 returns score 1 (strength override wins)")


def test_bs1_fires_strong_at_data_derived_fallback_for_uncovered_gene():
    """No VCEP rule and no mis_oe/MOI coverage → data-derived fallback: BS1 Strong at AF > 0.00929."""
    v = FakeVariant("UNKNOWN", gnomad_af=0.01, loeuf=1.5)
    score, rationale = get_bs1(v, None)
    assert score == 3, f"expected BS1 Strong (0.01 > 0.00929), got {score}"
    assert "fallback Strong threshold" in rationale
    print("PASS: BS1 fires Strong at the data-derived fallback when no VCEP/mis_oe rule applies")


# ---------- parser tests ----------

def test_parse_af_description_ba1():
    desc = ("Variants with a highest population minor allele frequency (MAF) "
            "≥0.007 (0.7%) in any continental population with >2000 alleles "
            "in gnomAD will meet BA1. Prevalence of 1:30,000, Penetrance 0.75.")
    parsed = parse_af_description(desc, "BA1")
    assert parsed is not None
    assert abs(parsed["threshold"] - 0.007) < 1e-9
    print("PASS: parses BA1 threshold 0.007 despite penetrance/prevalence distractors")


def test_parse_af_description_pm2_prefers_smallest():
    desc = ("Case AF cutoff ≤0.02 in cases; a threshold of ≤0.00004 in the "
            "subpopulation with the highest frequency activates this rule.")
    parsed = parse_af_description(desc, "PM2")
    assert parsed is not None
    assert abs(parsed["threshold"] - 0.00004) < 1e-9
    print("PASS: PM2 prefers smallest candidate")


def main():
    tests = [v for k, v in globals().items() if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
    print(f"\n{len(tests)} tests passed.")


if __name__ == "__main__":
    main()
