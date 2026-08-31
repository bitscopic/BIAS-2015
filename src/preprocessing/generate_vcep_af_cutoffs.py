"""
Generate a gene-level TSV of VCEP (Variant Curation Expert Panel) allele-frequency
cutoffs for BA1, BS1, and PM2 from the ClinGen Criteria Specification Registry (CSpec).

The registry publishes ~70 released Sequence Variant Interpretation (SVI) specifications,
each attached to one or more genes. Each specification recommends numeric AF thresholds
that override the ACMG/AMP 2015 defaults. BIAS ingests the TSV produced here and uses
these gene-specific thresholds in place of its LOEUF-tiered defaults whenever the gene
appears in the table.

The CSpec API is JSON-LD. Structural fields are machine-readable, but the numeric
thresholds are buried in free-text `description` strings on each `evidenceStrengths[]`
entry and must be regex-extracted. Rows that can't be parsed are written to a sidecar
_review.tsv so a human can hand-patch them; the script exits non-zero if the sidecar has
any rows unless --allow-unparseable is passed.

TSV columns:
    gene_symbol          HGNC gene symbol
    hgnc_id              CSpec-provided gene identifier URL (not a clean HGNC:xxxxx —
                         the API returns a genenames.org search URL and we pass it
                         through verbatim)
    acmg_code            BA1 | BS1 | PM2
    strength             StandAlone | Strong | Moderate | Supporting | VeryStrong
                         | NotApplicable
    threshold            numeric MAF; 0 for PM2 "absent from controls" rules
                         (classifier fires iff variant AF <= 0, i.e. absent from
                         gnomAD); blank when strength=NotApplicable
    comparator           >= | <= | > | <
    metric               raw_af | faf95 | grpmax | popmax
    population           free-text population qualifier (e.g., any_continental, AFR)
    min_allele_count     minimum allele count qualifier, blank if unspecified
    mode_of_inheritance  "" (any) | AR | AD | XL | Mito. Populated when a single
                         VCEP description encodes per-MoI thresholds (e.g., hearing
                         loss VCEPs specify distinct AF cutoffs for autosomal
                         recessive vs. dominant); one row emitted per applicable MoI.
    gn_id                source SVI id (e.g., GN021)
    version              SVI version string
    source_url           full URL to the CSpec doc for audit
    raw_description      unmodified description string, retained for review
"""
import argparse
import json
import logging
import re
import sys
import urllib.parse
import urllib.request


CSPEC_INDEX_URL = "https://cspec.genome.network/cspec/api/svis"
CSPEC_DOC_URL_TMPL = "https://cspec.genome.network/cspec/api/SequenceVariantInterpretation/id/{gn_id}"
CSPEC_UI_URL_TMPL = "https://cspec.genome.network/cspec/ui/svi/doc/{gn_id}?version={version}"

TARGET_CODES = {"BA1", "BS1", "PM2"}
DEFAULT_COMPARATOR = {"BA1": ">=", "BS1": ">=", "PM2": "<"}
RELEASED_STATUSES = {"Released", "Approved", "Approved For Release", "Released - Under Revision"}

TSV_COLUMNS = [
    "gene_symbol",
    "hgnc_id",
    "acmg_code",
    "strength",
    "threshold",
    "comparator",
    "metric",
    "population",
    "min_allele_count",
    "mode_of_inheritance",
    "gn_id",
    "version",
    "source_url",
    "raw_description",
]

REVIEW_COLUMNS = [
    "gn_id",
    "version",
    "gene_symbol",
    "acmg_code",
    "raw_description",
    "parse_failure_reason",
]

STRENGTH_LABEL_NORMALIZATION = {
    "stand alone": "StandAlone",
    "standalone": "StandAlone",
    "very strong": "VeryStrong",
    "strong": "Strong",
    "moderate": "Moderate",
    "supporting": "Supporting",
    "not applicable": "NotApplicable",
    "n/a": "NotApplicable",
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True,
                        help="Path to the primary TSV output.")
    parser.add_argument("--review-output", default=None,
                        help="Path to the review sidecar TSV (default: <output>_review.tsv).")
    parser.add_argument("--allow-unparseable", action="store_true",
                        help="Exit 0 even if the review sidecar has rows. Default: exit non-zero.")
    parser.add_argument("--gn-ids", default=None,
                        help="Comma-separated list of specific GN IDs to fetch (default: all released).")
    parser.add_argument("--verbose", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        default="INFO")
    options = parser.parse_args()
    logging.basicConfig(level=getattr(logging, options.verbose), format="%(message)s")
    if options.review_output is None:
        if options.output.endswith(".tsv"):
            options.review_output = options.output[:-4] + "_review.tsv"
        else:
            options.review_output = options.output + "_review.tsv"
    return options


def fetch_json(url):
    """Fetch a URL and return decoded JSON."""
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.loads(resp.read().decode("utf-8"))


def gn_id_from_url(url):
    """Extract 'GN021' from a URL tail like '.../GN021' or '.../GN021?version=1.0.0'."""
    tail = url.rstrip("/").split("/")[-1]
    return tail.split("?")[0]


def enumerate_gn_ids(explicit_list):
    """Fetch the SVI index and return [(gn_id, status), ...] sorted."""
    if explicit_list:
        return [(gn.strip(), "explicit") for gn in explicit_list.split(",") if gn.strip()]
    index = fetch_json(CSPEC_INDEX_URL)
    if isinstance(index, list):
        entries = index
    elif isinstance(index, dict):
        entries = index.get("data") or index.get("@graph") or []
    else:
        entries = []
    result = []
    for entry in entries:
        url = entry.get("url") or entry.get("@id", "")
        gn_id = gn_id_from_url(url)
        if not gn_id.startswith("GN"):
            continue
        status = entry.get("status", "")
        result.append((gn_id, status))
    return sorted(set(result))


def normalize_strength(label):
    """Map a strength label from CSpec into our canonical vocabulary."""
    if not label:
        return None
    key = label.strip().lower()
    return STRENGTH_LABEL_NORMALIZATION.get(key, label.strip())


def normalize_applicability(entry):
    """Return True if this evidence strength entry is applicable to this VCEP."""
    applicability = (entry.get("applicability") or "").strip().lower()
    if not applicability:
        return False
    if applicability.startswith("not applicable") or applicability in {"n/a", "no", "false"}:
        return False
    return True


# ---------------- description parsing helpers ----------------

# CSpec descriptions occasionally use HTML "x 10<sup>-N</sup>" notation or spell out
# small integers ("One out of 100,000"). Normalize those to standard `1.43e-7` and
# `1/100000` forms before running the numeric extractors.
_SCI_NOTATION_RE = re.compile(
    r"(?P<mant>\d+(?:\.\d+)?)\s*[xX×]\s*10\s*(?:<sup>)?\s*(?P<exp>-?\d+)\s*(?:</sup>)?"
)
_WORD_NUMS = {"one": 1, "two": 2, "three": 3, "four": 4, "five": 5,
              "six": 6, "seven": 7, "eight": 8, "nine": 9, "ten": 10}
_WORD_RATIO_RE = re.compile(
    r"\b(" + "|".join(_WORD_NUMS) + r")\s+out\s+of\s+(?P<denom>[\d,]+)",
    re.IGNORECASE,
)


def _normalize_notation(text):
    """Rewrite non-standard number formats to standard so the main regex catches them."""
    text = _SCI_NOTATION_RE.sub(lambda m: f"{m.group('mant')}e{m.group('exp')}", text)
    def _repl_word(m):
        n = _WORD_NUMS[m.group(1).lower()]
        return f"{n}/{m.group('denom')}"
    text = _WORD_RATIO_RE.sub(_repl_word, text)
    return text

# Fractions like "0.007", "7e-3", "0.0035" — require ≥1 non-zero digit after the dot
# so we don't accidentally match "0.0" out of "0.01%".
_FRACTION_RE = re.compile(
    r"(?<![\w.])(?P<frac>\d\.\d*[1-9]\d*(?:[eE][-+]?\d+)?|\d[eE][-+]?\d+)(?!\s*%)"
)
# Percents like "5%", "0.35%", ".5%" — allow a bare-decimal form (no leading 0).
_PERCENT_RE = re.compile(r"(?<![\w.])(?P<pct>\d*\.?\d+)\s*%")
# 1:N ratios like "1:10,000" or "1/10000" meaning 1/N.
_RATIO_RE = re.compile(r"(?<![\d.])1\s*[:/]\s*(?P<denom>\d[\d,]{2,})")
_COMPARATOR_RE = re.compile(r"(?P<cmp>>=|<=|>|<|≥|≤)")
# Anchors that indicate the following number is a threshold value. Includes symbolic
# comparators, word-form comparators, and topical keywords ("MAF cutoff of", "AF of").
# Order matters: longer phrases must come first so the alternation matches greedily.
_ANCHOR_ALT = (
    r"greater\s+than\s+or\s+equal\s+to"
    r"|greater\s+(?:than\s+)?or\s+equal\s+to"
    r"|less\s+than\s+or\s+equal\s+to"
    r"|less\s+(?:than\s+)?or\s+equal\s+to"
    r"|greater\s+than"
    r"|less\s+than"
    r"|at\s+least"
    r"|at\s+most"
    r"|no\s+more\s+than"
    r"|no\s+less\s+than"
    r"|above"
    r"|below"
    r"|of\s+at\s+least"
    r"|cutoff\s+of"
    r"|threshold\s+of"
    r"|maf\s+of"
    r"|af\s+of"
    r"|frequency\s+of"
    r"|frequency\s+is"
    r"|(?:maf|af|faf)\s*cutoff"
    r"|>=|<=|>|<|≥|≤"
)
# Number alternatives — percent variant is listed FIRST so "0.00333%" matches as
# 0.00333% (percent) rather than being partially consumed as the bare fraction 0.00333.
_CMP_NUM_RE = re.compile(
    r"(?P<cmp>" + _ANCHOR_ALT + r")\s*\**\s*"
    r"(?P<num>"
    r"\d*\.?\d+\s*%"
    r"|1\s*[:/]\s*\d[\d,]{2,}"
    r"|\d\.\d*[1-9]\d*(?:[eE][-+]?\d+)?"
    r"|\d[eE][-+]?\d+"
    r")",
    re.IGNORECASE,
)
_MIN_AC_RE = re.compile(
    r"(?:>|≥|>=|at\s+least|minimum(?:\s+of)?)\s*(?P<n>\d[\d,]*)\s+"
    r"(?:alleles|allele\s+count|observed\s+alleles)",
    re.IGNORECASE,
)
# Keywords that, when they appear shortly before an anchor+number match, indicate
# the number is NOT an allele-frequency threshold (e.g., "allele balance below 0.35",
# "penetrance at 30%", "prevalence of 1 in 5,000").
_NON_AF_CONTEXT_RE = re.compile(
    r"(?:allele\s+balance|variant\s+allele\s+fraction|penetrance|prevalence"
    r"|heterogeneity|allelic\s+heterogeneity|genetic\s+heterogeneity"
    r"|coverage\s+at|observed\s+alleles)",
    re.IGNORECASE,
)


def _to_fraction(raw):
    """Convert '0.007', '0.7%', '.5%', '1:10,000' to a float fraction, or None."""
    raw = raw.strip()
    if raw.endswith("%"):
        try:
            return float(raw[:-1].strip()) / 100.0
        except ValueError:
            return None
    if ":" in raw or "/" in raw:
        sep = ":" if ":" in raw else "/"
        left, right = raw.split(sep, 1)
        try:
            num = float(left.strip().replace(",", ""))
            denom = float(right.strip().replace(",", ""))
        except ValueError:
            return None
        if denom > 0:
            return num / denom
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def _metric_from_text(text):
    lower = text.lower()
    if "faf95" in lower or "filtering allele frequency" in lower or "faf " in lower or lower.startswith("faf"):
        return "faf95"
    if "grpmax" in lower or "group max" in lower:
        return "grpmax"
    if "popmax" in lower or "population max" in lower:
        return "popmax"
    return "raw_af"


def _population_from_text(text):
    lower = text.lower()
    if "any continental population" in lower or "continental" in lower:
        return "any_continental"
    m = re.search(r"\b(afr|amr|asj|eas|fin|nfe|sas|oth)\b", lower)
    if m:
        return m.group(1).upper()
    if "overall" in lower or "any population" in lower or "total" in lower:
        return "all"
    return ""


# Matches "≥0.005 (0.5%) for autosomal recessive" or "≤0.00002 for autosomal dominant"
# — the phrase-form that some VCEPs use to encode per-MoI thresholds within a single
# evidenceStrengths description.
_MOI_CONDITIONAL_RE = re.compile(
    r"(?P<cmp>>=|<=|>|<|≥|≤)?\s*\**\s*"
    r"(?P<num>\d*\.?\d+\s*%|\d\.\d+(?:[eE][-+]?\d+)?|\d[eE][-+]?\d+)"
    r"\s*(?:\([^)]*\))?\s*"
    r"(?:in|for)\s*(?:the\s+context\s+of\s+)?autosomal\s+(?P<moi>recessive|dominant)",
    re.IGNORECASE,
)


def parse_moi_conditional(description, code):
    """
    Detect per-MoI thresholds encoded in a single description string.
    Returns a list of (moi, parsed_dict) tuples, or [] if not MoI-conditional.
    """
    description = _normalize_notation(description)
    matches = list(_MOI_CONDITIONAL_RE.finditer(description))
    if len(matches) < 2:
        return []
    by_moi = {}
    for m in matches:
        moi_word = m.group("moi").lower()
        moi = "AR" if moi_word == "recessive" else "AD"
        v = _to_fraction(m.group("num"))
        if v is None or not (0 < v < 1):
            continue
        cmp_raw = m.group("cmp")
        if cmp_raw:
            comparator = {"≥": ">=", "≤": "<="}.get(cmp_raw, cmp_raw)
        else:
            comparator = DEFAULT_COMPARATOR.get(code, ">=")
        # Only keep the first (typical) threshold if the same MoI appears twice.
        if moi not in by_moi:
            by_moi[moi] = {
                "threshold": v,
                "comparator": comparator,
                "metric": _metric_from_text(description),
                "population": _population_from_text(description),
                "min_allele_count": _min_allele_count_from_text(description),
            }
    if len(by_moi) < 2:
        return []
    return sorted(by_moi.items())


def _min_allele_count_from_text(text):
    ac_match = _MIN_AC_RE.search(text)
    if ac_match:
        return ac_match.group("n").replace(",", "")
    return ""


# Recognizes PM2 descriptions that express "variant must be essentially absent from
# population controls" without giving a numeric AF threshold. Handled by encoding
# threshold=0 with comparator "<=" — a variant absent from gnomAD (AF=0) satisfies
# this while any observed variant (AF>0) does not.
_ABSENT_RE = re.compile(
    r"\b(?:absent(?:/rare)?|must\s+be\s+absent|one\s+or\s+fewer\s+alleles"
    r"|absent\s+in\s+gnomad)\b",
    re.IGNORECASE,
)


def parse_absent_rule(description, code):
    """
    For PM2 rules whose description is a semantic "absent from controls" statement
    rather than a numeric threshold, return a parsed dict with threshold=0 so the
    classifier fires PM2 iff the variant is absent from gnomAD.
    """
    if code != "PM2" or not description:
        return None
    if not _ABSENT_RE.search(description):
        return None
    return {
        "threshold": 0.0,
        "comparator": "<=",
        "metric": _metric_from_text(description),
        "population": _population_from_text(description),
        "min_allele_count": "",  # semantic rule, don't gate on AN
    }


def parse_af_description(description, code):
    """
    Extract (threshold, comparator, metric, population, min_allele_count) from
    a free-text description. Returns None if no numeric threshold can be found.
    """
    if not description:
        return None
    text = _normalize_notation(description.strip())

    # Prefer numbers that appear right after a comparator symbol — those are the
    # explicit thresholds. Skip matches where the preceding context indicates a
    # non-AF quantity (allele balance, penetrance, prevalence, etc.).
    candidates = []
    for m in _CMP_NUM_RE.finditer(text):
        # Check the 80 chars before this match for non-AF context keywords.
        prefix = text[max(0, m.start() - 80):m.start()]
        if _NON_AF_CONTEXT_RE.search(prefix):
            continue
        raw = m.group("num").strip()
        v = _to_fraction(raw)
        if v is not None and 0 < v < 1:
            candidates.append(v)

    if not candidates:
        return None

    # For BA1, filter out the generic ACMG default (0.05 = 5%) when VCEP-specific
    # candidates exist — some descriptions quote the ACMG default as context before
    # giving the VCEP override (e.g., "above 5% ... gnomAD >= 0.004").
    if code == "BA1" and len(candidates) > 1:
        vcep_candidates = [c for c in candidates if abs(c - 0.05) > 1e-9]
        if vcep_candidates:
            candidates = vcep_candidates

    # BS1 is the lower threshold in a range (">=X but <Y"): take the smallest.
    # BA1 is a stand-alone upper threshold: take the largest.
    # PM2 is a rarity threshold: take the smallest.
    if code == "BS1" or code == "PM2":
        threshold = min(candidates)
    else:
        threshold = max(candidates)

    cmp_match = _COMPARATOR_RE.search(text)
    if cmp_match:
        raw = cmp_match.group("cmp")
        comparator = {"≥": ">=", "≤": "<="}.get(raw, raw)
    else:
        comparator = DEFAULT_COMPARATOR.get(code, ">=")

    metric = _metric_from_text(text)
    population = _population_from_text(text)

    min_ac = _min_allele_count_from_text(text)

    return {
        "threshold": threshold,
        "comparator": comparator,
        "metric": metric,
        "population": population,
        "min_allele_count": min_ac,
    }


# ---------------- SVI walking ----------------

def _iter_criteria_codes(rule_set):
    for cc in rule_set.get("criteriaCodes", []) or []:
        yield cc


def _iter_evidence_strengths(criteria_code):
    for ev in criteria_code.get("evidenceStrengths", []) or []:
        yield ev


def _extract_genes(rule_set):
    """Return list of (gene_symbol, hgnc_id) tuples for a ruleSet."""
    genes = []
    for gene in rule_set.get("genes", []) or []:
        symbol = gene.get("label") or gene.get("symbol") or ""
        hgnc_id = gene.get("@id") or gene.get("hgnc_id") or ""
        if symbol:
            genes.append((symbol.strip(), hgnc_id.strip()))
    return genes


def walk_svi(gn_id, doc, out_rows, review_rows):
    """
    Walk one SVI document, appending TSV rows for BA1/BS1/PM2 rules and review
    rows for unparseable descriptions.
    """
    version = doc.get("version") or ""
    source_url = CSPEC_UI_URL_TMPL.format(
        gn_id=gn_id, version=urllib.parse.quote(str(version))
    )

    rule_sets = doc.get("ruleSets") or []
    if not rule_sets:
        logging.debug("%s: no ruleSets", gn_id)
        return

    for rule_set in rule_sets:
        genes = _extract_genes(rule_set)
        if not genes:
            logging.debug("%s: ruleSet has no genes; skipping", gn_id)
            continue

        for criteria_code in _iter_criteria_codes(rule_set):
            label = (criteria_code.get("label") or "").strip().upper()
            if label not in TARGET_CODES:
                continue

            applicable_entries = []
            for ev in _iter_evidence_strengths(criteria_code):
                if normalize_applicability(ev):
                    applicable_entries.append(ev)

            if not applicable_entries:
                # The whole code is not applicable to any strength for this VCEP.
                # Emit one NotApplicable row per gene so BIAS suppresses the code.
                for symbol, hgnc_id in genes:
                    _emit_row(out_rows, symbol, hgnc_id, label, "NotApplicable",
                              None, "", gn_id, version, source_url, "")
                continue

            for ev in applicable_entries:
                strength = normalize_strength(ev.get("label")) or ""
                description = (ev.get("description") or "").strip()
                rows = parse_description(description, label)
                if not rows:
                    for symbol, _ in genes:
                        review_rows.append({
                            "gn_id": gn_id, "version": version,
                            "gene_symbol": symbol, "acmg_code": label,
                            "raw_description": description,
                            "parse_failure_reason": "no_numeric_threshold_found",
                        })
                    continue
                for moi, parsed in rows:
                    for symbol, hgnc_id in genes:
                        _emit_row(out_rows, symbol, hgnc_id, label, strength,
                                  parsed, moi, gn_id, version, source_url,
                                  description)


def _format_threshold(value):
    """Format a float threshold. 0 → "0"; small values use sci-notation."""
    if value == 0:
        return "0"
    if value >= 1e-4:
        return f"{value:.10g}"
    return f"{value:.3e}"


def _emit_row(out_rows, symbol, hgnc_id, label, strength, parsed, moi,
              gn_id, version, source_url, description):
    """Append one TSV row. `parsed` is a dict from the parse functions, or None
    for NotApplicable rows (in which case all threshold fields are blank)."""
    if parsed is None:
        out_rows.append({
            "gene_symbol": symbol, "hgnc_id": hgnc_id, "acmg_code": label,
            "strength": strength, "threshold": "", "comparator": "",
            "metric": "", "population": "", "min_allele_count": "",
            "mode_of_inheritance": moi,
            "gn_id": gn_id, "version": version, "source_url": source_url,
            "raw_description": description,
        })
    else:
        out_rows.append({
            "gene_symbol": symbol, "hgnc_id": hgnc_id, "acmg_code": label,
            "strength": strength,
            "threshold": _format_threshold(parsed["threshold"]),
            "comparator": parsed["comparator"], "metric": parsed["metric"],
            "population": parsed["population"],
            "min_allele_count": parsed["min_allele_count"],
            "mode_of_inheritance": moi,
            "gn_id": gn_id, "version": version, "source_url": source_url,
            "raw_description": description,
        })


def parse_description(description, code):
    """
    Try all three parsing strategies in order and return a list of
    (moi, parsed_dict) pairs — one per row to emit. Returns [] if none apply.

    Strategy order:
      1. MoI-conditional  — description encodes distinct per-MoI thresholds; emit
                            one row per MoI.
      2. Numeric          — standard parser (comparator-anchored number extraction).
      3. Absent-rule      — PM2-only semantic fallback for "absent from controls".
    """
    moi_rows = parse_moi_conditional(description, code)
    if moi_rows:
        return moi_rows
    parsed = parse_af_description(description, code) or parse_absent_rule(description, code)
    if parsed is not None:
        return [("", parsed)]
    return []


def _sort_rows(rows):
    return sorted(rows, key=lambda r: (
        r.get("gn_id", ""),
        r.get("gene_symbol", ""),
        r.get("acmg_code", ""),
        r.get("strength", ""),
    ))


def _write_tsv(path, columns, rows):
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for row in rows:
            f.write("\t".join(_clean_field(row.get(c, "")) for c in columns) + "\n")


def _clean_field(value):
    if value is None:
        return ""
    return str(value).replace("\t", " ").replace("\n", " ").replace("\r", " ")


def main():
    options = parse_args()

    logging.info("Enumerating SVIs from %s", CSPEC_INDEX_URL)
    gn_entries = enumerate_gn_ids(options.gn_ids)
    logging.info("Discovered %d SVI entries", len(gn_entries))

    out_rows = []
    review_rows = []

    fetched = 0
    skipped_by_status = 0
    for gn_id, status in gn_entries:
        if not options.gn_ids and status and status not in RELEASED_STATUSES:
            logging.debug("Skipping %s: status=%s", gn_id, status)
            skipped_by_status += 1
            continue
        url = CSPEC_DOC_URL_TMPL.format(gn_id=gn_id)
        logging.info("Fetching %s (status=%s)", gn_id, status)
        try:
            doc = fetch_json(url)
        except Exception as exc:
            logging.warning("Failed to fetch %s: %s", gn_id, exc)
            review_rows.append({
                "gn_id": gn_id,
                "version": "",
                "gene_symbol": "",
                "acmg_code": "",
                "raw_description": "",
                "parse_failure_reason": f"fetch_failed: {exc}",
            })
            continue
        fetched += 1
        walk_svi(gn_id, doc, out_rows, review_rows)

    out_rows = _sort_rows(out_rows)
    review_rows = sorted(review_rows, key=lambda r: (
        r.get("gn_id", ""), r.get("gene_symbol", ""), r.get("acmg_code", "")
    ))

    _write_tsv(options.output, TSV_COLUMNS, out_rows)
    _write_tsv(options.review_output, REVIEW_COLUMNS, review_rows)

    logging.info("Fetched %d SVIs (%d skipped by status)", fetched, skipped_by_status)
    logging.info("Wrote %d rules to %s", len(out_rows), options.output)
    logging.info("Wrote %d review rows to %s", len(review_rows), options.review_output)

    if review_rows and not options.allow_unparseable:
        logging.error("Sidecar has %d unparseable rows. Fix regex or pass --allow-unparseable.",
                      len(review_rows))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
