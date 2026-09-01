"""
Missense o/e upper (`oe_mis_upper`) × MOI stratified AF cutoff resolver.

The `mis_oe_upper_binned_cutoffs.tsv` table gives BA1/BS1/PM2_Supporting AF cutoffs
binned by `oe_mis_upper` (gnomAD missense observed/expected upper bound) and stratified
by mode of inheritance (AD vs AR). This module owns:

  1) Loading the cutoff table into a nested dict.
  2) Loading the precomputed gene→MOI lookup (HPO > GenCC > ClinGen priority).
  3) Resolving a (variant, rule_label) request into a single cutoff by:
       - looking up MOI from the gene→MOI table, falling back to ClinGen validity,
       - clamping oe_mis_upper into the covered bin range,
       - returning the matching cutoff (or None if MOI is not AD/AR or oe_mis_upper is missing).
"""


def load_gene_to_moi(fp):
    """
    Load precomputed gene→MOI lookup from TSV (gene_symbol, inheritance_mode).
    Built from HPO > GenCC > ClinGen at calibration time.
    """
    gene_to_moi = {}
    with open(fp, "r", encoding="utf-8") as f:
        header = f.readline()
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                gene_to_moi[parts[0].strip()] = parts[1].strip()
    return gene_to_moi


def load_mis_oe_moi_af_cutoffs(fp):
    """
    Parse the mis_oe_upper × MOI cutoffs TSV.

    Returns:
        {moi: {rule: [(bin_lo, bin_hi, cutoff), ...]}} sorted by bin_lo.
    """
    table = {}
    with open(fp, "r", encoding="utf-8") as in_file:
        header = in_file.readline().rstrip("\n").split("\t")
        col = {name: header.index(name) for name in ["moi", "rule", "bin_lo", "bin_hi", "cutoff_value"]}
        for line in in_file:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(header):
                continue
            moi = fields[col["moi"]].strip()
            rule = fields[col["rule"]].strip()
            bin_lo = float(fields[col["bin_lo"]])
            bin_hi = float(fields[col["bin_hi"]])
            cutoff = float(fields[col["cutoff_value"]])
            table.setdefault(moi, {}).setdefault(rule, []).append((bin_lo, bin_hi, cutoff))

    for moi_rules in table.values():
        for rows in moi_rules.values():
            rows.sort(key=lambda row: row[0])
    return table


def resolve_moi(gene_name, gene_to_moi=None, clingen_gene_validity=None):
    """
    Pick a single MOI ("AD" or "AR") for a gene.

    Resolution order:
      1. Precomputed gene→MOI table (HPO > GenCC > ClinGen, built at calibration time).
         "mixed" → "AD" (conservative: AD cutoffs are tighter than AR).
         "XL" → "AD" (X-linked genes use AD cutoffs).
      2. ClinGen validity entries on the variant (runtime fallback).

    Returns "AD", "AR", or None.
    """
    # Level 1: precomputed lookup
    if gene_to_moi and gene_name:
        moi = gene_to_moi.get(gene_name)
        if moi == "AD" or moi == "mixed" or moi == "XL":
            return "AD"
        if moi == "AR":
            return "AR"

    # Level 2: ClinGen validity entries on the variant
    if not clingen_gene_validity:
        return None
    has_ad = False
    has_ar = False
    for entry in clingen_gene_validity:
        inheritance = (entry.get("inheritance") or "").lower()
        if "semi" in inheritance:
            continue
        if "autosomal dominant" in inheritance:
            has_ad = True
        elif "autosomal recessive" in inheritance:
            has_ar = True
    if has_ad:
        return "AD"
    if has_ar:
        return "AR"
    return None


def get_mis_oe_moi_cutoff(cutoff_table, variant, rule_label, gene_to_moi=None):
    """
    Look up the AF cutoff for a variant under the mis_oe × MOI table.

    Returns:
        (cutoff, moi, bin_lo, oe_mis_upper) on success,
        None if the variant's gene has no AD/AR MOI or no oe_mis_upper.

    Out-of-range oe_mis_upper is clamped to the lowest/highest available bin.
    """
    if not cutoff_table:
        return None

    gene_name = getattr(variant, "geneName", None)
    moi = resolve_moi(gene_name, gene_to_moi, getattr(variant, "clingen_gene_validity", None))
    if moi not in ("AD", "AR"):
        return None

    gene_gnomad = getattr(variant, "gene_gnomad", None) or {}
    oe_mis = gene_gnomad.get("oe_mis_upper")
    if oe_mis is None:
        return None

    rows = cutoff_table.get(moi, {}).get(rule_label)
    if not rows:
        return None

    if oe_mis < rows[0][0]:
        return rows[0][2], moi, rows[0][0], oe_mis
    if oe_mis >= rows[-1][1]:
        return rows[-1][2], moi, rows[-1][0], oe_mis
    for bin_lo, bin_hi, cutoff in rows:
        if bin_lo <= oe_mis < bin_hi:
            return cutoff, moi, bin_lo, oe_mis
    return rows[-1][2], moi, rows[-1][0], oe_mis
