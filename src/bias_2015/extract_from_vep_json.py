#!/usr/bin/env python3
# Chris Eisenhart - chris.eisenhart@bitscopic.com
"""
Extract information from a VEP (Variant Effect Predictor) JSON output file.

Each variant is stored in a VariantData class (shared with Nirvana ingestion).

The VEP JSON file is one variant per line, streamed through for memory efficiency.

VEP JSON structure (per line):
    - allele_string: "G/T" (ref/alt)
    - seq_region_name: chromosome
    - start/end: position
    - variant_class: "SNV", "insertion", "deletion", etc.
    - most_severe_consequence: top consequence term
    - transcript_consequences: array of transcript annotations
    - colocated_variants: array of overlapping known variants (dbSNP, gnomAD freqs, etc.)
"""
import json
import logging
import sys
import traceback
from src.bias_2015.variant_data import VariantData, open_file
from src.bias_2015 import bias_variant_classification


# Mapping table for three-letter to single-letter amino acids
AMINO_ACID_MAPPING = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q',
    'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F',
    'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Ter': '*'
}

# VEP variant_class to Bitscopic variantType mapping
VEP_VARIANT_CLASS_MAP = {
    "SNV": "SNV",
    "insertion": "insertion",
    "deletion": "deletion",
    "indel": "indel",
    "substitution": "MNV",
    "sequence_alteration": "indel",
}



def rank_transcript_entries(transcript_list, hgnc_to_gene_data):
    """
    Rank VEP transcript entries to select the best transcript.

    Criteria (highest to lowest priority):
    - Has AlphaMissense annotation: computational evidence is valuable
    - canonical: VEP marks the canonical transcript with canonical=1
    - protein_coding: biotype == "protein_coding"
    - Has coding info: cds_start/cds_end present
    - HGNC data available in gene dataset

    Args:
        transcript_list (list): VEP transcript_consequences array
        hgnc_to_gene_data (dict): Gene name to gene data mapping

    Returns:
        list: Transcript entries sorted by rank (best first)
    """
    rank_scores = {}
    for i, entry in enumerate(transcript_list):
        transcript_id = entry.get("transcript_id", f"unknown_{i}")
        rank_scores[transcript_id] = 0

        # Prioritize transcripts with AlphaMissense annotation
        if entry.get("alphamissense"):
            rank_scores[transcript_id] += 4

        if entry.get("canonical") == 1:
            rank_scores[transcript_id] += 3

        if entry.get("biotype") == "protein_coding":
            rank_scores[transcript_id] += 2

        if "cds_start" in entry or "protein_start" in entry:
            rank_scores[transcript_id] += 2

        if entry.get("gene_symbol") in hgnc_to_gene_data:
            rank_scores[transcript_id] += 1

    sorted_entries = sorted(
        transcript_list,
        key=lambda e: rank_scores.get(e.get("transcript_id", ""), 0),
        reverse=True
    )
    return sorted_entries


def convert_vep_mutation_format(hgvsp_str):
    """
    Convert VEP HGVS protein strings to shorthand format.

    VEP format: ENSP00000355651.4:p.Ala89Pro -> A89P
    VEP format: ENSP00000405344.2:p.Glu11Ter -> E11*

    Args:
        hgvsp_str (str): VEP hgvsp string

    Returns:
        str: Shorthand protein change (e.g., "A89P") or empty string if parsing fails
    """
    if not hgvsp_str or "?" in hgvsp_str:
        return ""

    try:
        # Extract the protein change portion after "p."
        if "p." not in hgvsp_str:
            return ""
        protein_change = hgvsp_str.split("p.")[1]

        # Handle VEP parenthesized format if present
        if protein_change.startswith("(") and protein_change.endswith(")"):
            protein_change = protein_change[1:-1]

        first_aa = AMINO_ACID_MAPPING.get(protein_change[:3])
        if not first_aa:
            return ""

        converted_str = first_aa
        if AMINO_ACID_MAPPING.get(protein_change[-3:]):
            base_number = protein_change[3:-3]
            second_aa = AMINO_ACID_MAPPING[protein_change[-3:]]
            converted_str = f"{first_aa}{base_number}{second_aa}"
        elif protein_change.endswith("="):
            base_number = protein_change[3:-1]
            converted_str = f"{first_aa}{base_number}*"
        elif "fs" in protein_change:
            if "Ter" in protein_change:
                term_num = protein_change.split("Ter")[1]
                pre_term_protein = protein_change.split("Ter")[0]
                base_number = pre_term_protein[3:-5]
                second_aa = AMINO_ACID_MAPPING.get(pre_term_protein[-5:-2], "")
                converted_str = f"{first_aa}{base_number}fs{second_aa}{term_num}*"
            else:
                # Frameshift without Ter
                base_number = protein_change[3:protein_change.index("fs")]
                converted_str = f"{first_aa}{base_number}fs"

        return converted_str
    except (IndexError, ValueError) as e:
        logging.debug("Could not parse VEP hgvsp: %s. Error: %s", hgvsp_str, e)
        return ""


def extract_gnomad_from_colocated(colocated_variants, alt_allele):
    """
    Extract gnomAD frequency data from VEP colocated_variants into a dict
    compatible with the Nirvana gnomad format used by BIAS classifiers.

    VEP stores gnomAD data in colocated_variants[].frequencies[alt_allele],
    with keys like gnomadg, gnomadg_afr, gnomadg_nfe, etc.

    Args:
        colocated_variants (list): VEP colocated_variants array
        alt_allele (str): The alternate allele to look up frequencies for

    Returns:
        dict: gnomAD data in Nirvana-compatible format, empty dict if not found
    """
    if not colocated_variants:
        return {}

    gnomad = {}
    for cv in colocated_variants:
        freqs = cv.get("frequencies", {})
        alt_freqs = freqs.get(alt_allele, {})
        if not alt_freqs:
            continue

        # Map VEP gnomAD fields to Nirvana gnomad format
        # VEP uses gnomadg (genomes) and gnomade (exomes)
        # Prefer genome data (gnomadg), fall back to exome (gnomade)
        all_af = alt_freqs.get("gnomadg", alt_freqs.get("gnomade", ""))
        if all_af != "":
            gnomad["allAf"] = all_af

        # Population-specific frequencies
        pop_mapping = {
            "afrAf": ["gnomadg_afr", "gnomade_afr"],
            "amrAf": ["gnomadg_amr", "gnomade_amr"],
            "easAf": ["gnomadg_eas", "gnomade_eas"],
            "nfeAf": ["gnomadg_nfe", "gnomade_nfe"],
            "sasAf": ["gnomadg_sas", "gnomade_sas"],
            "asjAf": ["gnomadg_asj", "gnomade_asj"],
            "finAf": ["gnomadg_fin", "gnomade_fin"],
            "midAf": ["gnomadg_mid", "gnomade_mid"],
            "amiAf": ["gnomadg_ami", "gnomade_ami"],
            "remainingAf": ["gnomadg_remaining", "gnomade_remaining"],
        }
        for nirvana_key, vep_keys in pop_mapping.items():
            for vep_key in vep_keys:
                val = alt_freqs.get(vep_key)
                if val is not None:
                    gnomad[nirvana_key] = val
                    break

        # If we found gnomAD data, no need to check more colocated variants
        if gnomad:
            break

    return gnomad


def extract_gnomad_from_custom_annotations(custom_annotations, alt_allele):
    """
    Extract gnomAD allele count and number data from VEP custom_annotations.

    VEP custom_annotations structure for gnomAD exomes:
        "custom_annotations": {
            "gnomAD_exomes": [{
                "allele": "T",
                "name": "rs765729544",
                "fields": {
                    "AN": 251494,
                    "AC": 4,
                    "AF": 1.5905e-05,
                    "AF_popmax": 0.000123047,
                    "FILTER": ["PASS"]
                }
            }]
        }

    These count fields (AC, AN) are not available from colocated_variants
    which only provide allele frequencies. The counts are needed for BS2
    classification (observed in healthy individuals).

    Args:
        custom_annotations (dict): VEP custom_annotations dictionary
        alt_allele (str): The alternate allele to match

    Returns:
        dict: gnomAD count data in Nirvana-compatible format, empty dict if not found
    """
    if not custom_annotations:
        return {}

    gnomad_entries = custom_annotations.get("gnomAD_exomes", [])
    if not gnomad_entries:
        return {}

    # Find the entry matching the alt allele
    best_entry = None
    for entry in gnomad_entries:
        if entry.get("allele", "") == alt_allele:
            best_entry = entry
            break

    # Fall back to first entry if no allele match
    if best_entry is None and gnomad_entries:
        best_entry = gnomad_entries[0]

    if best_entry is None:
        return {}

    fields = best_entry.get("fields", {})
    gnomad = {}

    ac = fields.get("AC")
    an = fields.get("AN")
    af = fields.get("AF")

    if ac is not None:
        gnomad["allAc"] = ac
        # Use total AC as a conservative proxy for controls AC since VEP
        # doesn't provide controls-specific counts. This is a slight
        # overestimate (includes non-control samples) but is directionally
        # correct for BS2 dominant evaluation.
        gnomad["controlsAllAc"] = ac
    if an is not None:
        gnomad["allAn"] = an
    if af is not None and "allAf" not in gnomad:
        gnomad["allAf"] = af
    if af is not None:
        # Estimate controlsAllAf from overall AF as an approximation
        gnomad["controlsAllAf"] = af

    return gnomad


def extract_onekg_from_colocated(colocated_variants, alt_allele):
    """
    Extract 1000 Genomes frequency data from VEP colocated_variants.

    VEP stores 1000 Genomes data in colocated_variants[].frequencies[alt_allele]
    with keys like af (global), eur, afr, amr, eas, sas.

    Args:
        colocated_variants (list): VEP colocated_variants array
        alt_allele (str): The alternate allele to look up frequencies for

    Returns:
        dict: 1000 Genomes data in Nirvana-compatible format, empty dict if not found
    """
    if not colocated_variants:
        return {}

    onekg = {}
    for cv in colocated_variants:
        freqs = cv.get("frequencies", {})
        alt_freqs = freqs.get(alt_allele, {})
        if not alt_freqs:
            continue

        # 1000 Genomes global frequency
        if "af" in alt_freqs:
            onekg["allAf"] = alt_freqs["af"]

        # Population-specific
        pop_mapping = {
            "afrAf": "afr",
            "amrAf": "amr",
            "eurAf": "eur",
            "easAf": "eas",
            "sasAf": "sas",
        }
        for nirvana_key, vep_key in pop_mapping.items():
            val = alt_freqs.get(vep_key)
            if val is not None:
                onekg[nirvana_key] = val

        if onekg:
            break

    return onekg


def extract_dbsnp_ids(colocated_variants):
    """
    Extract dbSNP IDs from VEP colocated_variants.

    Args:
        colocated_variants (list): VEP colocated_variants array

    Returns:
        str: Comma-separated dbSNP IDs, or empty string
    """
    if not colocated_variants:
        return ""
    rs_ids = []
    for cv in colocated_variants:
        cv_id = cv.get("id", "")
        if cv_id.startswith("rs"):
            rs_ids.append(cv_id)
    return ",".join(rs_ids)


def extract_clinvar_from_custom_annotations(custom_annotations, alt_allele):
    """
    Extract ClinVar information from VEP custom_annotations.

    VEP custom_annotations structure for ClinVar:
        "custom_annotations": {
            "ClinVar": [{
                "name": 132806,  # ClinVar Variation ID
                "allele": "A",
                "fields": {
                    "CLNSIG": "Likely_pathogenic",
                    "CLNREVSTAT": "reviewed_by_expert_panel",
                    "CLNDN": "Immunodeficiency_14",
                    "CLNVC": "single_nucleotide_variant"
                }
            }]
        }

    Args:
        custom_annotations (dict): VEP custom_annotations dictionary
        alt_allele (str): The alternate allele to match

    Returns:
        tuple: (clinvar_id, clinvar_significance, clinvar_review_status, clinvar_scv_count, clinvar_pathogenic_submitter_count)
            clinvar_scv_count is the number of SCV submissions (from CLNSIGSCV field)
            clinvar_pathogenic_submitter_count is the number of pathogenic/likely pathogenic submitters (from CLNSIGCONF field)
    """
    if not custom_annotations:
        return "", "", "", 0, 0

    clinvar_entries = custom_annotations.get("ClinVar", [])
    if not clinvar_entries:
        return "", "", "", 0, 0

    # Find the entry matching the alt allele
    best_entry = None
    for entry in clinvar_entries:
        entry_allele = entry.get("allele", "")
        if entry_allele == alt_allele:
            best_entry = entry
            break

    # If no allele match, use the first entry
    if best_entry is None and clinvar_entries:
        best_entry = clinvar_entries[0]

    if best_entry is None:
        return "", "", "", 0, 0

    # Extract the ClinVar ID (the "name" field is the variation ID as an integer)
    # Convert to VCV format to match Nirvana output (e.g., 132806 -> VCV000132806)
    raw_id = best_entry.get("name", "")
    if raw_id:
        clinvar_id = f"VCV{int(raw_id):09d}"
    else:
        clinvar_id = ""

    fields = best_entry.get("fields", {})

    # Extract and normalize significance
    # VEP format: "Likely_pathogenic" -> "likely pathogenic"
    raw_sig = fields.get("CLNSIG", "")
    clinvar_significance = raw_sig.replace("_", " ").lower() if raw_sig else ""

    # Extract and normalize review status
    # VEP format: "reviewed_by_expert_panel" -> "reviewed by expert panel"
    raw_review = fields.get("CLNREVSTAT", "")
    clinvar_review_status = raw_review.replace("_", " ").lower() if raw_review else ""

    # Extract SCV submission count from CLNSIGSCV field
    # CLNSIGSCV contains pipe-separated SCV IDs (e.g., "SCV000490382|SCV001135325|SCV001239025")
    # Each SCV represents an independent ClinVar submission
    clinvar_scv_count = 0
    raw_scv = fields.get("CLNSIGSCV", "")
    if raw_scv:
        clinvar_scv_count = len(raw_scv.split("|"))

    # Extract pathogenic submitter count from CLNSIGCONF field
    # CLNSIGCONF format: "Pathogenic(3)|Likely_benign(1)|Uncertain_significance(1)"
    # This gives us accurate per-classification counts even for conflicting variants
    clinvar_pathogenic_submitter_count = 0
    raw_conf = fields.get("CLNSIGCONF", "")
    if raw_conf:
        for entry in raw_conf.split("|"):
            entry_lower = entry.lower()
            if "pathogenic" in entry_lower:
                # Extract count from parentheses, e.g. "Pathogenic(3)" -> 3
                paren_start = entry.rfind("(")
                paren_end = entry.rfind(")")
                if paren_start != -1 and paren_end != -1:
                    try:
                        clinvar_pathogenic_submitter_count += int(entry[paren_start + 1:paren_end])
                    except ValueError:
                        pass

    return clinvar_id, clinvar_significance, clinvar_review_status, clinvar_scv_count, clinvar_pathogenic_submitter_count


def process_transcript(transcript, hgnc_to_gene_data):
    """
    Process a VEP transcript entry and extract gene-level information.

    Args:
        transcript (dict): A VEP transcript_consequences entry
        hgnc_to_gene_data (dict): Gene name to gene data mapping

    Returns:
        tuple: (gene_associated_disease, gene_name, consequence, clingen_gene_validity)
    """
    consequence = ",".join(transcript.get("consequence_terms", []))
    gene_name = transcript.get("gene_symbol", "n/a")

    gene_data = hgnc_to_gene_data.get(gene_name)
    if not gene_data:
        return "", gene_name, consequence, []

    clingen_gene_validity = gene_data.get("clingenGeneValidity", [])
    clingen_associated_disease = set()
    if clingen_gene_validity:
        for validity in clingen_gene_validity:
            disease = validity.get("disease", "").strip()
            if disease:
                clingen_associated_disease.add(disease)

    gene_associated_disease = sorted(list(clingen_associated_disease))
    return gene_associated_disease, gene_name, consequence, clingen_gene_validity


def process_variant(vep_entry, hgnc_to_gene_data):
    """
    Process a single VEP JSON entry and return a populated VariantData object.

    Args:
        vep_entry (dict): A parsed VEP JSON line
        hgnc_to_gene_data (dict): Gene name to gene data mapping

    Returns:
        VariantData: Populated variant object, or None if the entry cannot be parsed
    """
    # Parse alleles from allele_string (format: "G/T")
    allele_string = vep_entry.get("allele_string", "")
    alleles = allele_string.split("/")
    if len(alleles) < 2:
        logging.warning("Skipping VEP entry with unparseable allele_string: %s", allele_string)
        return None
    ref = alleles[0]
    alt = alleles[1]

    # Chromosome - ensure chr prefix
    chrom = vep_entry.get("seq_region_name", "")
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"

    position = str(vep_entry.get("start", ""))

    # Map VEP variant_class to Nirvana variantType
    vep_class = vep_entry.get("variant_class", "SNV")
    variant_type = VEP_VARIANT_CLASS_MAP.get(vep_class, vep_class)

    single_variant = VariantData(chrom, position, ref, alt, variant_type)

    # --- Colocated variants: gnomAD, 1000 Genomes, dbSNP ---
    colocated = vep_entry.get("colocated_variants", [])
    single_variant.gnomad = extract_gnomad_from_colocated(colocated, alt)
    single_variant.oneKg = extract_onekg_from_colocated(colocated, alt)
    single_variant.dbSnpIds = extract_dbsnp_ids(colocated)

    # Set default values for variants without transcript consequences
    single_variant.consequence = "n/a"
    single_variant.geneName = "n/a"

    # --- Transcript consequences ---
    transcript_list = vep_entry.get("transcript_consequences", [])
    if transcript_list:
        sorted_transcripts = rank_transcript_entries(transcript_list, hgnc_to_gene_data)
        single_variant.transcriptList = sorted_transcripts
        best = sorted_transcripts[0]
        single_variant.transcript = best.get("transcript_id", "")

        # HGVS notation
        if best.get("hgvsc"):
            single_variant.hgvsc = best["hgvsc"]
        if best.get("hgvsp"):
            single_variant.hgvsp = best["hgvsp"]
            single_variant.protein_variant = convert_vep_mutation_format(best["hgvsp"])

        # Gene and consequence from best transcript
        gene_associated_disease, gene_name, consequence, clingen_gene_validity = process_transcript(
            best, hgnc_to_gene_data
        )
        single_variant.geneName = gene_name
        single_variant.geneAssociatedDisease = gene_associated_disease
        single_variant.consequence = consequence
        single_variant.clingen_gene_validity = clingen_gene_validity

        # Gene-level gnomAD data (LOEUF, etc.)
        if hgnc_to_gene_data.get(gene_name):
            if hgnc_to_gene_data[gene_name].get("gnomAD"):
                single_variant.gene_gnomad = hgnc_to_gene_data[gene_name]["gnomAD"]

        # PolyPhen and SIFT predictions from transcript consequences
        if best.get("polyphen_score") is not None:
            single_variant.polyphen_score = best["polyphen_score"]
        if best.get("polyphen_prediction"):
            single_variant.polyphen_prediction = best["polyphen_prediction"]
        if best.get("sift_score") is not None:
            single_variant.sift_score = best["sift_score"]
        if best.get("sift_prediction"):
            single_variant.sift_prediction = best["sift_prediction"]

        # AlphaMissense predictions from VEP plugin
        alphamissense = best.get("alphamissense", {})
        if alphamissense:
            if alphamissense.get("am_pathogenicity") is not None:
                single_variant.alphamissense_score = alphamissense["am_pathogenicity"]
            if alphamissense.get("am_class"):
                single_variant.alphamissense_class = alphamissense["am_class"]

        # REVEL score from transcript consequences (VEP includes this with --everything flag)
        if best.get("revel") is not None:
            single_variant.revel = best["revel"]

    # --- Additional annotation fields ---
    # These may be present depending on VEP plugins/custom annotations.
    # Designed to be extended as new annotations are added to the VEP pipeline.
    # TODO: Implement this for VEP to get more BP4 calls
    #if vep_entry.get("phylopScore"):
    #    single_variant.phylopScore = vep_entry["phylopScore"]

    # --- Custom annotations: gnomAD counts and ClinVar ---
    custom_annotations = vep_entry.get("custom_annotations", {})
    if custom_annotations:
        # Extract gnomAD allele counts (AC, AN) from custom annotation
        # These are not available from colocated_variants which only has AFs
        gnomad_counts = extract_gnomad_from_custom_annotations(custom_annotations, alt)
        if gnomad_counts:
            # Merge count data into existing gnomad dict (which has AFs from colocated_variants)
            for key, value in gnomad_counts.items():
                if key not in single_variant.gnomad:
                    single_variant.gnomad[key] = value

        clinvar_id, clinvar_sig, clinvar_review, _, _ = extract_clinvar_from_custom_annotations(
            custom_annotations, alt
        )
        single_variant.clinvar_id = clinvar_id
        single_variant.clinvar_significance = clinvar_sig
        single_variant.clinvar_review_status = clinvar_review
        # PS4 pathogenic submitter counts are handled exclusively by the precomputed
        # file (loaded via bias_dataset_loader) and applied in get_ps4().

    return single_variant


def classify_variants(output_file, vep_json_file, hgnc_to_gene_data, name_to_dataset, variant_to_user_classification, skip_list):
    """
    Stream through a VEP JSON file, classify variants, and write results.

    VEP JSON files contain one JSON object per line, each representing a single variant.

    Args:
        output_file (str): Path to the output TSV file
        vep_json_file (str): Path to the VEP JSON input file
        hgnc_to_gene_data (dict): Gene-level data (ClinGen validity + gnomAD constraints) keyed by gene name
        name_to_dataset (dict): Dictionary containing all BIAS ACMG datasets
        variant_to_user_classification (dict): User-provided classifier overrides keyed by (chrom, pos, ref, alt)
        skip_list (list): List of ACMG codes to skip during classification

    Returns:
        int: Number of variants processed
    """
    v_count = 0
    try:
        with open(output_file, 'w') as o_file:
            headers = ["chromosome", "position", "refAllele", "altAllele", "variantType", "consequence",
                       "acmgClassification", "alleleFreq", "hgvsg", "hgvsc", "hgvsp", "aaChange", "geneName",
                       "pubmedIds", "associatedDiseases", "dbSnpids", "transcript", "rationale"]
            o_file.write("\t".join(headers) + "\n")

            with open_file(vep_json_file, "rt") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        vep_entry = json.loads(line)
                        single_variant = process_variant(vep_entry, hgnc_to_gene_data)
                        if single_variant is None:
                            continue

                        variant_key = (single_variant.chromosome, single_variant.position,
                                       single_variant.refAllele, single_variant.altAllele)
                        # Get supplemental user classification, if any
                        supplemental_codes = variant_to_user_classification.get(variant_key, {})

                        # Classify the variant
                        single_variant.significance, single_variant.justification = \
                            bias_variant_classification.get_variant_classification(
                                single_variant, name_to_dataset, supplemental_codes, skip_list)
                        v_count += 1
                        o_file.write(single_variant.to_tsv() + "\n")
                    except json.JSONDecodeError as e:
                        logging.error("Malformed JSON in VEP file: %s. Error: %s", line[:100], e)
                        continue
                    except KeyError as e:
                        logging.error("Missing key in VEP JSON data: %s. Error: %s", line[:100], e)
                        traceback.print_exc()
                        raise
                    except Exception as e:
                        logging.error("Unexpected error while processing VEP variant: %s. Error: %s", line[:100], e)
                        raise
                    if v_count % 25000 == 0:
                        logging.debug("Processed %d variants...", v_count)

        logging.info("Processed %d total variants from VEP file", v_count)
    except Exception as e:
        logging.critical("Error during VEP variant classification. Output file: %s. Error: %s", output_file, e)
        sys.exit(1)
    return v_count
