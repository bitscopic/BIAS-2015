#!/usr/bin/env python3
# Chris Eisenhart - chris.eisenhart@bitscopic.com
"""
Extract information from a variant line from the ICA/NIRVANA output

Each variant is stored in a VariantData class

The input file is streamed through as recommended by the NIRVANA team, each line is read and parsed individually
"""
import json
import logging
import sys
from src.bias_2015.constants import clinvar_review_status_to_level
from src.bias_2015.variant_data import VariantData, open_file
from src.bias_2015 import bias_variant_classification
from src.bias_2015.vcep_lookup import compute_popmax


def normalize_alleles(ref, alt):
    """
    Normalize VCF-style alleles to minimal representation.

    VCF uses padding (e.g., REF="TC", ALT="T" for a C deletion).
    Nirvana normalizes to minimal form (e.g., refAllele="C", altAllele="-").

    This function converts VCF-style to Nirvana-style for comparison.

    Args:
        ref (str): Reference allele
        alt (str): Alternate allele

    Returns:
        tuple: (normalized_ref, normalized_alt)
    """
    # Already normalized
    if ref == "-" or alt == "-":
        return ref, alt

    # Strip common prefix (VCF padding)
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]

    # Handle remaining padding - if one allele is now a single base that matches the start of the other
    if len(ref) > 1 and len(alt) == 1 and ref[0] == alt[0]:
        # Deletion: REF="TC" ALT="T" -> REF="C" ALT="-"
        ref = ref[1:]
        alt = "-"
    elif len(alt) > 1 and len(ref) == 1 and alt[0] == ref[0]:
        # Insertion: REF="T" ALT="TC" -> REF="-" ALT="C"
        alt = alt[1:]
        ref = "-"

    return ref, alt


def alleles_match(vcf_ref, vcf_alt, clinvar_ref, clinvar_alt):
    """
    Check if VCF alleles match ClinVar alleles, accounting for normalization differences.

    Args:
        vcf_ref (str): Reference allele from VCF/variant
        vcf_alt (str): Alternate allele from VCF/variant
        clinvar_ref (str): Reference allele from ClinVar entry
        clinvar_alt (str): Alternate allele from ClinVar entry

    Returns:
        bool: True if alleles match after normalization
    """
    # Direct match
    if vcf_ref == clinvar_ref and vcf_alt == clinvar_alt:
        return True

    # Normalize VCF alleles and compare
    norm_ref, norm_alt = normalize_alleles(vcf_ref, vcf_alt)
    if norm_ref == clinvar_ref and norm_alt == clinvar_alt:
        return True

    return False


def rank_clinvar_entries(entries):
    """
    Ranks a list of ClinVar entries based on predefined criteria.

    Args:
        entries (list): A list of JSON objects representing ClinVar entries.

    Returns:
        list: A new list of ClinVar entries sorted in descending order based on rank scores.

    Criteria:
    - Classification: Entries with pathogenic or likely pathogenic significance are ranked higher.
    - Submitters: Entries with multiple submitters are ranked higher.
    - Supporting PubMed IDs: Entries with more PubMed IDs are ranked higher.
    - Expert Panel review: Entries reviewed by ClinVar's Expert Panel are ranked higher.
    - Phenotypes: Entries with more associated phenotypes are ranked higher.
    - Relevant Phenotype: Entries associated with relevant phenotypes (e.g., Noonan syndrome) are ranked higher.
    - Last Updated Date: Entries with the most recent lastUpdatedDate are ranked higher.

    Note:
    - This function assumes that the provided JSON objects contain all the necessary fields.
    - The function calculates rank scores based on the defined criteria and returns a new list of entries sorted in descending order based on the rank scores.

    Example:
    >>> entries = [{'id': 'RCV000038275.6', 'variationId': 'VCV000045127.5', ...}, {...}, ...]
    >>> ref_allele = 'T'
    >>> alt_allele = 'A'
    >>> ranked_entries = rank_clinvar_entries(entries, ref_allele, alt_allele)
    """
    # Initialize a dictionary to store the rank scores for each entry
    rank_scores = {}

    # Iterate through each entry and assign a rank score based on the criteria
    for entry in entries:
        # Initialize the rank score for the current entry
        rank_scores[entry["id"]] = clinvar_review_status_to_level.get(entry["reviewStatus"], 0)

        # Increase rank score if the classification is likely pathogenic
        if "likely pathogenic" in entry["significance"]:
            rank_scores[entry["id"]] += 1

        # Increase rank score if the classification is pathogenic
        if "pathogenic" in entry["significance"]:
            rank_scores[entry["id"]] += 2

        # Increase rank score if there are supporting PubMed IDs
        if "pubMedIds" in entry and entry["pubMedIds"]:
            rank_scores[entry["id"]] += len(entry["pubMedIds"])

        # Increase rank score based on the number of phenotypes associated with the variant
        if "phenotypes" in entry and len(entry["phenotypes"]) > 0:
            rank_scores[entry["id"]] += len(entry["phenotypes"])

        # Increase rank score based on the relevance of the phenotype to your context (e.g., Noonan syndrome)
        if "phenotypes" in entry and "Noonan syndrome" in entry["phenotypes"]:
            rank_scores[entry["id"]] += 1

        # Increase rank score if the entry has the most recent lastUpdatedDate
        if entry["lastUpdatedDate"] == max(e["lastUpdatedDate"] for e in entries):
            rank_scores[entry["id"]] += 1

    # Sort the entries based on their rank scores (in descending order)
    sorted_entries = sorted(entries, key=lambda e: rank_scores[e["id"]], reverse=True)

    return sorted_entries


def return_clinvar_significance(clin_var_sig_list):
    """
    Return the most high rates (pathogenic) clinvar significance seen
    """
    valid_significance = set()
    for sig in clin_var_sig_list:
        valid_significance.add(sig)
    most_sig = ""
    for element in valid_significance:
        if element.lower() == "pathogenic":
            return "pathogenic"
        if element.lower() == "likely pathogenic":
            most_sig = "likely pathogenic"
        if element.lower() == "uncertain" and most_sig != "likely pathogenic":
            most_sig = "uncertain"
        if element.lower() == "likely benign" and most_sig != "uncertain" and  most_sig != "likely pathogenic":
            most_sig = "likely benign"
        if element.lower() == "benign" and most_sig != "likely benign" and most_sig != "uncertain" and  most_sig != "likely pathogenic":
            most_sig = "benign"
    return most_sig

def identify_clinvar_information(clin_var_list, ref_allele, alt_allele):
    """
    Identifies ClinVar information for a variant based on the reference allele and alternate allele.

    Args:
        clin_var_list (list): A list of ClinVar elements containing information about the variants.
        ref_allele (str): The reference allele of the variant.
        alt_allele (str): The alternate allele of the variant.

    Returns:
        tuple or "": A tuple containing the identified ClinVar information if it meets the criteria,
        or "" if there is no ClinVar data, or the variant does not meet the criteria.
        The tuple contains the following elements in order:
        - significance (str): The clinically relevant significances associated with the variant.
        - pubmed_joined_ids (str): The PubMed IDs associated with the variant and ClinVar entry.

    Raises:
        "".
    """
    # Go through the clinvar elements and remove variants that dont have the same ref and alt allele as the observed variant
    clean_clinvar_list = []
    clinvar_id = ""
    for clin_var in clin_var_list:
        if 'VCV' in clin_var['id']:
            clinvar_id = clin_var['id']
        if clin_var.get('variationId'):
            if 'VCV' in clin_var['variationId']:
                clinvar_id = clin_var['variationId']
        # Only consider variants that share the same reference base and alt allele
        # Use alleles_match to handle VCF padding vs Nirvana normalization differences
        clinvar_ref = clin_var.get("refAllele", "")
        clinvar_alt = clin_var.get("altAllele", "")
        if not alleles_match(ref_allele, alt_allele, clinvar_ref, clinvar_alt):
            continue
        # Skip entries without significance (e.g. "no assertion provided" submissions)
        if "significance" not in clin_var:
            continue
        clean_clinvar_list.append(clin_var)

    # Corner case where the none of the clinVar entries shared the same alt allele as current variant
    if not clean_clinvar_list:
        return "", "", "", "", clinvar_id

    # Identify the best clinvar entry using a ranking schema
    best_clin_var_entry = rank_clinvar_entries(clean_clinvar_list)[0]

    # Clinvar variant id with format 'RCV000005170.4'
    variant_id = best_clin_var_entry['id']

    # Identify the clinically relevant significances
    valid_significance = []
    ignored_significance_set = set(['benign', 'likely benign'])
    for sig in best_clin_var_entry["significance"]:
        if sig in ignored_significance_set:
            continue
        valid_significance.append(sig)

    # Harvest the clinVar significance
    significance = ""
    if best_clin_var_entry.get('significance'):
        significance = return_clinvar_significance(best_clin_var_entry["significance"])

    # Gather the pubmed ID's associated with this variant and clinvar entry
    pubmed_joined_ids = best_clin_var_entry.get("pubMedIds","")

    # Supporting evidence level
    review_status = best_clin_var_entry.get("reviewStatus", "")
    return variant_id, significance, pubmed_joined_ids, review_status, clinvar_id


def process_transcript(transcript, hgnc_to_gene_data):
    """
    Processes a transcript from a Nirvana JSON file and extracts relevant information.

    Args:
        transcript (dict): The transcript dictionary containing information about the transcript.
        hgnc_to_gene_data (dict): A dictionary mapping HGNC symbols to gene data.

    Returns:
        tuple or "": A tuple containing the processed transcript information if it meets the criteria,
        or "" if the transcript does not have sufficient annotations or does not meet the criteria.
        The tuple contains the following elements in order:
        - clingen_associated_disease (set): A set of unique diseases associated with the transcript.
        - valid_phenotype_set (set): A set of valid phenotypes associated with the transcript.
        - gene_name (str): The gene name associated with the transcript.
        - consequence (str): The consequences of the transcript.

    Raises:
        "".
    """
    # Gather the consequences as a string
    consequence = ",".join(transcript["consequence"])

    # Access gene data fields
    gene_name = transcript["hgnc"]

    # Verify that gene data was reported by Nirvana for this gene
    gene_data = hgnc_to_gene_data.get(gene_name)
    if not gene_data:
        return "", gene_name, consequence, []

    # The gene data must include clingenGeneValidity
    clingen_gene_validity = gene_data.get("clingenGeneValidity")
    # Gather the unique diseases seen in the clingenGeneValidity data
    clingen_associated_disease = set()
    if clingen_gene_validity:
        for validity in clingen_gene_validity:
            disease = validity["disease"].strip()
            clingen_associated_disease.add(disease)

    # Format for simplicity
    gene_associated_disease = sorted(list(clingen_associated_disease))

    return gene_associated_disease, gene_name, consequence, clingen_gene_validity


def convert_mutation_format(mutation_str):
    """
    Convert mutation strings to the shorthand version
        NM_000059.3:c.4563A>G(p.(Leu1521=))
        L1521*
        NP_000050.2:p.(Lys1691AsnfsTer15)
        K1691fsN15*
        NM_000059.3:c.6513G>C(p.(Val2171=))
        V2171*
    """
    # Throw away the exceptionally weird ones
    if "?" in mutation_str:
        return ""
    
    # Parse out the protein string within the parenthesis
    if "(p.(" in mutation_str:
        protein_change = mutation_str.split("(")[2][:-2]
    else:
        protein_change = mutation_str.split("(")[1][:-1]

    # Mapping table for three-letter to single-letter amino acids
    amino_acid_mapping = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q',
        'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F',
        'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Ter': '*'
    }

    # The first AA is always at the start of the protein string
    first_aa = amino_acid_mapping.get(protein_change[:3])

    if not first_aa:
        return ""

    # Go through the remaining protein string and parse out the final AA (if it exists)
    # and the number of bases
    converted_str = first_aa
    if amino_acid_mapping.get(protein_change[-3:]): # Basic change of one AA for another
        base_number = protein_change[3:-3]
        second_aa = amino_acid_mapping[protein_change[-3:]]
        converted_str = f"{first_aa}{base_number}{second_aa}"
    elif protein_change.endswith("="): # Stop codon
        base_number = protein_change[3:-1]
        converted_str = f"{first_aa}{base_number}*"
    elif "fs" in protein_change:
        term_num = protein_change.split("Ter")[1]
        pre_term_protein = protein_change.split("Ter")[0]
        base_number = pre_term_protein[3:-5]
        second_aa = amino_acid_mapping.get(pre_term_protein[-5:-2], "")
        converted_str = f"{first_aa}{base_number}fs{second_aa}{term_num}*"
    return converted_str


def rank_transcript_entries(transcript_list, hgnc_to_gene_data, transcript_database):
    """
    Identifies the ideal transcript from a list of transcripts based on predefined criteria.

    Args:
        transcript_list (list): A list of transcript dictionaries containing information about the transcripts.
        hgnc_to_gene_data (dict): A dictionary mapping HGNC symbols to gene data.
        transcript_database (str): The transcript database to use for identifying the ideal transcript.

    Returns:
        list: A new list of transcript dictionaries sorted in descending order based on rank scores.

    Criteria:
    - Is Canonical: Transcripts flagged as canonical are ranked higher.
    - Protein-Coding: Transcripts with bioType "mRNA" are ranked higher.
    - Has Coding Sequence Information: Transcripts with "cdsPos" or "proteinPos" are ranked higher.
    - Transcript Database: Transcripts from the specified transcript database are ranked higher.
    - HGNC Data: Transcripts with available HGNC data are ranked higher.

    Example:
    >>> transcript_list = [{'transcript': 'ENST000001', 'isCanonical': True, 'source': 'Ensembl', 'hgnc': 'GENE1'}, {...}, ...]
    >>> hgnc_to_gene_data = {'GENE1': {...}, ...}
    >>> transcript_database = 'Ensembl'
    >>> sorted_transcripts = rank_transcript_entries(transcript_list, hgnc_to_gene_data, transcript_database)
    """
    # Initialize a dictionary to store the rank scores for each entry
    rank_scores = {}

    # Iterate through each entry and assign a rank score based on the criteria
    for entry in transcript_list:
        transcript_id = entry["transcript"]
        rank_scores[transcript_id] = 0

        # Increase rank for canonical transcripts
        if entry.get("isCanonical"):
            rank_scores[transcript_id] += 3  # High priority

        # Prioritize protein-coding transcripts
        if entry.get("bioType") == "mRNA":
            rank_scores[transcript_id] += 2  # Medium priority

        # Prioritize transcripts with coding sequence information
        if "cdsPos" in entry or "proteinPos" in entry:
            rank_scores[transcript_id] += 2  # Medium priority

        # Increase rank if the transcript source matches the specified database
        if entry.get("source") in transcript_database:
            rank_scores[transcript_id] += 1

        # Increase rank if the transcript has HGNC mapping
        if entry.get("hgnc") in hgnc_to_gene_data:
            rank_scores[transcript_id] += 1

    # Sort the entries based on their rank scores (in descending order)
    sorted_entries = sorted(transcript_list, key=lambda e: rank_scores[e["transcript"]], reverse=True)
    return sorted_entries



def process_variant(variant, hgnc_to_gene_data, transcript_database, chrom, position, ref, alt):
    """
    Processes a variant from a Nirvana JSON file and extracts relevant information.

    Args:
        variant (dict): The variant dictionary containing information about the variant.
        hgnc_to_gene_data (dict): A dictionary mapping HGNC symbols to gene data.
        transcript_database (str): The transcript database to use for processing variants.
        chrom (str): The chromosome
        position (str): The position

    Returns:
        dict or "": A dictionary containing the processed variant information if it meets the criteria,
        or "" if the variant does not have sufficient annotations or does not meet the criteria.

    Raises:
        "".
    """
    # This sadly isnt standardized! Some ommit the chr, others require it. We ensure it is always there. 
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    
    # Assign essential variant elements
    single_variant = VariantData(chrom, position, ref, alt, variant["variantType"])

    if variant.get("phylopScore"):
        single_variant.phylopScore = variant["phylopScore"]

    # Grab the hgvsg name if it is present
    single_variant.hgvsg = variant.get("hgvsg", "")

    # Gather relvent dbsnp information on this variant
    if variant.get("dbsnp"):
        db_snp_ids = ",".join(variant["dbsnp"])
        single_variant.dbSnpIds = db_snp_ids

    if variant.get("gnomad"):
        single_variant.gnomad = variant.get("gnomad")
        # Compute popmax from per-population AFs so VCEP rules that specify
        # popmax/grpmax/FAF95 have a metric to evaluate against.
        popmax = compute_popmax(single_variant.gnomad)
        if popmax is not None:
            single_variant.gnomad["popmax"] = popmax
    else:
        single_variant.gnomad = {}

    if variant.get("oneKg"):
        single_variant.oneKg = variant['oneKg']
    else:
        single_variant.oneKg = {}

    if variant.get("revel"):
        single_variant.revel = variant['revel']['score']

    if variant.get("dannScore"):
        single_variant.dann = variant['dannScore']

    if variant.get("gerpScore"):
        single_variant.gerp = variant['gerpScore']

    # AlphaMissense score - verify alleles match before using
    if variant.get("AlphaMissense"):
        am_data = variant["AlphaMissense"]
        am_ref = am_data.get("refAllele", "")
        am_alt = am_data.get("altAllele", "")
        if am_ref == single_variant.refAllele and am_alt == single_variant.altAllele:
            single_variant.alphamissense_score = am_data.get("AM_score", "")

    # Gather the relevant clinvar information for this variant
    if variant.get("clinvar"):
        variant_id, significance, pubmed_joined_ids, review_status, clinvar_id = \
                identify_clinvar_information(variant["clinvar"], single_variant.refAllele, single_variant.altAllele)
        single_variant.variantId = variant_id
        single_variant.pubmedIds = pubmed_joined_ids
        single_variant.clinvar_significance = significance
        single_variant.clinvar_review_status = review_status
        single_variant.clinvar_id = clinvar_id
        # Count unique submitters classifying this variant as pathogenic/likely pathogenic
        # Each RCV entry represents one submitter's assertion
        patho_count = 0
        for cv in variant["clinvar"]:
            if not cv["id"].startswith("RCV"):
                continue
            cv_ref = cv.get("refAllele", "")
            cv_alt = cv.get("altAllele", "")
            if not alleles_match(single_variant.refAllele, single_variant.altAllele, cv_ref, cv_alt):
                continue
            for sig in cv.get("significance", []):
                if "pathogenic" in sig.lower():
                    patho_count += 1
                    break
        single_variant.clinvar_pathogenic_submitter_count = patho_count

    # Harvest transcript specific information
    transcript_list = variant.get("transcripts")
    if transcript_list:
        # Identify ideal transcript
        sorted_transcript_list = rank_transcript_entries(transcript_list, hgnc_to_gene_data, transcript_database)
        single_variant.transcriptList = sorted_transcript_list
        best_transcript = sorted_transcript_list[0] 
        single_variant.transcript = best_transcript['transcript']
        # Harvest information from the best transcript
        if best_transcript.get('hgvsc'):
            single_variant.hgvsc = best_transcript['hgvsc']
        if best_transcript.get('hgvsp'):
            single_variant.hgvsp = best_transcript['hgvsp']
            single_variant.protein_variant = convert_mutation_format(best_transcript['hgvsp'])
       
        gene_associated_disease, gene_name, consequence, clingen_gene_validity = process_transcript(
            best_transcript, hgnc_to_gene_data
        )
        single_variant.geneName = gene_name
        single_variant.geneAssociatedDisease = gene_associated_disease
        single_variant.consequence = consequence
        single_variant.clingen_gene_validity = clingen_gene_validity
        if hgnc_to_gene_data.get(gene_name):
            if hgnc_to_gene_data[gene_name].get('gnomAD'): 
                single_variant.gene_gnomad = hgnc_to_gene_data[gene_name]['gnomAD']
        if variant.get("topmed"):
            single_variant.topmed = variant['topmed']
    return single_variant


def classify_variants(output_file, nirvana_json_file, hgnc_to_gene_data, name_to_dataset, variant_to_user_classification, skip_list):
    """
    Stream through a Nirvana JSON file, classify variants, and write results.

    This function processes a Nirvana-annotated JSON file, applying ACMG classification
    to each variant and writing the results to a TSV output file.

    Args:
        output_file (str): Path to the output TSV file
        nirvana_json_file (str): Path to the Nirvana JSON input file
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

            with open_file(nirvana_json_file, "rt") as f:
                _ = f.readline()[10:-15]
                genes_section = False
                for line in f:
                    if genes_section:
                        break
                    if '"genes":[' in line:
                        genes_section = True
                        continue
                    if line.startswith("]}"):
                        continue
                    try:
                        data = json.loads(line[:-2] if line[-2] == "," else line)
                        chrom = data["chromosome"]
                        pos = str(data["position"])
                        ref = data["refAllele"]
                        if not data.get("altAlleles"):
                            continue  # Minor reference alleles - not supported currently
                        alt = data["altAlleles"][0]
                        for variant in data["variants"]:
                            single_variant = process_variant(
                                variant, hgnc_to_gene_data, "RefSeq", chrom, pos, ref, alt)
                            variant_key = (chrom, pos, single_variant.refAllele, single_variant.altAllele)
                            # Get supplemental user classification, if any
                            supplemental_codes = variant_to_user_classification.get(variant_key, {})

                            # Append the classification to the variant object (significance, justification)
                            single_variant.significance, single_variant.justification = \
                                bias_variant_classification.get_variant_classification(
                                    single_variant, name_to_dataset, supplemental_codes, skip_list)
                            v_count += 1
                            o_file.write(single_variant.to_tsv() + "\n")
                    except KeyError as e:
                        logging.error("Missing key in JSON data: %s. Error: %s", line.strip(), e)
                        raise
                    except Exception as e:
                        logging.error("Unexpected error while processing variant: %s. Error: %s", line.strip(), e)
                        raise
                    if v_count % 25000 == 0:
                        logging.debug("Processed %d variants...", v_count)

        logging.info("Processed %d total variants in %d genes", v_count, len(hgnc_to_gene_data))
    except Exception as e:
        logging.critical("Error during variant classification. Output file: %s. Error: %s", output_file, e)
        sys.exit(1)
    return v_count
