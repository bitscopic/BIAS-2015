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
import gzip
import os
import io
from src.bias_2015.constants import clinvar_review_status_to_level

class VariantData:
    """
    Helper class for storing annotated variant information for a single variant
    """
    def __init__(self, chromosome, position, ref_allele, alt_allele, variant_type):
        self.chromosome = chromosome
        self.position = position
        self.refAllele = ref_allele
        self.altAllele = alt_allele
        self.variantType = variant_type
        self.alleleFreq = ""
        self.variantReads = ""
        self.wildtypeReads = ""
        self.hgvsg = "n/a"
        self.hgvsc = "n/a"
        self.hgvsp = "n/a"
        self.protein_variant = "n/a" # Gene dependent, ex W515K
        self.geneName = "n/a" # The gene name
        self.consequence = "n/a" # Frameshift, missense, stop codon, etc
        self.significance = "uncertain" # The AMP/ACMG variant interpretation 
        self.justification = {} # AMP/ACMG Variant interpretation rationale and explanation
        self.clinvar_significance = ""
        self.pubmedIds = ""
        self.geneAssociatedDisease = ""
        self.dbSnpIds = ""
        self.transcript = ""
        self.variantId = ""
        self.transcriptList = ""
        self.gnomad = {}
        self.clinvar_review_status = ""
        self.oneKg = ""
        self.domain = ""
        self.gerp = ""
        self.dann = ""
        self.revel = ""
        self.topmed = {}
        self.phylopScore = ""
        self.clingen_gene_validity = []
        self.gene_gnomad = {}
        self.clinvar_id = ""

    def to_json(self):
        """
        Return a json format of the variant annotation class
        """
        return {
            "chromosome": self.chromosome,
            "position": self.position,
            "refAllele": self.refAllele,
            "altAllele": self.altAllele,
            "variantType": self.variantType,
            "alleleFrequency": self.alleleFreq,
            "variantReads": self.variantReads,
            "wildtypeReads": self.wildtypeReads,
            "hgvsg": self.hgvsg,
            "hgvsc": self.hgvsc,
            "hgvsp": self.hgvsp,
            "pdot": self.protein_variant,
            "geneName": self.geneName,
            "consequence": self.consequence,
            "significance": self.significance,
            "geneAssociatedDisease": self.geneAssociatedDisease,
            "dbSnpIds": self.dbSnpIds,
            "transcript": self.transcript,
            "variantId": self.variantId,
            "annotations": [
                {
                    "name": "ACMG Rationale",
                    "value": self.justification
                },
                {
                    "name": "pubmedIds",
                    "value": self.pubmedIds,
                },
                {
                    "name": "gnomad",
                    "value":
                        {
                        "alleleFrequency": self.gnomad.get('allAf', 0),
                        "coverage": self.gnomad.get('coverage', 0)
                        }
                },
                {
                    "name": "uniprot id",
                    "value": self.domain
                }
                ]
        }

    def to_tsv(self):
        """
        Return a tsv format of the variant annotation class
        """
        return "\t".join([self.chromosome,
            self.position,
            self.refAllele,
            self.altAllele,
            self.variantType,
            self.consequence,
            self.significance,
            self.alleleFreq,
            self.hgvsg,
            self.hgvsc,
            self.hgvsp,
            self.protein_variant,
            self.geneName,
            ",".join(self.pubmedIds),
            ",".join(self.geneAssociatedDisease),
            self.dbSnpIds,
            self.transcript,
            json.dumps(self.justification)])

    def __str__(self):
        return f"{self.chromosome}-{self.position}-{self.refAllele}-{self.altAllele}"

def open_file(file_path, mode):
    """
    Open either a normal or a .gz file
    """
    _, file_extension = os.path.splitext(file_path)
    if file_extension == ".gz":
        return gzip.open(file_path, mode)
    return io.open(file_path, mode, encoding="utf-8")

def load_nirvana_gene_information(nirvana_json_file):
    """
    NIRVANA gene information extraction with error handling.

    NIRVANA puts gene information at the end of the .json file, so you have to go through the
    file once to have it on hand for each variant
    """
    hgnc_to_gene_data = {}
    try:
        with open_file(nirvana_json_file, "rt") as f:
            _ = f.readline()[10:-15]
            genes_section = False
            for line in f:
                # Check if we have reached the end of the genes section
                if genes_section:
                    if line.strip() == "]}," or line.strip() == "]}":
                        break
                    try:
                        # Extract gene information
                        data = json.loads(line[:-2] if line[-2] == "," else line)
                        hgnc_id = data["name"]
                        hgnc_to_gene_data[hgnc_id] = data
                    except json.JSONDecodeError as e:
                        logging.error("Malformed JSON in line: %s. Error: %s", line.strip(), e)
                        continue
                # Check if we have reached the start of the genes section
                if '"genes":[' in line:
                    genes_section = True
    except Exception as e:
        logging.critical("Failed to load Nirvana gene information from file: %s. Error: %s", nirvana_json_file, e)
        sys.exit(1)
    return hgnc_to_gene_data


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
        if clin_var.get("refAllele", "") != ref_allele or clin_var.get("altAllele", "") != alt_allele:
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

    # Gather the relevant clinvar information for this variant
    if variant.get("clinvar"):
        variant_id, significance, pubmed_joined_ids, review_status, clinvar_id = \
                identify_clinvar_information(variant["clinvar"], single_variant.refAllele, single_variant.altAllele)
        single_variant.variantId = variant_id
        single_variant.pubmedIds = pubmed_joined_ids
        single_variant.clinvar_significance = significance
        single_variant.clinvar_review_status = review_status
        single_variant.clinvar_id = clinvar_id

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


def write_tsv(outFile, full_tsv_output, reference, creation_time):
    """
    Write out a tsv output with variant information
    """
    header = [
        'Chromosome',
        'Position',
        'Reference allele',
        'Alternative allele',
        'Variant type',
        'Allele freq',
        'Hgvsg',
        'Hgvsc',
        'Hgvsp',
        'Protein variant',
        'Gene',
        'Consequence',
        'Significance (ClinVar)',
        'Significance (Bitscopic)',
        'Consolidated Significance',
        'PubMed IDs',
        'Gene-Associated Disease',
        'Significant',
        'dbSNP IDs',
        'Transcript',
        'Variant ID',
    ]
    with open(outFile, "w") as oFile:
        oFile.write(f"#{reference}\n")
        oFile.write(f"#{creation_time}\n")
        oFile.write("#" + "\t".join(header) + "\n")
        oFile.write("\n".join(full_tsv_output))
