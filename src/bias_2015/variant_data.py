#!/usr/bin/env python3
# Chris Eisenhart - chris.eisenhart@bitscopic.com
"""
VariantData class for storing annotated variant information.

This module contains the core data structure used by BIAS-2015 to represent
a single variant with all its associated annotations. This class is shared
between different annotator ingestion modules (Nirvana, VEP, etc.).
"""
import json
import gzip
import os
import io


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
        self.polyphen_score = ""
        self.polyphen_prediction = ""
        self.sift_score = ""
        self.sift_prediction = ""
        self.topmed = {}
        self.phylopScore = ""
        self.clingen_gene_validity = []
        self.gene_gnomad = {}
        self.clinvar_id = ""
        self.clinvar_pathogenic_submitter_count = 0  # Number of ClinVar submitters classifying as pathogenic/likely pathogenic
        self.alphamissense_score = ""  # Continuous score 0-1
        self.alphamissense_class = ""  # benign, ambiguous, pathogenic

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

    Args:
        file_path (str): Path to the file to open
        mode (str): File mode (e.g., 'rt' for read text, 'wt' for write text)

    Returns:
        File handle for the opened file
    """
    _, file_extension = os.path.splitext(file_path)
    if file_extension == ".gz":
        return gzip.open(file_path, mode)
    return io.open(file_path, mode, encoding="utf-8")
