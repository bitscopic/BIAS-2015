"""
Unit testing for the Bitscopic Interpreting ACMG Standards (BIAS) pathogenic classifiers
"""
import os
import sys

# NOTE: This should be handled more elegantly for public facing code
# Some hacky path adjustments so users dont have to fiddle with env variables
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from src.bias_2015 import pathogenic_classifiers, extract_from_nirvana_json

clinvar_review_status_to_level = {
    "practice guideline": 5,
    "reviewed by expert panel": 5,
    "criteria provided, multiple submitters, no conflicts": 3,
    "criteria provided, single submitter": 2,
    "criteria provided, conflicting interpretations": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0
}

def test_get_pvs1():
    """
    Test the very strong pathogenic classifier reporting a PVS1
    """
    chromosome = 'chr2'
    position = '25470531'
    ref_allele = 'TC'
    alt_allele = 'AT'
    variant_type = 'MNV'
    variant = extract_from_nirvana_json.VariantData(chromosome, position, ref_allele, alt_allele, variant_type)
    variant.consequence = 'stop_gained'
    variant.geneName = 'DNMT3A'
    variant.transcriptList =  [{'transcript': 'NM_022552.5', 'source': 'RefSeq', 'bioType': 'mRNA', 'codons': 'tgGAtg/tgATtg',
                                'aminoAcids': 'WM/*', 'cdnaPos': '1219-1220', 'cdsPos': '942-943', 'exons': '8/23',
                                'proteinPos': '314-315', 'geneId': '1788', 'hgnc': 'DNMT3A', 'consequence': ['stop_gained'],
                                'hgvsc': 'NM_022552.5:c.942_943delinsAT', 'hgvsp': 'NP_072046.2:p.(Trp314Ter)', 
                                'isCanonical': True, 'proteinId': 'NP_072046.2'}, {'transcript': 'NR_135490.2', 'source':
                                'RefSeq', 'bioType': 'transcript', 'cdnaPos': '1173-1174', 'exons': '8/24', 'geneId':
                                '1788', 'hgnc': 'DNMT3A', 'consequence': ['stop_gained'], 'hgvsc':
                                'NR_135490.2:n.1173_1174delinsAT'}]
    lof_gene_to_pli = {'DNMT3A': '.99'}
    gene_name_to_3prime_region = {'DNMT3A': ('chr2', 25455829, 25455879)}
    pvs_score, pvs_rationale = pathogenic_classifiers.get_pvs(variant, lof_gene_to_pli, gene_name_to_3prime_region)
    assert pvs_score == 1
    assert pvs_rationale == [(1, 'PSV1 (1): Null variant consequence stop_gained in gene DNMT3A with GNOMAD PLI .99')]


def test_get_pvs1_caveat1():
    """
    Test the very strong pathogenic classifier caveat one dropout, when the variant consequence isn't significant enough
    to support loss of function.
    """
    chromosome = 'chr2'
    position = '25470531'
    ref_allele = 'TC'
    alt_allele = 'AT'
    variant_type = 'MNV'
    variant = extract_from_nirvana_json.VariantData(chromosome, position, ref_allele, alt_allele, variant_type)
    variant.consequence = 'synonymous'
    variant.geneName = 'DNMT3A'
    variant.transcriptList =  [{'transcript': 'NM_022552.5', 'source': 'RefSeq', 'bioType': 'mRNA', 'codons': 'tgGAtg/tgATtg',
                                'aminoAcids': 'WM/*', 'cdnaPos': '1219-1220', 'cdsPos': '942-943', 'exons': '8/23',
                                'proteinPos': '314-315', 'geneId': '1788', 'hgnc': 'DNMT3A', 'consequence': ['stop_gained'],
                                'hgvsc': 'NM_022552.5:c.942_943delinsAT', 'hgvsp': 'NP_072046.2:p.(Trp314Ter)', 
                                'isCanonical': True, 'proteinId': 'NP_072046.2'}, {'transcript': 'NR_135490.2', 'source':
                                'RefSeq', 'bioType': 'transcript', 'cdnaPos': '1173-1174', 'exons': '8/24', 'geneId':
                                '1788', 'hgnc': 'DNMT3A', 'consequence': ['stop_gained'], 'hgvsc':
                                'NR_135490.2:n.1173_1174delinsAT'}]
    lof_gene_list = ['DNMT3A']
    gene_name_to_3prime_region = {'DNMT3A': ('chr2', 25455829, 25455879)}
    pvs_score, pvs_rationale = pathogenic_classifiers.get_pvs(variant, lof_gene_list, gene_name_to_3prime_region)
    assert pvs_score == 0
    assert pvs_rationale == [(0, 'synonymous is not a null consequence')]


def test_get_pvs1_caveat2():
    """
    Test the very strong pathogenic classifier caveat two dropout, variant is at the extreme 3' end of a gene.

    This is testing that the region filtering is applied correctly, please see test_pathogenic_classifiers_dataset_loader
    to see the tests covering if the regions are generated correctly.
    """
    chromosome = 'chr2'
    position = '25455849'
    ref_allele = 'TC'
    alt_allele = 'AT'
    variant_type = 'MNV'
    variant = extract_from_nirvana_json.VariantData(chromosome, position, ref_allele, alt_allele, variant_type)
    variant.consequence = 'stop_gained'
    variant.geneName = 'DNMT3A'
    variant.transcriptList =  [{'transcript': 'NM_022552.5', 'source': 'RefSeq', 'bioType': 'mRNA', 'codons': 'tgGAtg/tgATtg',
                                'aminoAcids': 'WM/*', 'cdnaPos': '1219-1220', 'cdsPos': '942-943', 'exons': '8/23',
                                'proteinPos': '314-315', 'geneId': '1788', 'hgnc': 'DNMT3A', 'consequence': ['stop_gained'],
                                'hgvsc': 'NM_022552.5:c.942_943delinsAT', 'hgvsp': 'NP_072046.2:p.(Trp314Ter)', 
                                'isCanonical': True, 'proteinId': 'NP_072046.2'}, {'transcript': 'NR_135490.2', 'source':
                                'RefSeq', 'bioType': 'transcript', 'cdnaPos': '1173-1174', 'exons': '8/24', 'geneId':
                                '1788', 'hgnc': 'DNMT3A', 'consequence': ['stop_gained'], 'hgvsc':
                                'NR_135490.2:n.1173_1174delinsAT'}]
    lof_gene_list = ['DNMT3A']
    gene_name_to_3prime_region = {'DNMT3A': ('chr2', 25455829, 25455879)}
    pvs_score, pvs_rationale = pathogenic_classifiers.get_pvs(variant, lof_gene_list, gene_name_to_3prime_region)
    assert pvs_score == 0
    assert pvs_rationale == [(0, "Position 25455849 is in the extreme three prime region (25455829-25455879) of the last coding exon in gene DNMT3A")]


def test_get_pvs1_caveat4():
    """
    Test the very strong pathogenic classifier caveat four dropout, use caution when the variant is associated with multiple transcripts

    If there are multiple transcripts, ensure the variant consequence is seen more than one time. 
    """
    chromosome = 'chr2'
    position = '25470531'
    ref_allele = 'TC'
    alt_allele = 'AT'
    variant_type = 'MNV'
    variant = extract_from_nirvana_json.VariantData(chromosome, position, ref_allele, alt_allele, variant_type)
    variant.consequence = 'stop_gained'
    variant.geneName = 'DNMT3A'
    variant.transcriptList =  [{'transcript': 'NM_022552.5', 'source': 'RefSeq', 'bioType': 'mRNA', 'codons': 'tgGAtg/tgATtg',
                                'aminoAcids': 'WM/*', 'cdnaPos': '1219-1220', 'cdsPos': '942-943', 'exons': '8/23',
                                'proteinPos': '314-315', 'geneId': '1788', 'hgnc': 'DNMT3A', 'consequence': ['stop_gained'],
                                'hgvsc': 'NM_022552.5:c.942_943delinsAT', 'hgvsp': 'NP_072046.2:p.(Trp314Ter)', 
                                'isCanonical': True, 'proteinId': 'NP_072046.2'}, {'transcript': 'NR_135490.2', 'source':
                                'RefSeq', 'bioType': 'transcript', 'cdnaPos': '1173-1174', 'exons': '8/24', 'geneId':
                                '1788', 'hgnc': 'DNMT3A', 'consequence': ['synonymous'], 'hgvsc':
                                'NR_135490.2:n.1173_1174delinsAT'}]
    lof_gene_list = ['DNMT3A']
    gene_name_to_3prime_region = {'DNMT3A': ('chr2', 25455829, 25455879)}
    pvs_score, pvs_rationale = pathogenic_classifiers.get_pvs(variant, lof_gene_list, gene_name_to_3prime_region)
    assert pvs_score == 0
    assert pvs_rationale == [(0, "Out of 2 transcripts, variant consequence stop_gained in transcript NM_022552.5 was only seen once.")]


def test_get_ps1():
    """
    Test the strong pathogenic classifier, PS1: Same amino acid change as a previously established pathogenic variants
    regardless of nucleotide change
    """
    # Test data setup
    gene_mut_to_data = {
        ("SAMD11", "R793*"): ("761448939", "criteria provided, single submitter", "Pathogenic"),
    }
    gene_name = "SAMD11"
    aa_mut = "R793*"
    expected_score = 2
    expected_ps1 = ("PS1 (2): Amino acid change R793* in gene SAMD11 is associated with Clinvar Pathogenic variant "
                    "rs761448939 with review status criteria provided, single submitter")

    score, ps1 = pathogenic_classifiers.get_ps1(gene_name, aa_mut, gene_mut_to_data)
    assert (score, ps1) == (expected_score, expected_ps1)


def test_get_ps3_single_pubmed():
    """
    Test the get_ps3 function with a single pubmed ID supporting the variant.
    """
    # Initialize the variant data
    variant = extract_from_nirvana_json.VariantData("chr1", 1014142, "A", "T", "SNV")

    # Set up the test data
    gloc_to_pubmed_id_list = {
        ("chr1", 1014142): ["25307056"],
        ("chr1", 1014143): ["25307056"]
    }

    # Perform the test
    score, ps3 = pathogenic_classifiers.get_ps3(variant, gloc_to_pubmed_id_list)
    expected_score = 1
    expected_ps3 = "PS3 (1) - Variants pathogenicity is supported by pubmed study 25307056 as established by AVADA"
    assert (score, ps3) == (expected_score, expected_ps3)


def test_get_ps3_multiple_pubmeds():
    """
    Test the get_ps3 function with multiple pubmed IDs supporting the variant.
    """
    # Initialize the variant data
    variant = extract_from_nirvana_json.VariantData("chr1", 1014142, "A", "T", "SNV")

    # Set up the test data
    gloc_to_pubmed_id_list = {
        ("chr1", 1014142): ["25307056", "22859821"],
        ("chr1", 1014143): ["25307056"]
    }

    # Perform the test
    score, ps3 = pathogenic_classifiers.get_ps3(variant, gloc_to_pubmed_id_list)
    expected_score = 2
    expected_ps3 = "PS3 (2) - Variants pathogenicity is supported by multiple pubmed studies 22859821, 25307056 as established by AVADA"
    assert (score, ps3) == (expected_score, expected_ps3)

def test_get_ps4():
    """
    Test the get_ps4 function.
    """
    # Initialize the test data
    chrom_to_pos_to_gwas_data = {
        "chr1": {1014142:(6.0784, 0, 0, "30578418", "Pulse pressure", 7e-20, 'rs12345')},
    }

    # Initialize the variant data
    variant = extract_from_nirvana_json.VariantData(
        chromosome="chr1",
        position=1014142,
        ref_allele="A",
        alt_allele="T",
        variant_type="SNV"
    )
    variant.dbSnpIds = "rs537750491"

    # Perform the test
    score, ps4 = pathogenic_classifiers.get_ps4(variant, chrom_to_pos_to_gwas_data)
    expected_score = 3
    expected_ps4 = "PS4 (0): Genomic location chr1:1014142 is associated with trait Pulse pressure. It has OR 6.0784 " + \
            "with p-value 7e-20 and CI 0-0 as reported in PubMed 30578418"
    assert (score, ps4) == (expected_score, expected_ps4)

def test_get_pm1():
    """
    Test the get_pm1 function.
    """
    # Initialize the pathogenic domains data
    pathogenic_domains = [
        ("chr1", "17350491", "17354258", "P21912", 0.85091, 10),
    ]

    # Initialize the variant data
    variant = extract_from_nirvana_json.VariantData(
        chromosome="chr1",
        position=17350500,
        ref_allele="A",
        alt_allele="T",
        variant_type="SNV"
    )

    # Perform the test
    score, pm1 = pathogenic_classifiers.get_pm1(variant, pathogenic_domains)
    expected_score = 1
    expected_pm1 = "PM1 (1): Variant is in domain, P21912 where (0.85091 %) variants are associated with pathogenicity (var score 10)."
    assert (score, pm1) == (expected_score, expected_pm1)

def test_get_pm1_strong():
    """
    Test the get_pm1 function with very strong evidence
    """
    # Initialize the pathogenic domains data
    pathogenic_domains = [
        ("chr1", "17350491", "17354258", "P21912", 0.9090909091, 33),
    ]

    # Initialize the variant data
    variant = extract_from_nirvana_json.VariantData(
        chromosome="chr1",
        position=17350500,
        ref_allele="A",
        alt_allele="T",
        variant_type="SNV"
    )

    # Perform the test
    score, pm1 = pathogenic_classifiers.get_pm1(variant, pathogenic_domains)
    expected_score = 2
    expected_pm1 = "PM1 (2): Variant is in domain, P21912 where (0.9090909091 %) variants are associated with pathogenicity (var score 33)."
    assert (score, pm1) == (expected_score, expected_pm1)

def test_get_pm1_very_strong():
    """
    Test the get_pm1 function with very strong evidence
    """
    # Initialize the pathogenic domains data
    pathogenic_domains = [
        ("chr1", "17350491", "17354258", "P21912", 0.970909091, 75),
    ]

    # Initialize the variant data
    variant = extract_from_nirvana_json.VariantData(
        chromosome="chr1",
        position=17350500,
        ref_allele="A",
        alt_allele="T",
        variant_type="SNV"
    )

    # Perform the test
    score, pm1 = pathogenic_classifiers.get_pm1(variant, pathogenic_domains)
    expected_score = 3
    expected_pm1 = "PM1 (3): Variant is in domain, P21912 where (0.970909091 %) variants are associated with pathogenicity (var score 75)."
    assert (score, pm1) == (expected_score, expected_pm1)
