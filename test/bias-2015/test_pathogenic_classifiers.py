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

def test_get_pm2_absent():
    """
    Test the get_pm2 function when the variant is absent 
    """
    # Initialize the variant data
    variant = extract_from_nirvana_json.VariantData(
        chromosome="chr1",
        position=17350500,
        ref_allele="A",
        alt_allele="T",
        variant_type="SNV",
    )
    
    variant.gnomad = None
    variant.oneKg = None

    # Perform the test
    score, pm2 = pathogenic_classifiers.get_pm2(variant)
    expected_score = 1
    expected_pm2 = "PM2 (1): Variant is missing from GNOMAD and 1000 Genomes"
    assert (score, pm2) == (expected_score, expected_pm2)


def test_get_pm2_low_af():
    """
    Test the get_pm2 function when the variant allele frequency is low
    """
    # Initialize the variant data
    variant = extract_from_nirvana_json.VariantData(
        chromosome="chr1",
        position=17350500,
        ref_allele="A",
        alt_allele="T",
        variant_type="SNV",
    )

    # Set gnomad and oneKg attributes
    variant.gnomad = {"allAf": 0.005}
    variant.oneKg = {"allAf": 0.002}

    # Perform the test
    score, pm2 = pathogenic_classifiers.get_pm2(variant)
    expected_score = 1
    expected_pm2 = "PM2 (1): 1000 Genomes general population AF=0.50% is below the threshold of 1.0%"
    assert (score, pm2) == (expected_score, expected_pm2)

def test_get_aa_comparator_simple():
    """
    Test the get_aa_comparator function with a simple protein notation
    """
    p_notation = "Y524S"
    result = pathogenic_classifiers.get_aa_comparator(p_notation)
    expected = "Y524"
    assert result == expected

def test_get_aa_comparator_complex():
    """
    Test the get_aa_comparator function with a complex protein notation
    """
    p_notation = "V552fsS26*"
    result = pathogenic_classifiers.get_aa_comparator(p_notation)
    expected = "V552"
    assert result == expected

def test_get_pm4_non_repeat_region():
    """
    Test the get_pm4 function when the variant is in a non-repeat region.
    """
    # Initialize the chrom_to_repeat_regions data
    chrom_to_repeat_regions = {
        "chr1": [(100000, 200000, "repeat1"), (300000, 400000, "repeat2")]
    }

    # Initialize the variant data
    variant = extract_from_nirvana_json.VariantData(
        chromosome="chr1",
        position=250000,
        ref_allele="ATG",
        alt_allele="A",
        variant_type="del"
    )
    
    variant.consequence = "inframe_deletion"

    # Perform the test
    score, pm4 = pathogenic_classifiers.get_pm4(variant, chrom_to_repeat_regions)
    expected_score = 1
    expected_pm4 = "PM4 (1): inframe_deletion of ATG in a non-repeat region"
    assert (score, pm4) == (expected_score, expected_pm4)

def test_get_pm4_repeat_region(): #failing
    """ 
    Test the get_pm4 function when the variant is in a repeat region.
    """
    # Initialize the chrom_to_repeat_regions data
    chrom_to_repeat_regions = {
        "chr1": [(100000, 200000, "repeat1"), (300000, 400000, "repeat2")]
    }

    # Initialize the variant data
    variant = extract_from_nirvana_json.VariantData(
        chromosome="chr1",
        position=150000,
        ref_allele="A",
        alt_allele="ATG",
        variant_type="ins"
    )

    # Set the consequence attribute separately
    variant.consequence = "inframe_insertion"

    # Perform the test
    score, pm4 = pathogenic_classifiers.get_pm4(variant, chrom_to_repeat_regions)
    expected_score = 1
    expected_pm4 = "PM4 (1): inframe_insertion of ATG in a repeat region"
    assert (score, pm4) == (expected_score, expected_pm4)

def test_get_pm5_pathogenic(): #failing
    """
    Test the get_pm5 function when the variant matches a known pathogenic variant.
    """
    # Initialize the gene_aa_to_var_data and gene_mut_to_data
    gene_aa_to_var_data = {
        ("AGRN", "Y425"): ("Y425*", "1239736447", "criteria provided, single submitter", "Pathogenic")
    }
    gene_mut_to_data = {}

    # Initialize the variant data
    variant = extract_from_nirvana_json.VariantData(
        chromosome="chr17",
        position=41276045,
        ref_allele="T",
        alt_allele="C",
        variant_type="SNV",
    )
    geneName="AGRN",
    protein_variant="Y425*"

    # Perform the test
    score, pm5 = pathogenic_classifiers.get_pm5(variant, gene_aa_to_var_data, gene_mut_to_data)
    expected_score = 2
    expected_pm5 = "PM5 (2): Variant AA change Y425* is the same base AA as Y425* a Pathogenic variant with rs id 1239736447 and review status criteria provided, single submitter."
    assert (score, pm5) == (expected_score, expected_pm5)

def test_get_pm5_non_pathogenic(): 
    """
    Test the get_pm5 function when the variant does not match a known pathogenic variant.
    """
    # Initialize the gene_aa_to_var_data and gene_mut_to_data
    gene_aa_to_var_data = {}
    gene_mut_to_data = {}

    # Initialize the variant data
    variant = extract_from_nirvana_json.VariantData(
        chromosome="chr17",
        position=41276045,
        ref_allele="T",
        alt_allele="C",
        variant_type="SNV",
    )
    
    geneName="BRCA1",
    protein_variant="Y524S"

    # Perform the test
    score, pm5 = pathogenic_classifiers.get_pm5(variant, gene_aa_to_var_data, gene_mut_to_data)
    expected_score = 0
    expected_pm5 = ""
    assert (score, pm5) == (expected_score, expected_pm5)

def test_get_pp2(): #failing
    """
    Test the get_pp2 function with a known missense variant in a gene that has a low rate
    of benign missense variation.
    """
    # Initialize missense pathogenic genes
    missense_pathogenic_genes = {"AGRN"}

    # Initialize the variant data
    variant = extract_from_nirvana_json.VariantData(
        chromosome="chr17",
        position=41276045,
        ref_allele="T",
        alt_allele="C",
        variant_type="SNV",
    )
    geneName="AGRN"
    consequence="missense"

    # Perform the test
    score, pp2 = pathogenic_classifiers.get_pp2(variant, missense_pathogenic_genes)
    expected_score = 1
    expected_pp2 = "PP2 (1): Missense variant consequence missense in AGRN which was found to have a low rate of missense variation and where missense variants are a common mechanism of disease"
    assert (score, pp2) == (expected_score, expected_pp2)

def test_get_pp3(): #failing
    """
    Test the get_pp3 function with different levels of computational evidence.
    """
    # Test with various values for GERP, DANN, and REVEL
    variant_1 = extract_from_nirvana_json.VariantData(
        chromosome="chr17",
        position=41276045,
        ref_allele="T",
        alt_allele="C",
        variant_type="SNV",
    )
    revel=0.6
    dann=0.8
    gerp=2.5

    # Perform the test
    score, pp3 = pathogenic_classifiers.get_pp3(variant_1)
    expected_score = 3
    expected_pp3 = "PP3 (3): Computational evidence support a deleteterious effect; gerp 2.5 | dann 0.8 | revel 0.6"
    assert (score, pp3) == (expected_score, expected_pp3)

    # Test with only one positive computational evidence
    variant_2 = extract_from_nirvana_json.VariantData(
        chromosome="chr17",
        position=41276045,
        ref_allele="T",
        alt_allele="C",
        variant_type="SNV",
        revel=0.4,
        dann=0.8,
        gerp=1.5
    )

    # Perform the test
    score, pp3 = pathogenic_classifiers.get_pp3(variant_2)
    expected_score = 1
    expected_pp3 = "PP3 (1): Computational evidence support a deleteterious effect; gerp 1.5 | dann 0.8 | revel 0.4"
    assert (score, pp3) == (expected_score, expected_pp3)

    # Test with no positive computational evidence
    variant_3 = extract_from_nirvana_json.VariantData(
        chromosome="chr17",
        position=41276045,
        ref_allele="T",
        alt_allele="C",
        variant_type="SNV",
    )
    revel=0.3
    dann=0.6
    gerp=1.0

    # Perform the test
    score, pp3 = pathogenic_classifiers.get_pp3(variant_3)
    expected_score = 0
    expected_pp3 = ""
    assert (score, pp3) == (expected_score, expected_pp3)

def test_get_pp5(): #failing
    """
    Test the get_pp5 function with a variant reported as pathogenic in ClinVar.
    """
    # Initialize ClinVar review statuses and significance
    clinvar_review_status_to_level = {
        "criteria provided, single submitter": 1,
        "criteria provided, multiple submitters, no conflicts": 2
    }

    # Initialize the variant data
    variant = extract_from_nirvana_json.VariantData(
        chromosome="chr17",
        position=41276045,
        ref_allele="T",
        alt_allele="C",
        variant_type="SNV",
    )

    clinvar_review_status="criteria provided, single submitter",
    clinvar_significance="pathogenic"

    # Perform the test
    score, pp5 = pathogenic_classifiers.get_pp5(variant)
    expected_score = 1
    expected_pp5 = "PP5 (1): Variant was found reported as pathogenic in ClinVar with review status of criteria provided, single submitter."
    assert (score, pp5) == (expected_score, expected_pp5)

    # Test with a variant not reported as pathogenic
    variant_2 = extract_from_nirvana_json.VariantData(
        chromosome="chr17",
        position=41276045,
        ref_allele="T",
        alt_allele="C",
        variant_type="SNV",
        clinvar_review_status="criteria provided, single submitter",
        clinvar_significance="benign"
    )

    # Perform the test
    score, pp5 = pathogenic_classifiers.get_pp5(variant_2)
    expected_score = 0
    expected_pp5 = ""
    assert (score, pp5) == (expected_score, expected_pp5)

