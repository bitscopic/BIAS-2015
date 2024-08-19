"""
Unit testing for the Bitscopic Interpreting ACMG Standards (BIAS) dataset loader
"""
import os
import sys
import tempfile

# NOTE: This should be handled more elegantly for public facing code
# Some hacky path adjustments so users dont have to fiddle with env variables
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from src.bias_2015 import variant_interpretation_dataset_loader as vidl


# PVS
def test_get_lof_gene_to_pli():
    """
    opens the loss of function file and returns its contents as a list
    """
    test_data = "gene1\t1\ngene2\t.95\ngene3\t.9"

    # Create a temporary file
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt') as tmp_file:
        lof_gene_fp = tmp_file.name
        tmp_file.write(test_data)

    try:
        # Perform the test
        lof_gene_to_pli = vidl.get_lof_gene_to_pli(lof_gene_fp)
        assert lof_gene_to_pli == {"gene1": "1", "gene2":".95", "gene3":".9"}
    finally:
        # Clean up the temporary file
        os.remove(lof_gene_fp)


def test_calculate_3prime_region():
    """
    Test a default calculation of the 3' region with a positive strand
    """
    chrom = "chr1"
    strand = "+"
    exon_starts = [25455829, 25458575, 25459804, 25461998]
    exon_ends = [25457289, 25458694, 25459874, 25462084]
    exon_frames = ['0', '0', '0', '0']
    region = vidl.calculate_3prime_region(chrom, strand, exon_starts, exon_ends, exon_frames)
    assert region == ('chr1', 25462034, 25462084)


def test_calculate_3prime_region_neg():
    """
    Test a default calculation of the 3' region with a negative strand
    """
    chrom = "chr1"
    strand = "-"
    exon_starts = [25455829, 25458575, 25459804, 25461998]
    exon_ends = [25457289, 25458694, 25459874, 25462084]
    exon_frames = ['0', '0', '0', '0']
    region = vidl.calculate_3prime_region(chrom, strand, exon_starts, exon_ends, exon_frames)
    assert region == ('chr1', 25455829, 25455879)


def test_calculate_3prime_region_oof_exons():
    """
    Test a calculation of the 3' region with several out of frame (oof) exons
    that should be excluded
    """
    chrom = "chr1"
    strand = "+"
    exon_starts = [25455829, 25458575, 25459804, 25461998]
    exon_ends = [25457289, 25458694, 25459874, 25462084]
    exon_frames = ['0', '-1', '0', '-1']
    region = vidl.calculate_3prime_region(chrom, strand, exon_starts, exon_ends, exon_frames)
    assert region == ('chr1', 25459824, 25459874)


def test_get_gene_name_to_3prime_region():
    """
    Test the get_gene_name_to_3prime_region function.
    """
    test_data = "591\tNM_152486.4\tchr1\t+\t861110\t879954\t861321\t879533\t14\t861110,861301,865534,866418,871151,874419,874654,876523,877515,877789,877938," + \
                "878632,879077,879287,\t861180,861393,865716,866469,871276,874509,874840,876686,877631,877868,878438,878757,879188,879954,\t0\tSAMD11\tcmpl\t" + \
                "cmpl\t-1,0,0,2,2,1,1,1,2,1,2,1,0,0,\n" + \
                "591\tNM_198317.3\tchr1\t+\t895963\t901099\t896073\t900571\t12\t895963,896672,897008,897205,897734,898083,898488,898716,899299,899486,899728," + \
                "900342,\t896180,896932,897130,897427,897851,898297,898633,898884,899388,899560,899910,901099,\t0\tKLHL17\tcmpl\tcmpl\t0,2,1,0,0,0,1,2,2,1,0,2,"

    # Create a temporary file
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt') as tmp_file:
        ncbi_ref_seq_hgmd_fp = tmp_file.name
        tmp_file.write(test_data)

    try:
        # Perform the test
        gene_name_to_3prime_region = vidl.get_gene_name_to_3prime_region(ncbi_ref_seq_hgmd_fp)
        expected_output = {
            "SAMD11": ("chr1", 879904, 879954),
            "KLHL17": ("chr1", 901049, 901099)
        }
        assert gene_name_to_3prime_region == expected_output
    finally:
        # Clean up the temporary file
        os.remove(ncbi_ref_seq_hgmd_fp)

# PS1
def disabled_test_get_gene_mut_to_data():
    """
    NOTE: The implementation of this has changed and this test is now outdated. It needs to be re written. 
    Test the get_gene_mut_to_data function.
    """
    test_data = "SAMD11\tR793*\t761448939\tcriteria provided, single submitter\tPathogenic\n" \
                "ISG15\tE127*\t672601312\tcriteria provided, single submitter\tLikely_pathogenic\n" \
                "AGRN\tG76S\t756623659\tcriteria provided, single submitter\tLikely_pathogenic"

    # Create a temporary file
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt') as tmp_file:
        clinvar_pathogenic_aa_fp = tmp_file.name
        tmp_file.write(test_data)

    try:
        # Perform the test
        gene_mut_to_data = vidl.get_gene_mut_to_data(clinvar_pathogenic_aa_fp)
        expected_output = {
            ("SAMD11", "R793*"): ("761448939", "criteria provided, single submitter", "Pathogenic"),
            ("ISG15", "E127*"): ("672601312", "criteria provided, single submitter", "Likely_pathogenic"),
            ("AGRN", "G76S"): ("756623659", "criteria provided, single submitter", "Likely_pathogenic")
        }
        assert gene_mut_to_data == expected_output
    finally:
        # Clean up the temporary file
        os.remove(clinvar_pathogenic_aa_fp)

# PS3
def test_get_gloc_to_pubmed_id_list():
    """
    Test the get_gloc_to_pubmed_id_list function.
    """
    test_data = "chr1\t1014142\t1014143\t25307056\n" + \
                "chr1\t1014142\t1014143\t25307056\n" + \
                "chr1\t1014313\t1014314\t22859821"

    # Create a temporary file
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt') as tmp_file:
        literature_supported_variants_fp = tmp_file.name
        tmp_file.write(test_data)

    try:
        # Perform the test
        gloc_to_pubmed_id_list = vidl.get_gloc_to_pubmed_id_list(literature_supported_variants_fp)
        expected_output = {
            ("chr1", 1014142): ["25307056"],
            ("chr1", 1014313): ["22859821"]
        }
        assert gloc_to_pubmed_id_list == expected_output
    finally:
        # Clean up the temporary file
        os.remove(literature_supported_variants_fp)

# PS4
def disabled_test_get_dbsnpids_to_or():
    """
    NOTE: The implementation of this was updated and this test is outdated. It will need to be re written
    Test the get_dbsnpids_to_or function.
    """
    test_data = "30578418\tPulse pressure\trs537750491\t7e-20\t6.0784\t4.77\t-7.38\n" + \
                "30578418\tPulse pressure\trs530130707\t4e-19\t6.796\t5.3\t-8.29\n" + \
                "31490055\tCocaine use disorder x non-traditional parental care interaction\tchr1:15511771\t9e-10\t6.124\t0\t0"

    # Create a temporary file
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt') as tmp_file:
        gwas_dbsnp_fp = tmp_file.name
        tmp_file.write(test_data)

    try:
        # Perform the test
        dbsnpid_to_data = vidl.get_dbsnpids_to_or(gwas_dbsnp_fp)
        expected_output = {
            "rs537750491": (6.0784, 4.77, -7.38, "30578418", "Pulse pressure", 7e-20),
            "rs530130707": (6.796, 5.3, -8.29, "30578418", "Pulse pressure", 4e-19),
            "chr1:15511771": (6.124, 0.0, 0.0, "31490055", "Cocaine use disorder x non-traditional parental care interaction", 9e-10)
        }
        assert dbsnpid_to_data == expected_output
    finally:
        # Clean up the temporary file
        os.remove(gwas_dbsnp_fp)

# PM1
def test_get_pathogenic_domains():
    """
    Test the get_pathogenic_domains function.
    """
    test_data = "chr1\t11076971\t11082252\tQ13148\t1\t8\n" + \
                "chr1\t12052712\t12061881\tO95140\t1\t11\n" + \
                "chr1\t12067112\t12067115\tO95140\t1\t11\n" + \
                "chr1\t17350491\t17354258\tP21912\t0.9090909091\t10\n" + \
                "chr1\t17355118\t17371338\tP21912\t0.9090909091\t10"

    # Create a temporary file
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt') as tmp_file:
        pathogenic_domains_fp = tmp_file.name
        tmp_file.write(test_data)

    try:
        # Perform the test
        pathogenic_domains = vidl.get_pathogenic_domains(pathogenic_domains_fp)
        expected_output = [
            ("chr1", "11076971", "11082252", "Q13148", 1.0, 8),
            ("chr1", "12052712", "12061881", "O95140", 1.0, 11),
            ("chr1", "12067112", "12067115", "O95140", 1.0, 11),
            ("chr1", "17350491", "17354258", "P21912", 0.9090909091, 10),
            ("chr1", "17355118", "17371338", "P21912", 0.9090909091, 10)
        ]
        assert pathogenic_domains == expected_output
    finally:
        # Clean up the temporary file
        os.remove(pathogenic_domains_fp)

#PM4
def test_get_chrom_to_repeat_regions():
    """
    Test the get_chrom_to_repeat_regions function with three lines of data
    """
    # Data from hg19_coding_repeat_regions.tsv
    test_data = "chr1\t16777160\t16777470\t+\n" + \
                "chr1\t25165800\t25166089\t-\n" + \
                "chr1\t33553606\t33554646\t+\n"

    # Create a temporary file
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.tsv') as tmp_file:
        coding_repeat_region_fp = tmp_file.name
        tmp_file.write(test_data)

    try:
        # Perform the test
        chrom_to_repeat_regions = vidl.get_chrom_to_repeat_regions(coding_repeat_region_fp)
        expected_output = {
            "chr1": [
                (16777160, 16777470, "+"),
                (25165800, 25166089, "-"),
                (33553606, 33554646, "+")
            ]
        }
        assert chrom_to_repeat_regions == expected_output
    finally:
        # Clean up the temporary file
        os.remove(coding_repeat_region_fp)

#PP2
def test_get_missense_pathogenic_genes():
    """
    Test the get_missense_pathogenic_genes function
    """
    # Data from hg19_missense_pathogenic_genes.tsv
    test_data = "AAAS\t1.0\t0.0\n" + \
                "ABCC9\t1.0\t0.006745362563237774\n" + \
                "ABCD1\t1.0\t0.0\n"

    # Expected output
    expected_genes = ["AAAS", "ABCC9", "ABCD1"]

    # Create a temporary file
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt') as tmp_file:
        missense_pathogenic_genes_fp = tmp_file.name
        tmp_file.write(test_data)

    try:
        # Perform the test
        missense_pathogenic_genes = vidl.get_missense_pathogenic_genes(missense_pathogenic_genes_fp)
        assert missense_pathogenic_genes == expected_genes
    finally:
        # Clean up the temporary file
        os.remove(missense_pathogenic_genes_fp)

#BP1
def test_get_truncating_genes():
    """
    Test the get_truncating_genes function
    """
    # Data from hg19_truncating_genes.tsv
    test_data = "ABCA2\t1.0\n" + \
                "ABCC2\t1.0\n" + \
                "ABHD12\t1.0\n" + \
                "ADAM17\t0.8333333333333334\n" + \
                "ADAMTS13\t1.0\n"

    # Expected output
    expected_genes = ["ABCA2", "ABCC2", "ABHD12", "ADAM17", "ADAMTS13"]

    # Create a temporary file
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt') as tmp_file:
        truncating_genes_fp = tmp_file.name
        tmp_file.write(test_data)

    try:
        # Perform the test
        truncating_genes = vidl.get_truncating_genes(truncating_genes_fp)
        assert truncating_genes == expected_genes
    finally:
        # Clean up the temporary file
        os.remove(truncating_genes_fp)
