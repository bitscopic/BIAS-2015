"""
Load in the data files that will be used in the Bitscopic Interpreting ACMG Standards (BIAS) variant interpretation
process

The generation or aquisition of all input files i is provided in the doc directory under 'variant_interpretation.txt'
"""
import re

# PVS1
def get_lof_gene_to_pli(lof_gene_fp):
    """
    For PVS1, defining genes where loss of function is a mechanism of pathogenicity

    Loss of function genes (LOF) were defined as genes with a probability of Loss-of-function Intolerance (pLI) value
    greater than .9%.

    The pLI calculations and values are defined in the following GNOMAD/EXAC related publication.

    Karczewski, K. J., et al. (2020). The mutational constraint spectrum quantified from variation in 141,456 humans.
        Nature, 581(7809), 434–443. doi:10.1038/s41586-020-2308-7
    """
    lof_gene_to_pli = {}
    
    with open(lof_gene_fp) as in_file:
        for line in in_file:
            split_line = line.strip().split()
            lof_gene_to_pli[split_line[0]] = split_line[1]
    if not lof_gene_to_pli:
        print(f"File {lof_gene_fp} does not have any valid entries!")

    return lof_gene_to_pli


def calculate_3prime_region(chrom, strand, exon_starts, exon_ends, exon_frames):
    """
    Identify the final coding exon, then identify the last 50 bases. Return the
    region which spans these 50 bases as a tuple (chrom, start, stop).
    """
    # Only consider exons that are considered in-frame (exclude out of frame exons)
    coding_exon_starts = []
    coding_exon_ends = []
    exon_count = len(exon_starts)
    for x in range(0, exon_count):
        if exon_frames[x] == '-1':
            continue
        coding_exon_starts.append(exon_starts[x])
        coding_exon_ends.append(exon_ends[x])

    if not coding_exon_ends: return () 
    # For the positive strand, use the last exon
    if strand == "+":
        last_exon_end = coding_exon_ends[-1]
        start_3prime = last_exon_end - 50
        end_3prime = last_exon_end
    else: # For the negative strand, use the first exon
        first_exon_start = coding_exon_starts[0]
        end_3prime = first_exon_start + 50
        start_3prime = first_exon_start

    return chrom, start_3prime, end_3prime

# PVS1 Caveate 2
def get_gene_name_to_3prime_region(ncbi_ref_seq_hgmd_fp):
    """
    Load RefSeq hg19 data into a dictionary that maps gene name to the final 50 bases
    of coding region

    In general, nonsense-mediated decay NMD is not predicted to occur if the premature termination codon occurs in the 
    3′ most exon or within the 3′ most 50 nucleotides of the penultimate exon (Chang, Imam, & Wilkinson, 2007; Lewis,
    Green, & Brenner, 2003). When NMD is not predicted to occur, it is important to determine if the truncated or
    altered region is critical to protein function, often indicated by experimental or clinical evidence—such as
    pathogenic variants downstream of the new stop codon—supporting the biological relevance of the C-terminal region.

    Chang, Y. F., Imam, J. S., & Wilkinson, M. F. (2007). The nonsense-mediated decay RNA surveillance pathway. Annual 
        Review of Biochemistry, 76, 51–74. DOI: 10.1146/annurev.biochem.76.050106.093909
    
    Lewis, B. P., Green, R. E., & Brenner, S. E. (2003). Evidence for the widespread coupling of alternative splicing
        and nonsense-mediated mRNA decay in humans. Proceedings of the National Academy of Sciences, 100(1), 189-192.
        DOI: 10.1073/pnas.0136770100
    """
    gene_name_to_3prime_region = {}
    with open(ncbi_ref_seq_hgmd_fp, 'r') as in_file:
        for line in in_file:
            fields = line.strip().split("\t")
            chrom = fields[2]
            strand = fields[3]
            exon_starts = list(map(int, fields[9].strip()[:-1].split(",")))
            exon_ends = list(map(int, fields[10].strip()[:-1].split(",")))
            exon_frames = fields[-1].split(",")
            gene_name = fields[12]  # Extracting the gene name
            region = calculate_3prime_region(chrom, strand, exon_starts, exon_ends, exon_frames)
            if region:
                gene_name_to_3prime_region[gene_name] = region
    if not gene_name_to_3prime_region:
        print(f"File {ncbi_ref_seq_hgmd_fp} does not have any valid entries!")
    return gene_name_to_3prime_region


def get_aa_comparator(p_notation):
    """
    AA look like 'Y524S', 'E525K', 'V552fsS26*', we want to compare the first AA and position
    to the current variant
    """
    first_half = ""
    pos_read = True
    nums = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']
    for char in p_notation:
        if char in nums:
            pos_read = False
        if not pos_read and char not in nums:
            break
        first_half += char
    return first_half


# PS1
def get_gene_mut_to_data(clinvar_pathogenic_aa_fp):
    """
    Load in amino acid changes that are associated with known pathogenic and likely pathogenic variants.

    The method of using a variant database to define amino acid changes for the PS1 and PM5 classifiers is explored
    in this publication.

    Helbig, I., Nothnagel, M., May, P., Brünger, T., Ivaniuk, A., Pérez-Palma, E., Montanucci, L., Cohen, S., Smith, L.,
        Parthasarathy, S., & Lal, D. (2023). Pathogenic paralogous variants can be used to apply the ACMG PS1 and PM5
        variant interpretation criteria. medRxiv. https://doi.org/10.1101/2023.08.22.23294353.
    """
    gene_mut_to_data = {}
    gene_aa_to_var_data = {}
    with open(clinvar_pathogenic_aa_fp, 'r') as in_file:
        for line in in_file:
            split_line = line.strip().split("\t")
            if len(split_line) < 5: continue
            gene_name, p_notation, rs_id, criteria, signif = split_line
            aa_change = get_aa_comparator(p_notation)
            gene_mut = (gene_name, p_notation) # The full gene and mutation; ex ('SRY', 'V60L')
            gene_aa = (gene_name, aa_change) # The full gene and partial mutation; ex ('SRY', 'V60') with the final aa excluded
            gene_aa_to_var_data[gene_aa] = (p_notation, rs_id, criteria, signif)
            gene_mut_to_data[gene_mut] = (rs_id, criteria, signif) 
    if not gene_mut_to_data:
        print(f"File {clinvar_pathogenic_aa_fp} does not have any valid entries!")

    return gene_mut_to_data, gene_aa_to_var_data

# PS3
def get_gloc_to_pubmed_id_list(literature_supported_variants_fp):
    """
    AVADA uses AI to find pathogenic varaints from publications and link them to genomic coordinates.
    
    Johannes Birgmeier, Cole A. Deisseroth, Laura E. Hayward, Luisa M. T. Galhardo, Andrew P. Tierno, Karthik A.
        Jagadeesh, Peter D. Stenson, David N. Cooper, Jonathan A. Bernstein, Maximilian Haeussler, and Gill Bejerano.
        AVADA: Towards Automated Pathogenic Variant Evidence Retrieval Directly from the Full Text Literature.
        Genetics in Medicine. 2019. PMID: 31467448
    """
    # This will take lot of RAM but be very fast to lookup positions
    gloc_to_pubmed_id_list = {}

    with open(literature_supported_variants_fp, 'r') as in_file:
        for line in in_file:
            split_line = line.strip().split("\t")
            if len(split_line) < 4: continue
            chrom = split_line[0]
            start = split_line[1]
            stop = split_line[2]
            pubmed_id = split_line[3]
            pos = (chrom, start, stop) 
            if gloc_to_pubmed_id_list.get(pos):
                if pubmed_id not in gloc_to_pubmed_id_list[pos]: 
                    gloc_to_pubmed_id_list[pos].append(pubmed_id)
            else:
                gloc_to_pubmed_id_list[pos] = [pubmed_id]
    if not gloc_to_pubmed_id_list:
        print(f"File {literature_supported_variants_fp} does not have any valid entries!")
    return gloc_to_pubmed_id_list


def extract_floats(input_string):
    """
    The floats can be negative so we use a re 
    """
    # Define the regular expression pattern to match floats
    pattern = r'(-?\d+\.\d+)'
    
    # Find all matches in the input string
    matches = re.findall(pattern, input_string)
    
    first_float, last_float = 0, 0
    # Extract the first and last float
    if len(matches) >= 2:
        first_float = float(matches[0])
        last_float = float(matches[-1])

    return first_float, last_float


# PS4
def get_dbsnpids_to_or(gwas_dbsnp_fp):
    """
    The NHGRI-EBI GWAS Catalog is a curated collection of all published genome-wide association studies, produced by a 
    collaboration between EMBL-EBI and the National Human Genome Research Institute (NHGRI). It includes comprehensive 
    data on SNP-trait associations with genome-wide significance.

    Hindorff LA, Sethupathy P, Junkins HA, Ramos EM, Mehta JP, Collins FS, Manolio TA. Potential etiologic and
        functional implications of genome-wide association loci for human diseases and traits. Proc Natl Acad Sci U S A.
        2009 Jun 9;106(23):9362-7. PMID: 19474294; PMC: PMC2687147

    The NHGRI-EBI GWAS Catalog: a curated resource of SNP-trait associations. Buniello, A., MacArthur, J.A.L., 
        Cerezo, M., et al. Nucleic Acids Research, 2019, 47(D1): D1005-D1012. 
        https://academic.oup.com/nar/article/47/D1/D1005/5184712

    """
    # Map chrom to pos
    chrom_to_pos_to_gwas_data = {}
    with open(gwas_dbsnp_fp, 'r') as in_file:
        for line in in_file:
            split_line = line.strip().split("\t")
            if len(split_line) < 23: continue
            chrom, start, stop, rs_id, pubmed_id = split_line[1:6]
            trait = split_line[9]
            p_value = split_line[17]
            or_value = split_line[19]
            if not or_value: continue
            or_value = float(or_value)
            p_value = float(p_value)
            if or_value < 1.5: continue
            if p_value > .05: continue 
            ci_text = split_line[20]
            ci_lower, ci_upper = extract_floats(ci_text.split(" ")[0])
            data = (or_value, float(ci_lower), float(ci_upper), pubmed_id, trait, p_value, rs_id)
            for pos in range(int(start), int(stop)):
                if chrom_to_pos_to_gwas_data.get(chrom):
                    chrom_to_pos_to_gwas_data[chrom][pos] = data
                else:
                    chrom_to_pos_to_gwas_data[chrom] = {pos: data}
    if not chrom_to_pos_to_gwas_data:
        print(f"File {gwas_dbsnp_fp} does not have any valid entries!")
    return chrom_to_pos_to_gwas_data

# PM1
def get_pathogenic_domains(pathogenic_domains_fp):
    """
    Load in the pathogenic domains that are identified using Clinvar (see variant documentation) and the
    UCSC Genome Browser hg19 UniProt track

    UniProt Consortium. Reorganizing the protein space at the Universal Protein Resource (UniProt).
        Nucleic Acids Res. 2012 Jan;40(Database issue):D71-5. PMID: 22102590; PMC: PMC3245120

    Yip YL, Scheib H, Diemand AV, Gattiker A, Famiglietti LM, Gasteiger E, Bairoch A. The Swiss-Prot variant page and
        the ModSNP database: a resource for sequence and structure information on human protein variants.
        Hum Mutat. 2004 May;23(5):464-70. PMID: 15108278
    """
    pathogenic_domains = []
    with open(pathogenic_domains_fp, 'r') as in_file:
        for line in in_file:
            split_line = line.strip().split("\t")
            if len(split_line) < 6: continue
            chrom = split_line[0]
            start = split_line[1]
            stop = split_line[2]
            uniprot_acc = split_line[3]
            path_per = float(split_line[4])
            score = int(split_line[5])
            pathogenic_domains.append((chrom, start, stop, uniprot_acc, path_per, score))
    if not pathogenic_domains:
        print(f"File {pathogenic_domains_fp} does not have any valid entries!")
    return pathogenic_domains


# PM4
def get_chrom_to_repeat_regions(coding_repeat_region_fp):
    """
    UCSC repeat masker regions intersected with coding regions returned as a
    dict mapping chrom as a string to a list of repeat regions
    

    Repeat Masker:
    Smit AFA, Hubley R, Green P. RepeatMasker Open-3.0. http://www.repeatmasker.org. 1996-2010.

    Repbase Update is described in:

    Jurka J. Repbase Update: a database and an electronic journal of repetitive elements. Trends Genet.
        2000 Sep;16(9):418-420. PMID: 10973072

    For a discussion of repeats in mammalian genomes, see:

    Smit AF. Interspersed repeats and other mementos of transposable elements in mammalian genomes. Curr Opin Genet
        Dev. 1999 Dec;9(6):657-63. PMID: 10607616

    Smit AF. The origin of interspersed repeats in the human genome. Curr Opin Genet Dev. 1996 Dec;6(6):743-8.
        PMID: 8994846

    Consensus Coding:
    Hubbard T, Barker D, Birney E, Cameron G, Chen Y, Clark L, Cox T, Cuff J, Curwen V, Down T et al. The Ensembl
        genome database project. Nucleic Acids Res. 2002 Jan 1;30(1):38-41. PMID: 11752248; PMC: PMC99161

    Pruitt KD, Harrow J, Harte RA, Wallin C, Diekhans M, Maglott DR, Searle S, Farrell CM, Loveland JE, Ruef BJ et al. 
        The consensus coding sequence (CCDS) project: Identifying a common protein-coding gene set for the human and
        mouse genomes. Genome Res. 2009 Jul;19(7):1316-23. PMID: 19498102; PMC: PMC2704439

    Pruitt KD, Tatusova T, Maglott DR. NCBI Reference Sequence (RefSeq): a curated non-redundant sequence database of
        genomes, transcripts and proteins. Nucleic Acids Res. 2005 Jan 1;33(Database issue):D501-4. PMID: 15608248;
        PMC: PMC539979

    """
    chrom_to_repeat_regions = {}
    with open(coding_repeat_region_fp, 'r') as in_file:
        for line in in_file:
            split_line = line.strip().split("\t")
            if len(split_line) < 4: continue
            chrom, start, stop, strand = split_line
            region = (int(start), int(stop), strand)
            chrom_to_repeat_regions.setdefault(chrom, []).append(region)
    if not chrom_to_repeat_regions:
        print(f"File {coding_repeat_region_fp} does not have any valid entries!")
    return chrom_to_repeat_regions

# PP2
def get_missense_pathogenic_genes(missense_pathogenic_genes_fp):
    """
    Load in genes where missense is a common mechanism of pathogenicity
    """
    missense_pathogenic_genes = []

    with open(missense_pathogenic_genes_fp, 'r') as in_file:
        for line in in_file:
            split_line = line.strip().split("\t")
            if len(split_line) < 3: continue
            # NOTE: Pathogenic and benign percentage should be passed through and included in the classifier justification
            gene, _, _ = split_line
            missense_pathogenic_genes.append(gene)
    if not missense_pathogenic_genes:
        print(f"File {missense_pathogenic_genes_fp} does not have any valid entries!")
    return missense_pathogenic_genes

# BP1
def get_truncating_genes(truncating_genes_fp):
    """
    Genes identifies from ClinVar where over 80% of pathogenic variants are truncating
    variants
    """
    truncating_gene_list = []
    with open(truncating_genes_fp, 'r') as in_file:
        for line in in_file:
            split_line = line.strip().split("\t")
            if len(split_line) < 2: continue
            # NOTE: Pathogenic percentage should be passed through and included in the classifier justification
            gene, _ = split_line
            truncating_gene_list.append(gene)
    if not truncating_gene_list:
        print(f"File {truncating_genes_fp} does not have any valid entries!")
    return truncating_gene_list


def get_benign_domains(benign_domains_fp):
    """
    Load in the benign domains that are identified using Clinvar (see variant documentation)
    """
    benign_domains = []
    with open(benign_domains_fp, 'r') as in_file:
        for line in in_file:
            split_line = line.strip().split("\t")
            if len(split_line) < 4: continue
            chrom = split_line[0]
            start = split_line[1]
            stop = split_line[2]
            uniprot_acc = split_line[3]
            benign_domains.append((chrom, start, stop, uniprot_acc))
    if not benign_domains:
        print(f"File {benign_domains_fp} does not have any valid entries!")
    return benign_domains


def get_name_to_dataset(file_paths):
    """
    Gather the datasets and format them for easy python consumption
    """
    name_to_dataset = {}
    name_to_dataset['lof_gene_to_pli'] = get_lof_gene_to_pli(file_paths['lof_gene_fp']) # PVS1 
    name_to_dataset['gene_name_to_3prime_region'] = get_gene_name_to_3prime_region(file_paths['ncbi_ref_seq_hgmd_fp']) # PVS1 caveat
    name_to_dataset['gene_mut_to_data'], name_to_dataset['gene_aa_to_var_data'] = get_gene_mut_to_data(file_paths['clinvar_pathogenic_aa_fp']) # PS1/PM5
    name_to_dataset['chrom_to_pos_to_gwas_data'] = get_dbsnpids_to_or(file_paths['gwas_dbsnp_fp']) # PS4
    name_to_dataset['chrom_to_repeat_regions'] = get_chrom_to_repeat_regions(file_paths['coding_repeat_region_fp']) # PM4 & BP1
    name_to_dataset['missense_pathogenic_genes'] = get_missense_pathogenic_genes(file_paths['missense_pathogenic_genes_fp']) # PP2
    name_to_dataset['truncating_genes'] = get_truncating_genes(file_paths['truncating_genes_fp']) # BP1
    name_to_dataset['pathogenic_domains'] = get_pathogenic_domains(file_paths['pathogenic_domains_fp']) # PM1
    #name_to_dataset['benign_domains'] = get_benign_domains(file_paths['benign_domains_fp']) # TODO: BP1
    name_to_dataset['gloc_to_pubmed_id_list'] = get_gloc_to_pubmed_id_list(file_paths['literature_supported_variants_fp']) # PS3
    return name_to_dataset
