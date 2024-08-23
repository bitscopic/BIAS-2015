"""
BIAS-2015 implementation of the ACMG 2015 germline pathogenic classifiers
"""

# Map ClinVar reviewStatus to an integer value
# Weight the score based on the ClinVar review status. This is an application of the established fact that
# expertly reviewed entries should be considered more strongly than entries with limited supported evidence.
#T his concept and the concept of weighted scoring is explored in the following publications;

# Landrum, M. J., Lee, J. M., Benson, M., Brown, G. R., Chao, C., Chitipiralla, S., ... & Maglott, D. R. (2018). 
#    ClinVar: improving access to variant interpretations and supporting evidence. Nucleic Acids Research,
#    46(D1), D1062-D1067. https://doi.org/10.1093/nar/gkx1153.

# Li, Q., Wang, K., & Zhang, X. (2020). A scoring system integrating genetic and clinical data for the clinical
#    interpretation of variants in hereditary cancer genes. NPJ Genomic Medicine, 5(1), 1-12.
#    https://www.nature.com/articles/s41525-020-00146-4.
clinvar_review_status_to_level = {
    "practice guideline": 5,
    "reviewed by expert panel": 5,
    "criteria provided, multiple submitters, no conflicts": 3,
    "criteria provided, single submitter": 2,
    "criteria provided, conflicting interpretations": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0
}


def get_pvs(variant, lof_gene_to_pli, gene_name_to_3prime_region):
    """
    Bitscopic variant object
    List of gene names where loss of function is a known mechanism of disease
    Dict mapping gene names to the 3' region of their final exon 

    Very strong evidence of pathogenicity
        PVS1    Null variants (nonsense, frameshift, canonical +/−1 or 2 splice sites, initiation
                codon, single or multi-exon deletion) in a gene where loss of function (LOF)
                is a known mechanism of disease

        Caveats:
            Beware of genes where LOF is not a known disease mechanism (e.g. GFAP, MYH7)
            Use caution interpreting LOF variants at the extreme 3’ end of a gene
            Use caution with splice variants that are predicted to lead to exon skipping but leave the remainder of the protein intact
            Use caution in the presence of multiple transcripts
    """
    # Only consider variants with a null consequence 
    null_variant_consequence_list = ['frameshift_variant',
                                     'stop_gained',
                                     'start_lost',
                                     'splice_donor_variant',
                                     'splice_acceptor_variant'
                                     ]
    if variant.consequence not in null_variant_consequence_list:
        return 0, [(0, f"{variant.consequence} is not a null consequence")]

    # Caveat 1: Genes where loss of function is a not known mechanism of disease
    if variant.geneName not in lof_gene_to_pli:
        return 0, [(0, f"{variant.geneName} is not a recognized gene where lof is a known mechanism of disease")]

    # Caveat 2: LOF variant is at the extreme 3' end of a gene (final 50 of last coding exon)
    null_region = ('chr0', -10, -1)
    prime_region = gene_name_to_3prime_region.get(variant.geneName, null_region)
    _, prime_start, prime_end = prime_region
    if prime_start <= int(variant.position) <= prime_end:
        return 0, [(0, f"Position {variant.position} is in the extreme three prime region ({prime_start}-{prime_end}) of the last coding exon in gene {variant.geneName}")]

    # Caveat 3: Check splice variants that lead to exon skipping but leave the remainder of the protein intact
    # NOTE: This requires a splice site predictor to implement, see this article comparing the top algorithms
    # https://pubmed.ncbi.nlm.nih.gov/33942434/
    # One must also be cautious in assuming that a null variant will lead to disease if found in an exon where no
    # other pathogenic variants have been described given the possibility that the exon may be alternatively spliced.
    # This is particularly true if the predicted truncating variant is identified as an incidental finding (unrelated
    # to the indication for testing) given the low prior likelihood of finding a pathogenic variant in that setting.
    # if "splice" in variant.consequence:

    # Caveat 4: Use caution in the case of multiple transcripts
    # If there are multiple transcripts we ensure that the primary transcript consequence is seen in at least one other
    # transcript.
    same_cons = False
    if len(variant.transcriptList) < 2:
        same_cons = True
    for transcript in variant.transcriptList[1:]:
        cons = ",".join(transcript['consequence'])
        if variant.consequence == cons:
            same_cons = True
    if not same_cons:
        return 0, [(0, f"Out of {len(variant.transcriptList)} transcripts, variant consequence {variant.consequence} in " + \
                f"transcript {variant.transcriptList[0]['transcript']} was only seen once.")]

    # No caveat's apply - the very strong pathogenic classifier is applicable.
    score = 1
    psv1 = f'PSV1 ({score}): Null variant consequence {variant.consequence} in gene {variant.geneName} with GNOMAD PLI {lof_gene_to_pli[variant.geneName]}'
    return 1, [(score, psv1)]

def get_ps1(gene_name, aa_mut, gene_mut_to_data):
    """
    PS1    Same amino acid change as a previously established pathogenic variant
        regardless of nucleotide change
        Example:    Val->Leu caused by either G>C or G>T in the same codon
        Caveat:    Beware of changes that impact splicing rather than at the
        amino acid/protein level

    """
    gene_mut = (gene_name, aa_mut)
    score = 0
    ps1 = ""
    if gene_mut_to_data.get(gene_mut):
        rs_id, criteria, signif = gene_mut_to_data[gene_mut]
        score = 1
        if clinvar_review_status_to_level[criteria] > 3: # Score variants with high review levels higher
            score += 1
        if "likely" not in signif: # Score pathogenic variants higher
            score += 1
        ps1 = f'PS1 ({score}): Amino acid change {aa_mut} in gene {gene_name} is associated with Clinvar {signif} variant rs{rs_id} with review status {criteria}'
    return score, ps1


def get_ps2():
    """
    Future goal: This will require patient EHR integration to automate
    
    PS2     De novo (both maternity and paternity confirmed) in a patient with the
        disease and no family history

    This requires manual data injection, or pulling data from a patients EHR.
    """
    return ''


def get_ps3(variant, gloc_to_pubmed_id_list):
    """
    Functional study database integration

    PS3     Well-established in vitro or in vivo functional studies supportive of a
        damaging effect on the gene or gene product
    """
    score = 0
    ps3 = ""
    start_pos = int(variant.position)
    end_pos = start_pos + len(variant.refAllele)
    gloc = (variant.chromosome, str(start_pos), str(end_pos))
    all_pubmed_ids = gloc_to_pubmed_id_list.get(gloc, "")
    total_pubmed_count = len(all_pubmed_ids)
    sorted_pubmed_ids = sorted(list(all_pubmed_ids))
    if total_pubmed_count > 3: # NOTE: This ensures the variant is called pathogenic! 
        score = 2
        pubmed_str = ', '.join(sorted_pubmed_ids)
        ps3 = f"PS3 ({score}): Variants pathogenicity is supported by multiple pubmed studies {pubmed_str} as established by AVADA"
    elif total_pubmed_count > 0:
        score = 1
        pubmed = sorted_pubmed_ids[0]
        ps3 = f"PS3 ({score}): Variants pathogenicity is supported by pubmed study {pubmed} as established by AVADA"
    return score, ps3


def get_ps4(variant, chrom_to_pos_to_gwas_data):
    """
    PS4     The prevalence of the variant in affected individuals is significantly
        increased compared to the prevalence in controls

    Note 1: Relative risk or OR, as obtained from case–control studies, is >5.0, and the confidence interval around
            the estimate of relative risk or OR does not include 1.0.
    """
    score = 0
    ps4 = ""
    if chrom_to_pos_to_gwas_data.get(variant.chromosome):
        if chrom_to_pos_to_gwas_data[variant.chromosome].get(int(variant.position)):
            or_value, ci_lower, ci_upper, pubmed_id, trait, p_value, _ = chrom_to_pos_to_gwas_data[variant.chromosome][int(variant.position)]
            ps4_text = f"PS4 ({score}): Genomic location {variant.chromosome}:{variant.position} is associated with trait {trait}. It has OR {or_value} " \
                    + f"with p-value {p_value} and CI {ci_upper}-{ci_lower} as reported in PubMed {pubmed_id}"
            if or_value > 5.0 and p_value < 1e-4:
                score = 3
                ps4 = ps4_text
            elif 2.0 < or_value <= 5.0 and p_value < 0.05:
                score = 2
                ps4 = ps4_text
            elif 1.5 < or_value <= 2.0 and p_value < 0.05:
                score = 1
                ps4 = ps4_text
    return score, ps4


def get_ps(variant, gene_mut_to_data, chrom_to_pos_to_gwas_data, gloc_to_pubmed_id_list):
    """
    Strong evidence of pathogenicity
        PS1    Same amino acid change as a previously established pathogenic variants
            regardless of nucleotide change
            Example:    Val->Leu caused by either G>C or G>T in the same codon
            Caveat:    Beware of changes that impact splicing rather than at the
            amino acid/protein level
        PS2    De novo (both maternity and paternity confirmed) in a patient with the
            disease and no family history
            Note: Confirmation of paternity only is insufficient. Egg donation, surrogate
            motherhood, errors in embryo transfer, etc. can contribute to non-
            maternity
        PS3    Well-established in vitro or in vivo functional studies supportive of a
            damaging effect on the gene or gene product
            Note: Functional studies that have been validated and shown to be
            reproducible and robust in a clinical diagnostic laboratory setting are
            considered the most well-established
        PS4    The prevalence of the variants in affected individuals is significantly
            increased compared to the prevalence in controls
            Note 1: Relative risk (RR) or odds ratio (OR), as obtained from case-control
            studies, is >5.0 and the confidence interval around the estimate of RR or OR
            does not include 1.0. See manuscript for detailed guidance.
            Note 2: In instances of very rare variants where case-control studies may
            not reach statistical significance, the prior observation of the variants in
            multiple unrelated patients with the same phenotype, and its absence in
            controls, may be used as moderate level of evidence.
    """
    score_1, ps1 = get_ps1(variant.geneName, variant.protein_variant, gene_mut_to_data)
    ps2 = get_ps2()
    score_3, ps3 = get_ps3(variant, gloc_to_pubmed_id_list)
    score_4, ps4 = get_ps4(variant, chrom_to_pos_to_gwas_data)

    aggregate_score = score_1 + score_3 + score_4
    ps_rationale_list = [(score_1, ps1), (0, ps2), (score_3, ps3), (score_4, ps4)]
    return aggregate_score, ps_rationale_list


def get_pm1(variant, pathogenic_domains):
    """    
    PM1    Located in a mutational hot spot and/or critical and well-established
            functional domain (e.g. active site of an enzyme) without benign variation
    """
    score = 0
    pm1 = ""
    v_pos = int(variant.position)
    for domain in pathogenic_domains:
        d_chrom, d_start, d_end, d_name, d_percent, d_score = domain
        v_chrom = variant.chromosome
        if v_chrom == d_chrom and int(d_start) <= v_pos <= int(d_end):
            variant.domain = d_name
            score = 1
            # The variant score indicates how many pathogenic variants were found in this domain and their relative authority.
            # In domains with an overwhelming amount of evidence supporting pathogenicity, adjust the score accordingly
            if d_score > 30 and d_percent > .9:
                score = 2
            if d_score > 50 and d_percent > .95:
                score = 3
            pm1 = f"PM1 ({score}): Variant is in domain, {d_name} where ({d_percent} %) variants are associated with pathogenicity (var score {d_score})."
    return score, pm1


def get_pm2(variant, population = "all", af_threshold = .01):
    """
    PM2    Absent from controls (or at extremely low frequency if recessive) (see Table 6)
        in Exome Sequencing Project, 1000 Genomes or ExAC

        Caveat: Population data for indels may be poorly called by next generation
        sequencing

    GNOMAD data is annotated with ICA/NIRVANA, the raw data source was published under this citation

    Karczewski, K. J., Francioli, L. C., Tiao, G., Cummings, B. B., Alföldi, J., Wang, Q., ... & Daly, M. J. (2020).
        The mutational constraint spectrum quantified from variation in 141,456 humans. Nature, 581(7809), 434-443.
        DOI: 10.1038/s41586-020-2308-7
    """
    score = 0
    if not variant.gnomad and not variant.oneKg:
        score = 1
        pm2 = f"PM2 ({score}): Variant is missing from GNOMAD and 1000 Genomes"
        return score, pm2
    
    af_key = f'{population}Af'
    allele_frequency = variant.gnomad.get(af_key, -1)  # Default to -1 if not found
    pm2 = ""
    # Check if allele frequency data is available and assess against threshold
    if allele_frequency < 0:
        pm2 = ""
    elif allele_frequency < af_threshold:
        score = 1
        pm2 = f"PM2 ({score}): GNOMAD general population AF={allele_frequency*100:.2f}% is below the threshold of {af_threshold*100}%"
    if variant.oneKg:
        onekg_af = variant.oneKg.get('allAf')
        if onekg_af < af_threshold:
            score = 1
            pm2 = f"PM2 ({score}): 1000 Genomes general population AF={allele_frequency*100:.2f}% is below the threshold of {af_threshold*100}%"
    return score, pm2


def get_pm3():
    """
    Trans variant mapping, only applies when two variants are at the same position and both have ~.5 af, indicating
    that the patient has one variant on each allele. 
    """
    return ""


def get_pm4(variant, chrom_to_repeat_regions):
    """
    PM4    Protein length changes due to in-frame deletions/insertions in a non-repeat
            region or stop-loss variants 

    The deletion or insertion of one or more amino acids as well as the extension of a protein by changing the stop
    codon to an amino acid codon (e.g. a stop loss variant) is more likely to disrupt protein function as compared 
    to a missense change alone due to length changes in the protein. Therefore, in-frame deletions/insertions and
    stop losses are considered moderate evidence of pathogenicity. The larger the deletion, insertion or extension,
    and the more conserved the amino acids are in a deleted region, then the more substantial is the evidence to
    support pathogenicity. In contrast, small in-frame deletions/insertions in repetitive regions, or regions that
    are not well conserved in evolution, are less likely to be pathogenic.

    Longer inframe insertion and deletions are more likely to be pathogenic
    Cannon, S., Williams, M., Gunning, A. C., & Wright, C. F. (2023). Evaluation of in silico pathogenicity prediction
        tools for the classification of small in-frame indels. BMC Medical Genomics, 16, Article 36.
        https://doi.org/10.1186/s12920-023-01454-6
    """
    # Conservation regions be overlayed with the repeat regions to find highly conserved
    # non repeat regions. In frame insertion/deletions in these regions would be more likely to be pathogenic
    pm4 = ""
    score = 0
    
    # The variant should be an in frame insertion or deletion
    pm4_consequences = ['stop_lost', 'inframe_insertion', 'inframe_deletion']
    if not variant.consequence in pm4_consequences:
        return score, pm4

    # Check if the variant is in a repeat region
    repeat_region_list = chrom_to_repeat_regions[variant.chromosome]
    pos = int(variant.position)
    length = abs(len(variant.refAllele) - len(variant.altAllele))
    in_repeat_region = False
    for repeat_region in repeat_region_list:
        start, stop, _ = repeat_region
        if start <= (pos + length) <= stop:
            in_repeat_region = True
    
    # For non repeat regoin variants, update the score based on the INDEL length. Longer indels are more pathogenic
    if not in_repeat_region:
        if length < 10: # 3 amino acids
            score = 1
        elif length < 30: # 4-10 amino acides
            score = 2
        else: # > 10 amino acids
            score = 3
        pm4 = f"PM4 ({score}): {variant.consequence} in a non-repeat region"
        if "insertion" in variant.consequence:
            pm4 = f"PM4 ({score}): {variant.consequence} of {variant.altAllele} in a non-repeat region"
        if "deletion" in variant.consequence:
            pm4 = f"PM4 ({score}): {variant.consequence} of {variant.refAllele} in a non-repeat region"
    return score, pm4


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


def get_pm5(variant, gene_aa_to_var_data, gene_mut_to_data):
    """
    PM5     Novel missense change at an amino acid residue where a different
        missense change determined to be pathogenic has been seen before
    Example: Arg156His is pathogenic; now you observe Arg156Cys

    Caveat: Beware of changes that impact splicing rather than at the amino
        acid/protein level
    """
    pm5 = ""
    score = 0
    # Variant must be in a gene to have an AA change associated
    if not variant.geneName:
        return score, pm5

    # If the full variant change was seen then PS1 was applied and PM5 cannot be applied
    gene_mut = (variant.geneName, variant.protein_variant)
    if gene_mut_to_data.get(gene_mut):
        return score, pm5

    # Check if the base aa and position are associated with a known pathogenic variant
    variant_aa_comp = get_aa_comparator(variant.protein_variant)
    gene_aa = (variant.geneName, variant_aa_comp)
    if gene_aa_to_var_data.get(gene_aa):
        ref_mut, rs_id, criteria, signif = gene_aa_to_var_data[gene_aa]
        score = 1
        if clinvar_review_status_to_level[criteria] > 3: # Score variants with high review levels higher
            score += 1
        if "likely" not in signif: # Score pathogenic variants higher
            score += 1
        pm5 = f"PM5 ({score}): Variant AA change {variant.protein_variant} is the same base AA as {ref_mut} a {signif} variant with rs id {rs_id} and review status {criteria}."
    return score, pm5


def get_pm6():
    """
    Implied De Novo variant - difficult without phenotype and/or parental genotype/phenotype
    information
    """
    return ""


def get_pm(variant, chrom_to_repeat_regions, gene_aa_to_var_data, pathogenic_domains, gene_mut_to_data):
    """
    Moderate evidence of pathogenicity
        PM1    Located in a mutational hot spot and/or critical and well-established
            functional domain (e.g. active site of an enzyme) without benign variation
        PM2    Absent from controls (or at extremely low frequency if recessive) (see Table 6)
            in Exome Sequencing Project, 1000 Genomes or ExAC
            Caveat: Population data for indels may be poorly called by next generation
            sequencing
        PM3    For recessive disorders, detected in trans with a pathogenic variants 
            Note: This requires testing of parents (or offspring) to determine phase
        PM4    Protein length changes due to in-frame deletions/insertions in a non-repeat
            region or stop-loss variants 
        PM5    Novel missense change at an amino acid residue where a different
            missense change determined to be pathogenic has been seen before
            Example: Arg156His is pathogenic; now you observe Arg156Cys
            Caveat: Beware of changes that impact splicing rather than at the amino
            acid/protein level
        PM6    Assumed de novo, but without confirmation of paternity and maternity
    """

    score_1, pm1 = get_pm1(variant, pathogenic_domains)
    score_2, pm2 = get_pm2(variant)
    pm3 = get_pm3()
    score_4, pm4 = get_pm4(variant, chrom_to_repeat_regions)
    score_5, pm5 = get_pm5(variant, gene_aa_to_var_data, gene_mut_to_data)
    pm6 = get_pm6()

    aggregate_score = score_1 + score_2 + score_4 + score_5
    pm_rationale_list = [(score_1, pm1), (score_2, pm2), (0, pm3), (score_4, pm4), (score_5, pm5), (0, pm6)]
    return aggregate_score, pm_rationale_list


def get_pp1():
    """
    Future goal: This will require patient EHR integration to automate as it requires patient family genetic information

    PP1    Co-segregation with disease in multiple affected family members in a gene
        definitively known to cause the disease
        Note: May be used as stronger evidence with increasing segregation data
    """
    return ""


def get_pp2(variant, missense_pathogenic_genes):
    """
    PP2    Missense variants in a gene that has a low rate of benign missense variation
        and where missense variants are a common mechanism of disease
    """
    pp2 = ""
    score = 0
    if variant.geneName:
        if variant.geneName in missense_pathogenic_genes:
            if "missense" in variant.consequence:
                score = 1
                # TODO: These variant justifications should include the pathogenic and benign % calculated
                # for this gene
                pp2 = f"PP2 ({score}): Missense variant consequence {variant.consequence} in {variant.geneName} which was found to have a low rate of "\
                        + "missense variation and where missense variants are a common mechanism of disease" 
    return score, pp2


def get_pp3(variant):
    """
    PP3    Multiple lines of computational evidence support a deleterious effect on
            the gene or gene product (conservation, evolutionary, splicing impact, etc)
            Caveat: As many in silico algorithms use the same or very similar input for
            their predictions, each algorithm should not be counted as an independent
            criterion. PP3 can be used only once in any evaluation of a variant.

    GERP, REVEL and DANN data is annotated with ICA/NIRVANA, the raw data sources were published under these citations

    Davydov, E. V., et al. (2010). Identifying a high fraction of the human genome to be under selective constraint
        using GERP++. PLoS Computational Biology, 6(12), e1001025.

    Quang, D., Chen, Y., & Xie, X. (2015). DANN: a deep learning approach for annotating the pathogenicity of genetic
        variants. Bioinformatics, 31(5), 761-763.

    Ioannidis, N. M., et al. (2016). REVEL: An ensemble method for predicting the pathogenicity of rare missense
        variants. American Journal of Human Genetics, 99(4), 877-885.
    """
    # TODO In depth analysis of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9748256/pdf/main.pdf and update cutoffs and weights
    pp3 = ""
    score = 0
    revel = False
    if variant.revel:
        if variant.revel >= .5:
            score += 1
            revel = True
    dann = False
    if variant.dann:
        if variant.dann >= .7:
            score += 1
            dann = True
    gerp = False
    if variant.gerp:
        if variant.gerp >= 2.0:
            score += 1
            gerp = True

    if revel or dann or gerp:
        pp3 = f"PP3 ({score}): Computational evidence support a deleteterious effect; gerp {variant.gerp} | dann {variant.dann} | revel {variant.revel}" 
    return score, pp3


def get_pp4():
    """
    Future goal: This will require patient EHR integration to automate

    PP4    Patient’s phenotype or family history is highly specific for a disease with a
            single genetic etiology
    """
    return ""


def get_pp5(variant):
    """
    PP5    Reputable source recently reports variants as pathogenic but the evidence is
        not available to the laboratory to perform an independent evaluation
    """
    # Map clinvar review criteria to an integer
    pp5 = ""
    score = 0
    
    # Weight variants that were reviewed to a higher degree more heavily
    value = clinvar_review_status_to_level.get(variant.clinvar_review_status, 0)
    if value > 0 and "pathogenic" in variant.clinvar_significance:
        score = value
        pp5 = f"PP5 ({score}): Variant was found reported as {variant.clinvar_significance} in ClinVar with review status of {variant.clinvar_review_status}."
    return score, pp5


def get_pp(variant, missense_pathogenic_genes):
    """
    Supporting evidence of pathogenicity
        PP1    Co-segregation with disease in multiple affected family members in a gene
            definitively known to cause the disease
            Note: May be used as stronger evidence with increasing segregation data
        PP2    Missense variants in a gene that has a low rate of benign missense variation
            and where missense variants are a common mechanism of disease
        PP3    Multiple lines of computational evidence support a deleterious effect on
            the gene or gene product (conservation, evolutionary, splicing impact, etc)
            Caveat: As many in silico algorithms use the same or very similar input for
            their predictions, each algorithm should not be counted as an independent
            criterion. PP3 can be used only once in any evaluation of a variant.
        PP4    Patient’s phenotype or family history is highly specific for a disease with a
            single genetic etiology
        PP5    Reputable source recently reports variants as pathogenic but the evidence is
            not available to the laboratory to perform an independent evaluation
    """
    pp1 = get_pp1()
    score_2, pp2 = get_pp2(variant, missense_pathogenic_genes)
    score_3, pp3 = get_pp3(variant)
    pp4 = get_pp4()
    score_5, pp5 = get_pp5(variant)

    aggregate_score = score_2 + score_3 + score_5
    pp_rationale_list = [(0, pp1), (score_2, pp2), (score_3, pp3), (0, pp4), (score_5, pp5)]
    return aggregate_score, pp_rationale_list
