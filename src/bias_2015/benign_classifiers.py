"""
BIAS-2015 implementation of the ACMG 2015 germline benign classifiers
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


def get_ba1(variant):
    """
    Stand-Alone evidence of benign impact
        BA1    Allele frequency is above 5% in Exome Sequencing Project, 1000 Genomes,
        or ExAC
    
    ** ExAC has been replaced by Gnomad **
    """
    ba1 = ""
    score = 0
    if variant.oneKg:
        onekg_af = variant.oneKg.get('allAf')
        if onekg_af > 0.05:
            score = 1
            ba1 = f"BA1 ({score}): One Thousand Genome general population AF={onekg_af*100:.2f}% which is greater than 5%."
    elif variant.gnomad:
        gnomad_af = variant.gnomad.get('allAf', 0)
        if gnomad_af > 0.05:
            score = 1
            ba1 = f"BA1 ({score}): Gnomad general population AF={gnomad_af*100:.2f}% which is greater than 5%."
    return score, ba1


def get_ba(variant):
    """
    AMP/ACMG
    Stand-Alone evidence of benign impact
        BA1    Allele frequency is above 5% in Exome Sequencing Project, 1000 Genomes,
        or ExAC

    ** ExAC has been replaced by Gnomad **
    """
    ba1_score, ba1 = get_ba1(variant)
    aggregate_score = ba1_score
    ba_rationale_list = [(ba1_score, ba1)]
    return aggregate_score, ba_rationale_list


def get_bs1(variant):
    """
    BS1    Allele frequency is greater than expected for disorder (see table 6)
    """
    bs1 = ""
    score = 0
    if variant.oneKg:
        onekg_af = variant.oneKg.get('allAf')
        if 0.05 >= onekg_af > 0.01:
            score = 1
            bs1 = f"BS1 ({score}): One Thousand Genome general population AF={onekg_af*100:.2f}% which is between 1% and 5%."
    if variant.gnomad:
        gnomad_af = variant.gnomad.get('allAf', 0)
        if 0.05 >= gnomad_af > 0.01:
            score = 1
            bs1 = f"BS1 ({score}): Gnomad general population AF={gnomad_af*100:.2f}% which is between 1% and 5%."
    return score, bs1

# TODO
def get_bs2():
    """
    BS2    Observed in a healthy adult individual for a recessive (homozygous),
        dominant (heterozygous), or X-linked (hemizygous) disorder with full
        penetrance expected at an early age
    
    Looking at https://hpo.jax.org/ for disease:gene linkage and early onset genes
    and at https://search.clinicalgenome.org/kb/downloads for autosomal recessive genes

    Need to find a way to map diseases between the two tables (MONDO to OMIM?)
    """
    return ""


def get_bs3():
    """
    BS3    Well-established in vitro or in vivo functional studies shows no damaging
        effect on protein function or splicing

    NLP literature scraping (akin to AVADA) could be applied, however none was identified currently (July 2024)
    """
    return ""


def get_bs4():
    """
    BS4    Lack of segregation in affected members of a family
        Caveat: The presence of phenocopies for common phenotypes (i.e. cancer,
        epilepsy) can mimic lack of segregation among affected individuals. Also,
        families may have more than one pathogenic variants contributing to an
        autosomal dominant disorder, further confounding an apparent lack of
        segregation.
    """
    return ""


def get_bs(variant):
    """
    Strong evidence of benign impact
        BS1    Allele frequency is greater than expected for disorder (see table 6)
        BS2    Observed in a healthy adult individual for a recessive (homozygous),
            dominant (heterozygous), or X-linked (hemizygous) disorder with full
            penetrance expected at an early age
        BS3    Well-established in vitro or in vivo functional studies shows no damaging
            effect on protein function or splicing
        BS4    Lack of segregation in affected members of a family
            Caveat: The presence of phenocopies for common phenotypes (i.e. cancer,
            epilepsy) can mimic lack of segregation among affected individuals. Also,
            families may have more than one pathogenic variants contributing to an
            autosomal dominant disorder, further confounding an apparent lack of
            segregation.
    """

    bs1_score, bs1 = get_bs1(variant)
    bs2 = get_bs2()
    bs3 = get_bs3()
    bs4 = get_bs4()

    aggregate_score = bs1_score
    bs_rationale_list = [(bs1_score, bs1), (0, bs2), (0, bs3), (0, bs4)]
    return aggregate_score, bs_rationale_list


def get_bp1(variant, truncating_genes):
    """
    BP1    Missense variants in a gene for which primarily truncating variants are
        known to cause disease
    """
    score = 0
    bp1 = ""
    if variant.geneName and variant.consequence:
        if variant.geneName in truncating_genes:
            if "missense" in variant.consequence:
                score = 1
                # TODO: These variant justifications should include the pathogenic and benign % calculated
                # for this gene
                bp1 = f"BP1: Missense variant type {variant.consequence} in gene {variant.geneName} which has over 80% truncating pathogenic variants"  
    return score, bp1


def get_bp2():
    """    
    BP2    Observed in trans with a pathogenic variants for a fully penetrant dominant
            gene/disorder; or observed in cis with a pathogenic variants in any
            inheritance pattern
    """
    bp2 = ""
    return bp2


def get_bp3(variant, chrom_to_repeat_regions):
    """
    BP3    In-frame deletions/insertions in a repetitive region without a known
        function
    """
    score = 0
    bp3 = ""
    var_len = abs(len(variant.refAllele.replace("-", "")) - len(variant.altAllele.replace("-", "")))
    # Exclude missense and other SNP
    if var_len <3 and variant.consequence != 'stop_lost':
        return score, bp3
    # Exclude indels that are not in frame (length is divisible by 3)
    if var_len % 3 != 0 and variant.consequence != 'stop_lost':
        return score, bp3
   
    # Determine if the variant lands within a repeat region
    repeat_region_list = chrom_to_repeat_regions.get(variant.chromosome)
    if not repeat_region_list:
        return score, bp3
    pos = int(variant.position)
    in_repeat_region = False
    seen_repeat_region = ""
    for repeat_region in repeat_region_list:
        start, stop, _ = repeat_region
        if start <= pos <= stop:
            in_repeat_region = True
            seen_repeat_region = (start, stop)

    # TODO: Use conservation data if applicable, PolyPhen
    if in_repeat_region:
        if var_len % 3 == 0:
            if var_len < 60:  
                score = 1
                if var_len < 30:
                    score += 1
                if var_len < 5:
                    score += 1
                bp3 = f"BP3 ({score}): In-frame INDEL of length {var_len} in repeat region {variant.chromosome} {seen_repeat_region[0]}-{seen_repeat_region[1]}"
        else:
            bp3 = f"BP3 ({score}): Consequence {variant.consequence} in repeat region {variant.chromosome} {seen_repeat_region[0]}-{seen_repeat_region[1]}"
    return score, bp3


def get_bp4(variant):
    """    
    BP4    Multiple lines of computational evidence suggest no impact on gene or
            gene product (conservation, evolutionary, splicing impact, etc)
            Caveat: As many in silico algorithms use the same or very similar input for
            their predictions, each algorithm cannot be counted as an independent
            criterion. BP4 can be used only once in any evaluation of a variants.
    """
    bp4 = ""
    score = 0
    revel = False
    if variant.revel:
        if variant.revel <= .25:
            revel = True
            score += 1
    dann = False
    if variant.dann:
        if variant.dann <= .25:
            score += 1
            dann = True
    gerp = False
    if variant.gerp:
        if variant.gerp <= 2.0:
            score += 1
            gerp = True
    if revel or dann or gerp:
        bp4 = f"BP4 ({score}): Computational evidence support a benign effect; gerp {variant.gerp} | dann {variant.dann} | revel {variant.revel}" 
    return score, bp4


def get_bp5():
    """    
    BP5    Variant found in a case with an alternate molecular basis for disease
    """
    bp5 = ""
    return bp5


def get_bp6(variant):
    """
    BP6    Reputable source recently reports variants as benign but the evidence is not
        available to the laboratory to perform an independent evaluation
    """
    score = 0
    bp6 = ""
    value = clinvar_review_status_to_level.get(variant.clinvar_review_status, 0)
    if value > 0 and "benign" in variant.clinvar_significance:
        score = value
        bp6 = f"BP6 ({score}): Variant was found in ClinVar as {variant.clinvar_significance} with review status "+\
                f"of {variant.clinvar_review_status} and given a weighted PP5 value of {value}"
    return score, bp6


# TODO Splice and conservation
def get_bp7(variant):
    """
    BP7    A synonymous (silent) variants for which splicing prediction algorithms
            predict no impact to the splice consensus sequence nor the creation of a
            new splice site AND the nucleotide is not highly conserved
    """
    score = 0
    bp7 = ""
    if "synonymous" in variant.consequence and "splice" not in variant.consequence:
        score = 1
        bp7 = f"BP7 ({score}): Variant has synonymous associated consequence {variant.consequence}"
    elif "intron" in variant.consequence and "splice" not in variant.consequence:
        score = 1
        bp7 = f"BP7 ({score}): Variant has intronic associated consequence {variant.consequence}"
    elif variant.geneName == "intergenic":
        score = 1
        bp7 = f"BP7 ({score}): Intergenic variant" 

    return score, bp7


def get_bp(variant, truncating_genes, chrom_to_repeat_regions):
    """
    Supporting evidence of benign impact
        BP1    Missense variants in a gene for which primarily truncating variants are
            known to cause disease
        BP2    Observed in trans with a pathogenic variants for a fully penetrant dominant
            gene/disorder; or observed in cis with a pathogenic variants in any
            inheritance pattern
        BP3    In-frame deletions/insertions in a repetitive region without a known
            function
        BP4    Multiple lines of computational evidence suggest no impact on gene or
            gene product (conservation, evolutionary, splicing impact, etc)
            Caveat: As many in silico algorithms use the same or very similar input for
            their predictions, each algorithm cannot be counted as an independent
            criterion. BP4 can be used only once in any evaluation of a variants.
        BP5    Variant found in a case with an alternate molecular basis for disease
        BP6    Reputable source recently reports variants as benign but the evidence is not
            available to the laboratory to perform an independent evaluation
        BP7    A synonymous (silent) variants for which splicing prediction algorithms
            predict no impact to the splice consensus sequence nor the creation of a
            new splice site AND the nucleotide is not highly conserved
    """
    bp1_score, bp1 = get_bp1(variant, truncating_genes)
    bp2 = get_bp2()
    bp3_score, bp3 = get_bp3(variant, chrom_to_repeat_regions)
    bp4_score, bp4 = get_bp4(variant)
    bp5 = get_bp5()
    bp6_score, bp6 = get_bp6(variant)
    bp7_score, bp7 = get_bp7(variant)

    aggregate_score = bp1_score + bp3_score + bp4_score + bp6_score + bp7_score 
    bp_rationale_list = [(bp1_score, bp1), (0, bp2), (bp3_score, bp3), (bp4_score, bp4), (0, bp5), (bp6_score, bp6), (bp7_score, bp7)]
    return aggregate_score, bp_rationale_list
