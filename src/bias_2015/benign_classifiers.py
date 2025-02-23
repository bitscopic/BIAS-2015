"""
BIAS-2015 implementation of the ACMG 2015 germline benign classifiers
"""
from src.bias_2015.constants import clinvar_review_status_to_level, score_to_hum_readable, loeuf_thresholds, benign_thresholds

def get_ba1(variant):
    """
    BA1: Stand-Alone Evidence of Benign Impact
         A variant is classified as BA1 if its allele frequency (AF) is above a gene-specific threshold.
         The default threshold is 5%, but this is adjusted based on LOEUF from GnomAD.
    
    BIAS uses LOEUF modulated AF cutoff values, see constants.py for the values set. 

    ** ExAC has been replaced by GnomAD **
    """
    ba1 = ""
    score = 0

    # Retrieve allele frequencies safely
    onekg_af = variant.oneKg.get('allAf', 0) if variant.oneKg else 0
    gnomad_af = variant.gnomad.get('allAf', 0) if variant.gnomad else 0
    highest_af = max(onekg_af, gnomad_af)

    loeuf = variant.gene_gnomad.get('loeuf', 1) if hasattr(variant, 'gene_gnomad') and variant.gene_gnomad else 1

    # Determine BA1 cutoff based on LOEUF
    pop_thresholds = {}
    for loeuf_tier in loeuf_thresholds:
        if loeuf < loeuf_tier['max_loeuf']:
            pop_thresholds = loeuf_tier
            break
    ba1_cutoff = pop_thresholds['ba1_cutoff']
    
    loeuf_str = f"({variant.geneName}:{loeuf})"
    # Check if the variant exceeds the BA1 threshold
    if highest_af > ba1_cutoff:
        score = 5
        source = "both One Thousand Genomes & GnomAD" if onekg_af > ba1_cutoff and gnomad_af > ba1_cutoff else \
                 "One Thousand Genomes" if onekg_af > ba1_cutoff else "GnomAD"
        ba1 = f"BA1: {source} general population AF={highest_af*100:.5f}% exceeds " + \
                f"LOEUF{loeuf_str}-based threshold ({ba1_cutoff*100:.5f}%) based on LOEUF={loeuf:.5f}."

    return score, ba1

def get_ba(variant, skip_list):
    """
    AMP/ACMG
    Stand-Alone evidence of benign impact
        BA1    Allele frequency is above 5% in Exome Sequencing Project, 1000 Genomes,
        or ExAC

    ** ExAC has been replaced by Gnomad **
    """
    if 'ba1' not in skip_list:
        ba1_score, ba1 = get_ba1(variant)
    else:
        ba1_score = 0
        ba1 = ""
    ba_rationale_list = {'ba1': (ba1_score, ba1)}
    return ba_rationale_list

def get_bs1(variant):
    """
    BS1: Allele frequency is greater than expected for disorder.
    
    - BS1 applies if AF is above expected but below BA1 cutoff.
    - For very low cutoffs, a minimum allele count (AC) is also required.
    
    evRepo assigns the following weights:
    - BS1 (3) = Strong
    - BS1_Supporting (1) = Supporting
    """
    bs1 = ""
    score = 0

    # Retrieve allele frequencies and counts safely
    onekg_af = variant.oneKg.get('allAf', 0) if variant.oneKg else 0
    gnomad_af = variant.gnomad.get('allAf', 0) if variant.gnomad else 0
    highest_af = max(onekg_af, gnomad_af)

    # Retrieve LOEUF safely
    loeuf = variant.gene_gnomad.get('loeuf', 1.0) if hasattr(variant, 'gene_gnomad') and variant.gene_gnomad else 1.0

    # Determine BS1 cutoffs based on LOEUF
    pop_thresholds = {}
    for loeuf_tier in loeuf_thresholds:
        if loeuf < loeuf_tier['max_loeuf']:
            pop_thresholds = loeuf_tier
            break
    ba1_cutoff = pop_thresholds['ba1_cutoff']
    bs1_cutoff_strong = pop_thresholds['bs1_cutoff_strong']
    bs1_cutoff_supporting= pop_thresholds['bs1_cutoff']
    
    # Ensure a variant **cannot** qualify for both BA1 and BS1
    if highest_af >= ba1_cutoff:
        return 0, ""  # Variant already meets BA1 criteria, so BS1 should not be applied

    loeuf_str = f"({variant.geneName}:{loeuf})"
    # Apply BS1 based on gene-specific cutoffs and AC validation
    if highest_af > bs1_cutoff_strong:
        score = 3  # BS1 (Strong)
        source = "both One Thousand Genomes & GnomAD" if onekg_af > bs1_cutoff_supporting and gnomad_af > bs1_cutoff_supporting else \
                 "One Thousand Genomes" if onekg_af > bs1_cutoff_supporting else "GnomAD"
        bs1 = f"BS1: {source} general population AF={highest_af*100:.5f}% is between LOEUF{loeuf_str}-based thresholds " + \
                f"{bs1_cutoff_supporting*100:.5f}% and {bs1_cutoff_strong*100:.5f}%"

    elif highest_af > bs1_cutoff_supporting:
        score = 1  # BS1_Supporting
        source = "both One Thousand Genomes & GnomAD" if onekg_af > bs1_cutoff_supporting and gnomad_af > bs1_cutoff_supporting else \
                 "One Thousand Genomes" if onekg_af > bs1_cutoff_supporting else "GnomAD"
        bs1 = f"BS1_{score_to_hum_readable[score]}: {source} general population AF={highest_af*100:.5f}%" + \
                f" exceeds LOEUF{loeuf_str}-based threshold {bs1_cutoff_supporting*100:.5f}%"

    return score, bs1

def get_bs2(variant):
    """
    BS2: Observed in a healthy adult individual for a disorder with full penetrance expected at an early age.
    LOEUF-adjusted thresholds are applied based on gene constraint.
    """
    bs2 = ""
    score = 0

    if not variant.oneKg and not variant.gnomad:
        return score, bs2
    if not variant.clingen_gene_validity:
        return score, bs2

    controls_af, controls_ac, gnom_all_ac, all_homozygous, male_hc, female_ac, onekg_all_ac = 0, 0, 0, 0, 0, 0, 0
    if variant.gnomad:
        controls_ac = variant.gnomad.get("controlsAllAc", 0)
        controls_af = variant.gnomad.get("controlsAllAf", 0)
        gnom_all_ac = variant.gnomad.get("allAc", 0)
        all_homozygous = variant.gnomad.get("allHc", 0)
        male_hc = variant.gnomad.get("maleAc", 0)
        female_ac = variant.gnomad.get("femaleHc", 0)
    if variant.oneKg:
        onekg_all_ac = variant.oneKg.get("allAc", 0)

    loeuf = variant.gene_gnomad.get('loeuf', 1.0) if hasattr(variant, 'gene_gnomad') and variant.gene_gnomad else 1.0
    pop_thresholds = {}
    for loeuf_tier in loeuf_thresholds:
        if loeuf < loeuf_tier['max_loeuf']:
            pop_thresholds = loeuf_tier
            break
    recessive_homozygous_threshold = pop_thresholds["bs2_recessive_homozygous_threshold"]
    dominant_allele_threshold = pop_thresholds["bs2_dominant_allele_threshold"]
    xlinked_threshold = pop_thresholds["bs2_xlinked_threshold"]
    bs2_af_threshold = pop_thresholds["bs2_af_threshold"]
    loeuf_str = f"({variant.geneName}:{loeuf:.5f})"

    for validity_entry in variant.clingen_gene_validity:
        inheritance = validity_entry.get("inheritance", "").lower()
        disease = validity_entry.get("disease", "unknown disease")

        if "recessive" in inheritance and all_homozygous >= recessive_homozygous_threshold:
            score = 3 if all_homozygous > recessive_homozygous_threshold * 2 else 1
            bs2 = f"BS2_{score_to_hum_readable[score]}: Observed as homozygous ({all_homozygous}) in healthy individuals for autosomal recessive disease {disease} " \
                  f"exceeding LOEUF{loeuf_str}-based threshold ({recessive_homozygous_threshold})."

        elif "dominant" in inheritance:
            if loeuf > 0.7:
                continue
            if bs2_af_threshold and variant.gnomad and controls_af >= bs2_af_threshold:
                continue
            if (controls_ac >= dominant_allele_threshold or gnom_all_ac >= dominant_allele_threshold or onekg_all_ac >= dominant_allele_threshold):
                score = 3
                bs2 = f"BS2: Observed in {max(controls_ac, gnom_all_ac, onekg_all_ac)} healthy individuals for autosomal dominant disease {disease} " \
                      f"exceeding LOEUF{loeuf_str}-based threshold ({dominant_allele_threshold})."

        elif "x-linked" in inheritance and (male_hc >= xlinked_threshold or female_ac >= xlinked_threshold):
            score = 3
            bs2 = f"BS2: Observed in healthy males (hemizygous: {male_hc}) or females (homozygous: {female_ac}) for X-linked disease {disease} " \
                  f"exceeding LOEUF{loeuf_str}-based threshold ({xlinked_threshold})."

    return score, bs2

def get_bs3():
    """
    BS3    Well-established in vitro or in vivo functional studies shows no damaging
           effect on protein function or splicing

    NLP literature scraping (akin to AVADA) could be applied, however none was identified currently (July 2024) that
    can conclusively rule the literature derived variant as benign.
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


def get_bs(variant, skip_list):
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

    if "bs1" not in skip_list:
        bs1_score, bs1 = get_bs1(variant)
    else:
        bs1_score = 0
        bs1 = ""

    if "bs2" not in skip_list:
        bs2_score, bs2 = get_bs2(variant)
    else:
        bs2_score = 0
        bs2 = ""

    bs3 = get_bs3()
    bs4 = get_bs4()

    bs_rationale_list = {
        'bs1': (bs1_score, bs1),
        'bs2': (bs2_score, bs2),
        'bs3': (0, bs3),
        'bs4': (0, bs4)
    }
    
    return bs_rationale_list

def get_bp1(variant, truncating_gene_to_data):
    """
    BP1    Missense variants in a gene for which primarily truncating variants are
           known to cause disease
    
    Clinvar genes where >80% of pathogenic variants are truncating variants, and where
    <20% of benign variants are truncating variants were included. 
    """
    score = 0
    bp1 = ""
    if "missense" not in variant.consequence:
        return score, bp1

    data = truncating_gene_to_data.get(variant.geneName)
    if data:
        path_per, ben_per = data
        score = 1
        bp1 = f"BP1: Missense variant type {variant.consequence} in gene {variant.geneName} where {path_per} of pathogenic variants" + \
                f"are truncating variants, and {ben_per} of benign variants are truncating"

    return score, bp1


def get_bp2():
    """    
    BP2    Observed in trans with a pathogenic variants for a fully penetrant dominant
            gene/disorder; or observed in cis with a pathogenic variants in any
            inheritance pattern
    """
    bp2 = ""
    return bp2

# NOTE: This is duplicated in benign and pathogenic classifiers modules, it should be librarized if we keep
# seeing duplicated functions
def get_variant_indel_length(ref, alt):
    """
    Compute the true length of an indel by trimming common prefixes and suffixes.
    Returns the length of the actual insertion or deletion event.
    """
    # Trim common prefix
    while len(ref) > 0 and len(alt) > 0 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]

    # Trim common suffix
    while len(ref) > 0 and len(alt) > 0 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]

    # True length of the indel (in nucleotides)
    return max(len(ref), len(alt))

def get_bp3(variant, chrom_to_repeat_regions):
    """
    BP3    In-frame deletions/insertions in a repetitive region without a known
           function

    Repeat regions were idenitfied by intersecting UCSC repeat masker track with the UCSC
    consensus coding track. This results in repeat regions within coding regions.

    CCDS - Pruitt, K.D., Harrow, J., Harte, R.A., et al. (2009). The consensus coding sequence (CCDS) project: 
    Identifying a common protein-coding gene set for the human and mouse genomes. Genome Research, 19(7), 1316–1323.

    RepeatMasker - Smit, A.F.A., Hubley, R., & Green, P. RepeatMasker Open-4.0. 1996-2010.
    """
    score = 0
    bp3 = ""
    var_len = get_variant_indel_length(variant.refAllele, variant.altAllele)

    # Exclude missense and other SNP
    if var_len < 3 and variant.consequence != 'stop_lost':
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

    # Confirm that the variant is in a repeat region
    if in_repeat_region:
        # Confirm that the variant is an in-frame deletion
        if var_len % 3 == 0:
            if var_len < benign_thresholds['bp3_in_frame_max_length']: # 19 or less amino acids are deleted
                score = 1
                if var_len < benign_thresholds['bp3_in_frame_strong_length']: # less than 5 amino acids are deleted
                    score = 3
                if score == 1:
                    bp3 = f"BP3: In-frame INDEL of length {var_len} in repeat region " + \
                            f"{variant.chromosome} {seen_repeat_region[0]}-{seen_repeat_region[1]}"
                else:
                    bp3 = f"BP3_{score_to_hum_readable[score]}: In-frame INDEL of length {var_len} in repeat region " + \
                            f"{variant.chromosome} {seen_repeat_region[0]}-{seen_repeat_region[1]}"
    return score, bp3


def get_bp4(variant, chrom_to_pos_to_alt_to_splice_score):
    """    
    BP4    Multiple lines of computational evidence suggest no impact on gene or
           gene product (conservation, evolutionary, splicing impact, etc)

           Caveat: As many in silico algorithms use the same or very similar input for
                   their predictions, each algorithm cannot be counted as an independent
                   criterion. BP4 can be used only once in any evaluation of a variants.

    PhyloP, GERP, REVEL and DANN data are annotated with ICA/NIRVANA, the raw data sources were published under these citations

    Pollard KS, Hubisz MJ, Rosenbloom KR, Siepel A. "Detection of nonneutral substitution rates on mammalian phylogenies."
        Genome Research. 2010;20(1):110-121. doi: 10.1101/gr.097857.109.

    Davydov, E. V., et al. (2010). Identifying a high fraction of the human genome to be under selective constraint
        using GERP++. PLoS Computational Biology, 6(12), e1001025.

    Quang, D., Chen, Y., & Xie, X. (2015). DANN: a deep learning approach for annotating the pathogenicity of genetic
        variants. Bioinformatics, 31(5), 761-763.

    Ioannidis, N. M., et al. (2016). REVEL: An ensemble method for predicting the pathogenicity of rare missense
        variants. American Journal of Human Genetics, 99(4), 877-885.
    
    Cutoff values were used from "Calibration of computational tools for missense variant pathogenicity
    classification and ClinGen recommendations for PP3/BP4 criteria" - table 2
    https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/pdf/main.pdf

    ABSplice cutoff values were derived from UCSC values which the authors provide on github 
    Wagner N, Çelik MH, Hölzlwimmer FR, Mertes C, Prokisch H, Yépez VA, Gagneur J. Aberrant splicing prediction across 
        human tissues. Nat Genet. 2023 May;55(5):861-870. PMID: 37142848
    """
    bp4 = ""
    score = 0
    printout_text = []
    alg_to_score = {}

    if variant.phylopScore:
        if variant.phylopScore >= benign_thresholds['bp4_phylop_path_cutoff']:
            return score, bp4
        if variant.phylopScore <= benign_thresholds['bp4_phylop_low']:
            phylop_weight = "supporting"
            if variant.phylopScore <= benign_thresholds['bp4_phylop_very_low']:
                phylop_weight = "moderate"
                alg_to_score['phylop'] = 2
            else:
                alg_to_score['phylop'] = 1
            printout_text.append(f"{phylop_weight} phylop {variant.phylopScore}")
    
    if variant.revel: 
        if variant.revel >= benign_thresholds['bp4_revel_path_cutoff']: # Strong pathogenic 
            return score, bp4
        if variant.revel <= benign_thresholds['bp4_revel_supporting']:
            revel_weight = "supporting"
            if variant.revel <= benign_thresholds['bp4_revel_very_strong']: # Very strong benign
                alg_to_score['revel'] = 4
                revel_weight = "very strong"
            elif variant.revel <= benign_thresholds['bp4_revel_strong']: # strong
                alg_to_score['revel'] = 3
                revel_weight = "strong"
            elif variant.revel <= benign_thresholds['bp4_revel_moderate']: # moderate
                alg_to_score['revel'] = 2
                revel_weight = "moderate"
            else: # supporting
                alg_to_score['revel'] = 1 
            printout_text.append(f"{revel_weight} revel {variant.revel}")
    
    # The DANN (Deleterious Annotation of genetic variants using Neural Networks) score ranges from 0 to 1, where values
    # closer to 1 indicate a higher likelihood of the variant being deleterious.
    if variant.dann: #  
        if variant.dann <= benign_thresholds['bp4_dann_supporting']:
            dann_weight = "supporting"
            if variant.dann <= benign_thresholds['bp4_dann_strong']: # strong
                alg_to_score['dann'] = 3
                dann_weight = "strong"
            elif variant.dann <= benign_thresholds['bp4_dann_moderate']: # moderate
                alg_to_score['dann'] = 2
                dann_weight = "moderate"
            else: # supporting
                alg_to_score['dann'] = 1 
            printout_text.append(f"{dann_weight} dann {variant.dann}")

    # Conservation data, gerp
    if variant.gerp:
        if variant.gerp <= benign_thresholds['bp4_gerp_supporting']:
            gerp_weight = "supporting"
            if variant.gerp <= benign_thresholds['bp4_gerp_moderate']: # moderate
                alg_to_score['gerp'] = 2
                gerp_weight = "moderate"
            else: # supporting
                alg_to_score['gerp'] = 1
            printout_text.append(f"{gerp_weight} gerp {variant.gerp}")
    
    if "splice" in variant.consequence:
        # consider bases that are predicted to have a strong impact on splicing
        # Authors suggest .2 as a high end cutoff, .05 as a middle end, and .01 as a low https://github.com/gagneurlab/absplice?tab=readme-ov-file#output
        splice_score = None # 0 indicates no impact on splicing.
        if chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome):
            if chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome).get(variant.position):
                if chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome).get(variant.position).get(variant.altAllele):
                    splice_score = chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome).get(variant.position).get(variant.altAllele)
        if splice_score:
            if splice_score <= benign_thresholds['bp4_absplice_strong']:
                absplice_weight = "strong"
                alg_to_score['absplice'] = 3
                printout_text.append(f"{absplice_weight} ABSplice {splice_score}")
    
    if printout_text:
        formatted_printout_text = " | ".join(printout_text)
        lines_of_evidence = len(printout_text)
        # As per ACMG guidelines, only use one score - therefore average it and round down
        # Floor division is used to round down because we want to be conservative, increasing specificity
        if alg_to_score:
            score = 0
            best_score = 0
            best_algs = []
            # Identify the best score and its algorithm(s)
            for alg, a_score in alg_to_score.items():
                if a_score >= best_score:
                    best_algs.append(alg)
                    score = a_score
            # If multiple classifiers agree on the same max score (that isn't supporting), increment by one
            if len(best_algs) > 1 and best_score > 1:
                score += 1
            if len(alg_to_score) < 2 and score == 1:  # A single algorithm at supporting strength is not acceptable
                return 0, ""
            score = min(score, 4) # Do not weigh above very strong
        if score == 1:
            bp4 = f"BP4: {lines_of_evidence} line(s) of computational evidence support a benign effect; {formatted_printout_text}"
        else:
            if score == 2: # There is no moderate benign in default ACMG
                score = 3  # Make it strong
            bp4 = f"BP4_{score_to_hum_readable[score]}: {lines_of_evidence} line(s) of computational evidence support a benign effect; {formatted_printout_text}"
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

    ClinVar annotations are applied through NIRVANA, the most clinically relevant
    annotation is evaluated here.
    """
    score = 0
    bp6 = ""
    value = clinvar_review_status_to_level.get(variant.clinvar_review_status, 0)
    if value >= benign_thresholds['bp6_clinvar_minimum_review'] and "benign" in variant.clinvar_significance:
        score = value
        if score == 1:
            bp6 = f"BP6: Variant was found reported as {variant.clinvar_significance} "+\
                    f"in ClinVar as {variant.clinvar_id} with review status of {variant.clinvar_review_status}."
        else:
            bp6 = f"BP6_{score_to_hum_readable[score]}: Variant was found reported as {variant.clinvar_significance} "+\
                    f"in ClinVar as {variant.clinvar_id} with review status of {variant.clinvar_review_status}."
    return score, bp6


def get_bp7(variant, chrom_to_pos_to_alt_to_splice_score):
    """
    BP7    A synonymous (silent) variants for which splicing prediction algorithms
           predict no impact to the splice consensus sequence nor the creation of a
           new splice site AND the nucleotide is not highly conserved
   
    ABSplice - Wagner N, Çelik MH, Hölzlwimmer FR, Mertes C, Prokisch H, Yépez VA, Gagneur J. Aberrant splicing prediction across human tissues.

    Phylop - Pollard KS, Hubisz MJ, Rosenbloom KR, Siepel A. "Detection of nonneutral substitution rates on mammalian phylogenies."
        Genome Research. 2010;20(1):110-121. doi: 10.1101/gr.097857.109.

    Cutoff values were used from "Calibration of computational tools for missense variant pathogenicity
    classification and ClinGen recommendations for PP3/BP4 criteria" - table 2. When available. If unavailable, the
    cutoffs were taken from the base publication.
    https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/pdf/main.pdf
    

    evRepo assigns the following weights when applying BB7 depending on supporting evidence
    BP7_Strong (3), BP7 (1)
    """
    score = 0
    bp7 = ""

    # Conservation data (phylop) indicates the nucleotide is highly conserved
    if variant.phylopScore:
        if variant.phylopScore > benign_thresholds['bp7_phylop_cutoff']:
            return score, bp7

    if "splice" in variant.consequence:
        return score, bp7

    splice_score = None 
    # Consider bases that are predicted to have a strong impact on splicing
    # Authors believe .2 implies pathogenicity, .05 as a middle end, and below .01 as a benign indicator
    # https://github.com/gagneurlab/absplice?tab=readme-ov-file#output
    if chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome):
        if chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome).get(variant.position):
            if chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome).get(variant.position).get(variant.altAllele):
                splice_score = chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome).get(variant.position).get(variant.altAllele)
                if splice_score > benign_thresholds['bp7_splice_high']:
                    return score, bp7

    score = 1
    if splice_score:
        # Not conserved at all, nor associated with any significant splicing. 
        if splice_score < benign_thresholds['bp7_splice_low'] and variant.phylopScore <= benign_thresholds['bp7_phylop_supporting']:
            score = 3

    if "synonymous" in variant.consequence: 
        if score == 1:
            bp7 = f"BP7: Variant has synonymous associated consequence {variant.consequence} with ABSplice score {splice_score}"
        else:
            bp7 = f"BP7_{score_to_hum_readable[score]}: Variant has synonymous associated consequence {variant.consequence} with ABSplice score {splice_score}"
    elif "intron" in variant.consequence: 
        if score == 1:
            bp7 = f"BP7: Variant has intronic associated consequence {variant.consequence} with ABSplice score {splice_score}"
        else:
            bp7 = f"BP7_{score_to_hum_readable[score]}: Variant has intronic associated consequence {variant.consequence} with ABSplice score {splice_score}"
    elif variant.geneName == "intergenic":
        if score == 1:
            bp7 = f"BP7: Intergenic variant with ABSplice score {splice_score}"
        else:
            bp7 = f"BP7_{score_to_hum_readable[score]}: Intergenic variant with ABSplice score {splice_score}"

    else:
        return 0, ""

    return score, bp7


def get_bp(variant, truncating_gene_to_data, chrom_to_repeat_regions, chrom_to_pos_to_alt_to_splice_score, skip_list):
    """
    Supporting evidence of benign impact
    """
    
    if "bp1" not in skip_list:
        bp1_score, bp1 = get_bp1(variant, truncating_gene_to_data)
    else:
        bp1_score = 0
        bp1 = ""
    bp2 = get_bp2()
    if "bp3" not in skip_list:
        bp3_score, bp3 = get_bp3(variant, chrom_to_repeat_regions)
    else:
        bp3_score = 0
        bp3 = ""
    if "bp4" not in skip_list:
        bp4_score, bp4 = get_bp4(variant, chrom_to_pos_to_alt_to_splice_score)
    else:
        bp4_score = 0
        bp4 = ""
    bp5 = get_bp5()
    if "bp6" not in skip_list:
        bp6_score, bp6 = get_bp6(variant)
    else:
        bp6_score = 0
        bp6 = ""
    if "bp7" not in skip_list:
        bp7_score, bp7 = get_bp7(variant, chrom_to_pos_to_alt_to_splice_score)
    else:
        bp7_score = 0
        bp7 = ""

    bp_rationale_list = {'bp1': (bp1_score, bp1),
                         'bp2': (0, bp2),
                         'bp3': (bp3_score, bp3),
                         'bp4': (bp4_score, bp4),
                         'bp5': (0, bp5),
                         'bp6': (bp6_score, bp6),
                         'bp7': (bp7_score, bp7)
                         }
    return bp_rationale_list
