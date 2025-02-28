'''
Implementing the guidelines set forth in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/
ACMG/AMG 2015 interpretation guidelines. For each variant as many factors will be used as possible.
'''
from . import benign_classifiers, pathogenic_classifiers
from .constants import high_uncertain_threshold, low_uncertain_threshold 

def calc_pathogenic(pa, pvs, ps, pm, pp):
    """
    From table 5 of Richards et al 2015
    Pathogenic
        1. 1 Very Strong (PVS1) AND
                ≥1 Strong (PS1–PS4) OR
                ≥2 Moderate (PM1–PM6) OR
                1 Moderate (PM1–PM6) and 1 Supporting (PP1–PP5) OR
                ≥2 Supporting (PP1–PP5) 
        2. ≥2 Strong (PS1–PS4) OR
            1 Strong (PS1–PS4) AND
                ≥3 Moderate (PM1–PM6) OR
                2 Moderate (PM1–PM6) AND ≥2 Supporting (PP1–PP5) OR
                1 Moderate (PM1–PM6) AND ≥4 Supporting (PP1–PP5)

    Bitscopic Specific
        3. 1 or more Stand Alone (PA1)

        This weight is only applied when a variant is reported pathogenic or likely pathogenic 
        in a reliable database (PP5) with no conflicts and a review status of 'Practice Guidelines'
        or 'Reviewed by Expert Panel'. 
    """
    is_pathogenic = False
    if pa > 0: # One or more stand alone pathogenic - Bitscopic specific 
        is_pathogenic = True
    if pvs > 0: # At least one very strong pathogenic
        if ps > 0: # At least one strong pathogenic
            is_pathogenic = True
        if pm > 1: # More than one moderate pathogenic
            is_pathogenic = True
        if pm == 1 and pp > 0: # One moderate and one supporting
            is_pathogenic = True
        if pp > 1: # More than one supporting
            is_pathogenic = True
    if ps > 1: # Two or more strong criteria
        is_pathogenic = True
    if ps > 0: # One or more strong criteria
        if pm > 2: # At least three moderate
            is_pathogenic = True
        if pm > 1 and pp > 1: # At least two moderate and two supporting
            is_pathogenic = True
        if pm > 0 and pp > 3: # One moderate, four supporting
            is_pathogenic = True
    return is_pathogenic

def calc_likely_pathogenic(pvs, ps, pm, pp):
    """
    ACMG/AMP
    Likely Pathogenic
        1 Very Strong (PVS1) AND 1 Moderate (PM1–PM6) OR
        1 Strong (PS1–PS4) AND 1–2 Moderate (PM1–PM6) OR
        1 Strong (PS1–PS4) AND ≥2 Supporting (PP1–PP5) OR
        ≥3 Moderate (PM1–PM6) OR
        2 Moderate (PM1–PM6) AND ≥2 Supporting (PP1–PP5) OR
        1 Moderate (PM1–PM6) AND ≥4 Supporting (PP1–PP5)

    """
    is_likely_pathogenic = False
    if pvs > 0 and pm > 0: # At least one very strong and one moderate
        is_likely_pathogenic = True
    if ps > 0 and pm > 0: # At least one strong and one or more moderate
        is_likely_pathogenic = True
    if ps > 0 and pp > 1: # At least one strong and more than one supporting
        is_likely_pathogenic = True
    if pm > 2: # More than two moderate
        is_likely_pathogenic = True
    if pm > 1 and pp > 1: # More than one moderate and more than one supporting
        is_likely_pathogenic = True
    if pm > 0 and pp > 3: # At least one moderate and more than 3 supporting
        is_likely_pathogenic = True
    return is_likely_pathogenic

def calc_benign(ba, bvs, bs, bm, bp):
    """
    ACMG base publicatin did not consider benign very strong or benign moderate values for
    classifiers, however subsequent publications have applies these weights when evaluating
    the ACMG classifiers.  To ensure these weights can be considered in benign calculations
    we have adopted the logic from the pathogenic calculation, which already considers codes
    with weights of very strong and moderate.

    We note that when comparing the likely pathogenic and likely benign calculations that the
    pathogenic supporting codes of the same class are not evaluated the same as the supporting codes 
    in benign. This is by intent, as benign variants are more common and should not need as much
    supporting evidence to be labeled as such. To keep in line with this principle, we have lowered
    the threshold for moderate/supporting classification
    
    ACMG Benign definitions
        1 Stand-Alone (BA1) OR
        ≥2 Strong (BS1–BS4)
        
    Bitscopic Merged Definitions
        1. 1 Stand-Alone (BA1)
        2. 1 Very Strong (BVS1) AND
                ≥1 Strong (BS1–BS4) OR
                ≥2 Moderate () OR
                1 Moderate () and 1 Supporting (BP1–BP7) OR
                ≥2 Supporting (BP1–BP7) 
        3. ≥2 Strong (BS1–BS4) OR
            1 Strong (BS1–BS4) AND
                ≥2 Moderate () OR
                1 Moderate () AND ≥2 Supporting (BP1–BP7) OR
    """
    is_benign = False
    if ba > 0: # One stand alone
        is_benign = True
    if bvs > 0: # One very strong
        if bs > 0: # One strong
            is_benign = True
        if bm > 1: # Two moderate
            is_benign = True
        if bm > 0 and bp > 0: # One moderate and two supporting
            is_benign = True
        if bp > 1: # Two or more supporting
            is_benign = True
    if bs > 1: # More than one strong
        is_benign = True
    if bs > 0: # One strong
        if bm > 1: # Three moderate 
            is_benign = True
        if bm > 0 and bp > 1:
            is_benign = True

    return is_benign

def calc_likely_benign(bvs, bs, bm, bp):
    """
    ACMG base publicatin did not consider benign very strong or benign moderate values for
    AMP/ACMG Likely Benign
        1 Strong (BS1–BS4) and 1 Supporting (BP1–BP7) OR
        ≥2 Supporting (BP1–BP7)

    They did however, include values in the likely pathogenic 

    ACMG/AMP
    Likely Pathogenic
        1 Very Strong (PVS1) AND 1 Moderate (PM1–PM6) OR
        1 Strong (PS1–PS4) AND 1–2 Moderate (PM1–PM6) OR
        1 Strong (PS1–PS4) AND ≥2 Supporting (PP1–PP5) OR
        ≥3 Moderate (PM1–PM6) OR
        2 Moderate (PM1–PM6) AND ≥2 Supporting (PP1–PP5) OR
        1 Moderate (PM1–PM6) AND ≥4 Supporting (PP1–PP5)

    Where we can see they weigh supporting codes significantly less than in benign. This is likely
    by intent, as benign variants are more common and should not need as much supporting evidence
    to be labeled as such. 

    Bitscopic Likely Benign definitions
        1 Very-strong (BVS) OR
        1 Strong (BS1–BS4) and 1 Moderate (BM) OR
        1 Strong (BS1–BS4) and 1 Supporting (BP1-BP7) OR
        1 Moderate (BM) and >= 1 Supporting (BP1-BP7) OR
        >=2 Moderate (BM) OR
        >=2 Supporting (BP1-BP7)

    In the absense of BVS and BM scores, this will operate the same as the default ACMG combining
    criteria. BIAS does not assign BVS or BM scores currently, these weights can only be applied
    by the external call file. 
    """
    value = False
    if bvs > 0:
        value = True
    elif bs > 0:
        if bm > 0:
            value = True
        elif bp > 0:
            value = True
    elif bm > 0:
        if bp > 0:
            value = True
    elif bm > 1:
        value = True
    elif bp > 1:
        value = True
    return value


def calc_conflicting(pa, pvs, ps, pm, pp, ba, bvs, bs, bm, bp, threshold):
    """
    The 2015 ACMG paper states,
    
    'If a variant does not fulfill criteria using either of these sets (pathogenic or benign), or the evidence for
            benign and pathogenic is conflicting, the variant defaults to uncertain significance'
    """
    sa_weight = 10
    v_strong_weight = 7
    strong_weight = 5
    moderate_weight = 3
    supporting_weight = 1
    
    total_pathogenic = pa*sa_weight + pvs*v_strong_weight + ps*strong_weight + pm*moderate_weight + pp*supporting_weight
    total_benign = ba*sa_weight + bvs*v_strong_weight + bs*strong_weight + bm*moderate_weight + bp*supporting_weight

    if total_pathogenic >= threshold and total_benign >= threshold:
        return True
    return False



def merge_supplemental_codes_into_rationale_dict(rationale_dict, supplemental_codes):
    """
    Override the algorithmically determined classifiers with user provided codes and scores
    """
    merged_rationale_dict = {}
    for category, sup_code_to_rationale in supplemental_codes.items():
        code_to_rationale = rationale_dict[category]
        for code, sup_rationale in sup_code_to_rationale.items():
            if code_to_rationale[code] != sup_rationale:
                code_to_rationale[code] = sup_rationale
        merged_rationale_dict[category] = code_to_rationale
    return merged_rationale_dict


def recalculate_scores(rationale_dict):
    """
    Recalculate scores for each category after merging in user provided data
    """
    pvs_score, ps_score, pm_score, pp_score, ba_score, bs_score, bp_score = 0, 0, 0, 0, 0, 0, 0
    for category, code_list in rationale_dict.items():
        category_score = 0
        for code in code_list:
            category_score += code[0]
        if "pvs" in category:
            pvs_score = category_score
        elif "ps" in category:
            ps_score = category_score
        elif "pm" in category:
            pm_score = category_score
        elif "pp" in category:
            pp_score = category_score
        elif "ba" in category:
            ba_score = category_score
        elif "bs" in category:
            bs_score = category_score
        elif "bp" in category:
            bp_score = category_score

    return pvs_score, ps_score, pm_score, pp_score, ba_score, bs_score, bp_score
 

def get_weight_adjusted_scores(rationale_dict, b_or_p):
    """
    Each code has a weight that determines which class it should be applied too. This is somewhat confusing
    as the code's names imply a strength.  It is possible to have a PVS (pathogenic very strong) with a strength
    of supporting (1). Conversely, it is also possible to have a PP7 (pathogenic supporting #7)  code with 
    a strength of very strong (4).  

    Return the scores for each condition based on the weights applied to them
    """
    # These are lists as it was helpful for debugging. Technically we are only returning and using
    # the number of elements in these lists, so these could be integers.
    stand_alone = [] 
    very_strong = []
    strong = []
    moderate = []
    supporting = []
    for _, class_dict in rationale_dict.items():
        for code, evidence in class_dict.items():
            if code.startswith(b_or_p):
                code_score, _ = evidence # justification isn't needed
                if code_score >= 5: # Stand alone
                    stand_alone.append(code)
                if code_score >= 4: # Very strong
                    very_strong.append(code)
                elif code_score >= 3: # Strong
                    strong.append(code)
                elif code_score >= 2: # Moderate
                    moderate.append(code)
                elif code_score >= 1: # Supporting
                    supporting.append(code)
    return len(stand_alone), len(very_strong), len(strong), len(moderate), len(supporting)


def apply_ACMG_codes(variant, name_to_dataset, skip_list):
    """
    Apply all applicable ACMG codes to the variant and return them, along with detailed rationale
    as to why each code was applied, back in a rationale dictionary.
    """
    # Very strong pathogenic classifiers
    pvs_rationale = pathogenic_classifiers.get_pvs(variant, name_to_dataset['PVS1_gene_name_to_3prime_region'],
                                                   name_to_dataset['PVS1_PP3_BP4_BP7_splice_dict'], skip_list)

    # Strong pathogenic classifiers
    ps_rationale = pathogenic_classifiers.get_ps(variant, name_to_dataset['PS1_gene_mut_to_data'], name_to_dataset['PS3_lit_gene_mut_to_data'],
                                                 name_to_dataset['PS3_lit_variant_to_data'], name_to_dataset['PS4_chrom_to_pos_to_gwas_data'],
                                                 name_to_dataset['PVS1_PP3_BP4_BP7_splice_dict'], skip_list)
    
    # Moderate pathogenic classifiers
    pm_rationale = pathogenic_classifiers.get_pm(variant, name_to_dataset['PM4_BP3_chrom_to_repeat_regions'], name_to_dataset['PM5_gene_aa_to_var_data'],
                                                 name_to_dataset['PM1_chrom_to_pathogenic_domain_list'], name_to_dataset['PS1_gene_mut_to_data'], 
                                                 name_to_dataset['PVS1_PP3_BP4_BP7_splice_dict'], skip_list)

    # Supporting pathogenic classifiers
    pp_rationale = pathogenic_classifiers.get_pp(variant, name_to_dataset['PP2_missense_pathogenic_gene_to_region_list'],
                                                 name_to_dataset['PVS1_PP3_BP4_BP7_splice_dict'], skip_list)

    # Stand alone benign classifiers
    ba_rationale = benign_classifiers.get_ba(variant, skip_list)

    # Strong benign classifiers
    bs_rationale = benign_classifiers.get_bs(variant, skip_list)

    # Supporting benign classifiers
    bp_rationale = benign_classifiers.get_bp(variant, name_to_dataset['BP1_truncating_gene_to_data'], name_to_dataset['PM4_BP3_chrom_to_repeat_regions'], 
                                             name_to_dataset['PVS1_PP3_BP4_BP7_splice_dict'], skip_list)
    
    # Keep track of all the reasonings for each classifier
    rationale_dict = {
            "pvs": pvs_rationale,
            "ps": ps_rationale,
            "pm": pm_rationale,
            "pp": pp_rationale,
            "ba": ba_rationale,
            "bs": bs_rationale,
            "bp": bp_rationale
            }
    return rationale_dict


def get_variant_classification(variant, name_to_dataset, supplemental_codes, skip_list):
    """
    Calculate the pathogenicity evidence scores for the variants 
    then apply the classification system to determine pathogenicicty

    name_to_dataset maps strings (human readable) to python datastructs holding the defined data
    """
    # Apply the ACMG codes and store them in a rationale dictionary
    rationale_dict = apply_ACMG_codes(variant, name_to_dataset, skip_list)

    # The user has provided overriding input for some of the classifiers! Merge in their classifiers then
    # recalculate scores for each category
    if supplemental_codes:
        rationale_dict = merge_supplemental_codes_into_rationale_dict(rationale_dict, supplemental_codes)

    # Each code has a weight that determines which class it should be applied too. Calculate these to be
    # used with the ACMG combining criteria
    pa_score, pvs_score, ps_score, pm_score, pp_score = get_weight_adjusted_scores(rationale_dict, 'p')
    ba_score, bvs_score, bs_score, bm_score, bp_score = get_weight_adjusted_scores(rationale_dict, 'b')

    # Evaluate if the variant meets the ACMG combining criteria to be classified into a category. If the
    # variant cannot be categoriezed into Pathogenic, Likely Pathogenic, Benign or Likely Benign, then it
    # is reported as uncertain (equivalent to VUS)
    classification = 'uncertain'

    # Calculate if there the variant has conflicting classifiers with a high tolerance
    if calc_conflicting(pa_score, pvs_score, ps_score, pm_score, pp_score, ba_score, bvs_score, bs_score, bm_score, bp_score, high_uncertain_threshold):
        pass
    elif calc_pathogenic(pa_score, pvs_score, ps_score, pm_score, pp_score):
        classification = 'pathogenic'
    elif calc_benign(ba_score, bvs_score, bs_score, bm_score, bp_score):
        classification = 'benign'
    # Calculate if there the variant has conflicting classifiers with a low tolerance
    elif calc_conflicting(pa_score, pvs_score, ps_score, pm_score, pp_score, ba_score, bvs_score, bs_score, bm_score, bp_score, low_uncertain_threshold):
        pass
    elif calc_likely_pathogenic(pvs_score, ps_score, pm_score, pp_score):
        classification = 'likely pathogenic'
    elif calc_likely_benign(bvs_score, bs_score, bm_score, bp_score):
        classification = 'likely benign'
    return classification, rationale_dict
