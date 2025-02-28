'''
Implementing the guidelines set forth in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/
ACMG/AMG 2015 interpretation guidelines. For each variant as many factors will be used as possible.
'''

from . import benign_classifiers, pathogenic_classifiers

def calc_pathogenic(pvs, ps, pm, pp):
    """
    From table 5 of the publication
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
    """
    is_pathogenic = False
    if pvs > 0: # At least one very strong
        if ps > 0: # At least one strong
            is_pathogenic = True
        if pm > 1: # More than one moderate 
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

    Bitscopic
        ≥5 Supporting (PP1–PP5)
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
    if pp > 4: # Bitscopic analysis, greater than 4 supporting
        is_likely_pathogenic = True
    return is_likely_pathogenic

def calc_benign(ba, bs):
    """
    Benign
        1 Stand-Alone (BA1) OR
        ≥2 Strong (BS1–BS4)
    """
    if ba > 0: # One stand alone
        return True 
    if bs > 1: # More than one strong
        return True
    return False

def calc_likely_benign(bs, bp):
    """
    AMP/ACMG
    Likely Benign
        1 Strong (BS1–BS4) and 1 Supporting (BP1–BP7) OR
        ≥2 Supporting (BP1–BP7)
    """
    if bs > 0 and bp > 0:
        return True
    if bp > 1:
        return True 
    return False


def calc_conflicting(pvs_score, ps_score, pm_score, pp_score, ba_score, bs_score, bp_score):
    """
    The 2015 ACMG paper states,
    
    'If a variant does not fulfill criteria using either of these sets (pathogenic or benign), or the evidence for
            benign and pathogenic is conflicting, the variant defaults to uncertain significance'
    """
    total_pathogenic = pp_score + pm_score + ps_score + pvs_score
    total_benign = bp_score + bs_score + ba_score

    if total_pathogenic > 3 and total_benign > 3:
        return True
    return False


def merge_supplemental_codes_into_rationale_dict(rationale_dict, supplemental_codes):
    """
    Override the algorithmically determined classifiers with user provided codes and scores
    """
    merged_rationale_dict = {}
    for category, code_list in supplemental_codes.items():
        rationale_code_list = rationale_dict[category]
        index = 0
        new_code_list = []
        for code in code_list:
            if code:
                new_code_list.append(code)
            elif rationale_code_list[index]:
                new_code_list.append(rationale_code_list[index])
            index += 1
        merged_rationale_dict[category] = new_code_list
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


def get_variant_classification(variant, name_to_dataset, supplemental_codes):
    """
    Calculate the pathogenicity evidence scores for the variants 
    then apply the classification system to determine pathogenicicty

    name_to_dataset maps strings (human readable) to python datastructs holding the defined data
    """
    # Very strong pathogenic classifiers
    pvs_score, pvs_rationale = pathogenic_classifiers.get_pvs(variant, name_to_dataset['lof_gene_to_pli'], name_to_dataset['gene_name_to_3prime_region'])

    # Strong pathogenic classifiers
    ps_score, ps_rationale = pathogenic_classifiers.get_ps(variant, name_to_dataset['gene_mut_to_data'], name_to_dataset['chrom_to_pos_to_gwas_data'],
                                    name_to_dataset['gloc_to_pubmed_id_list'])
    
    # Moderate pathogenic classifiers
    pm_score, pm_rationale = pathogenic_classifiers.get_pm(variant, name_to_dataset['chrom_to_repeat_regions'], name_to_dataset['gene_aa_to_var_data'],
                                    name_to_dataset['pathogenic_domains'], name_to_dataset['gene_mut_to_data'])

    # Supporting pathogenic classifiers
    pp_score, pp_rationale = pathogenic_classifiers.get_pp(variant, name_to_dataset['missense_pathogenic_genes'])

    # Stand alone benign classifiers
    ba_score, ba_rationale = benign_classifiers.get_ba(variant)

    # Strong benign classifiers
    bs_score, bs_rationale = benign_classifiers.get_bs(variant)

    # Supporting benign classifiers
    bp_score, bp_rationale = benign_classifiers.get_bp(variant, name_to_dataset['truncating_genes'], name_to_dataset['chrom_to_repeat_regions'])

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

    # The user has provided overriding input for some of the classifiers! Merge in their classifiers then
    # recalculate scores for each category
    if supplemental_codes:
        rationale_dict = merge_supplemental_codes_into_rationale_dict(rationale_dict, supplemental_codes)
        pvs_score, ps_score, pm_score, pp_score, ba_score, bs_score, bp_score = recalculate_scores(rationale_dict)

    # Evaluate if the variant meets the vriteria to be clasified into any category
    if calc_pathogenic(pvs_score, ps_score, pm_score, pp_score):
        return 'pathogenic', rationale_dict
    if calc_likely_pathogenic(pvs_score, ps_score, pm_score, pp_score):
        return 'likely pathogenic', rationale_dict
    if calc_conflicting(pvs_score, ps_score, pm_score, pp_score, ba_score, bs_score, bp_score):
        return 'uncertain', rationale_dict
    if calc_benign(ba_score, bs_score):
        return 'benign', rationale_dict
    if calc_likely_benign(bs_score, bp_score):
        return 'likely benign', rationale_dict
    return "uncertain", rationale_dict
