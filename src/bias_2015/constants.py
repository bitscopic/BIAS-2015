"""
Constants for BIAS-2015. Users may change these at their own discretion. Bitscopic does not take responsibility
for any changes made by users. Re-validation is highly encouraged if any values are changed, the validation
code is provided in the scripts directory, see documentation for more details. 
"""

# The codes BIAS currently evaluates, please use the optional 'skip list' argument to ignore classifiers
# do not comment them out here. 
evaluated_codes = ['pvs1', 'ps1', 'ps3', 'ps4', 'pm1', 'pm2', 'pm4', 'pm5', 'pp2', 'pp3', 'pp5',
                   'ba1', 'bs1', 'bs2', 'bp1', 'bp3', 'bp4', 'bp6', 'bp7']

# Map weight scores to their human readable equivalents 
score_to_hum_readable = {
    5: 'stand-alone',
    4: 'very-strong',
    3: 'strong',
    2: 'moderate',
    1: 'supporting'
}

# Score thresholds for assigning uncertainty - if total path score AND total ben score are higher than the threshold
high_uncertain_threshold=10
low_uncertain_threshold=3

# Map ClinVar reviewStatus to an integer value
# Weight the score based on the ClinVar review status. This is an application of the established fact that
# expertly reviewed entries should be considered more strongly than entries with limited supported evidence.
# This concept and the concept of weighted scoring is explored in the following publications;

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
    "criteria provided, conflicting classifications": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classification provided": 0,
    "no classification for the single variant": 0,
    "no classifications from unflagged records": 0,
}

# Existing content remains unchanged, adding new constants below:

benign_thresholds = {
    # BA1, BS1, and BS2 have LOEUF dependent cutoff values, please see futher below

    # BP3 in-frame deletion/insertion length thresholds from ACMG guidelines and ClinGen expert panel standards
    "bp3_in_frame_max_length": 60,     # Maximum in-frame length for BP3 (19 amino acids)
    "bp3_in_frame_strong_length": 15,  # Strong BP3 threshold (5 amino acids or less)

    #Cutoff values were used from "Calibration of computational tools for missense variant pathogenicity
    #classification and ClinGen recommendations for PP3/BP4 criteria" - table 2
    # https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/pdf/main.pdf
    # BP4 thresholds for conservation and computational evidence, adapted from ClinGen's calibration paper
    "bp4_phylop_path_cutoff": 9.741,   # PnyloP score indicating a pathogenic signal - this disqualifies BP4
    "bp4_phylop_low": 1.879,           # PhyloP score indicating low conservation
    "bp4_phylop_very_low": 0.021,      # PhyloP very low conservation cutoff (high benignity support)
    "bp4_revel_path_cutoff": 0.932,    # Revel score indicating a pathogenic signal - this disqualifies BP4
    "bp4_revel_supporting": 0.290,     # REVEL supporting benign cutoff from ACMG calibration
    "bp4_revel_moderate": 0.183,       # REVEL moderate benign cutoff
    "bp4_revel_strong": 0.016,         # REVEL strong benign cutoff
    "bp4_revel_very_strong": 0.003,    # REVEL very strong benign cutoff
    "bp4_dann_supporting": 0.5,        # DANN supporting benign cutoff
    "bp4_dann_moderate": 0.35,         # DANN supporting benign cutoff
    "bp4_dann_strong": 0.2,            # DANN strong benign cutoff based on score calibration
    "bp4_gerp_supporting": 2.7,        # GERP supporting benign cutoff (low evolutionary constraint)
    "bp4_gerp_moderate": -4.54,        # GERP moderate benign cutoff from ClinGen recommendations
    "bp4_absplice_strong": .01,        # ABSplice cutoff for assigning a strong score

    # BP6 thresholds for using clinvar evidence
    "bp6_clinvar_minimum_review": 1,

    # BP7 thresholds for synonymous variants with no splicing impact, adapted from ABSplice authors' recommended cutoffs
    "bp7_phylop_cutoff": 7.367,        # PhyloP threshold for conserved regions
    "bp7_phylop_supporting": 0.21,     # PhyloP threshold for conserved regions
    "bp7_splice_high": 0.05,           # ABSplice high impact cutoff from Gagneur Lab
    "bp7_splice_low": 0.01             # ABSplice low impact benign cutoff from Gagneur Lab
}

pathogenic_thresholds = {
    # PVS1 thresholds for splice impact from ABSplice (Gagneur Lab)
    "pvs1_splice_moderate": 0.05,
    "pvs1_splice_strong": 0.2,

    # PS1 thresholds for ClinVar review weight
    # Source: ACMG guidelines - Higher review status from ClinVar indicates stronger evidence
    "ps1_strong_review": 4,         # Strong: well-reviewed ClinVar variants
    "ps1_moderate_review": 2,       # Moderate: limited review or lower confidence
    "ps1_supporting_review": 0,     # Supporting: minimal review
    "ps1_splice_moderate": 0.05,    # Splice variants with splice score higher than this are excluded from PS1 - as protein change is not a likely method of pathogenicity

    # PS3 thresholds for the number of supporting articles (AVADA)
    "ps3_supporting_articles": 0,   # Minimum articles for supporting
    "ps3_moderate_articles": 3,     # Minimum articles for moderate
    "ps3_strong_articles": 5,       # Minimum articles for strong

    # PS4 odds ratio (OR) and p-value thresholds from ACMG guidelines
    "ps4_or_strong": 3.0,
    "ps4_or_moderate": 2.0,
    "ps4_or_supporting": 1.5,
    "ps4_pvalue_strong": 0.05,
    "ps4_pvalue_moderate": 0.1,
    "ps4_pvalue_supporting": 0.2,
    # PS4 TOPMed fallback thresholds
    "ps4_topmed_ac_supporting_max": 4,  # Maximum allele count for supporting evidence
    "ps4_topmed_ac_supporting_min": 1,  # Maximum allele count for supporting evidence
    "ps4_topmed_hc_supporting_max": 9,  # Maximum homozygous count for supporting evidence
    "ps4_topmed_hc_supporting_min": 3,  # Maximum homozygous count for supporting evidence

    # PM1 thresholds for pathogenic domains
    # Domains with higher pathogenic-to-benign variant ratios provide stronger evidence of pathogenicity.
    "pm1_strong_ratio": 2.5,     # Strong: Pathogenic-to-benign ratio > 2.5
    "pm1_supporting_ratio": 1.5,  # Supporting: Pathogenic-to-benign ratio > 1.0
    "pm1_cutoff_ratio": 1.0,  # Variants which fall into regions with a ratio below this value will not have PM1 applied.

    # PM4 length thresholds for pathogenicity from Cannon et al. (2023)
    "pm4_strong_length": 30,
    "pm4_moderate_length": 10,
    "pm4_supporting_length": 0,
    "pm5_strong_review": 4,      # Strong: well-reviewed ClinVar variants
    "pm5_supporting_review": 1,  # Supporting: minimal review
    "pm5_review_cutoff": 0,      # Variants with a review status below this will not be used in PM5
    "pm5_splice_moderate": 0.05,


    # PP2 O/E thresholds and pathogenicity ratios from ClinVar and GnomAD RMC. A region can have ClinVar and/or Gnomad RMC annotation
    # all applicable cutoffs will be applied
    "pp2_oe_cutoff": 0.25,
    "pp2_benign_cutoff": 0.25,
    "pp2_path_cutoff": .5,
    "pp2_oe_very_low": 0.1,
    "pp2_pathogenic_high": 0.85,
    "pp2_benign_low": 0.15,

    # PP3 thresholds for computational predictions (PhyloP, REVEL, ABSplice)
    "pp3_phylop_moderate": 9.741,
    "pp3_phylop_supporting": 7.376,
    "pp3_revel_strong": 0.932,
    "pp3_revel_moderate": 0.773,
    "pp3_revel_supporting": 0.644,
    "pp3_absplice_strong": 0.2,
    
    # PP5 thresholds for using clinvar evidence
    "pp5_clinvar_minimum_review": 1,
}


# LOEUF (Loss-of-Function Observed/Expected Upper Bound Fraction) measures how tolerant a gene is to loss-of-function 
# variants, with lower values indicating stronger selection against LoF variants and higher values indicating greater tolerance.
#
# Highly constrained genes (low LOEUF) are more likely to have pathogenic variants → BA1 should be harder to trigger 
# (higher AF threshold) to avoid false benign calls, PM2 should be easier to trigger (higher AF threshold) to detect
# more variants.
#
# Unconstrained genes (high LOEUF) are less likely to have pathogenic variants → BA1 should be easier to trigger 
# (lower AF threshold) to avoid false benign calls, PM2 should be harder to trigger (low AF threshold) to detect
# less variants
loeuf_thresholds = [
    { # Very highly constrained genes
        "max_loeuf": 0.05, 
        "ba1_cutoff": 0.05, 
        "bs1_cutoff_strong": 0.02, 
        "bs1_cutoff": 0.01, 
        "pm2_cutoff": 0.001, 
        "pm2_cutoff_strong": 0.0005, 
        "bs2_recessive_homozygous_threshold": 8, 
        "bs2_dominant_allele_threshold": 8, 
        "bs2_xlinked_threshold": 8, 
        "bs2_af_threshold": 2e-5
    },
    { # Strongly constrained genes
        "max_loeuf": 0.2, 
        "ba1_cutoff": 0.01, 
        "bs1_cutoff_strong": 0.002, 
        "bs1_cutoff": 0.001, 
        "pm2_cutoff": 0.0005, 
        "pm2_cutoff_strong": 0.00025, 
        "bs2_recessive_homozygous_threshold": 6, 
        "bs2_dominant_allele_threshold": 6, 
        "bs2_xlinked_threshold": 6, 
        "bs2_af_threshold": 5e-5
    },
    { # Moderately constrained genes
        "max_loeuf": 0.5,
        "ba1_cutoff": 0.0025, 
        "bs1_cutoff_strong": 0.0002, 
        "bs1_cutoff": 0.0001, 
        "pm2_cutoff": 0.00001, 
        "pm2_cutoff_strong": 0.000001, 
        "bs2_recessive_homozygous_threshold": 5, 
        "bs2_dominant_allele_threshold": 5, 
        "bs2_xlinked_threshold": 5, 
        "bs2_af_threshold": 5e-4
    },
    { # Weakly constrained genes
        "max_loeuf": 1.0,
        "ba1_cutoff": 0.0025, 
        "bs1_cutoff_strong": 0.0004, 
        "bs1_cutoff": 0.0002, 
        "pm2_cutoff": 0.00001, 
        "pm2_cutoff_strong": 0.000001, 
        "bs2_recessive_homozygous_threshold": 3, 
        "bs2_xlinked_threshold": 3, 
        "bs2_af_threshold": 5e-4,
        "bs2_dominant_allele_threshold": 3
    },
    { # Tolerant genes (default)
        "max_loeuf": 100.0, 
        "ba1_cutoff": 0.0025, 
        "bs1_cutoff_strong": 0.0004, 
        "bs1_cutoff": 0.0002, 
        "pm2_cutoff": 0.00001, 
        "pm2_cutoff_strong": 0.000001, 
        "bs2_recessive_homozygous_threshold": 2, 
        "bs2_xlinked_threshold": 2, 
        "bs2_af_threshold": 5e-4, 
        "bs2_dominant_allele_threshold": 2
    }
]


# All ClinVar classifications. Each is mapped to one of three values; benign, uncertain and pathogenic
# the 'likely' values have been collapsed into their main category.  These are used to identify genomic
# regions with certain features, ex. regions with high rate of pathogenic missense variation. 
classification_mapping = {
    'Benign|confers_sensitivity': 'benign',
    'Conflicting_classifications_of_pathogenicity|protective': 'uncertain',
    'Benign': 'benign',
    'drug_response': 'uncertain',
    'Conflicting_classifications_of_pathogenicity|association': 'uncertain',
    'Likely_pathogenic|risk_factor': 'pathogenic',
    'Benign/Likely_benign|other': 'benign',
    'Likely_risk_allele': 'uncertain',
    'Affects|risk_factor': 'pathogenic',
    'Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance/Established_risk_allele': 'pathogenic',
    'Pathogenic/Likely_pathogenic/Likely_risk_allele': 'pathogenic',
    'not_provided': 'uncertain',
    'Conflicting_classifications_of_pathogenicity|other': 'uncertain',
    'Pathogenic/Pathogenic,_low_penetrance|other': 'pathogenic',
    'Likely_benign': 'benign',
    'Pathogenic/Likely_pathogenic|drug_response': 'pathogenic',
    'drug_response|risk_factor': 'uncertain',
    'Uncertain_significance': 'uncertain',
    'association': 'uncertain',
    'Likely_pathogenic|drug_response': 'pathogenic',
    'Affects': 'pathogenic',
    'Benign/Likely_benign|association': 'benign',
    'Benign|drug_response': 'benign',
    'Likely_benign|drug_response|other': 'benign',
    'drug_response|other': 'uncertain',
    'Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance|other': 'pathogenic',
    'Pathogenic|drug_response': 'pathogenic',
    'Pathogenic|risk_factor': 'pathogenic',
    'Conflicting_classifications_of_pathogenicity|risk_factor': 'uncertain',
    'Uncertain_risk_allele|risk_factor': 'uncertain',
    'Benign/Likely_benign': 'benign',
    'Benign|other': 'benign',
    'Pathogenic/Likely_pathogenic|other': 'pathogenic',
    'risk_factor': 'uncertain',
    'Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance': 'pathogenic',
    'Pathogenic|other': 'pathogenic',
    'Pathogenic|Affects': 'pathogenic',
    'no_classification_for_the_single_variant': 'uncertain',
    'Benign/Likely_benign|other|risk_factor': 'benign',
    'Pathogenic/Likely_risk_allele': 'pathogenic',
    'Conflicting_classifications_of_pathogenicity|association|risk_factor': 'uncertain',
    'Benign/Likely_benign|drug_response|other': 'benign',
    'Pathogenic/Likely_pathogenic|risk_factor': 'pathogenic',
    'protective': 'uncertain',
    'Likely_pathogenic': 'pathogenic',
    'Uncertain_significance|risk_factor': 'uncertain',
    'other': 'uncertain',
    'Likely_benign|other': 'benign',
    'Pathogenic/Pathogenic,_low_penetrance|other|risk_factor': 'pathogenic',
    'Benign/Likely_benign|drug_response': 'benign',
    'Likely_pathogenic/Likely_risk_allele': 'pathogenic',
    'Conflicting_classifications_of_pathogenicity|other|risk_factor': 'uncertain',
    'Uncertain_significance/Uncertain_risk_allele': 'uncertain',
    'Pathogenic': 'pathogenic',
    'Pathogenic/Likely_pathogenic': 'pathogenic',
    'Likely_pathogenic|other': 'pathogenic',
    'Uncertain_significance|drug_response': 'uncertain',
    'protective|risk_factor': 'uncertain',
    'Benign/Likely_benign|risk_factor': 'benign',
    'Conflicting_classifications_of_pathogenicity': 'uncertain',
    'Pathogenic/Likely_pathogenic|association': 'pathogenic',
    'Uncertain_significance|other': 'uncertain',
    'Benign|association': 'benign',
    'Likely_pathogenic|association': 'pathogenic',
    'Uncertain_risk_allele': 'uncertain',
    'Likely_benign|association': 'benign',
    'Conflicting_classifications_of_pathogenicity|drug_response': 'pathogenic'
}
