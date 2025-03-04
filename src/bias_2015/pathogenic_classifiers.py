"""
BIAS-2015 implementation of the ACMG 2015 germline pathogenic classifiers
"""
from .constants import clinvar_review_status_to_level, score_to_hum_readable, loeuf_thresholds, pathogenic_thresholds

def get_pvs(variant, gene_name_to_3prime_region, chrom_to_pos_to_alt_to_splice_score, skip_list):
    """
    Very strong evidence of pathogenicity
        PVS1    Null variants (nonsense, frameshift, canonical +/−1 or 2 splice sites, initiation
                codon, single or multi-exon deletion) in a gene where loss of function (LOF)
                is a known mechanism of disease

        Caveats:
            1. Beware of genes where LOF is not a known disease mechanism (e.g. GFAP, MYH7)
            2. Use caution interpreting LOF variants at the extreme 3’ end of a gene
            3. Use caution with splice variants that are predicted to lead to exon skipping but leave the remainder of the protein intact
            4. Use caution in the presence of multiple transcripts
    
    LOEUF values were identifies from GNOMAD, see preprocessing for details. 3' regions were identified
    using RefSeq, see preprocessing for details.

    evRepo assigns the following weights when applying PVS1 depending on supporting evidence
    PVS1 (4), PVS1_Strong (3),  PVS1_Moderate (2), PVS1_Supporting (1)
    """
    if 'pvs1' in skip_list:
        return {'pvs1': (0, "")}
    # Consider variants with a null consequence 
    null_variant_consequence_list = ['frameshift_variant',
                                     'stop_gained',
                                     'start_lost',
                                     'splice_donor_variant',
                                     'splice_acceptor_variant',
                                     'protein_altering_variant',
                                     ]
    valid_variant = False
    for null_cons in null_variant_consequence_list:
        if null_cons in variant.consequence:
            valid_variant = True
            break

    if not valid_variant:
        return {'pvs1': (0, f"{variant.consequence} is not a LoF consequence")}

    # Caveat 2: LOF variant is at the extreme 3' end of a gene (final 50 of last coding exon)
    null_region = ('chr0', -10, -1)
    prime_region = gene_name_to_3prime_region.get(variant.geneName, null_region)
    _, prime_start, prime_end = prime_region
    if prime_start <= int(variant.position) <= prime_end:
        return {'pvs1': (0, f"Position {variant.position} is in the extreme three prime region ({prime_start}-{prime_end}) of the last coding exon in gene {variant.geneName}")}

    # Caveat 3: Check splice variants that lead to exon skipping but leave the remainder of the protein intact
    # One must also be cautious in assuming that a null variant will lead to disease if found in an exon where no
    # other pathogenic variants have been described given the possibility that the exon may be alternatively spliced.
    # This is particularly true if the predicted truncating variant is identified as an incidental finding (unrelated
    # to the indication for testing) given the low prior likelihood of finding a pathogenic variant in that setting.

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
        return {'pvs1': (0, f"Out of {len(variant.transcriptList)} transcripts, variant consequence {variant.consequence} in " + \
                f"transcript {variant.transcriptList[0]['transcript']} was only seen once.")}
     
    # Consider bases that are predicted to have a strong impact on splicing
    # Authors suggest .2 as a high end cutoff https://github.com/gagneurlab/absplice?tab=readme-ov-file#output
    splice_score = None # 0 indicates no impact on splicing.
    if chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome):
        if chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome).get(variant.position):
            if chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome).get(variant.position).get(variant.altAllele):
                splice_score = chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome).get(variant.position).get(variant.altAllele)

    if 'splice' in variant.consequence:
        if splice_score:
            if splice_score < pathogenic_thresholds['pvs1_splice_moderate']:
                score = 2
            elif splice_score < pathogenic_thresholds['pvs1_splice_strong']:
                score = 3
            else:
                score = 4
        else: # Splice variant with no splice score, put it at moderate weight
            score = 2 
    else:
        score = 4
    loeuf = variant.gene_gnomad.get('loeuf', 1.0) if variant.gene_gnomad else 1.0
    if not variant.gene_gnomad:
        score = 4
        pvs1 = f"PVS1_{score_to_hum_readable[score]}: Null variant consequence {variant.consequence} in gene {variant.geneName}"
        if splice_score:
            pvs1 += f" with AbSplice impact {splice_score:.5f}."
        return {'pvs1': (score, pvs1)}

    if loeuf > 1.0:  # Gene not sensitive to mutations
        score = 3
        pvs1 = f"PVS1_{score_to_hum_readable[score]}: Null variant consequence {variant.consequence} in gene {variant.geneName} with gnomAD LOEUF={loeuf:.5f}."
        if splice_score:
            pvs1 += f" AbSplice impact={splice_score:.5f}."
        return {'pvs1': (score, pvs1)}

    # Default: No caveats apply
    score = 4
    pvs1 = f"PVS1: Null variant consequence {variant.consequence} in gene {variant.geneName} with gnomAD LOEUF={loeuf:.5f}."
    if splice_score:
        pvs1 += f" AbSplice impact={splice_score:.5f}."
    return {'pvs1': (score, pvs1)}

def check_valid_gene_mut_for_ps1(gene_mut_to_data, gene_mut, variant, chrom_to_pos_to_alt_to_splice_score):
    """
    Verify this gene mutation is applicable for PS1
    """
    if not gene_mut_to_data.get(gene_mut):
        return False
    _, _, _, clinvar_id = gene_mut_to_data[gene_mut][0]
    if clinvar_id in variant.clinvar_id: # Don't self cite.
        return False
    splice_without_score = False
    if "splice" in variant.consequence:
        if chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome):
            if chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome).get(variant.position):
                if chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome).get(variant.position).get(variant.altAllele):
                    splice_score = chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome).get(variant.position).get(variant.altAllele)
                    if splice_score: 
                        if splice_score > pathogenic_thresholds['ps1_splice_moderate']:
                            return False
                    else:
                        splice_without_score = True
                else:
                    splice_without_score = True
            else:
                splice_without_score = True
        else:
            splice_without_score = True
    if splice_without_score:
        return False
    return True

def get_ps1(variant, gene_name, aa_mut, gene_mut_to_data, chrom_to_pos_to_alt_to_splice_score):
    """
    PS1    Same amino acid change as a previously established pathogenic variants
           regardless of nucleotide change
           
           Example:    Val->Leu caused by either G>C or G>T in the same codon
           
           Caveat:    Beware of changes that impact splicing rather than at the
                      amino acid/protein level

    evRepo assigns the following weights when applying PS1 depending on supporting evidence
    PS1 (3), PS1_Moderate (2), PS1_Supporting (1)
    """
    gene_mut = (gene_name, aa_mut)
    score = 0
    ps1 = ""
    if "synonymous" in variant.consequence:
        return 0, ""

    if gene_mut_to_data.get(gene_mut):
        valid_entries = []
        for entry in gene_mut_to_data[gene_mut]:
            if check_valid_gene_mut_for_ps1(gene_mut_to_data, gene_mut, variant, chrom_to_pos_to_alt_to_splice_score):
                valid_entries.append(entry)

        if not valid_entries:
            return 0, ""

        # Find the best entry based on ClinVar review status
        best_entry = max(valid_entries, key=lambda x: clinvar_review_status_to_level.get(x[2], 0))
        _, criteria, signif, _ = best_entry

        review_value = clinvar_review_status_to_level.get(criteria, 0)
        score = 3
        if review_value < pathogenic_thresholds['ps1_strong_review']:  # Only consider well-defined and reviewed variants
            score = 2
        elif review_value < pathogenic_thresholds['ps1_moderate_review']:
            score = 1
        elif review_value < pathogenic_thresholds['ps1_supporting_review']:
            score = 1
        else:
            return 0, ""

        if score == 3:
            ps1 = f'PS1: Amino acid change {aa_mut} in gene {gene_name} ' \
                  + f'is associated with Clinvar {signif} variant {variant.clinvar_id} with review status {criteria}'
        else:
            ps1 = f'PS1_{score_to_hum_readable[score]}: Amino acid change {aa_mut} in gene {gene_name} ' \
                  + f'is associated with Clinvar {signif} variant {variant.clinvar_id} with review status {criteria}'
    return score, ps1


def get_ps2():
    """ 
    PS2    De novo (both maternity and paternity confirmed) in a patient with the
           disease and no family history
        
           Note: Confirmation of paternity only is insufficient. Egg donation, surrogate
                 motherhood, errors in embryo transfer, etc. can contribute to non-maternity

    This requires manual data injection, or pulling data from a patients EHR. This information
    can be provided to the algorithm via the external call file. It cannot be calculated internally
    at this point.
    """
    return ''


def get_ps3(variant, lit_gene_mut_to_data, lit_variant_to_data):
    """
    PS3    Well-established in vitro or in vivo functional studies supportive of a
           damaging effect on the gene or gene product
        
           Note: Functional studies that have been validated and shown to be
                 reproducible and robust in a clinical diagnostic laboratory setting are
                 considered the most well-established
    
    Literature information was derived from AVADA, see preprocessing for details.
    AVADA uses AI to find pathogenic varaints from publications and link them to genomic coordinates.
    
    evRepo assigns the following weights when applying PS3 depending on supporting evidence
    PS3 (3), PS3_Moderate (2), PS3_Supporting (1)
    """
    score = 0
    ps3 = ""

    # Identify pubmed articles that list this variant via nucleotide change
    lit_variant = (variant.chromosome, variant.position, variant.refAllele, variant.altAllele)
    variant_pubmed_set = lit_variant_to_data.get(lit_variant)
    
    # Identify pubmed articles that list this variant via protein change + gene name
    gene_mut = (variant.geneName, variant.protein_variant)
    gene_mut_pubmed_set = lit_gene_mut_to_data.get(gene_mut)
   
    if not variant_pubmed_set and not gene_mut_pubmed_set:
        return score, ps3

    
    # Identify the unique pubmeds for this variant
    if variant_pubmed_set and gene_mut_pubmed_set:
        unique_pubmed_ids = list(gene_mut_pubmed_set | variant_pubmed_set) 
    elif variant_pubmed_set:
        unique_pubmed_ids = sorted(list(variant_pubmed_set))
    else:
        unique_pubmed_ids = sorted(list(gene_mut_pubmed_set))

    if len(unique_pubmed_ids) > pathogenic_thresholds['ps3_supporting_articles']:
        score = 1
    if len(unique_pubmed_ids) > pathogenic_thresholds['ps3_moderate_articles']: # Variant was seen in at least three different papers
        score = 2
    if len(unique_pubmed_ids) > pathogenic_thresholds['ps3_strong_articles']: # Variant was seen in at least 5 papers - this is a very strong signal.
        score = 3
    if score == 0:
        return score, ps3

    pubmed_str = ", ".join(unique_pubmed_ids)
    if score == 3:
        ps3 = f"PS3: Variants pathogenicity is supported by multiple pubmed studies {pubmed_str} as established by AVADA"
    else:
        if len(unique_pubmed_ids) == 1:
            ps3 = f"PS3_{score_to_hum_readable[score]}: Variants pathogenicity is supported by a single pubmed study {pubmed_str} as established by AVADA"
        else:
            ps3 = f"PS3_{score_to_hum_readable[score]}: Variants pathogenicity is supported by multiple pubmed studies {pubmed_str} as established by AVADA"
    return score, ps3


def get_ps4(variant, chrom_to_pos_to_gwas_data):
    """
    PS4    The prevalence of the variants in affected individuals is significantly
           increased compared to the prevalence in controls

           Note 1: Relative risk (RR) or odds ratio (OR), as obtained from case-control
                   studies, is >5.0 and the confidence interval around the estimate of RR or OR
                   does not include 1.0. See manuscript for detailed guidance.
        
           Note 2: In instances of very rare variants where case-control studies may
                   not reach statistical significance, the prior observation of the variants in
                   multiple unrelated patients with the same phenotype, and its absence in
                   controls, may be used as moderate level of evidence.
    
    GWAS scores from the GWAS catalog hosted via UCSC Genome Browser were used. See preprocessing
    for more details.
    
    evRepo assigns the following weights when applying PS4 depending on supporting evidence
    PS4 (3), PS4_Moderate (2), PS4_Supporting (1)
    """
    score = 0
    ps4 = ""

    # Check GWAS data first
    if chrom_to_pos_to_gwas_data.get(variant.chromosome):
        if chrom_to_pos_to_gwas_data[variant.chromosome].get(int(variant.position)):
            or_value, ci_lower, ci_upper, pubmed_id, trait, p_value, _ = chrom_to_pos_to_gwas_data[variant.chromosome][int(variant.position)]
            if ci_lower > 1.0: # This is a firm requirement set by ACMG, confidence interval must not overlap 1. 
                if or_value > pathogenic_thresholds['ps4_or_strong'] and p_value < pathogenic_thresholds['ps4_pvalue_strong']:
                    score = 3
                elif or_value > pathogenic_thresholds['ps4_or_moderate'] and p_value < pathogenic_thresholds['ps4_pvalue_moderate']:
                    score = 2
                elif or_value > pathogenic_thresholds['ps4_or_supporting'] and p_value < pathogenic_thresholds['ps4_pvalue_supporting']:
                    score = 1
            if score > 0:
                prefix = "PS4" if score == 3 else f"PS4_{score_to_hum_readable[score]}"
                ps4 = (
                    f"{prefix}: Genomic location {variant.chromosome}:{variant.position} associated with {trait}. "
                    f"OR={or_value} (95% CI: {ci_lower}-{ci_upper}), p={p_value}, PubMed={pubmed_id}"
                )
        # If no GWAS data, fallback to TOPMed
        elif hasattr(variant, 'topmed'):
            topmed_af = variant.topmed.get('allAf', None)
            if topmed_af is not None:
                if variant.topmed.get('failedFilter', False):
                    return 0, ""
                if (pathogenic_thresholds['ps4_topmed_ac_supporting_min'] <= variant.topmed.get('allAc', 0) <= pathogenic_thresholds['ps4_topmed_ac_supporting_max']) \
                        or (pathogenic_thresholds['ps4_topmed_hc_supporting_min'] <= variant.topmed.get('allHc') <= pathogenic_thresholds['ps4_topmed_hc_supporting_max']):
                    score = 1
                    prefix = f"PS4_{score_to_hum_readable[score]}"
                    ps4 = (
                        f"{prefix}: No GWAS data found. "
                        f"TOPMed allele count {variant.topmed.get('allAc', 0)} and homozygous count {variant.topmed.get('allHc', 0)} "
                        f"indicate rarity consistent with pathogenicity."
                    )
    return score, ps4


def get_ps(variant, gene_mut_to_data, lit_gene_mut_to_data, lit_variant_to_data, chrom_to_pos_to_gwas_data, chrom_to_pos_to_alt_to_splice_score, skip_list):
    """
    Look for strong evidence of pathogenicity, apply the PS1, PS3, and PS4 classifiers.

    PS2 requires familial data and is not calculated currently. Values can be provided for it through the
    external call file.
    """

    if "ps1" not in skip_list:
        score_1, ps1 = get_ps1(variant, variant.geneName, variant.protein_variant, gene_mut_to_data, chrom_to_pos_to_alt_to_splice_score)
    else:
        score_1 = 0
        ps1 = ""

    ps2 = get_ps2()

    if "ps3" not in skip_list:
        score_3, ps3 = get_ps3(variant, lit_gene_mut_to_data, lit_variant_to_data)
    else:
        score_3 = 0
        ps3 = ""

    if "ps4" not in skip_list:
        score_4, ps4 = get_ps4(variant, chrom_to_pos_to_gwas_data)
    else:
        score_4 = 0
        ps4 = ""

    ps_rationale_list = {
        'ps1': (score_1, ps1),
        'ps2': (0, ps2),
        'ps3': (score_3, ps3),
        'ps4': (score_4, ps4)
    }

    return ps_rationale_list


def get_pm1(variant, chrom_to_pathogenic_domain):
    """    
    PM1    Located in a mutational hot spot and/or critical and well-established
           functional domain (e.g. active site of an enzyme) without benign variation

    UniProt domains were annotated with ClinVar data to identify domains with a high rate (>50%)
    of pathogenic variation and a low rate of benign variation (<50%)
    
    Landrum, M.J., Lee, J.M., Benson, M., et al. (2018). ClinVar: improving access to variant interpretations and
    supporting evidence. Nucleic Acids Research, 46(D1), D1062–D1067.

    UniProt Consortium. Reorganizing the protein space at the Universal Protein Resource (UniProt). Nucleic Acids Res. 
    2012 Jan;40(Database issue):D71-5. PMID: 22102590; PMC: PMC3245120
    
    evRepo assigns the following weights when applying PM1 depending on supporting evidence
    PM1_Strong (3), PM1 (2), PM1_Supporting (1)
    """
    score = 0
    pm1 = ""
    v_pos = int(variant.position)
    pathogenic_domain_list = chrom_to_pathogenic_domain.get(variant.chromosome)
    best_score = 0
    rationale = ""
    if pathogenic_domain_list:
        for domain in pathogenic_domain_list:
            d_start, d_end, d_name, d_path_percent, _, d_ben_percent, _, _, source, path_ben_ratio = domain
            if d_start <= v_pos <= d_end:
                # The variant score indicates how many pathogenic variants were found in this domain and their relative authority.
                # In domains with an overwhelming amount of evidence supporting pathogenicity, adjust the score accordingly
                if path_ben_ratio < pathogenic_thresholds['pm1_cutoff_ratio']:
                    score = 0, ""
                    continue
                variant.domain = d_name
                score = 2
                if path_ben_ratio < pathogenic_thresholds['pm1_supporting_ratio']: # Not that pathogenic of a domain
                    score = 1
                if path_ben_ratio > pathogenic_thresholds['pm1_strong_ratio']:
                    score = 3
                if score == 2:
                    pm1_code = "PM1"
                else:
                    pm1_code = f"PM1_{score_to_hum_readable[score]}"
                    
                pm1 = f"{pm1_code}: Variant is in {source} domain, {d_name} with observed ClinVar pathogenic rate of " +\
                            f"{100*d_path_percent:.5f}% and benign rate of {100*d_ben_percent:.5f}."
                if score > best_score:
                    best_score = score
                    rationale = pm1
              
    return best_score, rationale


def get_pm2(variant):
    """
    PM2: Absent from controls (or at extremely low frequency if recessive)
         in Exome Sequencing Project, 1000 Genomes, or GnomAD.
    
    Uses LOEUF to determine appropriate AF thresholds. See constants file for more documentation.

    PM2_Strong (3), PM2 (2), PM2_Supporting (1)
    """
    pm2 = ""
    score = 0

    # Retrieve allele frequencies safely
    gnomad_af = variant.gnomad.get('allAf', 0) if variant.gnomad else None
    onekg_af = variant.oneKg.get('allAf', 0) if variant.oneKg else None

    # Retrieve LOEUF safely, defaulting to 1.0 if missing
    loeuf = variant.gene_gnomad.get('loeuf', 1.0) if variant.gene_gnomad else 1.0

    # Determine PM2 cutoff based on LOEUF
    pop_thresholds = {}
    for loeuf_tier in loeuf_thresholds:
        if loeuf < loeuf_tier['max_loeuf']:
            pop_thresholds = loeuf_tier
            break
    pm2_cutoff_supporting = pop_thresholds['pm2_cutoff']
    pm2_cutoff_strong = pop_thresholds['pm2_cutoff_strong']
    
    # 1. **Variant is entirely missing in both GnomAD and 1000 Genomes**
    if not variant.gnomad and not variant.oneKg:
        score = 1
        pm2 = "PM2: Variant is absent from both GnomAD and 1000 Genomes."
        return score, pm2

    # 2. **Variant is found in both but at extremely low AF**
    loeuf_str = f"({variant.geneName}:{loeuf})"
    if variant.gnomad and variant.oneKg:
        max_af = max(gnomad_af, onekg_af)
        if max_af < pm2_cutoff_supporting:
            score = 2
            pm2 = f"PM2_{score_to_hum_readable[score]}: Both GnomAD {gnomad_af*100:.5f}% and 1000 Genomes {onekg_af*100:.5f}% are below " + \
                    f"LOEUF{loeuf_str}-based threshold {pm2_cutoff_supporting*100:.5f}%."
        if max_af < pm2_cutoff_strong: 
            score = 3
            pm2 = f"PM2: Both GnomAD {gnomad_af*100:.5f}% and 1000 Genomes {onekg_af*100:.5f}% are below " + \
                    f"LOEUF{loeuf_str}-based threshold {pm2_cutoff_strong*100:.5f}%."
    # 3. **Variant is found in 1k**
    elif variant.oneKg:
        if onekg_af < pm2_cutoff_supporting:
            score = 2
            pm2 = f"PM2_{score_to_hum_readable[score]}: 1000 Genomes {onekg_af*100:.5f}% is below LOEUF{loeuf_str}-based threshold {pm2_cutoff_supporting*100:.5f}%."
        if onekg_af < pm2_cutoff_strong: 
            score = 3
            pm2 = f"PM2: 1000 Genomes {onekg_af*100:.5f}% is below LOEUF{loeuf_str}-based threshold {pm2_cutoff_strong*100:.5f}%."
    # 4. **Variant is found in Gnomad**
    else:
        if gnomad_af < pm2_cutoff_supporting:
            score = 2
            pm2 = f"PM2_{score_to_hum_readable[score]}: gnomAD {gnomad_af*100:.5f}% is below LOEUF{loeuf_str}-based threshold {pm2_cutoff_supporting*100:.5f}%."
        if gnomad_af < pm2_cutoff_strong: 
            score = 3
            pm2 = f"PM2: gnomAD {gnomad_af*100:.5f}% is below LOEUF{loeuf_str}-based threshold {pm2_cutoff_strong*100:.5f}%."
    return score, pm2


def get_pm3():
    """
    Trans variant mapping, only applies when two variants are at the same position and both have ~.5 af, indicating
    that the patient has one variant on each allele. 
    """
    # NOTE: eRepo uses this flag quite a bit - However we don't see how a single variant can be 'trans' in orientation unless
    # its being evaluated against another variant. This is a bit out of scope for the algorithm as its currently
    # implemented, as it takes a single variant in at a time to classify it.  We could add neighboring variant
    # information in, however we would need to understand the context. 
    return ""


#NOTE: This is duplicated - find a better place for it
def get_variant_indel_length(ref, alt):
    """
    Compute the true length of an indel by trimming common prefixes and suffixes.
    Handles empty strings, complex indel sequences, and ensures correct length calculation.
    """
    # Edge case: if both are identical
    if ref == alt:
        return 0

    # Trim common prefix
    i = 0
    while i < len(ref) and i < len(alt) and ref[i] == alt[i]:
        i += 1
    ref = ref[i:]
    alt = alt[i:]

    # Trim common suffix
    j = 0
    while j < len(ref) and j < len(alt) and ref[-(j+1)] == alt[-(j+1)]:
        j += 1
    if j > 0:
        ref = ref[:-j]
        alt = alt[:-j]

    # True length is the maximum length of remaining sequences
    return max(len(ref), len(alt))

def get_pm4(variant, chrom_to_repeat_regions):
    """
    PM4    Protein length changes due to in-frame deletions/insertions in a non-repeat
            region or stop-loss variants 

    The deletion or insertion of one or more amino acids as well as the extension of a protein by changing the stop
    codon to an amino acid codon (e.g. a stop loss variant) is more likely to disrupt protein function as compared 
    to a missense change alone due to length changes in the protein. Therefore, in-frame deletions/insertions and
    stop lsses are considered moderate evidence of pathogenicity. The larger the deletion, insertion or extension,
    and the more conserved the amino acids are in a deleted region, then the more substantial is the evidence to
    support pathogenicity. In contrast, small in-frame deletions/insertions in repetitive regions, or regions that
    are not well conserved in evolution, are less likely to be pathogenic.

    Longer inframe insertion and deletions are more likely to be pathogenic
    Cannon, S., Williams, M., Gunning, A. C., & Wright, C. F. (2023). Evaluation of in silico pathogenicity prediction
        tools for the classification of small in-frame indels. BMC Medical Genomics, 16, Article 36.
        https://doi.org/10.1186/s12920-023-01454-6
    
    evRepo assigns the following weights when applying PM4 depending on supporting evidence
    PM4_Strong (3), PM4 (2), PM4_Supporting(1)
    """
    # Conservation regions be overlayed with the repeat regions to find highly conserved
    # non repeat regions. In frame insertion/deletions in these regions would be more likely to be pathogenic
    pm4 = ""
    score = 0
    pm4_consequences = ['stop_lost', 'inframe_insertion', 'inframe_deletion']
    if not variant.consequence in pm4_consequences:
        return score, pm4

    # Check if the variant is in a repeat region
    if chrom_to_repeat_regions.get(variant.chromosome):
        repeat_region_list = chrom_to_repeat_regions[variant.chromosome]
        pos = int(variant.position)
        length = get_variant_indel_length(variant.refAllele, variant.altAllele)
        for repeat_region in repeat_region_list:
            start, stop, _ = repeat_region
            if start <= pos <= stop or start <= (pos + length - 1) <= stop:
                return score, pm4

    # For non repeat regoin variants, update the score based on the INDEL length. Longer indels are more pathogenic
    if length > pathogenic_thresholds['pm4_strong_length']:
        score = 3
    elif length > pathogenic_thresholds['pm4_moderate_length']:
        score = 2
    elif length > pathogenic_thresholds['pm4_supporting_length']:
        score = 1
    else:
        score = 0
    
    if score == 2:
        pm4 = f"PM4: {variant.consequence} in a non-repeat region"
        if "insertion" in variant.consequence:
            pm4 = f"PM4: {variant.consequence} of {variant.altAllele} in a non-repeat region"
        if "deletion" in variant.consequence:
            pm4 = f"PM4: {variant.consequence} of {variant.refAllele} in a non-repeat region"
    else:
        pm4 = f"PM4_{score_to_hum_readable[score]}: {variant.consequence} in a non-repeat region"
        if "insertion" in variant.consequence:
            pm4 = f"PM4_{score_to_hum_readable[score]}: {variant.consequence} of {variant.altAllele} in a non-repeat region"
        if "deletion" in variant.consequence:
            pm4 = f"PM4_{score_to_hum_readable[score]}: {variant.consequence} of {variant.refAllele} in a non-repeat region"
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


def get_pm5(variant, gene_aa_to_var_data, gene_mut_to_data, chrom_to_pos_to_alt_to_splice_score):
    """
    PM5     Novel missense change at an amino acid residue where a different
            missense change determined to be pathogenic has been seen before
            
            Example: Arg156His is pathogenic; now you observe Arg156Cys

            Caveat: Beware of changes that impact splicing rather than at the amino
                    acid/protein level
    
    evRepo assigns the following weights when applying PM5 depending on supporting evidence
    PM5_Strong (3), PM5 (2), PM5_Supporting(1)
    """
    pm5 = ""
    score = 0
    # Variant must be in a gene to have an AA change associated
    if not variant.geneName:
        return score, pm5

    # If the full variant change was seen then PS1 was applied and PM5 cannot be applied
    gene_mut = (variant.geneName, variant.protein_variant)
    if check_valid_gene_mut_for_ps1(gene_mut_to_data, gene_mut, variant, chrom_to_pos_to_alt_to_splice_score):
        return 0, ""

    if "missense" not in variant.consequence:
        return score, pm5

    # Check if the base aa and position are associated with a known pathogenic variant
    variant_aa_comp = get_aa_comparator(variant.protein_variant)
    gene_aa = (variant.geneName, variant_aa_comp)
    if gene_aa_to_var_data.get(gene_aa):
        aa_data_list = gene_aa_to_var_data[gene_aa]
        ref_mut, rs_id, criteria, signif, clinvar_id = aa_data_list[0]
        not_valid = False
        if clinvar_id in variant.clinvar_id: # Don't self cite.
            not_valid = True
        if chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome):
            if chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome).get(variant.position):
                if chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome).get(variant.position).get(variant.altAllele):
                    splice_score = chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome).get(variant.position).get(variant.altAllele)
                    if splice_score and "splice" in variant.consequence:
                        if splice_score > pathogenic_thresholds['pm5_splice_moderate']:
                            not_valid = True
        if not_valid:
            return 0, ""
        # Find the data with the best score based on ClinVar review status
        best_data = max(aa_data_list, key=lambda x: clinvar_review_status_to_level[x[2]]) 

        ref_mut, rs_id, criteria, signif, clinvar_id = best_data
        score = 2
        if clinvar_review_status_to_level[criteria] < pathogenic_thresholds['pm5_review_cutoff']:
            return 0, ""
        if clinvar_review_status_to_level[criteria] < pathogenic_thresholds['pm5_supporting_review']:
            score = 1
        elif clinvar_review_status_to_level[criteria] > pathogenic_thresholds['pm5_strong_review']:  # Score variants with high review levels higher
            score = 3
        if score == 2:
            pm5_code = "PM5"
        else:
            pm5_code = f"PM5_{score_to_hum_readable[score]}"

        pm5 = f"{pm5_code}: Variant causes amino acid change {variant.protein_variant} which is " \
              + f"the same base amino acid as {ref_mut}, a {signif} variant with ClinVar ID {clinvar_id}, rs ID {rs_id}, and review status {criteria}."
        return score, pm5
    return score, pm5


def get_pm6():
    """
    Implied De Novo variant - difficult without phenotype and/or parental genotype/phenotype
    information
    """
    return ""


def get_pm(variant, chrom_to_repeat_regions, gene_aa_to_var_data, chrom_to_pathogenic_domain, gene_mut_to_data,
           chrom_to_pos_to_alt_to_splice_score, skip_list):
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

    if "pm1" not in skip_list:
        score_1, pm1 = get_pm1(variant, chrom_to_pathogenic_domain)
    else:
        score_1 = 0
        pm1 = ""

    if "pm2" not in skip_list:
        score_2, pm2 = get_pm2(variant)
    else:
        score_2 = 0
        pm2 = ""

    pm3 = get_pm3()

    if "pm4" not in skip_list:
        score_4, pm4 = get_pm4(variant, chrom_to_repeat_regions)
    else:
        score_4 = 0
        pm4 = ""

    if "pm5" not in skip_list:
        score_5, pm5 = get_pm5(variant, gene_aa_to_var_data, gene_mut_to_data, chrom_to_pos_to_alt_to_splice_score)
    else:
        score_5 = 0
        pm5 = ""

    pm6 = get_pm6()

    pm_rationale_list = {
        'pm1': (score_1, pm1),
        'pm2': (score_2, pm2),
        'pm3': (0, pm3),
        'pm4': (score_4, pm4),
        'pm5': (score_5, pm5),
        'pm6': (0, pm6)
    }
    
    return pm_rationale_list


def get_pp1():
    """
    Future goal: This will require patient EHR integration to automate as it requires patient family genetic information

    PP1    Co-segregation with disease in multiple affected family members in a gene
           definitively known to cause the disease
           Note: May be used as stronger evidence with increasing segregation data
    """
    return ""


def get_pp2(variant, missense_pathogenic_gene_to_region_list):
    """
    PP2    Missense variants in a gene that has a low rate of benign missense variation
           and where missense variants are a common mechanism of disease

    We consider two metrics for determining candidate reigons. GNOMAD RMC values, and
    regions identified in ClinVar as having high rates of missense pathogenicity and low
    rates of missense benign variants. 


    evRepo assigns the following weights when applying PP2 depending on supporting evidence
    PP2_Moderate (2), PP2 (1)
    """
    pp2 = ""
    score = 0
    if "missense" not in variant.consequence: # Only consider missense variants 
        return score, pp2

    if not variant.geneName: # Variant must be associated with a gene
        return score, pp2

    if not missense_pathogenic_gene_to_region_list.get(variant.geneName): # Gene is not associated with missense pathogenicity
        return score, pp2

    region_list = missense_pathogenic_gene_to_region_list[variant.geneName]
    full_gene = ""
    found_region = ""
    for region in region_list:
        start, end, oe, path_per, ben_per = region
        if start == 0 and end == 0: # Mapping was done for the full gene, this is a fallback if no precise gnomad RMC is found
            full_gene = region
            continue
        if start <= int(variant.position) <= end: # Variant is in a gnomad RMC region
            found_region = region
    if not found_region:
        found_region = full_gene
  
    # Variant was not in a region associated with missense pathogenicity 
    if not found_region:
        return score, pp2

    start, end, oe, path_per, ben_per = found_region
    failed = False
    if oe:
        if 1 > oe > pathogenic_thresholds['pp2_oe_cutoff']:
            failed = True
    if path_per < pathogenic_thresholds['pp2_path_cutoff'] or ben_per > pathogenic_thresholds['pp2_benign_cutoff'] or ben_per == 0 or path_per == 0:
        failed = True
    if failed:    
        return 0, ""

    score = 1
    if oe:
        if oe < pathogenic_thresholds['pp2_benign_low']: # A very low O/E ratio suggests this region is particularly intolerant to missense mutation
            score += 1
    if path_per and ben_per:
        if path_per > pathogenic_thresholds['pp2_pathogenic_high'] and ben_per < pathogenic_thresholds['pp2_benign_low']: # Highly pathogenic and rarely benign gene
            score += 1
    if score == 1:
        pp2_code = "PP2"
    else:
        pp2_code = f"PP2_{score_to_hum_readable[score]}"

    if path_per > 0 and oe: # RMC and region overlap
        pp2 = f"{pp2_code}: Missense variant consequence {variant.consequence} in {variant.geneName} where {100*path_per:.5f}% of pathogenic " + \
                f"variants are missense, {ben_per*100:.5f}% of missense variants are benign, and the region {start}-{end} has GNOMAD RMC OE value {oe:.5f}"
    elif path_per > 0: # Full gene overlap
        pp2 = f"{pp2_code}: Missense variant consequence {variant.consequence} in {variant.geneName} where {100*path_per:.5f}% of pathogenic " + \
                f"variants are missense, {ben_per*100:.5f}% of missense variants are benign"
    elif oe: # RMC overlap
        pp2 = f"{pp2_code}: Missense variant consequence {variant.consequence} in {variant.geneName} where " + \
                f"the region {start}-{end} has GNOMAD RMC OE value {oe:.5f}"

    return score, pp2

def get_pp3(variant, chrom_to_pos_to_alt_to_splice_score):
    """
    PP3    Multiple lines of computational evidence support a deleterious effect on
            the gene or gene product (conservation, evolutionary, splicing impact, etc)
            Caveat: As many in silico algorithms use the same or very similar input for
            their predictions, each algorithm should not be counted as an independent
            criterion. PP3 can be used only once in any evaluation of a variant.

    PhyloP, GERP, REVEL and DANN data is annotated with ICA/NIRVANA, the raw data sources were published under these citations

    Davydov, E. V., et al. (2010). Identifying a high fraction of the human genome to be under selective constraint
        using GERP++. PLoS Computational Biology, 6(12), e1001025.

    Quang, D., Chen, Y., & Xie, X. (2015). DANN: a deep learning approach for annotating the pathogenicity of genetic
        variants. Bioinformatics, 31(5), 761-763.

    Ioannidis, N. M., et al. (2016). REVEL: An ensemble method for predicting the pathogenicity of rare missense
        variants. American Journal of Human Genetics, 99(4), 877-885.
    
    Cutoff values were used from "Calibration of computational tools for missense variant pathogenicity
    classification and ClinGen recommendations for PP3/BP4 criteria" - table 2
    https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/pdf/main.pdf

    ABSplice cutoff values were derived from UCSC and the authors input on github 
    Wagner N, Çelik MH, Hölzlwimmer FR, Mertes C, Prokisch H, Yépez VA, Gagneur J. Aberrant splicing prediction across 
        human tissues. Nat Genet. 2023 May;55(5):861-870. PMID: 37142848

    evRepo assigns the following weights when applying PP3 depending on supporting evidence
    PP3_Strong (3), PP3_Moderate (2), PP3 (1)
    """
    pp3 = ""
    score = 0
    printout_text = []
    alg_to_score = {}

    if variant.phylopScore:
        if variant.phylopScore >= pathogenic_thresholds['pp3_phylop_supporting']:
            phylop_weight = "supporting"
            if variant.phylopScore >= pathogenic_thresholds['pp3_phylop_moderate']:
                phylop_weight = "moderate"
                alg_to_score['phylop'] = 2
            else:
                alg_to_score['phylop'] = 1
            printout_text.append(f"{phylop_weight} phylop {variant.phylopScore}")

    if variant.revel:
        if variant.revel >= pathogenic_thresholds['pp3_revel_supporting']:
            revel_weight = "supporting"
            if variant.revel >= pathogenic_thresholds['pp3_revel_strong']: # Strong pathogenic 
                alg_to_score['revel'] = 3
                revel_weight = "strong"
            elif variant.revel >= pathogenic_thresholds['pp3_revel_moderate']: # moderate
                alg_to_score['revel'] = 2
                revel_weight = "moderate"
            else: # supporting
                alg_to_score['revel'] = 1
            printout_text.append(f"{revel_weight} revel {variant.revel}")

    if 'splice' in variant.consequence: 
        # Also consider bases that are predicted to have a strong impact on splicing
        # Authors suggest .2 as a high end cutoff, .05 as a middle end, and .01 as a low https://github.com/gagneurlab/absplice?tab=readme-ov-file#output
        splice_score = 0 # 0 indicates no impact on splicing.
        if chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome):
            if chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome).get(variant.position):
                if chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome).get(variant.position).get(variant.altAllele):
                    splice_score = chrom_to_pos_to_alt_to_splice_score.get(variant.chromosome).get(variant.position).get(variant.altAllele)
        if splice_score >= pathogenic_thresholds['pp3_absplice_strong']:
            absplice_weight = "strong"
            alg_to_score['absplice'] = 3
            printout_text.append(f"{absplice_weight} ABSplice {splice_score}")
    
    if alg_to_score:
        formatted_printout_text = " | ".join(printout_text)
        lines_of_evidence = len(printout_text)
        score = 0
        best_score = 0
        best_algs = []
        # Identify the best score and its algorithm(s)
        for alg, a_score in alg_to_score.items():
            if a_score >= best_score:
                best_algs.append(alg)
                score = a_score
        # If multiple classifiers agree on the same max score, increment by one
        if len(best_algs) > 1:
            score += 1
        score = min(score, 3) # Do not weigh above strong
        if score == 1:
            pp3 = f"PP3: {lines_of_evidence} line(s) of computational evidence support a pathogenic effect; {formatted_printout_text}"
        else:
            pp3 = f"PP3_{score_to_hum_readable[score]}: {lines_of_evidence} line(s) of computational evidence support a pathogenic effect; {formatted_printout_text}"

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
    review_value = clinvar_review_status_to_level.get(variant.clinvar_review_status, 0)
    if review_value >= pathogenic_thresholds['pp5_clinvar_minimum_review'] and "pathogenic" in variant.clinvar_significance:
        score = review_value
        if score == 1:
            pp5 = f"PP5: Variant was found reported as {variant.clinvar_significance} " + \
                    f"in ClinVar as {variant.clinvar_id} with review status of {variant.clinvar_review_status}."
        else:
            pp5 = f"PP5_{score_to_hum_readable[score]}: Variant was found reported as {variant.clinvar_significance} " + \
                    f"in ClinVar as {variant.clinvar_id} with review status of {variant.clinvar_review_status}."
    return score, pp5

def get_pp(variant, missense_pathogenic_gene_to_region_list, chrom_to_pos_to_alt_to_splice_score, skip_list):
    """
    Supporting evidence of pathogenicity
    """
    pp1 = get_pp1()
    if "pp2" not in skip_list:
        score_2, pp2 = get_pp2(variant, missense_pathogenic_gene_to_region_list)
    else:
        score_2 = 0
        pp2 = ""
    if "pp3" not in skip_list:
        score_3, pp3 = get_pp3(variant, chrom_to_pos_to_alt_to_splice_score)
    else:
        score_3 = 0
        pp3 = ""
    pp4 = get_pp4()
    if "pp5" not in skip_list:
        score_5, pp5 = get_pp5(variant)
    else:
        score_5 = 0
        pp5 = ""
    pp_rationale_list = {
        'pp1': (0, pp1),
        'pp2': (score_2, pp2),
        'pp3': (score_3, pp3),
        'pp4': (0, pp4),
        'pp5': (score_5, pp5)
    }
    return pp_rationale_list
