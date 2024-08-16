"""
Identify the domains associated with pathogenic variants and domains associated with benign variants.

Write them to output files 'pathogenic_domains.tsv' and 'benign_domains.tsv'
"""
import sys
import os
import io
import gzip
import argparse


def parseArgs(): 
    """
    Parse the command line arguments into useful python objects.  '--' variables are optional
    set_defaults only applies when the argument is not provided (it won't override)
    """
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument("clinvar_vcf",
                        help = " Input file 1",
                        action = "store")
    parser.add_argument("output",
                        help = " Output file 1",
                        action = "store")
    parser.add_argument("--verbose",
                        help = " The verbosity level for stdout messages (default INFO)",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        action = "store")

    parser.set_defaults(verbose = "INFO")

    options = parser.parse_args()
    return options

def open_file(file_path, mode):
    """
    Open either a normal or a .gz file
    """
    _, file_extension = os.path.splitext(file_path)
    if file_extension == ".gz":
        return gzip.open(file_path, mode)
    return io.open(file_path, mode, encoding="utf-8")


def get_uniprot_accession_to_clinvar_significance(uniprot_accession_to_variant_list):
    """
    Take in a dict mapping uniprot accessions to a list of variants. For each uniprot accession (domain)
    calculate a score for each clinical significance category (benign, likely benign, uncertain, 
    likely pahogenic, pathogenic). The score is reflective of the quantity and quality of variants with
    that clinical significance.
    """
    # Map clinvar review status to a weighted integer. Weight more if the review is of higher authority.
    review_rankings = {
        'practice_guideline': 3,
        'reviewed_by_expert_panel': 3,
        'criteria_provided,_multiple_submitters,_no_conflicts': 2,
        'criteria_provided,_single_submitter': 1,
        'criteria_provided,_conflicting_classifications': 0,
        'no_assertion_criteria_provided': 0,
        'no_classification_provided': 0,
        'no_classification_for_the_single_variant': 0
    }
    uniprot_accession_to_clinvar_significance = {}
    for domain, variant_list in uniprot_accession_to_variant_list.items():
        for var in variant_list:
            data = var[4]
            if data.get('CLNSIG'):
                # Variants with more supporting evidence are weighted higher
                weighted_value = review_rankings[data['CLNREVSTAT']]
                if weighted_value > 0:
                    if uniprot_accession_to_clinvar_significance.get(domain):
                        if uniprot_accession_to_clinvar_significance[domain].get(data['CLNSIG']):
                            uniprot_accession_to_clinvar_significance[domain][data['CLNSIG']] += weighted_value
                        else:
                            uniprot_accession_to_clinvar_significance[domain][data['CLNSIG']] = weighted_value
                    else:
                        uniprot_accession_to_clinvar_significance[domain] = {data['CLNSIG']: weighted_value}
    return uniprot_accession_to_clinvar_significance


def get_uniprot_accession_to_variant_list(inVcf):
    """
    take in a clinvar vcf and generate a dictionary mapping the UniProt identifier
    to the variant information
    """
    uniprot_accession_to_variant_list = {} 
    with open(inVcf) as in_file:
        for line in in_file:
            if not line: continue
            if line.startswith("#"): continue
            split_line = line.strip().split()
            chrom = split_line[0]
            pos = split_line[1]
            ref = split_line[3]
            alt = split_line[4]
            if ref[0] == alt[0]:
                if len(ref) > 1:
                    ref = ref[1:]
                else:
                    ref = "-"
                if len(alt) > 1:
                    alt = alt[1:]
                else:
                    alt = "-"
            info_data = {item.split('=')[0]: item.split('=')[1] for item in split_line[7].split(";")}
            variant = (chrom, pos, ref, alt, info_data)
            if "CLNVI" in info_data:
                clin_data = {item.split(':')[0]: item.split(':')[1] for item in info_data['CLNVI'].split("|")}
                if "UniProtKB" in clin_data:
                    uniprot_accession = clin_data['UniProtKB'].split('#')[0]
                    if uniprot_accession_to_variant_list.get(uniprot_accession):
                        uniprot_accession_to_variant_list[uniprot_accession].append(variant)
                    else:
                        uniprot_accession_to_variant_list[uniprot_accession] = [variant]
    return uniprot_accession_to_variant_list


def identify_domains(uniprot_accession_to_clinvar_significance):
    """
    Identify the domains that are associated with pathogenic and benign variations
    """
    
    # All ClinVar classifications. Each is mapped to one of three values; benign, uncertain and pathogenic
    # the 'likely' values have been collapsed into their main category. 
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
        'Uncertain_significance|drug_response': 'uncertain',
        'protective|risk_factor': 'uncertain',
        'Benign/Likely_benign|risk_factor': 'benign',
        'Conflicting_classifications_of_pathogenicity': 'uncertain'
    }

    pathogenic_domains = [] # >80% pathogenic or likely pathogenic variants
    benign_domains = [] # >80% benign or likely benign variants
    for uniprot_access, clinvar_significance in uniprot_accession_to_clinvar_significance.items():
        total_score = 0
        pathogenic_score = 0
        benign_score = 0
        # Calculate totals across the broad categories of 'pathogenic' and 'benign'
        for sig, score in clinvar_significance.items():
            if classification_mapping[sig] == 'pathogenic':
                pathogenic_score += score
                total_score += score
            elif classification_mapping[sig] == 'benign':
                benign_score += score
                total_score += score
            
        # Require at least 5 score to classify a domain, ensure that domains with limited supporting
        # variants are not falsely included.
        if total_score < 5: continue

        path_per = float(pathogenic_score/total_score)
        benign_per = float(benign_score/total_score)
        if path_per > .8:
            pathogenic_domains.append((uniprot_access, path_per, pathogenic_score))
        if benign_per >.8:
            benign_domains.append((uniprot_access, benign_per, benign_score))
    return benign_domains, pathogenic_domains


def evaluate_clinvar_domains(in_vcf, pathogenic_output):
    """
    Evaluate the clinvar domain information. Identify the variants associated with each domain.
    
    For each domain, calculate a weighted score for each clinvar significance seen

    For each domain, identify the total pathogenic, uncertain and benign scores. If the domain
    has a pathogenic percentage greater than 80% it is added to the pathogenic domains list.
    If it has greater than 80% benign percentage, it is added to the benign domains list.
    """
    # A mapping of domains (uniprot id) to all variant in the domain
    uniprot_accession_to_variant_list = get_uniprot_accession_to_variant_list(in_vcf)

    # A mapping of domains (uniprot id) to significance (Ex 'benign|confers_sensitivity') to a weighted score
    uniprot_accession_to_clinvar_significance = get_uniprot_accession_to_clinvar_significance(uniprot_accession_to_variant_list)

    # Identify the domains associated with pathogenic and benign variations
    benign_domains, pathogenic_domains = identify_domains(uniprot_accession_to_clinvar_significance)

    with open("benign_domains.tsv", 'w') as o_file:
        for domain_data in sorted(benign_domains, key=lambda x: x[0]):
            o_file.write(f"{domain_data[0]}\t{domain_data[1]}\t{domain_data[2]}\n")

    with open(pathogenic_output, 'w') as o_file:
        for domain_data in sorted(pathogenic_domains, key=lambda x: x[0]):
            o_file.write(f"{domain_data[0]}\t{domain_data[1]}\t{domain_data[2]}\n")

def main():
    """
    main function, calls other functions
    """
    options = parseArgs()
    evaluate_clinvar_domains(options.clinvar_vcf, options.output)

if __name__ == "__main__":
    sys.exit(main())
