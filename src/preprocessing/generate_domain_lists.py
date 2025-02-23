"""
Identify the domains associated with pathogenic variants and domains associated with benign variants.

Uses UniProt domains from UCSC - see the UCSC uniprot domains track detail for full credits and details.

UniProt Consortium. Reorganizing the protein space at the Universal Protein Resource (UniProt). Nucleic Acids Res.
    2012 Jan;40(Database issue):D71-5. PMID: 22102590; PMC: PMC3245120
"""
import sys
import os
import io
import gzip
import argparse
import logging

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from src.bias_2015.constants import clinvar_review_status_to_level, classification_mapping

def parseArgs(): 
    """
    Parse the command line arguments into useful python objects.  '--' variables are optional
    set_defaults only applies when the argument is not provided (it won't override)
    """
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument("clinvar_vcf",
                        help = " Input file 1",
                        action = "store")
    parser.add_argument("uniprot_bed",
                        help = " Input file 2",
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
    logging.basicConfig(level=getattr(logging, options.verbose), format='%(message)s')
    return options

def open_file(file_path, mode):
    """
    Open either a normal or a .gz file
    """
    _, file_extension = os.path.splitext(file_path)
    if file_extension == ".gz":
        return gzip.open(file_path, mode)
    return io.open(file_path, mode, encoding="utf-8")

def get_uniprot_domain_to_region_list(unip_bed):
    """
    Create a python structure holding the SwissProt curated uniprot domain locational information
    
    SwissProt is a manually curated protein sequence database that provides high-quality annotations based on experimental
    evidence and expert review. Its rigorous curation process ensures reliable and well-established functional insights,
    making it a trusted resource for identifying critical protein domains and their biological significance.
    """
    chrom_to_uniprot_domain_to_region_list = {}
    domain_name_to_annotation_source = {}
    with open(unip_bed, 'r') as in_file:
        for line in in_file:
            split_line = line.strip().split("\t")
            chrom, start, end = split_line[:3]
            start = int(start)
            end = int(end)
            strand = split_line[5]
            block_count, block_sizes, chrom_starts = split_line[9:12]
            annotation_source = split_line[16]
            domain_full_name = split_line[22]
            domain_name_to_annotation_source[domain_full_name] = annotation_source
            if int(block_count) > 1:
                split_block_sizes = block_sizes.split(",")
                split_chrom_starts = chrom_starts.split(",")
                for block_size, chrom_start in zip(split_block_sizes, split_chrom_starts):
                    if strand == "+":
                        if chrom_to_uniprot_domain_to_region_list.get(chrom):
                            if chrom_to_uniprot_domain_to_region_list[chrom].get(domain_full_name):
                                chrom_to_uniprot_domain_to_region_list[chrom][domain_full_name].append(
                                        (chrom, start + int(chrom_start), start + int(chrom_start) + int(block_size)))
                            else:
                                chrom_to_uniprot_domain_to_region_list[chrom][domain_full_name] = [(chrom, start + int(chrom_start), start + int(chrom_start) + int(block_size))]
                        else:
                            chrom_to_uniprot_domain_to_region_list[chrom] = {domain_full_name: [(chrom, start, end)]}
                    else:
                        if chrom_to_uniprot_domain_to_region_list.get(chrom):
                            if chrom_to_uniprot_domain_to_region_list[chrom].get(domain_full_name):
                                chrom_to_uniprot_domain_to_region_list[chrom][domain_full_name].append((chrom,  end - int(chrom_start) - int(block_size), end - int(chrom_start)))
                            else:
                                chrom_to_uniprot_domain_to_region_list[chrom][domain_full_name] = [(chrom, end - int(chrom_start) - int(block_size), end - int(chrom_start))]
                        else:
                            chrom_to_uniprot_domain_to_region_list[chrom] = {domain_full_name: [(chrom, start, end)]}
            else:
                if chrom_to_uniprot_domain_to_region_list.get(chrom):
                    if chrom_to_uniprot_domain_to_region_list.get(domain_full_name):
                        chrom_to_uniprot_domain_to_region_list[chrom][domain_full_name].append((chrom, start, end))
                    else:
                        chrom_to_uniprot_domain_to_region_list[chrom][domain_full_name] = [(chrom, start, end)]
                else:
                    chrom_to_uniprot_domain_to_region_list[chrom] = {domain_full_name: [(chrom, start, end)]}
    return chrom_to_uniprot_domain_to_region_list, domain_name_to_annotation_source

def get_uniprot_accession_to_clinvar_significance(uniprot_accession_to_variant_list):
    """
    Take in a dict mapping uniprot accessions to a list of variants. For each uniprot accession (domain)
    calculate a score for each clinical significance category (benign, likely benign, uncertain, 
    likely pahogenic, pathogenic). The score is reflective of the quantity and quality of variants with
    that clinical significance.
    """
    # Map clinvar review status to a weighted integer. Weight more if the review is of higher authority.
    uniprot_accession_to_clinvar_significance = {}
    for domain, variant_list in uniprot_accession_to_variant_list.items():
        for var in variant_list:
            data = var[4]
            if data.get('CLNSIG'):
                # Variants with more supporting evidence are weighted higher
                weighted_value = clinvar_review_status_to_level[data['CLNREVSTAT'].replace("_", " ")]
                if weighted_value > 0:
                    if uniprot_accession_to_clinvar_significance.get(domain):
                        if uniprot_accession_to_clinvar_significance[domain].get(data['CLNSIG']):
                            uniprot_accession_to_clinvar_significance[domain][data['CLNSIG']] += weighted_value
                        else:
                            uniprot_accession_to_clinvar_significance[domain][data['CLNSIG']] = weighted_value
                    else:
                        uniprot_accession_to_clinvar_significance[domain] = {data['CLNSIG']: weighted_value}
    return uniprot_accession_to_clinvar_significance


def get_uniprot_domain_to_variant_list(in_vcf, chrom_to_uniprot_domain_to_region_list):
    """
    take in a clinvar vcf and generate a dictionary mapping the UniProt identifier
    to the variant information
    """
    uniprot_accession_to_variant_list = {}
    with open(in_vcf) as in_file:
        count = 0
        for line in in_file:
            if not line: continue
            if line.startswith("#"): continue
            count += 1
            if count % 100000 == 0:
                logging.debug("%i", count)
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
            pos = int(pos)
            uniprot_domain_to_region_list = chrom_to_uniprot_domain_to_region_list.get(chrom)
            if not uniprot_domain_to_region_list:
                continue
            for uniprot_domain, region_list in uniprot_domain_to_region_list.items():
                for region in region_list:
                    _, r_start, r_end = region
                    if r_start <= pos <= r_end:
                        if uniprot_accession_to_variant_list.get(uniprot_domain):
                            uniprot_accession_to_variant_list[uniprot_domain].append(variant)
                        else:
                            uniprot_accession_to_variant_list[uniprot_domain] = [variant]
    return uniprot_accession_to_variant_list


def identify_domains(uniprot_accession_to_clinvar_significance, chrom_to_uniprot_domain_to_region_list, dom_name_to_annot_source):
    """
    Identify the domains that are associated with pathogenic and benign variations
    """
    pathogenic_domains = [] # >80% pathogenic or likely pathogenic variants
    benign_domains = [] # >80% benign or likely benign variants
    missing_mappings = set()
    for uniprot_access, clinvar_significance in uniprot_accession_to_clinvar_significance.items():
        total_score = 0
        uncertain_vars, pathogenic_variants, benign_variants = 0, 0 ,0 
        pathogenic_score = 0
        benign_score = 0
        # Calculate totals across the broad categories of 'pathogenic' and 'benign'
        for sig, score in clinvar_significance.items():
            if not classification_mapping.get(sig):
                missing_mappings.add(sig)

            if classification_mapping.get(sig, 'uncertain') == 'pathogenic':
                pathogenic_score += score
                total_score += score
                pathogenic_variants += 1
            elif classification_mapping.get(sig, 'uncertain') == 'benign':
                benign_score += score
                benign_variants += 1
                total_score += score
            else:
                uncertain_vars += 1
                total_score += score
        path_per = float(pathogenic_variants/(pathogenic_variants + benign_variants + uncertain_vars)) 
        benign_per = float(benign_variants/(pathogenic_variants + benign_variants + uncertain_vars)) 
        path_ben_ratio = 0
        if benign_variants:
            path_ben_ratio = pathogenic_variants/benign_variants
        # Require at least 5 score to classify a domain, ensure that domains with limited supporting
        # variants are not falsely included.
        if total_score < 5: continue
        if pathogenic_variants < 1: continue
        if float(uncertain_vars/(pathogenic_variants + benign_variants + uncertain_vars)) > .5: # If the gene has over 50% uncertain variants, then disregard it
            continue
        if path_ben_ratio < 1 or benign_per > .33 or path_per > .66:
            continue

        region_list = []
        for chrom, uniprot_domain_to_region_list in chrom_to_uniprot_domain_to_region_list.items():
            if uniprot_domain_to_region_list.get(uniprot_access):
                region_list = uniprot_domain_to_region_list[uniprot_access]
        for region in region_list:
            chrom, start, end = region
            pathogenic_domains.append((chrom, str(start), str(end), uniprot_access, dom_name_to_annot_source[uniprot_access],
                                           str(path_per), str(pathogenic_score), str(benign_per), str(benign_score), str(total_score),
                                       str(path_ben_ratio)
                                       )
                                      )
        
        # NOTE: this is unused currently, the above filters may not make this applicable 
        if benign_per >.5 and path_per <.1:
            for region in region_list:
                chrom, start, end = region
                benign_domains.append((chrom, str(start), str(end), uniprot_access, dom_name_to_annot_source[uniprot_access], 
                                       str(path_per), str(pathogenic_score), str(benign_per), str(benign_score), str(total_score)))
    if len(missing_mappings) > 0:
        print(missing_mappings)
    return benign_domains, pathogenic_domains



def evaluate_clinvar_domains(in_vcf, uniprot_bed, pathogenic_output):
    """
    Evaluate the clinvar domain information. Identify the variants associated with each domain.
    
    For each domain, calculate a weighted score for each clinvar significance seen

    For each domain, identify the total pathogenic, uncertain and benign scores. If the domain
    has a pathogenic percentage greater than 80% it is added to the pathogenic domains list.
    If it has greater than 80% benign percentage, it is added to the benign domains list.
    """
    # A mapping of uniprot domains to precise regions associated with them
    chrom_to_uniprot_domain_to_region_list, dom_name_to_annot_source = get_uniprot_domain_to_region_list(uniprot_bed)

    # A mapping of domains (uniprot id) to all variant in the domain
    uniprot_domain_to_variant_list = get_uniprot_domain_to_variant_list(in_vcf, chrom_to_uniprot_domain_to_region_list)

    # A mapping of domains (uniprot id) to significance (Ex 'benign|confers_sensitivity') to a weighted score
    uniprot_domain_to_clinvar_significance = get_uniprot_accession_to_clinvar_significance(uniprot_domain_to_variant_list)

    # Identify the domains associated with pathogenic and benign variations
    benign_domains, pathogenic_domains = identify_domains(uniprot_domain_to_clinvar_significance,
                                                          chrom_to_uniprot_domain_to_region_list, dom_name_to_annot_source)

    with open(pathogenic_output, 'w') as o_file:
        for domain_data in sorted(pathogenic_domains, key=lambda x: x[0]):
            o_file.write("\t".join(domain_data) + "\n")

    with open("benign_domains.tsv", 'w') as o_file:
        for domain_data in sorted(benign_domains, key=lambda x: x[0]):
            o_file.write("\t".join(domain_data) + "\n")

def main():
    """
    main function, calls other functions
    """
    options = parseArgs()
    evaluate_clinvar_domains(options.clinvar_vcf, options.uniprot_bed, options.output)

if __name__ == "__main__":
    sys.exit(main())
