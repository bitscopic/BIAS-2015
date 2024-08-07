"""
Take a clinvar vcf nirvana json file and generate a tsv file with the pathogenic
aa changes
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
                        help = "Clinvar VCF",
                        action = "store")
    parser.add_argument("output_file",
                        help = " Missense pathogenic genes",
                        action = "store")
    parser.add_argument("output_file2",
                        help = " Truncating pathogenic genes",
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


def get_clinvar_to_data(inVcf):
    """
    take in a clinvar vcf and generate a dictionary mapping the UniProt identifier
    to the variant information
    """
    clinvar_to_data = {}
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
            variant = (chrom, pos, ref, alt)
            clinvar_to_data[variant] = info_data
    return clinvar_to_data


def get_gene_to_signif_to_variant_consequence_list(clinvar_to_data):
    """
    Extract the relevant fields from the clinvar data and format it to be namespaced under each gene

    ex; {BRCA: {'Pathogenic': ['missense', 'expert reviewed'], 'Benign':['missense', 'single reviewer']} ... }
    """
    gene_to_signif_to_variant_consequence_list = {}
    for _, data, in clinvar_to_data.items():
        gene_info = data.get('GENEINFO')
        if gene_info:
            signif = data.get('CLNSIG')
            consequence = data.get('MC')
            review_rank = data.get('CLNREVSTAT')
            if gene_to_signif_to_variant_consequence_list.get(gene_info):
                if gene_to_signif_to_variant_consequence_list[gene_info].get(signif):
                    gene_to_signif_to_variant_consequence_list[gene_info][signif].append((consequence, review_rank))
                else:
                    gene_to_signif_to_variant_consequence_list[gene_info][signif] = [(consequence, review_rank)]
            else:
                gene_to_signif_to_variant_consequence_list[gene_info] = {signif: [(consequence, review_rank)]}
    return gene_to_signif_to_variant_consequence_list


def get_genes_where_missense_variants_are_often_pathogenic(gene_to_signif_to_variant_consequence_list): 
    """
    PP2    Missense variants in a gene that has a low rate of benign missense variation
        and where missense variants are a common mechanism of disease 

    Looking for genes where over 80% of pathogenic variants are missense, and
    when looking at all missense variants in the gene, benign variants make up less than 10%
    """
    
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
   
    genes_where_missense_variants_are_often_pathogenic = set()
    for gene, signif_to_variant_consequence_list in gene_to_signif_to_variant_consequence_list.items():
        missense_variants = 0.0
        pathogenic_variants = 0.0
        pathogenic_missense_variants = 0.0
        benign_missense_variants = 0.0
        for signif, variant_consequence_list in signif_to_variant_consequence_list.items():
            for var_cons in variant_consequence_list:
                cons, review_rating = var_cons
                if not cons:
                    continue
                if not review_rankings.get(review_rating):
                    continue
                if review_rankings[review_rating] < 1: # Exclude variants that have limited evidence supporting them
                    continue
                if "missense" in cons:
                    missense_variants += 1
                if "benign" in signif.lower():
                    if "missense" in cons:
                        benign_missense_variants += 1
                if "pathogenic" in signif.lower():
                    pathogenic_variants += 1
                    if "missense" in cons:
                        pathogenic_missense_variants += 1
        if not pathogenic_variants: continue
        if not missense_variants: continue
        path_missense_percentage = pathogenic_missense_variants/pathogenic_variants
        ben_missense_percentage = benign_missense_variants/missense_variants
        if path_missense_percentage > .8 and ben_missense_percentage < .1:
            genes_where_missense_variants_are_often_pathogenic.add((gene.split(":")[0], str(path_missense_percentage), str(ben_missense_percentage)))
    return sorted(list(genes_where_missense_variants_are_often_pathogenic), key=lambda x: x[0])

def get_genes_where_most_pathogenic_variants_are_truncating(gene_to_signif_to_variant_consequence_list):
    """
    BP1    Missense variants in a gene for which primarily truncating variants are
        known to cause disease
    
    Looking for genes where over 80% of pathogenic variants are truncating variants (stop-gain, stop-loss, frameshift
    indel, or those disrupting splice sites.)
    """
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
   
    genes_where_most_pathogenic_variants_are_truncating = set()
    for gene, signif_to_variant_consequence_list in gene_to_signif_to_variant_consequence_list.items():
        pathogenic_variants = 0.0
        pathogenic_truncating_variants = 0.0
        for signif, variant_consequence_list in signif_to_variant_consequence_list.items():
            for var_cons in variant_consequence_list:
                cons, review_rating = var_cons
                if not cons:
                    continue
                if not review_rankings.get(review_rating):
                    continue
                if review_rankings[review_rating] < 1: # Exclude variants that have limited evidence supporting them
                    continue
                truncating_ontologies = ['stop_gained', 'splice_donor_variant', 'frameshift_variant', 'stop_lost', 'splice_acceptor_variant']
                if "pathogenic" in signif.lower():
                    pathogenic_variants += 1
                    ont = cons.split("|")[1]
                    if ont in truncating_ontologies:
                        pathogenic_truncating_variants += 1
        if not pathogenic_truncating_variants: continue
        path_trunc_per = pathogenic_truncating_variants/pathogenic_variants
        if path_trunc_per > .8:
            genes_where_most_pathogenic_variants_are_truncating.add((gene.split(":")[0], str(path_trunc_per)))
    return sorted(list(genes_where_most_pathogenic_variants_are_truncating), key=lambda x: x[0])


def main():
    """
    Load ClinVar data into memory, then perform two different operations on it. One to determine genes
    where missense variants are often pathogen (for PP2) and another to determine genes where most pathogenic
    variants are truncating (for BP1). 
    """
    options = parseArgs()

    # Load the clinvar data into memory
    clinvar_to_data = get_clinvar_to_data(options.clinvar_vcf)

    # Extract the relevant fields from the clinvar data and format it to be namespaced under each gene.
    gene_to_signif_to_variant_consequence_list = get_gene_to_signif_to_variant_consequence_list(clinvar_to_data)

    # PP2
    genes_where_missense_variants_are_often_pathogenic = \
            get_genes_where_missense_variants_are_often_pathogenic(gene_to_signif_to_variant_consequence_list) 

    # BP1
    genes_where_most_pathogenic_variants_are_truncating = \
            get_genes_where_most_pathogenic_variants_are_truncating(gene_to_signif_to_variant_consequence_list)

    with open(options.output_file, 'w') as o_file:
        for entry in genes_where_missense_variants_are_often_pathogenic:
            o_file.write("\t".join(entry) + "\n")

    with open(options.output_file2, 'w') as o_file:
        for entry in genes_where_most_pathogenic_variants_are_truncating:
            o_file.write("\t".join(entry) + "\n")


if __name__ == "__main__":
    sys.exit(main())
