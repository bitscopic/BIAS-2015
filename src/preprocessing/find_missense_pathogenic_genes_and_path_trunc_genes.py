"""
Take a clinvar vcf nirvana json file and generate a tsv file with the pathogenic
aa changes
"""
import sys
import os
import io
import gzip
import argparse
import logging

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from src.bias_2015.constants import clinvar_review_status_to_level

def parseArgs(): 
    """
    Parse the command line arguments into useful python objects.  '--' variables are optional
    set_defaults only applies when the argument is not provided (it won't override)
    """
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument("clinvar_vcf",
                        help = "Clinvar VCF",
                        action = "store")
    parser.add_argument("gnomad_rmc_file",
                        help = "GNOMAD regional constraint metrics",
                        action = "store")
    parser.add_argument("output_file",
                        help = " Missense pathogenic genes",
                        action = "store")
    parser.add_argument("output_file2",
                        help = " Truncating pathogenic genes",
                        action = "store")
    parser.add_argument("ref_b",
                        help = " The reference build",
                        choices = ['hg19', 'hg38'],
                        action = "store")
    parser.add_argument("--verbose",
                        help = " The verbosity level for stdout messages (default INFO)",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        action = "store")

    parser.set_defaults(verbose = "WARNING")
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

def overlap_gnomad_rmc(clinvar, chrom_to_gnomad_rmc_list):
    """
    Identify if a gnomad regional missense constraint region overlaps
    this variant
    """
    # Unpack clinvar
    clinvar_chrom, clinvar_pos, _, _ = clinvar
    clinvar_pos = int(clinvar_pos)

    # Loop through chrom_to_gnomad_rmc_list and check for overlap
    gnomad_rmc_list = chrom_to_gnomad_rmc_list.get(clinvar_chrom)
    if gnomad_rmc_list:
        for start, end, gene, oe in gnomad_rmc_list:
            if start <= clinvar_pos <= end:
                return gene, start, end, oe
    return "", "", "", ""

def get_gene_to_signif_to_variant_consequence_list(clinvar_to_data, chrom_to_gnomad_rmc_list):
    """
    Extract the relevant fields from the clinvar data and format it to be namespaced under each gene

    ex; {BRCA: {'Pathogenic': ['missense', 'expert reviewed'], 'Benign':['missense', 'single reviewer']} ... }
    """
    gene_to_signif_to_variant_consequence_list = {}
    for clinvar, data, in clinvar_to_data.items():
        
        gnomad_gene, gnomad_start, gnomad_end, gnomad_oe = overlap_gnomad_rmc(clinvar, chrom_to_gnomad_rmc_list)
        gene_info = data.get('GENEINFO')
        if gene_info:
            signif = data.get('CLNSIG')
            consequence = data.get('MC')
            review_rank = data.get('CLNREVSTAT')
            if gnomad_gene:
                gene_info = (gnomad_gene, gnomad_start, gnomad_end, gnomad_oe)
            if gene_to_signif_to_variant_consequence_list.get(gene_info):
                if gene_to_signif_to_variant_consequence_list[gene_info].get(signif):
                    gene_to_signif_to_variant_consequence_list[gene_info][signif].append((consequence, review_rank))
                else:
                    gene_to_signif_to_variant_consequence_list[gene_info][signif] = [(consequence, review_rank)]
            else:
                gene_to_signif_to_variant_consequence_list[gene_info] = {signif: [(consequence, review_rank)]}
    return gene_to_signif_to_variant_consequence_list


def get_genes_where_missense_variants_are_often_pathogenic(gene_to_signif_to_variant_consequence_list, chrom_to_gnomad_rmc_list): 
    """
    PP2    Missense variants in a gene that has a low rate of benign missense variation
           and where missense variants are a common mechanism of disease 

    Looking for genes where over 80% of pathogenic variants are missense, and
    when looking at all missense variants in the gene, benign variants make up less than 10%

    Also evaluate the GNOMAD regional mutation constraint obs/exp (oe) ratio. Regions that have an O/E less
    than 1 indicate an intolerance to missense mutation.
    """
    # Prepopulate with the RMC data, merging is handled at the end of line extraction
    genes_where_missense_variants_are_often_pathogenic = set()
    for _, gnomad_rmc_list in chrom_to_gnomad_rmc_list.items():
        for gnomad_rmc in gnomad_rmc_list:
            gnomad_start, gnomad_end, gene_name, gnomad_oe = gnomad_rmc
            genes_where_missense_variants_are_often_pathogenic.add((gene_name, str(gnomad_start), str(gnomad_end), 
                                                                    gnomad_oe, '0', '0'))
    for gene, signif_to_variant_consequence_list in gene_to_signif_to_variant_consequence_list.items():
        missense_variants = 0.0
        pathogenic_variants = 0.0
        benign_variants = 0
        pathogenic_missense_variants = 0.0
        benign_missense_variants = 0.0
        uncertain_vars = 0
        for signif, variant_consequence_list in signif_to_variant_consequence_list.items():
            for var_cons in variant_consequence_list:
                cons, review_rating = var_cons
                if review_rating:
                    review_rating = review_rating.replace("_", " ")
                if not cons:
                    continue
                if clinvar_review_status_to_level.get(review_rating, None) is None:
                    continue
                if clinvar_review_status_to_level[review_rating] < 1: # Exclude variants that have limited evidence supporting them
                    continue
                if "missense" in cons:
                    missense_variants += 1
                if "benign" in signif.lower():
                    benign_variants += 1
                    if "missense" in cons:
                        benign_missense_variants += 1
                elif "pathogenic" in signif.lower():
                    pathogenic_variants += 1
                    if "missense" in cons:
                        pathogenic_missense_variants += 1
                else:
                    uncertain_vars += 1
        # Calculate core stats
        path_missense_percentage = 0
        if pathogenic_variants:
            path_missense_percentage = pathogenic_missense_variants/pathogenic_variants
        ben_missense_percentage = 0
        if missense_variants:
            ben_missense_percentage = benign_missense_variants/missense_variants
        gnomad_start = 0
        gnomad_end = 0
        gnomad_oe = ""
        if len(gene) == 4:
            gene_name = gene[0]
            gnomad_start = gene[1]
            gnomad_end = gene[2]
            gnomad_oe = gene[3]
        else:
            gene_name = gene.split(":")[0]

        # Begin filtering, if we encounter a gnomad curated entry, we use OE
        if gnomad_oe:
            if float(gnomad_oe) >= .4: # Indicates that the region is not tolerant to missense mutation 
                continue
        else:
            if not missense_variants: # Require at least one missense_variant
                continue
            if pathogenic_missense_variants < 1:  # Require multiple pathogenic missense variants
                continue
            if float(uncertain_vars/(pathogenic_variants + benign_variants + uncertain_vars)) > .5: # If the gene has over 50% uncertain variants, then disregard it
                continue
            if (pathogenic_variants + benign_variants + uncertain_vars) < 5: # Need at least 5 variants to consider a gene
                continue 
            if ben_missense_percentage > .25: # Benign missense variants are too common
                continue
            if path_missense_percentage < .5: # Pathogenic missense variants are too rare
                continue
        # This handles merging the Gnomad RMC with our ClinVar derived regions
        if (gene_name, str(gnomad_start), str(gnomad_end), str(gnomad_oe), '0', '0') in genes_where_missense_variants_are_often_pathogenic:
            genes_where_missense_variants_are_often_pathogenic.remove((gene_name, str(gnomad_start), str(gnomad_end), str(gnomad_oe), '0', '0'))
        # Add the new entry
        genes_where_missense_variants_are_often_pathogenic.add((gene_name, str(gnomad_start), str(gnomad_end), 
                                                                    str(gnomad_oe), str(path_missense_percentage),
                                                                    str(ben_missense_percentage)))
    return sorted(list(genes_where_missense_variants_are_often_pathogenic), key=lambda x: x[0])

def get_genes_where_most_pathogenic_variants_are_truncating(gene_to_signif_to_variant_consequence_list):
    """
    BP1    Missense variants in a gene for which primarily truncating variants are
           known to cause disease
    
    Looking for genes where over 80% of pathogenic variants are truncating variants (stop-gain, stop-loss, frameshift
    indel, or those disrupting splice sites.)
    """
    genes_where_most_pathogenic_variants_are_truncating = set()
    for gene, signif_to_variant_consequence_list in gene_to_signif_to_variant_consequence_list.items():
        pathogenic_variants = 0.0
        pathogenic_truncating_variants = 0.0
        benign_variants = 0.0
        benign_truncating_variants = 0.0
        uncertain_vars = 0
        for signif, variant_consequence_list in signif_to_variant_consequence_list.items():
            for var_cons in variant_consequence_list:
                cons, review_rating = var_cons
                if review_rating:
                    review_rating = review_rating.replace("_", " ")
                if not cons:
                    continue
                if not clinvar_review_status_to_level.get(review_rating):
                    continue
                if clinvar_review_status_to_level[review_rating] < 1: # Exclude variants that have limited evidence supporting them
                    continue
                truncating_ontologies = ['stop_gained', 'splice_donor_variant', 'frameshift_variant', 'stop_lost'] 
                if "pathogenic" in signif.lower():
                    pathogenic_variants += 1
                    ont = cons.split("|")[1]
                    if ont in truncating_ontologies:
                        pathogenic_truncating_variants += 1
                elif "benign" in signif.lower():
                    benign_variants += 1
                    ont = cons.split("|")[1]
                    if ont in truncating_ontologies:
                        benign_truncating_variants += 1
                else:
                    uncertain_vars += 1
        path_trunc_per, ben_trunc_per = 0,0
        if pathogenic_variants:
            path_trunc_per = pathogenic_truncating_variants/pathogenic_variants
        if benign_variants:
            ben_trunc_per = benign_truncating_variants/benign_variants
        # Begin filtering
        if pathogenic_truncating_variants < 1: # We must see at least one pathogenic truncating variants
            continue
        if float(uncertain_vars/(pathogenic_variants + benign_variants + uncertain_vars)) > .5: # If the gene has over 50% uncertain variants, then disregard it
            continue
        if path_trunc_per < .75: # Truncating variants in this gene are not a majority of pathogenic variants
            continue
        if ben_trunc_per > .25: # Benign truncating variants are common in this gene
            continue
        if len(gene) == 4:
            gene_name = gene[0]
        else:
            gene_name = gene.split(":")[0]
        genes_where_most_pathogenic_variants_are_truncating.add((gene_name, str(path_trunc_per), str(ben_trunc_per)))
    return sorted(list(genes_where_most_pathogenic_variants_are_truncating), key=lambda x: x[0])


def get_chrom_to_gnomad_rmc_list(gnomad_rmc_file):
    """
    Takes in the UCSC Gnomad missense constraint browser file and extracts the useful information
    to a list.
    """
    chrom_to_gnomad_rmc_list = {}

    with open(gnomad_rmc_file, 'r') as in_file:
        for line in in_file:
            split_line = line.strip().split("\t")
            if len(split_line) < 15:
                continue
            chrom, start, end = split_line[:3]
            gene = split_line[12]
            oe_score = split_line[15]
            if chrom_to_gnomad_rmc_list.get(chrom):
                chrom_to_gnomad_rmc_list[chrom].append((int(start), int(end), gene, oe_score))
            else:
                chrom_to_gnomad_rmc_list[chrom] = [(int(start), int(end), gene, oe_score)]

    return chrom_to_gnomad_rmc_list


def find_missense_pathogenic_genes_and_path_trunc_genes(clinvar_vcf, gnomad_rmc_file, output_file, output_file_2, ref_b):
    """
    Find missense pathogenic_genes
    """
    # Load the clinvar data into memory
    clinvar_to_data = get_clinvar_to_data(clinvar_vcf)

    # The GNOMAD regional missense constraint information
    if ref_b == "hg19":
        chrom_to_gnomad_rmc_list = get_chrom_to_gnomad_rmc_list(gnomad_rmc_file)
    else:
        chrom_to_gnomad_rmc_list = {}


    # Extract the relevant fields from the clinvar data and format it to be namespaced under each gene.
    logging.info("Building gene to significanct to variant consequence lists dictionary") 
    gene_to_signif_to_variant_consequence_list = get_gene_to_signif_to_variant_consequence_list(clinvar_to_data, chrom_to_gnomad_rmc_list)

    # PP2
    genes_where_missense_variants_are_often_pathogenic = \
            get_genes_where_missense_variants_are_often_pathogenic(gene_to_signif_to_variant_consequence_list, chrom_to_gnomad_rmc_list) 

    # BP1
    genes_where_most_pathogenic_variants_are_truncating = \
            get_genes_where_most_pathogenic_variants_are_truncating(gene_to_signif_to_variant_consequence_list)

    with open(output_file, 'w') as o_file:
        for entry in genes_where_missense_variants_are_often_pathogenic:
            o_file.write("\t".join(entry) + "\n")

    with open(output_file_2, 'w') as o_file:
        for entry in genes_where_most_pathogenic_variants_are_truncating:
            o_file.write("\t".join(entry) + "\n")


def main():
    """
    Load ClinVar data into memory, then perform two different operations on it. One to determine genes
    where missense variants are often pathogen (for PP2) and another to determine genes where most pathogenic
    variants are truncating (for BP1). 
    """
    options = parseArgs()

    find_missense_pathogenic_genes_and_path_trunc_genes(options.clinvar_vcf, options.gnomad_rmc_file, options.output_file, options.output_file2, options.ref_b)


if __name__ == "__main__":
    sys.exit(main())
