"""
Parse loss of function information from genes from Gnomad data, we store the LOEUF value for genes with LOEUF < 2


LOEUF (Loss-of-Function Observed/Expected Upper Fraction) values from the gnomAD database quantify a gene's tolerance
to loss-of-function (LoF) variants, with lower values indicating higher intolerance. These values are used to prioritize
genes in rare disease studies and functional analyses, as they highlight genes under strong purifying selection. LOEUF
scores, available through gnomAD, range from near 0 (high intolerance) to >1 (greater tolerance), aiding gene
prioritization and population genetics studies.

"""
import argparse

def filter_lof_genes(input_file, output_file, ref_build, loeuf_threshold=2.0):
    """
    Identify genes that have a sufficient LOEUF to be considered genes where loss of function is a mechanism of disease

    "We also introduced the loss-of-function observed/expected upper bound fraction (LOEUF) metric, which uses the
     upper bound of the Poisson 90% confidence interval for the number of expected LoF variants in each gene."
    (Karczewski et al. Nature 581, 434–443 (2020))
    """
    # Open the input and output files
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Read the header line
        header = infile.readline().strip().split('\t')
        # Find the indices for gene and pLI columns
        gene_index = header.index('gene')
        if ref_build == "hg19":
            loeuf_index = header.index('oe_lof_upper')
        else:
            loeuf_index = header.index('lof.oe_ci.upper')

        # Write the header for the output file
        outfile.write('gene\tloeuf\n')
        
        # Process each line in the input file
        for line in infile:
            fields = line.strip().split('\t')
            gene = fields[gene_index]
            try:
                loeuf = float(fields[loeuf_index])
            except ValueError:
                continue

            # Check if the pLI score meets the threshold
            if loeuf <= loeuf_threshold:
                outfile.write(f'{gene}\t{loeuf}\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter genes based on pLI score.')
    parser.add_argument('input_file', type=str, help='Path to the input file.')
    parser.add_argument('output_file', type=str, help='Path to the output file.')
    parser.add_argument('ref_build', type=str, help='The reference genome build', choices=['hg19', 'hg38'])
    parser.add_argument('--loeuf_threshold', type=float, default=2.0, help='Threshold for LOEUF score (default: 2.0).')

    args = parser.parse_args()
    filter_lof_genes(args.input_file, args.output_file, args.ref_build, args.loeuf_threshold)
