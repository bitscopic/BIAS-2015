"""
Parse loss of function (lof) genes from Gnomad data
"""
import argparse

def filter_lof_genes(input_file, output_file, ref_build, pli_threshold=0.8):
    """
    Identify genes that have a sufficient pLi
    """
    # Open the input and output files
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Read the header line
        header = infile.readline().strip().split('\t')
        # Find the indices for gene and pLI columns
        gene_index = header.index('gene')
        if ref_build == "hg19":
            pli_index = header.index('pLI')
        else:
            pli_index = header.index('lof.pLI')

        # Write the header for the output file
        outfile.write('gene\tpLI\n')
        
        # Process each line in the input file
        for line in infile:
            fields = line.strip().split('\t')
            gene = fields[gene_index]
            try:
                pli = float(fields[pli_index])
            except ValueError:
                continue

            # Check if the pLI score meets the threshold
            if pli >= pli_threshold:
                outfile.write(f'{gene}\t{pli}\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter genes based on pLI score.')
    parser.add_argument('input_file', type=str, help='Path to the input file.')
    parser.add_argument('output_file', type=str, help='Path to the output file.')
    parser.add_argument('ref_build', type=str, help='The reference genome build', choices=['hg19', 'hg38'])
    parser.add_argument('--pli_threshold', type=float, default=0.8, help='Threshold for pLI score (default: 0.8).')

    args = parser.parse_args()
    filter_lof_genes(args.input_file, args.output_file, args.pli_threshold)
