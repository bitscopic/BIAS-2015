"""
Converts Erepo ACMG curated variant information in tsv format to hg19 DRAGEN vcf format

Erepo
https://erepo.clinicalgenome.org/evrepo/

"""
import pysam

def get_dna_sequence(chrom, start, end):
    """
    Return an hg19 dna sequence
    """
    reference_genome = "/Users/ceisenhart/Bioinformatics/data/reference_genomes/hg19.fa"
    # Open the reference genome file
    fasta = pysam.FastaFile(reference_genome)

    # Fetch the sequence
    sequence = fasta.fetch(chrom, start, end)

    # Close the fasta file
    fasta.close()

    return sequence.upper()

hg19_accessions = [
    "NC_000001.10",  # Chromosome 1
    "NC_000002.11",  # Chromosome 2
    "NC_000003.11",  # Chromosome 3
    "NC_000004.11",  # Chromosome 4
    "NC_000005.9",   # Chromosome 5
    "NC_000006.11",  # Chromosome 6
    "NC_000007.13",  # Chromosome 7
    "NC_000008.10",  # Chromosome 8
    "NC_000009.11",  # Chromosome 9
    "NC_000010.10",  # Chromosome 10
    "NC_000011.9",   # Chromosome 11
    "NC_000012.11",  # Chromosome 12
    "NC_000013.10",  # Chromosome 13
    "NC_000014.8",   # Chromosome 14
    "NC_000015.9",   # Chromosome 15
    "NC_000016.9",   # Chromosome 16
    "NC_000017.10",  # Chromosome 17
    "NC_000018.9",   # Chromosome 18
    "NC_000019.9",   # Chromosome 19
    "NC_000020.10",  # Chromosome 20
    "NC_000021.8",   # Chromosome 21
    "NC_000022.10",  # Chromosome 22
    "NC_000023.10",  # Chromosome X
    "NC_000024.9"    # Chromosome Y
]

def generate_vcf_data(erepo_filename, verbose = False):
    """

    NC_000010.9:g.89680808G>A,
    NC_000010.9:g.89682800_89682802dup,
    NC_000010.9:g.89701896_89701897delinsAT,
    NC_000010.9:g.89715049_89715051del,
    NC_000017.9:g.75696123_75696124insAGCGGGC,
    """
    # List of chrom, pos, ref, alt
    vcf_data = []
    with open(erepo_filename, 'r') as erepo_file:
        next(erepo_file)  # Skip header line
        int_str = [str(i) for i in range(10)]
        line_count = 0
        for line in erepo_file:
            line_count += 1
            columns = line.split('\t')
            no_hg19_coordinates = True
            for gen_id in columns[3].split(" "):
                if gen_id.startswith("NC_0000"):
                    components = gen_id.strip(",").split(":")
                    if "?" in gen_id: 
                        continue # Ignore super strange ones
                    if components[0] in hg19_accessions:
                        # Go through the second half of the g. notation and get the first location. For SNP is is the only location
                        no_hg19_coordinates = False
                        position = ""
                        end_position = ""
                        for char in components[1][2:]:
                            if char in int_str:
                                position += char
                            else:
                                break
                        chromosome = 'chr'+str(int(components[0].split(".")[0][7:]))
                        if chromosome == "chr23":
                            chromosome = "chrX"
                        if chromosome == "chr24":
                            chromosome = "chrY"
                        ref, alt = ".", "."
                        if ">" in components[1]: # SNP and MNP
                            ref, alt = components[1][len(position) + 2:].split(">")
                        elif "dup" in components[1]: # Duplications 
                            end_position = components[1][len(position) + 3:-3]
                            position = int(position) - 1
                            if not end_position:
                                end_position = int(position) + 1
                            ref = get_dna_sequence(chromosome, int(position), int(end_position))
                            alt = ref[0] + ref[1:] * 2
                            if len(alt) == 1:
                                alt = alt*2
                        elif "delins" in components[1]: # Insertion/deletion combination
                            if "_" in components[1]:
                                end_position, alt = components[1][len(position) + 3:].split("delins")
                            else:
                                alt = components[1].split("delins")[1]
                                end_position = int(position) + 1
                            position = int(position) - 1
                            ref = get_dna_sequence(chromosome, position, int(end_position))
                            alt = ref[0] + alt
                        elif "del" in components[1]: # Deletions
                            end_position = components[1][len(position) + 3:-3]
                            if not end_position:
                                end_position = int(position) + 1
                            ref = get_dna_sequence(chromosome, int(position) -1, int(end_position) + 1)
                            alt = ref[0]
                        elif "ins" in components[1]: # Insertions
                            alt = components[1][len(position) + 3:].split("ins")[1]
                            ref = get_dna_sequence(chromosome, int(position), int(position) + 1)
                            alt = ref + alt
                        position = str(position)
                        vcf_data.append((chromosome, position, ref, alt, line_count))
                        break 
            if no_hg19_coordinates and verbose:
                print(columns[3], line_count)
    print(f"Wrote vcf entries for {len(vcf_data)} out of {line_count} total Erepo variants")
    return vcf_data


def main():
    """
    Converting the erepo tsv file into a DRAGEN formated vcf file
    """
    # Generate VCF entries from the erepo file
    erepo_filename = 'erepo.tabbed.txt'
    vcf_data = generate_vcf_data(erepo_filename)

    def generate_vcf_entry(chrom, pos, ref, alt, line_count):
        # Create a template for the VCF entry with placeholders for chrom, pos, ref, and alt
        vcf_template = "{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tDP=562;MQ=247.50;FractionInformativeReads=0.966;AQ=100.00" + \
                        ";GermlineStatus=Germline_DB;EreppoLine={line_count}\tGT:SQ:AD:AF:F1R2:F2R1:DP:SB:MB\t0/1:11.80:434" + \
                        ",109:0.201:209,56:225,53:543:226,208,51,58:217,217,48,61"

        # Format the template with the provided chrom, pos, ref, and alt values
        vcf_entry = vcf_template.format(chrom=chrom, pos=pos, ref=ref, alt=alt, line_count=line_count)
        
        return vcf_entry

    # Output VCF entries
    output_filename = 'output.vcf'
    with open(output_filename, 'w') as vcf_file:
        with open("vcf_header.txt", 'r') as in_file:
            for line in in_file:
                vcf_file.write(line)
            for entry in vcf_data:
                vcf_entry = generate_vcf_entry(entry[0], entry[1], entry[2], entry[3], entry[4])

                vcf_file.write(vcf_entry + '\n')

    print(f"VCF entries have been written to {output_filename}")

if __name__ == "__main__" :
    main()
