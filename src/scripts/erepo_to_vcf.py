"""
Converts Erepo ACMG curated variant information in tsv format to hg19 DRAGEN vcf format

Erepo
https://erepo.clinicalgenome.org/evrepo/

Requires a reference genome fasta (hg19.fa or hg38.fa)
"""
import argparse
from datetime import datetime
from pysam import FastaFile 

# GenBank and RefSeq accessions from 
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
ref_to_accessions = {
    'hg19': {
        "CM000663.2": "chr1", "NC_000001.11": "chr1",
        "CM000664.2": "chr2", "NC_000002.12": "chr2",
        "CM000665.2": "chr3", "NC_000003.12": "chr3",
        "CM000666.2": "chr4", "NC_000004.12": "chr4",
        "CM000667.2": "chr5", "NC_000005.10": "chr5",
        "CM000668.2": "chr6", "NC_000006.12": "chr6",
        "CM000669.2": "chr7", "NC_000007.14": "chr7",
        "CM000670.2": "chr8", "NC_000008.11": "chr8",
        "CM000671.2": "chr9", "NC_000009.12": "chr9",
        "CM000672.2": "chr10", "NC_000010.11": "chr10",
        "CM000673.2": "chr11", "NC_000011.10": "chr11",
        "CM000674.2": "chr12", "NC_000012.12": "chr12",
        "CM000675.2": "chr13", "NC_000013.11": "chr13",
        "CM000676.2": "chr14", "NC_000014.9": "chr14",
        "CM000677.2": "chr15", "NC_000015.10": "chr15",
        "CM000678.2": "chr16", "NC_000016.10": "chr16",
        "CM000679.2": "chr17", "NC_000017.11": "chr17",
        "CM000680.2": "chr18", "NC_000018.10": "chr18",
        "CM000681.2": "chr19", "NC_000019.10": "chr19",
        "CM000682.2": "chr20", "NC_000020.11": "chr20",
        "CM000683.2": "chr21", "NC_000021.9": "chr21",
        "CM000684.2": "chr22", "NC_000022.11": "chr22",
        "CM000685.2": "chrX", "NC_000023.11": "chrX",
        "CM000686.2": "chrY", "NC_000024.10": "chrY",
    },
    'hg38': {
        "CM000663.2": "chr1", "NC_000001.11": "chr1",
        "CM000664.2": "chr2", "NC_000002.12": "chr2",
        "CM000665.2": "chr3", "NC_000003.12": "chr3",
        "CM000666.2": "chr4", "NC_000004.12": "chr4",
        "CM000667.2": "chr5", "NC_000005.10": "chr5",
        "CM000668.2": "chr6", "NC_000006.12": "chr6",
        "CM000669.2": "chr7", "NC_000007.14": "chr7",
        "CM000670.2": "chr8", "NC_000008.11": "chr8",
        "CM000671.2": "chr9", "NC_000009.12": "chr9",
        "CM000672.2": "chr10", "NC_000010.11": "chr10",
        "CM000673.2": "chr11", "NC_000011.10": "chr11",
        "CM000674.2": "chr12", "NC_000012.12": "chr12",
        "CM000675.2": "chr13", "NC_000013.11": "chr13",
        "CM000676.2": "chr14", "NC_000014.9": "chr14",
        "CM000677.2": "chr15", "NC_000015.10": "chr15",
        "CM000678.2": "chr16", "NC_000016.10": "chr16",
        "CM000679.2": "chr17", "NC_000017.11": "chr17",
        "CM000680.2": "chr18", "NC_000018.10": "chr18",
        "CM000681.2": "chr19", "NC_000019.10": "chr19",
        "CM000682.2": "chr20", "NC_000020.11": "chr20",
        "CM000683.2": "chr21", "NC_000021.9": "chr21",
        "CM000684.2": "chr22", "NC_000022.11": "chr22",
        "CM000685.2": "chrX", "NC_000023.11": "chrX",
        "CM000686.2": "chrY", "NC_000024.10": "chrY",
        }
}

def get_dna_sequence(reference_genome, chrom, start, end):
    """
    Return an hg19 dna sequence
    """
    
    fasta = FastaFile(reference_genome)
    sequence = fasta.fetch(chrom, start, end)
    fasta.close()
    return sequence.upper()


def is_before(date1, date2):
    """
    Check if the first date is before the second date, handling both "YYYY-MM-DD" and "YYYY/M/D" formats.
    """
    formats = ["%Y-%m-%d", "%Y/%m/%d"]

    def parse_date(date_str):
        for fmt in formats:
            try:
                return datetime.strptime(date_str, fmt)
            except ValueError:
                continue
        raise ValueError(f"Invalid date format: {date_str}")

    return parse_date(date1) < parse_date(date2)


def extract_ref_alt_position_from_g_coords(gen_coords, chromosome, reference_genome, components):
    """
    Take in genomic coordinates  g.102917130T>C and determine the ref and alt.

    Sometimes we need additional reference bases to accurately make the vcf, we use the reference genome 
    to find them in these cases. 
    """
    snp_count, dup_count, delins_count, del_count, ins_count = 0, 0, 0, 0, 0
    int_str = {str(i) for i in range(10)}
    if "_" in gen_coords: 
        position = "".join(c for c in gen_coords[2:].split("_")[0] if c in int_str)
    else:
        position = "".join(c for c in gen_coords[2:] if c in int_str)
    ref, alt = ".", "."
    if ">" in gen_coords:  # SNP and MNP
        ref, alt = gen_coords[len(position) + 2:].split(">")
        snp_count += 1
    elif "dup" in gen_coords:  # Duplications
        end_position = gen_coords[len(position) + 3:-3] or str(int(position) + 1)
        position = int(position) - 1
        if int(end_position) - position > 1:
            ref = get_dna_sequence(reference_genome, chromosome, int(position) -1, int(end_position))
            alt = ref[0] + ref[1:] * 2
        else:
            ref = get_dna_sequence(reference_genome, chromosome, int(position), int(end_position))
            alt = ref[0]*2 + ref[1:]
        dup_count += 1
    elif "delins" in gen_coords:  # Insertion/Deletion combination
        if "_" in gen_coords:
            end_position, alt = components[1][len(position) + 3:].split("delins")
        else:
            end_position = str(int(position) + 1)
            alt = gen_coords.split("delins")[1]
        position = int(position) - 1
        ref = get_dna_sequence(reference_genome, chromosome, position, int(end_position))
        alt = ref[0] + alt
        delins_count += 1
    elif "del" in gen_coords:  # Deletions
        end_position = gen_coords[len(position) + 3:-3] or str(int(position) + 1)
        ref = get_dna_sequence(reference_genome, chromosome, int(position) - 1, int(end_position) + 1)
        alt = ref[0]
        del_count += 1
    elif "ins" in gen_coords:  # Insertions
        alt = gen_coords[len(position) + 3:].split("ins")[1]
        ref = get_dna_sequence(reference_genome, chromosome, int(position), int(position) + 1)
        alt = ref + alt
        ins_count += 1
    v_type = (snp_count, dup_count, delins_count, del_count, ins_count)
    return position, ref, alt, v_type

def generate_vcf_data(erepo_filename, reference_genome, cutoff_date, ref_b, verbose=False):
    """
    Generate VCF data from an Erepo TSV file.

    NC_000010.9:g.89680808G>A,
    NC_000010.9:g.89682800_89682802dup,
    NC_000010.9:g.89701896_89701897delinsAT,
    NC_000010.9:g.89715049_89715051del,
    NC_000017.9:g.75696123_75696124insAGCGGGC,
    """
    vcf_data = []
    snp_count, dup_count, delins_count, del_count, ins_count = 0, 0, 0, 0, 0
    no_assembly_comp, past_cutoff_date = 0, 0
    with open(erepo_filename, 'r') as erepo_file:
        line_count = 0
        erepo_file.readline()
        for line in erepo_file:
            line_count += 1
            columns = line.strip().split('\t')
            hgvs_expressions = columns[3].strip().split(",")
            publish_date = columns[-4]
            if not is_before(publish_date, cutoff_date):
                past_cutoff_date += 1
                continue
            no_assembly_coordinates = True
            for hgvs in hgvs_expressions: # ex NM_000277.2:c.1A>G or NC_000012.12:g.102917130T>C
                if "?" in hgvs:
                    continue  # Ignore strange ones
                components = hgvs.strip(",").split(":") # The separate elements, ex NC_000012.12 and g.102917130T>C
                if len(components) < 2: continue
                reference_assembly = components[0] # NOTE These are chromosome specific!
                gen_coords = components[1]
                if reference_assembly in ref_to_accessions[ref_b]: # Only looking for standard hg19 chrosomes
                    no_assembly_coordinates = False
                    # The genomic coordinats ex g.102917130T>C start with 'g.'
                    chromosome = ref_to_accessions[ref_b][reference_assembly]
                    try:
                        position, ref, alt, v_type = extract_ref_alt_position_from_g_coords(gen_coords, chromosome, reference_genome, components)
                    except ValueError: # Complex variants that are difficult to accurately translate
                        continue
                    snp_count += v_type[0]
                    dup_count += v_type[1]
                    delins_count += v_type[2]
                    del_count += v_type[3]
                    ins_count += v_type[4]
                    vcf_data.append((chromosome, str(position), ref, alt, line_count))
                    break # Only one entry per variant
            if no_assembly_coordinates:
                no_assembly_comp += 1
                if verbose:
                    print(line, line_count)
    print(f"Skipped {past_cutoff_date} variants that were submitted past the cutoff date {cutoff_date}")
    print(f"Skipped {no_assembly_comp} variants where {ref} genomic coordinates were not found")
    print(f"Converted {snp_count} snp, {dup_count} dup, {delins_count} del-ins, {del_count} del, {ins_count} ins")
    print(f"Wrote VCF entries for {len(vcf_data)} out of {line_count} total Erepo variants")
    return vcf_data

def generate_vcf_entry(chrom, pos, ref, alt, line_count):
    """
    Format a VCF entry with ideal, high-quality variant calling values.
    """
    return (
        f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t"
        f"DP=10000;MQ=60.00;FractionInformativeReads=1.000;AQ=100.00;"
        f"GermlineStatus=Germline_DB;EreppoLine={line_count}\t"
        f"GT:SQ:AD:AF:F1R2:F2R1:DP:SB:MB\t0/1:99.99:5000,5000:0.500:2500,2500:2500,2500:10000:5000,5000,5000,5000:5000,5000,5000,5000\n"
    )

def chrom_sort_key(variant):
    """
    Helper function for sorting the vcf
    """
    chrom_order = {str(i): i for i in range(1, 23)}
    chrom_order.update({"X": 23, "Y": 24, "M": 25})  # Assigning numeric order to non-numeric chromosomes
    chrom, pos = variant[0].lstrip("chr"), int(variant[1])
    return (chrom_order.get(chrom, 26), pos)  # Default to 26 if an unexpected chromosome appears


def main():
    """
    Converting the erepo tsv file into a DRAGEN formatted vcf file
    """
    parser = argparse.ArgumentParser(description="Convert Erepo ACMG variant information to DRAGEN VCF format.")
    parser.add_argument("-i", "--input", required=True, help="Path to the Erepo TSV file")
    parser.add_argument("-o", "--output", required=True, help="Path to the output VCF file")
    parser.add_argument("-r", "--reference_fasta", required=True, help="Path to the reference genome FASTA file (e.g., hg19.fa or hg38.fa)")
    parser.add_argument("-H", "--header", required=True, help="Path to the VCF header file")
    parser.add_argument("-R", "--ref_build", required= True, action="store", choices=['hg19', 'hg38'], help="Available ref builds")
    parser.add_argument("-cd", "--cutoff_date", action="store", help="Provide a cutoff date")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose mode")
    parser.set_defaults(cutoff_date = "2025-01-31")
    parser.set_defaults(ref = 'hg19')
    args = parser.parse_args()

    vcf_data = generate_vcf_data(args.input, args.reference_fasta, args.cutoff_date, args.ref_build, args.verbose)
    sorted_vcf_data = sorted(vcf_data, key=chrom_sort_key)
    with open(args.output, 'w') as vcf_file:
        with open(args.header, 'r') as header_file:
            vcf_file.writelines(header_file.readlines())
        for entry in sorted_vcf_data:
            vcf_file.write(generate_vcf_entry(entry[0], entry[1], entry[2], entry[3], entry[4]))

    print(f"VCF entries have been written to {args.output}")

if __name__ == "__main__":
    main()
