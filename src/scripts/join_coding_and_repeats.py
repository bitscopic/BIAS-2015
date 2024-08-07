"""
Overlap consensus coding sequences with repeat regions

    Repeat Masker:
    Smit AFA, Hubley R, Green P. RepeatMasker Open-3.0. http://www.repeatmasker.org. 1996-2010.

    Repbase Update is described in:

    Jurka J. Repbase Update: a database and an electronic journal of repetitive elements. Trends Genet.
        2000 Sep;16(9):418-420. PMID: 10973072

    For a discussion of repeats in mammalian genomes, see:

    Smit AF. Interspersed repeats and other mementos of transposable elements in mammalian genomes. Curr Opin Genet
        Dev. 1999 Dec;9(6):657-63. PMID: 10607616

    Smit AF. The origin of interspersed repeats in the human genome. Curr Opin Genet Dev. 1996 Dec;6(6):743-8.
        PMID: 8994846

    Consensus Coding:
    Hubbard T, Barker D, Birney E, Cameron G, Chen Y, Clark L, Cox T, Cuff J, Curwen V, Down T et al. The Ensembl
        genome database project. Nucleic Acids Res. 2002 Jan 1;30(1):38-41. PMID: 11752248; PMC: PMC99161

    Pruitt KD, Harrow J, Harte RA, Wallin C, Diekhans M, Maglott DR, Searle S, Farrell CM, Loveland JE, Ruef BJ et al. 
        The consensus coding sequence (CCDS) project: Identifying a common protein-coding gene set for the human and
        mouse genomes. Genome Res. 2009 Jul;19(7):1316-23. PMID: 19498102; PMC: PMC2704439

    Pruitt KD, Tatusova T, Maglott DR. NCBI Reference Sequence (RefSeq): a curated non-redundant sequence database of
        genomes, transcripts and proteins. Nucleic Acids Res. 2005 Jan 1;33(Database issue):D501-4. PMID: 15608248;
        PMC: PMC539979

"""
import argparse
import csv
from intervaltree import IntervalTree


def parseArgs(): 
    """
    Parse the command line arguments into useful python objects.  '--' variables are optional
    set_defaults only applies when the argument is not provided (it won't override)
    """
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument("--inRepeatFile",
                        help = " The input repeat region file",
                        action = "store")
    parser.add_argument("--inConsensusCodingFile",
                        help = " The input consensus coding file",
                        action = "store")
    parser.add_argument("--outMergedFile",
                        help = " The output file",
                        action = "store")
    parser.add_argument("--verbose",
                        help = " The verbosity level for stdout messages (default INFO)",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        action = "store")

    parser.set_defaults(inRepeatFile = "INFO")
    parser.set_defaults(inConsensusCodingFile = "INFO")
    parser.set_defaults(outMergedFile = "INFO")
    parser.set_defaults(verbose = "INFO")

    options = parser.parse_args()
    return options

# Step 1: Parse the files
def parse_file(filepath, chrom_col, start_col, end_col, strand_col=None):
    """
    Load files into a list of tuples
    """
    intervals = {}
    with open(filepath, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            chrom = row[chrom_col]
            start = int(row[start_col])
            end = int(row[end_col])
            strand = row[strand_col] if strand_col else None
            if chrom not in intervals:
                intervals[chrom] = []
            intervals[chrom].append((start, end, strand))
    return intervals

def main():
    """
    Load the two files then find any repeat regions that overlap with the consensus coding regions
    """
    options = parseArgs()
    
    repeat_regions_intervals = parse_file(options.inRepeatFile, 'genoName', 'genoStart', 'genoEnd', 'strand')
    consensus_coding_intervals = parse_file(options.inConsensusCodingFile, 'chrom', 'txStart', 'txEnd', 'strand')

    # Step 2: Build an interval tree for the second file
    interval_trees = {}
    for chrom, intervals in consensus_coding_intervals.items():
        tree = IntervalTree()
        for start, end, _ in intervals:
            tree[start:end] = None
        interval_trees[chrom] = tree

    # Step 3: Query overlaps
    overlapping_intervals = []
    for chrom, intervals in repeat_regions_intervals.items():
        if chrom in interval_trees:
            tree = interval_trees[chrom]
            for start, end, strand in intervals:
                overlaps = tree[start:end]
                if overlaps:
                    overlapping_intervals.append((chrom, start, end, strand))

    # Step 4: Output the results
    with open(options.outMergedFile, 'w') as o_file:
        for interval in overlapping_intervals:
            o_file.write(f'{interval[0]}\t{interval[1]}\t{interval[2]}\t{interval[3]}\n')

if __name__ == "__main__" :
    main()
