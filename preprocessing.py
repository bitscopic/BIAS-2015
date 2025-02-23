"""
Handles all preprocessing for hg19 and hg38.

Tries to be intelligent when it can, this may result in very dumb behavior though (machines). Feel free
to edit it for your own purposes. We are not responsible for any effects of your editting. 
"""
import os
import subprocess
import argparse
import json
import shutil
import logging
from src.preprocessing import (
        generate_pathogenic_aa_list,
        extract_from_avada_track,
        join_coding_and_repeats,
        find_missense_pathogenic_genes_and_path_trunc_genes,
        generate_domain_lists
)

def parse_args():
    """
    Parse command-line arguments into usable Python objects.
    """
    parser = argparse.ArgumentParser(description="Preprocess data for variant interpretation on hg19/GRCh37.")
    parser.add_argument("--output_dir",
                        required=True,
                        help="Output directory where results will be written.")
    parser.add_argument("--reference_build",
                        required=True,
                        help="The reference build to use.",
                        choices=['hg19','hg38'])
    parser.add_argument("--verbose",
                        help = " The verbosity level for stdout messages (default INFO)",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        action = "store")    
    parser.add_argument("--os_type",
                        required=True,
                        help="The operating system type.",
                        choices=['linux', 'macos'])
    parser.add_argument("--nirvana_bin_dir",
                        required=True,
                        help="Directory containing Nirvana binaries.")
    parser.add_argument("--nirvana_data_dir",
                        required=True,
                        help="Directory containing Nirvana data.")

    parser.set_defaults(verbose = "WARNING")
    options = parser.parse_args()
    logging.basicConfig(level=getattr(logging, options.verbose), format='%(message)s')
    return options

def run_command(command):
    """
    Runs a shell command. If the current logging level is INFO or DEBUG,
    prints the subprocess output; otherwise (WARNING, ERROR, or CRITICAL),
    hides it.
    """
    current_level = logging.getLogger().getEffectiveLevel()
    
    if current_level <= logging.INFO:
        # Print subprocess output when level is DEBUG or INFO
        subprocess.run(command, shell=True, check=True)
    else:
        # Hide subprocess output when level is WARNING, ERROR, or CRITICAL
        subprocess.run(
            command,
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

def gunzip_to_file(input_file, output_file):
    """
    gunzip to a file name
    """
    logging.debug("Gunzipping %s to %s", input_file, output_file)
    with open(output_file, "wb") as out_fh:
        subprocess.run(["gunzip", "-c", input_file], stdout=out_fh, check=True)


def bgzip_decompress_to_file(input_file, output_file):
    """
    bgzip to a file name
    """
    new_name = input_file.replace(".bgz", ".gz")
    logging.debug("Renaming %s to %s", input_file, new_name)
    os.rename(input_file, new_name)
    gunzip_to_file(new_name, output_file)


def download_and_extract(url, output_file):
    """
    Download a file from a URL and extract it if compressed, skipping if the file or its extracted form already exists.

    The URL should be the compressed file (.gz, .bgz) whereas the output file will be the uncompressed data
    """
    download_file = url.split("/")[-1]
    base_file = download_file.rsplit('.', 1)[0]
    if os.path.exists(download_file) or os.path.exists(base_file): # This is a weird corner case, most likely only hit if something goes wrong
        logging.info("File %s already exists. Skipping download", download_file)
        if url.endswith(".gz"):
            gunzip_to_file(download_file, output_file)
        elif url.endswith(".bgz"):
            bgzip_decompress_to_file(download_file, output_file)
        elif url.endswith(".bb"):
            shutil.copy(download_file, output_file)
        return
   
    logging.debug("Downloading %s", url)
    command = ["wget", url]
    run_command(" ".join(command))

    if url.endswith(".gz"):
        gunzip_to_file(download_file, output_file)
    elif url.endswith(".bgz"):
        bgzip_decompress_to_file(download_file, output_file)
    elif url.endswith(".bb"):
        shutil.copy(download_file, output_file)


def prepend_chr_to_vcf(input_vcf, output_vcf):
    """
    Add 'chr' prefix to chromosome names in a VCF file.
    """
    logging.debug("prepend_chr %s %s", input_vcf, output_vcf)
    with open(input_vcf, "r") as infile, open(output_vcf, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
            else:
                fields = line.strip().split("\t")
                fields[0] = f"chr{fields[0]}"
                outfile.write("\t".join(fields) + "\n")

def filter_vcf(input_vcf, output_vcf):
    """
    Filter VCF file for pathogenic and likely pathogenic variants.
    """ 
    logging.debug("filter_vcf %s %s", input_vcf, output_vcf)
    with open(input_vcf, "r") as infile, open(output_vcf, "w") as outfile:
        for line in infile:
            if line.startswith("#") or "CLNSIG=Likely_pathogenic" in line or "CLNSIG=Pathogenic" in line or 'CLNGSIG=drug_response' in line:
                outfile.write(line)

def generate_output_json(paths, output_json):
    """
    Generate a JSON file listing all output file paths.
    """
    with open(output_json, "w") as json_file:
        json.dump(paths, json_file, indent=4)

def process_ncbi_ref_seq_hgmd(ncbi_file, ref_seq_download_file, ref_seq_file):
    """
    Process ncbiRefSeqHgmd data.
    """
    download_and_extract(ref_seq_download_file, ref_seq_file)
    logging.debug("cp %s %s", ref_seq_file, ncbi_file)
    shutil.copy(ref_seq_file, ncbi_file)

def annotate_clinvar_with_nirvana(nirvana_bin_dir, nirvana_data_dir, in_file, out_name, ref_b):
    """
    Annotate ClinVar VCF with Nirvana 
    """
    if ref_b == 'hg19':
        ref_b = 'GRCh37'
    else:
        ref_b = 'GRCh38'
    command = f"dotnet {nirvana_bin_dir}/Nirvana.dll --cache {nirvana_data_dir}/Cache/{ref_b}/Both --sd {nirvana_data_dir}/SupplementaryAnnotation/{ref_b} "+\
            f"--ref {nirvana_data_dir}/References/Homo_sapiens.{ref_b}.Nirvana.dat --in {in_file} --o {out_name}"
    logging.debug("%s", command)
    run_command(command)


def process_clinvar_vcf(clinvar_vcf_download_file, clinvar_vcf, clinvar_with_chr, clean_clinvar, aa_output,
                        nirvana_bin_dir, nirvana_data_dir, nirvana_output_tag, ref_b):
    """
    Process ClinVar VCF file. Isolate the variants which are pathogenic and likely pathogenic, then
    annotate them with NIRVANA/ICA and extract the amino acids changes identified. Store the
    gene:amino acid information for use in PS1
    """
    download_and_extract(clinvar_vcf_download_file, clinvar_vcf)
    prepend_chr_to_vcf(clinvar_vcf, clinvar_with_chr)
    filter_vcf(clinvar_with_chr, clean_clinvar)
    annotate_clinvar_with_nirvana(nirvana_bin_dir, nirvana_data_dir, clean_clinvar, nirvana_output_tag, ref_b)
    logging.debug("generate_pathogenic_aa_list %s %s %s", f"{nirvana_output_tag}.json.gz", clean_clinvar, aa_output)
    generate_pathogenic_aa_list.extract_aa_information(f"{nirvana_output_tag}.json.gz", clean_clinvar, aa_output)

def get_big_bed_executable(bigbed_download_file, bigbed_executable):
    """
    BigBed files are a UCSC compressed bed format. They can be converted to normal bed files
    using the UCSC executable bigBedToBed which this downloads and uses.

    BigBed files have a .bb extension.
    """
    download_and_extract(bigbed_download_file, bigbed_executable)
    run_command(f"chmod +x {bigbed_executable}")


def process_absplice(absplice_download_bb, output_dir, absplice_bb, absplice_bed, splice_output):
    """
    Process the absplice big bed file. PVS1, PP3, BP4, BP7 

    NOTE: The full absplice file is ~60 gb when extracted. To avoid this, we unpack one chromosome at a
    time, extract the information BIAS needs, merge outputs, and immediately delete intermediate files.
    """
    download_and_extract(absplice_download_bb, absplice_bb)

    # Create empty merged output files
    with open(splice_output, 'w') as o_file:
        for chromosome in [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]:
            # Modify filenames to be chromosome-specific
            chrom_absplice_bed = f"{absplice_bed}_{chromosome}"

            # Convert BigBed to Bed for the specific chromosome
            run_command(f"{output_dir}/bigBedToBed -chrom={chromosome} {absplice_bb} {chrom_absplice_bed}")

            # Extract essential data from the bed file
            with open(chrom_absplice_bed, 'r') as in_file:
                for line in in_file:
                    split_line = line.strip().split("\t")
                    chrom = split_line[0]
                    pos = split_line[1]
                    mut = split_line[3]
                    splice_score = split_line[11]
                    o_file.write(f"{chrom}\t{pos}\t{mut}\t{splice_score}\n")

            # Remove the chrom level splice bed file to save disk space
            os.remove(chrom_absplice_bed)


def process_avada_big_bed_file(avada_download_bb, output_dir, avada_bb, avada_bed, lit_gene_mut_output, lit_variant_output):
    """
    Process avada big bed file. Standardize gene:mut and variant coding (chrom, pos, ref, alt) definitions. Used in PS3
    """
    download_and_extract(avada_download_bb, avada_bb)
    run_command(f"{output_dir}/bigBedToBed {avada_bb} {avada_bed}")
    extract_from_avada_track.clean_and_extract_columns(avada_bed, lit_gene_mut_output, lit_variant_output)
    

def process_gwas_catalog(gwas_download_file, gwas_int_file, gwas_file):
    """
    Process GWAS catalog. PS4
    """
    download_and_extract(gwas_download_file, gwas_int_file) 
    logging.debug("cp %s %s", gwas_int_file, gwas_file)
    shutil.copy(gwas_int_file, gwas_file)

def process_pathogenic_domains(unip_domain_download_file, unip_bb, unip_bed, clinvar_vcf, pathogenic_output, output_dir):
    """
    PM1 - Generate lists of pathogenic domains using ClinVar data.

    """
    # Download and extract the UniProt domain file
    download_and_extract(unip_domain_download_file, unip_bb)
    
    # Convert BigBed to BED
    logging.debug("Converting BigBed to BED: %s -> %s", unip_bb, unip_bed)
    run_command(f"{output_dir}/bigBedToBed {unip_bb} {unip_bed}")
    
    # Generate ClinVar-based pathogenic domains
    logging.debug("generate_domain_lists %s %s %s", clinvar_vcf, unip_bed, pathogenic_output)
    generate_domain_lists.evaluate_clinvar_domains(clinvar_vcf, unip_bed, pathogenic_output)


def process_repeat_and_coding_gene_info(rmsk_download_file, ccds_download_file, rmsk_file, ccds_file, repeat_output):
    """
    Process repeat and coding gene information.
    """
    download_and_extract(rmsk_download_file, rmsk_file)
    download_and_extract(ccds_download_file, ccds_file)  
    logging.debug("join_coding_and_repeats %s %s %s", rmsk_file, ccds_file, repeat_output)
    join_coding_and_repeats.join_coding_and_repeats(rmsk_file, ccds_file, repeat_output)


def process_missense_pathogenic_genes_and_truncating_genes(full_clinvar, gnomad_rmc_download_bb, gnomad_rmc_bb, 
                                                           gnomad_rmc_bed, missense_output, truncating_output, output_dir, ref_b):
    """
    Find missense pathogenic genes.
    """
    if ref_b == 'hg19':
        download_and_extract(gnomad_rmc_download_bb, gnomad_rmc_bb)
        run_command(f"{output_dir}/bigBedToBed {gnomad_rmc_bb} {gnomad_rmc_bed}")
    logging.debug("find_missense_pathogenic_genes_and_path_trunc_genes %s %s %s %s", full_clinvar, gnomad_rmc_bed, missense_output, truncating_output)
    find_missense_pathogenic_genes_and_path_trunc_genes.find_missense_pathogenic_genes_and_path_trunc_genes(full_clinvar, gnomad_rmc_bed, missense_output, truncating_output, ref_b)


def main():
    """
    Main function to orchestrate preprocessing steps.
    """
    options = parse_args()
    output_dir = os.path.abspath(options.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Define reference build
    ref_b = options.reference_build
    
    # Definte the ICA/Nirvana path and the data path
    nirvana_bin_dir = options.nirvana_bin_dir
    nirvana_data_dir = options.nirvana_data_dir    
    
    # This is an executable used to convert bigBed (bb) files to bed files
    logging.basicConfig(level=getattr(logging, options.verbose), format='%(message)s')
    if options.os_type == 'linux':
        bigbed_download_file = "https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed"
    elif options.os_type == 'macos':
        bigbed_download_file = "https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/bigBedToBed"
    bigbed_executable = os.path.join(output_dir, "bigBedToBed")
    get_big_bed_executable(bigbed_download_file, bigbed_executable)
     
    # Download links
    if ref_b == 'hg19':
        ref_seq_download_file = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ncbiRefSeqHgmd.txt.gz"
        clinvar_vcf_download_file = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz"
        rmsk_download_file = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz"
        ccds_download_file = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ccdsGene.txt.gz"
        gwas_download_file = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gwasCatalog.txt.gz"    
        absplice_download_bb = "https://hgdownload.soe.ucsc.edu/gbdb/hg19/abSplice/AbSplice.bb"
        avada_download_bb = "http://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/avada.bb"
        unip_domain_download_file = "https://hgdownload.soe.ucsc.edu/gbdb/hg19/uniprot/unipDomain.bb"
        gnomad_rmc_download_bb = "https://hgdownload.soe.ucsc.edu/gbdb/hg19/gnomAD/missense/missenseConstrained.bb"
    elif ref_b == 'hg38':
        ref_seq_download_file = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqHgmd.txt.gz"
        clinvar_vcf_download_file = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
        rmsk_download_file = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
        ccds_download_file = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ccdsGene.txt.gz"
        gwas_download_file = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gwasCatalog.txt.gz"    
        absplice_download_bb = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/abSplice/AbSplice.bb"
        avada_download_bb = "http://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/avada.bb"
        unip_domain_download_file = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/uniprot/unipDomain.bb"
        gnomad_rmc_download_bb = "" # UCSC has not yet made this track for hg38

    # Intermediate files & names
    clinvar_vcf = os.path.join(output_dir, f"{ref_b}_clinvar.vcf")
    clinvar_with_chr = os.path.join(output_dir, f"{ref_b}_clinvar_with_chr.vcf")
    clean_clinvar = os.path.join(output_dir, f"{ref_b}_clean_clinvar.vcf")
    absplice_bb = os.path.join(output_dir, f"{ref_b}_absplice.bb")
    absplice_bed = os.path.join(output_dir, f"{ref_b}_absplice.bed")
    avada_bb = os.path.join(output_dir, f"{ref_b}_avada.bb")
    avada_bed = os.path.join(output_dir, f"{ref_b}_avada.bed")
    gnomad_rmc_bb = os.path.join(output_dir, f"{ref_b}_gnomad_rmc.bb")
    gnomad_rmc_bed = os.path.join(output_dir, f"{ref_b}_gnomad_rmc.bed")
    rmsk_file = os.path.join(output_dir, f"{ref_b}_rmsk.txt")
    ccds_file = os.path.join(output_dir, f"{ref_b}_ccdsGene.txt")
    unip_bb = os.path.join(output_dir, f"{ref_b}_unipDomain.bb")
    unip_bed = os.path.join(output_dir, f"{ref_b}_unipDomain.bed")
    ref_seq_file = os.path.join(output_dir, f"{ref_b}_ncbiRefSeqHgmd.txt")
    nirvana_output_tag = os.path.join(output_dir, f'{ref_b}_clean_clinvar_nirvana')
    gwas_int_file = os.path.join(output_dir, f"{ref_b}_gwasCatalog.txt")

    # Output files
    ncbi_file = os.path.join(output_dir, f"{ref_b}_PVS1_ncbiRefSeqHgmd.tsv")
    splice_output = os.path.join(output_dir, f"{ref_b}_PVS1_PP3_BP4_BP7_splice_data.tsv")
    aa_output = os.path.join(output_dir, f"{ref_b}_PS1_PM5_clinvar_pathogenic_aa.tsv")
    lit_gene_aa_output = os.path.join(output_dir, f"{ref_b}_PS3_lit_gene_aa.tsv")
    lit_variant_output = os.path.join(output_dir, f"{ref_b}_PS3_lit_variant.tsv")
    gwas_file = os.path.join(output_dir, f"{ref_b}_PS4_gwasCatalog.txt")
    pathogenic_domain_output = os.path.join(output_dir, f"{ref_b}_PM1_chrom_to_pathogenic_domain_list.tsv")
    repeat_output = os.path.join(output_dir, f"{ref_b}_PM4_BP3_coding_repeat_regions.tsv")
    missense_output = os.path.join(output_dir, f"{ref_b}_PP2_missense_pathogenic_genes.tsv")
    truncating_output = os.path.join(output_dir, f"{ref_b}_BP1_truncating_genes.tsv")
    
    
    logging.info("Processing NCBI RefSeq HGMD data... (PVS1 caveat 1)")
    process_ncbi_ref_seq_hgmd(ncbi_file, ref_seq_download_file, ref_seq_file) # PVS1 caveat 1
   
    logging.info("Processing ABSplice data... (PVS1 & PP3 & BP4 & BP7)") 
    process_absplice(absplice_download_bb, output_dir, absplice_bb, absplice_bed, splice_output) # Multiple codes

    logging.info("Processing ClinVar VCF file... (PS1)")
    process_clinvar_vcf(clinvar_vcf_download_file, clinvar_vcf, clinvar_with_chr, clean_clinvar, aa_output,
                        nirvana_bin_dir, nirvana_data_dir, nirvana_output_tag, ref_b) # PS1
    
    logging.info("Processing avada big bed file... (PS3)")
    process_avada_big_bed_file(avada_download_bb, output_dir, avada_bb, avada_bed, lit_gene_aa_output, lit_variant_output) # PS3
    
    logging.info("Processing GWAS catalog... (PS4)")
    process_gwas_catalog(gwas_download_file, gwas_int_file, gwas_file) #PS4
    
    logging.info("Processing pathogenic domains... (PM1)")
    process_pathogenic_domains(unip_domain_download_file, unip_bb, unip_bed, clinvar_with_chr, pathogenic_domain_output, output_dir) # PM1
    
    logging.info("Processing repeat and coding gene information... (PM4 & BP1)")
    process_repeat_and_coding_gene_info(rmsk_download_file, ccds_download_file, rmsk_file, ccds_file, repeat_output) # PM4 & BP1
    
    logging.info("Processing missense pathogenic genes and truncating genes... (PP2 & BP1)")
    process_missense_pathogenic_genes_and_truncating_genes(clinvar_with_chr, gnomad_rmc_download_bb, gnomad_rmc_bb,
                                                           gnomad_rmc_bed, missense_output, truncating_output, output_dir, ref_b) # PP2 & BP1

    # Write an output required_paths file
    output_json = os.path.join(output_dir, f"{ref_b}_required_paths.json")
    generate_output_json({
        "PVS1_ncbi_ref_seq_hgmd_fp": ncbi_file,
        "PVS1_PP3_BP4_BP7_splice_fp": splice_output,
        "PS1_PM5_clinvar_pathogenic_aa_fp": aa_output,
        "PS3_literature_gene_aa_fp": lit_gene_aa_output,
        "PS3_literature_variant_fp": lit_variant_output,
        "PS4_gwas_dbsnp_fp": gwas_file,
        "PM1_chrom_to_pathogenic_domain_list_fp": pathogenic_domain_output,
        "PM4_BP3_coding_repeat_region_fp": repeat_output,
        "PP2_missense_pathogenic_gene_to_region_list_fp": missense_output,
        "BP1_truncating_gene_to_data_fp": truncating_output,
    }, output_json)

    print(f"Setup complete. Output written to {output_json}")

if __name__ == "__main__":
    main()
