"""
Handles all preprocessing for hg19 and hg38.

Supports both Nirvana and VEP as the variant annotator. Use --annotator to select.

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
        generate_domain_lists,
        generate_vcep_af_cutoffs,
)
from src.preprocessing.generate_clinvar_submitter_counts import generate_submitter_counts

# Supported annotators
SUPPORTED_ANNOTATORS = ['nirvana', 'vep', 'all']


def skip_if_exists(output_files, description):
    """
    Check if output file(s) exist and log skip message if so.

    Args:
        output_files: Single file path (str) or list of file paths to check
        description: Description of the processing step for logging

    Returns:
        bool: True if ALL output files exist (should skip), False otherwise
    """
    if isinstance(output_files, str):
        output_files = [output_files]

    all_exist = all(os.path.exists(f) for f in output_files)
    if all_exist:
        logging.info("Skipping %s - output already exists: %s", description, ", ".join(output_files))
        return True
    return False

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
    parser.add_argument("--annotator",
                        help="The annotator to use for VCF annotation (default: nirvana)",
                        choices=SUPPORTED_ANNOTATORS,
                        default="nirvana")
    parser.add_argument("--verbose",
                        help=" The verbosity level for stdout messages (default INFO)",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        action="store")
    parser.add_argument("--os_type",
                        required=True,
                        help="The operating system type.",
                        choices=['linux', 'macos'])

    # Nirvana-specific arguments
    parser.add_argument("--nirvana_bin_dir",
                        help="Directory containing Nirvana binaries (required for --annotator nirvana).")
    parser.add_argument("--nirvana_data_dir",
                        help="Directory containing Nirvana data (required for --annotator nirvana).")

    # VEP-specific arguments
    parser.add_argument("--vep_cache_dir",
                        help="Directory containing VEP cache (homo_sapiens/). Required for --annotator vep or all.")

    parser.set_defaults(verbose="WARNING")
    options = parser.parse_args()
    logging.basicConfig(level=getattr(logging, options.verbose), format='%(message)s')

    # Validate annotator-specific arguments
    if options.annotator in ('nirvana', 'all'):
        if not options.nirvana_bin_dir or not options.nirvana_data_dir:
            parser.error("--nirvana_bin_dir and --nirvana_data_dir are required when using --annotator nirvana or all")
    if options.annotator in ('vep', 'all'):
        if not options.vep_cache_dir:
            parser.error("--vep_cache_dir is required when using --annotator vep or all")

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
    # Derive download filename from output_file to avoid collisions when
    # different builds share the same URL filename (e.g. both hg19 and hg38
    # ClinVar VCFs are named "clinvar.vcf.gz" at the source).
    url_filename = url.split("/")[-1]
    url_ext = ""
    if url_filename.endswith(".gz"):
        url_ext = ".gz"
    elif url_filename.endswith(".bgz"):
        url_ext = ".bgz"
    download_file = os.path.basename(output_file) + url_ext if url_ext else os.path.basename(output_file)
    base_file = download_file.rsplit('.', 1)[0]
    if os.path.exists(download_file):
        # The compressed/raw download file exists — decompress it
        logging.info("File %s already exists. Skipping download", download_file)
        if url.endswith(".gz"):
            gunzip_to_file(download_file, output_file)
        elif url.endswith(".bgz"):
            bgzip_decompress_to_file(download_file, output_file)
        else:
            shutil.copy(download_file, output_file)
        return
    if os.path.exists(base_file):
        # The extracted file already exists (from a previous run) — just copy it
        logging.info("Extracted file %s already exists. Skipping download", base_file)
        shutil.copy(base_file, output_file)
        return

    logging.debug("Downloading %s as %s", url, download_file)
    command = ["wget", "-O", download_file, url]
    run_command(" ".join(command))

    if url.endswith(".gz"):
        gunzip_to_file(download_file, output_file)
    elif url.endswith(".bgz"):
        bgzip_decompress_to_file(download_file, output_file)
    else:
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


def annotate_clinvar_with_vep(in_file, out_name, ref_b, output_dir, vep_cache_dir):
    """
    Annotate ClinVar VCF with VEP via Docker for preprocessing only.

    Only gene name, HGVSp (protein notation), and consequence are needed for PS1/PM5
    AA extraction, so we use a minimal VEP invocation with just --hgvs. No gnomAD or
    other custom annotations are required.

    Args:
        in_file (str): Input VCF file path (must be in output_dir)
        out_name (str): Output file path prefix (produces out_name.json)
        ref_b (str): Reference build ('hg19' or 'hg38')
        output_dir (str): Working directory (mounted as /work)
        vep_cache_dir (str): Host directory containing VEP cache (mounted as /cache)
    """
    if ref_b == 'hg19':
        assembly = 'GRCh37'
    elif ref_b == 'hg38':
        assembly = 'GRCh38'

    # Get relative paths for Docker
    in_filename = os.path.basename(in_file)
    out_filename = os.path.basename(out_name)

    command = (
        f'docker run --rm '
        f'-v "{output_dir}:/work" '
        f'-v "{vep_cache_dir}:/cache" '
        f'ensemblorg/ensembl-vep vep '
        f'--input_file /work/{in_filename} '
        f'--output_file /work/{out_filename}.json '
        f'--json --cache --offline '
        f'--assembly {assembly} '
        f'--dir_cache /cache '
        f'--hgvs --symbol --biotype '
        f'--force_overwrite'
    )
    logging.debug("%s", command)
    run_command(command)


def process_clinvar_aa_extraction(clean_clinvar, aa_output, annotator, annotator_output_tag, ref_b,
                                  nirvana_bin_dir=None, nirvana_data_dir=None, output_dir=None,
                                  vep_cache_dir=None):
    """
    Annotate a filtered ClinVar VCF with the specified annotator and extract amino acid changes for PS1/PM5.

    Args:
        clean_clinvar (str): Path to filtered (pathogenic/likely pathogenic) ClinVar VCF
        aa_output (str): Path for amino acid output file
        annotator (str): 'nirvana' or 'vep'
        annotator_output_tag (str): Output tag/path for annotator output
        ref_b (str): Reference build ('hg19' or 'hg38')
        nirvana_bin_dir (str): Nirvana binary directory (required for nirvana)
        nirvana_data_dir (str): Nirvana data directory (required for nirvana)
        output_dir (str): Working directory (required for vep)
    """
    if annotator == 'nirvana':
        annotated_json = f"{annotator_output_tag}.json.gz"
        if not os.path.exists(annotated_json):
            annotate_clinvar_with_nirvana(nirvana_bin_dir, nirvana_data_dir, clean_clinvar, annotator_output_tag, ref_b)
        else:
            logging.info("Skipping Nirvana annotation - file already exists: %s", annotated_json)
        logging.debug("generate_pathogenic_aa_list %s %s %s", annotated_json, clean_clinvar, aa_output)
        generate_pathogenic_aa_list.extract_aa_information(annotated_json, clean_clinvar, aa_output)
    elif annotator == 'vep':
        annotated_json = f"{annotator_output_tag}.json"
        if not os.path.exists(annotated_json):
            annotate_clinvar_with_vep(clean_clinvar, annotator_output_tag, ref_b, output_dir, vep_cache_dir)
        else:
            logging.info("Skipping VEP annotation - file already exists: %s", annotated_json)
        logging.debug("generate_pathogenic_aa_list (VEP) %s %s %s", annotated_json, clean_clinvar, aa_output)
        generate_pathogenic_aa_list.extract_aa_information_vep(annotated_json, clean_clinvar, aa_output)


def get_big_bed_executable(bigbed_download_file, bigbed_executable):
    """
    BigBed files are a UCSC compressed bed format. They can be converted to normal bed files
    using the UCSC executable bigBedToBed which this downloads and uses.

    BigBed files have a .bb extension.
    """
    download_and_extract(bigbed_download_file, bigbed_executable)
    os.chmod(bigbed_executable, 0o755)


def process_absplice(absplice_download_bb, output_dir, absplice_bb, absplice_bed, splice_output):
    """
    Process the absplice big bed file. PVS1, PP3, BP4, BP7

    NOTE: The full absplice file is ~60 gb when extracted. To avoid this, we unpack one chromosome at a
    time, extract the information BIAS needs, merge outputs, and immediately delete intermediate files.
    """
    download_and_extract(absplice_download_bb, absplice_bb)

    # Create empty merged output files
    with open(splice_output, 'w') as o_file:
        for chromosome in [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]:
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

    # Define reference build and annotator
    ref_b = options.reference_build
    annotator = options.annotator

    # Determine which annotators to run
    if annotator == 'all':
        annotators_to_run = ['nirvana', 'vep']
    else:
        annotators_to_run = [annotator]

    # Annotator-specific configuration
    nirvana_bin_dir = options.nirvana_bin_dir
    nirvana_data_dir = options.nirvana_data_dir
    vep_cache_dir = options.vep_cache_dir

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
    gwas_int_file = os.path.join(output_dir, f"{ref_b}_gwasCatalog.txt")

    # Gene-level data URLs (reference-build independent)
    clingen_download_url = "https://search.clinicalgenome.org/kb/gene-validity/download"
    gnomad_constraints_download_url = "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"

    # ClinVar submission-level data URLs (VEP-only, reference-build independent)
    clinvar_submission_summary_url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz"
    clinvar_variant_summary_url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"

    # Gene-level data output files
    clingen_file = os.path.join(output_dir, f"{ref_b}_clingen_gene_disease_validity.csv")
    gnomad_constraints_file = os.path.join(output_dir, f"{ref_b}_gnomad_gene_constraints.txt")

    # Output files
    ncbi_file = os.path.join(output_dir, f"{ref_b}_PVS1_ncbiRefSeqHgmd.tsv")
    splice_output = os.path.join(output_dir, f"{ref_b}_PVS1_PP3_BP4_BP7_splice_data.tsv")
    lit_gene_aa_output = os.path.join(output_dir, f"{ref_b}_PS3_lit_gene_aa.tsv")
    lit_variant_output = os.path.join(output_dir, f"{ref_b}_PS3_lit_variant.tsv")
    gwas_file = os.path.join(output_dir, f"{ref_b}_PS4_gwasCatalog.txt")
    pathogenic_domain_output = os.path.join(output_dir, f"{ref_b}_PM1_chrom_to_pathogenic_domain_list.tsv")
    repeat_output = os.path.join(output_dir, f"{ref_b}_PM4_BP3_coding_repeat_regions.tsv")
    missense_output = os.path.join(output_dir, f"{ref_b}_PP2_missense_pathogenic_genes.tsv")
    truncating_output = os.path.join(output_dir, f"{ref_b}_BP1_truncating_genes.tsv")
    # VCEP AF cutoffs are gene-only (build-agnostic) and committed to the repo at
    # data/vcep/. Rerunning preprocessing refreshes the committed snapshot in place.
    _repo_root = os.path.abspath(os.path.dirname(__file__))
    vcep_af_output = os.path.join(_repo_root, "data", "vcep", "vcep_af_cutoffs.tsv")
    vcep_af_review = os.path.join(_repo_root, "data", "vcep", "vcep_af_cutoffs_review.tsv")
    # Derived m o/e × MOI binned AF cutoffs and the MOI lookup are also gene-only
    # (build-agnostic) and ship committed under data/af_cutoffs/.
    mis_oe_moi_af_cutoffs_file = os.path.join(_repo_root, "data", "af_cutoffs", "mis_oe_upper_binned_cutoffs.tsv")
    gene_to_moi_file = os.path.join(_repo_root, "data", "af_cutoffs", "gene_to_moi.tsv")

    # VEP-only: ClinVar submission-level data for PS4 submitter counts
    clinvar_submission_summary_gz = os.path.join(output_dir, "submission_summary.txt.gz")
    clinvar_variant_summary_gz = os.path.join(output_dir, "variant_summary.txt.gz")
    clinvar_submitter_counts_output = os.path.join(output_dir, f"{ref_b}_PS4_clinvar_submitter_counts.tsv")

    # ── Prepare ClinVar VCF (shared by PM1, PP2, BP1, and annotator-specific PS1/PM5) ──
    if not os.path.exists(clinvar_vcf):
        download_and_extract(clinvar_vcf_download_file, clinvar_vcf)
    else:
        logging.info("Skipping download - ClinVar VCF already exists: %s", clinvar_vcf)
    if not os.path.exists(clinvar_with_chr):
        prepend_chr_to_vcf(clinvar_vcf, clinvar_with_chr)
    else:
        logging.info("Skipping chr prepend - file already exists: %s", clinvar_with_chr)
    if not os.path.exists(clean_clinvar):
        filter_vcf(clinvar_with_chr, clean_clinvar)
    else:
        logging.info("Skipping VCF filter - file already exists: %s", clean_clinvar)

    # Process NCBI RefSeq HGMD data (PVS1 caveat 1)
    if not skip_if_exists(ncbi_file, "NCBI RefSeq HGMD (PVS1)"):
        logging.info("Processing NCBI RefSeq HGMD data... (PVS1 caveat 1)")
        process_ncbi_ref_seq_hgmd(ncbi_file, ref_seq_download_file, ref_seq_file)

    # Process ABSplice data (PVS1 & PP3 & BP4 & BP7)
    if not skip_if_exists(splice_output, "ABSplice (PVS1 & PP3 & BP4 & BP7)"):
        logging.info("Processing ABSplice data... (PVS1 & PP3 & BP4 & BP7)")
        process_absplice(absplice_download_bb, output_dir, absplice_bb, absplice_bed, splice_output)

    # Process avada big bed file (PS3)
    if not skip_if_exists([lit_gene_aa_output, lit_variant_output], "Avada (PS3)"):
        logging.info("Processing avada big bed file... (PS3)")
        process_avada_big_bed_file(avada_download_bb, output_dir, avada_bb, avada_bed, lit_gene_aa_output, lit_variant_output)

    # Process GWAS catalog (PS4)
    if not skip_if_exists(gwas_file, "GWAS Catalog (PS4)"):
        logging.info("Processing GWAS catalog... (PS4)")
        process_gwas_catalog(gwas_download_file, gwas_int_file, gwas_file)

    # Process pathogenic domains (PM1)
    if not skip_if_exists(pathogenic_domain_output, "Pathogenic Domains (PM1)"):
        logging.info("Processing pathogenic domains... (PM1)")
        process_pathogenic_domains(unip_domain_download_file, unip_bb, unip_bed, clinvar_with_chr, pathogenic_domain_output, output_dir)

    # Process repeat and coding gene information (PM4 & BP3)
    if not skip_if_exists(repeat_output, "Repeat/Coding (PM4 & BP3)"):
        logging.info("Processing repeat and coding gene information... (PM4 & BP3)")
        process_repeat_and_coding_gene_info(rmsk_download_file, ccds_download_file, rmsk_file, ccds_file, repeat_output)

    # Process missense pathogenic genes and truncating genes (PP2 & BP1)
    if not skip_if_exists([missense_output, truncating_output], "Missense/Truncating (PP2 & BP1)"):
        logging.info("Processing missense pathogenic genes and truncating genes... (PP2 & BP1)")
        process_missense_pathogenic_genes_and_truncating_genes(clinvar_with_chr, gnomad_rmc_download_bb, gnomad_rmc_bb,
                                                               gnomad_rmc_bed, missense_output, truncating_output, output_dir, ref_b)

    # Download ClinGen Gene-Disease Validity (BS2, gene-level data)
    if not skip_if_exists(clingen_file, "ClinGen Gene-Disease Validity"):
        logging.info("Downloading ClinGen Gene-Disease Validity data... (BS2)")
        download_file = "Clingen-Gene-Disease-Summary.csv"
        if not os.path.exists(download_file):
            run_command(f"wget -O {download_file} {clingen_download_url}")
        shutil.copy(download_file, clingen_file)

    # Download gnomAD gene constraint metrics. Consumed by PVS1 (raw LOEUF drives
    # strength adjustment) and BS2 (LOEUF-tiered observed-count thresholds).
    if not skip_if_exists(gnomad_constraints_file, "gnomAD Gene Constraints"):
        logging.info("Downloading gnomAD gene constraint metrics... (PVS1 LOEUF, BS2 count tiers)")
        download_and_extract(gnomad_constraints_download_url, gnomad_constraints_file)

    # Fetch CSpec VCEP AF cutoffs (BA1, BS1, PM2 gene-specific overrides)
    if not skip_if_exists(vcep_af_output, "CSpec VCEP AF cutoffs (BA1/BS1/PM2)"):
        logging.info("Fetching VCEP AF cutoffs from ClinGen CSpec... (BA1/BS1/PM2)")
        import sys as _sys
        _orig_argv = _sys.argv
        _sys.argv = [
            "generate_vcep_af_cutoffs.py",
            "--output", vcep_af_output,
            "--review-output", vcep_af_review,
            "--allow-unparseable",
            "--verbose", "WARNING",
        ]
        try:
            generate_vcep_af_cutoffs.main()
        finally:
            _sys.argv = _orig_argv

    # ── Annotator-specific steps: PS1/PM5 AA extraction and PS4 submitter counts ──
    # These are the only steps that differ between annotators. When --annotator all,
    # we loop over both and generate separate output files + path JSONs.
    for ann in annotators_to_run:
        ann_tag = os.path.join(output_dir, f"{ref_b}_clean_clinvar_{ann}")
        aa_output = os.path.join(output_dir, f"{ref_b}_PS1_PM5_clinvar_pathogenic_aa_{ann}.tsv") if annotator == 'all' \
            else os.path.join(output_dir, f"{ref_b}_PS1_PM5_clinvar_pathogenic_aa.tsv")

        # Annotate ClinVar and extract AA changes (PS1/PM5)
        if not skip_if_exists(aa_output, f"ClinVar AA extraction - {ann} (PS1)"):
            logging.info("Extracting AA changes with %s annotator... (PS1/PM5)", ann)
            process_clinvar_aa_extraction(clean_clinvar, aa_output, ann, ann_tag, ref_b,
                                          nirvana_bin_dir=nirvana_bin_dir, nirvana_data_dir=nirvana_data_dir,
                                          output_dir=output_dir, vep_cache_dir=vep_cache_dir)

        # Process ClinVar submission-level pathogenic counts (PS4, VEP-only)
        if ann == 'vep':
            if not skip_if_exists(clinvar_submitter_counts_output, "ClinVar Submitter Counts (PS4)"):
                logging.info("Processing ClinVar submission-level pathogenic counts... (PS4)")
                if not os.path.exists(clinvar_submission_summary_gz):
                    run_command(f"wget -O {clinvar_submission_summary_gz} {clinvar_submission_summary_url}")
                if not os.path.exists(clinvar_variant_summary_gz):
                    run_command(f"wget -O {clinvar_variant_summary_gz} {clinvar_variant_summary_url}")
                assembly = "GRCh37" if ref_b == "hg19" else "GRCh38"
                generate_submitter_counts(clinvar_variant_summary_gz, clinvar_submission_summary_gz,
                                          clinvar_submitter_counts_output, assembly)

        # Write the output required_paths JSON for this annotator
        if annotator == 'all':
            output_json = os.path.join(output_dir, f"{ref_b}_{ann}_required_paths.json")
        else:
            output_json = os.path.join(output_dir, f"{ref_b}_required_paths.json")
        paths = {
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
            "BS2_clingen_gene_disease_validity_fp": clingen_file,
            "PVS1_PM2_BA1_BS1_BS2_gnomad_gene_constraints_fp": gnomad_constraints_file,
            "BA1_BS1_PM2_vcep_af_cutoffs_fp": vcep_af_output,
            "BA1_BS1_PM2_mis_oe_moi_af_cutoffs_fp": mis_oe_moi_af_cutoffs_file,
            "BA1_BS1_PM2_gene_to_moi_fp": gene_to_moi_file,
        }
        if ann == 'vep':
            paths["PS4_clinvar_submitter_counts_fp"] = clinvar_submitter_counts_output
        generate_output_json(paths, output_json)
        print(f"Setup complete for {ann}. Output written to {output_json}")

if __name__ == "__main__":
    main()
