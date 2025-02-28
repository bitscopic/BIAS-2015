#!/bin/bash

# Set paths for necessary directories and files
ICA_BIN_DIR=~/Bioinformatics/bin/IlluminaConnectedAnnotations
ICA_DATA_DIR=~/Bioinformatics/data
CONFIG_DIR=~/Bioinformatics/config
SCRIPTS_DIR=../src/scripts
OUTPUT_DIR=$(pwd)  # This gets the absolute path to where the script is being run

# Make sure to use full paths in your configuration file.
echo "Running variant interpretation setup for BIAS-2015 on hg38/GRCh38"

################### PVS #####################

# PVS1 - Download the gnomAD constraint metrics
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv
python3 $SCRIPTS_DIR/find_lof_genes.py gnomad.v4.1.constraint_metrics.tsv hg38_lof_genes.tsv hg38

# PVS1 Caveat 2 - Download ncbiRefSeqHgmd from UCSC
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqHgmd.txt.gz
gunzip ncbiRefSeqHgmd.txt.gz
mv ncbiRefSeqHgmd.txt hg38_ncbiRefSeqHgmd.tsv

################### PS ######################

# PS1 & PM5 - Download the latest version of ClinVar for hg38
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
gunzip clinvar.vcf.gz

# Prepend 'chr' onto each line without a '#' in ClinVar VCF
cat clinvar.vcf | awk '{if($0 !~ /^#/) print "chr"$0; else print $0;}' > clinvar_with_chr.vcf

# Filter pathogenic and likely pathogenic variants
grep '^#' clinvar_with_chr.vcf > clean_clinvar.vcf
grep -E 'CLNSIG=Likely_pathogenic|CLNSIG=Pathogenic' clinvar_with_chr.vcf >> clean_clinvar.vcf

# Annotate the ClinVar VCF with Illumina Connected Annotations
dotnet $ICA_BIN_DIR/Nirvana.dll --cache $ICA_DATA_DIR/Cache --sd $ICA_DATA_DIR/SupplementaryAnnotation/GRCh38 --ref $ICA_DATA_DIR/References/Homo_sapiens.GRCh38.Nirvana.dat --in clean_clinvar.vcf --o clean_clinvar_ica

# Generate data files for interpretation
python3 $SCRIPTS_DIR/generate_pathogenic_aa_list.py clean_clinvar_ica.json.gz clean_clinvar.vcf hg38_clinvar_pathogenic_aa.tsv

# PS3 - Download avada big bed file and convert it to bed format
wget http://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/avada.bb
wget https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/bigBedToBed
chmod +x bigBedToBed
./bigBedToBed avada.bb avada.bed

# Extract necessary columns from avada.bed
python3 $SCRIPTS_DIR/extract_from_avada_track.py avada.bed hg38_literature_supported_variants.tsv

# PS4 - Download GWAS catalog
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gwasCatalog.txt.gz
gunzip gwasCatalog.txt.gz
mv gwasCatalog.txt hg38GwasCatalog.txt

################### PM ######################

# PM1 - Generate lists of pathogenic domains using ClinVar data
python3 $SCRIPTS_DIR/generate_domain_lists.py clinvar.vcf pathogenic_domains.tsv

# Convert UCSC UniProt domains track to tsv and sort
wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/uniprot/unipDomain.bb
./bigBedToBed unipDomain.bb unipDomain.bed
cat unipDomain.bed | cut -f 1-3,28 | awk '{print $4"\t"$1"\t"$2"\t"$3}' | sort > sorted_unipDomain.bed

# Join ClinVar domains with UniProt domains
join sorted_unipDomain.bed pathogenic_domains.tsv | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$5"\t"$6}' > hg38_pathogenic_domains.tsv

# PM4 & BP3 - Download repeat and coding gene information
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ccdsGene.txt.gz
gunzip rmsk.txt.gz
gunzip ccdsGene.txt.gz
python3 $SCRIPTS_DIR/join_coding_and_repeats.py --inRepeatFile rmsk.txt --inConsensusCodingFile ccdsGene.txt --out hg38_coding_repeat_regions.tsv

################### PP ######################

# PP2 & BP1 - Find missense pathogenic genes
python3 $SCRIPTS_DIR/find_missense_pathogenic_genes.py clinvar.vcf hg38_missense_pathogenic_genes.tsv hg38_truncating_genes.tsv

# Create output JSON file with absolute paths
cat << EOF > $OUTPUT_DIR/hg38_required_paths.json
{
    "lof_gene_fp": "$OUTPUT_DIR/hg38_lof_genes.tsv",
    "ncbi_ref_seq_hgmd_fp": "$OUTPUT_DIR/hg38_ncbiRefSeqHgmd.tsv",
    "clinvar_pathogenic_aa_fp": "$OUTPUT_DIR/hg38_clinvar_pathogenic_aa.tsv",
    "gwas_dbsnp_fp": "$OUTPUT_DIR/hg38GwasCatalog.txt",
    "coding_repeat_region_fp": "$OUTPUT_DIR/hg38_coding_repeat_regions.tsv",
    "missense_pathogenic_genes_fp": "$OUTPUT_DIR/hg38_missense_pathogenic_genes.tsv",
    "truncating_genes_fp": "$OUTPUT_DIR/hg38_truncating_genes.tsv",
    "pathogenic_domains_fp": "$OUTPUT_DIR/hg38_pathogenic_domains.tsv",
    "literature_supported_variants_fp": "$OUTPUT_DIR/hg38_literature_supported_variants.tsv"
}
EOF

echo "Setup complete. Output written to $OUTPUT_DIR/hg38_required_paths.json."
