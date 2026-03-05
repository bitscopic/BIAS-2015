# VEP (Variant Effect Predictor) Setup for BIAS-2015 #

BIAS-2015 supports Ensembl's VEP as an alternative to Nirvana for variant annotation. This guide covers
everything needed to set up VEP and run BIAS-2015 with VEP-annotated variants. Both hg19/GRCh37 and
hg38/GRCh38 reference builds are supported.

**IMPORTANT: BIAS-2015 preprocessing data files are annotator-specific. You must use VEP preprocessing files
with VEP-annotated variants. Do not use Nirvana preprocessing files with VEP or vice versa — this will produce
incorrect classifications.**

## Prerequisites ##

Docker must be installed and running.
- macOS: https://docs.docker.com/desktop/install/mac-install/
- Linux: https://docs.docker.com/engine/install/

Disk space requirements:
- VEP cache: ~15 GB per assembly
- gnomAD v2.1.1 exomes VCF: ~60 GB
- ClinVar VCF: ~300 MB
- REVEL scores: ~3 GB
- AlphaMissense scores: ~700 MB
- **Total: ~80 GB minimum per assembly**

## VEP Annotation Files ##

VEP requires several annotation files that are used as custom annotations and plugins during variant annotation.
The files differ between hg19 and hg38 — use the section that matches your reference build.

### hg19/GRCh37 ###

#### ClinVar (from NCBI) ####

```
mkdir -p vep_data/hg19
cd vep_data/hg19

# Download ClinVar VCF + index (GRCh37) directly from NCBI
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi
```

#### REVEL and AlphaMissense (from Bitscopic S3) ####

These files are processed for VEP compatibility (source: [REVEL](https://sites.google.com/site/revelgenomics),
[AlphaMissense](https://zenodo.org/records/10813168)).

```
cd vep_data/hg19

# Download REVEL scores + index
aws s3 cp s3://bias-2015/vep/revel_final.tsv.gz . --no-sign-request
aws s3 cp s3://bias-2015/vep/revel_final.tsv.gz.tbi . --no-sign-request

# Download AlphaMissense scores + index (hg19)
aws s3 cp s3://bias-2015/vep/AlphaMissense_hg19.tsv.gz . --no-sign-request
aws s3 cp s3://bias-2015/vep/AlphaMissense_hg19.tsv.gz.tbi . --no-sign-request
```

#### From the Broad Institute (gnomAD v2.1.1) ####

gnomAD v2.1.1 is on GRCh37 and is hosted publicly by the Broad Institute.

```
mkdir -p vep_data/hg19/gnomad
cd vep_data/hg19/gnomad

# Download gnomAD v2.1.1 exomes VCF (~60 GB - this will take a while)
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi
```

Alternatively, gnomAD is also available on AWS:
```
aws s3 cp s3://gnomad-public-us-east-1/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz . --no-sign-request
aws s3 cp s3://gnomad-public-us-east-1/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi . --no-sign-request
```

### hg38/GRCh38 ###

#### ClinVar (from NCBI) ####

```
mkdir -p vep_data/hg38
cd vep_data/hg38

# Download ClinVar VCF + index (GRCh38) directly from NCBI
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
```

#### REVEL and AlphaMissense (from Bitscopic S3) ####

These files are processed for VEP compatibility (source: [REVEL](https://sites.google.com/site/revelgenomics),
[AlphaMissense](https://zenodo.org/records/10813168)).

```
cd vep_data/hg38

# Download REVEL scores + index (GRCh38)
aws s3 cp s3://bias-2015/vep/hg38/revel_grch38.tsv.gz . --no-sign-request
aws s3 cp s3://bias-2015/vep/hg38/revel_grch38.tsv.gz.tbi . --no-sign-request

# Download AlphaMissense scores + index (hg38)
aws s3 cp s3://bias-2015/vep/hg38/AlphaMissense_hg38.tsv.gz . --no-sign-request
aws s3 cp s3://bias-2015/vep/hg38/AlphaMissense_hg38.tsv.gz.tbi . --no-sign-request
```

#### From the Broad Institute (gnomAD v4.1) ####

gnomAD v4.1 is on GRCh38 and is hosted publicly by the Broad Institute. The exomes VCF is split by
chromosome, so you need to download each one individually.

```
mkdir -p vep_data/hg38/gnomad
cd vep_data/hg38/gnomad

# Download gnomAD v4.1 exomes VCFs by chromosome
for chr in $(seq 1 22) X Y; do
    wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr${chr}.vcf.bgz
    wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr${chr}.vcf.bgz.tbi
done
```

After downloading, concatenate into a single VCF for VEP:
```
# Concatenate per-chromosome VCFs into a single file
bcftools concat gnomad.exomes.v4.1.sites.chr{1..22}.vcf.bgz \
    gnomad.exomes.v4.1.sites.chrX.vcf.bgz \
    gnomad.exomes.v4.1.sites.chrY.vcf.bgz \
    -Oz -o gnomad.exomes.v4.1.sites.vcf.bgz
tabix -p vcf gnomad.exomes.v4.1.sites.vcf.bgz
```

## VEP Docker Image and Cache ##

Pull a pinned version of the VEP Docker image and install the cache for your reference build.

```
docker pull ensemblorg/ensembl-vep:release_115.0

mkdir -p ~/.vep
```

For **hg19/GRCh37**:
```
docker run --rm -v ~/.vep:/opt/vep/.vep ensemblorg/ensembl-vep:release_115.0 \
    perl INSTALL.pl --AUTO cf --ASSEMBLY GRCh37 --SPECIES homo_sapiens --NO_UPDATE
```

For **hg38/GRCh38**:
```
docker run --rm -v ~/.vep:/opt/vep/.vep ensemblorg/ensembl-vep:release_115.0 \
    perl INSTALL.pl --AUTO cf --ASSEMBLY GRCh38 --SPECIES homo_sapiens --NO_UPDATE
```

After setup, your VEP annotation files should be organized as follows:

**hg19:**
```
vep_data/hg19/
├── clinvar.vcf.gz
├── clinvar.vcf.gz.tbi
├── revel_final.tsv.gz
├── revel_final.tsv.gz.tbi
├── AlphaMissense_hg19.tsv.gz
├── AlphaMissense_hg19.tsv.gz.tbi
└── gnomad/
    ├── gnomad.exomes.r2.1.1.sites.vcf.bgz
    └── gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi

~/.vep/
└── homo_sapiens/
    └── 115_GRCh37/
        └── (cache files and FASTA)
```

**hg38:**
```
vep_data/hg38/
├── clinvar.vcf.gz
├── clinvar.vcf.gz.tbi
├── revel_grch38.tsv.gz
├── revel_grch38.tsv.gz.tbi
├── AlphaMissense_hg38.tsv.gz
├── AlphaMissense_hg38.tsv.gz.tbi
└── gnomad/
    ├── gnomad.exomes.v4.1.sites.vcf.bgz
    └── gnomad.exomes.v4.1.sites.vcf.bgz.tbi

~/.vep/
└── homo_sapiens/
    └── 115_GRCh38/
        └── (cache files and FASTA)
```

## Annotating Your VCF with VEP ##

Run VEP on your input VCF to produce a JSON annotation file that BIAS-2015 can consume. Replace
`/path/to/vep_data` with the absolute path to your `vep_data` directory, and update the input/output
file paths as needed.

### hg19/GRCh37 ###

```
docker run --rm \
    -v "$PWD:$PWD" \
    -v "/path/to/vep_data/hg19:/data" \
    -v "$HOME/.vep:/opt/vep/.vep" \
    ensemblorg/ensembl-vep:release_115.0 vep \
    --input_file $PWD/your_variants.vcf \
    --output_file $PWD/your_variants.vep.json \
    --json \
    --cache \
    --offline \
    --assembly GRCh37 \
    --dir_cache /opt/vep/.vep \
    --fasta /opt/vep/.vep/homo_sapiens/115_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
    --everything \
    --custom file=/data/gnomad/gnomad.exomes.r2.1.1.sites.vcf.bgz,short_name=gnomAD_exomes,format=vcf,type=exact,coords=0,fields=AF%AF_popmax%AC%AN \
    --custom file=/data/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNVC%CLNSIGSCV \
    --plugin AlphaMissense,file=/data/AlphaMissense_hg19.tsv.gz \
    --plugin REVEL,file=/data/revel_final.tsv.gz \
    --force_overwrite
```

### hg38/GRCh38 ###

```
docker run --rm \
    -v "$PWD:$PWD" \
    -v "/path/to/vep_data/hg38:/data" \
    -v "$HOME/.vep:/opt/vep/.vep" \
    ensemblorg/ensembl-vep:release_115.0 vep \
    --input_file $PWD/your_variants.vcf \
    --output_file $PWD/your_variants.vep.json \
    --json \
    --cache \
    --offline \
    --assembly GRCh38 \
    --dir_cache /opt/vep/.vep \
    --fasta /opt/vep/.vep/homo_sapiens/115_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    --everything \
    --custom file=/data/gnomad/gnomad.exomes.v4.1.sites.vcf.bgz,short_name=gnomAD_exomes,format=vcf,type=exact,coords=0,fields=AF%AF_joint%AC%AN \
    --custom file=/data/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNVC%CLNSIGSCV \
    --plugin AlphaMissense,file=/data/AlphaMissense_hg38.tsv.gz \
    --plugin REVEL,file=/data/revel_grch38.tsv.gz \
    --force_overwrite
```

The `--json` flag produces one JSON object per line, which is the format BIAS-2015 expects for VEP input.
The `--everything` flag enables all standard VEP annotations (SIFT, PolyPhen, regulatory, etc.).
The `--offline` flag ensures VEP uses only local cache data and does not make network requests.

**Note:** gnomAD v4.1 (hg38) uses `AF_joint` instead of `AF_popmax` for the maximum population allele
frequency. BIAS-2015 handles this field difference automatically.

## BIAS VEP Data Files & Preprocessing ##

BIAS-2015 requires preprocessed data files to run. These files are annotator-specific — you must use the
VEP version when running with `--annotator vep`.

### Downloading Pre-built Data Files (v3.0.0) ###

The required VEP BIAS data files for both hg19 and hg38 are available from AWS S3. Files for both
builds are in the same directory, prefixed with `hg19_` or `hg38_`.

```
# List available files
aws s3 ls s3://bias-2015/v3.0.0_datasets/2026.03.01/ --no-sign-request

# Download all hg19 VEP files
mkdir bias_hg19_vep_data_files
aws s3 cp s3://bias-2015/v3.0.0_datasets/2026.03.01/ bias_hg19_vep_data_files/ --no-sign-request --recursive --exclude "*" --include "hg19_*"

# Download all hg38 VEP files
mkdir bias_hg38_vep_data_files
aws s3 cp s3://bias-2015/v3.0.0_datasets/2026.03.01/ bias_hg38_vep_data_files/ --no-sign-request --recursive --exclude "*" --include "hg38_*"
```

The download includes 12 shared data files plus the VEP-specific files:
- `{ref}_PS1_PM5_clinvar_pathogenic_aa_vep.tsv` (PS1/PM5 — use the `_vep.tsv` variant, not `_nirvana.tsv`)
- `{ref}_PS4_clinvar_submitter_counts.tsv` (PS4 — VEP only)

The download also includes pre-built required_paths.json files (`hg19_vep_required_paths.json`,
`hg38_vep_required_paths.json`). You can use these directly — just update the file paths inside
to match your local directory. Alternatively, generate a new one:
```
# hg19
python3 src/scripts/create_new_required_paths_file.py bias_hg19_vep_data_files hg19 hg19_vep_required_paths.json --annotator vep

# hg38
python3 src/scripts/create_new_required_paths_file.py bias_hg38_vep_data_files hg38 hg38_vep_required_paths.json --annotator vep
```

### Running Preprocessing Locally ###

Users can generate the VEP data files by running the preprocessing locally for either build. This requires
the VEP Docker image and cache described above. This process takes multiple hours and will download
multiple GB of files. A populated required_paths.json file will be written at the end of preprocessing
and can be used immediately with BIAS.

**hg19:**
```
mkdir bias_hg19_vep_data_files
cd bias_hg19_vep_data_files
python3 ../preprocessing.py \
   --reference_build hg19 \
   --annotator vep \
   --output_dir . \
   --os_type linux \
   --vep_cache_dir ~/.vep \
   --verbose=DEBUG
```

**hg38:**
```
mkdir bias_hg38_vep_data_files
cd bias_hg38_vep_data_files
python3 ../preprocessing.py \
   --reference_build hg38 \
   --annotator vep \
   --output_dir . \
   --os_type linux \
   --vep_cache_dir ~/.vep \
   --verbose=DEBUG
```

**hg38 note:** The gnomAD regional missense constraint (RMC) track is not available for hg38 from UCSC.
This means PP2 classification will rely solely on ClinVar-derived gene-level missense statistics rather
than positional constraint data. All other ACMG criteria are fully supported.

## Running the Pipeline ##

Once you have the VEP annotation file and VEP-specific BIAS data files, run BIAS-2015 with the
`--annotator vep` flag.

### hg19 ###

We provide a pre-annotated VEP test file so you can verify your setup without running VEP yourself.
The test data set is the same 100 randomly selected eRepo variants used for the Nirvana test
(`test/data/bias-2015_test_file.vcf`), annotated with VEP.

```
python3 bias_2015.py test/data/bias-2015_test_file.vep.json hg19_vep_required_paths.json test_output.tsv --annotator vep
```

If you downloaded the BIAS VEP data files and correctly generated a required paths json, then this diff
will show no differences.
```
diff test_output.tsv test/data/bias-2015_test_file.vep.bias_output.tsv
```

### hg38 ###

```
python3 bias_2015.py your_variants.vep.json hg38_vep_required_paths.json output.tsv --annotator vep
```

You can view the output file manually or through any tsv reader (excel) to view the classification and
rationale assigned to each variant. Alternatively users can use the
[BIAS-2015-ui](https://github.com/bitscopic/BIAS-2015-ui) which simplifies viewing and manually
updating results.

## Troubleshooting ##

### Docker permission issues ###
If you get permission errors, ensure your user is in the docker group:
```
sudo usermod -aG docker $USER
# Log out and back in for changes to take effect
```

### VEP cache not found ###
Ensure the cache is installed at `~/.vep` and that the Docker command mounts it correctly
(`-v "$HOME/.vep:/opt/vep/.vep"`). The cache version (115) must match the Docker image version.

### gnomAD VCF indexing issues ###
If VEP reports index issues with gnomAD, regenerate the index:

**hg19:**
```
docker run --rm -v /path/to/vep_data/hg19:/data ensemblorg/ensembl-vep:release_115.0 \
    tabix -p vcf /data/gnomad/gnomad.exomes.r2.1.1.sites.vcf.bgz
```

**hg38:**
```
docker run --rm -v /path/to/vep_data/hg38:/data ensemblorg/ensembl-vep:release_115.0 \
    tabix -p vcf /data/gnomad/gnomad.exomes.v4.1.sites.vcf.bgz
```

### Memory issues ###
VEP with gnomAD custom annotations can be memory-intensive. Ensure Docker has at least 8 GB of
memory allocated (Docker Desktop > Settings > Resources).
