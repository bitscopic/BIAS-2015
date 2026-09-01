# BIAS-2015 #
Bitscopic Interpreting ACMG 2015 Standards
 
This software is dual licensed. Commercial use requires a license. [Contact us](mailto:bill@bitscopic.com) to discuss a commercial license to fit your needs. The software is free for academic use and is made available under the AGPL license for that purpose.

BIAS-2015 also has a graphical user interface [BIAS-2015-ui](https://github.com/bitscopic/BIAS-2015-ui) available to
view and modify classification results.

BIAS-2015 is the core variant classification engine for Bitscopic's commercial platform that supports tertiary analysis and
report generation. BIAS can be extended to support lab specific customizations. For a commercial demo please [contact us](mailto:bill@bitscopic.com). 

## Setup ##

### Required Libraries ###

**NOTE:** All commands should be run from the `bias-2015` root directory unless otherwise specified.

#### General Use ####
BIAS-2015 requires Nirvana to annotate VCF files—see below for installation instructions.
Nirvana requires [.NET](https://www.microsoft.com/net/download/core) as a dependency.

BIAS-2015 exclusively uses Python3 standard libraries, so no additional Python libraries are needed.

The machine must be able to read and unzip `.gz` files. This requires different libraries depending on the  
operating system. See the official documentation for more details: [gzip.org](https://www.gzip.org/).

The machine must have **AWS CLI** installed. 

Install: `brew install awscli`

Verify: `aws --version`

#### Preprocessing ####
The preprocessing requires an executable (`bigBedToBed`) that is only available on Mac and Linux machines. Additionally
the preprocessing uses system calls to `wget` to download remote files, so this must be installed and accessible to Python

#### Validation ####
The validation requires several python libraries defined in doc/validation_requirements.txt. These are not required
for generic use. 

### Nirvana ###
To process your own VCF file, you must first install Nirvana 3.18.1.

```
mkdir bin
cd bin
wget https://github.com/Illumina/Nirvana/releases/download/v3.18.1/Nirvana-3.18.1-net6.0.zip
unzip Nirvana-3.18.1-net6.0.zip
cd ../
mkdir data
cd data
mkdir GRCh37
dotnet ~/bin/Nirvana-v3.18.1/Downloader.dll --ga GRCh37 --out ~/data/GRCh37/
mkdir GRCh38
dotnet ~/bin/Nirvana-v3.18.1/Downloader.dll --ga GRCh38 --out ~/data/GRCh38/
```

Nirvana's supplementary annotation files are no longer actively maintained by Illumina. Bitscopic
hosts updated versions of two key files — **AlphaMissense** (new) and **ClinVar** (February 2026) —
that should replace the stale versions bundled with Nirvana 3.18.1.

After downloading the Nirvana data above, replace the outdated ClinVar file and add AlphaMissense.
The same process applies to GRCh38 (substitute the build name in the paths below).

```
cd ~/data/GRCh37/SupplementaryAnnotation/GRCh37/

# Remove the outdated ClinVar files
rm ClinVar_20230930.nsa ClinVar_20230930.nsa.idx ClinVar_20230930.nsi

# Download the updated ClinVar and new AlphaMissense files
aws s3 cp s3://bias-2015/nirvana_data/GRCh37_updated/ . --no-sign-request --recursive
```

This adds three files per annotation source: the `.nsa`, `.nsa.idx`, and `.nsa.schema` for both
ClinVar (February 2026) and AlphaMissense. Nirvana will automatically pick up the new files on the
next run.

Then run Nirvana on your VCF file to generate a .json output file.
```
dotnet ~/bin/Nirvana-v3.18.1/Nirvana.dll \
  --cache ~/data/GRCh37/Cache/Both \
  --sd ~/data/GRCh37/SupplementaryAnnotation/GRCh37 \
  --ref ~/data/GRCh37/References/Homo_sapiens.GRCh37.Nirvana.dat \
  --in test/data/bias-2015_test_file.vcf \
  --o bias-2015_test_file
```

If you have already have a Nirvana (or an ICA) .json file, you can provide it to BIAS directly. Please ensure that the
BIAS preprocessing data was generated with the same version of the annotator you use to make your json files. 

If you would like to test the  BIAS-2015 software without installing Nirvana, we have provided a test .json file in our
test/data directory.

### VEP (Alternative Annotator) ###
BIAS-2015 also supports Ensembl's VEP as an alternative annotator. VEP requires Docker and several additional
annotation files (gnomAD, ClinVar, REVEL, AlphaMissense). See [doc/vep_setup.md](doc/vep_setup.md) for full
setup and usage instructions.

### Data files & Preprocessing ###

**Quick start**: If you just want to run BIAS, download the prebuilt files from S3 (recommended).
Running preprocessing locally is an advanced feature and is not required.

BIAS-2015 v3.0.0 requires multiple data files to run. These are provided to the algorithm through a required_paths.json file
that lists the expected file and its path.

#### Downloading Pre-built Data Files ####

The required BIAS data files for hg19 and hg38 are available from AWS S3. Files for both builds
are in the same directory, prefixed with `hg19_` or `hg38_`.

**v3.0.0 (recommended)** — individual files:
```
# List available files
aws s3 ls s3://bias-2015/v3.0.0_datasets/2026.03.01/ --no-sign-request

# Download all hg19 files
mkdir bias_hg19_data_files
aws s3 cp s3://bias-2015/v3.0.0_datasets/2026.03.01/ bias_hg19_data_files/ --no-sign-request --recursive --exclude "*" --include "hg19_*"

# Download all hg38 files
mkdir bias_hg38_data_files
aws s3 cp s3://bias-2015/v3.0.0_datasets/2026.03.01/ bias_hg38_data_files/ --no-sign-request --recursive --exclude "*" --include "hg38_*"
```

**NOTE:** The PS1/PM5 file differs by annotator (`*_nirvana.tsv` vs `*_vep.tsv`). VEP users also
need the `PS4_clinvar_submitter_counts.tsv` file. The download includes both annotator variants;
extra files are harmless and will be ignored.

**BA1 / BS1 / PM2 allele-frequency cutoffs.** These codes resolve their AF thresholds in a
three-tier cascade: (1) a gene-specific VCEP rule from the ClinGen CSpec registry
(`data/vcep/vcep_af_cutoffs.tsv`) if one exists, (2) a data-derived cutoff stratified by gnomAD
missense O/E upper bound and mode of inheritance (`data/af_cutoffs/mis_oe_upper_binned_cutoffs.tsv`
+ `data/af_cutoffs/gene_to_moi.tsv`), and (3) a flat ACMG default (`ACMG_DEFAULT_AF_CUTOFFS` in
`src/bias_2015/constants.py`). All three data files ship committed in the repo; no download is
required for the AF cutoffs. LOEUF is no longer used for BA1/BS1/PM2; it remains in use for PVS1
(strength adjustment) and BS2 (observed-count tiering).

The download also includes pre-built required_paths.json files (`hg19_nirvana_required_paths.json`,
`hg19_vep_required_paths.json`, etc.). You can use these directly by updating the file paths
inside to match your local directory. However, generating a new required_paths.json using the script below is recommended, as it automatically sets the correct paths for your local machine.
```
# For Nirvana (default)
python3 src/scripts/create_new_required_paths_file.py bias_hg19_data_files hg19 hg19_required_paths.json

# For VEP
python3 src/scripts/create_new_required_paths_file.py bias_hg19_data_files hg19 hg19_vep_required_paths.json --annotator vep
```

**v2.#.# (legacy)** — zip archives for Nirvana only:
**NOTE:** This step can be skipped when using the recommended v3.0.0 data files.
```
aws s3 cp s3://bias-2015/v2.0.0_datasets/bias_v2.0.0_hg19_data_files.zip . --no-sign-request
aws s3 cp s3://bias-2015/v2.0.0_datasets/bias_v2.0.0_hg38_data_files.zip . --no-sign-request
```

#### Running Preprocessing Locally ####

Alternatively, users can generate their own BIAS data files by running the preprocessing locally. The entire preprocessing
pipeline has been included and can be run for hg19 or hg38 with a single command. This process takes multiple hours, will
download multiple GB of files to the running machine, and will use many GB of disk space as intermediate files. A populated
required_paths.json file will be written at the end of preprocessing and can be used immediately with BIAS. Please note that
hg38 requires significantly more disk space and downloads than hg19.

For **Nirvana** preprocessing (requires Nirvana installed):
```
mkdir bias_hg19_data_files
cd bias_hg19_data_files
python3 ../preprocessing.py \
   --reference_build hg19 \
   --annotator nirvana \
   --output_dir . \
   --os_type linux \
   --nirvana_bin_dir ~/bin/Nirvana-v3.18.1 \
   --nirvana_data_dir ~/data/GRCh37/ \
   --verbose=DEBUG
```

For **VEP** preprocessing (requires VEP Docker image and cache, see [doc/vep_setup.md](doc/vep_setup.md)):
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

## Running the Pipeline ##

A Nirvana .json output can be passed through BIAS-2015 with the following command structure. 

**NOTE:** Users must download or generate the BIAS-data and have a valid required_paths.json file before running BIAS! 

We provide a test data set comprised of 100 randomly selected eRepo variants at test/data/bias-2015_test_file.vcf. 
Also included is the Nirvana annotation file generated by processing the above VCF with Nirvana 3.18.1 and the static
Nirvana GRCh37 dataset.

You can run BIAS-2015 on the test data as follows. 
```
python3 bias_2015.py test/data/bias-2015_test_file.nirvana.json hg19_required_paths.json test_output.tsv
```

If you downloaded the v3.0.0 hg19 data files and correctly generated a required paths json, then
this diff will show no differences.
```
diff test_output.tsv test/data/bias-2015_test_file.nirvana.bias_output.tsv
```

You can view the output file manually or through any tsv reader (excel) to view the classification and rationale
assigned to each variant. Alternatively users can use the [BIAS-2015-ui](https://github.com/bitscopic/BIAS-2015-ui)
which simplifies viewing and manually updating results.

## eRepo VCF File ##
Our latest conversion of the eRepo to vcf is available here;

test/data/erepo_03.02.2026.vcf

We believe it is the best resource for benchmarking ACMG classifier algorithms as you can evaluate code F1 when
compared to the expert panels.  For more information on this process please refer to our Genome Medicine publication.

### Providing User Classifications ###
User classifications are the recommended way to populate the 9 ACMG codes that BIAS cannot automatically assign or to
override a classification with a users expert opinion. 

We have created the [BIAS-2015-ui](https://github.com/bitscopic/BIAS-2015-ui) which enables users to upload BIAS-2015
output files, view and modify them through the GUI, then save them as an user classifier file. If you would like to see
an example of a modified file, please use the UI to modify our test output, then download it and diff to see how your 
changes were stored.  

To include your own classifiers, please provide them in their own file using the optional --user_classifiers flag
```
python3 bias_2015.py test/data/bias-2015_test_file.json hg19_required_paths.json test_output.tsv --user_classifiers my_classifiers.tsv
```
Classifier files are the same format as BIAS-2015 output with the users updates included inline for each variant. When 
a user classifier file is provided, the user classifiers will override any algorithmic classifiers.  The algorithm will
attempt to assign any unassigned codes, then recalculate the combination criteria to determine the variants new 
classification. 

If you have many variants that need to be updated, or you wish to update your variant classifications with your own script, 
we recommend you run the pipeline first to make a classifier template then update the template. Once you have your 
classifiers ready in the template file, rerun the pipeline. 

You may update the template file via the GUI, by hand or programmatically.

Example
```
python3 bias_2015.py test/data/bias-2015_test_file.json hg19_required_paths.json test_output.tsv
mv test_output.tsv my_classifiers.tsv
```
Either manually or programmatically update my_classifiers.tsv to include your own ACMG classifiers. Then re run
```
python3 bias_2015.py test/data/bias-2015_test_file.json hg19_required_paths.json test_output.tsv --user_classifiers my_classifiers.tsv
```

## What do I cite? ##
Please use the following citation information when referencing BIAS 2015

Eisenhart, C., Brickey, R., Nadon, B. et al. Automating ACMG variant classifications with BIAS-2015 v2.1.1: algorithm 
analysis and benchmark against the FDA-approved eRepo dataset. Genome Med 17, 148 (2025).
https://doi.org/10.1186/s13073-025-01581-y


## Who do I talk to? ##
For algorithm questions please contact - 

Chris Eisenhart chris.eisenhart@bitscopic.com 

Rachel Brickey rachel@bitscopic.com

Joel Mewton joel@bitscopic.com
