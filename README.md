# BIAS-2015 #
Bitscopic Interpreting ACMG 2015 Standards
 
This software is dual licensed. You may choose the AGPL license or [contact us](mailto:bill@bitscopic.com) for a commercial license to fit your needs.
It is free for academic use.

BIAS-2015 also has a graphical user interface [BIAS-2015-ui](https://github.com/bitscopic/BIAS-2015-ui) available to
view and modify classification results.


## Setup ##

### Required Libraries ###

#### General Use ####
BIAS-2015 requires Nirvana to annotate VCF files—see below for installation instructions.
Nirvana requires [.NET](https://www.microsoft.com/net/download/core) as a dependency.

BIAS-2015 exclusively uses Python3 standard libraries, so no additional Python libraries are needed.

The machine must be able to read and unzip `.gz` files. This requires different libraries depending on the  
operating system. See the official documentation for more details: [gzip.org](https://www.gzip.org/).

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

Then run Nirvana on your VCF file to generate a .json output file.

```
dotnet ~/bin/Nirvana-v3.18.1/Nirvana.dll \
  --cache ~/data/GRCh37/Cache/Both \
  --sd ~/data/GRCh37/SupplementaryAnnotation/GRCh37 \
  --ref ~/data/GRCh37/References/Homo_sapiens.GRCh37.Nirvana.dat \
  --in test/data/bias-2015_test_file.vcf \
  --o bias-2015_test_file
```

If you have already have a Nirvana (or an ICA) .json file, you can provide it to BIAS directly. If you would like to test the 
BIAS-2015 software without installing Nirvana, we have provided a test .json file in our test/data directory.

### Data files & Preprocessing ###

BIAS-2015 v2.0.0 requires multiple data files to run. These are provided to the algorithm through a required_paths.json file
that lists the expected file and its path.

The required BIAS data files for hg19 or hg38 can be viewed and downloaded from AWS here. 
```
aws s3 ls s3://bias-2015 --no-sign-request
aws s3 cp s3://bias-2015/bias_v2.0.0_hg19_data_files.zip . --no-sign-request
aws s3 cp s3://bias-2015/bias_v2.0.0_hg38_data_files.zip . --no-sign-request
```

Only data sets v2.0.0+ are compatible with BIAS-2015 v2.0.0. Once the user has all the BIAS-2015 data files in a single directory,
they can generate a required_paths.json file by running the provided helper script.

```
unzip bias_v2.0.0_hg19_data_files.zip
python3 generate_required_paths_json.py bias_v2.0.0_hg19_data_files hg19 hg19_required_paths.json
```

Alternatively, users can generate their own BIAS data files by running the preprocessing locally.  The entire preprocessing
pipeline has been included and can be ruGn for hg19 or hg38 with a single command.  This process takes multiple hours, will
download multiple GB of files to the running machine, and will use many GB of disk space as intermediate files. Nirvana is 
required to run the preprocessing. A populated required_paths.json file will be written at the end of preprocessing and can
be used immediately with BIAS. Please note that hg38 requires significantly more disk space and downloads than hg19.

```
mkdir bias_hg19_data_files
cd bias_hg19_data_files 
python3 ../preprocessing.py \
   --reference_build hg19 \
   --out . \
   --os_type linux \
   --nirvana_bin_dir ~/bin/Nirvana-v3.18.1 \
   --nirvana_data_dir ~/data/GRCh37/ \
   --verbose=DEBUG
```

## Running the Pipeline ##

A Nirvana .json output can be passed through BIAS-2015 with the following command structure. 

NOTE - users must download or generate the BIAS-data and have a valid required_paths.json file before running BIAS! 

```
python bias-2015.py test/data/bias-2015_test_file.json hg19_required_paths.json test_output.tsv
```

If you downloaded the bias_v2.0.0_hg19_data_files.zip dataset, and correctly generated a required paths json, then
this diff should show no differences.
```
diff test_output.tsv test/data/bias-2015-expected_test_output.tsv
```

You can view the output file manually or through any tsv reader (excel) to view the classification and rationale
assigned to each variant. Alternatively users can use the [BIAS-2015-ui](https://github.com/bitscopic/BIAS-2015-ui)
which simplifies viewing and manually updating results.


### Providing User Classifications ###
To include your own classifiers, please provide the optional --user_classifiers flag
```
python bias-2015.py test/data/bias-2015_test_file.json hg19_required_paths.json test_output.tsv --user_classifiers my_classifiers.tsv
```

Classifier files are the same format as BIAS-2015 output. If you have many variants need to be updated, or you wish to update
your variant classifications with your own script, we recommend you run the pipeline first to make a classifier template then
update the template. Once you have your classifiers ready in the template file, rerun the pipeline. 

You may update the template file by hand or programmatically. Alternatively we have created the [BIAS-2015-ui](https://github.com/bitscopic/BIAS-2015-ui)
which enables users to upload BIAS-2015 output files, view and modify them through the GUI, then save them as an user classifier
file.

Example
```
python bias-2015.py test/data/bias-2015_test_file.json hg19_required_paths.json test_output.tsv
mv test_output.tsv my_classifiers.tsv
```
Either manually or programmatically update my_classifiers.tsv to include your own ACMG classifiers. Then re run
```
python bias-2015.py test/data/bias-2015_test_file.json hg19_required_paths.json test_output.tsv --user_classifiers my_classifiers.tsv
```

## Who do I talk to? ##

Chris Eisenhart chris.eisenhart@bitscopic.com
Rachel Brickey rachel@bitscopic.com
Brian Nadon brian@bitscopic.com
Joel Mewton joel@bitscopic.com
