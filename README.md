# BIAS-2015 #
Bitscopic Interpreting ACMG 2015 Standards
 
This software is dual licensed. You may choose the AGPL license or contact us for a commercial license to fit your needs.

## Setup ##

See doc/variant_interpretation.txt for generating the required data files. Once generated, update the paths in the file
config/requiredPaths.json to reflect the appropriate paths on your machine.

Alternatively, our version of the required data files can be downloaded here
https://github.com/bitscopic/BIAS-2015-data

Users will need the files in preprocessing_data/bias_2015_required_files/

BIAS-2015 exclusively uses Python3 standard libraries, you should not need to install any other dependencies.

To run the pipeline on a VCF file, you will first need to annotate it with Illumina Connected Annotations (ICA)

ICA (formerly Nirvana) is a dotnet software package that is free to Download online here
    https://support.illumina.com/downloads/illumina-connected-annotations.html

BIAS-2015 currently supports hg19/GRCh37 classifications, as such you must annotate with the ICA hg19/GRCh37 dataset.

```
mkdir data
dotnet IlluminaConnectedAnnotations/Downloader.dll --ga GRCh37 -o Data
```

We annotate our VCF with ICA/Nirvana in this manner
```
dotnet IlluminaConnectedAnnotations/Nirvana.dll --cache Data/Cache --sd Data/SupplementaryAnnotation/GRCh37 --ref Data/References/Homo_sapiens.GRCh37.Nirvana.dat --in test.vcf --o test_ica
```

If you are running DRAGEN >4.2, then you likely already have the output .json file in your results directory. You can
run BIAS-2015 directly on this json file. 

## Running the pipeline ##

Included in the repository is a set of test data. If the code is installed and setup correctly, the following commands should complete in <10 seconds.
```
python bias-2015.py test/data/bias-2015_test_file.json config/requiredPaths.json test_output.tsv
```

If you used the same data source files as the manuscript, you can verify the same output is generated with this command
```
diff test_output.tsv test/data/bias-2015-expected_test_output.tsv
```

To include your own classifiers, please provide the optional --user_classifiers flag
```
python bias-2015.py test/data/bias-2015_test_file.json config/requiredPaths.json test_output.tsv --user_classifiers my_classifiers.tsv
```

Classifier files are the same format as BIAS-2015 output! If you have many variants need to be updated, or you wish to update
your variant classifications with your own script, we recommend you run the pipeline first to make a classifier template then
update the template. Once you have your classifiers ready in the template file, rerun the pipeline. 

```
python bias-2015.py test/data/bias-2015_test_file.json config/requiredPaths.json test_output.tsv
mv test_output.tsv my_classifiers.tsv
```
 Either manually or programatically update my_classifiers.tsv to include your own ACMG classifiers. Then re run
```
python bias-2015.py test/data/bias-2015_test_file.json config/requiredPaths.json test_output.tsv --user_classifiers my_classifiers.tsv
```

### Unit tests ###
The unit tests are extremely helpful for developing or debugging specific portions of the code.
They are not required for running the pipeline nor do they provide any insight into the pipelines performance or
validity.

We use pytest as our unit test framework, users will need it installed and accessible to follow these testing commands.


To run the unit tests go to the test directory (some tests leverage local pathing)
```
cd test
```

Then use pytest 
```
pytest
```

If the unit tests are not working, there is likely something extremly wrong and we recommend re downloading the repository.

## Contribution guidelines ##

* All python code goes into the src directory.
* Every file in the src directory should have a paired file in the test directory.
* All python functions should have some degree of unit test coverage. This can be accomplished by either directly testing the
  function (preferred) or by testing a head function which calls other functions.
* Include publication citations directly in the code whenever possible.
* "To debug code you must be more clever than the person who wrote it, so don't be too clever when writing you own code".


## Who do I talk to? ##

Chris Eisenhart chris.eisenhart@bitscopic.com
Joel Mewton joel@bitscopic.com
