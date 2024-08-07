# BIAS-2015 #
Bitscopic Interpreting ACMG 2015 Standards
 
This software is dual licensed. You may choose the AGPL license or contact us for a commercial license to fit your needs.

## Setup ##

See doc/variant_interpretation.txt for generating the required data files. Once generated, update the paths in the file
config/requiredPaths.json to reflect the appropriate paths on your machine.

BIAS-2015 exclusively uses Python3 standard libraries, you should not need to install any other dependencies.

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
