# SMTa
Spatial meta-transcriptome analysis pipeline

## requirement
SMTa depends on several toolkits, please install/download following tool/database/R package in your linux and R respectively before using SMTa.

1. spaceranger1.0/2.0 (see https://support.10xgenomics.com/spatial-gene-expression/software/downloads/latest)
2. blastn 2.10.1+ and NT blast database (see https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
   if you have conda on server, you can install blast with desired version directly as following:
```Bash
conda install blast=2.10.1 -c bioconda
```
3. Seurat 4.0 in R (this require R to be updated to 4.0.0 or higher)
4. python 3.7.6 or higher

## install
```Bash
git clone https://github.com/affinis/SMTa.git
cd SMTa
bash install_SMTa.sh -n $ABSOLUTE_PATH_OF_NT_DATABASE -b $ABSOLUTE_PATH_OF_BLASTN -s $ABSOLUTE_PATH_OF_SPACERANGER
```
## test run
run test data under the directory of spaceranger software 
```Bash
./run_SMTa_test.sh
```

## run SMTa with your own data start from fastq
Two ways are provided for users to run their own data.
1. specify arguments with command line
```Bash
./run_SMTa.sh -d $ABSOLUTE_PATH_OF_SAMPLE_DIRECTORY -s $SAMPLE_ID -p $ABSOLUTE_PATH_OF_SAMPLE_IMAGE -r $ABSOLUTE_REFERENCE_PATH
```
2. specify arguments with the configuration file "sample.info", this way can run 4 sample (sample number with a typical slide) in order.\
Before that, you should modify "sample.info" in SMTa directory, remember all the path should be absolute path.
```Bash
./run_SMTa.sh -c sample.info
```

## run SMTa start from an existing directory of spaceranger output

## analysis with R

## trouble shooting
