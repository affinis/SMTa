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
Run test data under the directory of spaceranger software 
```Bash
./run_SMTa_test.sh
```

## run SMTa with your own data start from fastq
Two ways are provided for users to run their own data.
1. Specify arguments with command line
```Bash
./run_SMTa.sh -d $ABSOLUTE_PATH_OF_SAMPLE_DIRECTORY -s $SAMPLE_ID -p $ABSOLUTE_PATH_OF_SAMPLE_IMAGE -r $ABSOLUTE_REFERENCE_PATH
```
2. Specify arguments with the configuration file "sample.info", this way can run 4 sample (sample number with a typical slide) in order.\
Before that, you should modify "sample.info" in SMTa directory, remember all the path should be absolute path, here is an example:
```Bash
#!/bin/bash

#specify reference here
REF=/tmpdata/LyuLin/Tools/spaceranger-1.3.1/external/spaceranger_tiny_ref/1.0.0

#specify fastq path
SAMPLE_PATH=/tmpdata/LyuLin/Tools/spaceranger-1.3.1/external/spaceranger_tiny_inputs/fastqs

#first field of fastq file name, for example, "tinytest_S1_L001_R1.fastq.gz" and "tinytest_S1_L001_R2.fastq.gz" are 
#fastqs of one sample, the sample ID is "tinytest".
A="tinytest"
B=""
C=""
D=""

#specify image path, take care that for spaceranger2.0, it can reorient image automaticly while spaceranger1.0 not, 
#for spaceranger1.0 you should reorient the image to put the sandglass-like shape at upleft position.
A_image=/tmpdata/LyuLin/Tools/spaceranger-1.3.1/external/spaceranger_tiny_inputs/image/tinyimage.jpg
B_image=""
C_image=""
D_image=""
```
Then, just run SMTa using following command.
```Bash
./run_SMTa.sh -c sample.info
```

## run SMTa start from spaceranger output
If you have run spaceranegr before with your own visium data (or you have downloaded some public datasets from 10x genomics official website), SMTa can also be applied when bam files are available in these data.
```Bash
cd $SPACERANGER_OUTS_DIR
bash $PATH_TO_SMTa/src/microperc.sh -u 1 -b possorted_genome_bam.bam -d $ABSOLUTE_PATH_OF_NT_DATABASE
python3 $PATH_TO_SMTa/src/SMT_mtx_gener_blast.py > SMT_summary.txt
```

## analysis with R

## trouble shooting
