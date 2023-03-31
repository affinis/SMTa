# SMT
Spatial meta-transcriptome analysis pipeline \
publication: https://genome.cshlp.org/content/early/2023/03/27/gr.277178.122.full.pdf 

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
Final results should include 'matrix.mdx' and 'possorted_genome_bam.merged.microbe.tsv' and several tmp file except for spaceranger outputs.

## import microbial data into seurat object
It should be noticed that for spaceranger 2.x.x, there is a file called 'tissue_positions.csv' in dir 'spatial' instead of 'tissue_positions_list.csv' in the context of spaceranger 1.x.x, this will cause error while reading this file with 'CreateSeuratFromSMT', just use 'sed "1d" tissue_positions.csv > tissue_positions_list.csv' to solve it.<br>

Here's an example showing a basic usage of functions in 'SMT_core_functions.R' using a test data from our paper.<br>

First, read data under spaceranger outs dir including microbial data
```
setwd('$PATH_TO_SMTa')
source('./SMT_core_functions.R')
msrt<-CreateSeuratFromSMT('./output/testdata/')
```
Then, get an overview of microbe distribution if you like.
```
SMTPlotMod(msrt,'nCount_Microbe',log.scale=T)
```
![图片](https://user-images.githubusercontent.com/92193789/198481653-cf5332c1-895c-4f48-834e-75312ebcf5ea.png)

You can get composition of microbe-containing spots by use 'SMTStackPlot':
```
SMTStackPlot(msrt,min.nCount = 50,taxa.threshold = 0.1)
```
![图片](https://user-images.githubusercontent.com/92193789/198483824-1933d177-c796-4021-ab53-3b1a8b0b4895.png)

Or you want higher taxa level view:
```
msrt<-SMTlevel(msrt,"kingdom")
SMTStackPlot(msrt,min.nCount = 50,taxa.threshold = 0.1)
```
![图片](https://user-images.githubusercontent.com/92193789/198485553-304acc29-cd5e-405e-a8ac-ee3bfb13d436.png)

You can view certain taxon in feature plot.
```
msrt<-SMTlevel(msrt,"genus")
SMTPlotMod(msrt,'genus-Lactobacillus',log.scale=T)
```
![图片](https://user-images.githubusercontent.com/92193789/198486673-4bbb8152-9e17-4472-87c5-967d84cfd90c.png)

You can also view multiple taxa in a single plot (not recomended, just for showing):
```
SMTCoPlot(test,c("genus-Clostridium","genus-Duncaniella","genus-Staphylococcus"),tissue.non.tissue.colors = c("grey","white"),alpha = c(0.25,1),min.count = 8,log = T)
```
![图片](https://user-images.githubusercontent.com/92193789/198494247-bacfcc63-cddd-4679-893b-cfafa2bcb10a.png)



## trouble shooting
