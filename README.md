# SMTa
Spatial meta-transcriptome analysis pipeline

## requirement
SMTa depends on several toolkits, please install/download following tool/database/R package in your linux and R respectively before using SMTa, paths should be configured in 'bin/configs.conf' before using the pipeline.

1. spaceranger1.0/2.0
2. blastn 2.10.1+ and NT blast database
3. Seurat4.0 in R (this require R to be updated to 4.0.0 or higher)
4. python 3.7.6

## install
git clone https://github.com/affinis/SMTa.git
cd SMTa
bash install_SMTa.sh -n $PATH_OF_YOUR_NT_DATABASE

## test run

## run SMTa with your own data start from fastq

## run SMTa start from an existing directory of spaceranger output

## analysis with R
