#!/bin/bash

set -e
set -u
#set -o pipefail

SAMPLE=$1
INPUTFILE=$2
OUTPUTDIR=$3

source ./configs.conf

cd $OUTPUTDIR
zcat $INPUTFILE | head -n 400 | grep -A 1 '^@' | sed 's/^@/>/g' | sed '/^--$/d' > $OUTPUTDIR/$SAMPLE.sample.fasta
$BLASTN -query $OUTPUTDIR/$SAMPLE.sample.fasta -db $REF3 -max_target_seqs 5 -outfmt 6 -perc_identity 99 -qcov_hsp_perc 60 -out $OUTPUTDIR/Prealign.tsv -num_threads 32

MOUSE_COUNT=`cat $OUTPUTDIR/Prealign.tsv | grep 'Mouse_' | wc -l`
HUMAN_COUNT=`cat $OUTPUTDIR/Prealign.tsv | grep 'Human_' | wc -l`

#let MOUSE_COUNT=MOUSE_COUNT+1
#let HUMAN_COUNT=HUMAN_COUNT+1

#let HumanMouseRatio=`bc <<< "scale=5;HUMAN_COUNT/MOUSE_COUNT"`
#let MouseHumanRatio=HUMAN_COUNT/MOUSE_COUNT

if [ $MOUSE_COUNT -gt $HUMAN_COUNT ]
then
	echo 'Mouse'
elif [ $MOUSE_COUNT -lt $HUMAN_COUNT ]
then
	echo 'Human'
else
	echo 'Bad raw data'
fi
