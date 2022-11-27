#!/bin/bash

#set -e
#set -u
#set -o pipefail

SELF_PATH=`realpath $0`
. $3
. $4

SAMPLE=$1
QTAXA="`dirname $SELF`/querytaxa.py `dirname $SELF`"

if [ -d $INPUT_DIR/$SAMPLE ]
then
	rm -drf $INPUT_DIR/$SAMPLE
fi

if [ -d $OUTPUT_DIR/$SAMPLE ]
then
        rm -drf $OUTPUT_DIR/$SAMPLE
fi

echo "create dir $INPUT_DIR/$SAMPLE" 
mkdir -p $INPUT_DIR/$SAMPLE
echo "create dir $OUTPUT_DIR/$SAMPLE"
mkdir -p $OUTPUT_DIR/$SAMPLE

find $SAMPLE_PATH -name "*$SAMPLE*" | xargs -I {} ln -s {} $INPUT_DIR/$SAMPLE/
cd $INPUT_DIR/$SAMPLE/
$MODALLNAME

cd $SELF

REF_APP=$REF

echo `date` >> $SELF/run.info
echo "New sample started!" >> $SELF/run.info
echo -e "spaceranger arguments set:\n$SPACERANGER count --id=$SAMPLE --fastqs=$INPUT_DIR/$SAMPLE/ --sample=$SAMPLE --transcriptome=$REF_APP --localcores=$CORE --localmem=$MEM $SLIDE --image=$2" >> $SELF/run.info

cd $OUTPUT_DIR/$SAMPLE

$SPACERANGER count --id=$SAMPLE --fastqs=$INPUT_DIR/$SAMPLE/ --sample=$SAMPLE --transcriptome=$REF_APP --localcores=$CORE --localmem=$MEM $SLIDE --image=$2

if [ "$?" == 0 ]
then
	echo `date` >> $SELF/run.info
	echo "Successful" >> $SELF/run.info
else
	echo `date` >> $SELF/run.info
	echo "Failed" >> $SELF/run.info
fi

cd $SAMPLE/outs/

if [ "$MODE" == "SMT" ]
then
	$MIPMAIN -i possorted_genome_bam.bam -u 1 -f 0 -t $CORE -d $NT
	python3 $GENMATRIX > SMT_summary.txt
elif [ "$MODE" == "SMT.MRCA" ]
then
	$MIPMAIN -i possorted_genome_bam.bam -u 1 -f 0 -t $CORE -d $NT
	python3 ${GENMATRIX%.py}.MRCA.py > SMT_summary.txt
	cat taxaids.lst | $QTAXA > taxa.tsv
fi

#$PARSER

exit
#echo -e "spaceranger arguments set:\n$SPACERANGER count --id=$SAMPLE --fastqs=$INPUT_DIR/$SAMPLE/ --sample=$SAMPLE --transcriptome=$REF_APP --localcores=$CORE --localmem=$MEM --unknown-slide --image=$2" >> run.info
