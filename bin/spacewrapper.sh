#!/bin/bash

#set -e
#set -u
#set -o pipefail

source ./configs.conf

SAMPLE=$1

if [ -d $INPUT_DIR/$SAMPLE ]
then
	rm -drf $INPUT_DIR/$SAMPLE
fi

if [ -d $OUTPUT_DIR/$SAMPLE ]
then
        rm -drf $OUTPUT_DIR/$SAMPLE
fi
 
mkdir -p $INPUT_DIR/$SAMPLE
mkdir -p $OUTPUT_DIR/$SAMPLE

find $SAMPLE_PATH -name "*$SAMPLE*" | xargs -I {} ln -s {} $INPUT_DIR/$SAMPLE/
cd $INPUT_DIR/$SAMPLE/
$MODALLNAME

cd $SELF

REF_ANIMAL=`$SUBSAMPLE $SAMPLE $INPUT_DIR/$SAMPLE/*_R1_*.fastq.gz $OUTPUT_DIR/$SAMPLE`

if [ "$REF_ANIMAL" == "Human" ]
then
	REF_APP=$REF
elif [ "$REF_ANIMAL" == "Mouse" ]
then
	REF_APP=$REF2
else
	echo "$REF_ANIMAL, specify reference manually[hs/mm]:"
	read
	if [ $REPLY == "hs" ]
	then
		REF_APP=$REF
	elif [ $REPLY == "mm" ]
	then
		REF_APP=$REF2
	else
		echo "Invalid reference, stopping..."
	fi
fi
echo `date` >> $SELF/run.info
echo "New sample started!" >> $SELF/run.info
echo -e "spaceranger arguments set:\n$SPACERANGER count --id=$SAMPLE --fastqs=$INPUT_DIR/$SAMPLE/ --sample=$SAMPLE --transcriptome=$REF_APP --localcores=$CORE --localmem=$MEM --unknown-slide --image=$2" >> $SELF/run.info

cd $OUTPUT_DIR/$SAMPLE

$SPACERANGER count --id=$SAMPLE --fastqs=$INPUT_DIR/$SAMPLE/ --sample=$SAMPLE --transcriptome=$REF_APP --localcores=$CORE --localmem=$MEM --unknown-slide --image=$2

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
	$MIPMAIN -i possorted_genome_bam.bam -u 1 -f 0 -t $CORE
	python3 $GENMATRIX > SMT_summary.txt
fi

$PARSER

exit
#echo -e "spaceranger arguments set:\n$SPACERANGER count --id=$SAMPLE --fastqs=$INPUT_DIR/$SAMPLE/ --sample=$SAMPLE --transcriptome=$REF_APP --localcores=$CORE --localmem=$MEM --unknown-slide --image=$2" >> run.info
