#!/bin/bash

set -e
set -u

HELP=`cat <<EOF
        -h      print help.
        -d      absolute path of sample directory.
	-s	sample ID (first field of fastq filename).
	-p	absolute path of image file.
	-r	absolute reference path.
	-c      provide sample info file absolute path instead of pass by arguments
EOF
`

while getopts hd:s:p:r:c: options
do
	case $options in
	h)
		echo "$HELP"
		exit 0
		;;
	d)
		SAMPLE_DIR=$OPTARG
		;;
	s)
		SAMPLE=$OPTARG
		;;
	p)
		IMAGE=$OPTARG
		;;
	r)
		REF=$OPTARG
		;;
	c)
		SAMPLE_CONF=$OPTARG
		;;
	*)
		echo "$HELP"
		exit 0
		;;
	esac
done


SC_PATH=`realpath $0`
SC_PATH=${SC_PATH%/run_SMTa.sh}

if [ -z $SAMPLE_CONF ]
then
	sed -i "s|^SAPMLES=.*|SAPMLES=$SAMPLE_DIR|g" $SC_PATH/samples.info
	sed -i -e "s|^A=.*|A=tinytest|g" -e "s|^A_image=.*|A_image=|g" $SC_PATH/samples.info
	sed -i "s|^Reference=.*|Reference=$SPR_PATH/external/spaceranger_tiny_ref/1.0.0|g" $SC_PATH/samples.info
	bash $SC_PATH/src/main.sh $SC_PATH/samples.info $SC_PATH/src/configs.conf
else
	bash $SC_PATH/src/main.sh $SAMPLE_CONF $SC_PATH/src/configs.conf
fi
