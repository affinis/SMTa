#!/bin/bash

set -e
set -u

SC_PATH=`realpath $0`
SC_PATH=${SC_PATH%/run_SMTa_test.sh}

SPR_PATH=`grep '^SPACERANGER' $SC_PATH/src/configs.conf | cut -d '=' -f 2`
SPR_PATH=${SPR_PATH%/spaceranger}

echo "spaceranger path $SPR_PATH"

echo "testing SMTa via sample file......"

sed -i "s|^SAMPLE_PATH.*|SAMPLE_PATH=$SPR_PATH/external/spaceranger_tiny_inputs/fastqs|g" $SC_PATH/samples.info
sed -i -e "s|^A=.*|A=tinytest|g" -e "s|^A_image=.*|A_image=$SPR_PATH/external/spaceranger_tiny_inputs/image/tinyimage.jpg|g" $SC_PATH/samples.info

if [ "`grep -o -m 1 'spaceranger-[0-9]' <<< $SPR_PATH`" = "spaceranger-1" ]
then
	echo "Using spaceranger 1.0.0 ref"
	sed -i "s|^REF=.*|REF=$SPR_PATH/external/spaceranger_tiny_ref/1.0.0|g" $SC_PATH/samples.info

elif [ "`grep -o -m 1 'spaceranger-[0-9]' <<< $SPR_PATH`" = "spaceranger-2" ]
then
	echo "Using spaceranger 2.0.0 ref"
	sed -i "s|^REF.*|REF=$SPR_PATH/external/spaceranger_tiny_ref|g" $SC_PATH/samples.info
else
	echo "tiny reference of spaceranger not found, you may provide an incorrect spaceranger path while installing"
fi

bash src/main.sh $SC_PATH/samples.info $SC_PATH/src/configs.conf

sed -i "s|^SAMPLE_PATH=.*|SAMPLE_PATH=\"\"|g" $SC_PATH/samples.info
sed -i -e "s|^A=.*|A=\"\"|g" -e "s|^A_image=.*|A_image=\"\"|g" $SC_PATH/samples.info
sed -i "s|^REF=.*|REF=\"\"|g" $SC_PATH/samples.info
