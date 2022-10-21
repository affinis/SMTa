#!/bin/bash

set -e
set -u

echo "SMTa start"

M_PATH=`realpath $0`
source ${M_PATH%/main.sh}/configs.conf

if [ ! -z $A ]
then
	bash ${M_PATH%/main.sh}/spacewrapper.sh $A $A_image
else
	exit
fi

if [ ! -z $B ]
then
	bash ${M_PATH%/main.sh}/spacewrapper.sh $B $B_image
else
	exit
fi

if [ ! -z $C ]
then
	bash ${M_PATH%/main.sh}/spacewrapper.sh $C $C_image
else
	exit
fi

if [ ! -z $D ]
then
	bash ${M_PATH%/main.sh}/spacewrapper.sh $D $D_image
else
	exit
fi
