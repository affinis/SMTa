#!/bin/bash

set -e
set -u

echo "SMTa start"

M_PATH=`realpath $0`
source $1 $2

if [ ! -z $A ]
then
	bash ${M_PATH%/main.sh}/spacewrapper.sh $A $A_image $1 $2
else
	exit
fi

if [ ! -z $B ]
then
	bash ${M_PATH%/main.sh}/spacewrapper.sh $B $B_image $1 $2
else
	exit
fi

if [ ! -z $C ]
then
	bash ${M_PATH%/main.sh}/spacewrapper.sh $C $C_image $1 $2
else
	exit
fi

if [ ! -z $D ]
then
	bash ${M_PATH%/main.sh}/spacewrapper.sh $D $D_image $1 $2
else
	exit
fi
