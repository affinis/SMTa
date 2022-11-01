#!/bin/bash

set -e
set -u
#set -o pipefail

SC_PATH=`realpath $0`
SC_PATH=${SC_PATH%/modname.sh}
MODNAME=$SC_PATH/modname.py

date >> modname.log
newname=`$MODNAME cif $1`

if [ -e ./$newname ]
then
	echo "File $newname was found already exsiting while rename $1 to $newname"
else
	mv $1 `$MODNAME cif $1`
	echo "$1 > `$MODNAME cif $1`"
	echo "$1 > `$MODNAME cif $1`" >> modname.log
fi
