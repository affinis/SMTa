#!/bin/bash

set -e
set -u
#set -o pipefail

MODNAME=/data/LyuLin/Scripts/Modname/bin/modname.py
MODE="cif"

date >> modname.log
for file in `ls ./ | grep -E '\.fastq\.gz$|\.fastq$|\.fq\.gz$|\.fq$'`
do
	if [ `echo $file | grep '_[RrIi][12]' | wc -l` -eq 0 ]
	then
		echo "No obvious read info like 'r1,i1,r2,i2' or 'R1,I1,R2,I2' were identified in file name, it will be guessed"
	fi
	newname=`$MODNAME $MODE $file`
	if [ "$newname" == "In valid input file name! Please check!" ]
	then
		echo "File $file $newname"
		break
	fi
	#echo $newname
	if [ -e ./$newname ]
	then
		echo "File $newname was found already exsiting while rename $file to $newname"
		continue
	else
		mv $file `$MODNAME $MODE $file`
		if [ -f $file.md5 ]
                then
                        mv $file.md5 `$MODNAME $MODE $file`.md5
                        echo "$file.md5 > `$MODNAME $MODE $file`.md5"
                        echo "$file.md5 > `$MODNAME $MODE $file`.md5" >> modname.log
                fi
		echo "$file > `$MODNAME $MODE $file`"
		echo "$file > `$MODNAME $MODE $file`" >> modname.log
	fi
done
