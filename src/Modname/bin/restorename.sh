#!/bin/bash

set -e
set -u
#set -o pipefail


date >> modname.log
while read line
do
	historyname=`echo "$line" | awk 'BEGIN{FS=" > "}{print $1}'`
	currentname=`echo "$line" | awk 'BEGIN{FS=" > "}{print $2}'`
	if [ ! -e ./$currentname ]
	then
		echo "File $currentname not exist, skipping..."
		continue
	fi
	if [ -e ./$historyname ]
	then
		echo "File $historyname was found existing while convert current file $currentname to $historyname"
		continue
	else
		mv $currentname $historyname
		echo "$historyname < $currentname"
		echo "$historyname < $currentname" >> modname.log
	fi
done < <(cat modname.log | grep '>' | sort | uniq)
