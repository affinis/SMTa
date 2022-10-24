#!/bin/bash


HELP=`cat <<EOF
usage:
  bash install_SMTa.sh <opt> [args]

opts:
  -h	print help
  -n	absolute path of your nt database, eg: /data/Download/database/nt
  -s	absolute path of spaceranger, eg: /home/Foo/tools/spaceranger-2.0.0/spaceranger
  -b	absolute path of blastn, eg: /home/Foo/Download/ncbi-blast+/ncbi-blast-2.10.1+/bin/blastn
	
use 'which blastn' or 'which spaceranger' to get absolute path if you forget where they are.
EOF
`
if [ -z $1 ]
then
	echo "$HELP"
	exit
fi

while getopts hn:b:s: option
do
        case $option in
        h)
                echo "$HELP"
                exit
                ;;
        n)
                nt="$OPTARG"
                ;;
        b)
                blast="$OPTARG"
                ;;
	s)
		spaceranger="$OPTARG"
                ;;
        *)
                echo "$HELP"
                ;;
        esac
done

main_path=`realpath $0`
main_path=${main_path%/install_SMTa.sh}

sed -i "s|^MAINPATH=.*|MAINPATH=$main_path|g" src/configs.conf
sed -i "s|^NT=.*|NT=$nt|g" src/configs.conf
sed -i "s|^BLASTN=.*|BLASTN=$blast|g" src/configs.conf

if [ "`grep -o -m 1 'spaceranger-[0-9]' <<< $spaceranger`" = "spaceranger-1" ]
then
	sed -i "s|^SLIDE=.*|SLIDE='--unknown-slide'|g" src/configs.conf
	sed -i "s|^SPACERANGER=.*|SPACERANGER=$spaceranger|g" src/configs.conf
elif [ "`grep -o -m 1 'spaceranger-[0-9]' <<< $spaceranger`" = "spaceranger-2" ]
then
	sed -i "s|^SLIDE=.*|SLIDE='--unknown-slide visium-1'|g" src/configs.conf
	sed -i "s|^SPACERANGER=.*|SPACERANGER=$spaceranger|g" src/configs.conf
else
	echo "The spaceranegr path you provided is not valid, refer to -h for an example please"
	exit
fi
