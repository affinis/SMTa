#!/bin/bash

set -e
set -u

if [ -f ./unmapped.nonpolyB.blast.microbe.tsv ]
then
	rm -f ./unmapped.nonpolyB.blast.microbe.tsv ./unmapped.nonpolyB.blast.genus.list
fi


NT_BLASTDB_PATH="/data/LyuLin/Download/database"
MICROBE_GENUS_LIST="/data/LyuLin/Download/ncbi-taxon-names/MicroTaxidList"
THREADS=32
MAPPING="/data/LyuLin/Download/Mappings/taxid2taxon.tsv"
Rmerge=/data/LyuLin/Scripts/merge.R
FASTA=0

HELP=`cat <<EOF
    mip [opt] <args>
	-h	print help
	-i	input bam file
	-t	number of threads to use, default 32
	-u	umi availability, 1 or 0.
	-f	generate fasta only, if specified with 1.
EOF
`

while getopts hi:t:u:f: options
do
	case $options in
	h)
		echo "$HELP"
		exit 0
		;;
	i)
		BAM=$OPTARG
		;;
	t)
		THREADS=$OPTARG
		;;
	u)
		UMI=$OPTARG
		;;
	f)
		FASTA=$OPTARG
		;;
	*)
		echo "$HELP"
		exit 0
		;;
	esac
done

#extract unmmapped reads, delete reads with multiple ATCG(artificial reads) and custom to fasta
if [ ! -f ./${BAM%.bam}.unmapped.fasta ]
then
#	samtools view -@ $THREADS $BAM | awk 'match($0,/[ATCG]*-1/)&&$3=="*"{print ">"$1":"substr($0,RSTART,RLENGTH)"\n"$10}' | grep -E -v -B 1 "GGGGGGG*|CCCCCCC*|AAAAAAA*|TTTTTTT*|NNNNNNN*|^>" | sed '/^--$/d' > ${BAM%.bam}.unmapped.fasta
	if [ $UMI == "0" ]
	then
		samtools view -@ $THREADS $BAM | awk '$3=="*"{print ">"$1"|"b[0]"|"u[0]"\n"$10}' | grep -E -v -B 1 "GGGGGGGG*|CCCCCCCC*|AAAAAAAA*|TTTTTTTT*|NNNNNNNN*|^>" | sed '/^--$/d' > ${BAM%.bam}.unmapped.fasta
	else
		samtools view -@ $THREADS $BAM | awk 'match($0,/CB:Z:[ATCG]*-1/,b)&&match($0,/UB:Z:[ATCG]*/,u)&&$3=="*"{print ">"$1"|"b[0]"|"u[0]"\n"$10}' | grep -E -v -B 1 "GGGGGGGG*|CCCCCCCC*|AAAAAAAA*|TTTTTTTT*|NNNNNNNN*|^>" | sed '/^--$/d' > ${BAM%.bam}.unmapped.fasta
	fi
fi

if [ "$FASTA" == "1" ]
then
	exit
fi

#blast against nt database
blastn -query ${BAM%.bam}.unmapped.fasta -db $NT_BLASTDB_PATH/nt -qcov_hsp_perc 60 -perc_identity 80 -sorthits 4 -outfmt '6 qseqid staxid stitle length evalue pident qcovus' -max_target_seqs 1 -max_hsps 1 -culling_limit 1 -num_threads $THREADS > unmapped.nonpolyB.blast.tsv

#screen for microbe blast hits
#cat unmapped.nonpolyB.blast.tsv | awk -F $'\t' '$2!="10090"&&$2!="9606"&&$2!=9598&&$2!=10092{print $0}' > unmapped.nonpolyB.blast.cleaned.tsv
cat unmapped.nonpolyB.blast.tsv | cut -f 2 | qtaxa > ${BAM%.bam}.taxa.tsv
paste unmapped.nonpolyB.blast.tsv ${BAM%.bam}.taxa.tsv > ${BAM%.bam}.merged.tsv
#cat $MICROBE_GENUS_LIST | parallel -I {} awk\ \-F\ \$\'\\t\'\ \-v\ var\=\{\}\ \'\$2\=\=var\{print\ \$0\}\'\ unmapped.nonpolyB.blast.cleaned.tsv >> unmapped.nonpolyB.blast.microbe.tsv

\

if [ $UMI == "0" ]
then
	sed -i '1i\qseqid\ttaxid\tstitle\tlength\tevalue\tpident\tqcovus\ttaxid\tSuperkingdom\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' ${BAM%.bam}.merged.tsv
else
	sed -i 's/|[CU]B:Z:/\t/g' ${BAM%.bam}.merged.tsv
	sed -i '1i\qseqid\tbarcode\tUMI\ttaxid\tstitle\tlength\tevalue\tpident\tqcovus\ttaxid\tSuperkingdom\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' ${BAM%.bam}.merged.tsv
fi

cat ${BAM%.bam}.merged.tsv | grep -E 'kingdom_Fungi|superkingdom_Archaea|superkingdom_Bacteria|superkingdom_Viruses' > ${BAM%.bam}.merged.microbe.tsv.tmp
cat ${BAM%.bam}.merged.tsv | grep 'Eukaryota' | grep -v -E 'Metazoa|Viridiplantae|Fungi' >> ${BAM%.bam}.merged.microbe.tsv.tmp
cat ${BAM%.bam}.merged.microbe.tsv.tmp > ${BAM%.bam}.merged.microbe.tsv
rm ${BAM%.bam}.merged.microbe.tsv.tmp
#Rscript --vanilla $Rmerge -i unmapped.nonpolyB.blast.microbe.tsv -r $MAPPING -t tsv -f taxid

#genrate genus list
#cat $MICROBE_GENUS_LIST | parallel -I {} grep -w -o {} unmapped.nonpolyB.blast.tsv >> unmapped.nonpolyB.blast.genus.list
#cat unmapped.nonpolyB.blast.genus.list.0 | sort | uniq -c | sort -nr > unmapped.nonpolyB.blast.genus.list
