#!/bin/bash

set -e
set -u

cat possorted_genome_bam.merged.microbe.tsv | grep 'superkingdom_Bacteria' | cut -f 1 | parallel -I {} grep -w -A 1 {} possorted_genome_bam.unmapped.fasta >> bacterial.reads.fasta

blastn -query bacterial.reads.fasta -db /data/LyuLin/Ref/MergedReference/16s-23s_merge/16s.23s.merged -out rb.blast.out.tsv -outfmt '6 qaccver saccver pident length qstart qend sstart send evalue' -num_alignments 1 -qcov_hsp_perc 60 -perc_identity 90 -num_threads 96

bacterial_num=`cat bacterial.reads.fasta | grep '>' | wc -l`
n16s_num=`cat rb.blast.out.tsv | grep '16s_' | wc -l`
n23s_num=`cat rb.blast.out.tsv | grep '23s_' | wc -l`
n16s_rat=`bc <<< "scale=4; ${n16s_num}/${bacterial_num}"`
n23s_rat=`bc <<< "scale=4; ${n23s_num}/${bacterial_num}"`
ribosomal_rna_rat=`bc <<< "scale=4; ${n16s_rat}+${n23s_rat}"`

out=`cat << EOF
bacterial_read_num	${bacterial_num}
16s_read_num	${n16s_num}
23s_read_num	${n23s_num}
16s_read_rat	0${n16s_rat}
23s_read_rat	0${n23s_rat}
total_ribosomal_16s23s_rat	0${ribosomal_rna_rat}
EOF
`

echo "$out" > ribosomal_rna_summary.tsv
