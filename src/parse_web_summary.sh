#!/bin/bash

set -e
set -u

MeanReadsperSpotheader=`cat web_summary.html | grep -o '\["Mean Reads per Spot.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | head -1`
MeanReadsperSpotData=`cat web_summary.html | grep -o '\["Mean Reads per Spot.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | tail -1`
MedianGenesperSpotheader=`cat web_summary.html | grep -o '\["Median Genes per Spot.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | head -1`
MedianGenesperSpotData=`cat web_summary.html | grep -o '\["Median Genes per Spot.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | tail -1`
NumberofReadsheader=`cat web_summary.html | grep -o '\["Number of Reads.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | head -1`
NumberofReadsData=`cat web_summary.html | grep -o '\["Number of Reads.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | tail -1`
ValidBarcodesheader=`cat web_summary.html | grep -o '\["Valid Barcodes.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | head -1`
ValidBarcodesData=`cat web_summary.html | grep -o '\["Valid Barcodes.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | tail -1`
SequencingSaturationheader=`cat web_summary.html | grep -o '\["Sequencing Saturation.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | head -1`
SequencingSaturationData=`cat web_summary.html | grep -o '\["Sequencing Saturation.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | tail -1`
Q30BasesinBarcodeheader=`cat web_summary.html | grep -o '\["Q30 Bases in Barcode.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | head -1`
Q30BasesinBarcodeData=`cat web_summary.html | grep -o '\["Q30 Bases in Barcode.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | tail -1`
Q30BasesinRNAReadheader=`cat web_summary.html | grep -o '\["Q30 Bases in RNA Read.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | head -1`
Q30BasesinRNAReadData=`cat web_summary.html | grep -o '\["Q30 Bases in RNA Read.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | tail -1`
Q30BasesinUMIheader=`cat web_summary.html | grep -o '\["Q30 Bases in UMI.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | head -1`
Q30BasesinUMIData=`cat web_summary.html | grep -o '\["Q30 Bases in UMI.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | tail -1`
ReadsMappedtoGenomeheader=`cat web_summary.html | grep -o '\["Reads Mapped to Genome.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | head -1`
ReadsMappedtoGenomeData=`cat web_summary.html | grep -o '\["Reads Mapped to Genome.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | tail -1`
ReadsMappedConfidentlytoGenomeheader=`cat web_summary.html | grep -o '\["Reads Mapped Confidently to Genome.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | head -1`
ReadsMappedConfidentlytoGenomeData=`cat web_summary.html | grep -o '\["Reads Mapped Confidently to Genome.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | tail -1`
ReadsMappedConfidentlytoIntergenicRegionsheader=`cat web_summary.html | grep -o '\["Reads Mapped Confidently to Intergenic Regions.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | head -1`
ReadsMappedConfidentlytoIntergenicRegionsData=`cat web_summary.html | grep -o '\["Reads Mapped Confidently to Intergenic Regions.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | tail -1`
ReadsMappedConfidentlytoIntronicRegionsheader=`cat web_summary.html | grep -o '\["Reads Mapped Confidently to Intronic Regions.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | head -1`
ReadsMappedConfidentlytoIntronicRegionsData=`cat web_summary.html | grep -o '\["Reads Mapped Confidently to Intronic Regions.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | tail -1`
ReadsMappedConfidentlytoExonicRegionsheader=`cat web_summary.html | grep -o '\["Reads Mapped Confidently to Exonic Regions.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | head -1`
ReadsMappedConfidentlytoExonicRegionsData=`cat web_summary.html | grep -o '\["Reads Mapped Confidently to Exonic Regions.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | tail -1`
ReadsMappedConfidentlytoTranscriptomeheader=`cat web_summary.html | grep -o '\["Reads Mapped Confidently to Transcriptome.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | head -1`
ReadsMappedConfidentlytoTranscriptomeData=`cat web_summary.html | grep -o '\["Reads Mapped Confidently to Transcriptome.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | tail -1`
ReadsMappedAntisensetoGeneheader=`cat web_summary.html | grep -o '\["Reads Mapped Antisense to Gene.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | head -1`
ReadsMappedAntisensetoGeneData=`cat web_summary.html | grep -o '\["Reads Mapped Antisense to Gene.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | tail -1`
FractionReadsinSpotsUnderTissueheader=`cat web_summary.html | grep -o '\["Fraction Reads in Spots Under Tissue.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | head -1`
FractionReadsinSpotsUnderTissueData=`cat web_summary.html | grep -o '\["Fraction Reads in Spots Under Tissue.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | tail -1`
TotalGenesDetectedheader=`cat web_summary.html | grep -o '\["Total Genes Detected.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | head -1`
TotalGenesDetectedData=`cat web_summary.html | grep -o '\["Total Genes Detected.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | tail -1`
MedianUMICountsperSpotheader=`cat web_summary.html | grep -o '\["Median UMI Counts per Spot.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | head -1`
MedianUMICountsperSpotData=`cat web_summary.html | grep -o '\["Median UMI Counts per Spot.[^]]*"\]' | head -1 | sed 's/"//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/, /\n/g' | sed 's/,//g' | tail -1`

SUMMARY=`cat <<EOF
$NumberofReadsheader: $NumberofReadsData
$ReadsMappedtoGenomeheader: $ReadsMappedtoGenomeData
$ReadsMappedConfidentlytoGenomeheader: $ReadsMappedConfidentlytoGenomeData
$ReadsMappedConfidentlytoIntergenicRegionsheader: $ReadsMappedConfidentlytoIntergenicRegionsData
$ReadsMappedConfidentlytoIntronicRegionsheader: $ReadsMappedConfidentlytoIntronicRegionsData
$ReadsMappedConfidentlytoExonicRegionsheader: $ReadsMappedConfidentlytoExonicRegionsData
$ReadsMappedConfidentlytoTranscriptomeheader: $ReadsMappedConfidentlytoTranscriptomeData
$MeanReadsperSpotheader: $MeanReadsperSpotData
$MedianGenesperSpotheader: $MedianGenesperSpotData
$MedianUMICountsperSpotheader: $MedianUMICountsperSpotData
$ValidBarcodesheader: $ValidBarcodesData
$SequencingSaturationheader: $SequencingSaturationData
$Q30BasesinBarcodeheader: $Q30BasesinBarcodeData
$Q30BasesinRNAReadheader: $Q30BasesinRNAReadData
$Q30BasesinUMIheader: $Q30BasesinUMIData
$ReadsMappedAntisensetoGeneheader: $ReadsMappedAntisensetoGeneData
$FractionReadsinSpotsUnderTissueheader: $FractionReadsinSpotsUnderTissueData
$TotalGenesDetectedheader: $TotalGenesDetectedData
EOF`

echo "$SUMMARY" > websummary.par

if [ -f ./SMT_summary.txt ]
then
        MicrobialReadsHeader="Number of Microbial Reads"
        MicrobialReadsData=`cat SMT_summary.txt | grep 'number of reads processed' | grep -o [0-9]*`
	echo "Ratio of Microbial Reads: `bc -l <<< "scale=8;$MicrobialReadsData/$NumberofReadsData"`" >> websummary.par
	cat SMT_summary.txt >> websummary.par
fi

cat websummary.par > summary.txt
rm websummary.par

