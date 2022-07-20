import numpy as np
import pandas as pd
import sys
from pipetools import pipe
from functools import reduce

def calUMITaxa(df):
    #17 is the column of species
    species=list(np.unique(df[17]))
    if len(species)==1:
        #3 is the column of taxid
        taxid=list(df[3])[0]
    else:
        #17 is the column of species
        sums=pd.Series(df[17]).value_counts()
        if sums[0]==sums[1]:
            taxid="collision"
        else:
            #17 is the column of species, 3 is the column of taxid
            mostPreTaxa=sums.keys()[0]
            taxid=list(df[df[17]==mostPreTaxa][3])[0]
    return(taxid)

def calBarcodeTaxa(df):
    UMI=[]
    #1 is the column of barcode, 2 is the column of umi
    barcode=df[1][df.index[0]]
    barcode_reads=df
    n_reads=len(barcode_reads)
    umis=list(np.unique(barcode_reads[2]))
    n_umis=len(umis)
    for umi in umis:
#        print(umi)
        #2 is the column of UMI
        umi_reads=barcode_reads[barcode_reads[2]==umi]
        taxid=calUMITaxa(umi_reads)
#        print(taxid)
        UMI.append(taxid)
    occurences=pd.Series(UMI).value_counts()
    occurences=pd.DataFrame.from_dict(occurences,dtype='int64')
    occurences=occurences.rename(columns={0:str(barcode)})
    return(occurences)

try:
    df=pd.read_table('./possorted_genome_bam.merged.microbe.tsv',header=None)
#     df=pd.read_table('./possorted_genome_bam.merged.tsv',header=None)
except pd.errors.EmptyDataError:
    print("number of reads processed: 0")
    print("number of barcodes: 0")
    print("number of taxid: 0")
    print("number of valid microbe UMI: 0")
    print("number of microbe UMI collisions: 0")
    sys.exit()

#1 is the column of barcode
barcodes=list(np.unique(df[1]))
barcodes_counts=[]
for barcode in barcodes:
    #1 is the column of barcode
    barcode_count=calBarcodeTaxa(df[df[1]==barcode])
    barcodes_counts.append(barcode_count)
outdftotal=reduce(lambda  left,right: pd.merge(left,right,how='outer',left_index=True,right_index=True), barcodes_counts).fillna(0)

out=outdftotal[outdftotal.index!="collision"]
microbe_candidates_df=pd.read_table('./possorted_genome_bam.merged.microbe.tsv',header=None)
#3 is the column of taxid
microbe_candidates_taxids=list(np.unique(microbe_candidates_df[3]))
out=out[out.index.isin(microbe_candidates_taxids)]
num_of_reads=len(df)
num_of_candidate_microbe_reads=len(microbe_candidates_df)
num_of_umi=out.values.sum()
num_of_collision=outdftotal[outdftotal.index=="collision"].values.sum()
num_of_barcodes=len(out.keys())
num_of_taxid=len(out.index)
print("number of reads processed: "+str(num_of_reads))
print("number of candidate microbe reads: "+str(num_of_candidate_microbe_reads))
print("number of barcodes: "+str(num_of_barcodes))
print("number of taxid: "+str(num_of_taxid))
print("number of valid microbe UMI: "+str(num_of_umi))
print("number of microbe UMI collisions: "+str(num_of_collision))
print("top taxids: \n")
print(out.sum(axis=1).sort_values().tail(10))
out.to_csv('./matrix.mdx',sep="\t")
