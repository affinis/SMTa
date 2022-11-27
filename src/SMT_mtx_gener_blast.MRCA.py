import numpy as np
import pandas as pd
import os
import sys
from pipetools import pipe
from functools import reduce

def calUMITaxa(df):
    #get species list of this umi
    #17 is the column of species
    species=list(np.unique(df[17]))
    #get counts of these species
    sums=pd.Series(df[17]).value_counts()
    if len(species)==1 or sums[0]/sums[1]>=2:
        #3 is the column of taxid
        taxname=sums.keys()[0]
        if taxname=="species_":
          taxname="collision"
    else:
	#17 is the column of species
        current_taxa_level=17
        candidate1=sums.keys()[0].replace("species_","").replace("_"," ")
        candidate2=sums.keys()[1].replace("species_","").replace("_"," ")
        if candidate1=="" and candidate2!="":
          taxname=sums.keys()[1]
          return(taxname)
        if candidate2=="" and candidate1!="":
          taxname=sums.keys()[0]
          return(taxname)
        candidate1_id=name2id[candidate1]
        candidate2_id=name2id[candidate2]
        while candidate1_id!=candidate2_id and current_taxa_level!=11:
          candidate1_id=queryTaxon[candidate1_id]
          candidate2_id=queryTaxon[candidate2_id]
          current_taxa_level-=1
        if current_taxa_level==11:
          taxname="collision"
        else:
          taxname=id2name[candidate1_id]
          if current_taxa_level==16:
            taxname="genus_"+taxname
          elif current_taxa_level==15:
            taxname="family_"+taxname
          elif current_taxa_level==14:
            taxname="order_"+taxname
          elif current_taxa_level==13:
            taxname="class_"+taxname
          else:
            taxname="phylum_"+taxname
    return(taxname)

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

#generate taxa id to higher taxa id mapping
taxid2higherTaxidFile=open(os.path.dirname(sys.argv[0])+'/taxid2higherTaxid.tsv','r')
lines=taxid2higherTaxidFile.readlines()
queryTaxon={}
#print("Generating higher taxa hash.....")
for line in lines:
        line=line.rstrip()
        taxlow=line.split('\t')[0]
        taxhigh=line.split('\t')[1]
        queryTaxon[taxlow]=taxhigh

taxaName2taxaidFile=open(os.path.dirname(sys.argv[0])+'/taxid2taxname.tsv','r')
lines=taxaName2taxaidFile.readlines()
name2id={}
id2name={}
#print("Generating higher taxa hash.....")
for line in lines:
        line=line.rstrip()
        taxaName=line.split('\t')[1]
        taxaid=line.split('\t')[0]
        id2name[taxaid]=taxaName
        name2id[taxaName]=taxaid

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
i=0
new_index=[]
while i<len(out):
    new_index.append(name2id[out.index[i].replace("species_","").replace("genus_","").replace("family_","").replace("order_","").replace("class_","").replace("phylum_","").replace("_"," ")])
    i+=1
out.index=new_index
microbe_candidates_df=pd.read_table('./possorted_genome_bam.merged.microbe.tsv',header=None)
#3 is the column of taxid
#microbe_candidates_taxids=list(np.unique(microbe_candidates_df[3]))
#out=out[out.index.isin(microbe_candidates_taxids)]
num_of_reads=len(df)
num_of_candidate_microbe_reads=len(microbe_candidates_df)
num_of_umi=out.values.sum()
num_of_collision=outdftotal[outdftotal.index=="collision"].values.sum()
num_of_barcodes=len(out.keys())
num_of_taxa=len(out.index)
print("number of reads processed: "+str(num_of_reads))
print("number of candidate microbe reads: "+str(num_of_candidate_microbe_reads))
print("number of barcodes: "+str(num_of_barcodes))
print("number of taxid: "+str(num_of_taxa))
print("number of valid microbe UMI: "+str(num_of_umi))
print("number of microbe UMI collisions: "+str(num_of_collision))
print("top taxids: \n")
print(out.sum(axis=1).sort_values().tail(10))
out.to_csv('./matrix.mdx',sep="\t")
taxaidlist=open("./taxaids.lst","w")
for taxium in out.index:
#  taxium=taxium.replace("species_","").replace("genus_","").replace("family_","").replace("order_","").replace("class_","").replace("phylum_","").replace("_"," ")
#  taxaidlist.write(name2id[taxium]+"\n")
  taxaidlist.write(taxium+"\n")
taxaidlist.close()
