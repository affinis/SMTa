#!/home/vin/anaconda3/bin/python

import sys

#query=sys.argv[1]

taxid2higherTaxidFile=open(sys.argv[1]+'/taxid2higherTaxid.tsv','r')
lines=taxid2higherTaxidFile.readlines()
queryTaxon={}
#print("Generating higher taxa hash.....")
for line in lines:
	line=line.rstrip()
	taxlow=line.split('\t')[0]
	taxhigh=line.split('\t')[1]
	queryTaxon[taxlow]=taxhigh
#print("done")
taxid2taxnameFile=open(sys.argv[1]+'/taxid2taxname.tsv','r')
lines=taxid2taxnameFile.readlines()
queryTaxid={}

#print("Generating taxid-taxname hash.....")
for line in lines:
	line=line.rstrip()
	taxid=line.split('\t')[0]
	taxname=line.split('\t')[1]
	queryTaxid[taxid]=taxname
#print("done")

taxid2taxclassFile=open(sys.argv[1]+'/taxid2taxclass.tsv','r')
lines=taxid2taxclassFile.readlines()
taxidClass={}
for line in lines:
	line=line.rstrip()
	taxid=line.split('\t')[0]
	taxclass=line.split('\t')[1]
	taxidClass[taxid]=taxclass

#print("Start query..")
#line="taxid:"+query+"\t"+taxidClass[query]+":"+queryTaxid[query]


#handle=open('taxa.matrix.tsv','w')
#handle.write("taxid\tsuperkingdom\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n")
for query in sys.stdin:
	query=query.rstrip()
	taxinfo={"taxid":query,"species":"","genus":"","family":"","order":"","class":"","phylum":"","kingdom":"","superkingdom":""}
	while query in queryTaxon.keys() and taxidClass[query]!="no rank":
		if taxidClass[query] in ["species","genus","family","order","class","phylum","kingdom","superkingdom"]:
			taxinfo[taxidClass[query]]=queryTaxid[query]
			query=queryTaxon[query]
		else:
			query=queryTaxon[query]
	reply=""
	for key in ["taxid","superkingdom","kingdom","phylum","class","order","family","genus","species"]:
		reply=reply+key+"_"+taxinfo[key].replace(" ","_")+"\t"
		reply=reply.replace("taxid_","")
	reply=reply.rstrip('\t')
	print(reply)
#	handle.write(reply+"\n")
#handle.close()

