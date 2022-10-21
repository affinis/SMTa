#!/usr/bin/python3

import sys,re
sys.path.append(r"/data/LyuLin/Scripts/Modname/bin")
from FileName import *

option=sys.argv[1]

if option=='cif':
	name_to_mod=sys.argv[2]
	elements=re.split(r'[._]',name_to_mod)
	sampleid=elements[0]
	sample="_"
	lane="_"
	read="_"
	NO="_"
	
	for element in elements:
		if re.fullmatch(r'[Ss]{1,1}[\d]{1,2}',element):
			sample=element.upper()
		elif re.fullmatch(r'[Ll]{1,1}[\d]{3,3}',element):
			lane=element.upper()
		elif re.fullmatch(r'[RIri]{1,1}[12]{1,1}',element):
			read=element.upper()
		elif re.fullmatch(r'[\d]{3,3}',element):
			NO=element
		else:
			continue
	if read=="_":
#		print("No obvious read info like 'r1,i1,r2,i2' or 'R1,I1,R2,I2' were identified in file name, guessed")
		for element in elements:
			if re.fullmatch(r'[12]{1,1}',element):
				read=str("R"+element)
		if read=="_":
			print("In valid input file name! Please check!")
			exit()

	if (elements[-1]=="gz" and elements[-2]=="fastq") or (elements[-1]=="gz" and elements[-2]=="fq"):
		protoNameOb=CellRangerInputFileNameCompressed(sampleID=sampleid,sample=sample,lane=lane,read=read,NO=NO)
	elif elements[-1]=="fastq" or elements[-1]=="fq":
		protoNameOb=CellRangerInputFileName(sampleID=sampleid,sample=sample,lane=lane,read=read,NO=NO)
	else:
		print("In valid input file name! Please check!")
		exit()
	modNameOb=protoNameOb.stdbasename()
#	print(name_to_mod+" > "+protoNameOb.fullname+" > "+modNameOb.fullname)
	print(modNameOb.fullname)
