#!/usr/bin/python3

import copy

def fullname(FileNameObject):
	return(FileNameObject.basename+FileNameObject.suffix)

def basename(FileNameObject):
	return(FileNameObject.basename)

def suffix(FileNameObject):
	return(FileNameObject.suffix)

class FileName():
	def __init__(self,basename="basename",suffix=".suffix"):
		self.basename=basename
		self.suffix=suffix
		self.fullname=self.basename+self.suffix

	def getfullname(self):
		print(self.fullname)

	def getbasename(self):
		print(self.basename)

class CellRangerInputFileNameCompressed(FileName):
	
	def __init__(self,sampleID="UnknownID",sample="S1",lane="L001",read="R1",NO="001"):
		self.sampleID=sampleID
		self.sample=sample
		self.lane=lane
		self.read=read
		self.NO=NO
		self.basename=self.sampleID+"_"+self.sample+"_"+self.lane+"_"+self.read+"_"+self.NO
		self.suffix=".fastq.gz"
		self.fullname=self.basename+self.suffix

	def stdbasename(self):
		new=copy.deepcopy(self)
		if self.sampleID=="" or self.sampleID=="_":
			new.sampleID="UnknownID"
		if self.sample=="" or self.sample=="_":
			new.sample="S1"
		if self.lane=="" or self.lane=="_":
			new.lane="L001"
		if self.read=="" or self.read=="_":
			new.read="R1"
		if self.NO=="" or self.NO=="_":
			new.NO="001"
		new.basename=new.sampleID+"_"+new.sample+"_"+new.lane+"_"+new.read+"_"+new.NO
		new.fullname=new.basename+new.suffix
		return(new)

class CellRangerInputFileName(CellRangerInputFileNameCompressed):

	def __init__(self,sampleID="UnknownID",sample="S1",lane="L001",read="R1",NO="001"):
		CellRangerInputFileNameCompressed.__init__(self,sampleID,sample,lane,read,NO)
		self.suffix=".fastq"
		self.fullname=self.basename+self.suffix
