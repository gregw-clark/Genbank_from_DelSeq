#!/usr/bin/env python

import string, re, sys
from cPickle import load,dump
import multiprocessing
import itertools
import sys
import os
from collections import defaultdict
from itertools import izip
import glob
import re
from Bio.Seq import Seq
from operator import itemgetter
from APE import APEobj,_Exon,_Feature,_Other
import requests


def reverseSeq(aseq):
        lclSeq=Seq(aseq)
        revSeq=str(lclSeq.reverse_complement())
        return revSeq


def AnnotationRange(file):
	coords=filter(lambda k: k.startswith(">"),open(file,'r').readlines())
	glb=[]
	
	for c in coords:
		linedata=c.strip().split(":")[1].split("-")
		chromosome=c.strip().split("(")[0].lstrip(">chr")
	
		if len(linedata) == 2:
			c1,c2=int(linedata[0]),int(linedata[1])
			glb+=[c1,c2]
		elif len(linedata) == 1:
			c1=int(linedata[0])
			glb.append(c1)
		else:
			print c
			print "Exiting,no coordinates"
			sys.exit()
	minR=min(glb)-1000
	maxR=max(glb)+1000
	

	location=chromosome+":"+str(minR)+"-"+str(maxR)
	server= "https://rest.ensembl.org"
	ext = "/overlap/region/mouse/"+location+"?feature=transcript;feature=exon"

	try:
		r = requests.get(server+ext, headers={ "Content-Type" : "text/x-gff3"})
	
		if not r.ok:
			r.raise_for_status()
			sys.exit()

		gffdata=filter(lambda l: not l.startswith("#") and len(l),r.text.split("\n"))
	except requests.exceptions.HTTPError:
		print "WEB ERROR",
		print ext
		gffdata=[]
	
	return gffdata


def TargetedRange(file):
	coords=filter(lambda k: k.startswith(">"),open(file,'r').readlines())
	glb=[]
	server= "https://rest.ensembl.org"
	targets={}	
	for c in coords:
		linedata=c.strip().split(":")[1].split("-")
		chromosome=c.strip().split("(")[0].lstrip(">")
		if chromosome == "chr1":
			chromosome="1"
		if len(linedata) == 2:
			c1,c2=int(linedata[0]),int(linedata[1])
			location=chromosome+":"+str(c1)+"-"+str(c2)
		elif len(linedata) == 1:
			c1=int(linedata[0])
			location=chromosome+":"+str(c1)+"-"+str(c1)
		else:
			print c
			print "Exiting,no coordinates"
			sys.exit()

		ext = "/overlap/region/mouse/"+location+"?feature=gene;feature=CDS;feature=transcript;feature=exon"

		try:
			r = requests.get(server+ext, headers={ "Content-Type" : "text/x-gff3"})
		
			if not r.ok:
				r.raise_for_status()
				sys.exit()

			gffdata=filter(lambda l: not l.startswith("#") and len(l),r.text.split("\n"))
		except requests.exceptions.HTTPError:
			print "WEB ERROR",
			print ext
			gffdata=[]
	
		targets[c.lstrip(">").strip()]=gffdata
	return targets

def gff_process(gff):
	RExon=re.compile("exon_id=ENSMUSE[0-9]{3,15};")
	exons={}
	exonLine=filter(lambda p: RExon.search(p),gff)
	altFeatures=filter(lambda p: not RExon.search(p),gff)
	features={}
	for line in altFeatures:
		data=line.strip().split("\t")
		chromosome,anotBld,anotType,start,end=data[:5]
		strand=data[6]
		infoSp=data[8].strip().split(";")
		#if re.search("Parent=gene",data[8]):
		for sp in infoSp:
			if sp.startswith("Parent=gene"):
				genename=sp.split("gene:")[1]
				symbol=genename
		#print anotBld,anotType,data
		try:
			feature=filter(lambda l : re.search("ID=",l),infoSp)[0].split("=")[1]
			features[feature]=[feature.split(":")[1],chromosome,int(start),int(end),strand]
		except IndexError:
			pass
	
	for line in exonLine:
		data=line.strip().split("\t")
		chromosome,anotBld,anotType,start,end=data[:5]
		strand=data[6]
	
		exon=filter(lambda p: p.startswith("Name="),data[8].split(";"))[0].lstrip("Name=")
		rank=filter(lambda p: p.startswith("rank="),data[8].split(";"))[0].lstrip("rank=")
		phase=filter(lambda p: p.startswith("ensembl_phase="),data[8].split(";"))[0].lstrip("ensembl_phase=")
		end_phase=filter(lambda p: p.startswith("ensembl_end_phase="),data[8].split(";"))[0].lstrip("ensembl_end_phase=")
		transcript=filter(lambda p: p.startswith("Parent=transcript:"),data[8].split(";"))[0].lstrip("Parent=transcript:")
		if int(rank) < 3:
			if exon in exons:
				exons[exon][-1]+=[transcript]
			else:
				exons[exon]=[exon,chromosome,int(start),int(end),phase,end_phase,strand,[transcript]]

	return exons,features,genename


def processCoords(delfile):

	coords=filter(lambda d: d.startswith(">"),open(delfile,'r').readlines())
	positions=[]
	
	for c in coords:
		chromosome=c.strip().split(":")[0].lstrip(">chr").split("(")[0]
		posts=c.strip().split(":")[1].split("-")
		strand=c.strip().split("(")[1].split(")")[0]
		if len(posts) == 2:
			start,end=int(posts[0]),int(posts[1])
			positions+=[start,end]
		else:
			start=int(posts[0])
			positions.append(start)
	minP=min(positions)
	maxP=max(positions)
	server="http://rest.ensembl.org/sequence/region/mouse/"
	print "MIN" ,min(positions)
	print "MAX", max(positions)

	minP=101469587
	maxP=101468438

	if strand == "-":
	

		extFivePrime=str(chromosome)+":"+str(minP-1000)+".."+str(minP)+"?"+"content-type=text/plain;mask=soft"
		extThreePrime=str(chromosome)+":"+str(maxP)+".."+str(maxP+1000)+"?"+"content-type=text/plain;mask=soft"
	else:
		extFivePrime=str(chromosome)+":"+str(minP-1000)+".."+str(minP)+"?"+"content-type=text/plain;mask=soft"
		extThreePrime=str(chromosome)+":"+str(maxP)+".."+str(maxP+1000)+"?"+"content-type=text/plain;mask=soft"
	
	rFive = requests.get(server+extFivePrime)
	rThree = requests.get(server+extThreePrime)

	seqFiveP=rFive.text.strip()
	if seqFiveP.startswith("chrmsmGRCm"):
		seqFiveP=seqFiveP.lstrip("chrmsmGRCm")
	seqThreeP=rThree.text.strip()
	if seqThreeP.startswith("chrmsmGRCm"):
		seqThreeP=seqThreeP.lstrip("chrmsmGRCm")
	if strand == "-":
		print ">5'"
		print seqFiveP
		seqFiveP=str(Seq(seqFiveP).reverse_complement())
		print ">REV 5'"
		print seqFiveP,"\n\n"
		print ">3'"
		print seqThreeP
		seqThreeP=str(Seq(seqThreeP).reverse_complement())
		print ">REV 3'"
		print seqThreeP,"\n\n"
	return [[seqFiveP,minP-1000,minP],[seqThreeP,maxP,maxP+1000]]

delfiles=glob.glob("*.del")
comp=0
for file in delfiles:


	Mexons=[]
	fileinfo=file.split(".")[0]
	
	print fileinfo	
	FiveP,ThreeP=processCoords(file)


	initSeq=open(fileinfo+".txt",'r').readlines()
	#initSeq=open(fileinfo+".txt",'r').readlines()
	if len(initSeq) != 2:
		print "What is going on here"
		print fileinfo
		sys.exit()
	else:
		dSeq=initSeq[1].strip()
	print dSeq,"\n\n"

	print len(dSeq),"\n\n"


#	extendedDel=FiveP[0]+dSeq+ThreeP[0]

	print extendedDel
	sys.exit()

	### Find out if any of these deletions are overlapping with genomic features



	gff=AnnotationRange(file)
	exonA,featureA,symbol=gff_process(gff)

	if len(exonA):
		ordE=sorted(exonA.values(),key=itemgetter(3))#,reverse=True)
		for E in ordE:
			Mexons.append(E)
	### Create order of exons that appear in sequence

	deletion={}
	targff=TargetedRange(file)
	print "\n\n",5*"*","\t",symbol,"\t",fileinfo,"\t",5*"*"
	for t,gffT in targff.iteritems():
		Etar,Ftar,fum=gff_process(gffT)
		if len(Etar):
			delStart,delEnd=map(lambda d: int(d),t.split(":")[1].split("-"))
			for k,j in Etar.iteritems():
				ExonStart,ExonEnd=j[2],j[3]
				if delStart < ExonStart and delEnd > ExonEnd:
					deletion[k]="Full Deletion of exon %s, size %s NTs" % (k,j[3]-j[2])
				elif delStart < ExonStart and delEnd < ExonEnd:
					deletion[k]="5' deletion of %s NTs" % (delEnd-j[2])
				elif delStart > ExonStart and delEnd > ExonEnd:
					deletion[k]="3' deletion of %s NTs" % (ExonEnd-delStart)
				elif delStart > ExonStart and delEnd < ExonEnd:
					deletion[k]="IntraExon deletion of size %s,beginning %s NTs 5' and ending %s NTs 3'" % (delEnd-delStart,delStart-ExonStart,ExonEnd-delEnd)
				else:
					print k,"HUUUH",delStart,j[2],"\t\t",delEnd,j[3]
					sys.exit()
					
	server="http://rest.ensembl.org/sequence/region/mouse/"


	ColumnEx=zip(*Mexons)

	try:
		modStart=ColumnEx[2][0]-1000
		modEnd=ColumnEx[3][-1]+1000
		strand=ColumnEx[6][0]
		chromosome=ColumnEx[1][0]
	except IndexError:
		print Mexons
		sys.exit()
		
	ext=str(chromosome)+":"+str(modStart)+".."+str(modEnd)+":"+str(strand)+"?"+"content-type=text/plain;mask=soft"
	r = requests.get(server+ext)
	#r = requests.get("http://rest.ensembl.org/sequence/region/mouse/10:49098842..49784259:-1?content-type=text/x-fasta;mask=soft")
	 
	if not r.ok:
	  r.raise_for_status()
	  sys.exit()

	seqFlank=r.text.strip()
	if seqFlank.startswith("chrmsmGRCm"):
		seqFlank=seqFlank.lstrip("chrmsmGRCm")


	gDNAlen=abs(modStart-modEnd)



	jointname="_".join(ColumnEx[5])
	locl=APEobj(jointname)

	filetowrite=[]

	filetowrite.append(locl.LOCUS+symbol+"_gDNA\t\t"+str(gDNAlen)+" ds-DNA\t\tlinear\t\t"+time.strftime("%c"))
	filetowrite.append(locl.ACCESSION+symbol)
	filetowrite.append(locl.KEYWORDS +", ".join([symbol,fileinfo]))
	filetowrite.append(locl.SOURCE)
	filetowrite.append(locl.ORGANISM)
	filetowrite.append(locl.COMMENT_BLANK)
	filetowrite.append(locl.COMMENT_BLANK+ ">dna:chromosome:"+"GRCm38:" +chromosome+"("+str(strand)+")"+":"+str(modStart)+"..."+str(modEnd))
	filetowrite.append(locl.COMMENT_BLANK)
	filetowrite.append(locl.COMMENT_BLANK +symbol +" partial gDNA sequence based on "+fileinfo)
	filetowrite.append(locl.COMMENT_BLANK)
	for d,i in targff.iteritems():
		print d
		delInfo=d.split(":")[1].split("-")
		if len(delInfo) == 2:
			delSize=str(abs(int(delInfo[1])-int(delInfo[0])))
		else:
			delSize='1'
		filetowrite.append(locl.COMMENT_BLANK + "Deletion identified of size "+delSize+" : "+d)
	filetowrite.append(locl.COMMENT_BLANK)
	for d,i in Etar.iteritems():
		filetowrite.append(locl.COMMENT_BLANK + "EXON "+d+" implicated: Chr"+Etar[d][1]+":"+str(Etar[d][2])+"-"+str(Etar[d][3])+"...")
		filetowrite.append(locl.COMMENT_BLANK + "\t\tDel info: "+deletion[d])
		filetowrite.append(locl.COMMENT_BLANK)
	filetowrite.append(locl.COMMENT_BLANK)
	filetowrite.append(locl.COMMENT_BLANK +"Created by script on "+time.strftime("%c"))
	filetowrite.append(locl.COMMENT_BLANK)
	filetowrite.append(locl.FEATURES)

	wSeq=seqFlank

	mastLen=len(wSeq)

	masked_sequence=wSeq

	if len(masked_sequence) != len(wSeq):
		print "Error with masked sequence"
	locl.formatSeq(masked_sequence)



	compiled=[]
	for item in Mexons:
#	sys.exit()
		EnsId,chromosEx,e_start,e_end,startPhase,endPhase,e_strand,dummys=item
		if str(startPhase) == "-1":
			startPhase="-"
		if str(endPhase) == "-1":
			endPhase="-"
		EnsId=EnsId+" "+str(startPhase)+"/"+str(endPhase)
		lclinit = _Exon(EnsId)


		exonlength=abs(e_start-e_end)
		endDiff=abs(modEnd-e_end)

		
		a1=int(e_start)-modStart
		s1=mastLen
		s2=abs(e_end-e_start)+1
		if str(e_strand) == "-1":
			#print "ALSO NEG ONE"
		#	filetowrite.append(lclinit.exonl+str(endDiff+1)+".."+str(endDiff+exonlength+1))
			filetowrite.append(lclinit.exonl+str(s1-a1-s2+1)+".."+str(s1-a1))
			
		else:
			filetowrite.append(lclinit.exonl+str(a1+1)+".."+str(a1+s2))
			
		filetowrite.append(lclinit.label)
		filetowrite.append(lclinit.fwd)
		filetowrite.append(lclinit.rvs)
		filetowrite.append(lclinit.arw)


	for t,delz in targff.iteritems():
		lclO=_Other("Deletion:"+t)
		#print t,d
		delInfo=t.split(":")[1].split("-")
		if len(delInfo) == 2:
			delStart,delEnd=int(delInfo[0]),int(delInfo[1])
		else:
			delStart,delEnd=int(delInfo[0]),int(delInfo[0])
		a1=int(delStart)-modStart
		#s1=mastLen
		s2=abs(delEnd-delStart)+1
		filetowrite.append(lclO.exonl+str(a1+1)+".."+str(a1+s2))
		filetowrite.append(lclO.label)
		filetowrite.append(lclO.fwd)
		filetowrite.append(lclO.rvs)
		filetowrite.append(lclO.arw)

	filetowrite.append(locl.SEQUENCE)
	filetowrite.append(locl.formattedseq)

	pl=open(os.path.join("./APEdata",fileinfo+".ape"),'w')
	pl.write("\n".join(filetowrite))
	pl.close()
	comp+=1
	if comp > 5:
		sys.exit()
