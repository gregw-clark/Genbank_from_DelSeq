#!/usr/bin/env python

import string, re, sys
from cPickle import load,dump
import multiprocessing
import itertools
import sys
import os
from collections import defaultdict
sys.path.append('/home/clarkg/anaconda2/lib/python2.7/site-packages/')
sys.path.append('/var/www/musmusculims/musculus/C9Req/bin')
sys.path.append('/home/clarkg/_APE_/APEfiles/')
sys.path.append('/home/clarkg/')
from itertools import izip
import MySQLdb as mdb
import glob
import re
from Bio.Seq import Seq
import time
import os
from operator import itemgetter
from APE import APEobj,_Exon,_Feature
from RepeatMasker import repeat_mask
from cogent.db.ensembl import Genome
import requests


def reverseSeq(aseq):
        lclSeq=Seq(aseq)
        revSeq=str(lclSeq.reverse_complement())
        return revSeq

try:
	gene=sys.argv[1]
	genename=gene.capitalize()
	geneID=sys.argv[2]
except IndexError:
	print "No gene entered"
	gene=False
	sys.exit()



print "<html>"
print "<head>"
print "<style>\ntab1 { padding-left: 4em; }\ntab2 { padding-left: 8em; }\n"
print "body{font-family:helvetica} \ntab3 { padding-left: 12em; }\n p{color:#AACCFF}\n"
print "button {padding: 15px 32px;text-align: center; text-decoration: none;display: inline-block;font-size: 16px;}\n li {margin-top: 0px; margin-right: 5px;}\n</style>"
print "</head>"		##AACCFF
print "<body style=\"background-color:#222233;\" text=\"#FFFFA8\">"

mouse=Genome(Species='mouse',Release=87,account=None)
#coding=mouse.getGenesMatching(StableID=geneID)
#print coding
coding=mouse.getGenesMatching(Symbol=genename)

#print mouse.getGenesMatching()
#print dir(mouse.getGenesMatching())
#print mouse.getGenesMatching().__dict__

#print coding
#sys.exit()

if not coding:
	print "<p style=\"font-size: 55px;\">This is a fatal error. Can't find your gene<br>"
	print "</body>"
	sys.exit()
completed=[]
Mexons=[]
for g in coding:
	symbol=genename
	print "<p style=\"font-size: 12px;\">APE file will be generated from all protein-coding transcripts show below.<br>"
	print "Running script ./ENSEMBLget.py at %s.<br>" % time.strftime("%c")
	print "Files are created and stored on Development server,space permitting.</p>"
	print "<p style=\"font-size: 20px;\"><i>"+str(g.Symbol)+"<tab1></tab1>"+str(g.StableId)+"<tab1></tab1>"+str(g.Location)+"</i>"
	print "</p>"
	species,null,chromosome,locations,strand=str(g.Location).strip().split(":")
	giStart,giEnd=map(lambda i: int(i),locations.split("-"))
	if giStart > giEnd:
		seqStart=giEnd
		seqEnd=giStart
	else:
		seqStart=giStart
		seqEnd=giEnd
	#print chrom,geneStart,geneEnd,strand
	#sys.exit()	
	for t in g.Transcripts:
		trans_info=t._table_rows['transcript']#,#__dict__
		if trans_info[9] == "protein_coding":
			print "<b>"+t.StableId+"</b><br>"
			print "<ol style=\"list-style-type:square; font-size: 10px;\">"
			for e in t.Exons:
				#chromin,start,end,strand,sequence,EnsID,startPhase,endPhase,relStart,relEnd=item
				#e,e._table_rows['exon']
				try:

					print "<tab1><li>"+str(e)+"</tab1></li>"
					lclrow=e._table_rows['exon']
					ei_start,ei_end,e_strand,e_startPhase,e_endPhase,e_EnsId=itemgetter(2,3,4,5,6,9)(lclrow)
					if int(ei_start) > int(ei_end):
						e_start=int(ei_end)
						e_end=int(ei_start)
					else:
						e_start=int(ei_start)
						e_end=int(ei_end)
					if e_EnsId not in completed:
						Mexons.append([e_start,e_end,e_strand,e_startPhase,e_endPhase,e_EnsId])		
						completed.append(e_EnsId)
					#else:
					#	completed.append(e_EnsId)
				except KeyError:
					print "\t\tMISSING\t\n",e,e._table_rows,"\n\n"
			print "</ol>"
print "<br><br>"
server="http://rest.ensembl.org/sequence/region/mouse/"

Mexons=sorted(Mexons, key=lambda x: int(x[1]))
ColumnEx=zip(*Mexons)
strand=ColumnEx[2][0]



modStart=int(seqStart)-1000
modEnd=int(seqEnd)+1000
ext=str(chromosome)+":"+str(modStart)+".."+str(modEnd)+":"+str(strand)+"?"+"content-type=text/plain;mask=soft"
r = requests.get(server+ext)
#r = requests.get("http://rest.ensembl.org/sequence/region/mouse/10:49098842..49784259:-1?content-type=text/x-fasta;mask=soft")
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

seqFlank=r.text.strip()
if seqFlank.startswith("chrmsmGRCm"):
	seqFlank=seqFlank.lstrip("chrmsmGRCm")

Mexons=sorted(Mexons, key=lambda x: int(x[1]))

gDNAlen=abs(modStart-modEnd)


jointname="_".join(ColumnEx[5])
locl=APEobj(jointname)



filetowrite=[]

filetowrite.append(locl.LOCUS+symbol.capitalize()+"_gDNA\t\t"+str(gDNAlen)+" ds-DNA\t\tlinear\t\t"+time.strftime("%c"))
filetowrite.append(locl.ACCESSION+str(",".join(ColumnEx[5])))
filetowrite.append(locl.KEYWORDS +symbol)
filetowrite.append(locl.SOURCE)
filetowrite.append(locl.ORGANISM)
filetowrite.append(locl.COMMENT_BLANK+ ">dna:chromosome:"+"GRCm38:" +chromosome+"("+str(strand)+")"+":"+str(modStart)+"..."+str(modEnd))
filetowrite.append(locl.COMMENT_BLANK)
filetowrite.append(locl.COMMENT_BLANK +symbol.capitalize() +" gDNA and 1000bp upstream and downstream flank")
filetowrite.append(locl.COMMENT_BLANK)
filetowrite.append(locl.COMMENT_BLANK +"Created by script on "+time.strftime("%c"))
filetowrite.append(locl.COMMENT_BLANK)
filetowrite.append(locl.FEATURES)

wSeq=seqFlank

mastLen=len(wSeq)

#masked_sequence=repeat_mask(wSeq,genename)
masked_sequence=wSeq

if len(masked_sequence) != len(wSeq):
	print "Error with masked sequence"
locl.formatSeq(masked_sequence)



compiled=[]
for item in Mexons:
	e_start,e_end,e_strand,startPhase,endPhase,EnsId=item		
	if str(startPhase) == "-1":
		startPhase="-"
	if str(endPhase) == "-1":
		endPhase="-"
	EnsId=EnsId+" "+str(startPhase)+"/"+str(endPhase)
	lclinit = _Exon(EnsId)

#	print EnsId,abs(e_start-e_end)
#	print modStart,e_start,"\t",abs(modStart-e_start)
#	print modEnd,e_end,"\t",abs(modEnd-e_end)
#	print "\n\n"


	exonlength=abs(e_start-e_end)
	endDiff=abs(modEnd-e_end)

#	spacerFive=_Feature("Left Frame")
#	spacerThree=_Feature("Right Frame")
	
	
	a1=int(e_start)-modStart
	s1=mastLen
	s2=abs(e_end-e_start)+1
#	filetowrite.append(lclinit.exonl+str(endDiff+1)+".."+str(endDiff+exonlength+1))
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

filetowrite.append(locl.SEQUENCE)
filetowrite.append(locl.formattedseq)


pl=open(os.path.join("/var/www/musmusculims/musculus/C9Req/APEdata/",genename+".ape"),'w')
pl.write("\n".join(filetowrite))
pl.close()
print "<p><b>"
filename=genename+".ape"
print "File written to: %s" % filename
print "</b></p>"
print "<p>Button will trigger an automated download of file.</p>";  
print "<button type=\"button\" ;color:\"#900\"; onclick=\"location.href='Servefile.php?filename="+filename+"'\">Download APE file:<br><br><b>"+filename+"</b></button>";
print "</body>"
print "</html>"
