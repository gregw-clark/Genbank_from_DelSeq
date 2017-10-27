#!/usr/bin/env python

import string, re, sys
import operator
import os
import twobitreader as tbr
from Bio.Seq import Seq
import time
from collections import defaultdict
from itertools import izip
import glob
from operator import itemgetter
from GenBank import GENBANK,_Exon,_Feature,_Other,_OriginalSeq,_OriginalUpstream
import requests

class BlatInfo():
	"""
			This class carries data on Cas9 deletions on a single design basis. Some functions are 
			shared.

			Variables include:
				1.  ChromosomesDir (Static)= Hard-coded path to directory of chromosomes for running cas-offinder
				2.  blat (Static) - Hard-coded path to blat executable
				3.  Chromosome - chromosome for of gene/deletion sequence
				4,  tStart - Start of deletion from cas-offinder results (GRCm38) 
				5.  tEnd - End of deletion from cas-offinder results (GRCm38)
				6.  qStart - Start of deletion in query sequence
				7.  qEnd - End of deletion in query sequence
				8.  name - Name of deltetion. This will be filename suffix for *.class files
				9.  deletions - Nested list for GRCm38 co-ordinates of bp missing in genomic sequence...
					Each entry contains start and end co-ordinates
				10. deletionSequences - Dictionary of co-ordinate info/sequence pairs based on GRCm38... 
					(e.g. '>chr3(-):127679195-127681501':'ATCTCACTTCTCACAGGCAGGAAAT...')
	"""


	ChromosomesDir="/home/clarkg/cas-offinder-master/Chromosomes/"
	blat="/home/clarkg/SilicoPCR/blat"
	def __init__(self,name):
		self.name=name
		
	def group_consecutives(self,vals, step=1):
		run = []
		result = [run]
		expect = None
		for v in vals:
			if (v == expect) or (expect is None):
				run.append(v)
			else:
				run = [v]
				result.append(run)
			expect = v + step
		return result

	def missingRange(self):

		nts=[]
		#print "Start = %s, End =%s" % (self.tStart,self.tEnd)
		for l in self.genomicFlats:
			l.sort()
			nts+=range(l[0],l[1])
		nts=set(nts)
		
		fullrange=set(range(int(self.tStart),int(self.tEnd)))
		#print "Full Length = %s" % (len(fullrange))
		missingNT=list(fullrange.difference(nts))
		missingNT.sort()	
		spans=self.group_consecutives(missingNT)
		self.deletions=[]
		for g in spans:
			if len(g) > 1:
				self.deletions.append([g[0],g[-1]])
			elif len(g) == 1:
				self.deletions.append([g[0]])
				

	def pslReader(self,filename):
	
		io=open(filename,'r').readlines()
		if len(io):
			results=io#[5:]	##The initial 5 lines is some clunky header
		else:
			self.QCPass=False
			results=[]
		blats=[]
		for line in results:
			linedata=[]

			data=line.strip().split("\t")
			#local-hidden from class
			querySize=int(data[10])
			qStart,qEnd=int(data[11]),int(data[12])
			percentCoverage=1-(querySize - abs(qEnd-qStart))/float(querySize)
			querySize=int(data[10])
			qStart,self.qEnd=int(data[11]),int(data[12])
			Chromosome=data[13].strip()
			Insertion=int(data[5])
			qBaseInserts=int(data[5])
			qNumberInserts=int(data[4])
			tNumberInserts=int(data[6])
			tBaseInserts=int(data[7])
			### use locally held data to weed out nonsense results

			if percentCoverage < 0.50 :
				pass
			elif int(data[1])/float(int(data[0])) > 0.05:
				self.PercentageMismatch=int(data[1])/float(int(data[0]))
			else:


				qStarts=data[19]	##Cut Sites in Query
				qBlocks=data[18]	##SSize of Cuts
				tStart=data[15]
				tEnd=data[16]

				qStrt=map(lambda k: int(k),filter(lambda p: len(p),qStarts.split(",")))
				qBlk=map(lambda k: int(k),filter(lambda p: len(p),qBlocks.split(",")))
				tStart,tEnd=data[15],data[16]
				tCuts=data[20].strip().split(",")
				#linedata+=[querySize]
				if data[8] == "-":
					tStrand="-"
					qBlk.reverse()
					qCuts=[(querySize-k) if k > 0 else 0 for k in qStrt]
					if qCuts[0] == 0:
						qCuts.pop(0)
						qCuts.reverse()
						qCuts.insert(0,0)

					tCuts.reverse()
				else:
					tStrand="+"
					qCuts=qStrt

				linedata+=[qCuts]
				linedata+=[qBlk]

				linedata+=[tStart,tEnd,tStrand]
				linedata+=[filter(lambda l: int(len(l)),tCuts)]
				linedata+=[int(data[0])*percentCoverage,int(data[0]),percentCoverage,querySize,qStart,qEnd,Chromosome,qBaseInserts,tBaseInserts]
				blats.append(linedata)
		##Select the best based on matches and coverage	
		blats = sorted(blats, key = operator.itemgetter(6,7),reverse=True)

		if len(blats):
			##
			self.QCPass=True
			BestBlat=blats[0]	##We assume most plausible is based on best percent coverage/size
			##
			self.qCuts=BestBlat[0]
			self.qBlocks=BestBlat[1]
			self.tStart=BestBlat[2]
			self.tEnd=BestBlat[3]
			self.tStrand=BestBlat[4]
			self.tCuts=BestBlat[5]
			self.NumMisMatches=BestBlat[7]
			self.percentCoverage=BestBlat[8]
			self.querySize=BestBlat[9]
			self.qStart=BestBlat[10]
			self.qEnd=BestBlat[11]
			self.Chromosome=BestBlat[12]
			self.tBaseInserts=BestBlat[13]
			self.qBaseInserts=BestBlat[14]
		else:
			self.QCPass=False

	def _examinePSL(self,querySeq):

		chromfile=os.path.join(self.ChromosomesDir,self.Chromosome+".2bit")
		bitfile = tbr.TwoBitFile(chromfile)

		genomicSequence=bitfile[self.Chromosome][int(self.tStart):int(self.tEnd)]	
		self.genomicFlats=[]
		self.DelGenomicSeqs={}
		for x in range(len(self.qBlocks)):
			try:
				qB=int(self.qCuts[x])
				lenB=int(self.qBlocks[x])
				tS=int(self.tCuts[x])
				genomicCut=genomicSequence[(int(tS)-int(self.tStart)):(int(tS)-int(self.tStart)+lenB)]
				self.genomicFlats.append([int(tS),int(tS)+lenB])
				if self.tStrand == "-":
					Gseq=str(Seq(genomicCut).reverse_complement())
				else:
					Gseq=genomicCut
				Qseq=querySeq[qB:(qB+lenB)]
				location=self.Chromosome+"("+self.tStrand+"):"+str(tS)+"-"+str(int(tS)+lenB)
				self.DelGenomicSeqs[location]={'Del':Qseq,'Genomic':Gseq}
			except ValueError:
				pass
		##Send all information to function that will join overlapping ranges
		## And find the deletions in the genomic sequence
		self.missingRange()

	def processCoords(self):
		positions=[]

		##find the extended sequences flanking the deletion	
		for posts in self.deletions:
			if len(posts) == 2:
				start,end=int(posts[0]),int(posts[1])
				positions+=[start,end]
			else:
				start=int(posts[0])
				positions.append(start)

		self.minP=int(self.tStart)
		self.maxP=int(self.tEnd)
				
		server="http://rest.ensembl.org/sequence/region/mouse/"

		if self.Native ==True:
			OriginalSeq=str(self.Chromosome).lstrip("chr")+":"+str(self.minP-1500)+".."+str(self.maxP+1500)+"?"+"content-type=text/plain;mask=soft"
		if self.tStrand == "-":
			extFivePrime=str(self.Chromosome).lstrip("chr")+":"+str(self.minP-1500)+".."+str(self.minP)+"?"+"content-type=text/plain;mask=soft"
			extThreePrime=str(self.Chromosome).lstrip("chr")+":"+str(self.maxP+1)+".."+str(self.maxP+1500)+"?"+"content-type=text/plain;mask=soft"
			
		else:
			extFivePrime=str(self.Chromosome).lstrip("chr")+":"+str(self.minP-1500)+".."+str(self.minP)+"?"+"content-type=text/plain;mask=soft"
			extThreePrime=str(self.Chromosome).lstrip("chr")+":"+str(self.maxP+1)+".."+str(self.maxP+1500)+"?"+"content-type=text/plain;mask=soft"

		self.DelStart,self.DelEnd=self.minP-1500,self.maxP+1500
		
		rFive = requests.get(server+extFivePrime)
		time.sleep(1)

		rThree = requests.get(server+extThreePrime)
		time.sleep(1)

		NativeReturn=requests.get(server+OriginalSeq)
		time.sleep(1)

		NativeSeq=NativeReturn.text.strip()

		self.seqFiveP=rFive.text.strip()
		if self.seqFiveP.startswith("chrmsmGRCm"):
			self.seqFiveP=self.seqFiveP.lstrip("chrmsmGRCm")
		self.seqThreeP=rThree.text.strip()
		if self.seqThreeP.startswith("chrmsmGRCm"):
			self.seqThreeP=self.seqThreeP.lstrip("chrmsmGRCm")

		if self.tStrand == "-":
			self.seqFiveP=str(Seq(self.seqFiveP).reverse_complement())
			self.seqThreeP=str(Seq(self.seqThreeP).reverse_complement())
			#self.QuerySequence=str(Seq(self.QuerySequence).reverse_complement())
			self.FlankedDelSeq=self.seqThreeP+self.QuerySequence+self.seqFiveP
			self.FlankedSeq=str(Seq(NativeSeq).reverse_complement())
		else:
			self.FlankedDelSeq=self.seqFiveP+self.QuerySequence+self.seqThreeP
			self.FlankedSeq=NativeSeq

		print self.QuerySequence


	def RankControlFeatures(self,feature_deletion):
		ftrData=feature_deletion.values()
		#[u'ENSMUSE00000958228', u'10', 24892514, 24892585, 2, u'-1', u'-1', u'+', [u'ENSMUST00000175786']
		deleteList=[]
		for x in range(len(ftrData)):
				exonA,chrA,startA,endA,rankA=ftrData[x][:5]
				numTA=len(ftrData[x][8])
				for y in range(x+1,len(ftrData)):
					exonB,chrB,startB,endB,rankB=ftrData[y][:5]
					numTB=len(ftrData[y][8])
					overlap=abs(startA-startB)+abs(endA-endB)
					if overlap < 30:
						if rankA == rankB:
							if numTA >= numTB:
								deleteList.append(exonB)
							else:
								deleteList.append(exonB)
						elif rankA < rankB:
							deleteList.append(exonB)
						else:
							deleteList.append(exonA)
					elif startA >= startB and endA <= endB:
						if rankA < rankB:
							deleteList.append(exonB)
						else:
							deleteList.append(exonA)
					elif startB >= startA and endB <= endA:
						if rankB < rankA:
							deleteList.append(exonA)
						else:
							deleteList.append(exonB)
		deleteList=list(set(deleteList))
		return deleteList

	def RankControlExons(self,exons):
		ftrData=exons.values()
		#[u'ENSMUSE00000958228', u'10', 24892514, 24892585, 2, u'-1', u'-1', u'+', [u'ENSMUST00000175786']
		deleteList=[]
		for x in range(len(ftrData)):
				print ftrData[x]
				exonA,chrA,startA,endA,rankA=ftrData[x][:5]
				numTA=len(ftrData[x][8])
				for y in range(x+1,len(ftrData)):
					exonB,chrB,startB,endB,rankB=ftrData[y][:5]
					numTB=len(ftrData[y][8])
					overlap=abs(startA-startB)+abs(endA-endB)
					if overlap < 30:
						if rankA == rankB:
							if numTA >= numTB:
								deleteList.append(exonB)
							else:
								deleteList.append(exonB)
						elif rankA < rankB:
							deleteList.append(exonB)
						else:
							deleteList.append(exonA)
					elif startA >= startB and endA <= endB:
						if rankA < rankB:
							deleteList.append(exonB)
						else:
							deleteList.append(exonA)
					elif startB >= startA and endB <= endA:
						if rankB < rankA:
							deleteList.append(exonA)
						else:
							deleteList.append(exonB)
		deleteList=list(set(deleteList))
		return deleteList


	def AnnotationRange(self):
		#print self.deletions
		location=(self.Chromosome).lstrip("chr")+":"+str(self.DelStart)+"-"+str(self.DelEnd)
		server= "https://rest.ensembl.org"
		ext = "/overlap/region/mouse/"+location+"?feature=transcript;feature=exon"

		try:
			r = requests.get(server+ext, headers={ "Content-Type" : "text/x-gff3"})
			time.sleep(1)
			if not r.ok:
				r.raise_for_status()
				sys.exit()

			self.gffdata=filter(lambda l: not l.startswith("#") and len(l),r.text.split("\n"))
		except requests.exceptions.HTTPError:
			print "WEB ERROR",
			print ext
			print r.text()
			self.gffdata=[]

		RExon=re.compile("exon_id=ENSMUSE[0-9]{3,15};")
		exons={}
		exonLine=filter(lambda p: RExon.search(p),self.gffdata)
		altFeatures=filter(lambda p: not RExon.search(p),self.gffdata)
		features={}
		for line in altFeatures:
			data=line.strip().split("\t")
			chromosome,anotBld,anotType,start,end=data[:5]
			strand=data[6]
			infoSp=data[8].strip().split(";")
			for sp in infoSp:
				if sp.startswith("Parent=gene"):
					genename=sp.split("gene:")[1]
					symbol=genename
					self.symbol=symbol
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
			if int(rank) <= 10:
				if exon in exons:
					exons[exon][-1]+=[transcript]
				else:
					exons[exon]=[exon,chromosome,int(start),int(end),int(rank),phase,end_phase,strand,[transcript]]
		exonOrder=[]


		ordE=sorted(exons.values(),key=itemgetter(3))
		for E in ordE:	exonOrder.append(E[0])	##Lets order the exons for GENBANK creation

		#We store feature deletion locally for now
		#But once overlapping features are pruned for TSL significance and overlap, we define as part of self
		overlappingExons={}
		feature_deletion={}

		for deletedRange in self.deletions:
			if len(deletedRange) == 1: 
				delStart,delEnd=deletedRange,deletedRange
			else:
				delStart,delEnd=deletedRange
			for ensemblExon,ExonData in exons.iteritems():

				ExonStart,ExonEnd=ExonData[2],ExonData[3]
				print "General Start %s \t End %s" % (self.tStart,self.tEnd)
				print "Del Start = %s, Del End = %s\t%s" % (delStart,delEnd,deletedRange)
				print "EXON: %s\t\tExon Start %s\tExon End = %s\n\n" % (ExonData[0],ExonStart,ExonEnd)
				if delStart < ExonStart and delEnd > ExonEnd:
					feature_deletion[ensemblExon]="Full Deletion of exon %s, size %s NTs" % (ensemblExon,ExonEnd-ExonStart)
					overlappingExons[ensemblExon]=ExonData
				elif delStart < ExonStart and delEnd < ExonEnd and delEnd > ExonStart:
					feature_deletion[ensemblExon]="5' deletion of %s NTs" % (delEnd-ExonStart)
					overlappingExons[ensemblExon]=ExonData
				elif delStart > ExonStart and delStart < ExonEnd and delEnd > ExonEnd:
					feature_deletion[ensemblExon]="3' deletion of %s NTs" % (ExonEnd-delStart)
					overlappingExons[ensemblExon]=ExonData
				elif delStart > ExonStart and delEnd < ExonEnd:
					feature_deletion[ensemblExon]="IntraExon deletion of size %s,beginning %s NTs 5' and ending %s NTs 3'" % (delEnd-delStart,delStart-ExonStart,ExonEnd-delEnd)
					overlappingExons[ensemblExon]=ExonData
				else:
					##Exon does not overlap with deletion
					#print delStart,"\t\t",delEnd,"\n",ExonStart,"\t\t",ExonEnd
					pass
		#print "Original Exons",len(overlappingExons),overlappingExons.keys()
		deleteFeature=self.RankControlFeatures(overlappingExons)
		deleteExons=self.RankControlExons(exons)
		#print "\n".join(exonOrder)
		for d in list(set(deleteFeature+deleteExons)):
			try:
				del overlappingExons[d]
			except KeyError:
				pass
			try:
				del feature_deletion[d]
			except KeyError:
				pass
			try:
				del exons[d]
			except KeyError:
				pass
			try:
				exonOrder.remove(d)
			except ValueError:
				pass

		if not len(exons):
			print self.name,self.Chromosome,self.tStrand
			print self.deletions
			sys.exit()



		self.overlappingExons=overlappingExons
		self.exons=exons
		self.exonOrder=exonOrder
		self.features=feature_deletion
		self.genename=genename

	def APEdel(self):

		self.lclDelStart=int(self.DelStart)
		self.lclDelEnd=int(self.DelEnd)
		for delS in self.deletions:
		 		if len(delS) == 2:
		 			delSize=abs(delS[1]-delS[0])
		 		else:
		 			delSize=1
		 		if self.tStrand == "-":
		 			self.lclDelStart-=delSize
		 		else:
		 			self.lclDelEnd-=delSize

		gDNAlen=abs(int(self.lclDelEnd)-int(self.lclDelStart))
		locD=GENBANK(self.name+"_del")

		delfiletowrite=[]

		delfiletowrite.append(locD.LOCUS+self.genename+"_gDNA\t\t"+str(gDNAlen)+" ds-DNA\t\tlinear\t\t"+time.strftime("%c"))
		delfiletowrite.append(locD.ACCESSION+self.genename)
		delfiletowrite.append(locD.KEYWORDS +", ".join(["Deletion",self.symbol,self.genename,self.name]))
		delfiletowrite.append(locD.SOURCE)
		delfiletowrite.append(locD.ORGANISM)
		delfiletowrite.append(locD.COMMENT_BLANK+ ">dna:chromosome:"+"GRCm38:" +self.Chromosome+"("+str(self.tStrand)+")"+":"+str(self.DelStart)+"..."+str(self.DelEnd))
		delfiletowrite.append(locD.COMMENT_BLANK)
		delfiletowrite.append(locD.COMMENT_BLANK +self.symbol +" partial gDNA sequence with deletion removed based on sequencing data from "+self.name )
		delfiletowrite.append(locD.COMMENT_BLANK)
		delfiletowrite.append(locD.COMMENT_BLANK +"DELETION INFORMATION:")
		for n,d in enumerate(self.deletions):
			if len(d) == 2:
				delfiletowrite.append(locD.COMMENT_BLANK+"\t"+str(n+1)+". "+self.Chromosome+":"+str(d[0])+"-"+str(d[1])+". Deletion Size: "+str(abs(d[1]-d[0])))
			else:
				delfiletowrite.append(locD.COMMENT_BLANK+"\t"+str(n+1)+". "+self.Chromosome+":"+str(d[0])+". Deletion Size: 1")
		delfiletowrite.append(locD.COMMENT_BLANK)
		delfiletowrite.append(locD.COMMENT_BLANK+ 'TARGETED EXONS:')
		cnt=1
		for exn in self.exonOrder:
			try:
				print self.features[exn],exn
				delfiletowrite.append(locD.COMMENT_BLANK+"\t"+str(cnt)+"."+exn+": "+self.features[exn])
				cnt+=1
			except KeyError:
				pass

		delfiletowrite.append(locD.COMMENT_BLANK)
		delfiletowrite.append(locD.COMMENT_BLANK +"Created by script on "+time.strftime("%c"))
		delfiletowrite.append(locD.COMMENT_BLANK)
		delfiletowrite.append(locD.FEATURES)

	
		masked_sequence=self.FlankedDelSeq
		mastLen=len(masked_sequence)
		locD.formatSeq(masked_sequence)

		delShifts=[]
		side=1
		for cut,lencut in zip(self.tCuts,self.qBlocks):

			#print side%2 ==0,side,side%2
			for delcount,z in enumerate(self.deletions):
				if self.tStrand == "-":
					ds=int(self.DelEnd)-int(self.DelStart)

					dQ=abs(int(cut)-int(self.tEnd))+abs(int(self.tEnd)-int(self.DelEnd))+1
					##Thinking from a negative strand
				
					if len(z) == 2:
						endCUT,startCUT=z
					else:
						endCUT,startCUT=z[0],z[0]

					startCUT=abs(startCUT-int(self.DelEnd))
					endCUT=abs(endCUT-int(self.DelEnd))

					seqSTART=dQ-lencut-1
					seqEND=dQ

					if seqSTART > startCUT:
						seqSTART = startCUT +1
						seqEND = startCUT+lencut


					lclU=_OriginalUpstream(self.Chromosome+":"+cut+"-"+str(int(cut)+lencut))
					delfiletowrite.append(lclU.exonl+str(seqSTART)+".."+str(seqEND))
					delfiletowrite.append(lclU.label)

					#shifter=[seqSTART,seqEND]
					#shifter.sort()
					#delShifts.append(shifter)
					
				else:
					if len(z) == 2:
						startCUT,endCUT=z
					else:
						startCUT,endCUT=z[0],z[0]

					startCUT=abs(startCUT-int(self.DelStart))
					endCUT=abs(endCUT-int(self.DelStart))

					seqSTART=int(cut)-int(self.DelStart)+2
					seqEND=int(cut)-int(self.DelStart)+lencut

					if seqSTART > startCUT:
						seqSTART = startCUT +1
						seqEND = startCUT+lencut
					#print seqSTART+1,seqEND

					lclU=_OriginalUpstream(self.Chromosome+":"+cut+"-"+str(int(cut)+lencut))
					delfiletowrite.append(lclU.exonl+str(seqSTART)+".."+str(seqEND))
					delfiletowrite.append(lclU.label)

					#shifter=[seqSTART,seqEND]
					#shifter.sort()
					#delShifts.append(shifter)

					#shifter=[seqSTART,seqEND]
					#shifter.sort()
					#delShifts.append(shifter)
			side+=1
		#delShifts=sorted(delShifts,key=itemgetter(0))


		compiled=[]
		for item in self.exonOrder:

			EnsId,chromosEx,e_start,e_end,e_rank,startPhase,endPhase,e_strand,associatedTranscripts=self.exons[item]
			if str(startPhase) == "-1":
				startPhase="-"
			if str(endPhase) == "-1":
				endPhase="-"
			EnsId=EnsId+" "+str(startPhase)+"/"+str(endPhase)
			lclEdel = _Exon(EnsId)

			exonlength=abs(e_start-e_end)
		

			a1=int(e_start)-int(self.DelStart)
	
			# s1=mastLen
			# itemlength=abs(e_end-e_start)+1
			# if self.tStrand == "-":
			# 	a2=int(self.DelEnd)-int(e_start)+1
			# 	posA=a2-itemlength+1
			# 	posB=a2
			
			# else:
			# 	posA=a1+1
			# 	posB=a1+itemlength

			for dd in self.deletions:
				dd.sort()
				if self.tStrand == "-":
					dd.reverse()
					if len(dd)==2:
						#dd.sort()
						endCUT,startCUT=dd
						cutLength=abs(startCUT-endCUT)
					else:
						endCUT,startCUT=dd[0],dd[0]
						cutLength=1
				else:
					if len(dd)  == 2:
						startCUT,endCUT=dd
						cutLength=abs(startCUT-endCUT)
					else:
						startCUT,endCUT=dd[0],dd[0]
						cutLength=1

				print "DELETION:",startCUT,endCUT
				print EnsId,e_start,e_end
				intraExon=False
				if e_start >= startCUT and e_end <= endCUT:
					print "Scenario b - Exon within deletion - reference notes in book"
					e_start,e_end=0,0
				elif e_start <= startCUT and e_end <= endCUT and e_end >= endCUT:
					print "Scenario a - Exon spans upstream of deletion - reference notes in book"
					e_end=startCUT
				elif e_start >= startCUT and e_start <= endCUT and e_end >= startCUT:
					print "Scenario d - Exon spans downstream of deletion - reference notes in book"
					e_start,e_end=endCUT-cutLength,e_end - cutLength
				elif e_start <= startCUT and e_end >= endCUT:
					e_startA,e_endA=e_start,startCUT
					e_startB,e_endB=endCUT-cutLength,e_end-cutLength
					intraExon=True
					print "Scenario c - Exon spans past deletion either end- reference notes in book"
				elif e_end > endCUT and self.tStrand == "+":
					print "Scenario e - Exon exists entirely after "
					e_start,e_end=e_start - cutLength-1,e_end-cutLength
				print "\n"



			s1=mastLen
			itemlength=abs(e_end-e_start)+1
			if e_start or e_end:
				if self.tStrand == "-":
					if intraExon:
						itemlengthA=abs(e_endA-e_startA)
						a2=int(self.DelEnd)-int(e_startA)+1
						lclEdel = _Exon("Partial_"+EnsId[:-3])
						delfiletowrite.append(lclEdel.exonl+str(a2-itemlengthA+1)+".."+str(a2))
						delfiletowrite.append(lclEdel.label)


						itemlengthB=abs(e_endB-e_startB)
						a2=int(self.DelEnd)-int(e_startB)+1
						lclEdel = _Exon("Partial_"+EnsId[:-3])
						delfiletowrite.append(lclEdel.exonl+str(a2-itemlengthB+1)+".."+str(a2))
						delfiletowrite.append(lclEdel.label)
		

					else:
						a2=int(self.DelEnd)-int(e_start)+1
						delfiletowrite.append(lclEdel.exonl+str(a2-itemlength+1)+".."+str(a2))
						delfiletowrite.append(lclEdel.label)


				else:
					if intraExon:
						itemlengthA=abs(e_endA-e_startA)
						a1=int(e_startA)-int(self.lclDelStart)
						lclEdel = _Exon("Partial_"+EnsId[:-3])
						delfiletowrite.append(lclEdel.exonl+str(a1+1)+".."+str(a1+itemlengthA))
						delfiletowrite.append(lclEdel.label)


						itemlengthB=abs(e_endB-e_startB)
						a1=int(e_startB)-int(self.lclDelStart)
						lclEdel = _Exon("Partial_"+EnsId[:-3])
						delfiletowrite.append(lclEdel.exonl+str(a1+1)+".."+str(a1+itemlengthB))
						delfiletowrite.append(lclEdel.label)


					else:
						a1=int(e_start)-int(self.lclDelStart)
						delfiletowrite.append(lclEdel.exonl+str(a1+1)+".."+str(a1+itemlength))
						delfiletowrite.append(lclEdel.label)



		delfiletowrite.append(locD.SEQUENCE)
		delfiletowrite.append(locD.formattedseq)

		pl=open(os.path.join("/home/clarkg/IntentShare/APEfiles",self.name+"_del.ape"),'w')
		pl.write("\n".join(delfiletowrite).encode('UTF-8'))
		pl.close()

	def APEnative(self):
		gDNAlen=abs(int(self.DelEnd)-int(self.DelStart))

		locl=GENBANK(self.name)

		filetowrite=[]
		filetowrite.append(locl.LOCUS+self.genename+"_gDNA\t\t"+str(gDNAlen)+" ds-DNA\t\tlinear\t\t"+time.strftime("%c"))
		filetowrite.append(locl.ACCESSION+self.genename)
		filetowrite.append(locl.KEYWORDS +", ".join([self.symbol,self.genename,self.name]))
		filetowrite.append(locl.SOURCE)
		filetowrite.append(locl.ORGANISM)
		filetowrite.append(locl.COMMENT_BLANK+ ">dna:chromosome:"+"GRCm38:" +self.Chromosome+"("+str(self.tStrand)+")"+":"+str(self.DelStart)+"..."+str(self.DelEnd))
		filetowrite.append(locl.COMMENT_BLANK)
		filetowrite.append(locl.COMMENT_BLANK +self.symbol +" partial gDNA sequence based on sequencing data from "+self.name)
		filetowrite.append(locl.COMMENT_BLANK)
		filetowrite.append(locl.COMMENT_BLANK +"DELETION INFORMATION:")
		for n,d in enumerate(self.deletions):
			if len(d) == 2:
				filetowrite.append(locl.COMMENT_BLANK+"\t"+str(n+1)+". "+self.Chromosome+":"+str(d[0])+"-"+str(d[1])+". Deletion Size: "+str(abs(d[1]-d[0])))
			else:
				filetowrite.append(locl.COMMENT_BLANK+"\t"+str(n+1)+". "+self.Chromosome+":"+str(d[0])+". Deletion Size: 1")
		filetowrite.append(locl.COMMENT_BLANK)
		filetowrite.append(locl.COMMENT_BLANK+ 'TARGETED EXONS:')
		cnt=1
		for exn in self.exonOrder:
			try:
				filetowrite.append(locl.COMMENT_BLANK+"\t"+str(cnt)+"."+exn+": "+self.features[exn])
				cnt+=1
			except KeyError:
				pass

		filetowrite.append(locl.COMMENT_BLANK)
		filetowrite.append(locl.COMMENT_BLANK +"Created by script on "+time.strftime("%c"))
		filetowrite.append(locl.COMMENT_BLANK)
		filetowrite.append(locl.FEATURES)

	
		masked_sequence=self.FlankedSeq
		mastLen=len(masked_sequence)
		locl.formatSeq(masked_sequence)



		compiled=[]
		for item in self.exonOrder:
			#[exon,chromosome,int(start),int(end),int(rank),phase,end_phase,strand,[transcript]]
			EnsId,chromosEx,e_start,e_end,e_rank,startPhase,endPhase,e_strand,associatedTranscripts=self.exons[item]
			if str(startPhase) == "-1":
				startPhase="-"
			if str(endPhase) == "-1":
				endPhase="-"
			EnsId=EnsId+" "+str(startPhase)+"/"+str(endPhase)
			lclinit = _Exon(EnsId)


			exonlength=abs(e_start-e_end)

			a1=int(e_start)-int(self.DelStart)
			#print "ID,a1,\t,e_start,e_end\t,self.tStart"
			#print EnsId,a1,"\t",e_start,e_end,"\t",self.tStart
			s1=mastLen
			itemlength=abs(e_end-e_start)+1
			if self.tStrand == "-":
				a2=int(self.DelEnd)-int(e_start)+1
				filetowrite.append(lclinit.exonl+str(a2-itemlength+1)+".."+str(a2))
			else:
				filetowrite.append(lclinit.exonl+str(a1+1)+".."+str(a1+itemlength))
				
			filetowrite.append(lclinit.label)


		if len(self.deletions) > 1:
			pass
			#print "START",self.DelStart
			#print "END",self.DelEnd
			#print "\t".join(self.DelGenomicSeqs.keys())
			#print self.deletions

			#print self.qBlocks


		for cut,lencut in zip(self.tCuts,self.qBlocks):
			if self.tStrand == "-":
				ds=int(self.DelEnd)-int(self.DelStart)

				dQ=abs(int(cut)-int(self.tEnd))+abs(int(self.tEnd)-int(self.DelEnd))+1
				#print abs(int(cut)-int(self.tEnd))+abs(int(self.tEnd)-int(self.DelEnd))
				#print abs(int(cut)-int(self.tEnd))+abs(int(self.tEnd)-int(self.DelEnd)) - lencut
				location=[dQ-lencut,dQ-1]
				location.sort()
			else:

				location=[int(cut)-int(self.DelStart)+lencut,int(cut)-int(self.DelStart)]
				location.sort()

			lclN=_OriginalSeq(self.Chromosome+":"+str(location[0])+"-"+str(location[1]))
		
			filetowrite.append(lclN.exonl+str(location[0]+1)+".."+str(location[1]))

			filetowrite.append(lclN.label)


	
		for delInfo in self.deletions:
			lclO=_Other("Deletion:"+self.Chromosome+":"+"-".join(map(lambda s: str(s),delInfo)))
			if len(delInfo) == 2:
				delStart,delEnd=delInfo[0],delInfo[1]
			else:
				delStart,delEnd=delInfo[0],delInfo[0]
			#a1=int(delStart)-int(self.tStart)
			a1=int(delStart)-int(self.DelStart)
			#s1=mastLensudo mv
			itemlength=abs(delEnd-delStart)+1
			if self.tStrand == "-":
				a2=int(self.DelEnd)-int(delStart)+1
				#filetowrite.append(lclinit.exonl+str(a2-itemlength+1)+".."+str(a2))
				filetowrite.append(lclO.exonl+str(a2-itemlength+1)+".."+str(a2))			
			else:
				#filetowrite.append(lclinit.exonl+str(a1+1)+".."+str(a1+itemlength))
				filetowrite.append(lclO.exonl+str(a1+1)+".."+str(a1+itemlength))
			filetowrite.append(lclO.label)


		filetowrite.append(locl.SEQUENCE)
		filetowrite.append(locl.formattedseq)

		pl=open(os.path.join("/home/clarkg/IntentShare/APEfiles",self.name+".ape"),'w')
		pl.write("\n".join(filetowrite))
		pl.close()

			
