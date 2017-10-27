#!/usr/bin/env pyhton

import string, re, sys,os
from _dirset import dirset
from DELclass import BlatInfo
import time
from cPickle import dump,load
import subprocess
import twobitreader as tbr
from Bio.Seq import Seq
from Bio import pairwise2
import requests
import time
import operator

dirset=dirset()

def openLog():
	try:
		kl=load(open("./LocalLogs/GeneInfo.pkl",'rb'))
		dirset.geneinfo=kl
	except:
		dirset.geneinfo={}
def closeLog():
	io=open("./LocalLogs/GeneInfo.pkl",'wb')
	dump(dirset.geneinfo,io)
	io.close()

dna=re.compile("[ACT\/Gactg]{20}")



def EnsAPI(symbol):
	import requests, sys
	upperSymbol=symbol.upper()
	server = "https://rest.ensembl.org"
	ext = "/xrefs/symbol/mus_musculus/"+symbol.upper()+"?"

	if upperSymbol == "FAM188B":
		upperSymbol = "MINDY4"
	elif upperSymbol == "FAM188A":
		upperSymbol = "MINDY3"
	elif upperSymbol == "PTCHD2":
		upperSymbol = "DISP3"
	r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
	if not r.ok:
		r.raise_for_status()
		sys.exit()

	decoded = r.json()

	for g in range(len(decoded)):
		localID=decoded[g]
		ensID=localID['id']
		ext = "/lookup/id/"+ensID+"?"

		r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

		if not r.ok:
			r.raise_for_status()
			return False
			sys.exit()
		
		decod = r.json()
		try:
			if len(decoded) > 1:
				print "\n\n"
				print decoded,"\n\n"
				print symbol,decod['display_name']
			ensSymbol=str(decod['display_name'])
			ensUpper=ensSymbol.upper()
			if ensUpper == upperSymbol:#symbol.strip.upper():
				print symbol,ensID,localID,len(decoded)
				chromosome=decod['seq_region_name']
				start=decod['start']
				end=decod['end']
				strand=decod['strand']
				return [chromosome,strand,start,end]
		except KeyError:
			pass
	if upperSymbol == "2700054A10RIK":
		return ["17","-",13487021,13554094]
	return [False,False,False,False]

def _regionOverlap(chrom,start,end):
	location=chrom+":"+str(start)+"-"+str(end)
	server= "https://rest.ensembl.org"
	ext = "/overlap/region/mouse/"+location+"?feature=transcript;feature=exon"

	r = requests.get(server+ext, headers={ "Content-Type" : "text/x-gff3"})

	if not r.ok:
		r.raise_for_status()
		sys.exit()


	print r.text


def _seqIdentity(seqA,seqB):
	#"A" suffix is always our query, seqA,gapA,alnA
	#"B" suffix is the genomic
	errors=0
	length=0
	aligned=pairwise2.align.globalms(seqA,seqB,2, -.5, -1, -.6)
	alnA,alnB=aligned[0][:2]
	gapA,gapB=0,0
#	print ""
#	print ">ALN A","\n",alnA
#	print ">ALN B","\n",alnB
	for a,b in zip(alnA,alnB):
		if not re.search("-",a+b) and (a.upper() != b.upper()):
			errors+=1
		if a == "-":
			gapA+=1
		if b == "-":
			gapB+=1
		length+=1
	if errors > 10:
		print "TOO MANY ERRORS",errors
		sys.exit()
	elif length < 10:
		print "Missing sequence",len(seqA),len(seqB)
		sys.exit()
	return errors,length,gapA,gapB

def select(lst, *indices):
    return (lst[i] for i in indices)

def group_consecutives(vals, step=1):
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



def missingRange(genomicFlats,start,end):
	nts=[]
	#print "Start = %s, End =%s" % (start,end)
	for l in genomicFlats:
		l.sort()
		nts+=range(l[0],l[1])
	nts=set(nts)
	
	fullrange=set(range(int(start),int(end)))
	print "Full Length = %s" % (len(fullrange))
	missingNT=list(fullrange.difference(nts))
	missingNT.sort()	
	print "Number missing %s" % (len(missingNT))
	spans=group_consecutives(missingNT)
	deletions=[]
	print genomicFlats
	print spans,start,end
	#print len(spans)
	#if not len(spans):
	#	return []
	for g in spans:
		if len(g) > 1:
	#		print "MISSING SPAN:",g[0],g[-1],"\t\tLENGTH:",len(g)
			deletions.append([g[0],g[-1]])
		elif len(g) == 1:
			#print "INDEL",g,spans#int(g[0])-int(start)
			#print g,len(g)
			deletions.append([g[0]])
		else:
			print len(g),"Nothing here"
	#print deletions
	return deletions		

def _initAlign(querySeq,genomic,strand):
	if strand == "-":
		genomic=str(Seq(genomic).reverse_complement())

	#print ">FULL QUERY"
	#print querySeq,"\n"

	#print ">FULL GENOMIC"
	#print genomic,"\n"

	#print "Length of Query=%s, Length of Genomic=%s" % (len(querySeq),len(genomic))
	aligned=pairwise2.align.globalms(querySeq,genomic,2, -.5, -1, -.6)
	alnA,alnB=aligned[0][:2]
	if len(genomic) < len(querySeq):
		print "INSERTION OF LENGTH %s" % (len(querySeq)-len(genomic))
		print ">Full Query Sequence"
		print alnA

		print ">Genome Sequence"
		print alnB,"\nLength of   query = %s" % (len(querySeq))
		print        "Length of genomic = %s" % (len(genomic))
		#sys.exit()
		gapGenomic=[]
		lclCount=0
		for a,b in zip(alnA,alnB):
			if b=="-":
				gapGenomic.append(lclCount)
			lclCount+=1
		#gapSpan=group_consecutives(gapGenomic)
		return gapGenomic
	gapA,gapB=0,0
	length,errors=0,0
	for a,b in zip(alnA,alnB):
		if not re.search("-",a+b) and (a.upper() != b.upper()):
			errors+=1
		if a == "-":
			gapA+=1
		if b == "-":
			gapB+=1
		length+=1
	if abs(len(querySeq)-len(genomic)) == gapA:
		status="delQ"+str(gapA)
	elif abs(len(querySeq)-len(genomic)) == gapB:
		status="insQ"+str(gapB)
	elif (len(querySeq)+gapA) == (len(genomic)+gapB):
		status="unk"
	return False




def _examinePSL(querySeq,blatinfo):

	multiDels=[]
	for blatinfo in results:
		
		qsize,qCuts,qblocks,tStart,tEnd,tStrand,tCuts,matches,coverage,chromstr=blatinfo
		chromfile=os.path.join(dirset.chromosomes,chromstr+".2bit")
		bitfile = tbr.TwoBitFile(chromfile)

		genomicSequence=bitfile[chromstr][int(tStart):int(tEnd)]	
		#flagged=_initAlign(querySeq,genomicSequence,tStrand)
		#if flagged:
		#	genomeGaps=map(lambda l :int(tStart)+l,flagged)
		#	gapSpans=group_consecutives(genomeGaps)
		#	#print "INSERTIONS INTO QUERY at points:"
		#	for i,gp in enumerate(gapSpans):
		#		if len(gp) == 2:
		#			print "\t\t"+str(i+1)+".)\t"+chromstr+":"+"-".join([str(gp[0]),str(gp[1])])
		#		elif len(gp) == 1:
		#			print "\t\t"+str(i+1)+".)\t"+chromstr+":"+str(gp[0])
		#	qblocks=[]

		genomicFlats=[]
		for x in range(len(qblocks)):
			try:
				qB=int(qCuts[x])
				lenB=int(qblocks[x])
				tS=int(tCuts[x])
				genomicCut=genomicSequence[(int(tS)-int(tStart)):(int(tS)-int(tStart)+lenB)]
				genomicFlats.append([int(tS),int(tS)+lenB])
				location=chromstr+"("+tStrand+"):"+str(tS)
				if tStrand == "-":
					Gseq=str(Seq(genomicCut).reverse_complement())
				else:
					Gseq=genomicCut
				Qseq=querySeq[qB:(qB+lenB)]
			#	print ">Query\n",Qseq
			#	print ">Genom\n",Gseq
			except ValueError:
				pass

		deletions=missingRange(genomicFlats,tStart,tEnd)
		if len(deletions):
			deletions.insert(0,tStrand)
			multiDels.append(deletions)
#	print "\n\n"
#	print len(multiDels[0])	
#	time.sleep(5)
	if len(multiDels):
		return multiDels[0]
	else:
		return []

def _readFile():
	try:
		sys.argv[1]
		dirset.filename=sys.argv[1]
		io=open(dirset.filename,'r').readlines()
		delSeqs={}

		if re.search("Harwell",dirset.filename):
			for line in io:
				#line=	
				if ord(line[0]) == 62 or ord(line[0]) >122:
					if len(line.split()) == 2 and dna.match(line.split()[1].strip()):
						DelID=re.sub(r'\W+', '', line)+"_"+line.split("_")[-2].strip()
						delSeqs[DelID]=re.sub(r'\W+', '', line)
					elif len(line.split()) == 2 and not dna.match(line.split()[1].strip()):
						DelID=re.sub(r'\W+', '', line)+"_"+line.split("_")[-2].strip()
						delSeqs[DelID]=""
					elif len(line.split()) > 2:
						print "Please check on sequence %s, something is wrong" % line.split()[0]
						print "\n\nOutput of line below\n\n"
						print line.split()
					else:
						DelID=re.sub(r'\W+', '', line)+"_"+line.split("_")[-2].strip()
						delSeqs[DelID]=""
				elif  len(line) != 0:
					delSeqs[DelID]+=re.sub(r'\W+', '', line)
		elif re.search("UCD_GLT",dirset.filename):
			for line in io:
				if line.startswith(">"):
					DelID=re.sub(r'\W+', '', line)+"_"+line.split("_")[-1].strip()
					delSeqs[DelID]=""
				else:
					delSeqs[DelID]+=line.strip()	
		dirset.delSeqs=delSeqs
	except IndexError:
		print "No file has been selected"
		print "Please type \'#python findDels.py <filename>\'"
		sys.exit()
	


if __name__ == "__main__":
	openLog()
	_readFile()

	#dfunc.pslReader("stk.psl")
	#sys.exit()
	for k, j in dirset.delSeqs.iteritems():
		symbol=k.split("_")[-1]
		try:
			dirset.geneinfo[symbol]
			#print symbol,k.split("_"),dirset.geneinfo[symbol]
			#print j
		except KeyError:
			seqinfo=EnsAPI(symbol)
			if seqinfo[0]:
				dirset.geneinfo[symbol]=seqinfo
			else:
				print "PASSED", symbol,seqinfo
				sys.exit()
	closeLog()
	cnter=0
	for a,deletion_sequence in dirset.delSeqs.iteritems():
		gene_name=a.split("_")[-1]
		print a,gene_name#,dirset.geneinfo[gene_name][0]
		CHRFILE=os.path.join(dirset.chromosomes,"AllChromosomes.txt")
		try:
			found_chr=dirset.geneinfo[gene_name][0]
			chrpath=os.path.join(dirset.chromosomes,"chr"+found_chr+".2bit")
		except KeyError:
			##Can't find an associated gene name/chromosomes
			# so we run against all chromosomes
			chrpath=CHRFILE

		lclfile=a+".txt"
		outfile=a+".psl"
		if not os.path.exists(outfile):
				
			io=open(lclfile,'w')
			io.write(">"+a+"_Chr:"+found_chr+"\n"+deletion_sequence)
			io.close()

			process = subprocess.Popen([dirset.blat, chrpath,"-noHead",lclfile,outfile],stdout=subprocess.PIPE)
			stdoutdata, stderrdata = process.communicate()


		if re.search(a,"Ptchd2"):
			print results
			sys.exit()


		if not os.path.exists(a+".del"):	
			#binfo=BlatInfo(a)
			#binfo.pslReader(outfile)
			results=dirset.pslReader(outfile)
			if not len(results):
				##run a second time against all chromosomes
				process = subprocess.Popen([dirset.blat, CHRFILE,"-noHead",lclfile,outfile],stdout=subprocess.PIPE)
				stdoutdata, stderrdata = process.communicate()
				results=dirset.pslReader(outfile)
			if not len(results):
				print "CAUSING PROBLEMS",outfile

			else:
				delReturn=_examinePSL(deletion_sequence,results)
				delStrand=delReturn[0]
				deletions=delReturn[1:]

				chromstr=results[0][-1]
				chromfile=os.path.join(dirset.chromosomes,chromstr+".2bit")
				bitfile = tbr.TwoBitFile(chromfile)
				if len(deletions) > 1:
					#_mergeDeletion(deletions)
					deletions=sorted(deletions, key = operator.itemgetter(0))
					delfile=a+".del"
					io=open(delfile,'w')
					for delspot in deletions:
						#print delspot
						if len(delspot) == 2:
							startD,endD=delspot
							io.write(">"+chromstr+"("+delStrand+"):"+str(startD)+"-"+str(endD)+"\n")
							print chromstr,"("+delStrand+")",startD,endD
							print bitfile[chromstr][startD:(endD+1)]
							io.write(bitfile[chromstr][startD:(endD+1)]+"\n")
						else: 
							#pass
							io.write(">"+chromstr+"("+delStrand+"):"+str(delspot[0])+"\n")
							print chromstr,"("+delStrand+")",delspot[0]
							print bitfile[chromstr][delspot[0]]
							io.write(bitfile[chromstr][delspot[0]]+"\n")
					#io.close()
				else:
					#print "DELETIONS"
					#print deletions
					#delfile=a+".del"
					#io=open(delfile,'w')
					if len(deletions[0]) == 2:
						startD,endD=deletions[0]
						print chromstr,"("+delStrand+")",startD,endD
						print bitfile[chromstr][startD:(endD+1)]
						#io.write(">"+chromstr+"("+delStrand+"):"+str(startD)+"-"+str(endD)+"\n")
						#io.write(bitfile[chromstr][startD:(endD+1)]+"\n")
						#io.close()
					else:
						print deletions,len(deletions[0])
			#print "*****\t\t\t"+str(cnter)+"\t\t******\n\n\n"
			cnter+=1


