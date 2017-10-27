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

#dirset=dirset()

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
	time.sleep(1)
	if not r.ok:
		r.raise_for_status()
		sys.exit()

	decoded = r.json()

	for g in range(len(decoded)):
		localID=decoded[g]
		ensID=localID['id']
		ext = "/lookup/id/"+ensID+"?"

		r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
		time.sleep(1)

		if not r.ok:
			r.raise_for_status()
			return False
			sys.exit()
		
		decod = r.json()
		try:
			if len(decoded) > 1:
				print "\n\nEnsembl finds ambiguity from Gene Symbol %s..." % (symbol)
				print decoded,"\n"
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
	if upperSymbol =="SQRDL":		##Officially Symbol == Sqor
		return ["2","+" ,122765359,122809553]		#Chromosome 2: 122,765,237-122,809,569 
	if upperSymbol =="RMST":		##Officially Symbol == RMST_1
		return ["10","-" ,92071037,92165170]		#Chromosome 10: 92,071,037-92,165,170 
	if upperSymbol =="FAM178A":		##Officially Symbol == Slf2
		return ["19","+" ,44931119,44983787]		#Chromosome 19: 44,931,119-44,983,787 forward strand
	return [False,False,False,False]

def _regionOverlap(chrom,start,end):
	location=chrom.lstrip("chr")+":"+str(start)+"-"+str(end)
	server= "https://rest.ensembl.org"
	ext = "/overlap/region/mouse/"+location+"?feature=transcript;feature=exon"

	r = requests.get(server+ext, headers={ "Content-Type" : "text/x-gff3"})
	time.sleep(1)
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

def largestGap(cutlist):
	cutlist=list(set(cutlist))

	cutlist.sort()
	diffs=[(y-x,x) for x,y in zip(cutlist,cutlist[1:])]
	gap=sorted(diffs,key=operator.itemgetter(0),reverse=True)
	midgap=gap[0][1]
	gapsize=gap[0][0]
	indx=cutlist.index(midgap)
	mid=len(cutlist)/2
	if abs(indx-mid) > 2:
		print "WHAAAAT"
		print cutlist,gap[0]
		sys.exit()
	return midgap,gapsize
	


def _readFile():
	hargex=re.compile("_[0-9]{5}_[A-z0-9]{3,15}")
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
						tmpID=re.sub(r'\W+', '', line)+"_"+line.split("_")[-2].strip()
						#print tmpID
						DelID=hargex.findall(tmpID)
						#print DelID[0].split("_")
						dummy,internal,gene=DelID[0].split("_")
						delSeqs[tmpID]=re.sub(r'\W+', '', line)
					elif len(line.split()) == 2 and not dna.match(line.split()[1].strip()):
						tmpID=re.sub(r'\W+', '', line)+"_"+line.split("_")[-2].strip()
						#print tmpID
						DelID=hargex.findall(tmpID)
						#print DelID[0].split("_")
						dummy,internal,gene=DelID[0].split("_")
						delSeqs[tmpID]=""
					elif len(line.split()) > 2:
						print "Please check on sequence %s, something is wrong" % line.split()[0]
						print "\n\nOutput of line below\n\n"
						print line.split()
					else:
						tmpID=re.sub(r'\W+', '', line)+"_"+line.split("_")[-2].strip()
						#print tmpID
						DelID=hargex.findall(tmpID)
						#print DelID[0].split("_")
						dummy,internal,gene=DelID[0].split("_")
						delSeqs[tmpID]=""
				elif  len(line) != 0:
					delSeqs[tmpID]+=re.sub(r'\W+', '', line)
		elif re.search("UCD_GLT",dirset.filename):
			for line in io:
				if line.startswith(">"):
					DelID=re.sub(r'\W+', '', line)+"_"+line.split("_")[-1].strip()
					delSeqs[DelID]=""
				else:
					delSeqs[DelID]+=line.strip()	
		elif re.search("WTSI_2017",dirset.filename):
			header=io[0].strip().split("|")
			accessBYname=dict(zip(header,range(len(header))))
			accessBYindex=dict(zip(range(len(header)),header))
			#print header
			#['mi_attempt_link', 'consortium', 'production_centre', 'marker_symbol', 'mi_date', 'g0_obtained_date', 'genotype_confirmed_date', 'genes_inj', 'strain', 'primary_allele', 'secondary_allele', 'delivery_method', 'electroporation_voltage', 'electroporation_number_of_pulses', 'embryo_transfer_day', 'embryos_survived_to_2_cell', 'e_injected', 'embryos_survived', 'e_transferred', 'pups_born', 'cas9_d10a', 'mrna_protein', 'grnas_per_gene', 'templates_per_gene', 'template_prep', 'reagents', 'grna_sequences', 'grna_coordinates', 'grna_cleave_sites', 'mrna_nuclease_conc', 'protein_nuclease_conc', 'grna_conc', 'template_conc', 'reagents_conc', 'go_screened', 'total_count_of_mutagenized_g0_detected_by_screen', 'g0_with_nhej_mutation', 'g0_with_deletion_mutation', 'g0_with_hr_mutation', 'g0_with_hdr_mutation', 'g0_all_donors_inserted', 'g0_subset_of_donors_inserted', 'g0_bred', 'g0_glt', 'genotype_confirmed_f1s', 'experimental', 'repeated', 'comments', 'comments_abount_alleles', 'status', '1st Del Chr', '1st start bp', '1st end bp', '2nd Del Chr', '2nd start bp', '2nd end bp', '3rd Del Chr', '3rd start bp', '3rd end bp', 'DELETION ALLELE FOR IMITS']
			#sys.exit()
			io.pop(0)
			for line in io:
				data=line.strip().split("|")
				seq=data[accessBYname['DELETION ALLELE FOR IMITS']]
				DelID="_".join(	[
						data[accessBYname['mi_attempt_link']].split("/")[-1],
						data[accessBYname['production_centre']],
						data[accessBYname['marker_symbol']]
						])	

				grna_coords=data[accessBYname['grna_coordinates']]
				grna_cut_sites=map(lambda i: int(i),data[accessBYname['grna_cleave_sites']].split(","))
				grna_cut_sites.sort()
			
				if data[accessBYname['status']] == "Genotype confirmed":
				#	print ">"+DelID+"----"+data[accessBYname['status']]
				#	print seq
					midCut,smallGap=largestGap(grna_cut_sites)
					largeGap=max(grna_cut_sites)-min(grna_cut_sites)
				#	print grna_coords,"\n",smallGap,largeGap,"\n\n"
					if len(seq) > 0:
						delSeqs[DelID]=seq
					else:
						print "Missing INFO"
						print ">"+DelID+"----"+data[accessBYname['status']]
						print data,"\n\n"
						#print seq
						#print DelID
		dirset.delSeqs=delSeqs
	except IndexError:
		print "No file has been selected"
		print "Please type \'#python findDels.py <filename>\'"
		sys.exit()



if __name__ == "__main__":

	openLog()
	_readFile()

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
	for fileinfo,deletion_sequence in dirset.delSeqs.iteritems():
		gene_name=fileinfo.split("_")[-1]
	#	print a,gene_name#,dirset.geneinfo[gene_name][0]
		CHRFILE=os.path.join(dirset.chromosomes,"AllChromosomes.txt")
		try:
			found_chr=dirset.geneinfo[gene_name][0]
			chrpath=os.path.join(dirset.chromosomes,"chr"+found_chr+".2bit")
		except KeyError:
				##Can't find an associated gene name/chromosomes
				# so we run against all chromosomes
				chrpath=CHRFILE

		lclfile=fileinfo+".txt"
		outfile=fileinfo+".psl"
		if not os.path.exists(outfile):
				
			io=open(lclfile,'w')
			io.write(">"+fileinfo+"_Chr:"+found_chr+"\n"+deletion_sequence)
			io.close()

			process = subprocess.Popen([dirset.blat, chrpath,"-noHead",lclfile,outfile],stdout=subprocess.PIPE)
			stdoutdata, stderrdata = process.communicate()


		if re.search(fileinfo,"Ptchd2"):
			print results
			sys.exit()


		if not os.path.exists(fileinfo+".class"):	
			binfo=BlatInfo(fileinfo)
			binfo.QuerySequence=deletion_sequence
			binfo.pslReader(outfile)

		
			if not binfo.QCPass:
				##This means the expected chromosome was not found
				##run a second time against all chromosomes
				process = subprocess.Popen([dirset.blat, CHRFILE,"-noHead",lclfile,outfile],stdout=subprocess.PIPE)
				stdoutdata, stderrdata = process.communicate()
				binfo.pslReader(outfile)
			if not binfo.QCPass:
				print "CAUSING PROBLEMS",outfile
				print "Percent Coverage = %s" % (binfo.PercentCoverage)
				print "Percentage of Mismatches = %s" % (binfo.PercentageMismatch)

			else:
				binfo._examinePSL(deletion_sequence)
				
				chromfile=os.path.join(binfo.ChromosomesDir,binfo.Chromosome+".2bit")
				bitfile = tbr.TwoBitFile(chromfile)
				delfile=fileinfo+".class"
				io=open(delfile,'wb')
				print "File %s has %s deletions..." % (fileinfo,len(binfo.deletions))
				print "%s Strand" % (binfo.tStrand)
				binfo.deletionSequences={}
				for delRange in binfo.deletions:
					if len(delRange) == 2:
						startD,endD=delRange
						delLocation=">"+binfo.Chromosome+"("+binfo.tStrand+"):"+str(startD)+"-"+str(endD)
						delSequence=bitfile[binfo.Chromosome][startD:(endD+1)]
						binfo.deletionSequences[delLocation]=delSequence	
					elif len(delRange) == 1:
						delLocation=">"+binfo.Chromosome+"("+binfo.tStrand+"):"+str(delRange[0])
						delSequence=bitfile[binfo.Chromosome][delRange[0]]
						binfo.deletionSequences[delLocation]=delSequence	
					else:
						print "WWWWWWWWWWWWWWWWWWWWWWTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
						sys.exit()
				dump(binfo,io)
				io.close()
				if binfo.querySize != binfo.qEnd-binfo.qStart:
					print binfo.name,binfo.percentCoverage
					print binfo.querySize,"\t",binfo.qEnd-binfo.qStart

					print "Query Start:%s \t Query End: %s" % (binfo.qStart,binfo.qEnd)
					sys.exit()
				time.sleep(2)
				binfo.Native=True
				binfo.processCoords()
				binfo.AnnotationRange()
				binfo.APEnative()
				binfo.APEdel()	
				sys.exit()	

