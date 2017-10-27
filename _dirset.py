#!/usr/bin/env python

import re
import glob
import operator
###hard-coded directories
class dirset(object):

	chromosomes="/home/clarkg/cas-offinder-master/Chromosomes/"
	blat="/home/clarkg/SilicoPCR/blat"
	filename=""


	



	def pslReader(self,filename):
		##PSL format taken from web
#		0	matches - Number of matching bases that aren't repeats.
#		1	misMatches - Number of bases that don't match.
#		2	repMatches - Number of matching bases that are part of repeats.
#		3	nCount - Number of 'N' bases.
#		4	qNumInsert - Number of inserts in query.
#		5	qBaseInsert - Number of bases inserted into query.
#		6	tNumInsert - Number of inserts in target.
#		7	tBaseInsert - Number of bases inserted into target.
#		8	strand - defined as + (forward) or - (reverse) for query strand. In mouse, a second '+' or '-' indecates genomic strand.
#		9	qName - Query sequence name.
#		10	qSize - Query sequence size.
#		11	qStart - Alignment start position in query.
#		12	qEnd - Alignment end position in query.
#		13	tName - Target sequence name.
#		14	tSize - Target sequence size.
#		15	tStart - Alignment start position in target.
#		16	tEnd - Alignment end position in target.
#		17	blockCount - Number of blocks in the alignment.
#		18	blockSizes - Comma-separated list of sizes of each block.
#		19	qStarts - Comma-separated list of start position of each block in query.
#		20	tStarts - Comma-separated list of start position of each block in target.


		
			io=open(filename,'r').readlines()
			if len(io):
				results=io#[5:]	##The initial 5 lines is some clunky header
			else:
				results=[]
				return results
			blats=[]

			for line in results:
				linedata=[]

				data=line.strip().split("\t")
				querySize=int(data[10])
				qStart,qEnd=int(data[11]),int(data[12])
				percentCoverage=1-(querySize - abs(qEnd-qStart))/float(querySize)
				chromstr=data[13].strip()
				if percentCoverage < 0.50 :
					pass
				elif int(data[1])/float(int(data[0])) > 0.05:
					pass
				else:

					qStarts=data[19]	##Cut Sites in Query
					qBlocks=data[18]	##SSize of Cuts

					tStart=data[15]
					tEnd=data[16]

					qStrt=map(lambda k: int(k),filter(lambda p: len(p),qStarts.split(",")))
					qBlk=map(lambda k: int(k),filter(lambda p: len(p),qBlocks.split(",")))
					tStart,tEnd=data[15],data[16]
					tCuts=data[20].strip().split(",")
					linedata+=[querySize]
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
					linedata+=[int(data[0])*percentCoverage,int(data[0]),chromstr]
					blats.append(linedata)
			blats = sorted(blats, key = operator.itemgetter(7,8),reverse=True)
			##We select the best based on matches and coverage	
			if len(blats):
				return [blats[0]]
			else:
				return []
