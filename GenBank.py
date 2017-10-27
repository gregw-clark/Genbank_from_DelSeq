#!/usr/bin/env python

import string, re, sys


class GENBANK(object):
	def __init__(self,name):
		self.name = name
		self.formattedseq=""
		self.all_exons=[]
	LOCUS        ="LOCUS       "
	DEFINITION   ="DEFINITION  "
	BLANK        ="            "
	ACCESSION    ="ACCESSION   "
	KEYWORDS     ="KEYWORDS    ."
	SOURCE       ="SOURCE      house mouse"
	ORGANISM     ="  ORGANISM  Mus musculus Eukaryota; Opisthokonta; Metazoa; Eumetazoa;\n"\
		     +"            Bilateria; Deuterostomia; Chordata; Craniata; Vertebrata;\n"\
                     +"            Gnathostomata; Teleostomi; Euteleostomi; Sarcopterygii;\n" \
		     +"            Dipnotetrapodomorpha; Tetrapoda; Amniota; Mammalia; Theria;\n" \
                     +"            Eutheria; Boreoeutheria; Euarchontoglires; Glires; Rodentia;\n" \
                     +"            Sciurognathi; Muroidea; Muridae."	
	COMMENT_BLANK="COMMENT     "
	COMMENT_INIT ="This file was generated using Python. Lauryl Nutter, CMMR via gwc July 2016"
	FEATURES     ="FEATURES             Location/Qualifiers"
	SEQUENCE     ="ORIGIN"	

	def formatSeq(self,sequence):
		seqtmp="".join(map(lambda s: s.strip(), sequence.split()))
		splitseq=[seqtmp[i:(i+10)] for i in range(0,len(seqtmp),10)]
		cnt=1
		inc=1
		lclLine=""
		collected=""
		for k in splitseq:
			spaced=9-len(str(inc))
			if cnt == 1:
				lclLine=" "*spaced+str(inc)+" "+k
				inc+=60
			elif cnt > 1 and cnt < 6:
				lclLine+=" "+k
				
			elif cnt%6 == 0:
				if len(lclLine) and re.search("[0-9]{1,5}",lclLine):
					#if cnt  == 6:
					#	collo
					### If we are at the start
					collected+=lclLine+" "+k+"\n"
				lclLine=" "*spaced+str(inc)
				inc+=60
			else:
				lclLine+=" "+k
			cnt+=1
			
		collected+=lclLine+"\n//"
		
		self.formattedseq=collected
		
class _singleExon():
	
	def __init__(self,name,location):
		print name
		self.name=name
                self.exonl="     exon            "+location#37417..37995"
		#APEobj.label="                     /label="#ENSMUSE00000506880 -/2
		self.label="                     /label="+"\""+name+"\""#ENSMUSE00000506880 -/2
		self.note ="                     /note="+"\""+name+"\""#ENSMUSE00000506880.2\"


	def add_exon(self):
		return "\n".join([self.exonl,self.label,self.note,self.fwd,self.rvs,self.arw])

class _Exon(GENBANK):

	def __init__(self, name):
		GENBANK.__init__(self,name)
		#APEobj.name=name
		#print name, self.name
		self.exonl="     exon            "#37417..37995"
		self.label="                     /label=\""+self.name+"\""		#ENSMUSE00000506880 -/2
		self.note ="                     /note=\""+self.name+"\""	#ENSMUSE00000506880.2\"

	def add_exon(self,name):
		return self.name

class _Other(GENBANK):

	def __init__(self, name):
		GENBANK.__init__(self,name)
		#APEobj.name=name
		#print name, self.name
		self.exonl="     deletion        "#37417..37995"
		self.label="                     /label=\""+self.name+"\""		#ENSMUSE00000506880 -/2
		self.note ="                     /note=\""+self.name+"\""	#ENSMUSE00000506880.2\"

	def add_other(self,name):
		return self.name

class _OriginalSeq(GENBANK):

	def __init__(self, name):
		GENBANK.__init__(self,name)
		#APEobj.name=name
		#print name, self.name
		self.exonl="     input_seq    "#37417..37995"
		self.label="                     /label=\""+self.name+"\""		#ENSMUSE00000506880 -/2
		self.note ="                     /note=\""+self.name+"\""	#ENSMUSE00000506880.2\"

	def add_other(self,name):
		return self.name

class _OriginalUpstream(GENBANK):

	def __init__(self, name):
		GENBANK.__init__(self,name)
		#APEobj.name=name
		#print name, self.name
		self.exonl= '     input_seq       '#37417..37995"
		self.label= "                     /label=\""+self.name+"\""		#ENSMUSE00000506880 -/2
		self.note = "                     /note=\""+self.name+"\""	#ENSMUSE00000506880.2\"

	def add_other(self,name):
		return self.name

class _Feature(GENBANK):

	def __init__(self,name):
		GENBANK.__init__(self,name)
		self.feature="     misc_feature    "#37417..37995"
		self.label="                     /label=\""+self.name+"\""		#ENSMUSE00000506880 -/2
		self.note ="                     /note=\""+self.name+"\""	#ENSMUSE00000506880.2\"

