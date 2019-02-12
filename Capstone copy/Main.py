
from Bio.Seq import Seq #Seq allows for the creation of a sequence object to store the DNA string
from Bio.Blast import NCBIWWW #Allows access to NCBI wesbite
from Bio.Blast import NCBIXML #Allows access to BLAST database
import sys
from textwrap import wrap

def dnaSequence(userSequence):
	mySeq = Seq(userSequence)

	return mySeq

def splitSequence(myDNA, dnaSplitter):
	newDnaSeq = wrap(myDNA, dnaSplitter)
	return newDnaSeq


def checkNCBI(mySeq):
	result_handle = NCBIWWW.qblast("blastn", "nt", mySeq)  # Searches the NSCBI BLAST database with the nucleotide sequence # in my_seq
	with open("blast.xml", "w") as out_handle:
		out_handle.write(result_handle.read())
	result_handle.close()

def printResults():
	resultHandle = open("blast.xml")
	blast_records = NCBIXML.parse(resultHandle) #The returned information from result_handle is a messy XML file. This returns a cleaner version that is easier to read.
	E_VALUE_THRESH= 0.04

	for blast_record in blast_records:
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				if hsp.expect < E_VALUE_THRESH:
					# with open("test.txt", "w") as f:
					# 	f.write("**** Alignment ****\n")
					# 	f.write("Sequence: {}\n".format(alignment.title))
					# 	f.write("Length: {}\n".format(alignment.length))
					# 	f.write("e value: {}\n".format(hsp.expect))
					# 	f.write("********************\n ")
					# 	f.close()
					print("**** Alignment ****")
					print("Sequence: {}".format(alignment.title))
					print("Length: {}".format(alignment.length))
					print("e value: {}\n".format(hsp.expect))
					print(hsp.quary[0:75] + "...")
					print(hsp.match[0:75] + "...")
					print(hsp.sbjct[0:75] + "...")

	resultHandle.close()

def main():
	mySeq = raw_input("Enter Seqeunce: ")
	dnaSplitter = int(raw_input("Enter divider: "))
	dna = splitSequence(mySeq, dnaSplitter)

	for d in dna:
		print(d)
		mySeq = dnaSequence(d)
		#print("Made it to seqeunce")
		checkNCBI(mySeq)
		#print("made it to check")
		printResults()
		print("%%%%%%%%%%%%%%%%%%Done%%%%%%%%%%%%%%%%%%\n")

	print("\nCompleted")


main()