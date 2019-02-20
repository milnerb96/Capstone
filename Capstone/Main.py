#Blaze Milner
#DNA Sequence Project

from Bio.Seq import Seq #Seq allows for the creation of a sequence object to store the DNA string
from Bio.Blast import NCBIWWW #Allows access to NCBI wesbite
from Bio.Blast import NCBIXML #Allows access to BLAST database
import csv
from textwrap import wrap

def dnaSequence(userSequence):
	mySeq = Seq(userSequence)

	return mySeq

def splitSequence(myDNA, dnaSplitter):
	newDnaSeq = wrap(myDNA, dnaSplitter)
	return newDnaSeq


def checkNCBI(mySeq):
	result_handle = NCBIWWW.qblast("blastn", "nt", mySeq)  # Searches the NSCBI BLAST database with the nucleotide sequence # in my_seq

	return result_handle

def printResults(resultHandle):
	blast_records = NCBIXML.parse(resultHandle) #The returned information from result_handle is a messy XML file. This returns a cleaner version that is easier to read.
	E_VALUE_THRESH = 0.0

	for blast_record in blast_records:
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				if hsp.expect is not None:
					print("\n********************** Alignment **********************")
					print("Sequence: {}".format(alignment.title))
					print("Length: {}".format(alignment.length))
					print("e value: {}\n".format(hsp.expect))
					print("Start Value: {}".format(hsp.query_start)),
					print("End Value: {}".format((hsp.query_start + hsp.align_length) - 1))
					#print(hsp.query[0:75] + "...")
					#print(hsp.match[0:75] + "...")
					#print(hsp.sbjct[0:75] + "...")

					with open("DNA_file.csv", mode = "a") as DNA_file:
						DNA_writer = csv.writer(DNA_file, delimiter = ",")
						DNA_writer.writerow([alignment.title,hsp.query_start, (hsp.query_start + hsp.align_length) - 1])
				else:
					print("No matches!")

def main():
	mySeq = raw_input("Enter Seqeunce: ")
	dnaSplitter = int(raw_input("Enter divider: "))
	dna = splitSequence(mySeq, dnaSplitter)

	for d in dna:
		print(d)
		mySeq = dnaSequence(d)
		blastXML = checkNCBI(mySeq)
		printResults(blastXML)
		print("%%%%%%%%%%%%%%%%%%Done%%%%%%%%%%%%%%%%%%\n")

	print("\nCompleted")


main()