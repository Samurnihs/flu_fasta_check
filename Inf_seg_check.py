import argparse
import pandas as pd 
import sys
import os
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def make_kmers(st, length=10, number=500):
    """The function cuts the set of kmers from the given string.

    :param st: string for cutting k-mers
    :type st: str
    :param length: k-mer length (Default value = 10)
    :type length: int
    :param number: number of kmers (Default value = 500)
    :type number: int
    :returns: list of kmers


    """
    kmers = list()
    while len(kmers) < number: # works until number of kmers less thsn 500
        a = random.randint(0, len(st) - length) # random cutting kmer
        mer = st[a:a + length].upper() # bringing to upper register
        if 'N' not in mer: # only k-mer without N can be used
            kmers.append(mer)
    return kmers


def kmer_score(kmers, ref):
    """Calculates the fraction of k-mers found in sequence.

    :param kmers: list of k-mers. The single kmers has str type
    :param ref: the string in which to search for k-mers
    :returns: the fraction of k-mers found in sequence (float)

    """
    refu = ref.upper() # bringing to upper register
    n = 0
    for kmer in kmers:
        if kmer in refu:
            n+=1 # counting a found k-mer
    return n/len(kmers) # part of k-mers found in reference string



def seg_score(seq_record, refs, segnum):
    """Calculates the maximal measure of conformity between sequence from fasta file and some segment from reference sequnces.

    :param seq_record: a single record from fasta file
    :type seq_record: SeqRecord (Biopython)
    :param refs: list of reference records 
    :type refs: list of SeqRecords (Biopython)
    :param segnum: number of segment needed to be checked
    :type segnum: int
    :returns: maximal fraction of found k-mers (float) 

    """
    kmers = make_kmers(str(seq_record.seq)) # obtaining kmers from sequence
    cnds = [0] * len(refs) # empty list for scores
    for i in range(len(refs)):
        cnds[i] = kmer_score(kmers, str(refs[i][segnum-1].seq))
    return max(cnds) # returns maximal score


def isseg(seq_record, refs, segnum, threshold=0.1):
    """Calculates the verdict about confirmity of sequence from fasta file to some segment from reference sequnces.

    :param seq_record: a single record from fasta file
    :type seq_record: SeqRecord (Biopython)
    :param refs: list of reference records 
    :type refs: list of SeqRecords (Biopython)
    :param segnum: number of segment needed to be checked
    :type segnum: int
    :param threshold:  (Default value = 0.1)
    :type threshold: float
    :returns: True or False (bool)

    """
    if seg_score(seq_record, refs, segnum) >= threshold: # if score is higher than threshold, returns True, otherwise False
        return True
    else:
        return False


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='The script contains the set of functions needed to check if a sequence is the certain segment of the Influenza A/B genome. Use it with parameter -help for more detailed info.')
	parser.add_argument('input', type=str, nargs='+', help='Path to fasta file to analyze. In the example all sequences from the file will be checked.')
	parser.add_argument('type', type=str, help='Type of Influenza virus (A or B).', choices= ['A', 'B'])
	parser.add_argument('-seg', type=int, nargs='+', help='Number of the segement(s) to check the sequence for. \nPB2: 1\nPB1 and PB1-F2: 2\nPA: 3\nHA: 4\nNP: 5\nNA: 6\nM2 and M1: 7\nNS2 and NS1: 8', \
		choices=[1, 2, 3, 4, 5, 6, 7, 8], default=[1, 2, 3, 4, 5, 6, 7, 8])
	args = parser.parse_args()


	refpack = {'A': 
		list(map(lambda x: list(SeqIO.parse(x, "fasta")), list(map(lambda x: os.path.join(os.getcwd(),'InfRefs', \
		'Influenza_A', x), os.listdir(os.path.join(os.getcwd(),'InfRefs', 'Influenza_A')))))), \
		'B':
		list(map(lambda x: list(SeqIO.parse(x, "fasta")), list(map(lambda x: os.path.join(os.getcwd(),'InfRefs', \
		'Influenza_B', x), os.listdir(os.path.join(os.getcwd(),'InfRefs', 'Influenza_B'))))))} # dictionary with reference sequences

	df = pd.DataFrame(columns=['Sequence', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8']) # the output dataframe
	for file in args.input:
		seq = list(SeqIO.parse(file, "fasta"))
		oar = [None] * 9
		for i in range(len(seq)):
			oar[0] = seq[i].description
			for s in args.seg:
				oar[s] = isseg(seq[i], refpack[args.type], s)
			df.loc[len(df)] = oar

	df.to_csv(sys.stdout)
