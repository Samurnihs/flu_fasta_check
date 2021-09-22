import argparse
import pandas as pd 
import sys
import os
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import IUPACData

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

def check_empty(sequences, input_file):
    """If no fasta sequences were found, check that file is empty"""
    if not sequences:
        if os.stat(input_file).st_size != 0:
            return False

    return True

def check_header(segment, seq_ids):
    """Check fasta id correctness for each record in file"""
    if len(seq_ids) == 1:
        if int(seq_ids[0]) != segment:
            print(f"Invalid fasta header {seq_ids[0]}")
            return False
    elif len(seq_ids) > 1:
        for id in seq_ids:
            id_parsed = id.split('_')
            if len(id_parsed)!=2 or not id_parsed[1].isdigit():
                print(f"Invalid fasta header {id}")
                return False
    return True

def check_sequence(segment, sequences, max_N_ratio, min_chunk_length, max_segment_length):
    """Check sequence length, characters and N content"""
    summary_length = 0
    for record in sequences:

        # check minimum chunk length
        sequence_length= len(record.seq)
        if sequence_length < min_chunk_length:
            print(f"{record.id} fasta length is too small")
            return False
            
        summary_length+=sequence_length

        # check all letters are IUPAC nucleotides
        iupac_nucls = ''.join(list(IUPACData.ambiguous_dna_values.keys()))
        seq = record.seq.upper()
        if not seq.strip(iupac_nucls) == '':
            print(f"{record.id} fasta contains invalid nucleotide chars")
            return False

        # check N content
        N_num = seq.count('N')
        N_ratio = N_num/sequence_length

        if N_ratio > max_N_ratio:
            print(f"Too much N in {record.id} - {round(N_ratio*100,2)}%")
            return False
    
    # check maximum segment length
    if summary_length > max_segment_length:
        print(f"{segment} fasta summary length is too large")
        return False   

    return True

def check_segment(segment, sequences):
    """Check if sequences came from the corresponding segment"""
    is_seq_correct = []
    for record in sequences:
        if not isseg(record, reference_sequences, segment):
            is_seq_correct.append(False)
        else:
            is_seq_correct.append(True)

    any_seq_correct = any(is_seq_correct)
    all_seq_correct = all(is_seq_correct)

    # all chunks are valid
    if all_seq_correct:
        return 0

    # all chunks are invalid
    if not any_seq_correct:
        return 2
    # some chunks are invalid
    else:
        return 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='The script contains the set of functions needed to check if a sequence is the certain segment of the Influenza A/B genome. Use it with parameter --help for more detailed info.')
    parser.add_argument('input',
        type=str,
        nargs='+', 
        help='Path to 8 fasta files to analyze. Each file corresponds to the Influenza virus segment: segment1.fasta, segment2.fasta...segment8.fasta')
    parser.add_argument('type',
        type=str,
        help='Type of Influenza virus (A or B).',
        choices= ['A', 'B'])
    parser.add_argument('--max_N_ratio', 
        type=float,
        help='Maximum permissible content of N in each sequence , default: 0.5', 
        default=0.5)
    parser.add_argument('--min_chunk_length', 
        type=int,
        help='minimum length for each chunk of segment, default: 50', 
        default=50)
    parser.add_argument('--max_segment_length', 
        type=int,
        help='maximum length of segment, default: 3000', 
        default=3000)

    args = parser.parse_args()

    # Parse reference fasta files
    reference_folder = os.path.join(os.getcwd(),'InfRefs', 'Influenza_'+args.type)
    reference_files = list(map(lambda x: os.path.join(reference_folder, x), os.listdir(reference_folder)))
    reference_sequences = list(map(lambda x: list(SeqIO.parse(x, "fasta")), reference_files)) # list with reference sequences

    # Check if all 8 files for each segment were provided
    if len(args.input)!=8:
        sys.exit("8 fasta files should be provided")

    # Check if fasta headers and sequences are valid
    for i,input_file in enumerate(args.input):
        segment = i+1 # number of segment from 1 to 8
        sequences = list(SeqIO.parse(input_file, "fasta"))

        # If no fasta sequences were found, check that file is empty
        if not check_empty(sequences, input_file):
            sys.exit(f"Invalid fasta for segment {segment}")

        if not sequences:
            continue

        # Check fasta id correctness for each record in file
        seq_ids = [record.id for record in sequences]
        if not check_header(segment, seq_ids):
            sys.exit(f"Invalid header in fasta for segment {segment}")

        # Check sequence length, characters and N content
        if not check_sequence(segment, sequences, args.max_N_ratio, args.min_chunk_length, args.max_segment_length):
            sys.exit(f"Invalid sequence for segment {segment}")

        print(f"Segment {segment} fasta validation was passed!")
        
        # Check if sequences came from the corresponding segment
        status_info = {0:'correct', 1:'partially correct', 2:'invalid'}
        status = check_segment(segment, sequences)
        print(f"Sequences from file {input_file} correspond to segment {segment} reference: {status_info[status]}")