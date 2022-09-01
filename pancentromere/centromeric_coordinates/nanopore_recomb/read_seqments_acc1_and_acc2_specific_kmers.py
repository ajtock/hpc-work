#!/usr/bin/env python3

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 25/08/2022

# Usage:
# conda activate python_3.9.6
# ./read_seqments_acc1_and_acc2_specific_kmers.py -r Col_ler_f1_pollen_500bp_minq99 \
#  -a1 Col-0.ragtag_scaffolds_centromeres \
#  -a2 Ler-0_110x.ragtag_scaffolds_centromeres \ 
#  -k 24
# conda deactivate

# For each "hybrid" read containing acc1- AND acc2-specific k-mers,
# get the within-read locations of each matching k-mer, and output
# the read segments that span consecutive acc1-specific k-mers and, separately,
# the read segments that span consecutive acc2-specific k-mers in
# FASTQ format for alignment to the respective assemblies.

# ==== Import libraries
import sys
import os
import argparse
import pickle
import re

from Bio import SeqIO
from time import time, sleep
import timeit

outDir = "segments_fastq"

if not os.path.exists(outDir):
    os.makedirs(outDir)

# ==== Capture user input as command-line arguments
# https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
def create_parser():
    parser = argparse.ArgumentParser(description="Fasta filename variables.")
    #### Define command-line arguments
    parser.add_argument("-r", "--readsPrefix", type=str, default="Col_ler_f1_pollen_500bp_minq99",
                        help="The prefix of the FASTQ file name. Default: Col_ler_f1_pollen_500bp_minq99")
    parser.add_argument("-a1", "--acc1", type=str, default="Col-0.ragtag_scaffolds_centromeres",
                        help="The prefix of the first accession's sequences. Default: Col-0.ragtag_scaffolds_centromeres")
    parser.add_argument("-a2", "--acc2", type=str, default="Ler-0_110x.ragtag_scaffolds_centromeres",
                        help="The prefix of the second accession's sequences. Default: Ler-0_110x.ragtag_scaffolds_centromeres")
    parser.add_argument("-k", "--kmerSize", type=int, default="24",
                        help="The size of the k-mers to be found and counted in the FASTA file.")
    #### Create parser
    return parser

parser = create_parser().parse_args()
print(parser)


# Path to hybrid reads
input_fa = "fasta/" + parser.readsPrefix + \
    "_match_" + parser.acc1 + \
    "_specific_k" + str(parser.kmerSize) + \
    "_match_" + parser.acc2 + \
    "_specific_k" + str(parser.kmerSize) + \
    ".fa"

# Path to acc1-specific k-mers
acc1_fa = "fasta/" + \
    parser.acc1 + \
    "_specific_k" + str(parser.kmerSize) + \
    ".fa"
# Path to acc2-specific k-mers
acc2_fa = "fasta/" + \
    parser.acc2 + \
    "_specific_k" + str(parser.kmerSize) + \
    ".fa"

# Parse reads as FastaIterator
reads = list(SeqIO.parse(input_fa, "fasta")


acc1_kmers = list(SeqIO.parse(acc1_fa, "fasta"))


            outfile=outDir + "/" + parser.acc1c + "_specific_k" + str(parser.kmerSize) + ".fa")


# Find the 0-based start location in each read of all occurrences
# of each k-mer
print([match.start() for match in re.finditer(pattern, string)])

with open(




# Build a string translation table to enable
# forward-strand sequence to be translated into
# reverse-strand sequence
base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)


def count_kmer(h, k, seq):
    l = len(seq)
    if l < k: return
    for i in range(l - k + 1):
        kmer_for = seq[i:(i + k)]
        if "N" in kmer_for: continue
        kmer_rev = kmer_for.translate(comp_tab)[::-1]
        if kmer_for < kmer_rev:
            kmer = kmer_for
        else:
            kmer = kmer_rev
        if kmer in h:
            h[kmer] += 1
        else:
            h[kmer] = 1


def count_kmer_fa(k, fa_object):
    counter = {}
    seq = []
    lines = fa_object.readlines()
    for line in lines:
        # If line is FASTA header and
        # if seq contains sequence below previous FASTA header
        if line[0] == ">":
            if(len(seq) > 0):
                # Count kmer occurrences in that sequence
                count_kmer(counter, k, "".join(seq).upper())
                # Make sequence empty again so that next line in FASTA
                # can be processed in the same way
                seq = []
        # Add sequence under FASTA header to seq list
        else:
            seq.append(line[:-1])
    # Process final line of sequence in FASTA
    if len(seq) > 0:
        count_kmer(counter, k, "".join(seq).upper())
    return counter


if __name__ == "__main__":
    with open("fasta/" + parser.fasta, "r") as fa_object:
        tic = time()
        kmer_count_dict = count_kmer_fa(k=parser.kmerSize, fa_object=fa_object)
        print(f"Done in {time() - tic:.3f}s")

    with open(outDir + "/" + parser.fasta + "_k" + str(parser.kmerSize) + ".pickle", "wb") as handle:
        pickle.dump(kmer_count_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(outDir + "/" + parser.fasta + "_k" + str(parser.kmerSize) + ".pickle", "rb") as handle:
        kmer_count_dict_test = pickle.load(handle)

    print(kmer_count_dict == kmer_count_dict_test)
