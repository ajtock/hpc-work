#!/usr/bin/env python3

# Mainly based on and adapted from: https://github.com/lh3/kmer-cnt/blob/master/kc-py1.py
# Also see https://sourmash.readthedocs.io/en/latest/kmers-and-minhash.html
# Date: 08/08/2022

# Usage:
# conda activate python_3.9.6
# ./kmers_in_fasta.py -f Col-0.ragtag_scaffolds_centromere_Chr1.fa -k 24
# ./kmers_in_fasta.py -f Ler-0_110x.ragtag_scaffolds_centromere_Chr1.fa -k 24
# ./kmers_in_fasta.py -f Col-0.ragtag_scaffolds_not_centromere_Chr1.fa -k 24
# ./kmers_in_fasta.py -f Ler-0_110x.ragtag_scaffolds_not_centromere_Chr1.fa -k 24
# ./kmers_in_fasta.py -f t2t-col.20210610_ChrM_ChrC.fa -k 24
# conda deactivate

# Find and count all possible k-mers (substrings of length k)
# in the centromere or not_centromere of a given chromosome
# The bigger the k, the longer it will take

# ==== Import libraries
import sys
import os
import argparse
import pickle

from time import time, sleep

outDir = "kmers_in_fasta_dict"

if not os.path.exists(outDir):
    os.makedirs(outDir)


# ==== Capture user input as command-line arguments
# https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
def create_parser():
    parser = argparse.ArgumentParser(description="Fasta filename and k-mer size variables.")
    #### Define command-line arguments
    parser.add_argument("-f", "--fasta", type=str, default="Col-0.ragtag_scaffolds_centromere_Chr1.fa",
                        help="The filename of the fasta file to be analysed. Default: Col-0.ragtag_scaffolds_centromere_Chr1.fa")
    parser.add_argument("-k", "--kmerSize", type=int, default="24",
                        help="The size of the k-mers to be found and counted in the FASTA file. Default: 24")
    #### Create parser
    return parser

parser = create_parser().parse_args()
print(parser)


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
