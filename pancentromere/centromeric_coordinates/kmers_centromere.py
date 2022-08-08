#!/usr/bin/env python3

# Author: Andy Tock (ajt200@cam.ac.uk)
# Mainly based on and adapted from: https://github.com/lh3/kmer-cnt/blob/master/kc-py1.py
# and https://www.youtube.com/watch?v=dQG4-Gwo4BE
# and https://sourmash.readthedocs.io/en/latest/kmers-and-minhash.html
# Date: 08/08/2022

# Find and count all possible k-mers (substrings of length k)
# in the centromeres of a given genome
# The bigger the k, the longer it will take

# ==== Import libraries
import sys
import os
import argparse

from time import time, sleep

# ==== Capture user input as command-line arguments
# https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
def create_parser():
    parser = argparse.ArgumentParser(description="Fasta filename variables.")
    #### Define command-line arguments
    parser.add_argument("-f", "--fasta", type=str, default="Col-0.ragtag_scaffolds_centromeres.fa",
                        help="The filename of the fasta file to be analysed. Default: Col-0.ragtag_scaffolds_centromeres.fa")
    #### Create parser
    return parser

parser = create_parser().parse_args()


# Build a string translation table to enable
# forward-strand sequence to be translated into
# reverse-strand sequence
base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

#fa_object = open(parser.fasta)
#lines = fa_object.readlines()
#
#faRead = fa_object.read()
#faReadSplit = faRead.split()
#
#faReadSplitNoHeaders = []
#
#for i in range(len(faReadSplit)):
#    print(i)
#    if ">" not in faReadSplit[i]:
#        faReadSplitNoHeaders.append(faReadSplit[i])
#
#
#seq = faReadSplitNoHeaders[0]
#k = 178
#h = {}

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
        if line[0] == ">":
            if(len(seq) > 0):
                count_kmer(counter, k, "".join(seq).upper())
                seq = []
        else:
            seq.append(line[:-1])
    if len(seq) > 0:
        count_kmer(counter, k, "".join(seq).upper())
    return counter
 

fa_object = open(parser.fasta)
tic = time()
count_178mers = count_kmer_fa(k = 178, fa_object = fa_object)
print(f"Done in {time() - tic:.3f}s")
 
fa_object = open(parser.fasta)
tic = time()
count_89mers = count_kmer_fa(k = 89, fa_object = fa_object)
print(f"Done in {time() - tic:.3f}s")


kfreq = {}
kmers = []
n_kmers = len(sequence) - ksize + 1





