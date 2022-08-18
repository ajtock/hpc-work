#!/usr/bin/env python3

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 08/08/2022

# Usage:
# ./kmers_in_acc1_not_in_acc2.py -a1 Col-0.ragtag_scaffolds -a2 Ler-0_110x.ragtag_scaffolds -k 178 

# Find and count all possible k-mers (substrings of length k)
# in the centromeres of a given genome
# The bigger the k, the longer it will take

# ==== Import libraries
import sys
import os
import argparse
import pickle

from time import time, sleep
import timeit


# ==== Capture user input as command-line arguments
# https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
def create_parser():
    parser = argparse.ArgumentParser(description="Fasta filename variables.")
    #### Define command-line arguments
    parser.add_argument("-a1", "--acc1", type=str, default="Col-0.ragtag_scaffolds",
                        help="The prefix of the first accession. Default: Col-0.ragtag_scaffolds")
    parser.add_argument("-a2", "--acc2", type=str, default="Ler-0_110x.ragtag_scaffolds",
                        help="The prefix of the second accession. Default: Ler-0_110x.ragtag_scaffolds")
    parser.add_argument("-k", "--kmerSize", type=int, default="178",
                        help="The size of the k-mers to be found and counted in the FASTA file.")
    #### Create parser
    return parser

parser = create_parser().parse_args()
print(parser)


# Load k-mer count dictionaries (saved as pickle files)
with open(parser.acc1 + "_centromeres.fa_" + str(parser.kmerSize) + "mers.pickle", "rb") as handle:
  acc1_cen = pickle.load(handle)

with open(parser.acc2 + "_centromeres.fa_" + str(parser.kmerSize) + "mers.pickle", "rb") as handle:
  acc2_cen = pickle.load(handle)

with open(parser.acc1 + "_not_centromeres.fa_" + str(parser.kmerSize) + "mers.pickle", "rb") as handle:
  acc1_not_cen = pickle.load(handle)

with open(parser.acc2 + "_not_centromeres.fa_" + str(parser.kmerSize) + "mers.pickle", "rb") as handle:
  acc2_not_cen = pickle.load(handle)


# Make keys (k-mers) in each dictionary a set
acc1cen = set(list(acc1_cen.keys()))
acc2cen = set(list(acc2_cen.keys()))
acc1notcen = set(list(acc1_not_cen.keys()))
acc2notcen = set(list(acc2_not_cen.keys()))

print(len(acc1cen))
# 5376061
print(len(acc2cen))
# 5130730

# Define union of acc1notcen, acc2cen and acc2notcen
acc1notcen_acc2cen_acc2notcen_union = acc1notcen.union(acc2cen).union(acc2notcen)
acc1notcen_acc2cen_union = acc1notcen.union(acc2cen)
acc1notcen_acc2cen_acc2notcen_union2 = acc1notcen_acc2cen_union.union(acc2notcen)
acc1notcen_acc2cen_acc2notcen_union == acc1notcen_acc2cen_acc2notcen_union2
del acc1notcen_acc2cen_acc2notcen_union2

# Get k-mers unique to acc1cen by computing difference between
# acc1cen and the above-defined union
# (i.e., acc1cen \ acc1notcen_acc2cen_acc2notcen_union): acc1in
acc1in = acc1cen.difference(acc1notcen_acc2cen_acc2notcen_union)
print(len(acc1in))
# 5114056
print(len(acc1cen) - len(acc1in))
# 262005 
assert(len(acc1in) == len(list(acc1in)))


# Define union of acc2notcen, acc1cen and acc1notcen
acc2notcen_acc1cen_acc1notcen_union = acc2notcen.union(acc1cen).union(acc1notcen)
acc2notcen_acc1cen_union = acc2notcen.union(acc1cen)
acc2notcen_acc1cen_acc1notcen_union2 = acc2notcen_acc1cen_union.union(acc1notcen)
acc2notcen_acc1cen_acc1notcen_union == acc2notcen_acc1cen_acc1notcen_union2
del acc2notcen_acc1cen_acc1notcen_union2

# Get k-mers unique to acc2cen by computing difference between
# acc2cen and the second above-defined union
# (i.e., acc2cen \ acc2notcen_acc1cen_acc1notcen_union): acc2in
acc2in = acc2cen.difference(acc2notcen_acc1cen_acc1notcen_union)
print(len(acc2in))
# 4927462 
print(len(acc2cen) - len(acc2in))
# 203268
assert(len(acc2in) == len(list(acc2in)))

# Output as FASTA to supply to bbduk.sh as input k-mer database file
acc1in_list = 


def write_fasta(seq, outfile):
    out_fasta = open(outfile, "w")

    for s in 



def write_fasta(seq, outfile):
    out_fasta = open(outfile, "w")

    # Look through sequence ids (sorted alphabetically so output file is
    # reproducible).
    for s in sorted(seq.keys()):
        out_fasta.write(">" + s + "\n")
        out_fasta.write(seq[s] + "\n")

    out_fasta.close() 




#if __name__ == "__main__":
#    with open(parser.fasta, "r") as fa_object:
#        tic = time()
#        kmer_count_dict = count_kmer_fa(k=parser.kmerSize, fa_object=fa_object)
#        print(f"Done in {time() - tic:.3f}s")
#
#    with open(parser.fasta + "_" + str(parser.kmerSize) + "mers.pickle", "wb") as handle:
#        pickle.dump(kmer_count_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
#
#    with open(parser.fasta + "_" + str(parser.kmerSize) + "mers.pickle", "rb") as handle:
#        kmer_count_dict_test = pickle.load(handle)
#
#    print(kmer_count_dict == kmer_count_dict_test)
