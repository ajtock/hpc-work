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


# Function to define a list containing the union of
#  elements in an arbitrary number of lists 
def union_lists(*lists):
    """
    Get the union of elements in an arbitrary number of lists
    """
    return list(set.union(*map(set, lists)))


# Function to write dictionary of accession-specific centromeric k-mers
# to FASTA to supply to bbduk.sh as input k-mer database file
def write_fasta(kmer_dict, outfile):
    """
    Write dictionary of k-mers to FASTA
    """
    with open(outfile, "w") as fa_object:
        for s in kmer_dict.keys():
            fa_object.write(">" + str(s) + "\n")
            fa_object.write(kmer_dict[s] + "\n")


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


# Make keys (k-mers) in each dictionary a list
acc1cen = list(acc1_cen.keys())
acc2cen = list(acc2_cen.keys())
acc1notcen = list(acc1_not_cen.keys())
acc2notcen = list(acc2_not_cen.keys())

print(len(acc1cen))
# 5376061
print(len(acc2cen))
# 5130730


# Define union of acc1notcen, acc2cen and acc2notcen
acc1notcen_acc2cen_acc2notcen_union = union_lists(acc1notcen, acc2cen, acc2notcen)

# Get k-mers unique to acc1cen by computing difference between
# acc1cen and the above-defined union
# (i.e., acc1cen \ acc1notcen_acc2cen_acc2notcen_union): acc1in
acc1in = list(set(acc1cen).difference(set(acc1notcen_acc2cen_acc2notcen_union)))
acc1in = sorted(acc1in)
print(len(acc1in))
# 5114056
print(len(acc1cen) - len(acc1in))
# 262005

# Convert into dictionary
acc1in_1tolen = list(range(1, len(acc1in)+1)) 
acc1in_dict = dict(zip(acc1in_1tolen, acc1in))


# Define union of acc2notcen, acc1cen and acc1notcen 
acc2notcen_acc1cen_acc1notcen_union = union_lists(acc2notcen, acc1cen, acc1notcen)

# Get k-mers unique to acc2cen by computing difference between
# acc2cen and the second above-defined union
# (i.e., acc2cen \ acc2notcen_acc1cen_acc1notcen_union): acc2in
acc2in = list(set(acc2cen).difference(set(acc2notcen_acc1cen_acc1notcen_union)))
acc2in = sorted(acc2in)
print(len(acc2in))
# 4927462 
print(len(acc2cen) - len(acc2in))
# 203268

# Convert into dictionary
acc2in_1tolen = list(range(1, len(acc2in)+1)) 
acc2in_dict = dict(zip(acc2in_1tolen, acc2in))


# Write to FASTA
write_fasta(kmer_dict=acc1in_dict,
            outfile=parser.acc1 + "_centromere_specific_" + str(parser.kmerSize) + "mers.fasta")

write_fasta(kmer_dict=acc2in_dict,
            outfile=parser.acc2 + "_centromere_specific_" + str(parser.kmerSize) + "mers.fasta")
