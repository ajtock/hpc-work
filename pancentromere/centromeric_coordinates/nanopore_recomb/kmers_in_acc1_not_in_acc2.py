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

with open(parser.acc1 + "_centromeres.fa_" + str(parser.kmerSize) + "mers.pickle", "rb") as handle:
  acc1_cen = pickle.load(handle)

with open(parser.acc2 + "_centromeres.fa_" + str(parser.kmerSize) + "mers.pickle", "rb") as handle:
  acc2_cen = pickle.load(handle)

with open(parser.acc1 + "_not_centromeres.fa_" + str(parser.kmerSize) + "mers.pickle", "rb") as handle:
  acc1_not_cen = pickle.load(handle)

with open(parser.acc2 + "_not_centromeres.fa_" + str(parser.kmerSize) + "mers.pickle", "rb") as handle:
  acc2_not_cen = pickle.load(handle)


acc1_cen_set = set(list(acc1_cen.keys()))
acc2_cen_set = set(list(acc2_cen.keys()))


dataScientist = set(['Python', 'R', 'SQL', 'Git', 'Tableau', 'SAS'])
dataEngineer = set(['Python', 'Java', 'Scala', 'Git', 'SQL', 'Hadoop'])
graphicDesigner = {'InDesign', 'Photoshop', 'Acrobat', 'Premiere', 'Bridge'}
print(dataScientist)
print(dataEngineer)
print(graphicDesigner)
graphicDesigner.add("Illustrator")
print(graphicDesigner)

graphicDesigner.discard("Premiere")
print(graphicDesigner)

for skill in dataScientist:
    print(skill)

print(type(sorted(dataScientist, reverse=True)))


def remove_duplicates(original):
    unique = []
    [unique.append(n) for n in original if n not in unique]
    return unique

print(remove_duplicates([1, 2, 3, 1, 7]))

print(timeit.timeit("list(set([1, 2, 3, 1, 7]))", number=10000))
print(timeit.timeit("remove_duplicates([1, 2, 3, 1, 7])", globals=globals(), number=10000))


if __name__ == "__main__":
    with open(parser.fasta, "r") as fa_object:
        tic = time()
        kmer_count_dict = count_kmer_fa(k=parser.kmerSize, fa_object=fa_object)
        print(f"Done in {time() - tic:.3f}s")

    with open(parser.fasta + "_" + str(parser.kmerSize) + "mers.pickle", "wb") as handle:
        pickle.dump(kmer_count_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(parser.fasta + "_" + str(parser.kmerSize) + "mers.pickle", "rb") as handle:
        kmer_count_dict_test = pickle.load(handle)

    print(kmer_count_dict == kmer_count_dict_test)
