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
#import sys
import os
import argparse
#import pickle
import re
import numpy as np
import pandas as pd

from Bio import SeqIO
from pathlib import Path
#from time import time, sleep
#import timeit

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
                        help="The size of the k-mers to be found and counted in the FASTA file. Default: 24")
    parser.add_argument("-mh", "--minHits", type=int, default="3",
                        help="The minimum number of accession-specific k-mers found in a read. Default: 3")
    #### Create parser
    return parser

parser = create_parser().parse_args()
print(parser)

acc1 = parser.acc1.split(".")[0].split("_")[0]
acc2 = parser.acc2.split(".")[0].split("_")[0]

# Path to hybrid reads
input_fa = "fasta/" + parser.readsPrefix + \
    "_match_" + parser.acc1 + \
    "_specific_k" + str(parser.kmerSize) + \
    "_hits" + str(parser.minHits) + \
    "_match_" + parser.acc2 + \
    "_specific_k" + str(parser.kmerSize) + \
    "_hits" + str(parser.minHits) + \
    ".fa"
# File exists sanity check
Path(input_fa).resolve(strict=True)

# Path to acc1-specific k-mers
acc1_fa = "fasta/" + \
    parser.acc1 + \
    "_specific_k" + str(parser.kmerSize) + \
    ".fa"
# File exists sanity check
Path(acc1_fa).resolve(strict=True)

# Path to acc2-specific k-mers
acc2_fa = "fasta/" + \
    parser.acc2 + \
    "_specific_k" + str(parser.kmerSize) + \
    ".fa"
# File exists sanity check
Path(acc2_fa).resolve(strict=True)

# Parse reads as FastaIterator
reads = list(SeqIO.parse(input_fa, "fasta"))

# Parse acc1-specific k-mers as FastaIterator
acc1_kmers = list(SeqIO.parse(acc1_fa, "fasta"))

# Parse acc2-specific k-mers as FastaIterator
acc2_kmers = list(SeqIO.parse(acc2_fa, "fasta"))


# Build a string translation table to enable
# forward-strand sequence to be translated into
# reverse-strand sequence
base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)


# Function to define a list containing the union of
# elements in an arbitrary number of lists
def union_lists(*lists):
    """
    Get the union of elements in an arbitrary number of lists.
    """
    return list(set.union(*map(set, lists)))


# Within a read, find the 0-based start location of all occurrences
# of each accession-specific k-mer
def get_kmer_loc(kmers, read):
    """
    For a given read, get the within-read 0-based start locations of all k-mer matches.
    """
    kmer_loc_dict_list = []
    for j in range(len(kmers)):
        kmer_id = kmers[j].id
        kmer_acc = kmers[j].id.split("_", 1)[1]
        kmer_for = str(kmers[j].seq)
        kmer_rev = kmer_for.translate(comp_tab)[::-1] 
        kmer_for_matches = [match.start() for match in re.finditer(kmer_for, read)]
        kmer_rev_matches = [match.start() for match in re.finditer(kmer_rev, read)]
        kmer_matches = sorted(union_lists(kmer_for_matches, kmer_rev_matches))
        if kmer_for < kmer_rev:
            kmer = kmer_for
        else:
            kmer = kmer_rev
        if kmer not in kmer_loc_dict_list:
            if kmer_matches:
                for k in range(len(kmer_matches)):
                    kmer_loc_dict_list.append({"kmer": kmer,
                                               "id": kmer_id,
                                               "acc": kmer_acc,
                                               "hit_start": kmer_matches[k]})
        else:
            print("k-mer already present in object")
    #
    return pd.DataFrame(kmer_loc_dict_list)

# Get the within-read start locations of accession-specific k-mer matches
acc1_kmer_loc_df = get_kmer_loc(kmers=acc1_kmers, read=str(reads[0].seq)) 
acc2_kmer_loc_df = get_kmer_loc(kmers=acc2_kmers, read=str(reads[0].seq)) 

# Concatenate and sort by k-mer match start location in read
acc_kmer_loc_df = pd.concat(objs=[acc1_kmer_loc_df, acc2_kmer_loc_df],
                            axis=0,
                            ignore_index=True)
acc_kmer_loc_df_sort = acc_kmer_loc_df.sort_values(by="hit_start",
                                                   axis=0,
                                                   ascending=True,
                                                   kind="quicksort",
                                                   ignore_index=True)

# Get rows that correspond to accession-specific read segments and
# determine which segment from each accession is the largest
acc1_segs_counter = 0
acc2_segs_counter = 0
acc1_segs_lol = []
acc2_segs_lol = []
for rowtup in acc_kmer_loc_df_sort.itertuples():

rowtup = next(acc_kmer_loc_df_sort.itertuples())
if rowtup.acc == acc1:
    acc1_segs_lol.append(rowtup
else:
    rowtup_acc2 = rowtup


acc1_kmer_loc_df = pd.DataFrame.from_dict(acc1_kmer_loc)


def flatten(lol):
    """
    Use list comprehension to flatten a list of lists (lol) into a
    single list composed of all the elements in each sublist.
    """
    return [item for sublist in lol for item in sublist]


acc1_kmer_loc_values = sorted(flatten(list(acc1_kmer_loc.values())))
acc2_kmer_loc_values = sorted(flatten(list(acc2_kmer_loc.values())))



list(set.union(*map(set, lists)))
sorted(acc1_kmer_for_matches, acc1_kmer_rev_matches)
print([match.start() for match in re.finditer(pattern, string)])



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
