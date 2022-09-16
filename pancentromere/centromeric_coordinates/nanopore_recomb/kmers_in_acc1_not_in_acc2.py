#!/usr/bin/env python3

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 08/08/2022

# Usage:
# ./kmers_in_acc1_not_in_acc2.py -a1c Col-0.ragtag_scaffolds_centromeres -a2c Ler-0_110x.ragtag_scaffolds_centromeres -a1nc Col-0.ragtag_scaffolds_not_centromeres -a2nc Ler-0_110x.ragtag_scaffolds_not_centromeres -k 24 

# Write accession-specific centromeric k-mers
# to FASTA to supply as input file to bbduk.sh ("ref" parameter)
# containing query sequences to be matched against the FASTQ of ONT reads 

# ==== Import libraries
import sys
import os
import argparse
import screed
import pickle

from matplotlib import pyplot as plt
from matplotlib.pyplot import subplots
from venn import venn
from time import time, sleep
import timeit


# ==== Capture user input as command-line arguments
# https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
def create_parser():
    parser = argparse.ArgumentParser(description="Fasta filename variables.")
    #### Define command-line arguments
    parser.add_argument("-a1c", "--acc1c", type=str, default="Col-0.ragtag_scaffolds_centromeres",
                        help="The prefix of the first accession's centromeric sequences. Default: Col-0.ragtag_scaffolds_centromeres")
    parser.add_argument("-a2c", "--acc2c", type=str, default="Ler-0_110x.ragtag_scaffolds_centromeres",
                        help="The prefix of the second accession's centromeric sequences. Default: Ler-0_110x.ragtag_scaffolds_centromeres")
    parser.add_argument("-a1nc", "--acc1nc", type=str, default="Col-0.ragtag_scaffolds_not_centromeres",
                        help="The prefix of the first accession's non-centromeric sequences. Default: Col-0.ragtag_scaffolds_not_centromeres")
    parser.add_argument("-a2nc", "--acc2nc", type=str, default="Ler-0_110x.ragtag_scaffolds_not_centromeres",
                        help="The prefix of the second accession's non-centromeric sequences. Default: Ler-0_110x.ragtag_scaffolds_not_centromeres")
    parser.add_argument("-k", "--kmerSize", type=int, default="24",
                        help="The size of the k-mers to be found and counted in the FASTA file.")
    #### Create parser
    return parser

parser = create_parser().parse_args()
print(parser)


outDir = "fasta"
plotDir = outDir + "/plots"

if not os.path.exists(outDir):
    os.makedirs(outDir)

if not os.path.exists(plotDir):
    os.makedirs(plotDir)


# Define a list containing the union of
# elements in an arbitrary number of lists 
def union_lists(*lists):
    """
    Get the union of elements in an arbitrary number of lists
    """
    return list(set.union(*map(set, lists)))


# Write dictionary of accession-specific centromeric k-mers
# to FASTA to supply to bbduk.sh as input k-mer database file
def write_fasta(kmer_dict, acc_name, outfile):
    """
    Write dictionary of k-mers to FASTA
    """
    with open(outfile, "w") as fa_object:
        for s in kmer_dict.keys():
            fa_object.write(">" + str(s) + "_" + acc_name + "\n")
            fa_object.write(kmer_dict[s] + "\n")


# Convert a k-mer into a hash
# See https://sourmash.readthedocs.io/en/latest/kmers-and-minhash.html
# But using hash() instead of mmh3.hash64() because since version 3.4,
# Python hash() uses SipHash, which is more secure and less vulnerable
# to hash collision attacks
def hash_kmer(kmer):
    """
    Convert a k-mer into a hash.
    """
    # Get the reverse complement
    rc_kmer = screed.rc(kmer)
    #
    # Get the lesser of kmer and rc_kmer, based on lexicographical order
    if kmer < rc_kmer:
        canonical_kmer = kmer
    else:
        canonical_kmer = rc_kmer
    #
    # Calculate hash
    kmer_hash = hash(canonical_kmer)
    #
    return kmer_hash

# Convert a list of k-mers into a list of hashes
# See https://sourmash.readthedocs.io/en/latest/kmers-and-minhash.html
def hash_kmers(kmers):
    """
    Convert a list of k-mers into a list of hashes.
    """
    kmer_hashes = []
    for kmer in kmers:
        kmer_hashes.append(hash_kmer(kmer))
    #
    return kmer_hashes


# Load k-mer count dictionaries (saved as pickle files)
with open("kmers_in_fasta_dict/" + parser.acc1c + ".fa_k" + str(parser.kmerSize) + ".pickle", "rb") as handle:
  acc1cen_dict = pickle.load(handle)

with open("kmers_in_fasta_dict/" + parser.acc2c + ".fa_k" + str(parser.kmerSize) + ".pickle", "rb") as handle:
  acc2cen_dict = pickle.load(handle)

with open("kmers_in_fasta_dict/" + parser.acc1nc + ".fa_k" + str(parser.kmerSize) + ".pickle", "rb") as handle:
  acc1notcen_dict = pickle.load(handle)

with open("kmers_in_fasta_dict/" + parser.acc2nc + ".fa_k" + str(parser.kmerSize) + ".pickle", "rb") as handle:
  acc2notcen_dict = pickle.load(handle)

with open("kmers_in_fasta_dict/t2t-col.20210610_ChrM_ChrC.fa_k" + str(parser.kmerSize) + ".pickle", "rb") as handle:
  mitochloro_dict = pickle.load(handle)


# Make keys (k-mers) in each dictionary a list
acc1cen = list(acc1cen_dict.keys())
acc2cen = list(acc2cen_dict.keys())
acc1notcen = list(acc1notcen_dict.keys())
acc2notcen = list(acc2notcen_dict.keys())
mitochloro = list(mitochloro_dict.keys())

print(len(acc1cen))
# 660531
print(len(acc2cen))
# 659434 
print(len(acc1notcen))
# 109988424
print(len(acc2notcen))
# 110021023
print(len(mitochloro))
# 475021

## Convert each list of k-mers into a list of hashes
#acc1cen_hashes = hash_kmers(acc1cen)
#acc2cen_hashes = hash_kmers(acc2cen)
#acc1notcen_hashes = hash_kmers(acc1notcen)
#acc2notcen_hashes = hash_kmers(acc2notcen)
#mitochloro_hashes = hash_kmers(mitochloro)


# Make venn to show k-mer overlap
cmap = "plasma"

#hash_dict = {
#    "Col-0 cen": set(acc1cen_hashes),
#    "Ler-0 cen": set(acc2cen_hashes),
#    "Col-0 arm": set(acc1notcen_hashes),
#    "Ler-0 arm": set(acc2notcen_hashes)
#} 

dataset_dict = {
    "Col-0 cen": set(acc1cen),
    "Ler-0 cen": set(acc2cen),
    "Col-0 arm": set(acc1notcen),
    "Ler-0 arm": set(acc2notcen)
} 


venn_plot = venn(dataset_dict, cmap="plasma")
venn_fig = venn_plot.get_figure()
venn_fig.savefig(plotDir + "/venn.png")

#hash_venn_plot = venn(hash_dict, cmap="plasma")
#hash_venn_fig = venn_plot.get_figure()
#hash_venn_fig.savefig(plotDir + "/venn_hash.png")


# Define union of acc1notcen, acc2cen, acc2notcen and mitochloro
acc1notcen_acc2cen_acc2notcen_mitochloro_union = union_lists(acc1notcen, acc2cen, acc2notcen, mitochloro)

# Get k-mers unique to acc1cen by computing difference between
# acc1cen and the above-defined union
# (i.e., acc1cen \ acc1notcen_acc2cen_acc2notcen_mitochloro_union): acc1in
acc1in = list(set(acc1cen).difference(set(acc1notcen_acc2cen_acc2notcen_mitochloro_union)))
acc1in = sorted(acc1in)
print(len(acc1in))
# 309849
print(len(acc1cen) - len(acc1in))
# 350682

# Convert into dictionary
acc1in_1tolen = list(range(1, len(acc1in)+1)) 
acc1in_dict = dict(zip(acc1in_1tolen, acc1in))


# Define union of acc2notcen, acc1cen, acc1notcen and mitochloro
acc2notcen_acc1cen_acc1notcen_mitochloro_union = union_lists(acc2notcen, acc1cen, acc1notcen, mitochloro)

# Get k-mers unique to acc2cen by computing difference between
# acc2cen and the second above-defined union
# (i.e., acc2cen \ acc2notcen_acc1cen_acc1notcen_mitochloro_union): acc2in
acc2in = list(set(acc2cen).difference(set(acc2notcen_acc1cen_acc1notcen_mitochloro_union)))
acc2in = sorted(acc2in)
print(len(acc2in))
# 353341 
print(len(acc2cen) - len(acc2in))
# 306093

# Convert into dictionary
acc2in_1tolen = list(range(1, len(acc2in)+1)) 
acc2in_dict = dict(zip(acc2in_1tolen, acc2in))


# Write to FASTA
write_fasta(kmer_dict=acc1in_dict,
            acc_name=parser.acc1c[0:5],
            outfile=outDir + "/" + parser.acc1c + "_specific_k" + str(parser.kmerSize) + ".fa")

write_fasta(kmer_dict=acc2in_dict,
            acc_name=parser.acc2c[0:5],
            outfile=outDir + "/" + parser.acc2c + "_specific_k" + str(parser.kmerSize) + ".fa")
