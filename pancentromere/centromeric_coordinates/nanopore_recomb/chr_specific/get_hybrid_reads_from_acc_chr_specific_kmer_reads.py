#!/usr/bin/env python3

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 28/10/2022

# Usage:
# conda activate python_3.9.6
# ./get_hybrid_reads_from_acc_chr_specific_kmer_reads.py \
#  -r Col_Ler_F1_pollen_500bp_minq99 \
#  -a1 Col-0.ragtag_scaffolds_not_centromere \
#  -a2 Ler-0_110x.ragtag_scaffolds_not_centromere \
#  -c Chr1 \
#  -k 24 \
#  -op 0.9 \
#  -mh 10
# conda deactivate

# For each "hybrid" read containing acc1- AND acc2-specific k-mers,
# get the within-read locations of each matching k-mer, and output
# the read segments that span consecutive acc1-specific k-mers and, separately,
# the read segments that span consecutive acc2-specific k-mers in
# FASTQ format for alignment to the respective assemblies.

# ==== Import libraries
#import sys
import os
import glob
import argparse
##import pickle
import re
#import numpy as np
#import pandas as pd
import subprocess
#import screed

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser
from pathlib import Path
#from numba import jit
#from time import time, sleep
##import timeit

# ==== Capture user input as command-line arguments
# https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
def create_parser():
    parser = argparse.ArgumentParser(description="Fasta filename variables.")
    #### Define command-line arguments
    parser.add_argument("-r", "--readsPrefix", type=str, default="Col_Ler_F1_pollen_500bp_minq99",
                        help="The prefix of the FASTQ file name. Default: Col_Ler_F1_pollen_500bp_minq99")
    parser.add_argument("-a1", "--acc1", type=str, default="Col-0.ragtag_scaffolds_not_centromere",
                        help="The prefix of the first accession's sequences. Default: Col-0.ragtag_scaffolds_not_centromere")
    parser.add_argument("-a2", "--acc2", type=str, default="Ler-0_110x.ragtag_scaffolds_not_centromere",
                        help="The prefix of the second accession's sequences. Default: Ler-0_110x.ragtag_scaffolds_not_centromere")
    #parser.add_argument("-reg", "--region", type=str, default="not_centromere",
    #                    help="Region from which to get accession-specific, chromosome-specific k-mers. Default: not_centromere")
    parser.add_argument("-c", "--chrom", type=str, default="Chr1",
                        help="Name of chromosome from which to get accession-specific, chromosome-specific k-mers. Default: Chr1")
    parser.add_argument("-k", "--kmerSize", type=int, default="24",
                        help="The size of the k-mers to be found and counted in the FASTA file. Default: 24")
    parser.add_argument("-op", "--overlapProp", type=float, default="0.9",
                        help="The minimum proportion of an aligned k-mer's length that must overlap a genomic window for the aligned k-mer to be kept during downsampling of accession-specific k-mers. Default: 0.9")
    parser.add_argument("-mh", "--minHits", type=int, default="10",
                        help="The minimum number of accession-specific k-mers found in a read. Default: 10")
    return parser

parser = create_parser().parse_args()
print(parser)

region = re.sub(".+_scaffolds_", "", parser.acc1)
indir = "fasta"
outdir = region + "/" + parser.chrom + "/fasta"

if not os.path.exists(outdir):
    os.makedirs(outdir)


# Existing paths to reads containing accession-specific, chromosome-specific k-mers
acc1_fa = indir + "/" + parser.readsPrefix + \
    "_match_" + parser.acc1 + \
    "_" + parser.chrom + \
    "_specific_k" + str(parser.kmerSize) + \
    "_downsampled_op" + str(parser.overlapProp) + \
    "_hits" + str(parser.minHits) + \
    ".fa"
# File exists sanity check
Path(acc1_fa).resolve(strict=True)

acc2_fa = indir + "/" + parser.readsPrefix + \
    "_match_" + parser.acc2 + \
    "_" + parser.chrom + \
    "_specific_k" + str(parser.kmerSize) + \
    "_downsampled_op" + str(parser.overlapProp) + \
    "_hits" + str(parser.minHits) + \
    ".fa"
# File exists sanity check
Path(acc2_fa).resolve(strict=True)

# To-be-created path to hybrid reads containing accession-specific, chromosome-specific
# k-mers from both accessions
# NOTE: remember to symlink indir location to outdir
hybrid_fa = indir + "/" + parser.readsPrefix + \
    "_match_" + parser.acc1 + \
    "_" + parser.chrom + \
    "_specific_k" + str(parser.kmerSize) + \
    "_downsampled_op" + str(parser.overlapProp) + \
    "_hits" + str(parser.minHits) + \
    "_match_" + parser.acc2 + \
    "_" + parser.chrom + \
    "_specific_k" + str(parser.kmerSize) + \
    "_downsampled_op" + str(parser.overlapProp) + \
    "_hits" + str(parser.minHits) + \
    "_test.fa"


# Make a list containing the intersection of
# elements in an arbitrary number of lists
def intersection_lists(*lists):
    """
    Get the intersection of elements in an arbitrary number of lists.
    """
    return list(set.intersection(*map(set, lists)))


acc1_fa_dict = SeqIO.index(acc1_fa, "fasta")
acc2_fa_dict = SeqIO.index(acc2_fa, "fasta")

acc1_fa_id_list = [v.id for i, v in enumerate(acc1_fa_dict.values())]
acc2_fa_id_list = [v.id for i, v in enumerate(acc2_fa_dict.values())]

hybrid_fa_id_list = intersection_lists(acc1_fa_id_list, acc2_fa_id_list)
hybrid_fa_gen = (v for i, v in enumerate(acc1_fa_dict.values()) if v.id in hybrid_fa_id_list)

with open(hybrid_fa, "w") as hybrid_fa_handle:
    for record in hybrid_fa_gen:
        SeqIO.write(record, hybrid_fa_handle, "fasta")

os.chdir(outdir)
subprocess.run(["ln", "-s", "../../../" + hybrid_fa, "."])
os.chdir("../../../")
