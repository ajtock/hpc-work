#!/usr/bin/env python3

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 16/10/2022

# Usage:
# conda activate python_3.9.6
# nohup ./kmers_in_acc1_not_in_acc2_dedup_kmers_check_members.py -a Col-0.ragtag_scaffolds_not_centromeres -k 24 -op 0.9 &> nohup1.out&
# conda deactivate

# Write accession-specific centromeric k-mers
# to FASTA to supply as input file to bbduk.sh ("ref" parameter)
# containing query sequences to be matched against the FASTQ of ONT reads 

# ==== Import libraries
import sys
import os
import argparse
import screed
import pickle
import re
import subprocess
import pandas as pd
import seaborn as sns

from Bio import SeqIO
from pybedtools import BedTool
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
    parser.add_argument("-a", "--acc", type=str, default="Col-0.ragtag_scaffolds_not_centromeres",
                        help="The prefix of the first accession's centromeric sequences. Default: Col-0.ragtag_scaffolds_not_centromeres")
    parser.add_argument("-k", "--kmerSize", type=int, default="24",
                        help="The size of the k-mers to be found and counted in the FASTA file. Default: 24")
    parser.add_argument("-op", "--overlapProp", type=float, default="0.9",
                        help="The minimum proportion of an aligned k-mer's length that must overlap a genomic window for the aligned k-mer to be kept during downsampling of accession-specific k-mers. Default: 0.9")
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


# Make a list containing the union of
# elements in an arbitrary number of lists 
def union_lists(*lists):
    """
    Get the union of elements in an arbitrary number of lists.
    """
    return list(set.union(*map(set, lists)))


# Make a list containing the intersection of
# elements in an arbitrary number of lists
def intersection_lists(*lists):
    """
    Get the intersection of elements in an arbitrary number of lists.
    """
    return list(set.intersection(*map(set, lists)))


# Write dictionary of accession-specific centromeric k-mers
# to FASTA to supply to bbduk.sh as input k-mer database file
def write_fasta(kmer_dict, acc_name, outfile):
    """
    Write dictionary of k-mers to FASTA.
    """
    with open(outfile, "w") as fa_object:
        for s in kmer_dict.keys():
            fa_object.write(">" + str(s) + "_" + acc_name + "\n")
            fa_object.write(kmer_dict[s] + "\n")


# Build a string translation table to enable
# forward-strand sequence to be translated into
# reverse-strand sequence
base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)


# Deduplicate downsampled accession-specific k-mers, and
# keep the strand representation of each k-mer that is
# lexicographically smallest, as was done for full k-mer set,
# enabling subsequent test for membership of full set
def dedup_kmers_fa(kmers_fa_noheaders):
    #kmers_fa_noheaders=outDir + "/" + \
    #    parser.acc + "_specific_k" + \
    #    str(parser.kmerSize) + "_bowtie_sorted_intersect_op" + \
    #    str(parser.overlapProp) + "_merge_omg_noheaders.fa"
    """
    Deduplicate downsampled accession-specific k-mers,
    keeping the lexicographically smallest strand representation.
    """
    kmers_list = []
    with open(kmers_fa_noheaders, "r") as kmers_fa_noheaders_handle:
        for line in kmers_fa_noheaders_handle:
            #if line[0] == ">": continue
            kmer_for = line[:-1]
            #kmer_rev = screed.rc(kmer_for)
            kmer_rev = kmer_for.translate(comp_tab)[::-1]
            if kmer_for < kmer_rev:
                kmer = kmer_for
            else:
                kmer = kmer_rev
            if kmer not in kmers_list:
                kmers_list.append(kmer)
    #
    return kmers_list


#kmers_fa_noheaders=outDir + "/" + \
#    parser.acc + "_specific_k" + \
#    str(parser.kmerSize) + "_bowtie_sorted_intersect_op" + \
#    str(parser.overlapProp) + "_merge_omg_head200000.fa"
#start = time()
#acc_kmers_ds = dedup_kmers_fa(kmers_fa_noheaders=kmers_fa_noheaders)
#print(f"Done in {time() - start:.3f}s")
#
#kmers_fa=outDir + "/" + \
#    parser.acc + "_specific_k" + \
#    str(parser.kmerSize) + "_bowtie_sorted_intersect_op" + \
#    str(parser.overlapProp) + "_merge_omg_head200000.fa"
#start = time()
#acc_kmers_ds = dedup_kmers_fa(kmers_fa=kmers_fa)
#print(f"Done in {time() - start:.3f}s")


# Remove headers from full accession-specific k-mers FASTA
def remove_headers(kmers_fa):
    #kmers_fa=outDir + "/" + \
    #    parser.acc + "_specific_k" + \
    #    str(parser.kmerSize) + ".fa"
    """
    Remove headers from full accession-specific k-mers FASTA.
    """
    out_fa_noheaders = re.sub(".fa", "_noheaders.fa", kmers_fa)
    out_fa_noheaders_err = re.sub(".fa", "_noheaders.err", kmers_fa)
    noheaders_cmd = ["grep", "-v"] + \
                    ["^>", kmers_fa]
    with open(out_fa_noheaders, "w") as out_fa_noheaders_handle, \
        open(out_fa_noheaders_err, "w") as out_fa_noheaders_err_handle:
        subprocess.run(noheaders_cmd, stdout=out_fa_noheaders_handle, stderr=out_fa_noheaders_err_handle)
        # Delete empty error files
        if os.stat(out_fa_noheaders_err).st_size == 0:
            subprocess.run(["rm", out_fa_noheaders_err])


# Check the downsampled (ds) kmers in the list output
# from dedup_kmers_fa() (e.g., acc_kmers_ds) for membership of
# the corresponding full accession-specific k-mer set (e.g., acc_kmers),
# returning members as a sorted list
def get_members(ds_kmers_list, full_kmers_list):
    #ds_kmers_list=acc_kmers_ds
    #full_kmers_list=acc_kmers
    """
    Make a sorted list containing the downsampled k-mers in the list output
    from dedup_kmers_fa() that are members of the corresponding full
    accession-specific k-mer set.
    """
    print("k-mers in ds_kmers_list: " + str(len(ds_kmers_list)))
    members_ds_kmers_list = sorted(intersection_lists(ds_kmers_list, full_kmers_list))
    print("k-mers in members_ds_kmers_list: " + str(len(members_ds_kmers_list)))
    #
    return members_ds_kmers_list



def main():
    """
    Get accession-specific and downsampled accession-specific k-mers
    and write to FASTA files.
    """
    
    # Deduplicate downsampled accession-specific k-mers,
    # keeping the lexicographically smallest strand representation,
    # returned as a list
    # NOTE: long run time
    acc_kmers_ds = dedup_kmers_fa(
        kmers_fa_noheaders=outDir + "/" + \
            parser.acc + "_specific_k" + \
            str(parser.kmerSize) + "_bowtie_sorted_intersect_op" + \
            str(parser.overlapProp) + "_merge_omg_noheaders.fa")
    
    ## Get full list of accession-specific k-mers    
    full_kmers_fa = outDir + "/" + \
            parser.acc + "_specific_k" + \
            str(parser.kmerSize) + ".fa"
    remove_headers(
        kmers_fa=full_kmers_fa)
    full_kmers_fa_noheaders = re.sub(".fa", "_noheaders.fa", full_kmers_fa)
    with open(full_kmers_fa_noheaders, "r") as full_kmers_fa_noheaders_handle:
        acc_kmers = full_kmers_fa_noheaders_handle.read().splitlines()
    del full_kmers_fa_noheaders_handle, full_kmers_fa_noheaders, full_kmers_fa
    
    # Make a sorted list containing the downsampled k-mers in the list output
    # from dedup_kmers_fa() that are members of the corresponding full
    # accession-specific k-mer set
    acc_kmers_ds_members = get_members(
        ds_kmers_list=acc_kmers_ds,
        full_kmers_list=acc_kmers)
    del acc_kmers, acc_kmers_ds
    # Convert into dictionary
    acc_kmers_ds_members_1tolen = list(range(1, len(acc_kmers_ds_members)+1))
    acc_kmers_ds_members_dict = dict(zip(acc_kmers_ds_members_1tolen, acc_kmers_ds_members))
    del acc_kmers_ds_members, acc_kmers_ds_members_1tolen
    # Write to FASTA
    write_fasta(
        kmer_dict=acc_kmers_ds_members_dict,
        acc_name=parser.acc[0:5],
        outfile=outDir + "/" + \
            parser.acc + "_specific_k" + \
            str(parser.kmerSize) + "_downsampled_op" + \
            str(parser.overlapProp) + ".fa")
    del acc_kmers_ds_members_dict



if __name__ == "__main__":
    main()
