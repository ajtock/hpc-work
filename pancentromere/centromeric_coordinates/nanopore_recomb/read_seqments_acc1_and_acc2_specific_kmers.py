#!/usr/bin/env python3

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 25/08/2022

# Usage:
# conda activate python_3.9.6
# ./read_seqments_acc1_and_acc2_specific_kmers.py \
#  -r Col_ler_f1_pollen_500bp_minq99 \
#  -a1 Col-0.ragtag_scaffolds_centromeres \
#  -a2 Ler-0_110x.ragtag_scaffolds_centromeres \ 
#  -k 24 \
#  -mh 3 \
#  -hr 2
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
import subprocess

from Bio import SeqIO
from pathlib import Path
#from time import time, sleep
#import timeit

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
    parser.add_argument("-hr", "--hybReadNo", type=int, default="2",
                        help="The hybrid read number, defined according to the order it appears in the hybrid reads FASTA file. Default: 2")
    return parser

parser = create_parser().parse_args()
print(parser)

acc1_name = parser.acc1.split(".")[0].split("_")[0]
acc2_name = parser.acc2.split(".")[0].split("_")[0]

acc1_outdir_co = "segments/" + acc1_name + "/co"
acc2_outdir_co = "segments/" + acc2_name + "/co"
acc1_outdir_nco = "segments/" + acc1_name + "/nco"
acc2_outdir_nco = "segments/" + acc2_name + "/nco"

if not os.path.exists(acc1_outdir_co):
    os.makedirs(acc1_outdir_co)

if not os.path.exists(acc2_outdir_co):
    os.makedirs(acc2_outdir_co)

if not os.path.exists(acc1_outdir_nco):
    os.makedirs(acc1_outdir_nco)

if not os.path.exists(acc2_outdir_nco):
    os.makedirs(acc2_outdir_nco)


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


# Make a list containing the union of
# elements in an arbitrary number of lists
def union_lists(*lists):
    """
    Get the union of elements in an arbitrary number of lists.
    """
    return list(set.union(*map(set, lists)))


# Make a single list that contains all
# the elements of each sublist of a list
def flatten(lol):
    """
    Use list comprehension to flatten a list of lists (lol) into a
    single list composed of all the elements in each sublist.
    """
    return [item for sublist in lol for item in sublist]


# Within a read, find the 0-based start location of all occurrences
# of each accession-specific k-mer
def get_kmer_loc(kmers, read):
    """
    For a given read, get the within-read 0-based start locations of all k-mer matches.
    """
    kmer_loc_dict_list = []
    for h in range(len(kmers)):
        kmer_id = kmers[h].id
        kmer_acc = kmers[h].id.split("_", 1)[1]
        kmer_for = str(kmers[h].seq)
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
                                               "hit_start": kmer_matches[k],
                                               "hit_end": kmer_matches[k] + parser.kmerSize})
        else:
            print("k-mer already present in object")
    #
    return pd.DataFrame(kmer_loc_dict_list)


# Remove rows in acc1_kmer_loc_df whose k-mer match coordinate ranges
# overlap any of those in acc2_kmer_loc_df, and vice versa
# NOTE: requires two reciprocal function calls
def remove_overlaps(DF1, DF2):
    """
    Remove rows in DF1 whose k-mer match coordinate ranges
    overlap any of those in DF2.
    """
    DF1_no_overlaps =  pd.DataFrame()
    for h in range(len(DF1)):
        range_h = range(DF1["hit_start"][h], DF1["hit_end"][h])
        range_h_overlaps_counter = 0
        for j in range(len(DF2)):
            range_j = range(DF2["hit_start"][j], DF2["hit_end"][j])
            overlap_range_hj = range(max(range_h[0], range_j[0]),
                                     min(range_h[-1], range_j[-1])+1)
            if len(overlap_range_hj) > 0:
                range_h_overlaps_counter += 1
        #
        #print(range_h_overlaps_counter)
        if range_h_overlaps_counter == 0:
            DF1_no_overlaps = pd.concat(objs=[DF1_no_overlaps, DF1.iloc[h:h+1,:]],
                                        axis=0,
                                        ignore_index=True)
    #
    return DF1_no_overlaps


# For a given read, get accession-specific read segments
def get_read_segments(kmer_loc_df_sort):
    """
    For a given pandas DataFrame representing the sorted start locations of all
    accession-specific k-mers in a read, make a list in which each element is
    a subset of the DataFrame that represents an accession-specific read segment.
    This equates to extracting separate DataFrames where consecutive rows have the
    same value in the "acc" column.
    This function should be called twice:
    1. The first function call will exclude segments with < parser.minHits consecutive
       accession-specific k-mers, and the resulting list elements (retained segments)
       should be concatenated into a pandas DataFrame (sorted by k-mer hit_start location)
       separately (using concat_DF_list() on output), to be provided as the input to the second call.
    2. The second function call will be applied to the concatenated DataFrame consisting
       of filtered segments, in order to extract each segment DataFrame as a list element
       for segment length calculations. This second call will combine accession-specific
       segments into one extended segment where, following the first function call,
       there are no intervening short segments representing the other accession
       (excluded segments with < parser.minHits consecutive accession-specific k-mers).
    """
    segments_list = []
    segment = pd.DataFrame()
    for h in range(len(kmer_loc_df_sort)-1):
        if len(segment) == 0:
            segment = pd.concat(objs=[segment, kmer_loc_df_sort.iloc[h:h+1,:]],
                                axis=0,
                                ignore_index=True)
        if kmer_loc_df_sort.iloc[h+1].acc == kmer_loc_df_sort.iloc[h].acc:
            segment = pd.concat(objs=[segment, kmer_loc_df_sort.iloc[h+1:h+2,:]],
                                axis=0,
                                ignore_index=True)
        elif len(segment) >= parser.minHits:
            segments_list.append(segment)
            segment = pd.DataFrame()
        else:
            segment = pd.DataFrame()
    # Handle final segment in read
    if len(segment) >= parser.minHits:
        segments_list.append(segment)
    #
    return segments_list


# Concatenate a list of pandas DataFrames (corresponding to
# accession-specific read segments) into a single DataFrame,
# sorted by k-mer match start location in read
def concat_DF_list(DF_list):
    """
    Concatenate a list of pandas DataFrames into a single DataFrame,
    and sort by k-mer match start location in read.
    """
    concat_DF = pd.DataFrame()
    for h in range(len(DF_list)):
        concat_DF = pd.concat(objs=[concat_DF, DF_list[h]],
                              axis=0,
                              ignore_index=True)
    #
    concat_DF_sort = concat_DF.sort_values(by="hit_start",
                                           axis=0,
                                           ascending=True,
                                           kind="quicksort",
                                           ignore_index=True)
    #
    return concat_DF_sort


# Get the longest accession-specific read segment
def get_longest_read_segment(accspec_read_segments_list):
    """
    Get the longest accession-specific read segment from a list of
    pandas DataFrames containing the within-read start coordinates of
    accession-specific k-mer matches.
    """
    read_segments_length_dict = {}
    for h in range(len(accspec_read_segments_list)):
        read_segments_length_dict[h] = ( accspec_read_segments_list[h]["hit_start"][
            len(accspec_read_segments_list[h]) - 1
        ] + parser.kmerSize ) \
        - accspec_read_segments_list[h]["hit_start"][0]
    #
    read_segments_length_vals_max = max(read_segments_length_dict.values())
    for seg_key, seg_val in read_segments_length_dict.items():
        if seg_val == read_segments_length_vals_max:
            return accspec_read_segments_list[seg_key]


# Write accession-specific read segment to FASTA
# to supply to alignment software as input read
def write_fasta_from_SeqRecord(read, segment, outfile):
    """
    Extract read segment from read and write to FASTA.
    """
    record_segment = read[segment["hit_start"].iloc[0] :
                          segment["hit_end"].iloc[-1]]
    with open(outfile, "w") as output_handle:
        SeqIO.write(record_segment, output_handle, "fasta")


# Align accession-specific read segments to respective genome
def align_read_segment_wm_ont(segment_fasta, genome):
    """
    Align read in FASTA format to genome using winnowmap map-ont.
    """
    aln_cmd = ["winnowmap"] + \
              ["-W", "index/" + genome + "_repetitive_k15.txt"] + \
              ["-x", "map-ont"] + \
              ["-t", "32"] + \
              ["-p", "1.0"] + \
              ["-N", "10"] + \
              ["index/" + genome + ".fa"] + \
              [segment_fasta]
    outpaf = re.sub(".fasta", "_wm_ont.paf", segment_fasta)
    outerr = re.sub(".fasta", "_wm_ont.err", segment_fasta)
    with open(outpaf, "w") as outfile_handle, open(outerr, "w") as outerr_handle:
        subprocess.run(aln_cmd, stdout=outfile_handle, stderr=outerr_handle) 
    # Delete file(s) if unmapped
    if os.stat(outpaf).st_size == 0:
        subprocess.run(["rm", outpaf, outerr])
#    # Delete file(s) if unmapped
#    sam_view_cmd = ["samtools"] + \
#                   ["view", outsam]
#    sam_record = subprocess.Popen(sam_view_cmd, stdout=subprocess.PIPE)
#    sam_flag = int( subprocess.check_output(["cut", "-f2"], stdin=sam_record.stdout) )
#    sam_rm_cmd = ["rm"] + \
#                 [outsam, outerr]
#    if sam_flag == 4:
#        subprocess.run(sam_rm_cmd)


# Align accession-specific read segments to respective genome
def align_read_segment_mm_ont(segment_fasta, genome):
    """
    Align read in FASTA format to genome using minimap2 map-ont.
    """
    aln_cmd = ["minimap2"] + \
              ["-x", "map-ont"] + \
              ["-t", "32"] + \
              ["-p", "1.0"] + \
              ["-N", "10"] + \
              ["index/" + genome + ".fa"] + \
              [segment_fasta]
    outpaf = re.sub(".fasta", "_mm_ont.paf", segment_fasta)
    outerr = re.sub(".fasta", "_mm_ont.err", segment_fasta)
    with open(outpaf, "w") as outfile_handle, open(outerr, "w") as outerr_handle:
        subprocess.run(aln_cmd, stdout=outfile_handle, stderr=outerr_handle)
    # Delete file(s) if unmapped
    if os.stat(outpaf).st_size == 0:
        subprocess.run(["rm", outpaf, outerr])


# Align accession-specific read segments to respective genome
def align_read_segment_mm_sr(segment_fasta, genome):
    """
    Align read in FASTA format to genome using minimap2 sr.
    """
    aln_cmd = ["minimap2"] + \
              ["-x", "sr"] + \
              ["-t", "32"] + \
              ["-p", "1.0"] + \
              ["-N", "10"] + \
              ["index/" + genome + ".fa"] + \
              [segment_fasta]
    outpaf = re.sub(".fasta", "_mm_sr.paf", segment_fasta)
    outerr = re.sub(".fasta", "_mm_sr.err", segment_fasta)
    with open(outpaf, "w") as outfile_handle, open(outerr, "w") as outerr_handle:
        subprocess.run(aln_cmd, stdout=outfile_handle, stderr=outerr_handle)
    # Delete file(s) if unmapped
    if os.stat(outpaf).st_size == 0:
        subprocess.run(["rm", outpaf, outerr])


def main(read):
    """
    Execute functions on a given read to extract and map
    the longest accession-specific read segment.
    """
    #
    #
    # Get the within-read start locations of accession-specific k-mer matches
    acc1_kmer_loc_df_tmp = get_kmer_loc(kmers=acc1_kmers, read=str(read.seq)) 
    acc2_kmer_loc_df_tmp = get_kmer_loc(kmers=acc2_kmers, read=str(read.seq)) 
    #
    #
    # Remove rows in acc1_kmer_loc_df whose k-mer match coordinate ranges
    # overlap any of those in acc2_kmer_loc_df, and vice versa
    # NOTE: requires two reciprocal function calls
    acc1_kmer_loc_df = remove_overlaps(DF1=acc1_kmer_loc_df_tmp,
                                       DF2=acc2_kmer_loc_df_tmp)
    acc2_kmer_loc_df = remove_overlaps(DF1=acc2_kmer_loc_df_tmp,
                                       DF2=acc1_kmer_loc_df_tmp)
    del acc1_kmer_loc_df_tmp, acc2_kmer_loc_df_tmp
    #
    #
    # Stop main execution if acc1_kmer_loc_df or acc2_kmer_loc_df have < parser.minHits
    if len(acc1_kmer_loc_df) < parser.minHits or len(acc2_kmer_loc_df) < parser.minHits:
        print("Stopping for read " + read.id + " because\n" +
              "acc1_kmer_loc_df or acc2_kmer_loc_df has < " + str(parser.minHits) + " accession-specific k-mers")
        return
    #
    #
    # TODO: Use pyranges to remove rows in acc1_kmer_loc_df and acc2_kmer_loc_df
    # whose k-mer match coordinate ranges overlap between accessions
    # E.g.:
    # https://stackoverflow.com/questions/57032580/finding-overlaps-between-millions-of-ranges-intervals
    # https://stackoverflow.com/questions/49118347/pythonic-equivalent-to-reduce-in-r-granges-how-to-collapse-ranged-data
    # NOTE: Not needed given remove_overlaps() function using range() above
    #
    #
    # Concatenate and sort by k-mer match start location in read
    acc_kmer_loc_df = pd.concat(objs=[acc1_kmer_loc_df, acc2_kmer_loc_df],
                                axis=0,
                                ignore_index=True)
    acc_kmer_loc_df_sort_tmp = acc_kmer_loc_df.sort_values(by="hit_start",
                                                           axis=0,
                                                           ascending=True,
                                                           kind="quicksort",
                                                           ignore_index=True)
    #
    #
    # For a given read, get accession-specific read segments
    # get_read_segments function call 1:
    # The first function call will exclude segments with < parser.minHits consecutive
    # accession-specific k-mers, and the resulting list elements (retained segments)
    # should be concatenated into a pandas DataFrame (sorted by k-mer hit_start location)
    # separately (using concat_DF_list() on output), to be provided as the input to the second call
    acc_read_segments_list_tmp = get_read_segments(kmer_loc_df_sort=acc_kmer_loc_df_sort_tmp)
    del acc_kmer_loc_df_sort_tmp, acc_kmer_loc_df
    #
    #
    # Concatenate acc_read_segments_list_tmp into a single DataFrame,
    # sorted by k-mer hit_start location in read
    acc_kmer_loc_df_sort = concat_DF_list(DF_list=acc_read_segments_list_tmp)
    del acc_read_segments_list_tmp
    #
    #
    # For a given read, get accession-specific read segments
    # get_read_segments function call 2:
    # The second function call will be applied to the concatenated DataFrame consisting
    # of filtered segments, in order to extract each segment DataFrame as a list element
    # for segment length calculations. This second call will combine accession-specific
    # segments into one extended segment where, following the first function call,
    # there are no intervening short segments representing the other accession
    # (excluded segments with < parser.minHits consecutive accession-specific k-mers)
    acc_read_segments_list = get_read_segments(kmer_loc_df_sort=acc_kmer_loc_df_sort)
    #
    #
    # Get per-accession read segments lists
    acc1_read_segments_list = []
    acc2_read_segments_list = []
    for i in range(len(acc_read_segments_list)):
        if acc_read_segments_list[i]["acc"][0] == acc1_name and len(acc_read_segments_list[i]) >= parser.minHits:
            acc1_read_segments_list.append(acc_read_segments_list[i])
        elif acc_read_segments_list[i]["acc"][0] == acc2_name and len(acc_read_segments_list[i]) >= parser.minHits:
            acc2_read_segments_list.append(acc_read_segments_list[i])
    #
    #
    # Determine whether hybrid read represents a putative crossover or noncrossover
    if len(acc1_read_segments_list) > 1 or len(acc2_read_segments_list) > 1:
        acc1_outdir = acc1_outdir_nco
        acc2_outdir = acc2_outdir_nco
    elif len(acc1_read_segments_list) == 1 and len(acc2_read_segments_list) == 1:
        acc1_outdir = acc1_outdir_co
        acc2_outdir = acc2_outdir_co
    #
    #
    # Get the longest accession-specific read segment for each accession 
    acc1_longest_read_segment = get_longest_read_segment(accspec_read_segments_list=acc1_read_segments_list)
    acc2_longest_read_segment = get_longest_read_segment(accspec_read_segments_list=acc2_read_segments_list)
    #
    #
    # Define output FASTA file names for writing read segments 
    acc1_outfile = acc1_outdir + "/" + read.id + "__" + acc1_longest_read_segment["acc"][0] + ".fasta"  
    acc2_outfile = acc2_outdir + "/" + read.id + "__" + acc2_longest_read_segment["acc"][0] + ".fasta"  
    #
    #
    # Write accession-specific read segment to FASTA
    # to supply to alignment software as input read
    write_fasta_from_SeqRecord(read=read,
                               segment=acc1_longest_read_segment,
                               outfile=acc1_outfile)
    write_fasta_from_SeqRecord(read=read,
                               segment=acc2_longest_read_segment,
                               outfile=acc2_outfile)
    #
    #
    # Align accession-specific read segments to respective genome
    align_read_segment_wm_ont(segment_fasta=acc1_outfile,
                              genome=re.sub("_centromeres", "", parser.acc1))
    align_read_segment_wm_ont(segment_fasta=acc2_outfile,
                              genome=re.sub("_centromeres", "", parser.acc2))
    #
    align_read_segment_mm_ont(segment_fasta=acc1_outfile,
                              genome=re.sub("_centromeres", "", parser.acc1))
    align_read_segment_mm_ont(segment_fasta=acc2_outfile,
                              genome=re.sub("_centromeres", "", parser.acc2))
    #
    align_read_segment_mm_sr(segment_fasta=acc1_outfile,
                             genome=re.sub("_centromeres", "", parser.acc1))
    align_read_segment_mm_sr(segment_fasta=acc2_outfile,
                             genome=re.sub("_centromeres", "", parser.acc2))



## TODO: make Slurm script to parallelise main function over reads


if __name__ == "__main__":
    print("Hybrid read number: " + str(parser.hybReadNo))
    main(read=reads[parser.hybReadNo])
