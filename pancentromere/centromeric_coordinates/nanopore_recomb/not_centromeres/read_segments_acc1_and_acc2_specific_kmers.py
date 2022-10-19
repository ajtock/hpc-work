#!/usr/bin/env python3

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 14/10/2022

# Usage (sbatch read_segments_acc1_and_acc2_specific_kmers_py_icelake_slurm):
# conda activate python_3.9.6
# ./read_segments_acc1_and_acc2_specific_kmers.py \
#  -r Col_ler_f1_pollen_500bp_minq99 \
#  -a1 Col-0.ragtag_scaffolds_not_centromeres \
#  -a2 Ler-0_110x.ragtag_scaffolds_not_centromeres \ 
#  -k 24 \
#  -op 0.9 \
#  -mh 10 \
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
import glob
import argparse
#import pickle
import re
import numpy as np
import pandas as pd
import subprocess
import screed

from Bio import SeqIO
from pathlib import Path
from numba import jit
from time import time, sleep
#import timeit

# ==== Capture user input as command-line arguments
# https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
def create_parser():
    parser = argparse.ArgumentParser(description="Fasta filename variables.")
    #### Define command-line arguments
    parser.add_argument("-r", "--readsPrefix", type=str, default="Col_ler_f1_pollen_500bp_minq99",
                        help="The prefix of the FASTQ file name. Default: Col_ler_f1_pollen_500bp_minq99")
    parser.add_argument("-a1", "--acc1", type=str, default="Col-0.ragtag_scaffolds_not_centromeres",
                        help="The prefix of the first accession's sequences. Default: Col-0.ragtag_scaffolds_not_centromeres")
    parser.add_argument("-a2", "--acc2", type=str, default="Ler-0_110x.ragtag_scaffolds_not_centromeres",
                        help="The prefix of the second accession's sequences. Default: Ler-0_110x.ragtag_scaffolds_not_centromeres")
    parser.add_argument("-k", "--kmerSize", type=int, default="24",
                        help="The size of the k-mers to be found and counted in the FASTA file. Default: 24")
    parser.add_argument("-op", "--overlapProp", type=float, default="0.9",
                        help="The minimum proportion of an aligned k-mer's length that must overlap a genomic window for the aligned k-mer to be kept during downsampling of accession-specific k-mers. Default: 0.9")
    parser.add_argument("-mh", "--minHits", type=int, default="10",
                        help="The minimum number of accession-specific k-mers found in a read. Default: 10")
    parser.add_argument("-hr", "--hybReadNo", type=int, default="0",
                        help="The hybrid read number, defined according to the order it appears in the hybrid reads FASTA file. Default: 0")
    return parser

parser = create_parser().parse_args()
print(parser)


acc1_name = parser.acc1.split(".")[0].split("_")[0]
acc2_name = parser.acc2.split(".")[0].split("_")[0]

acc1_outdir_co = "segments/" + acc1_name + "/co"
acc2_outdir_co = "segments/" + acc2_name + "/co"
acc1_outdir_nco = "segments/" + acc1_name + "/nco"
acc2_outdir_nco = "segments/" + acc2_name + "/nco"

#if not os.path.exists(acc1_outdir_co):
#    os.makedirs(acc1_outdir_co)
#
#if not os.path.exists(acc2_outdir_co):
#    os.makedirs(acc2_outdir_co)
#
#if not os.path.exists(acc1_outdir_nco):
#    os.makedirs(acc1_outdir_nco)
#
#if not os.path.exists(acc2_outdir_nco):
#    os.makedirs(acc2_outdir_nco)


# Path to hybrid reads
input_fa = "fasta/" + parser.readsPrefix + \
    "_match_" + parser.acc1 + \
    "_specific_k" + str(parser.kmerSize) + \
    "_downsampled_op" + str(parser.overlapProp) + \
    "_hits" + str(parser.minHits) + \
    "_match_" + parser.acc2 + \
    "_specific_k" + str(parser.kmerSize) + \
    "_downsampled_op" + str(parser.overlapProp) + \
    "_hits" + str(parser.minHits) + \
    ".fa"
# File exists sanity check
Path(input_fa).resolve(strict=True)

# Path to acc1-specific k-mers
acc1_fa = "fasta/" + \
    parser.acc1 + "_specific_k" + \
    str(parser.kmerSize) + "_downsampled_op" + \
    str(parser.overlapProp) + "_noheaders.fa"
# File exists sanity check
Path(acc1_fa).resolve(strict=True)

# Path to acc2-specific k-mers
acc2_fa = "fasta/" + \
    parser.acc2 + "_specific_k" + \
    str(parser.kmerSize) + "_downsampled_op" + \
    str(parser.overlapProp) + "_noheaders.fa"
# File exists sanity check
Path(acc2_fa).resolve(strict=True)


# Parse reads and get read corresponding to index parser.hybReadNo
#reads = list(SeqIO.parse(input_fa, "fasta"))
#reads_iter = SeqIO.parse(input_fa, "fasta")
# Dictionary approach assumes use of Python >= 3.7, because in previous
# Python versions dictionaries were inherently unordered
reads_dict = SeqIO.index(input_fa, "fasta")
# Commented-out approach to indexed read extraction from dictionary takes longer, although simpler
#read = list(reads_dict.values())[parser.hybReadNo]
read = next(v for i, v in enumerate(reads_dict.values()) if i == parser.hybReadNo)
#read_test = [v for i, v in enumerate(reads_dict.values()) if i == parser.hybReadNo or i == 10749]
del reads_dict
print("Hybrid read number: " + str(parser.hybReadNo))


## Parse acc1-specific k-mers as FastaIterator
#acc1_kmers_iter = SeqIO.parse(acc1_fa, "fasta")
##acc1_kmers = list(SeqIO.parse(acc1_fa, "fasta"))
#
## Parse acc2-specific k-mers as FastaIterator
#acc2_kmers_iter = SeqIO.parse(acc2_fa, "fasta")
##acc2_kmers = list(SeqIO.parse(acc2_fa, "fasta"))


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


## From the list acc1_kmers or acc2_kmers, define a new list of kmers
## (of class 'Bio.SeqRecord.SeqRecord') that is the subset found in read
#def get_kmer_subset(kmers, read):
#    """
#    From the list acc1_kmers or acc2_kmers, define a new list of kmers
#    (of class 'Bio.SeqRecord.SeqRecord') that is the subset found in read.
#    """
#seqrecord_list2 = [record for record in acc1_kmers_iter if re.search(str(record.seq), str(read.seq)) or re.search(screed.rc(str(record.seq)), str(read.seq))]
#seqrecord_list = []
#
#SeqIO.write(seqrecord_iterator, acc1_subset_fa, "fasta")
#
#for record in seqrecord_iterator:
#    seqrecord_list.append(record)
#
#    kmer_seqrecord_list = []
#    for h in range(len(kmers)):
#        kmer_for = str(kmers[h].seq)
#        kmer_rev = screed.rc(kmer_for)
#        kmer_for_match = re.search(kmer_for, read)
#        kmer_rev_match = re.search(kmer_rev, read)
#        if kmer_for_match or kmer_rev_match:
#            kmer_seqrecord_list.append(kmers[h])
#    #
#    return kmer_seqrecord_list
#
#acc1_kmers_read_subset = get_kmer_subset(kmers=acc1_kmers, read=str(read.seq))

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
    with open(kmers_fa_noheaders, "r") as kmers_fa_noheaders_handle:
        lines = kmers_fa_noheaders_handle.read().splitlines()
        lines_rc = list(map(screed.rc, lines))
        lines_tuple = list(zip(lines, lines_rc))
        kmers_list = list(map(min, lines_tuple))
    # Dictionary approach to duplicate k-mer removal
    # (1 representative k-mer retained where duplicates exist)
    # will maintain the insertion order of the k-mers in the list
    # (assumes use of Python >= 3.7; in previous Python versions,
    # dictionaries were inherently unordered)
    return list(dict.fromkeys(kmers_list))



# Within a read, find the 0-based start location of all occurrences
# of each accession-specific k-mer
def get_kmer_loc_map(kmers_fa_noheaders, read_seq):
kmers_fa_noheaders=acc1_fa
kmers_fa=re.sub("_noheaders", "", kmers_fa_noheaders)
read_seq=str(read.seq)
"""
For a given read, get the within-read 0-based start locations of all k-mer matches.
"""
kmers_fa_noheaders_handle = open(kmers_fa_noheaders, "r")
kmers = kmers_fa_noheaders_handle.read().splitlines()

def kmer_in_read_search(kmer_x):
    if re.search(kmer_x, read_seq) or re.search(screed.rc(kmer_x), read_seq):
        return True
    else:
        return False

kmers_in_read = [kmer for kmer in kmers if
                 re.search(kmer, read_seq) or
                 re.search(screed.rc(kmer), read_seq)]
kmers_in_read_v2 = list(map(lambda kmer: re.search(kmer, read_seq), kmers))
kmers_in_read_v3 = filter(kmer_in_read_search, kmers)
kmers_in_read_v3_list = []
for x in kmers_in_read_v3:
    kmers_in_read_v3_list.append(x)
 

tmp = list(map re.search(lines[0], read_seq)

ages = [5, 12, 17, 18, 24, 32]

def myFunc(x):
  if x < 18:
    return False
  else:
    return True

adults = filter(myFunc, ages)

adults_list = []
for x in adults:
  adults_list.append(x)


kmers_iter = SeqIO.parse(kmers_fa, "fasta")
kmers = [record for record in kmers_iter if
         re.search(str(record.seq), read_seq) or
         re.search(screed.rc(str(record.seq)), read_seq)]
del kmers_iter
kmer_loc_dict_list = []
for h in range(len(kmers)):
    kmer_id = kmers[h].id
    kmer_acc = kmers[h].id.split("_", 1)[1]
    kmer_for = str(kmers[h].seq)
    #kmer_rev = kmer_for.translate(comp_tab)[::-1] 
    kmer_rev = screed.rc(kmer_for)
    kmer_for_matches = [match.start() for match in re.finditer(kmer_for, read_seq)]
    kmer_rev_matches = [match.start() for match in re.finditer(kmer_rev, read_seq)]
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


# Within a read, find the 0-based start location of all occurrences
# of each accession-specific k-mer
def get_kmer_loc(kmers_fa, read_seq):
    """
    For a given read, get the within-read 0-based start locations of all k-mer matches.
    """
    kmers_iter = SeqIO.parse(kmers_fa, "fasta")
    kmers = [record for record in kmers_iter if
             re.search(str(record.seq), read_seq) or
             re.search(screed.rc(str(record.seq)), read_seq)]
    del kmers_iter
    kmer_loc_dict_list = []
    for h in range(len(kmers)):
        kmer_id = kmers[h].id
        kmer_acc = kmers[h].id.split("_", 1)[1]
        kmer_for = str(kmers[h].seq)
        #kmer_rev = kmer_for.translate(comp_tab)[::-1] 
        kmer_rev = screed.rc(kmer_for)
        kmer_for_matches = [match.start() for match in re.finditer(kmer_for, read_seq)]
        kmer_rev_matches = [match.start() for match in re.finditer(kmer_rev, read_seq)]
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


# Within a read, find the 0-based start location of all occurrences
# of each accession-specific k-mer
def get_kmer_loc_ori(kmers, read):
    """
    For a given read, get the within-read 0-based start locations of all k-mer matches.
    """
    kmer_loc_dict_list = []
    for h in range(len(kmers)):
        kmer_id = kmers[h].id
        kmer_acc = kmers[h].id.split("_", 1)[1]
        kmer_for = str(kmers[h].seq)
        #kmer_rev = kmer_for.translate(comp_tab)[::-1] 
        kmer_rev = screed.rc(kmer_for)
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
              ["-t", "1"] + \
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
              ["-t", "1"] + \
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
              ["-t", "1"] + \
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


# Delete accession-specific segment alignment file if the equivalent
# file for the other accession doesn't exist (indicating an unmapped segment)
def delete_alignment(alignment_prefix1, alignment_prefix2):
    """
    Delete accession-specific segment alignment file if the equivalent
    file for the other accession doesn't exist (indicating an unmapped segment).
    """
    if not glob.glob(alignment_prefix2 + "*"):
        file_list = glob.glob(alignment_prefix1 + "*", recursive=True)
        for file_path in file_list:
            try:
                os.remove(file_path)
            except OSError:
                print("Error while deleting file " + file_path)



def main():
    """
    Execute functions on a given read to extract and map
    the longest accession-specific read segment.
    """
    
    
    # Get the within-read start locations of accession-specific k-mer matches
    acc1_kmer_loc_df_tmp = get_kmer_loc(kmers_fa=acc1_fa, read_seq=str(read.seq)) 
    acc2_kmer_loc_df_tmp = get_kmer_loc(kmers_fa=acc2_fa, read_seq=str(read.seq)) 
    
    
    # Remove rows in acc1_kmer_loc_df whose k-mer match coordinate ranges
    # overlap any of those in acc2_kmer_loc_df, and vice versa
    # NOTE: requires two reciprocal function calls
    acc1_kmer_loc_df = remove_overlaps(DF1=acc1_kmer_loc_df_tmp,
                                       DF2=acc2_kmer_loc_df_tmp)
    acc2_kmer_loc_df = remove_overlaps(DF1=acc2_kmer_loc_df_tmp,
                                       DF2=acc1_kmer_loc_df_tmp)
    del acc1_kmer_loc_df_tmp, acc2_kmer_loc_df_tmp
    
    
    # Stop main execution if acc1_kmer_loc_df or acc2_kmer_loc_df have < parser.minHits
    if len(acc1_kmer_loc_df) < parser.minHits or len(acc2_kmer_loc_df) < parser.minHits or not "acc1_kmer_loc_df" in locals() or not "acc1_kmer_loc_df" in locals():
        print("Stopping for read " + str(parser.hybReadNo) + ": " + read.id + " because\n" +
              "acc1_kmer_loc_df or acc2_kmer_loc_df has < " + str(parser.minHits) + " accession-specific k-mers")
        return
    
    
    # TODO: Use pyranges to remove rows in acc1_kmer_loc_df and acc2_kmer_loc_df
    # whose k-mer match coordinate ranges overlap between accessions
    # E.g.:
    # https://stackoverflow.com/questions/57032580/finding-overlaps-between-millions-of-ranges-intervals
    # https://stackoverflow.com/questions/49118347/pythonic-equivalent-to-reduce-in-r-granges-how-to-collapse-ranged-data
    # NOTE: Not needed given remove_overlaps() function using range() above
    
    
    # Concatenate and sort by k-mer match start location in read
    acc_kmer_loc_df = pd.concat(objs=[acc1_kmer_loc_df, acc2_kmer_loc_df],
                                axis=0,
                                ignore_index=True)
    acc_kmer_loc_df_sort_tmp = acc_kmer_loc_df.sort_values(by="hit_start",
                                                           axis=0,
                                                           ascending=True,
                                                           kind="quicksort",
                                                           ignore_index=True)
    
    
    # For a given read, get accession-specific read segments
    # get_read_segments function call 1:
    # The first function call will exclude segments with < parser.minHits consecutive
    # accession-specific k-mers, and the resulting list elements (retained segments)
    # should be concatenated into a pandas DataFrame (sorted by k-mer hit_start location)
    # separately (using concat_DF_list() on output), to be provided as the input to the second call
    acc_read_segments_list_tmp = get_read_segments(kmer_loc_df_sort=acc_kmer_loc_df_sort_tmp)
    del acc_kmer_loc_df_sort_tmp, acc_kmer_loc_df
    
    
    # Concatenate acc_read_segments_list_tmp into a single DataFrame,
    # sorted by k-mer hit_start location in read
    acc_kmer_loc_df_sort = concat_DF_list(DF_list=acc_read_segments_list_tmp)
    del acc_read_segments_list_tmp
    
    
    # For a given read, get accession-specific read segments
    # get_read_segments function call 2:
    # The second function call will be applied to the concatenated DataFrame consisting
    # of filtered segments, in order to extract each segment DataFrame as a list element
    # for segment length calculations. This second call will combine accession-specific
    # segments into one extended segment where, following the first function call,
    # there are no intervening short segments representing the other accession
    # (excluded segments with < parser.minHits consecutive accession-specific k-mers)
    acc_read_segments_list = get_read_segments(kmer_loc_df_sort=acc_kmer_loc_df_sort)
    
    
    # Stop main execution if acc_read_segments_list has < 2 elements
    # (accession-specific segments)
    if len(acc_read_segments_list) < 2:
        print("Stopping for read " + str(parser.hybReadNo) + ": " + read.id + " because\n" +
              "acc_read_segments_list has < 2 elements (accession-specific segments)")
        return
    
    
    # Get per-accession read segments lists
    acc1_read_segments_list = []
    acc2_read_segments_list = []
    for i in range(len(acc_read_segments_list)):
        if acc_read_segments_list[i]["acc"][0] == acc1_name and len(acc_read_segments_list[i]) >= parser.minHits:
            acc1_read_segments_list.append(acc_read_segments_list[i])
        elif acc_read_segments_list[i]["acc"][0] == acc2_name and len(acc_read_segments_list[i]) >= parser.minHits:
            acc2_read_segments_list.append(acc_read_segments_list[i])
    
    
    # Stop main execution if acc1_read_segments_list or acc2_read_segments_list
    # is empty (without accession-specific segments)
    if not acc1_read_segments_list or not acc2_read_segments_list:
        print("Stopping for read " + str(parser.hybReadNo) + ": " + read.id + " because\n" +
              "acc1_read_segments_list or acc2_read_segments_list is empty (without accession-specific segments)")
        return
    
    
    # Determine whether hybrid read represents a putative crossover or noncrossover
    if len(acc1_read_segments_list) > 1 or len(acc2_read_segments_list) > 1:
        acc1_outdir = acc1_outdir_nco
        acc2_outdir = acc2_outdir_nco
    elif len(acc1_read_segments_list) == 1 and len(acc2_read_segments_list) == 1:
        acc1_outdir = acc1_outdir_co
        acc2_outdir = acc2_outdir_co
    
    
    # Stop main execution if acc1_outdir or acc2_outdir is not defined
    if not "acc1_outdir" in locals() or not "acc2_outdir" in locals():
        print("Stopping for read " + str(parser.hybReadNo) + ": " + read.id + " because\n" +
              "acc1_outdir or acc2_outdir is not defined")
        return
    
    
    # Get the longest accession-specific read segment for each accession 
    acc1_longest_read_segment = get_longest_read_segment(accspec_read_segments_list=acc1_read_segments_list)
    acc2_longest_read_segment = get_longest_read_segment(accspec_read_segments_list=acc2_read_segments_list)
    
    
    # Stop main execution if acc1_longest_read_segment or acc2_longest_read_segment
    # is not defined
    if not "acc1_longest_read_segment" in locals() or not "acc2_longest_read_segment" in locals():
        print("Stopping for read " + str(parser.hybReadNo) + ": " + read.id + " because\n" +
              "acc1_longest_read_segment or acc2_longest_read_segment is not defined")
        return
    
    
    # Define output FASTA file names for writing read segments 
    acc1_outfile = acc1_outdir + "/" + read.id + "__" + acc1_longest_read_segment["acc"][0] + ".fasta"
    acc2_outfile = acc2_outdir + "/" + read.id + "__" + acc2_longest_read_segment["acc"][0] + ".fasta"
    
    
    # Write accession-specific read segment to FASTA
    # to supply to alignment software as input read
    write_fasta_from_SeqRecord(read=read,
                               segment=acc1_longest_read_segment,
                               outfile=acc1_outfile)
    write_fasta_from_SeqRecord(read=read,
                               segment=acc2_longest_read_segment,
                               outfile=acc2_outfile)
    
    
    # Align accession-specific read segments to respective genome
    align_read_segment_wm_ont(segment_fasta=acc1_outfile,
                              genome=re.sub(r"(_scaffolds)_.+", r"\1", parser.acc1))
    align_read_segment_wm_ont(segment_fasta=acc2_outfile,
                              genome=re.sub(r"(_scaffolds)_.+", r"\1", parser.acc2))
    
    align_read_segment_mm_ont(segment_fasta=acc1_outfile,
                              genome=re.sub(r"(_scaffolds)_.+", r"\1", parser.acc1))
    align_read_segment_mm_ont(segment_fasta=acc2_outfile,
                              genome=re.sub(r"(_scaffolds)_.+", r"\1", parser.acc2))
    
    align_read_segment_mm_sr(segment_fasta=acc1_outfile,
                             genome=re.sub(r"(_scaffolds)_.+", r"\1", parser.acc1))
    align_read_segment_mm_sr(segment_fasta=acc2_outfile,
                             genome=re.sub(r"(_scaffolds)_.+", r"\1", parser.acc2))
    
    
    # Delete accession-specific segment alignment file if the equivalent
    # file for the other accession doesn't exist (indicating an unmapped segment)
    acc1_alignment_prefix = acc1_outdir + "/" + read.id + "__" + acc1_name + "_"
    acc2_alignment_prefix = acc2_outdir + "/" + read.id + "__" + acc2_name + "_"
    
    delete_alignment(alignment_prefix1=acc1_alignment_prefix,
                     alignment_prefix2=acc2_alignment_prefix)
    delete_alignment(alignment_prefix1=acc2_alignment_prefix,
                     alignment_prefix2=acc1_alignment_prefix)



if __name__ == "__main__":
    main()
