#!/usr/bin/env python3

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 18/11/2022

# Extract putative recombination events detected in Col-0/Ler-0
# hybrid ONT reads identified on the basis that they contain >=
# a threshold number of accession-specific, chromosome-specific k-mers

# Usage:
# conda activate python_3.9.6
# ./get_chr_specific_segment_pairs_alnToDifferent.py \
#  -r ColLerF1pollen_1000bp_minq90 \
#  -a1 Col-0.ragtag_scaffolds \
#  -a2 Ler-0_110x.ragtag_scaffolds \
#  -k 24 \
#  -op 0.9 \
#  -mh 11 \
#  -md 1e6 \
#  -aq 0.90 \
#  -nq 0.90 \
#  -rt co \
#  -reg centromere \
#  -c 'Chr1'
# conda deactivate


# ==== Import libraries
import os
import argparse
import re
import gc
import pandas as pd
import numpy as np
import math
import subprocess
import glob

from pathlib import Path
from Bio import SeqIO


# ==== Capture user input as command-line arguments
# https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
def create_parser():
    parser = argparse.ArgumentParser(description="Fasta filename variables.")
    #### Define command-line arguments
    parser.add_argument("-r", "--readsPrefix", type=str, default="ColLerF1pollen_1000bp_minq90",
                        help="The prefix of the FASTQ file name. Default: ColLerF1pollen_1000bp_minq90")
    parser.add_argument("-a1", "--acc1", type=str, default="Col-0.ragtag_scaffolds",
                        help="The prefix of the first accession's sequences. Default: Col-0.ragtag_scaffolds")
    parser.add_argument("-a2", "--acc2", type=str, default="Ler-0_110x.ragtag_scaffolds",
                        help="The prefix of the second accession's sequences. Default: Ler-0_110x.ragtag_scaffolds")
    parser.add_argument("-k", "--kmerSize", type=int, default="24",
                        help="The size of the k-mers to be found and counted in the FASTA file. Default: 24")
    parser.add_argument("-op", "--overlapProp", type=float, default="0.9",
                        help="The minimum proportion of an aligned k-mer's length that must overlap a genomic window for the aligned k-mer to be kept during downsampling of accession-specific k-mers. Default: 0.9")
    parser.add_argument("-mh", "--minHits", type=int, default="11",
                        help="The minimum number of accession-specific k-mers found in a read. Default: 11")
    parser.add_argument("-mmq", "--minMAPQ", type=int, default="2",
                        help="The minimum alignment MAPQ score allowed across all aligners. Default: 2")
    parser.add_argument("-aq", "--alenTOqlen", type=float, default="0.90",
                        help="The minimum ratio of the read segment alignment length to the read segment length. Default: 0.90")
    parser.add_argument("-nq", "--nmatchTOqlen", type=float, default="0.90",
                        help="The minimum ratio of the read segment alignment number of matching bases to the read segment length. Default: 0.90")
    parser.add_argument("-rt", "--recombType", type=str, default="co",
                        help="The type/pattern of the recombination event identified based on the sequence of accession-specific read segments. Default: co")
    parser.add_argument("-reg", "--region", type=str, default="centromere",
                        help="The chromosome for which accession-specific, chromosome-specific read segments have been extracted and aligned. Default: centromere")
    parser.add_argument("-c", "--chrom", type=str, default="Chr1",
                        help="The chromosome for which accession-specific, chromosome-specific read segments have been extracted and aligned. Default: Chr1")
    return parser

parser = create_parser().parse_args()
print(parser)

chrom = parser.chrom.split(",")

#maxDist = int(parser.maxDist)
#
#if math.floor(math.log10(maxDist)) + 1 < 4:
#    maxDistName = str(int(maxDist)) + "bp"
#elif math.floor(math.log10(maxDist)) + 1 >= 4 and \
#     math.floor(math.log10(maxDist)) + 1 <= 6:
#    maxDistName = str(int(maxDist / 1e3)) + "kb"
#elif math.floor(math.log10(maxDist)) + 1 >= 7:
#    maxDistName = str(int(maxDist / 1e6)) + "Mb"

outdir = parser.region + "/segment_pairs/" + parser.recombType

if not os.path.exists(outdir):
    os.makedirs(outdir)

# Accession names
acc1_name = parser.acc1.split(".")[0].split("_")[0]
acc2_name = parser.acc2.split(".")[0].split("_")[0]

# Directories containing read segment alignment files
acc1_indir_list = [parser.region + "/" +  x + "/segments/" + acc1_name + "/" + parser.recombType for x in chrom]
acc2_indir_list = [parser.region + "/" +  x + "/segments/" + acc2_name + "/" + parser.recombType for x in chrom]


# CEN coordinates
CEN = pd.read_csv("/rds/project/rds-O5Ty9yVfQKg/pancentromere/centromeric_coordinates/" + \
                  "centromere_manual_EDTA4_fa.csv")
CEN["fasta.name"].replace(to_replace="\.fa", value="", regex=True, inplace=True)


# Genomic definitions
# acc1
acc1_fai = pd.read_csv("index/" + parser.acc1 + ".fa.fai",
                       sep="\t", header=None)
acc1_fai = acc1_fai.loc[acc1_fai[0].isin(chrom)]
#acc1_chrs = acc1_fai[0].to_string(index=False).split()
acc1_chrs = list(acc1_fai[0])
acc1_chrLens = list(acc1_fai[1])

acc1_CEN = CEN[CEN["fasta.name"] == parser.acc1]
acc1_CEN = acc1_CEN.loc[:, acc1_CEN.columns.isin(["chr", "start", "end"])]
acc1_CEN_new = pd.DataFrame()
for i in range(0, len(acc1_chrs)):
    acc1_CEN_chr = acc1_CEN[acc1_CEN["chr"] == acc1_chrs[i]]
    if acc1_CEN_chr.shape[0] > 1:
        acc1_CEN_chr = pd.DataFrame({ "chr":   [ acc1_CEN_chr["chr"].iloc[0] ],
                                      "start": [ acc1_CEN_chr["start"].iloc[0] ],
                                      "end":   [ acc1_CEN_chr["end"].iloc[-1] ] })
    acc1_CEN_new = pd.concat(objs=[acc1_CEN_new, acc1_CEN_chr],
                             axis=0,
                             ignore_index=True)

acc1_CEN = acc1_CEN_new
acc1_CENstart = acc1_CEN["start"]
acc1_CENend = acc1_CEN["end"]
acc1_chrs = [acc1_name + "_" + x for x in acc1_chrs]
acc1_CENsize = list(acc1_CEN["end"] - acc1_CEN["start"] + 1)


# acc2
acc2_fai = pd.read_csv("index/" + parser.acc2 + ".fa.fai",
                       sep="\t", header=None)
acc2_fai = acc2_fai.loc[acc2_fai[0].isin(chrom)]
#acc2_chrs = acc2_fai[0].to_string(index=False).split()
acc2_chrs = list(acc2_fai[0])
acc2_chrLens = list(acc2_fai[1])

acc2_CEN = CEN[CEN["fasta.name"] == parser.acc2]
acc2_CEN = acc2_CEN.loc[:, acc2_CEN.columns.isin(["chr", "start", "end"])]
acc2_CEN_new = pd.DataFrame()
for i in range(0, len(acc2_chrs)):
    acc2_CEN_chr = acc2_CEN[acc2_CEN["chr"] == acc2_chrs[i]]
    if acc2_CEN_chr.shape[0] > 1:
        acc2_CEN_chr = pd.DataFrame({ "chr":   [ acc2_CEN_chr["chr"].iloc[0] ],
                                      "start": [ acc2_CEN_chr["start"].iloc[0] ],
                                      "end":   [ acc2_CEN_chr["end"].iloc[-1] ] })
    acc2_CEN_new = pd.concat(objs=[acc2_CEN_new, acc2_CEN_chr],
                             axis=0,
                             ignore_index=True)

acc2_CEN = acc2_CEN_new
acc2_CENstart = acc2_CEN["start"]
acc2_CENend = acc2_CEN["end"]
acc2_chrs = [acc2_name + "_" + x for x in acc2_chrs]
acc2_CENsize = list(acc2_CEN["end"] - acc2_CEN["start"] + 1)


CENsize_nparray = np.array([acc1_CENsize,
                            acc2_CENsize])
min_CENsize_list = list(CENsize_nparray.min(axis=0))
max_CENsize_list = list(CENsize_nparray.max(axis=0))


# Concatenate all alignment files for the given chromosome,
# accession and aligner
# NOTE: for some reason the "find ..." approach below doesn't work
# from within python, so need to create and run an equivalent bash script:
#cat_cmd = ["find"] + \
#          [indir + "/"] + \
#          ["-mindepth", "1"] + \
#          ["-maxdepth", "1"] + \
#          ["-type", "f"] + \
#          ["-name", "*" + acc_name + suffix] + \
#          ["-exec cat {} + >> cat.paf"]
#subprocess.run(cat_cmd)
# For details on -exec, see https://stackoverflow.com/questions/2961673/find-missing-argument-to-exec
def cat_pafs(indir, acc_name, suffix):
    ##indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc1_indir_list[1]
    #indir=acc1_indir_list[1]
    #acc_name=acc1_name
    #suffix="_alnTo_" + parser.acc1 + "_Chr_mm_ont.paf"
    """
    Concatenate all alignment files for the given chromosome,
    accession and aligner.
    """
    outdir_paf = indir + "/cat_paf"
    if not os.path.exists(outdir_paf):
        os.makedirs(outdir_paf)
    #
    out_paf = outdir_paf + "/all_segments__" + acc_name + suffix
    if os.path.exists(out_paf):
        print("Concatenated alignment file " + out_paf + " already exists!")
        return
    else:
        cat_pafs_script = indir + "/find_cat_pafs.sh"
        with open(cat_pafs_script, "w") as cat_pafs_script_handle:
            cat_pafs_script_handle.write("#!/bin/bash\n\n" + \
                                         "find " + indir + "/ \\\n" + \
                                         "  -mindepth 1 \\\n" +\
                                         "  -maxdepth 1 \\\n" +\
                                         "  -type f \\\n" + \
                                         "  -name '*" + acc_name + suffix + "' \\\n" + \
                                         "  -exec cat {} + >> " + out_paf)
        #
        subprocess.run(["bash", cat_pafs_script])
        return


for x in range(0, len(acc1_indir_list)):
    print(acc1_indir_list[x])
    cat_pafs(indir=acc1_indir_list[x],
             acc_name=acc1_name,
             suffix="_alnTo_" + parser.acc1 + "_Chr_wm_ont.paf")
    cat_pafs(indir=acc1_indir_list[x],
             acc_name=acc1_name,
             suffix="_alnTo_" + parser.acc1 + "_Chr_mm_ont.paf")
    cat_pafs(indir=acc1_indir_list[x],
             acc_name=acc1_name,
             suffix="_alnTo_" + parser.acc1 + "_Chr_mm_sr.paf")
    cat_pafs(indir=acc1_indir_list[x],
             acc_name=acc1_name,
             suffix="_alnTo_" + parser.acc1 + "_Chr_bt2.paf")

for x in range(0, len(acc2_indir_list)):
    print(acc2_indir_list[x])
    cat_pafs(indir=acc2_indir_list[x],
             acc_name=acc2_name,
             suffix="_alnTo_" + parser.acc2 + "_Chr_wm_ont.paf")
    cat_pafs(indir=acc2_indir_list[x],
             acc_name=acc2_name,
             suffix="_alnTo_" + parser.acc2 + "_Chr_mm_ont.paf")
    cat_pafs(indir=acc2_indir_list[x],
             acc_name=acc2_name,
             suffix="_alnTo_" + parser.acc2 + "_Chr_mm_sr.paf")
    cat_pafs(indir=acc2_indir_list[x],
             acc_name=acc2_name,
             suffix="_alnTo_" + parser.acc2 + "_Chr_bt2.paf")


# Load concatenated read segment alignment file as a DataFrame
def load_cat_paf(indir, acc_name, suffix, aligner):
    ##indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc1_indir_list[1]
    #indir=acc1_indir_list[1]
    #acc_name=acc1_name
    #suffix="_alnTo_" + parser.acc1 + "_Chr_mm_ont.paf"
    #aligner="mm"
    """
    Load concatenated read segment alignment file as a DataFrame.
    """
    #cat_pafs(indir=indir, acc_name=acc_name, suffix=suffix)
    indir_paf = indir + "/cat_paf"
    in_paf = indir_paf + "/all_segments__" + acc_name + suffix
    aln_DF = pd.read_csv(in_paf,
                         sep="\t", header=None, usecols=list(range(0, 13)))
    aln_DF["aligner"] = aligner
    aln_DF.columns = ["qname", "qlen", "qstart0", "qend",
                      "strand", "tname", "tlen", "tstart0", "tend",
                      "nmatch", "alen", "mapq", "atype", "aligner"]
    #
    return aln_DF


# Load read segment alignment files as a combined DataFrame
def load_pafs_slowly(indir, acc_name, suffix, aligner):
    #indir=acc1_indir_list[0]
    #acc_name=acc1_name
    #suffix="_alnTo_" + parser.acc1 + "_Chr_mm_ont.paf"
    #aligner="mm"
    find_cmd = ["find"] + \
               [indir + "/"] + \
               ["-type", "f"] + \
               ["-name", "*" + acc_name + suffix] + \
               ["-print"]
    find_output = subprocess.run(find_cmd, capture_output=True)
    files_bytes = find_output.stdout.splitlines()
    files = [x.decode("utf-8") for x in files_bytes]
    #
    if len(files) > 0:
        aln_DF = pd.DataFrame()
        for h in range(0, len(files)):
            aln = pd.read_csv(files[h],
                              sep="\t", header=None, usecols=list(range(0, 13)))
            aln_DF = pd.concat(objs=[aln_DF, aln],
                               axis=0,
                               ignore_index=True)
        #aln_DF = aln_DF.iloc[:, :13]
        aln_DF["aligner"] = aligner
        aln_DF.columns = ["qname", "qlen", "qstart0", "qend0",
                          "strand", "tname", "tlen", "tstart0", "tend",
                          "nmatch", "alen", "mapq", "atype", "aligner"]
        #
        return aln_DF


# wm alignments
acc1_wm_list = []
for x in range(0, len(acc1_indir_list)):
    acc1_wm_Chr = load_cat_paf(indir=acc1_indir_list[x],
                               acc_name=acc1_name,
                               suffix="_alnTo_" + parser.acc1 + "_Chr_wm_ont.paf",
                               aligner="wm")
    acc1_wm_list.append(acc1_wm_Chr)
    del acc1_wm_Chr
    gc.collect()

acc2_wm_list = []
for x in range(0, len(acc2_indir_list)):
    acc2_wm_Chr = load_cat_paf(indir=acc2_indir_list[x],
                               acc_name=acc2_name,
                               suffix="_alnTo_" + parser.acc2 + "_Chr_wm_ont.paf",
                               aligner="wm")
    acc2_wm_list.append(acc2_wm_Chr)
    del acc2_wm_Chr
    gc.collect()

# mm alignments
acc1_mm_list = []
for x in range(0, len(acc1_indir_list)):
    acc1_mm_Chr = load_cat_paf(indir=acc1_indir_list[x],
                               acc_name=acc1_name,
                               suffix="_alnTo_" + parser.acc1 + "_Chr_mm_ont.paf",
                               aligner="mm")
    acc1_mm_list.append(acc1_mm_Chr)
    del acc1_mm_Chr
    gc.collect()

acc2_mm_list = []
for x in range(0, len(acc2_indir_list)):
    acc2_mm_Chr = load_cat_paf(indir=acc2_indir_list[x],
                               acc_name=acc2_name,
                               suffix="_alnTo_" + parser.acc2 + "_Chr_mm_ont.paf",
                               aligner="mm")
    acc2_mm_list.append(acc2_mm_Chr)
    del acc2_mm_Chr
    gc.collect()

# sr alignments
acc1_sr_list = []
for x in range(0, len(acc1_indir_list)):
    acc1_sr_Chr = load_cat_paf(indir=acc1_indir_list[x],
                               acc_name=acc1_name,
                               suffix="_alnTo_" + parser.acc1 + "_Chr_mm_sr.paf",
                               aligner="sr")
    acc1_sr_list.append(acc1_sr_Chr)
    del acc1_sr_Chr
    gc.collect()

acc2_sr_list = []
for x in range(0, len(acc2_indir_list)):
    acc2_sr_Chr = load_cat_paf(indir=acc2_indir_list[x],
                               acc_name=acc2_name,
                               suffix="_alnTo_" + parser.acc2 + "_Chr_mm_sr.paf",
                               aligner="sr")
    acc2_sr_list.append(acc2_sr_Chr)
    del acc2_sr_Chr
    gc.collect()

# bt2 alignments
acc1_bt2_list = []
for x in range(0, len(acc1_indir_list)):
    acc1_bt2_Chr = load_cat_paf(indir=acc1_indir_list[x],
                                acc_name=acc1_name,
                                suffix="_alnTo_" + parser.acc1 + "_Chr_bt2.paf",
                                aligner="bt2")
    acc1_bt2_list.append(acc1_bt2_Chr)
    del acc1_bt2_Chr
    gc.collect()

acc2_bt2_list = []
for x in range(0, len(acc2_indir_list)):
    acc2_bt2_Chr = load_cat_paf(indir=acc2_indir_list[x],
                                acc_name=acc2_name,
                                suffix="_alnTo_" + parser.acc2 + "_Chr_bt2.paf",
                                aligner="bt2")
    acc2_bt2_list.append(acc2_bt2_Chr)
    del acc2_bt2_Chr
    gc.collect()


# acc alignments - chromosome list of aligner lists
acc1_aln_chr_nested_list = []
for x in range(0, len(acc1_mm_list)):
    acc1_aln_chr_nested_list.append([acc1_wm_list[x], acc1_mm_list[x], acc1_sr_list[x]])

acc2_aln_chr_nested_list = []
for x in range(0, len(acc2_mm_list)):
    acc2_aln_chr_nested_list.append([acc2_wm_list[x], acc2_mm_list[x], acc2_sr_list[x]])


# Get best pair of acc1 and acc2 read segment alignments
def aln_best_pair(acc1_aln_DF_list, acc2_aln_DF_list):
    #acc1_aln_DF_list=acc1_aln_chr_nested_list[0]
    #acc2_aln_DF_list=acc2_aln_chr_nested_list[0]
    """
    Get the best pair of acc1 and acc2 read segment alignments, based on:
    1. The alignment number of matching bases (nmatch)
    2. The alignment length (alen)
    3. The alignment quality (mapq)
    4. The alignment type (atype)
    #. The alignment strand
    """
    # Each of the 2 list elements in acc1_aln_DF_list and acc2_aln_DF_list is
    # a DataFrame of alignments done by mm_ont or mm_sr
    acc1_aln_DF_concat = pd.concat(objs=acc1_aln_DF_list, axis=0, ignore_index=True)
    acc2_aln_DF_concat = pd.concat(objs=acc2_aln_DF_list, axis=0, ignore_index=True)
    #
    # Keep alignments with MAPQ >= parser.minMAPQ
    acc1_aln_DF_concat = acc1_aln_DF_concat.loc[acc1_aln_DF_concat["mapq"] >= parser.minMAPQ]
    acc2_aln_DF_concat = acc2_aln_DF_concat.loc[acc2_aln_DF_concat["mapq"] >= parser.minMAPQ]
    # 
    # For each read ID, get the best alignment from each of acc1_aln_DF_concat and
    # and acc2_aln_DF_concat
    # acc1
    acc1_aln_DF_concat_sort = acc1_aln_DF_concat.sort_values(by=["qname", "nmatch", "alen", "mapq", "atype"],
                                                             axis=0,
                                                             ascending=[True, False, False, False, True],
                                                             kind="quicksort",
                                                             ignore_index=True)
    acc1_aln_DF_concat_sort_list = list(acc1_aln_DF_concat_sort.groupby("qname"))
    acc1_aln_DF_best = pd.DataFrame()
    for read_tuple in acc1_aln_DF_concat_sort_list:
    #    if read_tuple[1].loc[read_tuple[1]["aligner"] == "wm"].shape[0] < 10 and \
    #       read_tuple[1].loc[read_tuple[1]["aligner"] == "mm"].shape[0] < 10 and \
    #       read_tuple[1].loc[read_tuple[1]["aligner"] == "sr"].shape[0] < 10:
        read_tuple_aln_DF_best = read_tuple[1].iloc[[0]]
        acc1_aln_DF_best = pd.concat(objs=[acc1_aln_DF_best, read_tuple_aln_DF_best],
                                     axis=0,
                                     ignore_index=True)
    del acc1_aln_DF_concat, acc1_aln_DF_concat_sort, acc1_aln_DF_concat_sort_list, read_tuple_aln_DF_best
    gc.collect()
    #
    # acc2
    acc2_aln_DF_best = pd.DataFrame()
    for read_id in list(acc1_aln_DF_best["qname"]):
        #print(read_id)
        acc1_aln_DF_read_id = acc1_aln_DF_best[acc1_aln_DF_best["qname"] == read_id] 
        acc2_aln_DF_read_id = acc2_aln_DF_concat[acc2_aln_DF_concat["qname"] == read_id] 
        if len(acc2_aln_DF_read_id) > 0:
            acc2_aln_DF_read_id_sort = acc2_aln_DF_read_id.sort_values(by=["nmatch", "alen", "mapq", "atype"],
                                                                       axis=0,
                                                                       ascending=[False, False, False, True],
                                                                       kind="quicksort",
                                                                       ignore_index=True)
            #if acc2_aln_DF_read_id_sort.loc[acc2_aln_DF_read_id_sort["aligner"] == "wm"].shape[0] < 10 and \
            #   acc2_aln_DF_read_id_sort.loc[acc2_aln_DF_read_id_sort["aligner"] == "mm"].shape[0] < 10 and \
            #   acc2_aln_DF_read_id_sort.loc[acc2_aln_DF_read_id_sort["aligner"] == "sr"].shape[0] < 10:
            acc2_aln_DF_read_id_sort_strand = acc2_aln_DF_read_id_sort[acc2_aln_DF_read_id_sort["strand"] == acc1_aln_DF_read_id["strand"].iloc[0]]
            if acc2_aln_DF_read_id_sort_strand.shape[0] > 0:
                #acc2_aln_DF_read_id_sort_select = acc2_aln_DF_read_id_sort_strand.iloc[[0]]
                # If "nmatch" and "alen" are to be given priority over finding pairs where both alignments are to the same strand:
                if acc2_aln_DF_read_id_sort_strand.iloc[0]["nmatch"] == acc2_aln_DF_read_id_sort.iloc[0]["nmatch"] and \
                   acc2_aln_DF_read_id_sort_strand.iloc[0]["alen"] == acc2_aln_DF_read_id_sort.iloc[0]["alen"]:
                    acc2_aln_DF_read_id_sort_select = acc2_aln_DF_read_id_sort_strand.iloc[[0]]
                else:
                    acc2_aln_DF_read_id_sort_select = acc2_aln_DF_read_id_sort.iloc[[0]]
            else:
                acc2_aln_DF_read_id_sort_select = acc2_aln_DF_read_id_sort.iloc[[0]]
            acc2_aln_DF_best = pd.concat(objs=[acc2_aln_DF_best, acc2_aln_DF_read_id_sort_select],
                                         axis=0,
                                         ignore_index=True)
    del acc2_aln_DF_concat, acc1_aln_DF_read_id, acc2_aln_DF_read_id, acc2_aln_DF_read_id_sort
    gc.collect()
    #
    acc1_aln_DF_best.columns = "acc1_" + acc1_aln_DF_best.columns 
    acc2_aln_DF_best.columns = "acc2_" + acc2_aln_DF_best.columns 
    #
    ## Stop if read ID order differs between acc1_aln_DF_best and acc2_aln_DF_best
    ##list(acc1_aln_DF_best["acc1_qname"]) == list(acc2_aln_DF_best["acc2_qname"])
    #if not acc1_aln_DF_best["acc1_qname"].equals(acc2_aln_DF_best["acc2_qname"]):
    #    print("Stopping because read ID order differs between acc1_aln_DF_best and acc2_aln_DF_best")
    #    return
    #
    aln_best_pair_DF = pd.merge(left=acc1_aln_DF_best, right=acc2_aln_DF_best,
                                how="inner", left_on="acc1_qname", right_on="acc2_qname")
    aln_best_pair_DF = aln_best_pair_DF.rename(columns = {"acc1_qname":"qname"})
    aln_best_pair_DF = aln_best_pair_DF.drop(columns="acc2_qname")
    #
    aln_best_pair_DF_sort = aln_best_pair_DF.sort_values(by=["acc1_tname", "acc1_tstart0", "acc1_tend"],
                                                         axis=0,
                                                         ascending=[True, True, True],
                                                         kind="quicksort",
                                                         ignore_index=True)
    del acc1_aln_DF_best, acc2_aln_DF_best, aln_best_pair_DF
    gc.collect()
    #
    return aln_best_pair_DF_sort


# Get best pair of aligned read segments for each read
aln_best_pair_DF_list = []
for x in range(0, len(acc1_aln_chr_nested_list)):
    aln_best_pair_DF_x = aln_best_pair(acc1_aln_DF_list=acc1_aln_chr_nested_list[x],
                                       acc2_aln_DF_list=acc2_aln_chr_nested_list[x])
    aln_best_pair_DF_list.append(aln_best_pair_DF_x)

# Concatenate list elements (per-chromosome pd.DataFrames) into one pd.DataFrame
if len(aln_best_pair_DF_list) > 1:
    aln_best_pair_DF = pd.concat(objs=aln_best_pair_DF_list, axis=0, ignore_index=True)
else:
    aln_best_pair_DF = aln_best_pair_DF_list[0]


# Report total
str(aln_best_pair_DF.shape[0]) + " validly aligning hybrid read segment pairs where unaligned segments are of '" + parser.recombType + " type', \n" + \
    "based on the sequence of accession-specific, chromosome-specific k-mers"

# Filter to retain hybrid read segments pairs where each segment aligns to the same chromosome
aln_best_pair_hom_DF = aln_best_pair_DF.loc[aln_best_pair_DF["acc1_tname"] == aln_best_pair_DF["acc2_tname"]]
str(aln_best_pair_hom_DF.shape[0]) + " '" + parser.recombType + "-type' hybrid read segments align to the same chromosome"
str( round( aln_best_pair_hom_DF.shape[0] / aln_best_pair_DF.shape[0], 4 ) * 100 ) + "% of '" + parser.recombType + "-type' hybrid read segment pairs align to the same chromosome"


# Get number of hybrid reads and read lengths for hybrid read IDs in aln_best_pair_hom_DF["qname"]
hybrid_reads_counter = 0
hybrid_read_lengths_DF_list = []
for x in range(0, len(chrom)):
    reads_fa = "fasta/" + parser.readsPrefix + \
        "_match_" + parser.acc1 + "_" + parser.region + "_" + chrom[x] + \
        "_specific_k" + str(parser.kmerSize) + "_downsampled_op" + str(parser.overlapProp) + "_hits" + str(parser.minHits) + \
        "_match_" + parser.acc2 + "_" + parser.region + "_" + chrom[x] + \
        "_specific_k" + str(parser.kmerSize) + "_downsampled_op" + str(parser.overlapProp) + "_hits" + str(parser.minHits) + \
        ".fa"
    reads_dict = SeqIO.index(reads_fa, "fasta")
    hybrid_reads_counter += len(reads_dict)
    reads = [v for i, v in enumerate(reads_dict.values()) if v.id in list(aln_best_pair_hom_DF["qname"])]
    #reads_iter = SeqIO.parse(reads_fa, "fasta")
    #reads2 = [x for x in reads_iter if x.id in list(aln_best_pair_hom_DF["qname"])]
    read_ids = [x.id for x in reads]
    read_lens = [len(x) for x in reads]
    hybrid_read_lengths_DF = pd.DataFrame({ "read_id": read_ids,
                                            "read_len": read_lens })
    hybrid_read_lengths_DF_list.append(hybrid_read_lengths_DF)

del hybrid_read_lengths_DF
gc.collect()

# Concatenate list elements (per-chromosome pd.DataFrames) into one pd.DataFrame
if len(hybrid_read_lengths_DF_list) > 1:
    hybrid_read_lengths_DF = pd.concat(objs=hybrid_read_lengths_DF_list, axis=0, ignore_index=True)
else:
    hybrid_read_lengths_DF = hybrid_read_lengths_DF_list[0]

aln_best_pair_hom_DF_read_lens = pd.merge(left=aln_best_pair_hom_DF, right=hybrid_read_lengths_DF,
                                          how="inner", left_on="qname", right_on="read_id")
del aln_best_pair_hom_DF
gc.collect()
aln_best_pair_hom_DF = aln_best_pair_hom_DF_read_lens.drop(columns="read_id")


# Filter to retain hybrid read segments pairs where the per-accession read segments align to within
# the centromere coordinates of the respective genome assembly
aln_best_pair_hom_maxDist_DF = pd.DataFrame()
for x in range(0, len(chrom)):
    print(chrom[x])
    acc1_CEN_chrom = acc1_CEN.loc[acc1_CEN["chr"] == chrom[x]]
    acc2_CEN_chrom = acc2_CEN.loc[acc2_CEN["chr"] == chrom[x]]
    aln_best_pair_hom_maxDist_DF_chrom = aln_best_pair_hom_DF.loc[ \
                                                                   (aln_best_pair_hom_DF["acc1_tname"] == chrom[x]) & \
                                                                   (aln_best_pair_hom_DF["acc1_tstart0"] + 1 <= int(acc1_CEN_chrom["end"])) & \
                                                                   (aln_best_pair_hom_DF["acc1_tend"] >= int(acc1_CEN_chrom["start"])) & \
                                                                   (aln_best_pair_hom_DF["acc2_tname"] == chrom[x]) & \
                                                                   (aln_best_pair_hom_DF["acc2_tstart0"] + 1 <= int(acc2_CEN_chrom["end"])) & \
                                                                   (aln_best_pair_hom_DF["acc2_tend"] >= int(acc2_CEN_chrom["start"])) \
                                                                  ] 
    if aln_best_pair_hom_maxDist_DF_chrom.shape[0] > 0:
        aln_best_pair_hom_maxDist_DF = pd.concat(objs=[aln_best_pair_hom_maxDist_DF, aln_best_pair_hom_maxDist_DF_chrom],
                                                 axis=0, ignore_index=True)


# Filter to retain putative recombination events between homologous chromosomes where
# the per-accession read segments align to within maxDist bp of each other in the same reference assembly, AND
# where the per-accession alignment length is >= parser.alenTOqlen of the segment length
aln_best_pair_hom_maxDist_alenTOqlen_DF = aln_best_pair_hom_maxDist_DF.loc[ ( aln_best_pair_hom_maxDist_DF["acc1_alen"] / \
                                                                              aln_best_pair_hom_maxDist_DF["acc1_qlen"] >= parser.alenTOqlen ) & \
                                                                            ( aln_best_pair_hom_maxDist_DF["acc2_alen"] / \
                                                                              aln_best_pair_hom_maxDist_DF["acc2_qlen"] >= parser.alenTOqlen ) ]


# Filter to retain putative recombination events between homologous chromosomes where
# the per-accession read segments align to within maxDist bp of each other in the same reference assembly, AND
# where the per-accession alignment number of matching bases is >= parser.alenTOqlen of the segment length
aln_best_pair_hom_maxDist_nmatchTOqlen_DF = aln_best_pair_hom_maxDist_DF.loc[ ( aln_best_pair_hom_maxDist_DF["acc1_nmatch"] / \
                                                                                aln_best_pair_hom_maxDist_DF["acc1_qlen"] >= parser.nmatchTOqlen ) & \
                                                                              ( aln_best_pair_hom_maxDist_DF["acc2_nmatch"] / \
                                                                                aln_best_pair_hom_maxDist_DF["acc2_qlen"] >= parser.nmatchTOqlen ) ]


# Make summary of counts obtained at each filtering stage
summary_DF = pd.DataFrame({ "Filter": [ "hybrid_reads",
                                        "hybrid_read_segment_pairs_aligned",
                                        "hybrid_read_segment_pairs_aligned_to_same_chr",
                                        "hybrid_read_segment_pairs_aligned_to_within_maxDist",
                                        "hybrid_read_segment_pairs_aligned_chr_maxDist_alenTOqlen",
                                        "hybrid_read_segment_pairs_aligned_chr_maxDist_nmatchTOqlen" ],
                            "Count": [ hybrid_reads_counter,
                                       aln_best_pair_DF.shape[0],
                                       aln_best_pair_hom_DF.shape[0],
                                       aln_best_pair_hom_maxDist_DF.shape[0],
                                       aln_best_pair_hom_maxDist_alenTOqlen_DF.shape[0],
                                       aln_best_pair_hom_maxDist_nmatchTOqlen_DF.shape[0] ],
                            "Proportion": [ hybrid_reads_counter / hybrid_reads_counter,
                                            aln_best_pair_DF.shape[0] / hybrid_reads_counter,
                                            aln_best_pair_hom_DF.shape[0] / aln_best_pair_DF.shape[0],
                                            aln_best_pair_hom_maxDist_DF.shape[0] / aln_best_pair_hom_DF.shape[0],
                                            aln_best_pair_hom_maxDist_alenTOqlen_DF.shape[0] / aln_best_pair_hom_maxDist_DF.shape[0],
                                            aln_best_pair_hom_maxDist_nmatchTOqlen_DF.shape[0] / aln_best_pair_hom_maxDist_DF.shape[0] ] })


# Write to TSV
summary_DF_filename = outdir + "/" + parser.readsPrefix + \
    "_" + parser.acc1 + "_" + parser.acc2 + \
    "_k" + str(parser.kmerSize) + "_op" + str(parser.overlapProp) + "_h" + str(parser.minHits) + \
    "_" + parser.recombType + \
    "_" +  \
    re.sub(",", "_", parser.chrom) + "_count_summary.tsv"
summary_DF.to_csv(summary_DF_filename, sep="\t", header=True, index=False)

aln_best_pair_DF_filename = outdir + "/" + parser.readsPrefix + \
    "_" + parser.acc1 + "_" + parser.acc2 + \
    "_k" + str(parser.kmerSize) + "_op" + str(parser.overlapProp) + "_h" + str(parser.minHits) + \
    "_" + parser.recombType + \
    "_" +  \
    re.sub(",", "_", parser.chrom) + ".tsv"
aln_best_pair_DF.to_csv(aln_best_pair_DF_filename, sep="\t", header=True, index=False)

aln_best_pair_hom_DF_filename = outdir + "/" + parser.readsPrefix + \
    "_" + parser.acc1 + "_" + parser.acc2 + \
    "_k" + str(parser.kmerSize) + "_op" + str(parser.overlapProp) + "_h" + str(parser.minHits) + \
    "_hom_" + parser.recombType + \
    "_" +  \
    re.sub(",", "_", parser.chrom) + ".tsv"
aln_best_pair_hom_DF.to_csv(aln_best_pair_hom_DF_filename, sep="\t", header=True, index=False)

aln_best_pair_hom_maxDist_DF_filename = outdir + "/" + parser.readsPrefix + \
    "_" + parser.acc1 + "_" + parser.acc2 + \
    "_k" + str(parser.kmerSize) + "_op" + str(parser.overlapProp) + "_h" + str(parser.minHits) + \
    "_hom_maxDist_" + parser.recombType + \
    "_" +  \
    re.sub(",", "_", parser.chrom) + ".tsv"
aln_best_pair_hom_maxDist_DF.to_csv(aln_best_pair_hom_maxDist_DF_filename, sep="\t", header=True, index=False)

aln_best_pair_hom_maxDist_alenTOqlen_DF_filename = outdir + "/" + parser.readsPrefix + \
    "_" + parser.acc1 + "_" + parser.acc2 + \
    "_k" + str(parser.kmerSize) + "_op" + str(parser.overlapProp) + "_h" + str(parser.minHits) + \
    "_hom_maxDist_aTOq" + str(parser.alenTOqlen) + "_" + parser.recombType + \
    "_" +  \
    re.sub(",", "_", parser.chrom) + ".tsv"
aln_best_pair_hom_maxDist_alenTOqlen_DF.to_csv(aln_best_pair_hom_maxDist_alenTOqlen_DF_filename, sep="\t", header=True, index=False)

aln_best_pair_hom_maxDist_nmatchTOqlen_DF_filename = outdir + "/" + parser.readsPrefix + \
    "_" + parser.acc1 + "_" + parser.acc2 + \
    "_k" + str(parser.kmerSize) + "_op" + str(parser.overlapProp) + "_h" + str(parser.minHits) + \
    "_hom_maxDist_nTOq" + str(parser.nmatchTOqlen) + "_" + parser.recombType + \
    "_" +  \
    re.sub(",", "_", parser.chrom) + ".tsv"
aln_best_pair_hom_maxDist_nmatchTOqlen_DF.to_csv(aln_best_pair_hom_maxDist_nmatchTOqlen_DF_filename, sep="\t", header=True, index=False)


# Write filtered hybrid reads to FASTA
for x in range(0, len(chrom)):
    # Input FASTA path to all hybrid reads for the given chromosome
    reads_fa = "fasta/" + parser.readsPrefix + \
        "_match_" + parser.acc1 + "_" + parser.region + "_" + chrom[x] + \
        "_specific_k" + str(parser.kmerSize) + "_downsampled_op" + str(parser.overlapProp) + "_hits" + str(parser.minHits) + \
        "_match_" + parser.acc2 + "_" + parser.region + "_" + chrom[x] + \
        "_specific_k" + str(parser.kmerSize) + "_downsampled_op" + str(parser.overlapProp) + "_hits" + str(parser.minHits) + \
        ".fa"
    reads_dict = SeqIO.index(reads_fa, "fasta")
    # Get the reads where the segments align to within the given distance of each other
    reads_maxDist = [v for i, v in enumerate(reads_dict.values()) if v.id in list(aln_best_pair_hom_maxDist_DF["qname"])]
    # Get the reads where the segments align to within the given distance of each other, and
    # where the alen / qlen ratio >= parser.alenTOqlen
    reads_maxDist_alenTOqlen = [v for i, v in enumerate(reads_dict.values()) if v.id in list(aln_best_pair_hom_maxDist_alenTOqlen_DF["qname"])]
    # Get the reads where the segments align to within the given distance of each other, and
    # where the nmatch / qlen ratio >= parser.nmatchTOqlen
    reads_maxDist_nmatchTOqlen = [v for i, v in enumerate(reads_dict.values()) if v.id in list(aln_best_pair_hom_maxDist_nmatchTOqlen_DF["qname"])]
    # Output FASTA path to maxDist-filtered hybrid reads for the given chromosome
    reads_maxDist_fa = outdir + "/" + parser.readsPrefix + \
        "_" + parser.acc1 + "_" + parser.acc2 + \
        "_k" + str(parser.kmerSize) + "_op" + str(parser.overlapProp) + "_h" + str(parser.minHits) + \
        "_hom_maxDist_" + parser.recombType + \
        "_" +  \
        chrom[x] + "_hybrid_reads.fa"
    # Output FASTA path to maxDist_alenTOqlen-filtered hybrid reads for the given chromosome
    reads_maxDist_alenTOqlen_fa = outdir + "/" + parser.readsPrefix + \
        "_" + parser.acc1 + "_" + parser.acc2 + \
        "_k" + str(parser.kmerSize) + "_op" + str(parser.overlapProp) + "_h" + str(parser.minHits) + \
        "_hom_maxDist_aTOq" + str(parser.alenTOqlen) + "_" + parser.recombType + \
        "_" +  \
        chrom[x] + "_hybrid_reads.fa"
    # Output FASTA path to maxDist_nmatchTOqlen-filtered hybrid reads for the given chromosome
    reads_maxDist_nmatchTOqlen_fa = outdir + "/" + parser.readsPrefix + \
        "_" + parser.acc1 + "_" + parser.acc2 + \
        "_k" + str(parser.kmerSize) + "_op" + str(parser.overlapProp) + "_h" + str(parser.minHits) + \
        "_hom_maxDist_nTOq" + str(parser.nmatchTOqlen) + "_" + parser.recombType + \
        "_" +  \
        chrom[x] + "_hybrid_reads.fa"
    # Write outputs
    with open(reads_maxDist_fa, "w") as reads_maxDist_fa_handle:
        SeqIO.write(reads_maxDist, reads_maxDist_fa_handle, "fasta")
    with open(reads_maxDist_alenTOqlen_fa, "w") as reads_maxDist_alenTOqlen_fa_handle:
        SeqIO.write(reads_maxDist_alenTOqlen, reads_maxDist_alenTOqlen_fa_handle, "fasta")
    with open(reads_maxDist_nmatchTOqlen_fa, "w") as reads_maxDist_nmatchTOqlen_fa_handle:
        SeqIO.write(reads_maxDist_nmatchTOqlen, reads_maxDist_nmatchTOqlen_fa_handle, "fasta")

