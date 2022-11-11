#!/usr/bin/env Rscript

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 09/11/2022

# Extract putative recombination events detected in Col-0/Ler-0
# hybrid ONT reads identified on the basis that they contain >=
# a threshold number of accession-specific, chromosome-specific k-mers

# Usage:
# conda activate python_3.9.6
# ./get_chr_specific_segment_pairs_alnToSame.R \
#  -r ColLerF1pollen_1000bp_minq90 \
#  -a1 Col-0.ragtag_scaffolds \
#  -a2 Ler-0_110x.ragtag_scaffolds \
#  -k 24 \
#  -op 0.9 \
#  -mh 11 \
#  -at Col-0.ragtag_scaffolds_Chr \
#  -aq 0.90 \
#  -rt co \
#  -reg not_centromere \
#  -c 'Chr5'
# conda deactivate


# ==== Import libraries
import os
import argparse
import re
import gc
import pandas as pd
import subprocess
import glob

from pathlib import Path


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
    parser.add_argument("-at", "--alnTo", type=str, default="Col-0.ragtag_scaffolds_Chr",
                        help="The prefix of the assembly to be used for read segment alignment. Default: Col-0.ragtag_scaffolds_Chr")
    parser.add_argument("-aq", "--alenTOqlen", type=float, default="0.90",
                        help="The minimum ratio of the read segment alignment length to the read segment length. Default: 0.90")
    parser.add_argument("-rt", "--recombType", type=str, default="co",
                        help="The type/pattern of the recombination event identified based on the sequence of accession-specific read segments. Default: co")
    parser.add_argument("-reg", "--region", type=str, default="not_centromere",
                        help="The chromosome for which accession-specific, chromosome-specific read segments have been extracted and aligned. Default: not_centromere")
    parser.add_argument("-c", "--chrom", type=str, default="Chr5",
                        help="The chromosome for which accession-specific, chromosome-specific read segments have been extracted and aligned. Default: Chr5")
    return parser

parser = create_parser().parse_args()
print(parser)

chrom = parser.chrom.split(",")

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


# Concatenate all alignment files for the given chromosome,
# accession and aligner
# NOTE: for some reason the "find ..." approach below doesn't work
# from within python, so need to create and run an equivalent bash script:
#cat_cmd = ["find"] + \
#          ["/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + indir + "/"] + \
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
    #suffix="_alnTo_" + parser.alnTo + "_mm_ont.paf"
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
    cat_pafs(indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc1_indir_list[x],
             acc_name=acc1_name,
             suffix="_alnTo_" + parser.alnTo + "_mm_ont.paf")
    cat_pafs(indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc1_indir_list[x],
             acc_name=acc1_name,
             suffix="_alnTo_" + parser.alnTo + "_mm_sr.paf")

for x in range(0, len(acc2_indir_list)):
    print(acc2_indir_list[x])
    cat_pafs(indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc2_indir_list[x],
             acc_name=acc2_name,
             suffix="_alnTo_" + parser.alnTo + "_mm_ont.paf")
    cat_pafs(indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc2_indir_list[x],
             acc_name=acc2_name,
             suffix="_alnTo_" + parser.alnTo + "_mm_sr.paf")


# Load concatenated read segment alignment file as a DataFrame
def load_cat_paf(indir, acc_name, suffix, aligner):
    ##indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc1_indir_list[1]
    #indir=acc1_indir_list[1]
    #acc_name=acc1_name
    #suffix="_alnTo_" + parser.alnTo + "_mm_ont.paf"
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
    aln_DF.columns = ["qname", "qlen", "qstart0", "qend0",
                      "strand", "tname", "tlen", "tstart", "tend",
                      "nmatch", "alen", "mapq", "atype", "aligner"]
    #
    return aln_DF


# Load read segment alignment files as a combined DataFrame
def load_pafs_slowly(indir, acc_name, suffix, aligner):
    #indir=acc1_indir_list[0]
    #acc_name=acc1_name
    #suffix="_alnTo_" + parser.alnTo + "_mm_ont.paf"
    #aligner="mm"
    find_cmd = ["find"] + \
               ["/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + indir + "/"] + \
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
                          "strand", "tname", "tlen", "tstart", "tend",
                          "nmatch", "alen", "mapq", "atype", "aligner"]
        #
        return aln_DF


# mm alignments
acc1_mm_list = []
for x in range(0, len(acc1_indir_list)):
    acc1_mm_Chr = load_cat_paf(indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc1_indir_list[x],
                               acc_name=acc1_name,
                               suffix="_alnTo_" + parser.alnTo + "_mm_ont.paf",
                               aligner="mm")
    acc1_mm_list.append(acc1_mm_Chr)
    del acc1_mm_Chr
    gc.collect()

acc2_mm_list = []
for x in range(0, len(acc2_indir_list)):
    acc2_mm_Chr = load_cat_paf(indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc2_indir_list[x],
                               acc_name=acc2_name,
                               suffix="_alnTo_" + parser.alnTo + "_mm_ont.paf",
                               aligner="mm")
    acc2_mm_list.append(acc2_mm_Chr)
    del acc2_mm_Chr
    gc.collect()

# sr alignments
acc1_sr_list = []
for x in range(0, len(acc1_indir_list)):
    acc1_sr_Chr = load_cat_paf(indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc1_indir_list[x],
                               acc_name=acc1_name,
                               suffix="_alnTo_" + parser.alnTo + "_mm_sr.paf",
                               aligner="sr")
    acc1_sr_list.append(acc1_sr_Chr)
    del acc1_sr_Chr
    gc.collect()

acc2_sr_list = []
for x in range(0, len(acc2_indir_list)):
    acc2_sr_Chr = load_cat_paf(indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc2_indir_list[x],
                               acc_name=acc2_name,
                               suffix="_alnTo_" + parser.alnTo + "_mm_sr.paf",
                               aligner="sr")
    acc2_sr_list.append(acc2_sr_Chr)
    del acc2_sr_Chr
    gc.collect()


# acc alignments - chromosome list of aligner lists
acc1_aln_chr_nested_list = []
for x in range(0, len(acc1_mm_list)):
    acc1_aln_chr_nested_list.append([acc1_mm_list[x], acc1_sr_list[x]])

acc2_aln_chr_nested_list = []
for x in range(0, len(acc2_mm_list)):
    acc2_aln_chr_nested_list.append([acc2_mm_list[x], acc2_sr_list[x]])


# Get best pair of acc1 and acc2 read segment alignments, based on:
# 1. The alignment length (alen)
# 2. The alignment number of matching bases (nmatch)
# 3. The alignment strand
def aln_best_pair(acc1_aln_DF_list, acc2_aln_DF_list):
acc1_aln_DF_list=acc1_aln_chr_nested_list[0]
acc2_aln_DF_list=acc2_aln_chr_nested_list[0]
"""
Get the best pair of acc1 and acc2 read segment alignments, based on:
1. The alignment length (alen)
2. The alignment number of matching bases (nmatch)
3. The alignment strand
"""
# Each of the 2 list elements in acc1_aln_DF_list and acc2_aln_DF_list is
# a DataFrame of alignments done by mm_ont or mm_sr
acc1_aln_DF_concat = pd.concat(objs=acc1_aln_DF_list, axis=0, ignore_index=True)
acc2_aln_DF_concat = pd.concat(objs=acc2_aln_DF_list, axis=0, ignore_index=True)
#
# For each read ID, get the best alignment from each of acc1_aln_DF_concat and
# and acc2_aln_DF_concat
acc1_aln_DF_concat_sort = acc1_aln_DF_concat.sort_values(by=["qname", "alen", "nmatch"],
                                                         axis=0,
                                                         ascending=[True, False, False],
                                                         kind="quicksort",
                                                         ignore_index=False)
acc1_aln_DF_concat_sort_list = list(acc1_aln_DF_concat_sort.groupby("qname"))

    kmers_DF_sort_list = list(kmers_DF_sort.groupby(["win_chr", "win_start0", "win_end"]))


acc1_aln_DF_best = pd.DataFrame()
acc2_aln_DF_best = pd.DataFrame()
for read_id in list(set(acc1_aln_DF_concat["qname"])):
acc1_aln_DF_read_id = acc1_aln_DF_concat[acc1_aln_DF_concat["qname"] == read_id] 
acc1_aln_DF_read_id_sort = acc1_aln_DF_read_id.sort_values(by=["alen", "nmatch"],
                                                           axis=0,
                                                           ascending=[False, False],
                                                           kind="quicksort",
                                                           ignore_index=False)
acc1_aln_DF_read_id_sort_select = acc1_aln_DF_read_id_sort.iloc[[0]]
acc1_aln_DF_best = pd.concat(objs=[acc1_aln_DF_best, acc1_aln_DF_read_id_sort_select],
                             axis=0,
                             ignore_index=True)
  
    kmers_DF_sort = kmers_DF.sort_values(by=["win_chr", "win_start0", "overlap_bp"],
                                         axis=0,
                                         ascending=[True, True, False],
                                         kind="quicksort",
                                         ignore_index=True)
    kmers_DF_sort_list = list(kmers_DF_sort.groupby(["win_chr", "win_start0", "win_end"]))
    kmers_DF_greatest = pd.DataFrame()
    for window_tuple in kmers_DF_sort_list:
        window_DF_greatest = window_tuple[1].iloc[[0]]
        kmers_DF_greatest = pd.concat(objs=[kmers_DF_greatest, window_DF_greatest],
                                      axis = 0,
                                      ignore_index=True)
    kmers_DF_greatest_sort = kmers_DF_greatest.iloc[:,3:10].sort_values(by=["kmer_chr", "kmer_start0"],
                                                                        axis=0,
                                                                        ascending=[True, True],
                                                                        kind="quicksort",
                                                                        ignore_index=True)
    kmers_DF_greatest_sort.to_csv(out_bed, sep="\t", header=False, index=False)








    for(read_id in unique(acc1_aln_DF_concat$qname)) {
        acc1_aln_DF_read_id = acc1_aln_DF_concat %>%
            dplyr::filter(qname == read_id)
        acc1_aln_DF_read_id = acc1_aln_DF_read_id[ with(acc1_aln_DF_read_id,
                                                        order(alen, nmatch, decreasing=T)), ]
        acc1_aln_DF_read_id_select = acc1_aln_DF_read_id[1, ]
        acc1_aln_DF_best = rbind(acc1_aln_DF_best, acc1_aln_DF_read_id_select)

        acc2_aln_DF_read_id = acc2_aln_DF_concat %>%
            dplyr::filter(qname == read_id)
        acc2_aln_DF_read_id = acc2_aln_DF_read_id[ with(acc2_aln_DF_read_id,
                                                        order(alen, nmatch, decreasing=T)), ]
        acc2_aln_DF_read_id_strand = acc2_aln_DF_read_id[ which(acc2_aln_DF_read_id$strand == acc1_aln_DF_read_id_select$strand), ]
        if(nrow(acc2_aln_DF_read_id_strand) > 0) {
            if(acc2_aln_DF_read_id_strand[1, ]$alen == acc2_aln_DF_read_id[1, ]$alen &
               acc2_aln_DF_read_id_strand[1, ]$nmatch == acc2_aln_DF_read_id[1, ]$nmatch) {
                acc2_aln_DF_read_id_select = acc2_aln_DF_read_id_strand[1, ]
            } else {
                acc2_aln_DF_read_id_select = acc2_aln_DF_read_id[1, ]
            }
        } else {
            acc2_aln_DF_read_id_select = acc2_aln_DF_read_id[1, ]
        }
        acc2_aln_DF_best = rbind(acc2_aln_DF_best, acc2_aln_DF_read_id_select)
    }

    colnames(acc1_aln_DF_best) = paste0("acc1_", colnames(acc1_aln_DF_best))
    colnames(acc2_aln_DF_best) = paste0("acc2_", colnames(acc2_aln_DF_best))
    stopifnot(identical(acc1_aln_DF_best$acc1_qname, acc2_aln_DF_best$acc2_qname))

    #aln_best_pair_DF = cbind(acc1_aln_DF_best, acc2_aln_DF_best)
    #aln_best_pair_DF = aln_best_pair_DF[ , -which(colnames(aln_best_pair_DF) == "acc2_qname") ]

    aln_best_pair_DF = base::merge(x = acc1_aln_DF_best,
                                   y = acc2_aln_DF_best,
                                   by.x = "acc1_qname",
                                   by.y = "acc2_qname")

    colnames(aln_best_pair_DF)[which(colnames(aln_best_pair_DF) == "acc1_qname")] = "qname"
    aln_best_pair_DF = aln_best_pair_DF[ with(aln_best_pair_DF,
                                              order(acc1_tname, acc1_tstart, acc1_tend, decreasing=F)), ]

    return(aln_best_pair_DF)
}


# Get best pair of aligned read segments for each read
aln_best_pair_DF = dplyr::bind_rows(
    mclapply(1:length(acc1_aln_chr_nested_list), function(x) {
        aln_best_pair(acc1_aln_DF_list=acc1_aln_chr_nested_list[[x]], acc2_aln_DF_list=acc2_aln_chr_nested_list[[x]])
    }, mc.preschedule=F, mc.cores=length(acc1_aln_chr_nested_list))
)

#aln_best_pair_DF = aln_best_pair_DF %>%
#    dplyr::filter(acc1_alen >= 200) %>%
#    dplyr::filter(acc2_alen >= 200)

print(paste0(nrow(aln_best_pair_DF), " putative ", recombType, " events"))

# Filter to retain putative recombination events between homologous chromosomes only
aln_best_pair_hom_DF = aln_best_pair_DF[ which(aln_best_pair_DF$acc1_tname == aln_best_pair_DF$acc2_tname), ]

print(paste0(nrow(aln_best_pair_hom_DF), " putative ", recombType, " events are between homologous chromosomes"))

print(paste0( round( ( nrow(aln_best_pair_hom_DF) / nrow(aln_best_pair_DF) ), 4 ) * 100, "% of putative ", recombType, " events are between homologous chromosomes"))

# Filter to retain putative recombination events between homologous chromosomes where
# the per-accession read segments align to within maxDist bp of each other in the same reference assembly
aln_dist_acc1_tstart_acc2_tstart = abs(aln_best_pair_hom_DF$acc1_tstart - aln_best_pair_hom_DF$acc2_tstart) + 1
aln_dist_acc1_tstart_acc2_tend = abs(aln_best_pair_hom_DF$acc1_tstart - aln_best_pair_hom_DF$acc2_tend) + 1
aln_dist_acc1_tend_acc2_tstart = abs(aln_best_pair_hom_DF$acc1_tend - aln_best_pair_hom_DF$acc2_tstart) + 1
aln_dist_acc1_tend_acc2_tend = abs(aln_best_pair_hom_DF$acc1_tend - aln_best_pair_hom_DF$acc2_tend) + 1

aln_dist_min = pmin(aln_dist_acc1_tstart_acc2_tstart, aln_dist_acc1_tstart_acc2_tend,
                    aln_dist_acc1_tend_acc2_tstart, aln_dist_acc1_tend_acc2_tend, na.rm = T)
aln_dist_max = pmax(aln_dist_acc1_tstart_acc2_tstart, aln_dist_acc1_tstart_acc2_tend,
                    aln_dist_acc1_tend_acc2_tstart, aln_dist_acc1_tend_acc2_tend, na.rm = T)

event_start = pmin(aln_best_pair_hom_DF$acc1_tstart, aln_best_pair_hom_DF$acc1_tend,
                   aln_best_pair_hom_DF$acc2_tstart, aln_best_pair_hom_DF$acc2_tend, na.rm = T)
event_end = pmax(aln_best_pair_hom_DF$acc1_tstart, aln_best_pair_hom_DF$acc1_tend,
                 aln_best_pair_hom_DF$acc2_tstart, aln_best_pair_hom_DF$acc2_tend, na.rm = T)
event_midpoint = round((event_start + event_end) / 2)

aln_best_pair_hom_DF$aln_dist_min = aln_dist_min
aln_best_pair_hom_DF$aln_dist_max = aln_dist_max
aln_best_pair_hom_DF$event_start = event_start
aln_best_pair_hom_DF$event_end = event_end
aln_best_pair_hom_DF$event_midpoint = event_midpoint

# Get read lengths for hybrid read IDs in aln_best_pair_hom_DF$qname
# to be used for retaining alignment pairs where the Col and Ler
# read segments align to within a given distance of each other
# (e.g., the given hybrid read length) in the same assembly
hybrid_read_lengths = distinct(dplyr::bind_rows(
    mclapply(1:length(chrName), function(x) {
        tmp_list = read.fasta(paste0("fasta/", readsPrefix,
                                     "_match_", acc1, "_", region, "_", chrName[x],
                                     "_specific_k", kmerSize, "_downsampled_op", overlapProp, "_hits", minHits,
                                     "_match_", acc2, "_", region, "_", chrName[x],
                                     "_specific_k", kmerSize, "_downsampled_op", overlapProp, "_hits", minHits,
                                     ".fa"),
                              forceDNAtolower=F)
        read_length_DF = data.frame()
        for(i in 1:length(tmp_list)) {
            if( attr(tmp_list[[i]], which="name", exact=T) %in% aln_best_pair_hom_DF$qname ) {
                read_length_DF_i = data.frame(read_id = attr(tmp_list[[i]], which="name", exact=T),
                                              read_len = length(getSequence(tmp_list[[i]])))
                read_length_DF = rbind(read_length_DF, read_length_DF_i)
            } 
        }                                                              
        read_length_DF
    }, mc.preschedule=F, mc.cores=length(chrName))
))

aln_best_pair_hom_DF = base::merge(x = aln_best_pair_hom_DF,
                                   y = hybrid_read_lengths,
                                   by.x = "qname",
                                   by.y = "read_id",
                                   sort = F)

aln_best_pair_hom_maxDist_DF = data.frame()
for(x in 1:nrow(aln_best_pair_hom_DF)) {
    if(aln_best_pair_hom_DF[x, ]$aln_dist_min <= aln_best_pair_hom_DF[x, ]$read_len * 2) {
        aln_best_pair_hom_maxDist_DF = rbind(aln_best_pair_hom_maxDist_DF, aln_best_pair_hom_DF[x, ])
    }
}

print(paste0(nrow(aln_best_pair_hom_maxDist_DF), " putative ", recombType, " events are between homologous chromosomes where the per-accession read segments align to within hybrid-read-length bp of each other in the same reference assembly"))

print(paste0( round( ( nrow(aln_best_pair_hom_maxDist_DF) / nrow(aln_best_pair_DF) ), 4 ) * 100, "% of putative ", recombType, " events are between homologous chromosomes where the per-accession read segments align to within hybrid-read-length bp of each other in the same reference assembly"))

print(paste0( round( ( nrow(aln_best_pair_hom_maxDist_DF) / nrow(aln_best_pair_hom_DF) ), 4 ) * 100, "% of putative ", recombType, " events that are between homologous chromosomes where the per-accession read segments align to within hybrid-read-length bp of each other in the same reference assembly"))



# Filter to retain putative recombination events between homologous chromosomes where
# the per-accession read segments align to within maxDist bp of each other in the same reference assembly, AND
# where the per-accession alignment length is >= alenTOqlen of the segment length
aln_best_pair_hom_maxDist_alenTOqlen_DF = aln_best_pair_hom_maxDist_DF[ which(aln_best_pair_hom_maxDist_DF$acc1_alen /
                                                                              aln_best_pair_hom_maxDist_DF$acc1_qlen >= alenTOqlen &
                                                                              aln_best_pair_hom_maxDist_DF$acc2_alen /
                                                                              aln_best_pair_hom_maxDist_DF$acc2_qlen >= alenTOqlen), ]

print(paste0(nrow(aln_best_pair_hom_maxDist_alenTOqlen_DF), " putative ", recombType, " events are between homologous chromosomes where the per-accession read segments align to within hybrid-read-length bp of each other in the same reference assembly, and where the per-accession alignment length is >= ", alenTOqlen * 100, "% of the segment length"))

print(paste0( round( ( nrow(aln_best_pair_hom_maxDist_alenTOqlen_DF) / nrow(aln_best_pair_DF) ), 4 ) * 100, "% of putative ", recombType, " events are between homologous chromosomes where the per-accession read segments align to within hybrid-read-length bp of each other in the same reference assembly, and where the per-accession alignment length is >= ", alenTOqlen * 100, "% of the segment length"))

print(paste0( round( ( nrow(aln_best_pair_hom_maxDist_alenTOqlen_DF) / nrow(aln_best_pair_hom_maxDist_DF) ), 4 ) * 100, "% of putative ", recombType, " events that are between homologous chromosomes and where the per-accession read segments align to within hybrid-read-length bp of each other in the same reference assembly, are those where the per-accession alignment length is >= ", alenTOqlen * 100, "% of the segment length"))



#write.table(aln_best_pair_DF,
#            paste0(outdir, readsPrefix,
#                   "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
#                   "_", recombType,
#                   "_alnTo_", alnTo, "_",
#                   paste0(chrName, collapse="_"), ".tsv"),
#            quote=F, sep="\t", col.names=T, row.names=F)
#write.table(aln_best_pair_hom_DF,
#            paste0(outdir, readsPrefix,
#                   "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
#                   "_hom_", recombType,
#                   "_alnTo_", alnTo, "_",
#                   paste0(chrName, collapse="_"), ".tsv"),
#            quote=F, sep="\t", col.names=T, row.names=F)
write.table(aln_best_pair_hom_maxDist_DF,
            paste0(outdir, readsPrefix,
                   "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                   "_hom_maxDist_", recombType,
                   "_alnTo_", alnTo, "_",
                   paste0(chrName, collapse="_"), ".tsv"),
            quote=F, sep="\t", col.names=T, row.names=F)
write.table(aln_best_pair_hom_maxDist_alenTOqlen_DF,
            paste0(outdir, readsPrefix,
                   "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                   "_hom_maxDist_aTOq", alenTOqlen, "_", recombType,
                   "_alnTo_", alnTo, "_",
                   paste0(chrName, collapse="_"), ".tsv"),
            quote=F, sep="\t", col.names=T, row.names=F)


#aln_best_pair_DF = read.table(paste0(outdir, readsPrefix,
#                                     "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
#                                     "_", recombType,
#                                     "_alnTo_", alnTo, "_",
#                                     paste0(chrName, collapse="_"), ".tsv"),
#                              header=T)
#aln_best_pair_hom_DF = read.table(paste0(outdir, readsPrefix,
#                                         "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
#                                         "_hom_", recombType,
#                                         "_alnTo_", alnTo, "_",
#                                         paste0(chrName, collapse="_"), ".tsv"),
#                                  header=T)
aln_best_pair_hom_maxDist_DF = read.table(paste0(outdir, readsPrefix,
                                                 "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                                                 "_hom_maxDist_", recombType,
                                                 "_alnTo_", alnTo, "_",
                                                 paste0(chrName, collapse="_"), ".tsv"),
                                          header=T)
aln_best_pair_hom_maxDist_alenTOqlen_DF = read.table(paste0(outdir, readsPrefix,
                                                            "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                                                            "_hom_maxDist_aTOq", alenTOqlen, "_", recombType,
                                                            "_alnTo_", alnTo, "_",
                                                            paste0(chrName, collapse="_"), ".tsv"),
                                                     header=T)
