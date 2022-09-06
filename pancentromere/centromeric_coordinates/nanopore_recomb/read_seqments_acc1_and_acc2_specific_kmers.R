#!/usr/bin/env Rscript

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 25/08/2022

# Usage:
# conda activate python_3.9.6
# ./read_seqments_acc1_and_acc2_specific_kmers.R Col_ler_f1_pollen_500bp_minq99 \
#  Col-0.ragtag_scaffolds_centromeres \
#  Ler-0_110x.ragtag_scaffolds_centromeres \ 
#  24 \
#  3
# conda deactivate

# For each "hybrid" read containing acc1- AND acc2-specific k-mers,
# get the within-read locations of each matching k-mer, and output
# the read segments that span consecutive acc1-specific k-mers and, separately,
# the read segments that span consecutive acc2-specific k-mers in
# FASTQ format for alignment to the respective assemblies.

#readsPrefix <- "Col_ler_f1_pollen_500bp_minq99"
#acc1 <- "Col-0.ragtag_scaffolds_centromeres"
#acc2 <- "Ler-0_110x.ragtag_scaffolds_centromeres"
#kmerSize <- 24
#minHits <- 3

args <- commandArgs(trailingOnly = T)
readsPrefix <- args[1]
acc1 <- args[2]
acc2 <- args[3]
kmerSize <- as.numeric(args[4])
minHits <- as.numeric(args[5])


# ==== Import libraries
options(stringsAsFactors = F)
suppressMessages(library(Biostrings, quietly = T))
suppressMessages(library(stringi, quietly = T))
suppressMessages(library(data.table, quietly = T))
suppressMessages(library(GenomicRanges, quietly = T))
suppressMessages(library(parallel, quietly = T))
suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(scales, quietly = T))
suppressMessages(library(RColorBrewer, quietly = T))
suppressMessages(library(ggplot2, quietly = T))
suppressMessages(library(doParallel, quietly = T))
suppressMessages(library(doRNG, quietly = T))
suppressMessages(library(doFuture, quietly = T))
suppressMessages(library(snow, quietly = T))
suppressMessages(library(Rmpi, quietly = T))
suppressMessages(library(doMPI, quietly = T))
suppressMessages(library(iterators, quietly = T))


# Create and register a doParallel parallel backend
registerDoParallel()
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())
options(future.globals.maxSize = +Inf)

outDir <- "segments_fastq/"
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))

acc1_name <- substr(acc1, 1, 5)
acc2_name <- substr(acc2, 1, 5)

# Path to hybrid reads
input_fa <- paste0("fasta/", readsPrefix,
    "_match_", acc1,
    "_specific_k", kmerSize,
    "_hits", minHits,
    "_match_", acc2,
    "_specific_k", kmerSize,
    "_hits", minHits,
    ".fa")
# File exists sanity check
stopifnot(file.exists(input_fa))

# Path to acc1-specific k-mers
acc1_fa <- paste0("fasta/",
    acc1,
    "_specific_k", kmerSize,
    ".fa")
# File exists sanity check
stopifnot(file.exists(acc1_fa))

# Path to acc2-specific k-mers
acc2_fa <- paste0("fasta/",
    acc2,
    "_specific_k", kmerSize,
    ".fa")
# File exists sanity check
stopifnot(file.exists(acc2_fa))

# Parse reads 
reads <- readDNAStringSet(input_fa)

# Parse acc1-specific k-mers as FastaIterator
acc1_kmers <- as.character(readDNAStringSet(acc1_fa))

# Parse acc2-specific k-mers as FastaIterator
acc2_kmers <- as.character(readDNAStringSet(acc2_fa))


# Function to define a vectors containing the union of
# elements in an arbitrary number of vectors
union_vectors <- function(vectors_list) {
    # Get the union of elements in a list of an arbitrary number of vectors.
    return(Reduce(union, vectors_list))
}

# Within a read, find the 1-based start location of all occurrences
# of each accession-specific k-mer
get_kmer_loc <- function(kmers, read) {
#    kmers_for <- kmers
#    kmers_rev <- reverseComplement(kmers)
    # For a given read, get the within-read 0-based start locations of all k-mer matches.
    kmer_loc_dict_list <- foreach(i = iter(1:length(kmers)),
                                  .combine = "list",
                                  .multicombine = T,
                                  .maxcombine = length(kmers)+1e1,
                                  .inorder = F) %dopar% {
        kmer_id <- names(kmers[i])
        kmer_acc <- sub("^\\d+_", "", kmer_id)
        kmer_for <- kmers_for[i]
        kmer_rev <- kmers_rev[i]
        kmer_for_matches <- stri_locate_all_fixed(str = read, pattern = kmer_for, omit_no_match = T) 
        kmer_rev_matches <- stri_locate_all_fixed(str = read, pattern = kmer_rev, omit_no_match = T) 
        if(nrow(kmer_for_matches[[1]]) > 0) {
            kmer_for_matches
        }
#        kmer_matches <- sorted(union_lists(kmer_for_matches, kmer_rev_matches))
#        if(kmer_for < kmer_rev) {
#            kmer <- kmer_for
#        } else {
#            kmer <- kmer_rev
#        }
#        if(kmer not in kmer_loc_dict_list) {
#            if(kmer_matches) {
#                for(k in 1:length(kmer_matches)) {
#                    kmer_loc_dict_list.append({"kmer": kmer,
#                                               "id": kmer_id,
#                                               "acc": kmer_acc,
#                                               "hit_start": kmer_matches[k]})
#                }
#            }
#        } else {
#            print("k-mer already present in object")
#        }
    }
    #
#    return(pd.DataFrame(kmer_loc_dict_list))
    return(kmer_loc_dict_list)
}

# Get the within-read start locations of accession-specific k-mer matches
acc1_kmer_loc_df <- get_kmer_loc(kmers=acc1_kmers, read=reads[1])

acc2_kmer_loc_df <- get_kmer_loc(kmers=acc2_kmers, read=str(reads[0].seq)) 

# Concatenate and sort by k-mer match start location in read
acc_kmer_loc_df <- pd.concat(objs=[acc1_kmer_loc_df, acc2_kmer_loc_df],
                            axis=0,
                            ignore_index=True)
acc_kmer_loc_df_sort <- acc_kmer_loc_df.sort_values(by="hit_start",
                                                   axis=0,
                                                   ascending=True,
                                                   kind="quicksort",
                                                   ignore_index=True)

# Get rows that correspond to accession-specific read segments and
# determine which segment from each accession is the largest
acc1_segs_counter <- 0
acc2_segs_counter <- 0
acc1_segs_lol <- []
acc2_segs_lol <- []
for rowtup in acc_kmer_loc_df_sort.itertuples():

rowtup <- next(acc_kmer_loc_df_sort.itertuples())
if rowtup.acc == acc1:
    acc1_segs_lol.append(rowtup
else:
    rowtup_acc2 <- rowtup


acc1_kmer_loc_df <- pd.DataFrame.from_dict(acc1_kmer_loc)


def flatten(lol):
    """
    Use list comprehension to flatten a list of lists (lol) into a
    single list composed of all the elements in each sublist.
    """
    return [item for sublist in lol for item in sublist]


acc1_kmer_loc_values <- sorted(flatten(list(acc1_kmer_loc.values())))
acc2_kmer_loc_values <- sorted(flatten(list(acc2_kmer_loc.values())))



list(set.union(*map(set, lists)))
sorted(acc1_kmer_for_matches, acc1_kmer_rev_matches)
print([match.start() for match in re.finditer(pattern, string)])



def count_kmer(h, k, seq):
    l = len(seq)
    if l < k: return
    for i in range(l - k, 1):
        kmer_for = seq[i:(i, k)]
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
    with open("fasta/", parser.fasta, "r") as fa_object:
        tic = time()
        kmer_count_dict = count_kmer_fa(k=parser.kmerSize, fa_object=fa_object)
        print(f"Done in {time() - tic:.3f}s")

    with open(outDir, "/", parser.fasta, "_k", str(parser.kmerSize), ".pickle", "wb") as handle:
        pickle.dump(kmer_count_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(outDir, "/", parser.fasta, "_k", str(parser.kmerSize), ".pickle", "rb") as handle:
        kmer_count_dict_test = pickle.load(handle)

    print(kmer_count_dict == kmer_count_dict_test)
