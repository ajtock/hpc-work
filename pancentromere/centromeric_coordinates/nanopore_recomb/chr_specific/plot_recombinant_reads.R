#!/usr/bin/env Rscript

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 24/10/2022

# Extract putative recombination events detected in Col-0/Ler-0
# hybrid ONT reads identified on the basis that they contain >=
# a threshold number of accession-specific, chromosome-specific k-mers

# Usage:
# conda activate python_3.9.6
# ./plot_recombinant_reads.R Col_Ler_F1_pollen_500bp_minq99 Col-0.ragtag_scaffolds Ler-0_110x.ragtag_scaffolds Col-0.ragtag_scaffolds_Chr 24 0.9 10 30000 0.90 co 'Chr1,Chr2,Chr3,Chr4,Chr5'
# conda deactivate

#readsPrefix = "Col_Ler_F1_pollen_500bp_minq99"
#acc1 = "Col-0.ragtag_scaffolds"
#acc2 = "Ler-0_110x.ragtag_scaffolds"
#alnTo = "Col-0.ragtag_scaffolds_Chr"
#kmerSize = 24
#overlapProp = 0.9
#minHits = 10
#maxDist = 30000
#alenTOqlen = 0.90
#recombType = "co"
#chrName = unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split=","))

args = commandArgs(trailingOnly=T)
readsPrefix = args[1]
acc1 = args[2]
acc2 = args[3]
alnTo = args[4]
kmerSize = as.integer(args[5])
overlapProp = as.numeric(args[6])
minHits = as.integer(args[7])
maxDist = as.integer(args[8])
alenTOqlen = as.numeric(args[9])
recombType = args[10]
chrName = unlist(strsplit(args[11],
                          split=","))


if(floor(log10(maxDist)) + 1 < 4) {
    maxDistName = paste0(maxDist, "bp")
} else if(floor(log10(maxDist)) + 1 >= 4 &
          floor(log10(maxDist)) + 1 <= 6) {
    maxDistName = paste0(maxDist/1e3, "kb")
} else if(floor(log10(maxDist)) + 1 >= 7) {
    maxDistName = paste0(maxDist/1e6, "Mb")
}

options(stringsAsFactors=F)
options(scipen=999)
#library(parallel)
#library(GenomicRanges)
#library(circlize)
library(ComplexHeatmap)
library(Cairo)
#library(gridBase)
#library(viridis)
#library(colorspace)
#library(data.table)
#library(dplyr)
#library(tidyr)
#library(ggplot2)
library(seqinr)


outDir = paste0("not_centromere/segment_pairs/", recombType, "/")
plotDir = paste0(outDir, "recombinant_read_plots/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# Accession names
acc1_name = strsplit( strsplit(acc1, split="\\.")[[1]][1],
                      split="_" )[[1]][1]
acc2_name = strsplit( strsplit(acc2, split="\\.")[[1]][1],
                      split="_" )[[1]][1]


aln_best_pair_hom_maxDist_DF = read.table(paste0(outDir, readsPrefix,
                                                  "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                                                  "_hom_maxDist_", recombType,
                                                  "_alnTo_", alnTo, "_",
                                                  paste0(chrName, collapse="_"), ".tsv"),
                                           header=T)
#aln_best_pair_hom_maxDist_alenTOqlen_DF = read.table(paste0(outDir, readsPrefix,
#                                                             "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
#                                                             "_hom_maxDist_aTOq", alenTOqlen, "_", recombType,
#                                                             "_alnTo_", alnTo, "_",
#                                                             paste0(chrName, collapse="_"), ".tsv"),
#                                                      header=T)

# Create list of all hybrid reads
hybrid_reads_list = list()
for(x in 1:length(chrName)) {
    hybrid_reads_list_chr = read.fasta(paste0("fasta/", readsPrefix,
                                              "_match_", acc1, "_not_centromere_", chrName[x],
                                              "_specific_k", kmerSize, "_downsampled_op", overlapProp, "_hits", minHits,
                                              "_match_", acc2, "_not_centromere_", chrName[x],
                                              "_specific_k", kmerSize, "_downsampled_op", overlapProp, "_hits", minHits,
                                              ".fa"),
                                       forceDNAtolower=F)
    hybrid_reads_list = c(hybrid_reads_list, hybrid_reads_list_chr)
}

# Make haplotype data.frames of reads based on k-mer locations
make_haplo_DF_list = function(aln_pair_DF, hybrid_reads_list) {
    #aln_pair_DF = aln_best_pair_hom_maxDist_DF
    #hybrid_reads_list = hybrid_reads_list
    DF_list = lapply(1:nrow(aln_pair_DF), function(x) {
        print(x)
        read_id = aln_pair_DF[x, ]$qname
        read_chr = aln_pair_DF[x, ]$acc1_tname
        read_rfa = hybrid_reads_list[[which(names(hybrid_reads_list) == read_id)]]
        read_seq = getSequence(read_rfa)
        read_len = length(read_seq)
        read_kmer_loc = read.table(paste0("not_centromere/", read_chr, "/kmer_loc_tsv/",
                                          read_id, "__kmer_loc.tsv"), header=T)
        read_haplo_DF = data.frame()
        for(h in 1:read_len) {
            if(h %in% (read_kmer_loc$hit_start + 1)) {
                pos_kmer_loc = read_kmer_loc[which(read_kmer_loc$hit_start + 1 == h),]
                pos_geno_DF = data.frame(pos = h,
                                         acc = pos_kmer_loc$acc)
            } else {
                pos_geno_DF = data.frame(pos = h,
                                         acc = "X")
            }
            read_haplo_DF = rbind(read_haplo_DF, pos_geno_DF)
        }
        read_haplo_DF = data.frame(pos = read_haplo_DF$pos,
                                   seq = read_seq,
                                   acc = read_haplo_DF$acc)
        read_haplo_DF
    })

    return(DF_list)
}

haplo_DF_list = make_haplo_DF_list(aln_pair_DF=aln_best_pair_hom_maxDist_DF,
                                   hybrid_reads_list=hybrid_reads_list) 

# Plot a haplotype heat map for each read
for(x in 1:length(haplo_DF_list)) {
  htmp = Heatmap(t(haplo_DF_list[[x]][ , 2:3]),
                 name = "Accession",
                 col = c("Col-0"="dodgerblue3", "Ler-0"="darkgoldenrod", "X"="grey60",
                         "A"="firebrick3", "T"="forestgreen", "G"="darkgoldenrod1", "C"="blue3"),
                 show_row_names = F,
                 column_title = paste0(aln_best_pair_hom_maxDist_DF$acc1_tname[x], " read ID: ", aln_best_pair_hom_maxDist_DF$qname[x]),
                 column_title_gp = gpar(fontsize = 20, fontface = "bold"),                                                                
                 show_column_names = T,
                 column_names_side = "bottom",
                 column_names_gp = gpar(fontsize = 16),
                 heatmap_legend_param = list(title = bquote(bold("Accession-specific, chromosome-specific") ~ bold(.(kmerSize)) * bold("-mer")),
                                             title_gp = gpar(fontsize = 16),
                                             title_position = "topleft",
                                             grid_height = unit(6, "mm"),
                                             grid_width = unit(10, "mm"),
                                             at = c("Col-0", "Ler-0", "X",
                                                    "A", "T", "G", "C"),
                                             labels = c("Col-0", "Ler-0", "Nonspecific",
                                                        "A", "T", "G", "C"),
                                             labels_gp = gpar(fontsize = 16),
                                             ncol = 7, by_row = T),
                 raster_device = "CairoPNG"
                )
  pdf(paste0(plotDir, aln_best_pair_hom_maxDist_DF$qname[x], "__", aln_best_pair_hom_maxDist_DF$acc1_tname[x],
             "_alnTo_", alnTo, "_haplo_heatmap.pdf"),
      height = 2, width = 100)
  draw(htmp,
       heatmap_legend_side = "bottom")
  dev.off()
}

