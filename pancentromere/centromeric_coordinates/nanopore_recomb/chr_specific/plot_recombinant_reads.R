#!/usr/bin/env Rscript

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 24/10/2022

# Extract putative recombination events detected in Col-0/Ler-0
# hybrid ONT reads identified on the basis that they contain >=
# a threshold number of accession-specific, chromosome-specific k-mers

# Usage:
# conda activate python_3.9.6
# ./plot_recombinant_reads.R ColLerF1pollen_1000bp_minq90 Col-0.ragtag_scaffolds Ler-0_110x.ragtag_scaffolds not_centromere Col-0.ragtag_scaffolds_Chr 24 0.9 11 30000 0.90 co 'Chr1,Chr2,Chr3,Chr4,Chr5'
# conda deactivate

#readsPrefix = "ColLerF1pollen_1000bp_minq90"
#acc1 = "Col-0.ragtag_scaffolds"
#acc2 = "Ler-0_110x.ragtag_scaffolds"
#region = "not_centromere"
#alnTo = "Col-0.ragtag_scaffolds_Chr"
#kmerSize = 24
#overlapProp = 0.9
#minHits = 11
#maxDist = 30000
#alenTOqlen = 0.90
#recombType = "co"
#chrName = unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split=","))

args = commandArgs(trailingOnly=T)
readsPrefix = args[1]
acc1 = args[2]
acc2 = args[3]
region = args[4]
alnTo = args[5]
kmerSize = as.integer(args[6])
overlapProp = as.numeric(args[7])
minHits = as.integer(args[8])
maxDist = as.integer(args[9])
alenTOqlen = as.numeric(args[10])
recombType = args[11]
chrName = unlist(strsplit(args[12],
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
library(ComplexHeatmap)
library(Cairo)
library(seqinr)
library(dplyr)
library(parallel)
library(doParallel)
# Create and register a doParallel parallel backend
registerDoParallel()
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())


outDir = paste0(region, "/segment_pairs/", recombType, "/")
plotDir = paste0(outDir, "recombinant_read_plots/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# Accession names
acc1_name = strsplit( strsplit(acc1, split="\\.")[[1]][1],
                      split="_" )[[1]][1]
acc2_name = strsplit( strsplit(acc2, split="\\.")[[1]][1],
                      split="_" )[[1]][1]


aln_best_pair_DF = dplyr::bind_rows(
    mclapply(1:length(acc1_aln_chr_list_of_lists), function(x) {
        aln_best_pair(acc1_aln_DF_list=acc1_aln_chr_list_of_lists[[x]], acc2_aln_DF_list=acc2_aln_chr_list_of_lists[[x]])
    }, mc.preschedule=F, mc.cores=length(acc1_aln_chr_list_of_lists))
)


aln_best_pair_hom_maxDist_DF_list = lapply(1:length(chrName), function(x) {
    read.table(paste0(outDir, readsPrefix,
                      "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                      "_hom_maxDist_", recombType,
                      "_alnTo_", alnTo, "_",
                      chrName[x], ".tsv"),
               header=T)
})
if(length(aln_best_pair_hom_maxDist_DF_list) > 1) {
    aln_best_pair_hom_maxDist_DF = dplyr::bind_rows(aln_best_pair_hom_maxDist_DF_list)
    } else {
    aln_best_pair_hom_maxDist_DF = aln_best_pair_hom_maxDist_DF[[1]]
}


aln_best_pair_hom_maxDist_alenTOqlen_DF_list = lapply(1:length(chrName), function(x) {
    read.table(paste0(outDir, readsPrefix,
                      "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                      "_hom_maxDist_aTOq", alenTOqlen, "_", recombType,
                      "_alnTo_", alnTo, "_",
                      chrName[x], ".tsv"),
               header=T)
})
if(length(aln_best_pair_hom_maxDist_alenTOqlen_DF_list) > 1) {
    aln_best_pair_hom_maxDist_alenTOqlen_DF = dplyr::bind_rows(aln_best_pair_hom_maxDist_alenTOqlen_DF_list)
    } else {
    aln_best_pair_hom_maxDist_alenTOqlen_DF = aln_best_pair_hom_maxDist_alenTOqlen_DF[[1]]
}


aln_DF = aln_best_pair_hom_maxDist_alenTOqlen_DF

# Create list of all hybrid reads
hybrid_reads_list = list()
for(x in 1:length(chrName)) {
    hybrid_reads_list_chr = read.fasta(paste0(outDir, readsPrefix,
                                              "_", acc1, "_", acc2,
                                              "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                                              "_hom_maxDist_aTOq", alenTOqlen, "_", recombType,
                                              "_alnTo_", alnTo, "_",
                                              chrName[x], "_hybrid_reads.fa"),
                                       forceDNAtolower=F)
    hybrid_reads_list = c(hybrid_reads_list, hybrid_reads_list_chr)
}

# Make haplotype data.frames of reads based on k-mer locations
make_haplo_DF_list = function(aln_pair_DF, hybrid_reads_list) {
    DF_list = mclapply(1:nrow(aln_pair_DF), function(x) {
        print(x)
        read_id = aln_pair_DF[x, ]$qname
        read_chr = aln_pair_DF[x, ]$acc1_tname
        read_rfa = hybrid_reads_list[[which(names(hybrid_reads_list) == read_id)]]
        read_seq = getSequence(read_rfa)
        read_len = length(read_seq)
        read_kmer_tsv = system(paste0("ls ", region, "/", read_chr, "/kmer_loc_tsv/",
                                      read_id, "__*", "_alnTo_", alnTo, "_kmer_loc.tsv"), intern=T)
        if(length(read_kmer_tsv) > 0) {
            read_kmer_loc = read.table(read_kmer_tsv, header=T)
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
        }
    }, mc.cores=detectCores(), mc.preschedule=T)

    return(DF_list)
}

haplo_DF_list = make_haplo_DF_list(aln_pair_DF = aln_DF,
                                   hybrid_reads_list = hybrid_reads_list) 

# Plot haplotype of recombinant read
haplo_heatmap = function(haplo_DF, aln_DF) {
  haplo_mat = t(haplo_DF[ , 2:3])
  colnames(haplo_mat) = 1:ncol(haplo_mat)
  Heatmap(
          haplo_mat,
          col = c("Col-0"="dodgerblue3", "Ler-0"="darkgoldenrod", "X"="grey70",
                  "A"="firebrick3", "T"="forestgreen", "G"="darkgoldenrod1", "C"="blue3"),
          show_row_names = F,
          column_title = paste0(aln_DF$acc1_tname[x], " read ID: ", aln_DF$qname[x]),
          column_title_gp = gpar(fontsize = 16, fontface = "bold"),
          column_labels = colnames(haplo_mat),
          column_names_rot = 90,
          column_names_centered = F,
          show_column_names = T,
          column_names_side = "top",
          column_names_gp = gpar(fontsize = 4),
          heatmap_legend_param = list(title = bquote(bold("Accession-specific, chromosome-specific" ~ .(as.character(kmerSize)) * "-mer")),
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
          #height = unit(2, "npc"),
          #width = unit(100, "npc"),
          rect_gp = gpar(col = "white", lwd = 0.01),
          raster_device = "CairoPNG"
         )
}

foreach(x = 1:length(haplo_DF_list)) %dopar% {
    if(!is.null(haplo_DF_list[[x]][[1]])) {
        haplo_htmp = haplo_heatmap(haplo_DF = haplo_DF_list[[x]], aln_DF=aln_DF)
        pdf(paste0(plotDir, aln_DF$acc1_tname[x], "_", aln_DF$qname[x],
                   "_alnTo_", alnTo, "_haplo_heatmap.pdf"),
            height = 2, width = 0.05 * nrow(haplo_DF_list[[x]]))
        draw(haplo_htmp,
             heatmap_legend_side = "bottom")
        dev.off()
    }
}
