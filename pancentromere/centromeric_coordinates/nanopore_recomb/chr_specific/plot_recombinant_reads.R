#!/usr/bin/env Rscript

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 24/10/2022

# Extract putative recombination events detected in Col-0/Ler-0
# hybrid ONT reads identified on the basis that they contain >=
# a threshold number of accession-specific, chromosome-specific k-mers

# Usage:
# conda activate python_3.9.6
# ./get_chr_specific_segment_pairs_alnToSame.R Col_Ler_F1_pollen_500bp_minq99 Col-0.ragtag_scaffolds Ler-0_110x.ragtag_scaffolds Col-0.ragtag_scaffolds_Chr 24 0.9 10 30000 0.90 co 'Chr1,Chr2,Chr3,Chr4,Chr5'
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
#library(seqinr)


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


make_haplo_mat_list = function(aln_pair_DF) {
    #aln_pair_DF = aln_best_pair_hom_maxDist_DF
    mat_list = lapply(1:nrow(aln_pair_DF), function(x) {
        read_id = aln_pair_DF[x, ]$qname
        read_chr = aln_pair_DF[x, ]$acc1_tname
        read_len = aln_pair_DF[x, ]$read_len
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
        matrix(read_haplo_DF$acc,
               ncol=nrow(read_haplo_DF))
    })

    return(mat_list)
}
     
haplo_mat_list = make_haplo_mat_list(aln_pair_DF=aln_best_pair_hom_maxDist_DF) 





# Plot a haplotype heat map for matches to each COGC haplotype
#for(x in seq_along(ABAmat_list)) {
for(x in 1:13) {
  ABAhtmp <- Heatmap(ABAmat_list[[x]][ ,1:(dim(tplpHapVar)[2]-1)],
                     name = "Allele",
                     col = c("A" = "red", "B" = "blue",
                             "X" = "goldenrod1", "I" = "green2", "N" = "black", "." = "grey60"),                                              
                     row_split = factor(data.frame(ABAmat_list[[x]], stringsAsFactors = F)$hap,                                               
                                        levels = unique(as.character(data.frame(ABAmat_list[[x]], stringsAsFactors = F)$hap))),               
                     row_gap = unit(1.5, "mm"),
                     row_title = paste0(sprintf('%2.0f', ABAmat_list_frequencies[[x]])),                                                      
                     row_title_rot = 0,
                     row_title_gp = gpar(fontsize = 10),
                     row_order = ABAmat_list_row_order[[x]],
                     show_row_names = F,
                     #column_split = colnames(mat1[ ,1:(dim(tplpHap_quant)[2]-1)]),                                                           
                     #column_gap = unit(1.0, "mm"),
                     #column_title = rep("", length(colnames(mat1[ ,1:(dim(tplpHap_quant)[2]-1)]))),                                          
                     column_title = paste0(ABAmat_list_patterns[x],
                                           " haplotype matches in ", sample, " ONT (ratio = ",                                                
                                           round(ABAmat_list_ratios[x], digits = 2), ")"),                                                    
                     column_title_gp = gpar(fontsize = 20, fontface = "bold"),                                                                
                     show_column_names = T,
                     column_names_side = "bottom",
                     column_names_gp = gpar(fontsize = 16),
                     heatmap_legend_param = list(title = bquote(bolditalic("3a") ~ bold("marker allele")),                                    
                                                 title_gp = gpar(fontsize = 16),                                                              
                                                 title_position = "topcenter",                                                                
                                                 grid_height = unit(6, "mm"),                                                                 
                                                 grid_width = unit(10, "mm"),                                                                 
                                                 at = c("A", "B", "X", "."),                                                                  
#                                                        "X", "I", "N", "."),                                                                 
                                                 labels = c("Col-0", "Ws-4", "Nonparental SNV", "Any"),                                       
#                                                            "Nonparental SNV", "Nonparental indel", "Missing", "Any"),                       
                                                 labels_gp = gpar(fontsize = 16),                                                             
                                                 ncol = 4, by_row = T),
                     raster_device = "CairoPNG"
                    )
  pdf(paste0(plotDir, sample, "_ONT_ABA_patterns_", ABAmat_list_patterns[x], "_matches_recombo_heatmap_v140120.pdf"), height = 18, width = 10)
  draw(ABAhtmp,
       heatmap_legend_side = "bottom")
  dev.off()
}

