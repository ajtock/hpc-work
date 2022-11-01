#!/usr/bin/env Rscript

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 09/09/2022

# Plot putative crossover events detected in Col-0/Ler-0
# hybrid ONT reads


# Usage:
# conda activate python_3.9.6
# ./circlize_co_nco.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 10000 Col-0.ragtag_scaffolds Ler-0_110x.ragtag_scaffolds co 201022
# conda deactivate

#chrName = unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split=","))
#genomeBinSize = 10000
#acc1 = "Col-0.ragtag_scaffolds"
#acc2 = "Ler-0_110x.ragtag_scaffolds"
#recombType = "nco"
#date = "201022"

args = commandArgs(trailingOnly=T)
chrName = unlist(strsplit(args[1],
                           split=","))
genomeBinSize = as.integer(args[2])
acc1 = args[3]
acc2 = args[4]
recombType = args[5]
date = as.character(args[6])

if(floor(log10(genomeBinSize)) + 1 < 4) {
    genomeBinName = paste0(genomeBinSize, "bp")
} else if(floor(log10(genomeBinSize)) + 1 >= 4 &
          floor(log10(genomeBinSize)) + 1 <= 6) {
    genomeBinName = paste0(genomeBinSize/1e3, "kb")
} else if(floor(log10(genomeBinSize)) + 1 >= 7) {
    genomeBinName = paste0(genomeBinSize/1e6, "Mb")
}

options(stringsAsFactors=F)
options(scipen=999)
library(parallel)
library(GenomicRanges)
library(circlize)
library(ComplexHeatmap)
library(gridBase)
library(viridis)
library(colorspace)
library(data.table)
library(dplyr)

plotDir = paste0("plots/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# Accession names
acc1_name = strsplit( strsplit(acc1, split="\\.")[[1]][1],
                      split="_" )[[1]][1]
acc2_name = strsplit( strsplit(acc2, split="\\.")[[1]][1],
                      split="_" )[[1]][1]

# Directories containing read segment alignment files
acc1_indir = paste0("segments/", acc1_name, "/", recombType, "/")
acc2_indir = paste0("segments/", acc2_name, "/", recombType, "/")

# CEN coordinates
CEN = read.csv(paste0("/home/ajt200/rds/hpc-work/pancentromere/centromeric_coordinates/",
                      "centromere_manual_EDTA4_fa.csv"),
               header=T)
CEN$fasta.name = gsub(".fa", "", CEN$fasta.name)

# Genomic definitions

#acc1
acc1_fai = read.table(paste0("index/", acc1, ".fa.fai"), header=F)
acc1_chrs = acc1_fai[which(acc1_fai$V1 %in% chrName),]$V1
acc1_chrLens = acc1_fai[which(acc1_fai$V1 %in% chrName),]$V2

acc1_CEN = CEN[grep(acc1, CEN$fasta.name),]
acc1_CEN = acc1_CEN[,which(colnames(acc1_CEN) %in% c("chr", "start", "end"))]
acc1_CEN_new = data.frame()
for(i in 1:length(acc1_chrs)) {
    acc1_CEN_chr = acc1_CEN[which(acc1_CEN$chr == acc1_chrs[i]),]
    if(nrow(acc1_CEN_chr) > 1) {
        acc1_CEN_chr = data.frame(chr=acc1_CEN_chr$chr[1],
                                  start=acc1_CEN_chr$start[1],
                                  end=acc1_CEN_chr$end[nrow(acc1_CEN_chr)])
    }
    acc1_CEN_new = rbind(acc1_CEN_new, acc1_CEN_chr)
}
acc1_CEN = acc1_CEN_new
acc1_CENstart = acc1_CEN$start
acc1_CENend = acc1_CEN$end
acc1_chrs = paste0(acc1_name, "_", acc1_chrs)

#acc2
acc2_fai = read.table(paste0("index/", acc2, ".fa.fai"), header=F)
acc2_chrs = acc2_fai[which(acc2_fai$V1 %in% chrName),]$V1
acc2_chrLens = acc2_fai[which(acc2_fai$V1 %in% chrName),]$V2

acc2_CEN = CEN[grep(acc2, CEN$fasta.name),]
acc2_CEN = acc2_CEN[,which(colnames(acc2_CEN) %in% c("chr", "start", "end"))]
acc2_CEN_new = data.frame()
for(i in 1:length(acc2_chrs)) {
    acc2_CEN_chr = acc2_CEN[which(acc2_CEN$chr == acc2_chrs[i]),]
    if(nrow(acc2_CEN_chr) > 1) {
        acc2_CEN_chr = data.frame(chr=acc2_CEN_chr$chr[1],
                                  start=acc2_CEN_chr$start[1],
                                  end=acc2_CEN_chr$end[nrow(acc2_CEN_chr)])
    }
    acc2_CEN_new = rbind(acc2_CEN_new, acc2_CEN_chr)
}
acc2_CEN = acc2_CEN_new
acc2_CENstart = acc2_CEN$start
acc2_CENend = acc2_CEN$end
acc2_chrs = paste0(acc2_name, "_", acc2_chrs)


# Load read segment alignment files as a combined data.frame
load_pafs = function(indir, acc_name, suffix, aligner) {
    #indir=acc1_indir
    #acc_name=acc1_name
    #suffix="_wm_ont.paf"
    #aligner="wm"
    files = system(paste0("ls -1 ", indir, "*", acc_name, suffix), intern=T)
    aln_DF = data.frame()
    for(h in 1:length(files)) {
        aln = fread(files[h],
                    header=F, fill=T, sep="\t", data.table=F)[,1:13]
        aln_DF = rbind(aln_DF, aln)
    }
    aln_DF$aligner = aligner
    colnames(aln_DF) = c("qname", "qlen", "qstart0", "qend0",
                         "strand", "tname", "tlen", "tstart", "tend",
                         "nmatch", "alen", "mapq", "atype", "aligner")
 
    return(aln_DF)
}


# wm alignments
acc1_wm = load_pafs(indir=acc1_indir, acc_name=acc1_name, suffix="_wm_ont.paf", aligner="wm")
acc2_wm = load_pafs(indir=acc2_indir, acc_name=acc2_name, suffix="_wm_ont.paf", aligner="wm")

# mm alignments
acc1_mm = load_pafs(indir=acc1_indir, acc_name=acc1_name, suffix="_mm_ont.paf", aligner="mm")
acc2_mm = load_pafs(indir=acc2_indir, acc_name=acc2_name, suffix="_mm_ont.paf", aligner="mm")

# sr alignments
acc1_sr = load_pafs(indir=acc1_indir, acc_name=acc1_name, suffix="_mm_sr.paf", aligner="sr")
acc2_sr = load_pafs(indir=acc2_indir, acc_name=acc2_name, suffix="_mm_sr.paf", aligner="sr")

# acc alignments list
acc1_aln_list = list(acc1_wm, acc1_mm, acc1_sr)
acc2_aln_list = list(acc2_wm, acc2_mm, acc2_sr)
names(acc1_aln_list) = c("wm", "mm", "sr")
names(acc2_aln_list) = c("wm", "mm", "sr")

# Get best pair of acc1 and acc2 read segment alignments, based on:
# 1. The alignment chromosome
# 2. The aligner:
#   Prioritise:
#     1. acc1_wm : acc2_wm, 2. acc1_wm : acc2_mm, 3. acc1_mm : acc2_wm,
#     4. acc1_wm : acc2_sr, 5. acc1_sr : acc2_wm, 6. acc1_mm : acc2_mm,
#     7. acc1_mm : acc2_sr, 8. acc1_sr : acc2_mm, 9. acc1_sr : acc2_sr
# 3. The alignment type (atype)
# 4. The alignment mapq score (mapq)
# 5. The alignment length (alen)
# 6. The alignment number of matching bases (nmatch)
aln_best_pair = function(acc1_aln_DF_list, acc2_aln_DF_list) {
    #acc1_aln_DF_list=acc1_aln_list
    #acc2_aln_DF_list=acc2_aln_list
    # Each of the 3 list elements in acc1_aln_DF_list is a
    # a data.frame of alignments done by wm_ont, mm_ont or mm_sr.
    aligner_list = lapply(1:length(acc1_aln_DF_list), function(l) {

        acc1_aln_DF_list_l = acc1_aln_DF_list[[l]]

        # For each read segment alignment in data.frame acc1_aln_DF_list_l,
        # get the best read segment alignment from acc2 corresponding to
        # the same read
        # For many read IDs, this will give more than 1 pair of
        # read segment alignments to subsequently select 1 from in the next loop
        acc1_aln_DF = data.frame()
        acc2_aln_DF = data.frame()
        for(m in 1:nrow(acc1_aln_DF_list_l)) {

            print(m)
            acc1_aln_m = acc1_aln_DF_list_l[m, ]
            acc1_aln_DF = rbind(acc1_aln_DF, acc1_aln_m)

            qname_match_acc2_aln_DF_list = lapply(1:length(acc2_aln_DF_list), function(x) {

                tmp = acc2_aln_DF_list[[x]] %>%
                    dplyr::filter(qname == acc1_aln_m$qname)
                tmp = tmp[ with(tmp,
                                order(mapq, alen, nmatch, decreasing=T)), ]
                tmp_chr = tmp[which(tmp$tname == acc1_aln_m$tname),]
                if(nrow(tmp_chr) > 0) {
                    if("tp:A:P" %in% tmp_chr$atype) {
                        tmp_select = tmp_chr[ which(tmp_chr$atype == "tp:A:P"), ][1, ]
                    } else {
                        tmp_select = tmp_chr[1, ]
                    }
                } else if(nrow(tmp) > 0) {
                    if("tp:A:P" %in% tmp$atype) {
                        tmp_select = tmp[ which(tmp$atype == "tp:A:P"), ][1, ]
                    } else {
                        tmp_select = tmp[1, ]
                    }
                } else {
                    tmp_select = tmp
                }
 
                tmp_select

            })

            acc2_aln_m_DF = dplyr::bind_rows(qname_match_acc2_aln_DF_list)
            acc2_aln_m_DF = acc2_aln_m_DF[ with(acc2_aln_m_DF,
                                                order(mapq, alen, nmatch, decreasing=T)), ]

            if(acc1_aln_m$tname %in% acc2_aln_m_DF$tname) {
                acc2_aln_m_DF = acc2_aln_m_DF[ which(acc2_aln_m_DF$tname == acc1_aln_m$tname), ]
            }

            if("wm" %in% acc2_aln_m_DF$aligner) {
                acc2_aln_m = acc2_aln_m_DF[ which(acc2_aln_m_DF$aligner == "wm"), ][1, ]
            } else if("mm" %in% acc2_aln_m_DF$aligner) {
                acc2_aln_m = acc2_aln_m_DF[ which(acc2_aln_m_DF$aligner == "mm"), ][1, ]
            } else {
                acc2_aln_m = acc2_aln_m_DF[ which(acc2_aln_m_DF$aligner == "sr"), ][1, ]
            }
 
            acc2_aln_DF = rbind(acc2_aln_DF, acc2_aln_m)

        }

        # Get the best alignment from acc1_aln_DF and corresponding row from acc2_aln_DF
        acc1_aln_DF_best = data.frame()
        acc2_aln_DF_best = data.frame()
        for(read_id in unique(acc1_aln_DF$qname)) {

            acc1_aln_DF_read_id = acc1_aln_DF[ which(acc1_aln_DF$qname == read_id), ]
            acc2_aln_DF_read_id = acc2_aln_DF[ which(acc2_aln_DF$qname == read_id), ]

            acc1_aln_DF_read_id_order =  with(acc1_aln_DF_read_id,
                                              order(mapq, alen, nmatch, decreasing=T))

            acc1_aln_DF_read_id = acc1_aln_DF_read_id[acc1_aln_DF_read_id_order, ]
            acc2_aln_DF_read_id = acc2_aln_DF_read_id[acc1_aln_DF_read_id_order, ]

            tname_match_row_idx = which(acc1_aln_DF_read_id$tname %in%
                                        acc2_aln_DF_read_id$tname)
            if(length(tname_match_row_idx) > 0) {
                acc1_aln_DF_read_id = acc1_aln_DF_read_id[ tname_match_row_idx, ]
                acc2_aln_DF_read_id = acc2_aln_DF_read_id[ tname_match_row_idx, ]
            }

            atype_match_row_idx = which(acc1_aln_DF_read_id$atype == "tp:A:P")
            if(length(atype_match_row_idx) > 0) {
                acc1_aln_DF_read_id = acc1_aln_DF_read_id[ atype_match_row_idx, ][1, ]
                acc2_aln_DF_read_id = acc2_aln_DF_read_id[ atype_match_row_idx, ][1, ]
            } else {
                acc1_aln_DF_read_id = acc1_aln_DF_read_id[1, ]
                acc2_aln_DF_read_id = acc2_aln_DF_read_id[1, ]
            }

            acc1_aln_DF_best = rbind(acc1_aln_DF_best, acc1_aln_DF_read_id)
            acc2_aln_DF_best = rbind(acc2_aln_DF_best, acc2_aln_DF_read_id)

        }

        colnames(acc1_aln_DF_best) = paste0("acc1_", colnames(acc1_aln_DF_best))
        colnames(acc2_aln_DF_best) = paste0("acc2_", colnames(acc2_aln_DF_best))
        stopifnot(identical(acc1_aln_DF_best$acc1_qname, acc2_aln_DF_best$acc2_qname))

        cbind(acc1_aln_DF_best, acc2_aln_DF_best)

    })

    names(aligner_list) = c("wm", "mm", "sr")

    aligner_bind_rows = dplyr::bind_rows(aligner_list, .id="id")
    stopifnot(identical(aligner_bind_rows$id, aligner_bind_rows$acc1_aligner))

    # Up to 3 read segment alignment pairs have been selected for each read
    # (potentially one pair for each of the 3 aligners used to align the acc1 segment)
    # Select the best pair based on the criteria described at the top of the function definition
    aln_best_pair_DF = data.frame()
    for(uniq_qname in unique(aligner_bind_rows$acc1_qname)) {

        tmp = aligner_bind_rows %>%
            dplyr::filter(acc1_qname == uniq_qname)
        tmp = tmp[ with(tmp,
                        order(acc1_mapq, acc2_mapq,
                              acc1_alen, acc2_alen,
                              acc1_nmatch, acc2_nmatch,
                              decreasing=T)), ]

        # Prioritise alignment pairs where both read segments align to the same chromosome
        tmp_chr = tmp[ which(tmp$acc1_tname == tmp$acc2_tname), ]
        if(nrow(tmp_chr) > 0) {
            tmp_chr_wm_wm = tmp_chr[ which(tmp_chr$acc1_aligner == "wm" &
                                           tmp_chr$acc2_aligner == "wm"), ]
            tmp_chr_wm = tmp_chr[ which(tmp_chr$acc1_aligner == "wm" |
                                        tmp_chr$acc2_aligner == "wm"), ]
            tmp_chr_mm_mm = tmp_chr[ which(tmp_chr$acc1_aligner == "mm" &
                                           tmp_chr$acc2_aligner == "mm"), ]
            tmp_chr_mm = tmp_chr[ which(tmp_chr$acc1_aligner == "mm" |
                                        tmp_chr$acc2_aligner == "mm"), ]
            if(nrow(tmp_chr_wm_wm) > 0) {
                tmp = tmp_chr_wm_wm[1, ]
            } else if(nrow(tmp_chr_wm) > 0) {
                tmp = tmp_chr_wm[1, ]
            } else if(nrow(tmp_chr_mm_mm) > 0) {
                tmp = tmp_chr_mm_mm[1, ]
            } else if(nrow(tmp_chr_mm) > 0) {
                tmp = tmp_chr_mm[1, ]
            } else {
                tmp = tmp_chr[1, ]
            }
        } else {
            tmp_wm_wm = tmp[ which(tmp$acc1_aligner == "wm" &
                                   tmp$acc2_aligner == "wm"), ]
            tmp_wm = tmp[ which(tmp$acc1_aligner == "wm" |
                                tmp$acc2_aligner == "wm"), ]
            tmp_mm_mm = tmp[ which(tmp$acc1_aligner == "mm" &
                                   tmp$acc2_aligner == "mm"), ]
            tmp_mm = tmp[ which(tmp$acc1_aligner == "mm" |
                                tmp$acc2_aligner == "mm"), ]
            if(nrow(tmp_wm_wm) > 0) {
                tmp = tmp_wm_wm[1, ]
            } else if(nrow(tmp_wm) > 0) {
                tmp = tmp_wm[1, ]
            } else if(nrow(tmp_mm_mm) > 0) {
                tmp = tmp_mm_mm[1, ]
            } else if(nrow(tmp_mm) > 0) {
                tmp = tmp_mm[1, ]
            } else {
                tmp = tmp[1, ]
            }
        }

        aln_best_pair_DF = rbind(aln_best_pair_DF, tmp)

    }

    return(aln_best_pair_DF)

}


# Get best pair of aligned read segments for each read
aln_best_pair_DF = aln_best_pair(acc1_aln_DF_list=acc1_aln_list, acc2_aln_DF_list=acc2_aln_list)
stopifnot(identical(aln_best_pair_DF$acc1_qname, aln_best_pair_DF$acc2_qname))

#aln_best_pair_DF_ori = aln_best_pair_DF
#
#aln_best_pair_DF = aln_best_pair_DF %>%
#    dplyr::filter(acc1_alen >= 200) %>%
#    dplyr::filter(acc2_alen >= 200)

print(paste0(nrow(aln_best_pair_DF), " putative ", recombType, " events"))
#[1] "213 putative co events"
#[1] "966 putative nco events"

# Filter to retain putative recombination events between homologous chromosomes only
aln_best_pair_hom_DF = aln_best_pair_DF[ which(aln_best_pair_DF$acc1_tname == aln_best_pair_DF$acc2_tname), ]

print(paste0(nrow(aln_best_pair_hom_DF), " putative ", recombType, " events between homologous chromosomes"))
#[1] "135 putative co events between homologous chromosomes"
#[1] "559 putative nco events between homologous chromosomes"

print(paste0( round( ( nrow(aln_best_pair_hom_DF) / nrow(aln_best_pair_DF) ), 2 ) * 100, "% of putative ", recombType, " events between homologous chromosomes"))
#[1] "63% of putative co events between homologous chromosomes"
#[1] "58% of putative nco events between homologous chromosomes"


# circlize

# Define "links" between the two accessions' chromosomes to represent
# recombination events in circlize plot

acc1_all_bed = data.frame(chr = paste0(acc1_name, "_", aln_best_pair_DF$acc1_tname),
                          start = aln_best_pair_DF$acc1_tstart-1,
                          end = aln_best_pair_DF$acc1_tend,
                          aligner = aln_best_pair_DF$acc1_aligner,
                          atype = aln_best_pair_DF$acc1_atype,
                          mapq = aln_best_pair_DF$acc1_mapq,
                          alen = aln_best_pair_DF$acc1_alen,
                          nmatch = aln_best_pair_DF$acc1_nmatch)

acc2_all_bed = data.frame(chr = paste0(acc2_name, "_", aln_best_pair_DF$acc2_tname),
                          start = aln_best_pair_DF$acc2_tstart-1,
                          end = aln_best_pair_DF$acc2_tend,
                          aligner = aln_best_pair_DF$acc2_aligner,
                          atype = aln_best_pair_DF$acc2_atype,
                          mapq = aln_best_pair_DF$acc2_mapq,
                          alen = aln_best_pair_DF$acc2_alen,
                          nmatch = aln_best_pair_DF$acc2_nmatch)

acc1_hom_bed = data.frame(chr = paste0(acc1_name, "_", aln_best_pair_hom_DF$acc1_tname),
                          start = aln_best_pair_hom_DF$acc1_tstart-1,
                          end = aln_best_pair_hom_DF$acc1_tend,
                          aligner = aln_best_pair_hom_DF$acc1_aligner,
                          atype = aln_best_pair_hom_DF$acc1_atype,
                          mapq = aln_best_pair_hom_DF$acc1_mapq,
                          alen = aln_best_pair_hom_DF$acc1_alen,
                          nmatch = aln_best_pair_hom_DF$acc1_nmatch)

acc2_hom_bed = data.frame(chr = paste0(acc2_name, "_", aln_best_pair_hom_DF$acc2_tname),
                          start = aln_best_pair_hom_DF$acc2_tstart-1,
                          end = aln_best_pair_hom_DF$acc2_tend,
                          aligner = aln_best_pair_hom_DF$acc2_aligner,
                          atype = aln_best_pair_hom_DF$acc2_atype,
                          mapq = aln_best_pair_hom_DF$acc2_mapq,
                          alen = aln_best_pair_hom_DF$acc2_alen,
                          nmatch = aln_best_pair_hom_DF$acc2_nmatch)


# Define genome data.frame for circlize
acc1_genome_DF = data.frame(chr = acc1_chrs,
                            start = rep(0, length(acc1_chrs)),
                            end = acc1_chrLens)
acc2_genome_DF = data.frame(chr = acc2_chrs,
                            start = rep(0, length(acc2_chrs)),
                            end = acc2_chrLens)
genome_DF = rbind(acc1_genome_DF, acc2_genome_DF)
chr_index = c(rev(acc2_chrs), acc1_chrs)
genome_DF[, 1] = factor(genome_DF[, 1], levels = chr_index)

chrs = c(acc1_chrs, acc2_chrs)
CENstart = c(acc1_CENstart, acc2_CENstart)
CENend = c(acc1_CENend, acc2_CENend)

CEN_DF = data.frame(chr = chrs,
                    start = CENstart,
                    end = CENend)
CEN_DF[, 1] = factor(CEN_DF[, 1], levels = chr_index)

CEN_DF_extend = CEN_DF
CEN_DF_extend$start = CEN_DF$start - 2e6
CEN_DF_extend$end = CEN_DF$end - 2e6

# Reverse orientation of given sectors (chromosomes)
rev_x = function(x, xrange = CELL_META$xlim) {
    xrange[2] - x + xrange[1]
}


aligner_col_fun = c("wm" = "#009E73", "mm" = "#0072B2", "sr" = "#CC79A7")
atype_col_fun = structure(c("#332288", "#117733"), names = c("tp:A:P", "tp:A:S"))
mapq_col_fun = colorRamp2(quantile(c(aln_best_pair_DF$acc1_mapq, aln_best_pair_DF$acc2_mapq),
                                   c(0.05, 0.20, 0.40, 0.60, 0.80, 0.95),
                                   na.rm = T),
                          viridis(6))
alen_col_fun = colorRamp2(quantile(c(aln_best_pair_DF$acc1_alen, aln_best_pair_DF$acc2_alen),
                                   c(0.05, 0.20, 0.40, 0.60, 0.80, 0.95),
                                   na.rm = T),
                          plasma(6))
nmatch_col_fun = colorRamp2(quantile(c(aln_best_pair_DF$acc1_nmatch, aln_best_pair_DF$acc2_nmatch),
                                     c(0.05, 0.20, 0.40, 0.60, 0.80, 0.95),
                                     na.rm = T),
                            heat.colors(6))

# Define corresponding heatmap legends
lgd_aligner = Legend(at = c("wm", "mm", "sr"), type = "grid", legend_gp = gpar(col = aligner_col_fun), background = NULL, title = "Aligner", title_gp = gpar(fontface = "bold"), title_position = "leftcenter-rot")   
lgd_atype = Legend(at = c("tp:A:P", "tp:A:S"), type = "grid", legend_gp = gpar(col = atype_col_fun), background = NULL, title = "Type", title_gp = gpar(fontface = "bold"), title_position = "leftcenter-rot") 
lgd_mapq = Legend(col_fun = mapq_col_fun, title = "MAPQ", title_gp = gpar(fontface = "bold"), title_position = "leftcenter-rot")
lgd_alen = Legend(col_fun = alen_col_fun, title = "Length", title_gp = gpar(fontface = "bold"), title_position = "leftcenter-rot")
lgd_nmatch = Legend(col_fun = nmatch_col_fun, title = "Matches", title_gp = gpar(fontface = "bold"), title_position = "leftcenter-rot")
#lgd_list1 <- packLegend(lgd_aligner, lgd_atype, lgd_mapq, lgd_alen, lgd_nmatch)
lgd_list1 <- packLegend(lgd_mapq, lgd_alen)




# Initialize circular layout
circlize_plot = function(acc1_bed, acc2_bed, genome_DF) {
 
    circos.par(
               gap.after = c(rep(1, length(acc1_chrs)-1), 5, rep(1, length(acc2_chrs)-1), 5),
               track.height = 0.15
              )

    circos.genomicInitialize(data = genome_DF,
                             plotType = NULL,
                             tickLabelsStartFromZero = True)
    circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    }, track.height = mm_h(1), cell.padding = c(0, 0 , 0, 0), bg.border = NA)
    highlight.chromosome(acc1_chrs, col = "dodgerblue3", track.index = 1)
    highlight.chromosome(acc2_chrs, col = "darkgoldenrod", track.index = 1)

    set_track_gap(gap = 0.1)

    circos.track(ylim = c(0, 1),
                 bg.col = "grey70",
                 bg.border = NA,
                 track.height = 0.05,
                 panel.fun = function(x, y) {
                     major.by = seq(from = genome_DF[ which(genome_DF$chr == CELL_META$sector.index), 2],
                                    to   = genome_DF[ which(genome_DF$chr == CELL_META$sector.index), 3],
                                    by   = 4e6)
                     if(CELL_META$sector.index %in% acc2_chrs) {
                         circos.genomicAxis(major.at = rev_x(major.by), labels = paste0((major.by/1e6), " Mb"),
                                            minor.ticks = 1, labels.facing = "clockwise")
                     } else {
                         circos.genomicAxis(major.at = major.by, labels = paste0((major.by/1e6), " Mb"),
                                            minor.ticks = 1, labels.facing = "clockwise")
                     }
                 })

    # Reverse centromere coordinates for acc2
    CEN_DF_rev = CEN_DF
    for(l in which(CEN_DF_rev$chr %in% acc2_chrs)) {
        CEN_DF_rev$end[l] = rev_x(CEN_DF$start[l], c(genome_DF[l, 2], genome_DF[l, 3]))
        CEN_DF_rev$start[l] = rev_x(CEN_DF$end[l], c(genome_DF[l, 2], genome_DF[l, 3]))
    } 
    sapply(1:nrow(CEN_DF_rev), function(x) {
        circos.text(x = (CEN_DF_rev$start[x]+CEN_DF_rev$end[x])/2,
                    y = 0.5,
                    sector.index = CEN_DF_rev$chr[x],
                    track.index = 2,
                    labels = paste0("CEN", gsub(".*Chr", "", CEN_DF_rev$chr[x])),
                    niceFacing = TRUE,
                    cex = 0.8,
                    col = "black",
                    font = 4)
    })

    acc2_bed_rev = acc2_bed
    for(l in which(acc2_bed_rev$chr %in% acc2_chrs)) {
        chr_l = acc2_bed_rev$chr[l]
        acc2_bed_rev$end[l] = rev_x(acc2_bed$start[l], c(genome_DF[ which(genome_DF$chr == chr_l), 2], genome_DF[ which(genome_DF$chr == chr_l), 3]))
        acc2_bed_rev$start[l] = rev_x(acc2_bed$end[l], c(genome_DF[ which(genome_DF$chr == chr_l), 2], genome_DF[ which(genome_DF$chr == chr_l), 3]))
    } 

    set_track_gap(gap = 0.005)

    # mapq heatmap
    acc1_bed_mapq = acc1_bed[, c(1:3, which(colnames(acc1_bed) == "mapq"))]
    acc2_bed_mapq = acc2_bed_rev[, c(1:3, which(colnames(acc2_bed_rev) == "mapq"))]
    acc_bed_mapq = rbind(acc1_bed_mapq, acc2_bed_mapq)
    circos.genomicHeatmap(bed = do.call(rbind, lapply(seq_along(chrs), function(x) {
        acc_bed_mapq[ which( acc_bed_mapq$chr == genome_DF$chr[x] &
                             acc_bed_mapq$start >= genome_DF[ which(genome_DF$chr == chrs[x]), 2] &
                             acc_bed_mapq$end <= genome_DF[ which(genome_DF$chr == chrs[x]), 3] ), ] })),
        col = mapq_col_fun,
        border = NA,
        side = "inside",
        heatmap_height = 0.05,
        connection_height = mm_h(8))

    # alen heatmap
    acc1_bed_alen = acc1_bed[, c(1:3, which(colnames(acc1_bed) == "alen"))]
    acc2_bed_alen = acc2_bed_rev[, c(1:3, which(colnames(acc2_bed_rev) == "alen"))]
    acc_bed_alen = rbind(acc1_bed_alen, acc2_bed_alen)
    circos.genomicHeatmap(bed = do.call(rbind, lapply(seq_along(chrs), function(x) {
        acc_bed_alen[ which( acc_bed_alen$chr == genome_DF$chr[x] &
                             acc_bed_alen$start >= genome_DF[ which(genome_DF$chr == chrs[x]), 2] &
                             acc_bed_alen$end <= genome_DF[ which(genome_DF$chr == chrs[x]), 3] ), ] })),
        col = alen_col_fun,
        border = NA,
        side = "inside",
        heatmap_height = 0.05,
        connection_height = mm_h(8))

#    # nmatch heatmap
#    acc1_bed_nmatch = acc1_bed[, c(1:3, which(colnames(acc1_bed) == "nmatch"))]
#    acc2_bed_nmatch = acc2_bed_rev[, c(1:3, which(colnames(acc2_bed_rev) == "nmatch"))]
#    acc_bed_nmatch = rbind(acc1_bed_nmatch, acc2_bed_nmatch)
#    circos.genomicHeatmap(bed = do.call(rbind, lapply(seq_along(chrs), function(x) {
#        acc_bed_nmatch[ which( acc_bed_nmatch$chr == genome_DF$chr[x] &
#                             acc_bed_nmatch$start >= genome_DF[ which(genome_DF$chr == chrs[x]), 2] &
#                             acc_bed_nmatch$end <= genome_DF[ which(genome_DF$chr == chrs[x]), 3] ), ] })),
#        col = nmatch_col_fun,
#        border = NA,
#        side = "inside",
#        heatmap_height = 0.05,
#        connection_height = mm_h(8))

    set_track_gap(gap = 0.005)

    # Links between acc1 and acc2 aligned read segment pairs
    circos.genomicLink(acc1_bed, acc2_bed_rev, col = rand_color(nrow(acc1_bed)))
#    circos.genomicLink(acc1_bed, acc2_bed_rev,
#                       col = mapq_col_fun(sapply(1:nrow(acc1_bed_mapq), function(x) mean(acc1_bed_mapq$mapq[x], acc2_bed_mapq$mapq[x]))))
#    circos.genomicLink(acc1_bed, acc2_bed_rev,
#                       col = alen_col_fun(sapply(1:nrow(acc1_bed_alen), function(x) mean(acc1_bed_alen$alen[x], acc2_bed_alen$alen[x]))))

    # Reset graphic parameters and internal variables
    circos.clear()
 
    text(-0.9, -0.8, paste0(acc2_name, "\ngenome"))
    text(0.9, 0.8, paste0(acc1_name, "\ngenome"))

}

pdf(paste0(plotDir,
           acc1, "_", acc2, "_all_",
           "putative_not_centromeric_", recombType, "_",
           paste0(chrName, collapse = "_"), "_circlize_v", date, ".pdf"))                                                                
circlize_plot(acc1_bed = acc1_all_bed,
              acc2_bed = acc2_all_bed,
              genome_DF = genome_DF)
draw(lgd_list1, x = unit(1, "npc") - unit(2, "mm"), y = unit(4, "mm"), just = c("right", "bottom"))
dev.off()

pdf(paste0(plotDir,
           acc1, "_", acc2, "_hom_",
           "putative_not_centromeric_", recombType, "_",
           paste0(chrName, collapse = "_"), "_circlize_v", date, ".pdf"))                                                                
circlize_plot(acc1_bed = acc1_hom_bed,
              acc2_bed = acc2_hom_bed,
              genome_DF = genome_DF)
draw(lgd_list1, x = unit(1, "npc") - unit(2, "mm"), y = unit(4, "mm"), just = c("right", "bottom"))
dev.off()

#pdf(paste0(plotDir,
#           acc1, "_", acc2, "_all_",
#           "putative_not_centromeric_", recombType, "_",
#           paste0(chrName, collapse = "_"), "_circlize_zoom_v", date, ".pdf"))                                                                
#circlize_plot(acc1_bed = acc1_all_bed,
#              acc2_bed = acc2_all_bed,
#              genome_DF = CEN_DF_extend)
#draw(lgd_list1, x = unit(1, "npc") - unit(2, "mm"), y = unit(4, "mm"), just = c("right", "bottom"))
#dev.off()
#
#pdf(paste0(plotDir,
#           acc1, "_", acc2, "_hom_",
#           "putative_not_centromeric_", recombType, "_",
#           paste0(chrName, collapse = "_"), "_circlize_zoom_v", date, ".pdf"))                                                                
#circlize_plot(acc1_bed = acc1_hom_bed,
#              acc2_bed = acc2_hom_bed,
#              genome_DF = CEN_DF_extend)
#draw(lgd_list1, x = unit(1, "npc") - unit(2, "mm"), y = unit(4, "mm"), just = c("right", "bottom"))
#dev.off()




#  # Gypsy heatmap
#  circos.genomicHeatmap(bed = do.call(rbind, lapply(seq_along(chrs), function(x) {
#    Gypsy_bed[Gypsy_bed$chr == chrs[x] &
#              Gypsy_bed$start >= genomeDF$start[x] &
#              Gypsy_bed$end <= genomeDF$end[x],] } )),
#    col = Gypsy_col_fun,
#    border = NA,
#    side = "inside",
#    heatmap_height = 0.05,
#    connection_height = NULL)
#  # CEN180 heatmap
#  circos.genomicHeatmap(bed = do.call(rbind, lapply(seq_along(chrs), function(x) {
#    CEN180freq_bed[CEN180freq_bed$chr == chrs[x] &
#                   CEN180freq_bed$start >= genomeDF$start[x] &
#                   CEN180freq_bed$end <= genomeDF$end[x],] } )),
#    col = CEN180freq_col_fun,
#    border = NA,
#    side = "inside",
#    heatmap_height = 0.05,
#    connection_height = NULL)
#  # CENAthila rainfall plot
#  circos.genomicRainfall(data = do.call(rbind, lapply(seq_along(chrs), function(x) {
#    CENAthila_bed[CENAthila_bed$chr == chrs[x] &
#                  CENAthila_bed$start >= genomeDF$start[x] &
#                  CENAthila_bed$end <= genomeDF$end[x],] } )),
#                         bg.border = NA,
#                         track.height = 0.05,
#                         pch = 16,
#                         cex = 0.4,
#                         col = c("#0000FF80"))
#
#  set_track_gap(gap = 0.005)
#
#  # Plot windowed CEN180 frequency for each quantile
#  circos.genomicTrack(data = do.call(rbind, lapply(seq_along(chrs), function(x) {
#    CENH3_in_bodies_CEN180_bed[CENH3_in_bodies_CEN180_bed$chr == chrs[x] &                                                                    
#                               CENH3_in_bodies_CEN180_bed$start >= genomeDF$start[x] &                                                        
#                               CENH3_in_bodies_CEN180_bed$end <= genomeDF$end[x],] } )),                                                      
#                      panel.fun = function(region, value, ...) {
#                        circos.genomicLines(region,
#                                            value,
#                                            col = quantileColours,
#                                            lwd = 1.5,
#                                            lty = 1,
#                                            area = FALSE,
#                                            ...)
#                      }, bg.border = NA)
#  circos.yaxis(side = "left", sector.index = get.all.sector.index()[1], labels.cex = 0.3, tick.length = convert_x(0.5, "mm"), lwd = 0.5)      
#
#  set_track_gap(gap = 0.005)
#
#  circos.genomicTrack(data = do.call(rbind, lapply(seq_along(chrs), function(x) {                                                             
#    HORlengthsSum_CEN180_bed[HORlengthsSum_CEN180_bed$chr == chrs[x] &
#                             HORlengthsSum_CEN180_bed$start >= genomeDF$start[x] &                                                            
#                             HORlengthsSum_CEN180_bed$end <= genomeDF$end[x],] } )),                                                          
#                      panel.fun = function(region, value, ...) {
#                        circos.genomicLines(region,
#                                            value,
#                                            col = quantileColours,
#                                            lwd = 1.5,
#                                            lty = 1,
#                                            area = FALSE,
#                                            ...)
#                      }, bg.border = NA)
#  circos.yaxis(side = "left", sector.index = get.all.sector.index()[1], labels.cex = 0.3, tick.length = convert_x(0.5, "mm"), lwd = 0.5)      
#
#  set_track_gap(gap = 0.005)
#
#  circos.genomicTrack(data = do.call(rbind, lapply(seq_along(chrs), function(x) {
#    wSNV_CEN180_bed[wSNV_CEN180_bed$chr == chrs[x] &
#                    wSNV_CEN180_bed$start >= genomeDF$start[x] &
#                    wSNV_CEN180_bed$end <= genomeDF$end[x],] } )),
#                      panel.fun = function(region, value, ...) {
#                        circos.genomicLines(region,
#                                            value,
#                                            col = quantileColours,
#                                            lwd = 1.5,
#                                            lty = 1,
#                                            area = FALSE,
#                                            ...)
#                      }, bg.border = NA)
#  circos.yaxis(side = "left", sector.index = get.all.sector.index()[1], labels.cex = 0.3, tick.length = convert_x(0.5, "mm"), lwd = 0.5)
#
#  set_track_gap(gap = 0.005)
#
#  circos.genomicTrack(data = do.call(rbind, lapply(seq_along(chrs), function(x) {
#    CENH3_bed[CENH3_bed$chr == chrs[x] &
#              CENH3_bed$start >= genomeDF$start[x] &
#              CENH3_bed$end <= genomeDF$end[x],] } )),
#                      panel.fun = function(region, value, ...) {
#                        circos.genomicLines(region,
#                                            value,
#                                            col = "purple",
#                                            area = TRUE,
#                                            baseline = 0,
#                                            border = NA,
#                                            ...)
#                      }, bg.border = NA)
#  circos.yaxis(side = "left", at = seq(0, 4, by = 2), sector.index = get.all.sector.index()[1], labels.cex = 0.3, tick.length = convert_x(0.5, "mm"), lwd = 0.5)
#
#  # Reset graphic parameters and internal variables
#  circos.clear()
#}
#
#
#pdf(paste0(plotDir,
#           "CEN180_frequency_per_", genomeBinName,
#           "_", quantileDef, "_", quantiles, "quantiles",
#           "_of_CEN180_in_t2t-col.20210610_",
#           paste0(chrName, collapse = "_"), "_circlize_zoom_v", date, ".pdf"))                                                                
#circlize_plot()
#draw(lgd_list2, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))                                                             
#draw(lgd_list1, x = unit(1, "npc") - unit(2, "mm"), y = unit(4, "mm"), just = c("right", "bottom"))                                           
#dev.off()
#
