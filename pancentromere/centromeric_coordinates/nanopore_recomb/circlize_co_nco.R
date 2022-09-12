#!/usr/bin/env Rscript

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 09/09/2022

# Plot putative crossover events detected in Col-0/Ler-0
# hybrid ONT reads


# Usage:
# ./circlize_co_nco.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 10000 Col-0.ragtag_scaffolds Ler-0_110x.ragtag_scaffolds co 090922

#chrName = unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split=","))
#genomeBinSize = 10000
#acc1 = "Col-0.ragtag_scaffolds"
#acc2 = "Ler-0_110x.ragtag_scaffolds"
#recombType = "co"
#date = "090922"

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
acc1_chrs = paste0(acc1_chrs, "_", acc1_name)

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
acc2_chrs = paste0(acc2_chrs, "_", acc2_name)


CENstart = c(acc1_CENstart, acc2_CENstart)
CENend = c(acc1_CENend, acc2_CENend)


# Load read segment alignment files as a combined data.frame
load_pafs = function(indir, acc_name, suffix, aligner) {
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

# acc1 alignments list
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



bed1=generateRandomBed(nr=100)
head(bed1)
bed1=bed1[sample(nrow(bed1), 20),]
head(bed1)




# circlize

# Define genome data.frame for circlize

acc1_genome_DF = data.frame(chr = acc1_chrs,
                            start = rep(0, length(acc1_chrs)),
                            end = acc1_chrLens)
acc2_genome_DF = data.frame(chr = acc2_chrs,
                            start = rep(0, length(acc2_chrs)),
                            end = acc2_chrLens)
genome_DF = rbind(acc1_genome_DF, acc2_genome_DF)

chr_index = c(acc1_chrs, rev(acc2_chrs))









# Initialize circular layout
circlize_plot = function() {
  gapDegree = 6
  circos.par(
             gap.after = c(rep(1, length(acc1_chrs)), 5, rep(1, length(acc2_chrs))),
             track.height = 0.15,
             canvas.xlim = c(-1.1, 1.1),
             canvas.ylim = c(-1.1, 1.1),
             gap.degree = c(rep(1, length(chrs)-1), gapDegree),
             start.degree = 90-(gapDegree/2)
             )
  circos.genomicInitialize(data = genome_DF,
                           plotType = NULL,
                           tickLabelsStartFromZero = FALSE)
  circos.track(ylim = c(0, 1),
               bg.col = "grey70",
               bg.border = NA,
               track.height = 0.05,
               panel.fun = function(x, y) {
                 circos.genomicAxis(h = "top",
                                    labels.facing = "clockwise",
                                    tickLabelsStartFromZero = FALSE)
               })
  sapply(seq_along(chrs), function(x) {
    circos.text(x = (CENstart[x]+CENend[x])/2,
                y = 0.5,
                sector.index = chrs[x],
                track.index = 1,
                labels = paste0("CEN", x),
                niceFacing = TRUE,
                cex = 0.8,
                col = "black",
                font = 4)
  })

  # Gypsy heatmap
  circos.genomicHeatmap(bed = do.call(rbind, lapply(seq_along(chrs), function(x) {
    Gypsy_bed[Gypsy_bed$chr == chrs[x] &
              Gypsy_bed$start >= genomeDF$start[x] &
              Gypsy_bed$end <= genomeDF$end[x],] } )),
    col = Gypsy_col_fun,
    border = NA,
    side = "inside",
    heatmap_height = 0.05,
    connection_height = NULL)
  # CEN180 heatmap
  circos.genomicHeatmap(bed = do.call(rbind, lapply(seq_along(chrs), function(x) {
    CEN180freq_bed[CEN180freq_bed$chr == chrs[x] &
                   CEN180freq_bed$start >= genomeDF$start[x] &
                   CEN180freq_bed$end <= genomeDF$end[x],] } )),
    col = CEN180freq_col_fun,
    border = NA,
    side = "inside",
    heatmap_height = 0.05,
    connection_height = NULL)
  # CENAthila rainfall plot
  circos.genomicRainfall(data = do.call(rbind, lapply(seq_along(chrs), function(x) {
    CENAthila_bed[CENAthila_bed$chr == chrs[x] &
                  CENAthila_bed$start >= genomeDF$start[x] &
                  CENAthila_bed$end <= genomeDF$end[x],] } )),
                         bg.border = NA,
                         track.height = 0.05,
                         pch = 16,
                         cex = 0.4,
                         col = c("#0000FF80"))

  set_track_gap(gap = 0.005)

  # Plot windowed CEN180 frequency for each quantile
  circos.genomicTrack(data = do.call(rbind, lapply(seq_along(chrs), function(x) {
    CENH3_in_bodies_CEN180_bed[CENH3_in_bodies_CEN180_bed$chr == chrs[x] &                                                                    
                               CENH3_in_bodies_CEN180_bed$start >= genomeDF$start[x] &                                                        
                               CENH3_in_bodies_CEN180_bed$end <= genomeDF$end[x],] } )),                                                      
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region,
                                            value,
                                            col = quantileColours,
                                            lwd = 1.5,
                                            lty = 1,
                                            area = FALSE,
                                            ...)
                      }, bg.border = NA)
  circos.yaxis(side = "left", sector.index = get.all.sector.index()[1], labels.cex = 0.3, tick.length = convert_x(0.5, "mm"), lwd = 0.5)      

  set_track_gap(gap = 0.005)

  circos.genomicTrack(data = do.call(rbind, lapply(seq_along(chrs), function(x) {                                                             
    HORlengthsSum_CEN180_bed[HORlengthsSum_CEN180_bed$chr == chrs[x] &
                             HORlengthsSum_CEN180_bed$start >= genomeDF$start[x] &                                                            
                             HORlengthsSum_CEN180_bed$end <= genomeDF$end[x],] } )),                                                          
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region,
                                            value,
                                            col = quantileColours,
                                            lwd = 1.5,
                                            lty = 1,
                                            area = FALSE,
                                            ...)
                      }, bg.border = NA)
  circos.yaxis(side = "left", sector.index = get.all.sector.index()[1], labels.cex = 0.3, tick.length = convert_x(0.5, "mm"), lwd = 0.5)      

  set_track_gap(gap = 0.005)

  circos.genomicTrack(data = do.call(rbind, lapply(seq_along(chrs), function(x) {
    wSNV_CEN180_bed[wSNV_CEN180_bed$chr == chrs[x] &
                    wSNV_CEN180_bed$start >= genomeDF$start[x] &
                    wSNV_CEN180_bed$end <= genomeDF$end[x],] } )),
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region,
                                            value,
                                            col = quantileColours,
                                            lwd = 1.5,
                                            lty = 1,
                                            area = FALSE,
                                            ...)
                      }, bg.border = NA)
  circos.yaxis(side = "left", sector.index = get.all.sector.index()[1], labels.cex = 0.3, tick.length = convert_x(0.5, "mm"), lwd = 0.5)

  set_track_gap(gap = 0.005)

  circos.genomicTrack(data = do.call(rbind, lapply(seq_along(chrs), function(x) {
    CENH3_bed[CENH3_bed$chr == chrs[x] &
              CENH3_bed$start >= genomeDF$start[x] &
              CENH3_bed$end <= genomeDF$end[x],] } )),
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region,
                                            value,
                                            col = "purple",
                                            area = TRUE,
                                            baseline = 0,
                                            border = NA,
                                            ...)
                      }, bg.border = NA)
  circos.yaxis(side = "left", at = seq(0, 4, by = 2), sector.index = get.all.sector.index()[1], labels.cex = 0.3, tick.length = convert_x(0.5, "mm"), lwd = 0.5)

  # Reset graphic parameters and internal variables
  circos.clear()
}


pdf(paste0(plotDir,
           "CEN180_frequency_per_", genomeBinName,
           "_", quantileDef, "_", quantiles, "quantiles",
           "_of_CEN180_in_t2t-col.20210610_",
           paste0(chrName, collapse = "_"), "_circlize_zoom_v", date, ".pdf"))                                                                
circlize_plot()
draw(lgd_list2, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))                                                             
draw(lgd_list1, x = unit(1, "npc") - unit(2, "mm"), y = unit(4, "mm"), just = c("right", "bottom"))                                           
dev.off()

