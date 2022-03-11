#!/applications/R/R-4.0.0/bin/Rscript

# Analyse correlations between CEN180 distance to the nearest CENATHILA
# and CEN180 HOR and divergence metrics

# Usage:
# ./58Atha_CEN180_vs_CENATHILA_v170222.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))

options(stringsAsFactors = F)
library(GenomicRanges)
library(parallel)
library(dplyr)
library(scales)
library(ggplot2)
library(cowplot)
library(viridis)

plotDir <- paste0("plots/")
plotDirHORlengthsSum <- paste0(plotDir, "HORlengthsSum/")
plotDirHORcount <- paste0(plotDir, "HORcount/")
plotDirWeightedConsensusScore <- paste0(plotDir, "WeightedConsensusScore/")
plotDirEditDistance <- paste0(plotDir, "EditDistance/")
plotDirAllMetrics <- paste0(plotDir, "AllMetrics/")
plotDirAllAccessions <- paste0(plotDir, "AllAccessions/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))
system(paste0("[ -d ", plotDirHORlengthsSum, " ] || mkdir -p ", plotDirHORlengthsSum))
system(paste0("[ -d ", plotDirHORcount, " ] || mkdir -p ", plotDirHORcount))
system(paste0("[ -d ", plotDirWeightedConsensusScore, " ] || mkdir -p ", plotDirWeightedConsensusScore))
system(paste0("[ -d ", plotDirEditDistance, " ] || mkdir -p ", plotDirEditDistance))
system(paste0("[ -d ", plotDirAllMetrics, " ] || mkdir -p ", plotDirAllMetrics))
system(paste0("[ -d ", plotDirAllAccessions, " ] || mkdir -p ", plotDirAllAccessions))

acc_full <- system("ls /home/ajt200/analysis/pancentromere/annotation/CEN180/repeats/*.fa*", intern = T)
acc_full <- gsub("/home/ajt200/analysis/pancentromere/annotation/CEN180/repeats/cen180.consensus.repetitiveness", "", acc_full)
acc_uniq <- unique( gsub("\\.fa\\..+", "", acc_full) )
acc_uniq_len <- NULL
for(i in acc_uniq) {
  acc_uniq_len <- c(acc_uniq_len, length( grep(i , acc_full)) )
}
acc <- acc_uniq[which(acc_uniq_len == length(chrName))]

CEN180_list <- mclapply(1:length(acc), function(x) {
  tab_list <- lapply(1:length(chrName), function(y) {
    read.csv(paste0("/home/ajt200/analysis/pancentromere/annotation/CEN180/repeats/",
                    "cen180.consensus.repetitiveness", acc[x], ".fa.", chrName[y], ".csv"),
             header = T)
  })
  if(length(chrName) > 1) {
    tab <- dplyr::bind_rows(tab_list)
  } else {
    tab <- tab_list[[1]]
  }
  tab <- tab[tab$class == "cen180",]
  colnames(tab)[9] <- "chr"
  tab$chr <- gsub("_RagTag_RagTag", "", tab$chr)
  tab$chr <- gsub("chr", "Chr", tab$chr)
  tab$fasta.file.name <- gsub("\\.fa.+", "", tab$fasta.file.name)
  tab <- tab[
             with( tab, order(chr, start, end) ),
            ]
  tab
}, mc.cores = length(acc), mc.preschedule = F)

CENATHILA_list <- mclapply(1:length(acc), function(x) {
  tab <- read.table(paste0("/home/ajt200/analysis/pancentromere/annotation/ATHILA/ATHILA/",
                           acc[x], "/CENATHILA_in_", acc[x], "_",
                           paste0(chrName, collapse = "_"), ".bed"),
                    header = F)
  tab$V2 <- tab$V2+1
  colnames(tab) <- c("chr", "start", "end", "name", "phylo", "strand")
  tab <- tab[
             with( tab, order(chr, start, end) ),
            ]
  tab
}, mc.cores = length(acc), mc.preschedule = F)

CENATHILA_DF_phylo <- dplyr::bind_rows(CENATHILA_list) 
phylo <- sort(unique(CENATHILA_DF_phylo$phylo))

CENranLoc_list <- mclapply(1:length(acc), function(x) {
  tab <- read.table(paste0("/home/ajt200/analysis/pancentromere/annotation/ATHILA/ATHILA/",
                           acc[x], "/CENATHILA_in_", acc[x], "_",
                           paste0(chrName, collapse = "_"), "_CENrandomLoci.bed"),
                    header = F)
  tab$V2 <- tab$V2+1
  colnames(tab) <- c("chr", "start", "end", "name", "phylo", "strand")
  tab <- tab[
             with( tab, order(chr, start, end) ),
            ]
  tab
}, mc.cores = length(acc), mc.preschedule = F)

# Function to get distance between each CEN180 and the nearest CENATHILA
CEN180distToCENATHILA <- function(CEN180, CENATHILA, CENranLoc) {
  CEN180_distTo_CENATHILA <- data.frame()
  for(i in 1:length(chrName)) {
    print(chrName[i])

    CEN180_chr <- CEN180[CEN180$chr == chrName[i],]
    CENATHILA_chr <- CENATHILA[CENATHILA$chr == chrName[i],]
    CENranLoc_chr <- CENranLoc[CENranLoc$chr == chrName[i],]

    if(nrow(CENATHILA_chr) > 0) {
      # Calculate distances from the start and end coordinates of each CEN180 and CENATHILA
      CEN180start_vs_CENATHILAstart <- mclapply(1:length(CEN180_chr$start), function(x) {
        abs(CEN180_chr$start[x] - CENATHILA_chr$start)
      }, mc.cores = detectCores(), mc.preschedule = T) 
      CEN180end_vs_CENATHILAend <- mclapply(1:length(CEN180_chr$end), function(x) {
        abs(CEN180_chr$end[x] - CENATHILA_chr$end)
      }, mc.cores = detectCores(), mc.preschedule = T) 
      CEN180start_vs_CENATHILAend <- mclapply(1:length(CEN180_chr$start), function(x) {
        abs(CEN180_chr$start[x] - CENATHILA_chr$end)
      }, mc.cores = detectCores(), mc.preschedule = T) 
      CEN180end_vs_CENATHILAstart <- mclapply(1:length(CEN180_chr$end), function(x) {
        abs(CEN180_chr$end[x] - CENATHILA_chr$start)
      }, mc.cores = detectCores(), mc.preschedule = T) 

      # Get distance between each CEN180 and the
      # nearest CENATHILA
      minDistToCENATHILA <- unlist(mclapply(1:length(CEN180_chr$start), function(x) {
        min(c(CEN180start_vs_CENATHILAstart[[x]],
              CEN180end_vs_CENATHILAend[[x]],
              CEN180start_vs_CENATHILAend[[x]],
              CEN180end_vs_CENATHILAstart[[x]]), na.rm = T)
      }, mc.cores = detectCores(), mc.preschedule = T))
      stopifnot(nrow(CEN180_chr) == length(minDistToCENATHILA))

      #phylo_chr <- unique(CENATHILA_chr$phylo)
      minDistToCENATHILA_phylo_list <- lapply(1:length(phylo), function(x) {
        CENATHILA_chr_phylo <- CENATHILA_chr[CENATHILA_chr$phylo == phylo[x],]

        if(nrow(CENATHILA_chr_phylo) > 0) {
          # Calculate distances from the start and end coordinates of each CEN180 and CENATHILA
          CEN180start_vs_CENATHILAstart <- mclapply(1:length(CEN180_chr$start), function(x) {
            abs(CEN180_chr$start[x] - CENATHILA_chr_phylo$start)
          }, mc.cores = detectCores(), mc.preschedule = T) 
          CEN180end_vs_CENATHILAend <- mclapply(1:length(CEN180_chr$end), function(x) {
            abs(CEN180_chr$end[x] - CENATHILA_chr_phylo$end)
          }, mc.cores = detectCores(), mc.preschedule = T) 
          CEN180start_vs_CENATHILAend <- mclapply(1:length(CEN180_chr$start), function(x) {
            abs(CEN180_chr$start[x] - CENATHILA_chr_phylo$end)
          }, mc.cores = detectCores(), mc.preschedule = T) 
          CEN180end_vs_CENATHILAstart <- mclapply(1:length(CEN180_chr$end), function(x) {
            abs(CEN180_chr$end[x] - CENATHILA_chr_phylo$start)
          }, mc.cores = detectCores(), mc.preschedule = T) 

          # Get distance between each CEN180 and the
          # nearest CENATHILA
          minDistToCENATHILA <- unlist(mclapply(1:length(CEN180_chr$start), function(x) {
            min(c(CEN180start_vs_CENATHILAstart[[x]],
                  CEN180end_vs_CENATHILAend[[x]],
                  CEN180start_vs_CENATHILAend[[x]],
                  CEN180end_vs_CENATHILAstart[[x]]), na.rm = T)
          }, mc.cores = detectCores(), mc.preschedule = T))
          stopifnot(nrow(CEN180_chr) == length(minDistToCENATHILA))
          minDistToCENATHILA
        } else {
          minDistToCENATHILA <- rep(Inf, length(CEN180_chr$start))
          minDistToCENATHILA
        } 
      })
    } else {
      minDistToCENATHILA <- rep(Inf, length(CEN180_chr$start))
      minDistToCENATHILA_phylo_list <- lapply(1:length(phylo), function(x) { rep(Inf, length(CEN180_chr$start)) } )
    } 

    minDistToCENATHILA_phylo_DF <- as.data.frame(dplyr::bind_cols(minDistToCENATHILA_phylo_list))
    minDistToCENATHILA_DF <- cbind(minDistToCENATHILA, minDistToCENATHILA_phylo_DF)
    colnames(minDistToCENATHILA_DF) <- c("ATHILA", phylo)
    colnames(minDistToCENATHILA_DF) <- gsub("ATHILA", "minDistToCENATHILA", colnames(minDistToCENATHILA_DF))

    if(nrow(CENranLoc_chr) > 0) {
      # Calculate distances from the start and end coordinates of each CEN180 and CENranLoc
      CEN180start_vs_CENranLocstart <- mclapply(1:length(CEN180_chr$start), function(x) {
        abs(CEN180_chr$start[x] - CENranLoc_chr$start)
      }, mc.cores = detectCores(), mc.preschedule = T) 
      CEN180end_vs_CENranLocend <- mclapply(1:length(CEN180_chr$end), function(x) {
        abs(CEN180_chr$end[x] - CENranLoc_chr$end)
      }, mc.cores = detectCores(), mc.preschedule = T) 
      CEN180start_vs_CENranLocend <- mclapply(1:length(CEN180_chr$start), function(x) {
        abs(CEN180_chr$start[x] - CENranLoc_chr$end)
      }, mc.cores = detectCores(), mc.preschedule = T) 
      CEN180end_vs_CENranLocstart <- mclapply(1:length(CEN180_chr$end), function(x) {
        abs(CEN180_chr$end[x] - CENranLoc_chr$start)
      }, mc.cores = detectCores(), mc.preschedule = T) 

      # Get distance between each CEN180 and the
      # nearest CENranLoc
      minDistToCENranLoc <- unlist(mclapply(1:length(CEN180_chr$start), function(x) {
        min(c(CEN180start_vs_CENranLocstart[[x]],
              CEN180end_vs_CENranLocend[[x]],
              CEN180start_vs_CENranLocend[[x]],
              CEN180end_vs_CENranLocstart[[x]]), na.rm = T)
      }, mc.cores = detectCores(), mc.preschedule = T))
      stopifnot(nrow(CEN180_chr) == length(minDistToCENranLoc))

      #phylo_chr <- unique(CENranLoc_chr$phylo)
      minDistToCENranLoc_phylo_list <- lapply(1:length(phylo), function(x) {
        CENranLoc_chr_phylo <- CENranLoc_chr[CENranLoc_chr$phylo == phylo[x],]

        if(nrow(CENranLoc_chr_phylo) > 0) {
          # Calculate distances from the start and end coordinates of each CEN180 and CENranLoc
          CEN180start_vs_CENranLocstart <- mclapply(1:length(CEN180_chr$start), function(x) {
            abs(CEN180_chr$start[x] - CENranLoc_chr_phylo$start)
          }, mc.cores = detectCores(), mc.preschedule = T) 
          CEN180end_vs_CENranLocend <- mclapply(1:length(CEN180_chr$end), function(x) {
            abs(CEN180_chr$end[x] - CENranLoc_chr_phylo$end)
          }, mc.cores = detectCores(), mc.preschedule = T) 
          CEN180start_vs_CENranLocend <- mclapply(1:length(CEN180_chr$start), function(x) {
            abs(CEN180_chr$start[x] - CENranLoc_chr_phylo$end)
          }, mc.cores = detectCores(), mc.preschedule = T) 
          CEN180end_vs_CENranLocstart <- mclapply(1:length(CEN180_chr$end), function(x) {
            abs(CEN180_chr$end[x] - CENranLoc_chr_phylo$start)
          }, mc.cores = detectCores(), mc.preschedule = T) 

          # Get distance between each CEN180 and the
          # nearest CENranLoc
          minDistToCENranLoc <- unlist(mclapply(1:length(CEN180_chr$start), function(x) {
            min(c(CEN180start_vs_CENranLocstart[[x]],
                  CEN180end_vs_CENranLocend[[x]],
                  CEN180start_vs_CENranLocend[[x]],
                  CEN180end_vs_CENranLocstart[[x]]), na.rm = T)
          }, mc.cores = detectCores(), mc.preschedule = T))
          stopifnot(nrow(CEN180_chr) == length(minDistToCENranLoc))
          minDistToCENranLoc
        } else {
          minDistToCENranLoc <- rep(Inf, length(CEN180_chr$start))
          minDistToCENranLoc
        } 
      })
    } else {
      minDistToCENranLoc <- rep(Inf, length(CEN180_chr$start))
      minDistToCENranLoc_phylo_list <- lapply(1:length(phylo), function(x) { rep(Inf, length(CEN180_chr$start)) } )
    } 

    minDistToCENranLoc_phylo_DF <- as.data.frame(dplyr::bind_cols(minDistToCENranLoc_phylo_list))
    minDistToCENranLoc_DF <- cbind(minDistToCENranLoc, minDistToCENranLoc_phylo_DF)
    colnames(minDistToCENranLoc_DF) <- c("ranLoc_ATHILA", paste0("ranLoc_", phylo))
    colnames(minDistToCENranLoc_DF) <- gsub("ranLoc", "minDistToCENranLoc", colnames(minDistToCENranLoc_DF))

    CEN180_distTo_CENATHILA_chr <- data.frame(CEN180_chr,
                                              minDistToCENATHILA_DF,
                                              minDistToCENranLoc_DF)
    CEN180_distTo_CENATHILA <- rbind(CEN180_distTo_CENATHILA, CEN180_distTo_CENATHILA_chr) 
  }
  CEN180_distTo_CENATHILA
}


# For each acc, get distance between each CEN180 and the nearest CENATHILA
CEN180_list_dist <- lapply(1:length(acc), function(x) {
  print(acc[x])
  CEN180distToCENATHILA(CEN180 = CEN180_list[[x]],
                        CENATHILA = CENATHILA_list[[x]],
                        CENranLoc = CENranLoc_list[[x]])
})

rm(CEN180_list, CENATHILA_list, CENranLoc_list); gc()


# Plot relationships and define groups
trendPlot <- function(acc_id, dataFrame, mapping, xvar, yvar, xlab, ylab, xaxtrans, yaxtrans, xbreaks, ybreaks, xlabels, ylabels) {
  xvar <- enquo(xvar)
  yvar <- enquo(yvar)
  ggplot(data = dataFrame,
         mapping = mapping) +
  geom_hex(bins = 60) +
  scale_x_continuous(trans = xaxtrans,
                     breaks = xbreaks,
                     labels = xlabels) +
  scale_y_continuous(trans = yaxtrans,
                     breaks = ybreaks,
                     labels = ylabels) +
  annotation_logticks(sides = "lb") +
  scale_fill_viridis() +
  geom_smooth(colour = "red", fill = "grey70", alpha = 0.9,
              method = "gam", formula = y ~ s(x, bs = "cs")) +
  labs(x = xlab,
       y = ylab) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 0.5, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        axis.line = element_line(size = 1.0, colour = "black"),
        panel.background = element_blank(),
        panel.border = element_blank(),
#        panel.border = element_rect(size = 1.0, colour = "black"),
#        panel.grid = element_blank(),
        strip.text.x = element_text(size = 20, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"),
        plot.margin = unit(c(0.3,1.2,0.3,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ggtitle(bquote(
                 .(acc_id) ~
                 "			" ~
                 "Chr1" ~ italic(r[s]) ~ "=" ~
                 .(round(cor.test(select(dataFrame[dataFrame[,9] == "Chr1",], !!enquo(xvar))[,1], select(dataFrame[dataFrame[,9] == "Chr1",], !!enquo(yvar))[,1], method = "spearman", use = "pairwise.complete.obs")$estimate[[1]],
                         digits = 2)) *
                 ";" ~ italic(P) ~ "=" ~
                 .(round(min(0.5, cor.test(select(dataFrame[dataFrame[,9] == "Chr1",], !!enquo(xvar))[,1], select(dataFrame[dataFrame[,9] == "Chr1",], !!enquo(yvar))[,1], method = "spearman", use = "pairwise.complete.obs")$p.value * sqrt( (dim(dataFrame[dataFrame[,9] == "Chr1",])[1]/100) )),
                         digits = 5)) ~
                 "			" ~
                 "Chr2" ~ italic(r[s]) ~ "=" ~
                 .(round(cor.test(select(dataFrame[dataFrame[,9] == "Chr2",], !!enquo(xvar))[,1], select(dataFrame[dataFrame[,9] == "Chr2",], !!enquo(yvar))[,1], method = "spearman", use = "pairwise.complete.obs")$estimate[[1]],
                         digits = 2)) *
                 ";" ~ italic(P) ~ "=" ~
                 .(round(min(0.5, cor.test(select(dataFrame[dataFrame[,9] == "Chr2",], !!enquo(xvar))[,1], select(dataFrame[dataFrame[,9] == "Chr2",], !!enquo(yvar))[,1], method = "spearman", use = "pairwise.complete.obs")$p.value * sqrt( (dim(dataFrame[dataFrame[,9] == "Chr2",])[1]/100) )),
                         digits = 5)) ~
                 "			" ~
                 "Chr3" ~ italic(r[s]) ~ "=" ~
                 .(round(cor.test(select(dataFrame[dataFrame[,9] == "Chr3",], !!enquo(xvar))[,1], select(dataFrame[dataFrame[,9] == "Chr3",], !!enquo(yvar))[,1], method = "spearman", use = "pairwise.complete.obs")$estimate[[1]],
                         digits = 2)) *
                 ";" ~ italic(P) ~ "=" ~
                 .(round(min(0.5, cor.test(select(dataFrame[dataFrame[,9] == "Chr3",], !!enquo(xvar))[,1], select(dataFrame[dataFrame[,9] == "Chr3",], !!enquo(yvar))[,1], method = "spearman", use = "pairwise.complete.obs")$p.value * sqrt( (dim(dataFrame[dataFrame[,9] == "Chr3",])[1]/100) )),
                         digits = 5)) ~
                 "			" ~
                 "Chr4" ~ italic(r[s]) ~ "=" ~
                 .(round(cor.test(select(dataFrame[dataFrame[,9] == "Chr4",], !!enquo(xvar))[,1], select(dataFrame[dataFrame[,9] == "Chr4",], !!enquo(yvar))[,1], method = "spearman", use = "pairwise.complete.obs")$estimate[[1]],
                         digits = 2)) *
                 ";" ~ italic(P) ~ "=" ~
                 .(round(min(0.5, cor.test(select(dataFrame[dataFrame[,9] == "Chr4",], !!enquo(xvar))[,1], select(dataFrame[dataFrame[,9] == "Chr4",], !!enquo(yvar))[,1], method = "spearman", use = "pairwise.complete.obs")$p.value * sqrt( (dim(dataFrame[dataFrame[,9] == "Chr4",])[1]/100) )),
                         digits = 5)) ~
                 "			" ~
                 "Chr5" ~ italic(r[s]) ~ "=" ~
                 .(round(cor.test(select(dataFrame[dataFrame[,9] == "Chr5",], !!enquo(xvar))[,1], select(dataFrame[dataFrame[,9] == "Chr5",], !!enquo(yvar))[,1], method = "spearman", use = "pairwise.complete.obs")$estimate[[1]],
                         digits = 2)) *
                 ";" ~ italic(P) ~ "=" ~
                 .(round(min(0.5, cor.test(select(dataFrame[dataFrame[,9] == "Chr5",], !!enquo(xvar))[,1], select(dataFrame[dataFrame[,9] == "Chr5",], !!enquo(yvar))[,1], method = "spearman", use = "pairwise.complete.obs")$p.value * sqrt( (dim(dataFrame[dataFrame[,9] == "Chr5",])[1]/100) )),
                         digits = 5)) ~
                 "			" ~
                 "All" ~ italic(r[s]) ~ "=" ~
                 .(round(cor.test(select(dataFrame, !!enquo(xvar))[,1], select(dataFrame, !!enquo(yvar))[,1], method = "spearman", use = "pairwise.complete.obs")$estimate[[1]],
                         digits = 2)) *
                 ";" ~ italic(P) ~ "=" ~
                 .(round(min(0.5, cor.test(select(dataFrame, !!enquo(xvar))[,1], select(dataFrame, !!enquo(yvar))[,1], method = "spearman", use = "pairwise.complete.obs")$p.value * sqrt( (dim(dataFrame)[1]/100) )),
                         digits = 5))
                )
         )
}


# Changes to each data.frame for plotting
CEN180_list_dist_tmp <- lapply(1:length(acc), function(x) {
  colnames(CEN180_list_dist[[x]]) <- gsub("minDistTo", "", colnames(CEN180_list_dist[[x]]))
  # Add offset (+ 1) for plotting on log scales
  CEN180_list_dist[[x]][,12:ncol(CEN180_list_dist[[x]])] <- CEN180_list_dist[[x]][,12:ncol(CEN180_list_dist[[x]])] + 1
  CEN180_list_dist[[x]]$est.chromosome.no <- NA
  CEN180_list_dist[[x]]
})

CEN180_list_dist_tmp <- lapply(1:length(acc), function(x) {
  CEN180_list_dist_tmp[[x]]$est.chromosome.no <- NA
  CEN180_list_dist_tmp[[x]]
})

# Bind rows of per-accession data.frames
CEN180_list_dist_tmp_allacc <- dplyr::bind_rows(CEN180_list_dist_tmp)

#rm(CEN180_list_dist); gc()

phylo_ext <- paste0("CEN", c("ATHILA", phylo))

#                    mapping = aes(x = CEN180_list_dist_tmp[[i]][,which(names(CEN180_list_dist_tmp[[i]]) == phylo_ext[j])] + 1,
#                                  y = HORlengthsSum + 1),
#                    xvar = as.name(names(CEN180_list_dist_tmp[[i]])[which(names(CEN180_list_dist_tmp[[i]]) == phylo_ext[j])]),

ggTrend_minDistToCENATHILA_HORlengthsSum_listOlists <- lapply(1:length(acc), function(i) {
  lapply(1:length(phylo_ext), function(j) {
    tP <- trendPlot(acc_id = acc[i],
                    dataFrame = CEN180_list_dist_tmp[[i]],
                    mapping = aes_(x = as.name(phylo_ext[j]),
                                   y = as.name("HORlengthsSum")),
                    xvar = as.name(phylo_ext[j]),
                    yvar = HORlengthsSum,
                    xlab = bquote("Distance to nearest" ~ italic(.(phylo_ext[j])) ~ "(bp)"),
                    ylab = bquote(italic("CEN180") ~ "repetitiveness"),
                    xaxtrans = log10_trans(),
                    yaxtrans = log10_trans(),
                    xbreaks = trans_breaks("log10", function(x) 10^x),
                    ybreaks = trans_breaks("log10", function(x) 10^x),
                    xlabels = trans_format("log10", math_format(10^.x)),
                    ylabels = trans_format("log10", math_format(10^.x)))
    tP <- tP +
      facet_grid(cols = vars(chr), margins = "chr", scales = "free_x")
    tP
  })
})
#}, mc.cores = detectCores(), mc.preschedule = F)

lapply(1:length(acc), function(i) {
  gg_cow_list <- ggTrend_minDistToCENATHILA_HORlengthsSum_listOlists[[i]]
  gg_cow <- plot_grid(plotlist = gg_cow_list,
                      labels = "AUTO", label_size = 30,
                      align = "hv",
                      axis = "l",
                      nrow = length(gg_cow_list), ncol = 1)
  ggsave(paste0(plotDirHORlengthsSum,
                "CEN180_MinDistToCENATHILA_vs_HORlengthsSum_trendPlot_",
                paste0(chrName, collapse = "_"),
                "_", acc[i], ".pdf"),
         plot = gg_cow,
         height = 5.5*length(gg_cow_list), width = 5*(length(chrName)+1), limitsize = F)
})
#}, mc.cores = detectCores(), mc.preschedule = F)


ggTrend_minDistToCENATHILA_HORcount_listOlists <- lapply(1:length(acc), function(i) {
  lapply(1:length(phylo_ext), function(j) {
    tP <- trendPlot(acc_id = acc[i],
                    dataFrame = CEN180_list_dist_tmp[[i]],
                    mapping = aes_(x = as.name(phylo_ext[j]),
                                   y = as.name("HORcount")),
                    xvar = as.name(phylo_ext[j]),
                    yvar = HORcount,
                    xlab = bquote("Distance to nearest" ~ italic(.(phylo_ext[j])) ~ "(bp)"),
                    ylab = bquote(italic("CEN180") ~ "HOR count"),
                    xaxtrans = log10_trans(),
                    yaxtrans = log10_trans(),
                    xbreaks = trans_breaks("log10", function(x) 10^x),
                    ybreaks = trans_breaks("log10", function(x) 10^x),
                    xlabels = trans_format("log10", math_format(10^.x)),
                    ylabels = trans_format("log10", math_format(10^.x)))
    tP <- tP +
      facet_grid(cols = vars(chr), margins = "chr", scales = "free_x")
    tP
  })
})

lapply(1:length(acc), function(i) {
  gg_cow_list <- ggTrend_minDistToCENATHILA_HORcount_listOlists[[i]]
  gg_cow <- plot_grid(plotlist = gg_cow_list,
                      labels = "AUTO", label_size = 30,
                      align = "hv",
                      axis = "l",
                      nrow = length(gg_cow_list), ncol = 1)
  ggsave(paste0(plotDirHORcount,
                "CEN180_MinDistToCENATHILA_vs_HORcount_trendPlot_",
                paste0(chrName, collapse = "_"),
                "_", acc[i], ".pdf"),
         plot = gg_cow,
         height = 5.5*length(gg_cow_list), width = 5*(length(chrName)+1), limitsize = F)
})


ggTrend_minDistToCENATHILA_WeightedConsensusScore_listOlists <- lapply(1:length(acc), function(i) {
  lapply(1:length(phylo_ext), function(j) {
    tP <- trendPlot(acc_id = acc[i],
                    dataFrame = CEN180_list_dist_tmp[[i]],
                    mapping = aes_(x = as.name(phylo_ext[j]),
                                   y = as.name("weighted.consensus.score")),
                    xvar = as.name(phylo_ext[j]),
                    yvar = weighted.consensus.score,
                    xlab = bquote("Distance to nearest" ~ italic(.(phylo_ext[j])) ~ "(bp)"),
                    ylab = bquote(italic("CEN180") ~ "consensus score"),
                    xaxtrans = log10_trans(),
                    yaxtrans = log2_trans(),
                    xbreaks = trans_breaks("log10", function(x) 10^x),
                    ybreaks = trans_breaks("log2", function(x) 2^x),
                    xlabels = trans_format("log10", math_format(10^.x)),
                    ylabels = trans_format("log2", math_format(2^.x)))
    tP <- tP +
      facet_grid(cols = vars(chr), margins = "chr", scales = "free_x")
    tP
  })
})

lapply(1:length(acc), function(i) {
  gg_cow_list <- ggTrend_minDistToCENATHILA_WeightedConsensusScore_listOlists[[i]]
  gg_cow <- plot_grid(plotlist = gg_cow_list,
                      labels = "AUTO", label_size = 30,
                      align = "hv",
                      axis = "l",
                      nrow = length(gg_cow_list), ncol = 1)
  ggsave(paste0(plotDirWeightedConsensusScore,
                "CEN180_MinDistToCENATHILA_vs_WeightedConsensusScore_trendPlot_",
                paste0(chrName, collapse = "_"),
                "_", acc[i], ".pdf"),
         plot = gg_cow,
         height = 5.5*length(gg_cow_list), width = 5*(length(chrName)+1), limitsize = F)
})


ggTrend_minDistToCENATHILA_EditDistance_listOlists <- lapply(1:length(acc), function(i) {
  lapply(1:length(phylo_ext), function(j) {
    tP <- trendPlot(acc_id = acc[i],
                    dataFrame = CEN180_list_dist_tmp[[i]],
                    mapping = aes_(x = as.name(phylo_ext[j]),
                                   y = as.name("edit.distance")),
                    xvar = as.name(phylo_ext[j]),
                    yvar = edit.distance,
                    xlab = bquote("Distance to nearest" ~ italic(.(phylo_ext[j])) ~ "(bp)"),
                    ylab = bquote(italic("CEN180") ~ "edit distance"),
                    xaxtrans = log10_trans(),
                    yaxtrans = log2_trans(),
                    xbreaks = trans_breaks("log10", function(x) 10^x),
                    ybreaks = trans_breaks("log2", function(x) 2^x),
                    xlabels = trans_format("log10", math_format(10^.x)),
                    ylabels = trans_format("log2", math_format(2^.x)))
    tP <- tP +
      facet_grid(cols = vars(chr), margins = "chr", scales = "free_x")
    tP
  })
})

lapply(1:length(acc), function(i) {
  gg_cow_list <- ggTrend_minDistToCENATHILA_EditDistance_listOlists[[i]]
  gg_cow <- plot_grid(plotlist = gg_cow_list,
                      labels = "AUTO", label_size = 30,
                      align = "hv",
                      axis = "l",
                      nrow = length(gg_cow_list), ncol = 1)
  ggsave(paste0(plotDirEditDistance,
                "CEN180_MinDistToCENATHILA_vs_EditDistance_trendPlot_",
                paste0(chrName, collapse = "_"),
                "_", acc[i], ".pdf"),
         plot = gg_cow,
         height = 5.5*length(gg_cow_list), width = 5*(length(chrName)+1), limitsize = F)
})


# Plot all CEN180 stats comparisons with distance to nearest CENATHILA,
# but not including subfamily-specific relationships
lapply(1:length(acc), function(i) {
  gg_cow_list <- list(
                      ggTrend_minDistToCENATHILA_HORlengthsSum_listOlists[[i]][[1]],
                      ggTrend_minDistToCENATHILA_HORcount_listOlists[[i]][[1]],
                      ggTrend_minDistToCENATHILA_WeightedConsensusScore_listOlists[[i]][[1]],
                      ggTrend_minDistToCENATHILA_EditDistance_listOlists[[i]][[1]]
                     )
  gg_cow <- plot_grid(plotlist = gg_cow_list,
                      labels = "AUTO", label_size = 30,
                      align = "hv",
                      axis = "l",
                      nrow = length(gg_cow_list), ncol = 1)
  ggsave(paste0(plotDirAllMetrics,
                "CEN180_MinDistToCENATHILA_vs_AllMetrics_trendPlot_",
                paste0(chrName, collapse = "_"),
                "_", acc[i], ".pdf"),
         plot = gg_cow,
         height = 5.5*length(gg_cow_list), width = 5*(length(chrName)+1), limitsize = F)
})



# Plot correlations across all accessions,
# including subfamily-specific relationships
ggTrend_minDistToCENATHILA_HORlengthsSum_list <- lapply(1:length(phylo_ext), function(j) {
  tP <- trendPlot(acc_id = paste0(length(unique(CEN180_list_dist_tmp_allacc$fasta.file.name)), " accessions"),
                  dataFrame = CEN180_list_dist_tmp_allacc,
                  mapping = aes_(x = as.name(phylo_ext[j]),
                                 y = as.name("HORlengthsSum")),
                  xvar = as.name(phylo_ext[j]),
                  yvar = HORlengthsSum,
                  xlab = bquote("Distance to nearest" ~ italic(.(phylo_ext[j])) ~ "(bp)"),
                  ylab = bquote(italic("CEN180") ~ "repetitiveness"),
                  xaxtrans = log10_trans(),
                  yaxtrans = log10_trans(),
                  xbreaks = trans_breaks("log10", function(x) 10^x),
                  ybreaks = trans_breaks("log10", function(x) 10^x),
                  xlabels = trans_format("log10", math_format(10^.x)),
                  ylabels = trans_format("log10", math_format(10^.x)))
  tP <- tP +
    facet_grid(cols = vars(chr), margins = "chr", scales = "free_x")
  tP
})

gg_cow_list <- ggTrend_minDistToCENATHILA_HORlengthsSum_list
gg_cow <- plot_grid(plotlist = gg_cow_list,
                    labels = "AUTO", label_size = 30,
                    align = "hv",
                    axis = "l",
                    nrow = length(gg_cow_list), ncol = 1)
ggsave(paste0(plotDirAllAccessions,
              "CEN180_MinDistToCENATHILA_vs_HORlengthsSum_trendPlot_",
              paste0(chrName, collapse = "_"),
              "_", length(unique(CEN180_list_dist_tmp_allacc$fasta.file.name)), "accessions.pdf"),
       plot = gg_cow,
       height = 5.5*length(gg_cow_list), width = 5*(length(chrName)+1), limitsize = F)


ggTrend_minDistToCENATHILA_HORcount_list <- lapply(1:length(phylo_ext), function(j) {
  tP <- trendPlot(acc_id = paste0(length(unique(CEN180_list_dist_tmp_allacc$fasta.file.name)), " accessions"),
                  dataFrame = CEN180_list_dist_tmp_allacc,
                  mapping = aes_(x = as.name(phylo_ext[j]),
                                 y = as.name("HORcount")),
                  xvar = as.name(phylo_ext[j]),
                  yvar = HORcount,
                  xlab = bquote("Distance to nearest" ~ italic(.(phylo_ext[j])) ~ "(bp)"),
                  ylab = bquote(italic("CEN180") ~ "HOR count"),
                  xaxtrans = log10_trans(),
                  yaxtrans = log10_trans(),
                  xbreaks = trans_breaks("log10", function(x) 10^x),
                  ybreaks = trans_breaks("log10", function(x) 10^x),
                  xlabels = trans_format("log10", math_format(10^.x)),
                  ylabels = trans_format("log10", math_format(10^.x)))
  tP <- tP +
    facet_grid(cols = vars(chr), margins = "chr", scales = "free_x")
  tP
})

gg_cow_list <- ggTrend_minDistToCENATHILA_HORcount_list
gg_cow <- plot_grid(plotlist = gg_cow_list,
                    labels = "AUTO", label_size = 30,
                    align = "hv",
                    axis = "l",
                    nrow = length(gg_cow_list), ncol = 1)
ggsave(paste0(plotDirAllAccessions,
              "CEN180_MinDistToCENATHILA_vs_HORcount_trendPlot_",
              paste0(chrName, collapse = "_"),
              "_", length(unique(CEN180_list_dist_tmp_allacc$fasta.file.name)), "accessions.pdf"),
       plot = gg_cow,
       height = 5.5*length(gg_cow_list), width = 5*(length(chrName)+1), limitsize = F)


ggTrend_minDistToCENATHILA_WeightedConsensusScore_list <- lapply(1:length(phylo_ext), function(j) {
  tP <- trendPlot(acc_id = paste0(length(unique(CEN180_list_dist_tmp_allacc$fasta.file.name)), " accessions"),
                  dataFrame = CEN180_list_dist_tmp_allacc,
                  mapping = aes_(x = as.name(phylo_ext[j]),
                                 y = as.name("weighted.consensus.score")),
                  xvar = as.name(phylo_ext[j]),
                  yvar = weighted.consensus.score,
                  xlab = bquote("Distance to nearest" ~ italic(.(phylo_ext[j])) ~ "(bp)"),
                  ylab = bquote(italic("CEN180") ~ "consensus score"),
                  xaxtrans = log10_trans(),
                  yaxtrans = log2_trans(),
                  xbreaks = trans_breaks("log10", function(x) 10^x),
                  ybreaks = trans_breaks("log2", function(x) 2^x),
                  xlabels = trans_format("log10", math_format(10^.x)),
                  ylabels = trans_format("log2", math_format(2^.x)))
  tP <- tP +
    facet_grid(cols = vars(chr), margins = "chr", scales = "free_x")
  tP
})

gg_cow_list <- ggTrend_minDistToCENATHILA_WeightedConsensusScore_list
gg_cow <- plot_grid(plotlist = gg_cow_list,
                    labels = "AUTO", label_size = 30,
                    align = "hv",
                    axis = "l",
                    nrow = length(gg_cow_list), ncol = 1)
ggsave(paste0(plotDirAllAccessions,
              "CEN180_MinDistToCENATHILA_vs_WeightedConsensusScore_trendPlot_",
              paste0(chrName, collapse = "_"),
              "_", length(unique(CEN180_list_dist_tmp_allacc$fasta.file.name)), "accessions.pdf"),
       plot = gg_cow,
       height = 5.5*length(gg_cow_list), width = 5*(length(chrName)+1), limitsize = F)


ggTrend_minDistToCENATHILA_EditDistance_list <- lapply(1:length(phylo_ext), function(j) {
  tP <- trendPlot(acc_id = paste0(length(unique(CEN180_list_dist_tmp_allacc$fasta.file.name)), " accessions"),
                  dataFrame = CEN180_list_dist_tmp_allacc,
                  mapping = aes_(x = as.name(phylo_ext[j]),
                                 y = as.name("edit.distance")),
                  xvar = as.name(phylo_ext[j]),
                  yvar = edit.distance,
                  xlab = bquote("Distance to nearest" ~ italic(.(phylo_ext[j])) ~ "(bp)"),
                  ylab = bquote(italic("CEN180") ~ "edit distance"),
                  xaxtrans = log10_trans(),
                  yaxtrans = log2_trans(),
                  xbreaks = trans_breaks("log10", function(x) 10^x),
                  ybreaks = trans_breaks("log2", function(x) 2^x),
                  xlabels = trans_format("log10", math_format(10^.x)),
                  ylabels = trans_format("log2", math_format(2^.x)))
  tP <- tP +
    facet_grid(cols = vars(chr), margins = "chr", scales = "free_x")
  tP
})

gg_cow_list <- ggTrend_minDistToCENATHILA_EditDistance_list
gg_cow <- plot_grid(plotlist = gg_cow_list,
                    labels = "AUTO", label_size = 30,
                    align = "hv",
                    axis = "l",
                    nrow = length(gg_cow_list), ncol = 1)
ggsave(paste0(plotDirAllAccessions,
              "CEN180_MinDistToCENATHILA_vs_EditDistance_trendPlot_",
              paste0(chrName, collapse = "_"),
              "_", length(unique(CEN180_list_dist_tmp_allacc$fasta.file.name)), "accessions.pdf"),
       plot = gg_cow,
       height = 5.5*length(gg_cow_list), width = 5*(length(chrName)+1), limitsize = F)


# Accross all accessions, plot all CEN180 stats comparisons with distance to nearest CENATHILA,
# but not including subfamily-specific relationships
gg_cow_list <- list(
                    ggTrend_minDistToCENATHILA_HORlengthsSum_list[[1]],
                    ggTrend_minDistToCENATHILA_HORcount_list[[1]],
                    ggTrend_minDistToCENATHILA_WeightedConsensusScore_list[[1]],
                    ggTrend_minDistToCENATHILA_EditDistance_list[[1]]
                   )
gg_cow <- plot_grid(plotlist = gg_cow_list,
                    labels = "AUTO", label_size = 30,
                    align = "hv",
                    axis = "l",
                    nrow = length(gg_cow_list), ncol = 1)
ggsave(paste0(plotDirAllAccessions,
              "CEN180_MinDistToCENATHILA_vs_AllMetrics_trendPlot_",
              paste0(chrName, collapse = "_"),
              "_", length(unique(CEN180_list_dist_tmp_allacc$fasta.file.name)), "accessions.pdf"),
       plot = gg_cow,
       height = 5.5*length(gg_cow_list), width = 5*(length(chrName)+1), limitsize = F)

