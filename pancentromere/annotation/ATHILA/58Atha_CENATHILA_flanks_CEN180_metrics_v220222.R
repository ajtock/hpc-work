#!/applications/R/R-4.0.0/bin/Rscript

# Compare average CEN180 metrics (HOR membership and divergence) in regions flanking
# centromeric ATHILA and matched centromeric random loci

# Usage:
# ./58Atha_CENATHILA_flanks_CEN180_metrics_v220222.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 500

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#flankSize <- 500

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
flankSize <- as.numeric(args[2])

if(floor(log10(flankSize)) + 1 < 4) {
  flankName <- paste0(flankSize, "bp")
} else if(floor(log10(flankSize)) + 1 >= 4 &
          floor(log10(flankSize)) + 1 <= 6) {
  flankName <- paste0(flankSize/1e3, "kb")
} else if(floor(log10(flankSize)) + 1 >= 7) {
  flankName <- paste0(flankSize/1e6, "Mb")
}
flankNamePlot <- paste0(c(strsplit(flankName, split = "")[[1]][1:(length(strsplit(flankName, split = "")[[1]])-2)],
                          "-",
                          strsplit(flankName, split = "")[[1]][(length(strsplit(flankName, split = "")[[1]])-1):(length(strsplit(flankName, split = "")[[1]]))]),
                          collapse = "")

options(stringsAsFactors = F)
source("/projects/meiosis/ajt200/Rfunctions/TTSplus.R")
library(GenomicRanges)
library(segmentSeq)
library(parallel)
library(dplyr)
library(scales)
library(ggplot2)
library(cowplot)
library(ggbeeswarm)
library(ggthemes)
library(viridis)

plotDir <- paste0("ATHILA/plots/")
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


# Get ranges corresponding to featRegion
if(featRegion == "bodies") {
  featGR <- featGR
} else if(featRegion == "promoters") {
  # Obtain 1000 bp upstream of start coordinates
  featGR <- promoters(featGR, upstream = 1000, downstream = 0)
} else if(featRegion == "terminators") {
  # Obtain 1000 bp downstream of end coordinates
  source("/projects/meiosis/ajt200/Rfunctions/TTSplus.R")
  featGR <- TTSplus(featGR, upstream = -1, downstream = 1000)
} else if(featRegion == "regions") {
  featGR <- GRanges(seqnames = seqnames(featGR),
                    ranges = IRanges(start = start(featGR)-1000,
                                     end = end(featGR)+1000),
                    strand = strand(featGR),
                    name = featGR$name,
                    score = featGR$score)
} else {
  stop("featRegion is none of bodies, promoters, terminators or regions")
}


# Function to calculate mean CEN180 metrics for regions
# upstream and downstream of CENATHILA and CENranLoc bodies
CEN180metricsAtCENATHILA <- function(CEN180, CENATHILA, featureName) {

  CEN180_GR <- GRanges(seqnames = as.character(CEN180$chr),
                       ranges = IRanges(start = as.integer(CEN180$start),
                                        end = as.integer(CEN180$end)),
                       strand = as.character(CEN180$strand),
                       acc = as.character(CEN180$fasta.file.name),
                       HORlengthsSum = as.numeric(CEN180$HORlengthsSum),
                       HORcount = as.numeric(CEN180$HORcount),
                       WeightedConsensusScore = as.numeric(CEN180$weighted.consensus.score),
                       EditDistance = as.numeric(CEN180$edit.distance))
  CENATHILA_GR <- GRanges(seqnames = as.character(CENATHILA$chr),
                          ranges = IRanges(start = as.integer(CENATHILA$start+1),
                                           end = as.integer(CENATHILA$end)),
                          strand = as.character(CENATHILA$strand),
                          phylo = as.character(CENATHILA$phylo))
 
  # Get flankSize bp upstream of start coordinates
  CENATHILA_up_GR <- promoters(CENATHILA_GR, upstream = flankSize, downstream = 0)
  CENATHILA_down_GR <- TTSplus(CENATHILA_GR, upstream = -1, downstream = flankSize)
  stopifnot(identical(CENATHILA_GR$phylo, CENATHILA_up_GR$phylo))
  stopifnot(identical(CENATHILA_GR$phylo, CENATHILA_down_GR$phylo))
  stopifnot(identical(seqnames(CENATHILA_GR), seqnames(CENATHILA_up_GR)))
  stopifnot(identical(seqnames(CENATHILA_GR), seqnames(CENATHILA_down_GR)))
  stopifnot(identical(strand(CENATHILA_GR), strand(CENATHILA_up_GR)))
  stopifnot(identical(strand(CENATHILA_GR), strand(CENATHILA_down_GR)))

  CENATHILA_CEN180_metrics <- data.frame()
  for(i in 1:length(chrName)) {
    print(chrName[i])

    CEN180_GR_chr <- CEN180_GR[seqnames(CEN180_GR) == chrName[i]]
    CENATHILA_GR_chr <- CENATHILA_GR[seqnames(CENATHILA_GR) == chrName[i]]
    CENATHILA_up_GR_chr <- CENATHILA_up_GR[seqnames(CENATHILA_up_GR) == chrName[i]]
    CENATHILA_down_GR_chr <- CENATHILA_down_GR[seqnames(CENATHILA_down_GR) == chrName[i]]
    stopifnot(identical(CENATHILA_GR_chr$phylo, CENATHILA_up_GR_chr$phylo))
    stopifnot(identical(CENATHILA_GR_chr$phylo, CENATHILA_down_GR_chr$phylo))
    stopifnot(identical(strand(CENATHILA_GR_chr), strand(CENATHILA_up_GR_chr)))
    stopifnot(identical(strand(CENATHILA_GR_chr), strand(CENATHILA_down_GR_chr)))

    if(length(CENATHILA_GR_chr) > 0) {
#     ## Note: findOverlaps() approach does not work where a query does not overlap
#     ##       any positions in subjects
#     CENATHILA_up_CEN180_foverlaps <- findOverlaps(query = CENATHILA_up_GR_chr,
#                                                   subject = CEN180_GR_chr,
#                                                   type = "any",
#                                                   select = "all",
#                                                   ignore.strand = TRUE)
#     # Convert CENATHILA_up_CEN180_foverlaps into list object equivalent to that
#     # generated by segmentSeq::getOverlaps(), in which each
#     # list element corresponds to a vector of sequentially numbered indices of
#     # CEN180_GR_chr that overlap a feature in CENATHILA_up_GR_chr
#     CENATHILA_up_CEN180 <- mclapply(seq_along(unique(queryHits(CENATHILA_up_CEN180_foverlaps))),
#                            function(x) {
#                              subjectHits(CENATHILA_up_CEN180_foverlaps[queryHits(CENATHILA_up_CEN180_foverlaps) == x])
#                            }, mc.cores = detectCores(), mc.preschedule = F)
      CENATHILA_up_CEN180 <- getOverlaps(coordinates = CENATHILA_up_GR_chr,
                                         segments = CEN180_GR_chr,
                                         overlapType = "overlapping",
                                         whichOverlaps = TRUE,
                                         ignoreStrand = TRUE)
      CENATHILA_up_CEN180_HORlengthsSum <- sapply(1:length(CENATHILA_up_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_up_CEN180[[x]]]$HORlengthsSum, na.rm = T)
      })
      CENATHILA_up_CEN180_HORcount <- sapply(1:length(CENATHILA_up_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_up_CEN180[[x]]]$HORcount, na.rm = T)
      })
      CENATHILA_up_CEN180_WeightedConsensusScore <- sapply(1:length(CENATHILA_up_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_up_CEN180[[x]]]$WeightedConsensusScore, na.rm = T)
      })
      CENATHILA_up_CEN180_EditDistance <- sapply(1:length(CENATHILA_up_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_up_CEN180[[x]]]$EditDistance, na.rm = T)
      })

      CENATHILA_down_CEN180 <- getOverlaps(coordinates = CENATHILA_down_GR_chr,
                                           segments = CEN180_GR_chr,
                                           overlapType = "overlapping",
                                           whichOverlaps = TRUE,
                                           ignoreStrand = TRUE)
      CENATHILA_down_CEN180_HORlengthsSum <- sapply(1:length(CENATHILA_down_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_down_CEN180[[x]]]$HORlengthsSum, na.rm = T)
      })
      CENATHILA_down_CEN180_HORcount <- sapply(1:length(CENATHILA_down_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_down_CEN180[[x]]]$HORcount, na.rm = T)
      })
      CENATHILA_down_CEN180_WeightedConsensusScore <- sapply(1:length(CENATHILA_down_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_down_CEN180[[x]]]$WeightedConsensusScore, na.rm = T)
      })
      CENATHILA_down_CEN180_EditDistance <- sapply(1:length(CENATHILA_down_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_down_CEN180[[x]]]$EditDistance, na.rm = T)
      })

      #CENATHILA_reg_GR_chr <- c(CENATHILA_up_GR_chr, CENATHILA_down_GR_chr)
      CENATHILA_reg_GR_chr <- c(CENATHILA_GR_chr, CENATHILA_GR_chr)
      CENATHILA_chr <- data.frame(CENATHILA_reg_GR_chr,
                                  feature = rep(featureName, length(CENATHILA_reg_GR_chr)),
                                  accession = rep(CEN180_GR_chr$acc[1], length(CENATHILA_reg_GR_chr)),
                                  region = rep(c("Upstream", "Downstream"), each = length(CENATHILA_GR_chr)),
                                  HORlengthsSum = c(CENATHILA_up_CEN180_HORlengthsSum, CENATHILA_down_CEN180_HORlengthsSum),
                                  HORcount = c(CENATHILA_up_CEN180_HORcount, CENATHILA_down_CEN180_HORcount),
                                  WeightedConsensusScore = c(CENATHILA_up_CEN180_WeightedConsensusScore, CENATHILA_down_CEN180_WeightedConsensusScore),
                                  EditDistance = c(CENATHILA_up_CEN180_EditDistance, CENATHILA_down_CEN180_EditDistance))
      colnames(CENATHILA_chr)[1] <- "chr"
      CENATHILA_chr$chr <- as.character(CENATHILA_chr$chr)
      CENATHILA_chr$strand <- as.character(CENATHILA_chr$strand)
      colnames(CENATHILA_chr)[which(colnames(CENATHILA_chr) == "phylo")] <- "Family"
      colnames(CENATHILA_chr)[which(colnames(CENATHILA_chr) == "feature")] <- "Feature"
      colnames(CENATHILA_chr)[which(colnames(CENATHILA_chr) == "region")] <- "Region"

      CENATHILA_CEN180_metrics <- rbind(CENATHILA_CEN180_metrics, CENATHILA_chr)
    } 
  }
  CENATHILA_CEN180_metrics
}

# For each acc, apply CEN180metricsAtCENATHILA() function
CENATHILA_CEN180_metrics_list <- mclapply(1:length(acc), function(x) {
  print(acc[x])
  CEN180metricsAtCENATHILA(CEN180 = CEN180_list[[x]],
                           CENATHILA = CENATHILA_list[[x]],
                           featureName = "CENATHILA")
}, mc.cores = detectCores(), mc.preschedule = F)

CENranLoc_CEN180_metrics_list <- mclapply(1:length(acc), function(x) {
  print(acc[x])
  CEN180metricsAtCENATHILA(CEN180 = CEN180_list[[x]],
                           CENATHILA = CENranLoc_list[[x]],
                           featureName = "CENranLoc")
}, mc.cores = detectCores(), mc.preschedule = F)

#rm(CEN180_list, CENATHILA_list, CENranLoc_list); gc()

CENfeats_CEN180_metrics_list <- lapply(1:length(acc), function(x) {
  rbind(CENATHILA_CEN180_metrics_list[[x]],
        CENranLoc_CEN180_metrics_list[[x]])
})

# Bind rows of per-accession data.frames
CENfeats_CEN180_metrics_allacc <- dplyr::bind_rows(CENfeats_CEN180_metrics_list)


# Function to make boxplot of CEN180 metrics overlapping
# regions flanking CENATHILA and CENranLoc
pointPlot <- function(acc_id, dataFrame, mapping, box_mapping, pvals, xlab, ylab, yaxtrans, ybreaks, ylabels) {
  ggplot(data = dataFrame,
         mapping = mapping) +
  geom_boxplot(inherit.aes = F,
               mapping = box_mapping,
               colour = "grey60") +
  geom_beeswarm(cex = 3,
                size = 3) +
  scale_y_continuous(trans = yaxtrans,
                     breaks = ybreaks,
                     labels = ylabels) +
  labs(x = xlab,
       y = ylab) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 0.5, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 16, colour = "black", angle = 45, vjust = 1.0, hjust = 1.0, face = "italic"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        axis.line = element_line(size = 1.0, colour = "black"),
        legend.text = element_text(face = "italic"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        strip.text.x = element_text(size = 20, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"),
        plot.margin = unit(c(0.3,1.2,0.3,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ggtitle(bquote(
                 .(acc_id) ~
                 "                      " ~
                 .(pvals)
                )
         )
}

# Function to make boxplot of CEN180 metrics overlapping
# regions flanking CENATHILA and CENranLoc
violinPlot <- function(acc_id, dataFrame, mapping, pvals, xlab, ylab, yaxtrans, ybreaks, ylabels) {
  ggplot(data = dataFrame,
         mapping = mapping) +
  geom_violin(scale = "count",
              trim = T,
              draw_quantiles = c(0.25, 0.50, 0.75)) +
  scale_y_continuous(trans = yaxtrans,
                     breaks = ybreaks,
                     labels = ylabels) +
  labs(x = xlab,
       y = ylab) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 0.5, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 16, colour = "black", angle = 45, vjust = 1.0, hjust = 1.0, face = "italic"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        axis.line = element_line(size = 1.0, colour = "black"),
        legend.text = element_text(face = "italic"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        strip.text.x = element_text(size = 20, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"),
        plot.margin = unit(c(0.3,1.2,0.3,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ggtitle(bquote(
                 .(acc_id) ~
                 "                      " ~
                 .(pvals)
                )
         )
}


# Disable scientific notation (e.g., 0.0001 rather than 1e-04)
options(scipen = 100)

ggPoint_HORlengthsSum_list <- lapply(1:length(acc), function(i) {
  print(i)
  CENfeats_CEN180_metrics_list[[i]]$Region <- factor(CENfeats_CEN180_metrics_list[[i]]$Region,
                                                     levels = sort(levels(factor(CENfeats_CEN180_metrics_list[[i]]$Region)), decreasing = T))

  acc_chrName <- c(sort(unique(CENfeats_CEN180_metrics_list[[i]]$chr)), "All")
  Utest_Pvals <- as.numeric(sapply(acc_chrName, function(z) {
    print(z)
    if(z != "All") {
      tt <- tryCatch(
                     wilcox.test(x = select(CENfeats_CEN180_metrics_list[[i]][CENfeats_CEN180_metrics_list[[i]]$chr == z &
                                                                              CENfeats_CEN180_metrics_list[[i]]$Feature == "CENATHILA",],
                                            "HORlengthsSum")[,1],
                                 y = select(CENfeats_CEN180_metrics_list[[i]][CENfeats_CEN180_metrics_list[[i]]$chr == z &
                                                                              CENfeats_CEN180_metrics_list[[i]]$Feature == "CENranLoc",],
                                            "HORlengthsSum")[,1],
                                 alternative = "two.sided")$p.value,
                     error = function(e) e
                    )
      if(is(tt, "error")) {
        NA
      } else {
        tt
      }
    } else {
      tt <- tryCatch(
                     wilcox.test(x = select(CENfeats_CEN180_metrics_list[[i]][CENfeats_CEN180_metrics_list[[i]]$Feature == "CENATHILA",],
                                            "HORlengthsSum")[,1],
                                 y = select(CENfeats_CEN180_metrics_list[[i]][CENfeats_CEN180_metrics_list[[i]]$Feature == "CENranLoc",],
                                            "HORlengthsSum")[,1],
                                 alternative = "two.sided")$p.value,
                     error = function(e) e
                    )
      if(is(tt, "error")) {
        NA
      } else {
        tt
      }
    }
  }))
  #Utest_Pvals <- sapply(1:length(Utest_list), function(z) { Utest_list[[z]]$p.value } )
  Utest_PvalsChar <- sapply(1:length(Utest_Pvals), function(z) {
    if(is.na(Utest_Pvals[z])) {
      paste0(acc_chrName[z], " MWW P = NA")
    } else if(Utest_Pvals[z] < 0.0001) {
      paste0(acc_chrName[z], " MWW P < 0.0001")
    } else {
      paste0(acc_chrName[z], " MWW P = ", as.character(round(Utest_Pvals[z], digits = 4)))
    }
  })

  tP <- pointPlot(acc_id = acc[i],
                  dataFrame = CENfeats_CEN180_metrics_list[[i]],
                  mapping = aes(x = Feature,
                                y = HORlengthsSum + 1,
                                shape = Region,
                                colour = Family),
                  box_mapping = aes(x = Feature,
                                    y = HORlengthsSum + 1),
                  pvals = paste0(Utest_PvalsChar, collapse = "                "),
                  xlab = bquote(.(flankNamePlot) ~ "regions flanking centromeric" ~ italic("ATHILA") ~ "and random loci"),
                  ylab = bquote(italic("CEN180") ~ "repetitiveness"),
                  yaxtrans = log2_trans(),
                  ybreaks = trans_breaks("log2", function(x) 2^x),
                  ylabels = trans_format("log2", math_format(2^.x)))
  tP <- tP +
    facet_grid(cols = vars(chr), margins = "chr", scales = "fixed")
  tP
})

ggPoint_HORcount_list <- lapply(1:length(acc), function(i) {
  print(i)
  CENfeats_CEN180_metrics_list[[i]]$Region <- factor(CENfeats_CEN180_metrics_list[[i]]$Region,
                                                     levels = sort(levels(factor(CENfeats_CEN180_metrics_list[[i]]$Region)), decreasing = T))

  acc_chrName <- c(sort(unique(CENfeats_CEN180_metrics_list[[i]]$chr)), "All")
  Utest_Pvals <- as.numeric(sapply(acc_chrName, function(z) {
    print(z)
    if(z != "All") {
      tt <- tryCatch(
                     wilcox.test(x = select(CENfeats_CEN180_metrics_list[[i]][CENfeats_CEN180_metrics_list[[i]]$chr == z &
                                                                              CENfeats_CEN180_metrics_list[[i]]$Feature == "CENATHILA",],
                                            "HORcount")[,1],
                                 y = select(CENfeats_CEN180_metrics_list[[i]][CENfeats_CEN180_metrics_list[[i]]$chr == z &
                                                                              CENfeats_CEN180_metrics_list[[i]]$Feature == "CENranLoc",],
                                            "HORcount")[,1],
                                 alternative = "two.sided")$p.value,
                     error = function(e) e
                    )
      if(is(tt, "error")) {
        NA
      } else {
        tt
      }
    } else {
      tt <- tryCatch(
                     wilcox.test(x = select(CENfeats_CEN180_metrics_list[[i]][CENfeats_CEN180_metrics_list[[i]]$Feature == "CENATHILA",],
                                            "HORcount")[,1],
                                 y = select(CENfeats_CEN180_metrics_list[[i]][CENfeats_CEN180_metrics_list[[i]]$Feature == "CENranLoc",],
                                            "HORcount")[,1],
                                 alternative = "two.sided")$p.value,
                     error = function(e) e
                    )
      if(is(tt, "error")) {
        NA
      } else {
        tt
      }
    }
  }))
  #Utest_Pvals <- sapply(1:length(Utest_list), function(z) { Utest_list[[z]]$p.value } )
  Utest_PvalsChar <- sapply(1:length(Utest_Pvals), function(z) {
    if(is.na(Utest_Pvals[z])) {
      paste0(acc_chrName[z], " MWW P = NA")
    } else if(Utest_Pvals[z] < 0.0001) {
      paste0(acc_chrName[z], " MWW P < 0.0001")
    } else {
      paste0(acc_chrName[z], " MWW P = ", as.character(round(Utest_Pvals[z], digits = 4)))
    }
  })

  tP <- pointPlot(acc_id = acc[i],
                  dataFrame = CENfeats_CEN180_metrics_list[[i]],
                  mapping = aes(x = Feature,
                                y = HORcount + 1,
                                shape = Region,
                                colour = Family),
                  box_mapping = aes(x = Feature,
                                    y = HORcount + 1),
                  pvals = paste0(Utest_PvalsChar, collapse = "                "),
                  xlab = bquote(.(flankNamePlot) ~ "regions flanking centromeric" ~ italic("ATHILA") ~ "and random loci"),
                  ylab = bquote(italic("CEN180") ~ "HOR count"),
                  yaxtrans = log2_trans(),
                  ybreaks = trans_breaks("log2", function(x) 2^x),
                  ylabels = trans_format("log2", math_format(2^.x)))
  tP <- tP +
    facet_grid(cols = vars(chr), margins = "chr", scales = "fixed")
  tP
})

ggPoint_WeightedConsensusScore_list <- lapply(1:length(acc), function(i) {
  print(i)
  CENfeats_CEN180_metrics_list[[i]]$Region <- factor(CENfeats_CEN180_metrics_list[[i]]$Region,
                                                     levels = sort(levels(factor(CENfeats_CEN180_metrics_list[[i]]$Region)), decreasing = T))

  acc_chrName <- c(sort(unique(CENfeats_CEN180_metrics_list[[i]]$chr)), "All")
  Utest_Pvals <- as.numeric(sapply(acc_chrName, function(z) {
    print(z)
    if(z != "All") {
      tt <- tryCatch(
                     wilcox.test(x = select(CENfeats_CEN180_metrics_list[[i]][CENfeats_CEN180_metrics_list[[i]]$chr == z &
                                                                              CENfeats_CEN180_metrics_list[[i]]$Feature == "CENATHILA",],
                                            "WeightedConsensusScore")[,1],
                                 y = select(CENfeats_CEN180_metrics_list[[i]][CENfeats_CEN180_metrics_list[[i]]$chr == z &
                                                                              CENfeats_CEN180_metrics_list[[i]]$Feature == "CENranLoc",],
                                            "WeightedConsensusScore")[,1],
                                 alternative = "two.sided")$p.value,
                     error = function(e) e
                    )
      if(is(tt, "error")) {
        NA
      } else {
        tt
      }
    } else {
      tt <- tryCatch(
                     wilcox.test(x = select(CENfeats_CEN180_metrics_list[[i]][CENfeats_CEN180_metrics_list[[i]]$Feature == "CENATHILA",],
                                            "WeightedConsensusScore")[,1],
                                 y = select(CENfeats_CEN180_metrics_list[[i]][CENfeats_CEN180_metrics_list[[i]]$Feature == "CENranLoc",],
                                            "WeightedConsensusScore")[,1],
                                 alternative = "two.sided")$p.value,
                     error = function(e) e
                    )
      if(is(tt, "error")) {
        NA
      } else {
        tt
      }
    }
  }))
  #Utest_Pvals <- sapply(1:length(Utest_list), function(z) { Utest_list[[z]]$p.value } )
  Utest_PvalsChar <- sapply(1:length(Utest_Pvals), function(z) {
    if(is.na(Utest_Pvals[z])) {
      paste0(acc_chrName[z], " MWW P = NA")
    } else if(Utest_Pvals[z] < 0.0001) {
      paste0(acc_chrName[z], " MWW P < 0.0001")
    } else {
      paste0(acc_chrName[z], " MWW P = ", as.character(round(Utest_Pvals[z], digits = 4)))
    }
  })

  tP <- pointPlot(acc_id = acc[i],
                  dataFrame = CENfeats_CEN180_metrics_list[[i]],
                  mapping = aes(x = Feature,
                                y = WeightedConsensusScore + 1,
                                shape = Region,
                                colour = Family),
                  box_mapping = aes(x = Feature,
                                    y = WeightedConsensusScore + 1),
                  pvals = paste0(Utest_PvalsChar, collapse = "                "),
                  xlab = bquote(.(flankNamePlot) ~ "regions flanking centromeric" ~ italic("ATHILA") ~ "and random loci"),
                  ylab = bquote(italic("CEN180") ~ "consensus score"),
                  yaxtrans = log2_trans(),
                  ybreaks = trans_breaks("log2", function(x) 2^x),
                  ylabels = trans_format("log2", math_format(2^.x)))
  tP <- tP +
    facet_grid(cols = vars(chr), margins = "chr", scales = "fixed")
  tP
})

ggPoint_EditDistance_list <- lapply(1:length(acc), function(i) {
  print(i)
  CENfeats_CEN180_metrics_list[[i]]$Region <- factor(CENfeats_CEN180_metrics_list[[i]]$Region,
                                                     levels = sort(levels(factor(CENfeats_CEN180_metrics_list[[i]]$Region)), decreasing = T))

  acc_chrName <- c(sort(unique(CENfeats_CEN180_metrics_list[[i]]$chr)), "All")
  Utest_Pvals <- as.numeric(sapply(acc_chrName, function(z) {
    print(z)
    if(z != "All") {
      tt <- tryCatch(
                     wilcox.test(x = select(CENfeats_CEN180_metrics_list[[i]][CENfeats_CEN180_metrics_list[[i]]$chr == z &
                                                                              CENfeats_CEN180_metrics_list[[i]]$Feature == "CENATHILA",],
                                            "EditDistance")[,1],
                                 y = select(CENfeats_CEN180_metrics_list[[i]][CENfeats_CEN180_metrics_list[[i]]$chr == z &
                                                                              CENfeats_CEN180_metrics_list[[i]]$Feature == "CENranLoc",],
                                            "EditDistance")[,1],
                                 alternative = "two.sided")$p.value,
                     error = function(e) e
                    )
      if(is(tt, "error")) {
        NA
      } else {
        tt
      }
    } else {
      tt <- tryCatch(
                     wilcox.test(x = select(CENfeats_CEN180_metrics_list[[i]][CENfeats_CEN180_metrics_list[[i]]$Feature == "CENATHILA",],
                                            "EditDistance")[,1],
                                 y = select(CENfeats_CEN180_metrics_list[[i]][CENfeats_CEN180_metrics_list[[i]]$Feature == "CENranLoc",],
                                            "EditDistance")[,1],
                                 alternative = "two.sided")$p.value,
                     error = function(e) e
                    )
      if(is(tt, "error")) {
        NA
      } else {
        tt
      }
    }
  }))
  #Utest_Pvals <- sapply(1:length(Utest_list), function(z) { Utest_list[[z]]$p.value } )
  Utest_PvalsChar <- sapply(1:length(Utest_Pvals), function(z) {
    if(is.na(Utest_Pvals[z])) {
      paste0(acc_chrName[z], " MWW P = NA")
    } else if(Utest_Pvals[z] < 0.0001) {
      paste0(acc_chrName[z], " MWW P < 0.0001")
    } else {
      paste0(acc_chrName[z], " MWW P = ", as.character(round(Utest_Pvals[z], digits = 4)))
    }
  })

  tP <- pointPlot(acc_id = acc[i],
                  dataFrame = CENfeats_CEN180_metrics_list[[i]],
                  mapping = aes(x = Feature,
                                y = EditDistance + 1,
                                shape = Region,
                                colour = Family),
                  box_mapping = aes(x = Feature,
                                    y = EditDistance + 1),
                  pvals = paste0(Utest_PvalsChar, collapse = "                "),
                  xlab = bquote(.(flankNamePlot) ~ "regions flanking centromeric" ~ italic("ATHILA") ~ "and random loci"),
                  ylab = bquote(italic("CEN180") ~ "edit distance"),
                  yaxtrans = log2_trans(),
                  ybreaks = trans_breaks("log2", function(x) 2^x),
                  ylabels = trans_format("log2", math_format(2^.x)))
  tP <- tP +
    facet_grid(cols = vars(chr), margins = "chr", scales = "fixed")
  tP
})

# Per-accession plots
mclapply(1:length(acc), function(i) {
  acc_chrName <- c(sort(unique(CENfeats_CEN180_metrics_list[[i]]$chr)), "All")
  gg_cow_list <- list(
                      ggPoint_HORlengthsSum_list[[i]],
                      ggPoint_HORcount_list[[i]],
                      ggPoint_WeightedConsensusScore_list[[i]],
                      ggPoint_EditDistance_list[[i]]
                     )
  gg_cow <- plot_grid(plotlist = gg_cow_list,
                      labels = "AUTO", label_size = 30,
                      align = "hv",
                      axis = "l",
                      nrow = length(gg_cow_list), ncol = 1)
  ggsave(paste0(plotDirAllMetrics,
                "CENATHILA_CEN180_", flankName, "_flanks_AllMetrics_pointPlot_",
                paste0(chrName, collapse = "_"),
                "_", acc[i], ".pdf"),
         plot = gg_cow,
         height = 5.5*length(gg_cow_list), width = 5*(length(acc_chrName)), limitsize = F)
}, mc.cores = detectCores(), mc.preschedule = F)




# Plot across all accessions, with facet rows as Family 

CENfeats_CEN180_metrics_allacc$Region <- factor(CENfeats_CEN180_metrics_allacc$Region,
                                                levels = sort(levels(factor(CENfeats_CEN180_metrics_allacc$Region)), decreasing = T))
acc_chrName <- c(sort(unique(CENfeats_CEN180_metrics_allacc$chr)), "All")

# HORlengthsSum
Utest_Pvals_HORlengthsSum <- as.numeric(sapply(acc_chrName, function(z) {
  print(z)
  if(z != "All") {
    tt <- tryCatch(
                   wilcox.test(x = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$chr == z &
                                                                         CENfeats_CEN180_metrics_allacc$Feature == "CENATHILA",],
                                          "HORlengthsSum")[,1],
                               y = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$chr == z &
                                                                         CENfeats_CEN180_metrics_allacc$Feature == "CENranLoc",],
                                          "HORlengthsSum")[,1],
                               alternative = "two.sided")$p.value,
                   error = function(e) e
                  )
    if(is(tt, "error")) {
      NA
    } else {
      tt
    }
  } else {
    tt <- tryCatch(
                   wilcox.test(x = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Feature == "CENATHILA",],
                                          "HORlengthsSum")[,1],
                               y = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Feature == "CENranLoc",],
                                          "HORlengthsSum")[,1],
                               alternative = "two.sided")$p.value,
                   error = function(e) e
                  )
    if(is(tt, "error")) {
      NA
    } else {
      tt
    }
  }
}))
#Utest_Pvals_HORlengthsSum <- sapply(1:length(Utest_list_HORlengthsSum), function(z) { Utest_list_HORlengthsSum[[z]]$p.value } )
Utest_PvalsChar_HORlengthsSum <- sapply(1:length(Utest_Pvals_HORlengthsSum), function(z) {
  if(is.na(Utest_Pvals_HORlengthsSum[z])) {
    paste0(acc_chrName[z], " MWW P = NA")
  } else if(Utest_Pvals_HORlengthsSum[z] < 0.0001) {
    paste0(acc_chrName[z], " MWW P < 0.0001")
  } else {
    paste0(acc_chrName[z], " MWW P = ", as.character(round(Utest_Pvals_HORlengthsSum[z], digits = 4)))
  }
})

ggPoint_HORlengthsSum_allacc <- violinPlot(acc_id = paste0(length(unique(CENfeats_CEN180_metrics_allacc$accession)), " accessions"),
                                           dataFrame = CENfeats_CEN180_metrics_allacc,
                                           mapping = aes(x = Feature,
                                                         y = HORlengthsSum + 1,
                                                         fill = Family),
                                           pvals = paste0(Utest_PvalsChar_HORlengthsSum, collapse = "                "),
                                           xlab = bquote(.(flankNamePlot) ~ "regions flanking centromeric" ~ italic("ATHILA") ~ "and random loci"),
                                           ylab = bquote(italic("CEN180") ~ "repetitiveness"),
                                           yaxtrans = log2_trans(),
                                           ybreaks = trans_breaks("log2", function(x) 2^x),
                                           ylabels = trans_format("log2", math_format(2^.x)))
ggPoint_HORlengthsSum_allacc <- ggPoint_HORlengthsSum_allacc +
  facet_grid(cols = vars(chr), rows = vars(Family), margins = "chr", scales = "fixed")


# HORcount
Utest_Pvals_HORcount <- as.numeric(sapply(acc_chrName, function(z) {
  print(z)
  if(z != "All") {
    tt <- tryCatch(
                   wilcox.test(x = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$chr == z &
                                                                         CENfeats_CEN180_metrics_allacc$Feature == "CENATHILA",],
                                          "HORcount")[,1],
                               y = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$chr == z &
                                                                         CENfeats_CEN180_metrics_allacc$Feature == "CENranLoc",],
                                          "HORcount")[,1],
                               alternative = "two.sided")$p.value,
                   error = function(e) e
                  )
    if(is(tt, "error")) {
      NA
    } else {
      tt
    }
  } else {
    tt <- tryCatch(
                   wilcox.test(x = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Feature == "CENATHILA",],
                                          "HORcount")[,1],
                               y = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Feature == "CENranLoc",],
                                          "HORcount")[,1],
                               alternative = "two.sided")$p.value,
                   error = function(e) e
                  )
    if(is(tt, "error")) {
      NA
    } else {
      tt
    }
  }
}))
#Utest_Pvals_HORcount <- sapply(1:length(Utest_list_HORcount), function(z) { Utest_list_HORcount[[z]]$p.value } )
Utest_PvalsChar_HORcount <- sapply(1:length(Utest_Pvals_HORcount), function(z) {
  if(is.na(Utest_Pvals_HORcount[z])) {
    paste0(acc_chrName[z], " MWW P = NA")
  } else if(Utest_Pvals_HORcount[z] < 0.0001) {
    paste0(acc_chrName[z], " MWW P < 0.0001")
  } else {
    paste0(acc_chrName[z], " MWW P = ", as.character(round(Utest_Pvals_HORcount[z], digits = 4)))
  }
})

ggPoint_HORcount_allacc <- violinPlot(acc_id = paste0(length(unique(CENfeats_CEN180_metrics_allacc$accession)), " accessions"),
                                      dataFrame = CENfeats_CEN180_metrics_allacc,
                                      mapping = aes(x = Feature,
                                                    y = HORcount + 1,
                                                    fill = Family),
                                      pvals = paste0(Utest_PvalsChar_HORcount, collapse = "                "),
                                      xlab = bquote(.(flankNamePlot) ~ "regions flanking centromeric" ~ italic("ATHILA") ~ "and random loci"),
                                      ylab = bquote(italic("CEN180") ~ "HOR count"),
                                      yaxtrans = log2_trans(),
                                      ybreaks = trans_breaks("log2", function(x) 2^x),
                                      ylabels = trans_format("log2", math_format(2^.x)))
ggPoint_HORcount_allacc <- ggPoint_HORcount_allacc +
  facet_grid(cols = vars(chr), rows = vars(Family), margins = "chr", scales = "fixed")


# WeightedConsensusScore
Utest_Pvals_WeightedConsensusScore <- as.numeric(sapply(acc_chrName, function(z) {
  print(z)
  if(z != "All") {
    tt <- tryCatch(
                   wilcox.test(x = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$chr == z &
                                                                         CENfeats_CEN180_metrics_allacc$Feature == "CENATHILA",],
                                          "WeightedConsensusScore")[,1],
                               y = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$chr == z &
                                                                         CENfeats_CEN180_metrics_allacc$Feature == "CENranLoc",],
                                          "WeightedConsensusScore")[,1],
                               alternative = "two.sided")$p.value,
                   error = function(e) e
                  )
    if(is(tt, "error")) {
      NA
    } else {
      tt
    }
  } else {
    tt <- tryCatch(
                   wilcox.test(x = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Feature == "CENATHILA",],
                                          "WeightedConsensusScore")[,1],
                               y = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Feature == "CENranLoc",],
                                          "WeightedConsensusScore")[,1],
                               alternative = "two.sided")$p.value,
                   error = function(e) e
                  )
    if(is(tt, "error")) {
      NA
    } else {
      tt
    }
  }
}))
#Utest_Pvals_WeightedConsensusScore <- sapply(1:length(Utest_list_WeightedConsensusScore), function(z) { Utest_list_WeightedConsensusScore[[z]]$p.value } )
Utest_PvalsChar_WeightedConsensusScore <- sapply(1:length(Utest_Pvals_WeightedConsensusScore), function(z) {
  if(is.na(Utest_Pvals_WeightedConsensusScore[z])) {
    paste0(acc_chrName[z], " MWW P = NA")
  } else if(Utest_Pvals_WeightedConsensusScore[z] < 0.0001) {
    paste0(acc_chrName[z], " MWW P < 0.0001")
  } else {
    paste0(acc_chrName[z], " MWW P = ", as.character(round(Utest_Pvals_WeightedConsensusScore[z], digits = 4)))
  }
})

ggPoint_WeightedConsensusScore_allacc <- violinPlot(acc_id = paste0(length(unique(CENfeats_CEN180_metrics_allacc$accession)), " accessions"),
                                                    dataFrame = CENfeats_CEN180_metrics_allacc,
                                                    mapping = aes(x = Feature,
                                                                  y = WeightedConsensusScore + 1,
                                                                  fill = Family),
                                                    pvals = paste0(Utest_PvalsChar_WeightedConsensusScore, collapse = "                "),
                                                    xlab = bquote(.(flankNamePlot) ~ "regions flanking centromeric" ~ italic("ATHILA") ~ "and random loci"),
                                                    ylab = bquote(italic("CEN180") ~ "consensus score"),
                                                    yaxtrans = log2_trans(),
                                                    ybreaks = trans_breaks("log2", function(x) 2^x),
                                                    ylabels = trans_format("log2", math_format(2^.x)))
ggPoint_WeightedConsensusScore_allacc <- ggPoint_WeightedConsensusScore_allacc +
  facet_grid(cols = vars(chr), rows = vars(Family), margins = "chr", scales = "fixed")


# EditDistance
Utest_Pvals_EditDistance <- as.numeric(sapply(acc_chrName, function(z) {
  print(z)
  if(z != "All") {
    tt <- tryCatch(
                   wilcox.test(x = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$chr == z &
                                                                         CENfeats_CEN180_metrics_allacc$Feature == "CENATHILA",],
                                          "EditDistance")[,1],
                               y = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$chr == z &
                                                                         CENfeats_CEN180_metrics_allacc$Feature == "CENranLoc",],
                                          "EditDistance")[,1],
                               alternative = "two.sided")$p.value,
                   error = function(e) e
                  )
    if(is(tt, "error")) {
      NA
    } else {
      tt
    }
  } else {
    tt <- tryCatch(
                   wilcox.test(x = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Feature == "CENATHILA",],
                                          "EditDistance")[,1],
                               y = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Feature == "CENranLoc",],
                                          "EditDistance")[,1],
                               alternative = "two.sided")$p.value,
                   error = function(e) e
                  )
    if(is(tt, "error")) {
      NA
    } else {
      tt
    }
  }
}))
#Utest_Pvals_EditDistance <- sapply(1:length(Utest_list_EditDistance), function(z) { Utest_list_EditDistance[[z]]$p.value } )
Utest_PvalsChar_EditDistance <- sapply(1:length(Utest_Pvals_EditDistance), function(z) {
  if(is.na(Utest_Pvals_EditDistance[z])) {
    paste0(acc_chrName[z], " MWW P = NA")
  } else if(Utest_Pvals_EditDistance[z] < 0.0001) {
    paste0(acc_chrName[z], " MWW P < 0.0001")
  } else {
    paste0(acc_chrName[z], " MWW P = ", as.character(round(Utest_Pvals_EditDistance[z], digits = 4)))
  }
})

ggPoint_EditDistance_allacc <- violinPlot(acc_id = paste0(length(unique(CENfeats_CEN180_metrics_allacc$accession)), " accessions"),
                                          dataFrame = CENfeats_CEN180_metrics_allacc,
                                          mapping = aes(x = Feature,
                                                        y = EditDistance + 1,
                                                        fill = Family),
                                          pvals = paste0(Utest_PvalsChar_EditDistance, collapse = "                "),
                                          xlab = bquote(.(flankNamePlot) ~ "regions flanking centromeric" ~ italic("ATHILA") ~ "and random loci"),
                                          ylab = bquote(italic("CEN180") ~ "edit distance"),
                                          yaxtrans = log2_trans(),
                                          ybreaks = trans_breaks("log2", function(x) 2^x),
                                          ylabels = trans_format("log2", math_format(2^.x)))
ggPoint_EditDistance_allacc <- ggPoint_EditDistance_allacc +
  facet_grid(cols = vars(chr), rows = vars(Family), margins = "chr", scales = "fixed")

# Plot as one with facet rows as Family
gg_cow_list <- list(
                    ggPoint_HORlengthsSum_allacc,
                    ggPoint_HORcount_allacc,
                    ggPoint_WeightedConsensusScore_allacc,
                    ggPoint_EditDistance_allacc
                   )
gg_cow <- plot_grid(plotlist = gg_cow_list,
                    labels = "AUTO", label_size = 30,
                    align = "hv",
                    axis = "l",
                    nrow = length(gg_cow_list), ncol = 1)
ggsave(paste0(plotDirAllAccessions,
              "CENATHILA_CEN180_", flankName, "_flanks_AllMetrics_violinPlot_",
              paste0(chrName, collapse = "_"),
              "_", length(unique(CENfeats_CEN180_metrics_allacc$accession)), "accessions.pdf"),
       plot = gg_cow,
       height = (2*length(phylo))*length(gg_cow_list), width = 5*(length(acc_chrName)), limitsize = F)



# Plot across all accessions, with gg_cow rows as Family 

# HORlengthsSum
Utest_Pvals_HORlengthsSum_phylo <- lapply(1:length(phylo), function(x) {
  as.numeric(sapply(acc_chrName, function(z) {
    print(z)
    if(z != "All") {
      tt <- tryCatch(
                     wilcox.test(x = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$chr == z &
                                                                           CENfeats_CEN180_metrics_allacc$Family == phylo[x] &
                                                                           CENfeats_CEN180_metrics_allacc$Feature == "CENATHILA",],
                                            "HORlengthsSum")[,1],
                                 y = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$chr == z &
                                                                           CENfeats_CEN180_metrics_allacc$Family == phylo[x] &
                                                                           CENfeats_CEN180_metrics_allacc$Feature == "CENranLoc",],
                                            "HORlengthsSum")[,1],
                                 alternative = "two.sided")$p.value,
                     error = function(e) e
                    )
      if(is(tt, "error")) {
        NA
      } else {
        tt
      }
    } else {
      tt <- tryCatch(
                     wilcox.test(x = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Family == phylo[x] &
                                                                           CENfeats_CEN180_metrics_allacc$Feature == "CENATHILA",],
                                            "HORlengthsSum")[,1],
                                 y = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Family == phylo[x] &
                                                                           CENfeats_CEN180_metrics_allacc$Feature == "CENranLoc",],
                                            "HORlengthsSum")[,1],
                                 alternative = "two.sided")$p.value,
                     error = function(e) e
                    )
      if(is(tt, "error")) {
        NA
      } else {
        tt
      }
    }
  }))
})
Utest_PvalsChar_HORlengthsSum_phylo <- lapply(1:length(phylo), function(x) {
  sapply(1:length(Utest_Pvals_HORlengthsSum_phylo[[x]]), function(z) {
    if(is.na(Utest_Pvals_HORlengthsSum_phylo[[x]][z])) {
      paste0(acc_chrName[z], " MWW P = NA")
    } else if(Utest_Pvals_HORlengthsSum_phylo[[x]][z] < 0.0001) {
      paste0(acc_chrName[z], " MWW P < 0.0001")
    } else {
      paste0(acc_chrName[z], " MWW P = ", as.character(round(Utest_Pvals_HORlengthsSum_phylo[[x]][z], digits = 4)))
    }
  })
})

ggPoint_HORlengthsSum_allacc_phylo <- lapply(1:length(phylo), function(x) {
  tP <- violinPlot(acc_id = paste0(length(unique(CENfeats_CEN180_metrics_allacc$accession)), " accessions"),
                   dataFrame = CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Family == phylo[x],],
                   mapping = aes(x = Feature,
                                 y = HORlengthsSum + 1,
                                 fill = Feature),
                   pvals = paste0(Utest_PvalsChar_HORlengthsSum_phylo[[x]], collapse = "                "),
                   xlab = bquote(.(flankNamePlot) ~ "regions flanking centromeric" ~ italic(.(phylo[x])) ~ "and random loci"),
                   ylab = bquote(italic("CEN180") ~ "repetitiveness"),
                   yaxtrans = log2_trans(),
                   ybreaks = trans_breaks("log2", function(x) 2^x),
                   ylabels = trans_format("log2", math_format(2^.x)))
  tP <- tP +
    facet_grid(cols = vars(chr), margins = "chr", scales = "fixed")
})

# Plot with gg_cow rows as Family
gg_cow_list <- ggPoint_HORlengthsSum_allacc_phylo
gg_cow <- plot_grid(plotlist = gg_cow_list,
                    labels = "AUTO", label_size = 30,
                    align = "hv",
                    axis = "l",
                    nrow = length(gg_cow_list), ncol = 1)
ggsave(paste0(plotDirAllAccessions,
              "CENATHILA_CEN180_", flankName, "_flanks_HORlengthsSum_violinPlot_",
              paste0(chrName, collapse = "_"),
              "_", length(unique(CENfeats_CEN180_metrics_allacc$accession)), "accessions.pdf"),
       plot = gg_cow,
       height = 5.5*length(gg_cow_list), width = 5*(length(acc_chrName)), limitsize = F)


# HORcount
Utest_Pvals_HORcount_phylo <- lapply(1:length(phylo), function(x) {
  as.numeric(sapply(acc_chrName, function(z) {
    print(z)
    if(z != "All") {
      tt <- tryCatch(
                     wilcox.test(x = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$chr == z &
                                                                           CENfeats_CEN180_metrics_allacc$Family == phylo[x] &
                                                                           CENfeats_CEN180_metrics_allacc$Feature == "CENATHILA",],
                                            "HORcount")[,1],
                                 y = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$chr == z &
                                                                           CENfeats_CEN180_metrics_allacc$Family == phylo[x] &
                                                                           CENfeats_CEN180_metrics_allacc$Feature == "CENranLoc",],
                                            "HORcount")[,1],
                                 alternative = "two.sided")$p.value,
                     error = function(e) e
                    )
      if(is(tt, "error")) {
        NA
      } else {
        tt
      }
    } else {
      tt <- tryCatch(
                     wilcox.test(x = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Family == phylo[x] &
                                                                           CENfeats_CEN180_metrics_allacc$Feature == "CENATHILA",],
                                            "HORcount")[,1],
                                 y = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Family == phylo[x] &
                                                                           CENfeats_CEN180_metrics_allacc$Feature == "CENranLoc",],
                                            "HORcount")[,1],
                                 alternative = "two.sided")$p.value,
                     error = function(e) e
                    )
      if(is(tt, "error")) {
        NA
      } else {
        tt
      }
    }
  }))
})
Utest_PvalsChar_HORcount_phylo <- lapply(1:length(phylo), function(x) {
  sapply(1:length(Utest_Pvals_HORcount_phylo[[x]]), function(z) {
    if(is.na(Utest_Pvals_HORcount_phylo[[x]][z])) {
      paste0(acc_chrName[z], " MWW P = NA")
    } else if(Utest_Pvals_HORcount_phylo[[x]][z] < 0.0001) {
      paste0(acc_chrName[z], " MWW P < 0.0001")
    } else {
      paste0(acc_chrName[z], " MWW P = ", as.character(round(Utest_Pvals_HORcount_phylo[[x]][z], digits = 4)))
    }
  })
})

ggPoint_HORcount_allacc_phylo <- lapply(1:length(phylo), function(x) {
  tP <- violinPlot(acc_id = paste0(length(unique(CENfeats_CEN180_metrics_allacc$accession)), " accessions"),
                   dataFrame = CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Family == phylo[x],],
                   mapping = aes(x = Feature,
                                 y = HORcount + 1,
                                 fill = Feature),
                   pvals = paste0(Utest_PvalsChar_HORcount_phylo[[x]], collapse = "                "),
                   xlab = bquote(.(flankNamePlot) ~ "regions flanking centromeric" ~ italic(.(phylo[x])) ~ "and random loci"),
                   ylab = bquote(italic("CEN180") ~ "HOR count"),
                   yaxtrans = log2_trans(),
                   ybreaks = trans_breaks("log2", function(x) 2^x),
                   ylabels = trans_format("log2", math_format(2^.x)))
  tP <- tP +
    facet_grid(cols = vars(chr), margins = "chr", scales = "fixed")
})

# Plot with gg_cow rows as Family
gg_cow_list <- ggPoint_HORcount_allacc_phylo
gg_cow <- plot_grid(plotlist = gg_cow_list,
                    labels = "AUTO", label_size = 30,
                    align = "hv",
                    axis = "l",
                    nrow = length(gg_cow_list), ncol = 1)
ggsave(paste0(plotDirAllAccessions,
              "CENATHILA_CEN180_", flankName, "_flanks_HORcount_violinPlot_",
              paste0(chrName, collapse = "_"),
              "_", length(unique(CENfeats_CEN180_metrics_allacc$accession)), "accessions.pdf"),
       plot = gg_cow,
       height = 5.5*length(gg_cow_list), width = 5*(length(acc_chrName)), limitsize = F)


# WeightedConsensusScore
Utest_Pvals_WeightedConsensusScore_phylo <- lapply(1:length(phylo), function(x) {
  as.numeric(sapply(acc_chrName, function(z) {
    print(z)
    if(z != "All") {
      tt <- tryCatch(
                     wilcox.test(x = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$chr == z &
                                                                           CENfeats_CEN180_metrics_allacc$Family == phylo[x] &
                                                                           CENfeats_CEN180_metrics_allacc$Feature == "CENATHILA",],
                                            "WeightedConsensusScore")[,1],
                                 y = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$chr == z &
                                                                           CENfeats_CEN180_metrics_allacc$Family == phylo[x] &
                                                                           CENfeats_CEN180_metrics_allacc$Feature == "CENranLoc",],
                                            "WeightedConsensusScore")[,1],
                                 alternative = "two.sided")$p.value,
                     error = function(e) e
                    )
      if(is(tt, "error")) {
        NA
      } else {
        tt
      }
    } else {
      tt <- tryCatch(
                     wilcox.test(x = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Family == phylo[x] &
                                                                           CENfeats_CEN180_metrics_allacc$Feature == "CENATHILA",],
                                            "WeightedConsensusScore")[,1],
                                 y = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Family == phylo[x] &
                                                                           CENfeats_CEN180_metrics_allacc$Feature == "CENranLoc",],
                                            "WeightedConsensusScore")[,1],
                                 alternative = "two.sided")$p.value,
                     error = function(e) e
                    )
      if(is(tt, "error")) {
        NA
      } else {
        tt
      }
    }
  }))
})
Utest_PvalsChar_WeightedConsensusScore_phylo <- lapply(1:length(phylo), function(x) {
  sapply(1:length(Utest_Pvals_WeightedConsensusScore_phylo[[x]]), function(z) {
    if(is.na(Utest_Pvals_WeightedConsensusScore_phylo[[x]][z])) {
      paste0(acc_chrName[z], " MWW P = NA")
    } else if(Utest_Pvals_WeightedConsensusScore_phylo[[x]][z] < 0.0001) {
      paste0(acc_chrName[z], " MWW P < 0.0001")
    } else {
      paste0(acc_chrName[z], " MWW P = ", as.character(round(Utest_Pvals_WeightedConsensusScore_phylo[[x]][z], digits = 4)))
    }
  })
})

ggPoint_WeightedConsensusScore_allacc_phylo <- lapply(1:length(phylo), function(x) {
  tP <- violinPlot(acc_id = paste0(length(unique(CENfeats_CEN180_metrics_allacc$accession)), " accessions"),
                   dataFrame = CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Family == phylo[x],],
                   mapping = aes(x = Feature,
                                 y = WeightedConsensusScore + 1,
                                 fill = Feature),
                   pvals = paste0(Utest_PvalsChar_WeightedConsensusScore_phylo[[x]], collapse = "                "),
                   xlab = bquote(.(flankNamePlot) ~ "regions flanking centromeric" ~ italic(.(phylo[x])) ~ "and random loci"),
                   ylab = bquote(italic("CEN180") ~ "consensus score"),
                   yaxtrans = log2_trans(),
                   ybreaks = trans_breaks("log2", function(x) 2^x),
                   ylabels = trans_format("log2", math_format(2^.x)))
  tP <- tP +
    facet_grid(cols = vars(chr), margins = "chr", scales = "fixed")
})

# Plot with gg_cow rows as Family
gg_cow_list <- ggPoint_WeightedConsensusScore_allacc_phylo
gg_cow <- plot_grid(plotlist = gg_cow_list,
                    labels = "AUTO", label_size = 30,
                    align = "hv",
                    axis = "l",
                    nrow = length(gg_cow_list), ncol = 1)
ggsave(paste0(plotDirAllAccessions,
              "CENATHILA_CEN180_", flankName, "_flanks_WeightedConsensusScore_violinPlot_",
              paste0(chrName, collapse = "_"),
              "_", length(unique(CENfeats_CEN180_metrics_allacc$accession)), "accessions.pdf"),
       plot = gg_cow,
       height = 5.5*length(gg_cow_list), width = 5*(length(acc_chrName)), limitsize = F)


# EditDistance
Utest_Pvals_EditDistance_phylo <- lapply(1:length(phylo), function(x) {
  as.numeric(sapply(acc_chrName, function(z) {
    print(z)
    if(z != "All") {
      tt <- tryCatch(
                     wilcox.test(x = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$chr == z &
                                                                           CENfeats_CEN180_metrics_allacc$Family == phylo[x] &
                                                                           CENfeats_CEN180_metrics_allacc$Feature == "CENATHILA",],
                                            "EditDistance")[,1],
                                 y = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$chr == z &
                                                                           CENfeats_CEN180_metrics_allacc$Family == phylo[x] &
                                                                           CENfeats_CEN180_metrics_allacc$Feature == "CENranLoc",],
                                            "EditDistance")[,1],
                                 alternative = "two.sided")$p.value,
                     error = function(e) e
                    )
      if(is(tt, "error")) {
        NA
      } else {
        tt
      }
    } else {
      tt <- tryCatch(
                     wilcox.test(x = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Family == phylo[x] &
                                                                           CENfeats_CEN180_metrics_allacc$Feature == "CENATHILA",],
                                            "EditDistance")[,1],
                                 y = select(CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Family == phylo[x] &
                                                                           CENfeats_CEN180_metrics_allacc$Feature == "CENranLoc",],
                                            "EditDistance")[,1],
                                 alternative = "two.sided")$p.value,
                     error = function(e) e
                    )
      if(is(tt, "error")) {
        NA
      } else {
        tt
      }
    }
  }))
})
Utest_PvalsChar_EditDistance_phylo <- lapply(1:length(phylo), function(x) {
  sapply(1:length(Utest_Pvals_EditDistance_phylo[[x]]), function(z) {
    if(is.na(Utest_Pvals_EditDistance_phylo[[x]][z])) {
      paste0(acc_chrName[z], " MWW P = NA")
    } else if(Utest_Pvals_EditDistance_phylo[[x]][z] < 0.0001) {
      paste0(acc_chrName[z], " MWW P < 0.0001")
    } else {
      paste0(acc_chrName[z], " MWW P = ", as.character(round(Utest_Pvals_EditDistance_phylo[[x]][z], digits = 4)))
    }
  })
})

ggPoint_EditDistance_allacc_phylo <- lapply(1:length(phylo), function(x) {
  tP <- violinPlot(acc_id = paste0(length(unique(CENfeats_CEN180_metrics_allacc$accession)), " accessions"),
                   dataFrame = CENfeats_CEN180_metrics_allacc[CENfeats_CEN180_metrics_allacc$Family == phylo[x],],
                   mapping = aes(x = Feature,
                                 y = EditDistance + 1,
                                 fill = Feature),
                   pvals = paste0(Utest_PvalsChar_EditDistance_phylo[[x]], collapse = "                "),
                   xlab = bquote(.(flankNamePlot) ~ "regions flanking centromeric" ~ italic(.(phylo[x])) ~ "and random loci"),
                   ylab = bquote(italic("CEN180") ~ "edit distance"),
                   yaxtrans = log2_trans(),
                   ybreaks = trans_breaks("log2", function(x) 2^x),
                   ylabels = trans_format("log2", math_format(2^.x)))
  tP <- tP +
    facet_grid(cols = vars(chr), margins = "chr", scales = "fixed")
})

# Plot with gg_cow rows as Family
gg_cow_list <- ggPoint_EditDistance_allacc_phylo
gg_cow <- plot_grid(plotlist = gg_cow_list,
                    labels = "AUTO", label_size = 30,
                    align = "hv",
                    axis = "l",
                    nrow = length(gg_cow_list), ncol = 1)
ggsave(paste0(plotDirAllAccessions,
              "CENATHILA_CEN180_", flankName, "_flanks_EditDistance_violinPlot_",
              paste0(chrName, collapse = "_"),
              "_", length(unique(CENfeats_CEN180_metrics_allacc$accession)), "accessions.pdf"),
       plot = gg_cow,
       height = 5.5*length(gg_cow_list), width = 5*(length(acc_chrName)), limitsize = F)

