#!/usr/bin/env Rscript

# Plot across all accessions comparison of average CEN180 metrics (HOR membership and divergence)
# in regions flanking centromeric ATHILA and matched randomly positioned centromeric loci

# Usage:
# conda activate R-4.1.2
# ./66Atha_CENATHILA_flanks_CEN180_metrics_v260522_hpc_PLOT_ONLY.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 1000 1e4 Flanks
# conda deactivate

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#flankSize <- 1000
#perms <- 1e4
#regionName <- "Flanks"

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
flankSize <- as.numeric(args[2])
perms <- as.numeric(args[3])
regionName <- args[4]

# Set minimum possible P-value for permutation test result with
# perms sets of random loci
minPval <- 1 - ( (perms - 1) / perms)

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
suppressMessages(library(data.table, quietly = T))
suppressMessages(library(GenomicRanges, quietly = T))
suppressMessages(library(segmentSeq, quietly = T))
suppressMessages(library(parallel, quietly = T))
suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(scales, quietly = T))
suppressMessages(library(pals, quietly = T))
suppressMessages(library(RColorBrewer, quietly = T))
suppressMessages(library(ggplot2, quietly = T))
suppressMessages(library(cowplot, quietly = T))
suppressMessages(library(ggbeeswarm, quietly = T))
suppressMessages(library(viridis, quietly = T))

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

combined <- fread(file = paste0(plotDirAllMetrics,
                                "CENATHILA_", flankName, "_regions_CEN180_AllMetrics_combined_",
                                paste0(chrName, collapse = "_"), "_",
                                perms, "perms_",
                                regionName,
                                ".tsv"),
                  data.table = F)

combined$Accession <- factor(combined$Accession,
                             levels = c(paste0(length(unique(combined$Accession))-1, " accessions"),
                                        sort(unique(combined$Accession))[-grep(" accessions", sort(unique(combined$Accession)))]))
combined$Metric <- factor(combined$Metric,
                          levels = unique(combined$Metric))
combined$Region <- factor(combined$Region,
                          levels = unique(combined$Region))
combined$Family <- factor(combined$Family,
                          levels = sort(unique(combined$Family)))

# Per-accession plots (including All accessions plot)
bp_per <- ggplot(data = combined,
                 mapping = aes(x = Family,
                               y = log2obsexp,
                               fill = Metric)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
#  scale_fill_discrete(name = "") +
  scale_fill_brewer(name = "",
                    palette = "Dark2") +
  geom_point(mapping = aes(x = Family,
                           y = log2alpha),
             position = position_dodge(0.9),
             shape = "-", colour  = "grey70", size = 12) +
  xlab(bquote(italic("CENATHILA") ~ .(flankNamePlot) ~ "flanking regions")) +
  ylab(bquote("Log"[2] * "(observed/expected)" ~ italic("CEN178") ~ "metric")) +
#  scale_y_continuous(limits = c(-4.0, 4.0)) +
  scale_x_discrete(position = "bottom") +
  guides(fill = guide_legend(direction = "vertical",
                             label.position = "right",
                             label.theme = element_text(size = 16, hjust = 0, vjust = 0.5, angle = 0),
                             nrow = length(unique(combined$Family)),
                             byrow = TRUE)) +

  theme_bw() +
  theme(axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks.y = element_line(size = 0.5, colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black", hjust = 0.5, vjust = 0.5, angle = 90),
        axis.title.y = element_text(size = 20, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 20, colour = "black", hjust = 0.5, vjust = 0.5, angle = 45),
        axis.title.x = element_text(size = 20, colour = "black", margin = margin(t = 20)),
        strip.text.x = element_text(size = 20, colour = "white"),
        strip.text.y = element_text(size = 20, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"),
#        panel.grid = element_blank(),
#        panel.border = element_blank(),
#        panel.background = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        plot.title = element_text(size = 20, colour = "black", hjust = 0.5)) +
  ggtitle(bquote("Permutations =" ~
                 .(prettyNum(perms,
                             big.mark = ",",
                             trim = T)) ~ "sets of randomly positioned centromeric loci"))
bp_per <- bp_per +
  facet_grid(cols = vars(Accession), rows = vars(Region), scales = "fixed")
ggsave(paste0(plotDirAllMetrics,
              "CENATHILA_", flankName, "_regions_CEN180_AllMetrics_combined_bargraph_",
              paste0(chrName, collapse = "_"), "_",
              perms, "perms_",
              paste0(unique(as.vector(combined$Region)), collapse = "_"),
              ".pdf"),
       plot = bp_per,
       width = 10*length(unique(as.vector(combined$Accession))), height = 8, limitsize = F)

# All accessions plot (excluding per-accession plots)
all_accessions <- combined[grep(" accessions", combined$Accession),] 
bp_all <- ggplot(data = all_accessions,
                 mapping = aes(x = Family,
                               y = log2obsexp,
                               fill = Metric)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
#  scale_fill_discrete(name = "") +
  scale_fill_brewer(name = "",
                    palette = "Dark2") +
  geom_point(mapping = aes(x = Family,
                           y = log2alpha),
             position = position_dodge(0.9),
             shape = "-", colour  = "grey70", size = 12) +
  xlab(bquote(italic("CENATHILA") ~ .(flankNamePlot) ~ "flanking regions")) +
  ylab(bquote("Log"[2] * "(observed/expected)" ~ italic("CEN178") ~ "metric")) +
#  scale_y_continuous(limits = c(-4.0, 4.0)) +
  scale_x_discrete(position = "bottom") +
  guides(fill = guide_legend(direction = "vertical",
                             label.position = "right",
                             label.theme = element_text(size = 16, hjust = 0, vjust = 0.5, angle = 0),
                             nrow = length(unique(all_accessions$Family)),
                             byrow = TRUE)) +

  theme_bw() +
  theme(axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks.y = element_line(size = 0.5, colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black", hjust = 0.5, vjust = 0.5, angle = 90),
        axis.title.y = element_text(size = 20, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 20, colour = "black", hjust = 0.5, vjust = 0.5, angle = 45),
        axis.title.x = element_text(size = 20, colour = "black", margin = margin(t = 20)),
        strip.text.x = element_text(size = 20, colour = "white"),
        strip.text.y = element_text(size = 20, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"),
#        panel.grid = element_blank(),
#        panel.border = element_blank(),
#        panel.background = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        plot.title = element_text(size = 20, colour = "black", hjust = 0.5)) +
  ggtitle(bquote("Permutations =" ~
                 .(prettyNum(perms,
                             big.mark = ",",
                             trim = T)) ~ "sets of randomly positioned centromeric loci"))
bp_all <- bp_all +
  facet_grid(cols = vars(Accession), scales = "fixed")
ggsave(paste0(plotDirAllMetrics,
              "CENATHILA_", flankName, "_regions_CEN180_AllMetrics_combined_all_accessions_bargraph_",
              paste0(chrName, collapse = "_"), "_",
              perms, "perms_",
              paste0(unique(as.vector(all_accessions$Region)), collapse = "_"),
              ".pdf"),
       plot = bp_all,
       width = 14*length(unique(as.vector(all_accessions$Accession))), height = 8, limitsize = F)
