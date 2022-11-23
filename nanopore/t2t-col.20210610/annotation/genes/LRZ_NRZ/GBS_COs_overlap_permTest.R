#!/usr/bin/env Rscript

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 23/11/2022

# Do permutation tests to evaluate overlap between WGS-derived crossovers
# and features of interest

# Usage:
# conda activate python_3.9.6
# ./GBS_COs_overlap_permTest.R t2t-col.20210610 'Chr1,Chr2,Chr3,Chr4,Chr5' 10000 LRZs
# conda deactivate

args = commandArgs(trailingOnly=T)
alnTo = args[1]
chrom = unlist(strsplit(args[2], split=","))
perms = as.integer(args[3])
region = args[4]

#alnTo = "t2t-col.20210610"
#chrom = unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                        split=","))
#perms = as.integer(10000)
#region = "LRZs"

# Set minimum possible P-value for permutation test result with
# perms sets of random loci
min_pval = 1 - ( (perms - 1) / perms)


library(regioneR)
library(rtracklayer)
library(dplyr)
library(ggplot2)


outdir = "perm_tests/"
plotdir = paste0(outdir, "plots/")
system(paste0("[ -d ", plotdir, " ] || mkdir -p ", plotdir))



# Genomic definitions
fai = read.table(paste0("/rds/user/ajt200/hpc-work/nanopore/", alnTo, "/", alnTo, ".fa.fai"), header=F)
chrs = fai[ which(fai$V1 %in% chrom), ]$V1
chrLens = fai[ which(fai$V1 %in% chrom), ]$V2

genome_GR = GRanges(seqnames=chrs,
                    ranges=IRanges(start=rep(1, length(chrs)),
                                   end=chrLens),
                    strand="*")

# Region definitions
regions_DF = read.csv("Table_S1_LRZ_NRZ.csv", header=T)[,-1]
regions_DF = regions_DF[ which(regions_DF$Chr %in% chrs), ]
regions_DF = regions_DF[ with(regions_DF, order(Chr, decreasing=F)), ]

LRZstart = regions_DF$LRZ.left
LRZend = regions_DF$LRZ.right
NRZstart = regions_DF$NRZ.left
NRZend = regions_DF$NRZ.right

LRZ_GR = GRanges(seqnames=rep(regions_DF$Chr, 2),
                 ranges=IRanges(start=c(LRZstart, NRZend+1),
                                end=c(NRZstart-1, LRZend)),
                 strand="*")

NRZ_GR = GRanges(seqnames=regions_DF$Chr,
                 ranges=IRanges(start=NRZstart,
                                end=NRZend),
                 strand="*")

arm_GR = GRanges(seqnames=rep(regions_DF$Chr, 2),
                 ranges=IRanges(start=c(rep(1, nrow(regions_DF)),
                                      LRZend+1),
                                end=c(LRZstart-1,
                                      chrLens)),
                 strand="*")


# COs
COs_DF = read.table("cos.all.ajt.txt", header=T)
COs_GR = GRanges(seqnames=COs_DF$Chr_id,
                 ranges=IRanges(start=COs_DF$COs_start,
                                end=COs_DF$COs_stop),
                 strand="*")
COs_GR = sortSeqlevels(COs_GR)
COs_GR = sort(COs_GR)


# Define region to be masked out of analysis
if(region == "LRZs") {
    region_GR = LRZ_GR
    mask_GR = c(arm_GR, NRZ_GR)
    region_plot = "LRZ"
} else if(region == "arms") {
    region_GR = arm_GR
    mask_GR = c(LRZ_GR, NRZ_GR)
    region_plot = "Arm"
} else if(region == "genomewide") {
    region = genome_GR()
    mask_GR = GRanges()
    region_plot = "Genome-wide"
} else {
    stop("region is not 'LRZs' 'arms' or 'genomewide'")
}


# Subset to include only those not overlapping masked region (e.g., centromere)
mask_COs_overlap = findOverlaps(query=mask_GR,
                                subject=COs_GR,
                                type = "any",
                                select = "all",
                                ignore.strand = TRUE)
if(length(mask_COs_overlap) > 0) {
    COs_GR = COs_GR[-subjectHits(mask_COs_overlap)]
}
print("Candidate recombination events:")
print(COs_GR)


# Load table of gene coordinates in t2t-col.20210610 (GFF3 derived from Liftoff tool)
genes = readGFF(paste0("/rds/user/ajt200/hpc-work/nanopore/", alnTo, "/annotation/genes/", alnTo, ".genes.gff3"))
genes = genes[ which(genes$seqid %in% chrs), ]
genes = genes[ which(genes$type == "gene"), ]
print(dim(genes))
#[1] 28504    20

genes_GR = GRanges(seqnames=genes$seqid,
                   ranges=IRanges(start=genes$start,
                                  end=genes$end),
                   strand=genes$strand)

# Subset to include only those not overlapping masked region
mask_genes_overlap = findOverlaps(query = mask_GR,
                                  subject = genes_GR,
                                  type = "any",
                                  select = "all",
                                  ignore.strand = TRUE)
if(length(mask_genes_overlap) > 0) {
    genes_GR = genes_GR[-subjectHits(mask_genes_overlap)]
}

genes_GR = unique(genes_GR)

# NOTE: Retain strand information until after obtaining promoters, etc.
# Overlap analysis should be strand-unaware

# Obtain 1000-bp gene promoters
promoters_GR = promoters(genes_GR, upstream=1000, downstream=0)
promoters_GR = unique(promoters_GR)
strand(promoters_GR) = "*"
print(promoters_GR)

# Obtain regions immediately downstream of gene TSSs (gene 5' ends: TSS to TSS+499 bp)
g5ends_GR = promoters(genes_GR, upstream=0, downstream=500)
g5ends_GR = unique(g5ends_GR)
strand(g5ends_GR) = "*"
print(g5ends_GR)

# Obtain regions relative to TTS using TTSplus()
TTSplus = function(x, upstream=100, downstream=1000, ...) {
    if(!isSingleNumber(upstream))
        stop("'upstream' must be a single integer")
    if(!is.integer(upstream))
        upstream = as.numeric(upstream)
    if(!isSingleNumber(downstream))
        stop("'downstream' must be a single integer")
    if(!is.integer(downstream))
        downstream = as.numeric(downstream)
    if(downstream < 0)
        stop("'downstream' must be an integer >= 0")
    if(any(strand(x) == "*"))
        warning("'*' ranges were treated as '+'")
    on_plus = which(strand(x) == "+" | strand(x) == "*")
    on_plus_TTS = end(x)[on_plus]
    start(x)[on_plus] = on_plus_TTS - upstream
    end(x)[on_plus] = on_plus_TTS + downstream
    on_minus = which(strand(x) == "-")
    on_minus_TTS = start(x)[on_minus]
    end(x)[on_minus] = on_minus_TTS + upstream
    start(x)[on_minus] = on_minus_TTS - downstream
    return(x)
}

# Obtain regions immediately upstream of gene TTSs (gene 3' ends: TTS to TTS-499 bp)
g3ends_GR = TTSplus(genes_GR, upstream=499, downstream=0)
g3ends_GR = unique(g3ends_GR)
strand(g3ends_GR) = "*"
print(g3ends_GR)

# Obtain 1000-bp gene terminators
terminators_GR = TTSplus(genes_GR, upstream=-1, downstream=1000)
terminators_GR = unique(terminators_GR)
strand(terminators_GR) = "*"
print(terminators_GR)

# Remove strand information from genes_GR
strand(genes_GR) = "*"
print(genes_GR)

# Feature names
features_names = c(
                   "Genes",
                   "1 kb upstream of TSSs",                  
                   "500 bp downstream of TSSs",
                   "500 bp upstream of TTSs",
                   "1 kb downstream of TTSs"
                  )
# GRanges list
features_GR_list = c(
                     "Genes"=genes_GR,
                     "1 kb upstream of TSSs"=promoters_GR,
                     "500 bp downstream of TSSs"=g5ends_GR,
                     "500 bp upstream of TTSs"=g3ends_GR,
                     "1 kb downstream of TTSs"=terminators_GR
                    )


# Perform permutation tests with randomized regions generated on a per chromosome basis
set.seed(47393573)
pt_COs_vs_features = lapply(1:length(features_GR_list), function(x) {
    permTest(A=COs_GR,
             B=features_GR_list[[x]],
             alternative="auto",
             ntimes=perms,
             randomize.function=randomizeRegions,
             genome=genome_GR,
             mask=mask_GR,
             allow.overlaps=TRUE,
             per.chromosome=TRUE,
             evaluate.function=numOverlaps,
             count.once=TRUE,
             mc.set.seed=FALSE,
             mc.cores=detectCores())
})



pt_dist_DF = data.frame(
                        Region = unlist(lapply(1:length(pt_COs_vs_features), function(x) {
                            rep(region_plot,
                                times=length(pt_COs_vs_features[[x]]$numOverlaps$permuted))
                        })),
                        Feature = unlist(lapply(1:length(pt_COs_vs_features), function(x) {
                            rep(names(features_GR_list)[x],
                                times=length(pt_COs_vs_features[[x]]$numOverlaps$permuted))
                        })),
                        Permuted = unlist(lapply(1:length(pt_COs_vs_features), function(x) {
                            pt_COs_vs_features[[x]]$numOverlaps$permuted
                        }))
                      )
                        

pt_DF = data.frame(
                   Region = rep(region_plot, times=length(names(features_GR_list))),
                   Feature = names(features_GR_list),
                   Number_of_features = unlist(lapply(1:length(features_GR_list), function(x) {
                       length(features_GR_list[[x]])
                   })),
                   Observed = unlist(lapply(1:length(pt_COs_vs_features), function(x) {
                       pt_COs_vs_features[[x]]$numOverlaps$observed
                   })),
                   Expected = unlist(lapply(1:length(pt_COs_vs_features), function(x) {
                       mean(pt_COs_vs_features[[x]]$numOverlaps$permuted, na.rm=T)
                   })),
                   AlphaThreshold = unlist(lapply(1:length(pt_COs_vs_features), function(x) {
                       if(pt_COs_vs_features[[x]]$numOverlaps$alternative == "less") {
                           quantile(pt_COs_vs_features[[x]]$numOverlaps$permuted, probs=0.05, na.rm=T)[[1]]
                       } else if(pt_COs_vs_features[[x]]$numOverlaps$alternative == "greater") {
                           quantile(pt_COs_vs_features[[x]]$numOverlaps$permuted, probs=0.95, na.rm=T)[[1]]
                       } else {
                           stop(paste0(pt_COs_vs_features[[x]]$numOverlaps$alternative, " is not 'less' or 'greater'"))
                       }
                   })),
                   Pvalue = unlist(lapply(1:length(pt_COs_vs_features), function(x) {
                       pt_COs_vs_features[[x]]$numOverlaps$pval
                   })),
                   Zscore = unlist(lapply(1:length(pt_COs_vs_features), function(x) {
                       pt_COs_vs_features[[x]]$numOverlaps$zscore
                   })),
                   AlternativeHypothesis = unlist(lapply(1:length(pt_COs_vs_features), function(x) {
                       pt_COs_vs_features[[x]]$numOverlaps$alternative
                   })),
                   Permutations = unlist(lapply(1:length(pt_COs_vs_features), function(x) {
                       pt_COs_vs_features[[x]]$numOverlaps$ntimes
                   }))
                  )
                   


write.table(pt_dist_DF,
            file=paste0(outdir, region, "_GBS_COs",
                        "_alnTo_", alnTo, "_",
                        paste0(chrs, collapse="_"), "_genes_perm_test_distribution.tsv"),
            quote=F, sep="\t", row.names=F, col.names=T)
write.table(pt_DF,
            file = paste0(outdir, region, "_GBS_COs",
                          "_alnTo_", alnTo, "_",
                          paste0(chrs, collapse="_"), "_genes_perm_test_summary.tsv"),
            quote=F, sep="\t", row.names=F, col.names=T)

pt_dist_DF = read.table(paste0(outdir, region, "_GBS_COs",
                               "_alnTo_", alnTo, "_",
                               paste0(chrs, collapse="_"), "_genes_perm_test_distribution.tsv"),
                        sep="\t", header=T)
pt_DF = read.table(paste0(outdir, region, "_GBS_COs",
                          "_alnTo_", alnTo, "_",
                          paste0(chrs, collapse="_"), "_genes_perm_test_summary.tsv"),
                   sep="\t", header=T)

            
pt_DF$Feature = factor(pt_DF$Feature,
                       levels = rev(unique(pt_DF$Feature)))
pt_dist_DF$Feature = factor(pt_dist_DF$Feature,
                            levels = rev(unique(pt_dist_DF$Feature)))

pt_DF$Region = factor(pt_DF$Region,
                      levels = rev(unique(pt_DF$Region)))
pt_dist_DF$Region = factor(pt_dist_DF$Region,
                           levels = rev(unique(pt_dist_DF$Region)))


# Define function to make colours transparent,
# to aid visibility where points overlap
makeTransparent = function(thisColour, alpha = 230)
{
  newColour = col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}


vp_all = ggplot(data = pt_dist_DF,
                mapping = aes(x = Feature,
                              y = Permuted)) +
  xlab("Feature category") +
  ylab(bquote(.(region_plot) ~ "COs vs random overlaps with features")) +
  geom_violin(trim = F,
              scale = "count",
              colour = "grey70",
              fill = "grey70") +
  geom_point(data = pt_DF,
             mapping = aes(x = Feature,
                           y = AlphaThreshold),
             position = position_nudge(x = -0.04),
             shape = 124, colour = makeTransparent("darkorange1"), size = 12) +
  geom_point(data = pt_DF,
             mapping = aes(x = Feature,
                           y = Observed),
             position = position_nudge(x = -0.04),
             shape = 124, colour = makeTransparent("dodgerblue2"), size = 12) +
  geom_point(data = pt_DF,
             mapping = aes(x = Feature,
                           y = Expected),
             position = position_nudge(x = -0.04),
             shape = 124, colour = makeTransparent("black"), size = 12) +

  coord_flip() +
  theme_bw() +
  theme(
#        axis.line.y = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_line(linewidth = 0.5, colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black", hjust = 0.5, vjust = 0.5, angle = 0),
        axis.title.y = element_text(size = 20, colour = "black"),
#        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.x = element_line(linewidth = 0.5, colour = "black"),
        axis.text.x = element_text(size = 18, colour = "black", hjust = 0.5, vjust = 0.5, angle = 0),
        axis.title.x = element_text(size = 20, colour = "black", margin = margin(t = 20)),
        strip.text.x = element_text(size = 20, colour = "black"),
        strip.text.y = element_text(size = 20, colour = "black"),
        strip.background = element_rect(fill = "white", colour = "white"),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        plot.title = element_text(size = 16, colour = "black", hjust = 0.5)) +
  ggtitle(bquote("Permutations =" ~
                 .(prettyNum(perms,
                             big.mark = ",",
                             trim = T)) ~
                 "sets of randomly positioned loci"))

ggsave(paste0(plotdir, region, "_GBS_COs",
              "_alnTo_", alnTo, "_",
              paste0(chrs, collapse="_"), "_genes_perm_test_violin.pdf"),
       plot = vp_all,
       width = 12, height = 8, limitsize = F)
