#!/usr/bin/env Rscript

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 17/11/2022

# Do permutation tests to evaluate overlap between candidate recombinant
# read segment pair alignment intervals and features of interest

# Usage:
# conda activate python_3.9.6
# ./recomb_overlap_permTest.R ColLerF1pollen_1000bp_minq90 Col-0.ragtag_scaffolds Ler-0_110x.ragtag_scaffolds 24 0.9 11 Col-0.ragtag_scaffolds_Chr 0.9 not_centromere 'Chr1,Chr2,Chr3,Chr4,Chr5' 10000 COs_and_NCOs
# conda deactivate

args = commandArgs(trailingOnly=T)
readsPrefix = args[1]
acc1 = args[2]
acc2 = args[3]
kmerSize = as.integer(args[4])
overlapProp = as.numeric(args[5])
minHits = as.integer(args[6])
alnTo = args[7]
alenTOqlen = args[8]
region = args[9]
chrom = unlist(strsplit(args[10], split=","))
perms = as.integer(args[11])
recombType = args[12] 

#readsPrefix = "ColLerF1pollen_1000bp_minq90"
#acc1 = "Col-0.ragtag_scaffolds"
#acc2 = "Ler-0_110x.ragtag_scaffolds"
#kmerSize = as.integer(24)
#overlapProp = 0.9
#minHits = as.integer(11)
#alnTo = "Col-0.ragtag_scaffolds_Chr"
#alenTOqlen = 0.9
#region = "not_centromere"
#chrom = unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                        split=","))
#perms = as.integer(10000)
#recombType = "COs_and_NCOs"

# Set minimum possible P-value for permutation test result with
# perms sets of random loci
min_pval = 1 - ( (perms - 1) / perms)
region_plot = gsub("not_", "non-", region)
region_plot = gsub("centromere", "centromeric", region_plot)

library(regioneR)
library(rtracklayer)
library(dplyr)
library(ggplot2)


indir = paste0(region, "/segment_pairs/")
outdir = paste0(indir, "perm_tests/")
plotdir = paste0(outdir, "plots/")
system(paste0("[ -d ", plotdir, " ] || mkdir -p ", plotdir))


# Accession names
acc1_name = strsplit( strsplit(acc1, split="\\.")[[1]][1],
                      split="_" )[[1]][1]
acc2_name = strsplit( strsplit(acc2, split="\\.")[[1]][1],
                      split="_" )[[1]][1]
alnTo_name =  strsplit( strsplit(alnTo, split="\\.")[[1]][1],
                        split="_" )[[1]][1] 


# Filtered read segment alignment pairs TSV (candidate recombinant reads)
# COs
COs_DF_list = lapply(1:length(chrom), function(x) {
    COs_tsv = paste0(indir, "co/", readsPrefix,
                     "_", acc1, "_", acc2,
                     "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                     "_hom_maxDist_aTOq", alenTOqlen, "_co",
                     "_alnTo_", alnTo, "_",
                     chrom[x], ".tsv")
    # File exists sanity check
    COs_DF_x = read.table(COs_tsv, header=T)
    COs_DF_x_filt = COs_DF_x[which(COs_DF_x$acc1_tname == chrom[x]),]
    # Destack
    COs_DF_x_filt_destack = data.frame()
    for(y in 1:nrow(COs_DF_x_filt)) {
        COs_DF_x_filt_y = COs_DF_x_filt[y,]
        COs_DF_x_filt_y_start_dups = COs_DF_x_filt[ which( (COs_DF_x_filt$event_start == COs_DF_x_filt_y$event_start &
                                                            COs_DF_x_filt$acc1_tname == COs_DF_x_filt_y$acc1_tname) |
                                                           (COs_DF_x_filt$event_end == COs_DF_x_filt_y$event_end &
                                                            COs_DF_x_filt$acc1_tname == COs_DF_x_filt_y$acc1_tname) ), ]
        if(nrow(COs_DF_x_filt_y_start_dups) == 1) {
            COs_DF_x_filt_destack = rbind(COs_DF_x_filt_destack, COs_DF_x_filt_y)
        }
    }
    COs_DF_x_filt_destack
})
if(length(COs_DF_list) > 1) {
    COs_DF = dplyr::bind_rows(COs_DF_list)
} else {
    COs_DF = COs_DF_list[[1]]
}

write.table(COs_DF,
            file=paste0(indir, "co/", readsPrefix,
                        "_", acc1, "_", acc2,
                        "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                        "_hom_maxDist_aTOq", alenTOqlen, "_co",
                        "_alnTo_", alnTo, "_",
                        "destacked.tsv"),
            quote=F, sep="\t", row.names=F, col.names=T)

COs_GR = GRanges(seqnames=COs_DF$acc1_tname,
                 ranges=IRanges(start=COs_DF$event_start,
                                end=COs_DF$event_end),
                 strand="*")


# NCOs
NCOs_DF_list = lapply(1:length(chrom), function(x) {
    NCOs_tsv = paste0(indir, "nco/", readsPrefix,
                      "_", acc1, "_", acc2,
                      "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                      "_hom_maxDist_aTOq", alenTOqlen, "_nco",
                      "_alnTo_", alnTo, "_",
                      chrom[x], ".tsv")
    # File exists sanity check
    NCOs_DF_x = read.table(NCOs_tsv, header=T)
    NCOs_DF_x_filt = NCOs_DF_x[which(NCOs_DF_x$acc1_tname == chrom[x]),]
    # Destack
    NCOs_DF_x_filt_destack = data.frame()
    for(y in 1:nrow(NCOs_DF_x_filt)) {
        NCOs_DF_x_filt_y = NCOs_DF_x_filt[y,]
        NCOs_DF_x_filt_y_start_dups = NCOs_DF_x_filt[ which( (NCOs_DF_x_filt$event_start == NCOs_DF_x_filt_y$event_start &
                                                              NCOs_DF_x_filt$acc1_tname == NCOs_DF_x_filt_y$acc1_tname) |
                                                             (NCOs_DF_x_filt$event_end == NCOs_DF_x_filt_y$event_end &
                                                              NCOs_DF_x_filt$acc1_tname == NCOs_DF_x_filt_y$acc1_tname) ), ]
        if(nrow(NCOs_DF_x_filt_y_start_dups) == 1) {
            NCOs_DF_x_filt_destack = rbind(NCOs_DF_x_filt_destack, NCOs_DF_x_filt_y)
        }
    }
    NCOs_DF_x_filt_destack
})
if(length(NCOs_DF_list) > 1) {
    NCOs_DF = dplyr::bind_rows(NCOs_DF_list)
} else {
    NCOs_DF = NCOs_DF_list[[1]]
}

write.table(NCOs_DF,
            file=paste0(indir, "nco/", readsPrefix,
                        "_", acc1, "_", acc2,
                        "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                        "_hom_maxDist_aTOq", alenTOqlen, "_nco",
                        "_alnTo_", alnTo, "_",
                        "destacked.tsv"),
            quote=F, sep="\t", row.names=F, col.names=T)

NCOs_GR = GRanges(seqnames=NCOs_DF$acc1_tname,
                  ranges=IRanges(start=NCOs_DF$event_start,
                                 end=NCOs_DF$event_end),
                  strand="*")


if(recombType == "COs_and_NCOs") {
   recomb_GR = c(COs_GR, NCOs_GR) 
} else if(recombType == "COs") {
   recomb_GR = COs_GR
} else if(recombType == "NCOs") {
   recomb_GR = NCOs_GR
} else {
   stop("recombType is not 'COs_and_NCOs', 'COs' or 'NCOs'")
}


# Genomic definitions

# CEN coordinates
CEN = read.csv(paste0("/rds/project/rds-O5Ty9yVfQKg/pancentromere/centromeric_coordinates/",
                      "centromere_manual_EDTA4_fa.csv"),
               header=T)
CEN$fasta.name = gsub(".fa", "", CEN$fasta.name)

# alnTo accession
alnTo_fai = read.table(paste0("index/", alnTo, ".fa.fai"),
                       header=F)
alnTo_fai =  alnTo_fai[which(alnTo_fai$V1 %in% chrom),]
alnTo_chrs = alnTo_fai$V1
alnTo_chrLens = alnTo_fai$V2

alnTo_CEN = CEN[which(CEN$fasta.name == gsub("_Chr", "", alnTo)),]
alnTo_CEN = alnTo_CEN[,which(colnames(alnTo_CEN) %in% c("chr", "start", "end"))]
alnTo_CEN_new = data.frame()
for(i in 1:length(alnTo_chrs)) {
  alnTo_CEN_chr = alnTo_CEN[which(alnTo_CEN$chr == alnTo_chrs[i]),]
  if(nrow(alnTo_CEN_chr) > 1) {
    alnTo_CEN_chr = data.frame(chr=alnTo_CEN_chr$chr[1],
                               start=alnTo_CEN_chr$start[1],
                               end=alnTo_CEN_chr$end[nrow(alnTo_CEN_chr)])
  }
  alnTo_CEN_new = rbind(alnTo_CEN_new, alnTo_CEN_chr)
}
alnTo_CEN = alnTo_CEN_new

alnTo_CEN_GR = GRanges(seqnames=alnTo_CEN$chr,
                       ranges=IRanges(start=alnTo_CEN$start,
                                      end=alnTo_CEN$end),
                       strand="*")

alnTo_nonCEN = data.frame(chr=rep(alnTo_chrs, 2),
                          start=c( rep(1, length(alnTo_chrs)), (alnTo_CEN$end + 1) ),
                          end=c( (alnTo_CEN$start - 1), alnTo_chrLens))

alnTo_nonCEN_GR = GRanges(seqnames=alnTo_nonCEN$chr,
                          ranges=IRanges(start=alnTo_nonCEN$start,
                                         end=alnTo_nonCEN$end),
                          strand="*")

# Define region to be analysed
if(region == "not_centromere") {
    region_GR = alnTo_nonCEN_GR
    mask_GR = alnTo_CEN_GR
} else if(region == "centromere") {
    region_GR = alnTo_CEN_GR
    mask_GR = alnTo_nonCEN_GR
} else {
    stop("region is not 'not_centromere' or 'centromere'")
}

genome_GR = GRanges(seqnames=alnTo_chrs,
                    ranges=IRanges(start=rep(1, length(alnTo_chrs)),
                                   end=alnTo_chrLens),
                    strand="*")


# Subset to include only those not overlapping masked region (e.g., centromere)                                                          
mask_recomb_overlap = findOverlaps(query=mask_GR,
                                   subject=recomb_GR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
if(length(mask_recomb_overlap) > 0) {
    recomb_GR = recomb_GR[-subjectHits(mask_recomb_overlap)]
}
print("Candidate recombination events:")
print(recomb_GR)


# Load genes
genes = readGFF(paste0("/rds/project/rds-O5Ty9yVfQKg/pancentromere/annotation/genes/",
                        alnTo, "/", alnTo, ".genes.gff3"))
genes = genes[which(genes$seqid %in% alnTo_chrs),]
genes = genes[which(genes$type == "mRNA"),]
print(dim(genes))

genes_GR = GRanges(seqnames=genes$seqid,
                   ranges=IRanges(start=genes$start,
                                  end=genes$end),
                   strand=genes$strand)

# Subset to include only those not overlapping masked region
mask_genes_overlap <- findOverlaps(query = mask_GR,
                                   subject = genes_GR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
genes_GR <- genes_GR[-subjectHits(mask_genes_overlap)]

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
pt_recomb_vs_features = lapply(1:length(features_GR_list), function(x) {
    permTest(A=recomb_GR,
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
                        Feature = unlist(lapply(1:length(pt_recomb_vs_features), function(x) {
                            rep(names(features_GR_list)[x],
                                times=length(pt_recomb_vs_features[[x]]$numOverlaps$permuted))
                        })),
                        Permuted = unlist(lapply(1:length(pt_recomb_vs_features), function(x) {
                            pt_recomb_vs_features[[x]]$numOverlaps$permuted
                        }))
                      )
                        

pt_DF = data.frame(
                   Feature = names(features_GR_list),
                   Number_of_features = unlist(lapply(1:length(features_GR_list), function(x) {
                       length(features_GR_list[[x]])
                   })),
                   Observed = unlist(lapply(1:length(pt_recomb_vs_features), function(x) {
                       pt_recomb_vs_features[[x]]$numOverlaps$observed
                   })),
                   Expected = unlist(lapply(1:length(pt_recomb_vs_features), function(x) {
                       mean(pt_recomb_vs_features[[x]]$numOverlaps$permuted, na.rm=T)
                   })),
                   AlphaThreshold = unlist(lapply(1:length(pt_recomb_vs_features), function(x) {
                       if(pt_recomb_vs_features[[x]]$numOverlaps$alternative == "less") {
                           quantile(pt_recomb_vs_features[[x]]$numOverlaps$permuted, probs=0.05, na.rm=T)[[1]]
                       } else if(pt_recomb_vs_features[[x]]$numOverlaps$alternative == "greater") {
                           quantile(pt_recomb_vs_features[[x]]$numOverlaps$permuted, probs=0.95, na.rm=T)[[1]]
                       } else {
                           stop(paste0(pt_recomb_vs_features[[x]]$numOverlaps$alternative, " is not 'less' or 'greater'"))
                       }
                   })),
                   Pvalue = unlist(lapply(1:length(pt_recomb_vs_features), function(x) {
                       pt_recomb_vs_features[[x]]$numOverlaps$pval
                   })),
                   Zscore = unlist(lapply(1:length(pt_recomb_vs_features), function(x) {
                       pt_recomb_vs_features[[x]]$numOverlaps$zscore
                   })),
                   AlternativeHypothesis = unlist(lapply(1:length(pt_recomb_vs_features), function(x) {
                       pt_recomb_vs_features[[x]]$numOverlaps$alternative
                   })),
                   Permutations = unlist(lapply(1:length(pt_recomb_vs_features), function(x) {
                       pt_recomb_vs_features[[x]]$numOverlaps$ntimes
                   }))
                  )
                   


write.table(pt_dist_DF,
            file=paste0(outdir, readsPrefix,
                          "_", acc1, "_", acc2,
                          "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                          "_hom_maxDist_aTOq", alenTOqlen, "_", recombType,
                          "_alnTo_", alnTo, "_",
                          paste0(alnTo_chrs, collapse="_"), "_genes_perm_test_distribution.tsv"),
            quote=F, sep="\t", row.names=F, col.names=T)
write.table(pt_DF,
            file = paste0(outdir, readsPrefix,
                          "_", acc1, "_", acc2,
                          "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                          "_hom_maxDist_aTOq", alenTOqlen, "_", recombType,
                          "_alnTo_", alnTo, "_",
                          paste0(alnTo_chrs, collapse="_"), "_genes_perm_test_summary.tsv"),
            quote=F, sep="\t", row.names=F, col.names=T)

pt_dist_DF = read.table(paste0(outdir, readsPrefix,
                               "_", acc1, "_", acc2,
                               "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                               "_hom_maxDist_aTOq", alenTOqlen, "_", recombType,
                               "_alnTo_", alnTo, "_",
                               paste0(alnTo_chrs, collapse="_"), "_genes_perm_test_distribution.tsv"),
                        sep="\t", header=T)
pt_DF = read.table(paste0(outdir, readsPrefix,
                          "_", acc1, "_", acc2,
                          "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                          "_hom_maxDist_aTOq", alenTOqlen, "_", recombType,
                          "_alnTo_", alnTo, "_",
                          paste0(alnTo_chrs, collapse="_"), "_genes_perm_test_summary.tsv"),
                   sep="\t", header=T)

            
pt_DF$Feature <- factor(pt_DF$Feature,
                        levels = rev(unique(pt_DF$Feature)))
pt_dist_DF$Feature <- factor(pt_dist_DF$Feature,
                             levels = rev(unique(pt_dist_DF$Feature)))


# Define function to make colours transparent,
# to aid visibility where points overlap
makeTransparent <- function(thisColour, alpha = 230)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}


vp_all <- ggplot(data = pt_dist_DF,
                 mapping = aes(x = Feature,
                               y = Permuted)) +
  xlab("Feature category") +
  ylab(bquote(.(gsub("_", " ", recombType)) ~ "vs random overlaps with features")) +
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
                 "sets of randomly positioned non-centromeric loci"))

ggsave(paste0(plotdir, readsPrefix,
              "_", acc1, "_", acc2,
              "_k", kmerSize, "_op", overlapProp, "_h", minHits,
              "_hom_maxDist_aTOq", alenTOqlen, "_", recombType,
              "_alnTo_", alnTo, "_",
              paste0(alnTo_chrs, collapse="_"), "_genes_perm_test_violin.pdf"),
       plot = vp_all,
       width = 12, height = 8, limitsize = F)
