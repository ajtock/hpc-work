#!/usr/bin/env Rscript

# Use DMRcaller v1.24.0 to identify hypomethylated DMRs between two conditions
# (e.g. mutant and wild type)

# Usage:
# source ~/.bashrc
# conda activate python_3.9.6
# ./DMRcaller_hypo_nomerge.R --condition1 'WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013,WT_BSseq_Rep1_2014' \
#                            --condition2 'cmt2_BSseq_Rep1' \
#                            --refbase t2t-col.20210610 \
#                            --chrName 'Chr1,Chr2,Chr3,Chr4,Chr5' \
#                            --genomeRegion genomewide \
#                            --quantiles 6 \
#                            --context CHG
# conda deactivate

library(argparse)
library(DMRcaller)
library(betareg) # required by DMRcaller::computeDMRsReplicates 
library(rtracklayer)

# Create parser object
parser <- ArgumentParser()

# Specify arguments
# ArgumentParser will add a help option by default
parser$add_argument("--condition1", type = "character",
                    help="Sample replicate prefixes for first condition for inclusion in DMRcaller contrast.")
parser$add_argument("--condition2", type = "character",
                    help="Sample replicate prefixes for second condition for inclusion in DMRcaller contrast.")
parser$add_argument("--refbase", type = "character", default = "t2t-col.20210610",
                    help="Reference genome base name. Default: %(default)s")
parser$add_argument("--chrName", type = "character",
                    help="Reference genome chromosome names for inclusion in DMRcaller contrast.")
parser$add_argument("--genomeRegion", type = "character", default = "genomewide",
                    help="Reference genome region from which to extract and export DMRs. Default: %(default)s")
parser$add_argument("--quantiles", type = "double", default = 6,
                    help="Number of groups between which to divide DMRs based on percentile rank. Default: %(default)s")
parser$add_argument("--context", type = "character",
                    help="cytosine methylation context.")

# parse arguments
args <- parser$parse_args()
print(args)

#args_file <- "tempArgsObjectFile.rds"
#saveRDS(args, args_file); print(args)
##system("./DMRcaller_hypo_nomerge.R --condition1 'WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013,WT_BSseq_Rep1_2014' --condition2 'cmt2_BSseq_Rep1' --refbase 't2t-col.20210610' --chrName 'Chr1,Chr2,Chr3,Chr4,Chr5' --genomeRegion genomewide --quantiles 6 --context CHG")
#args <- readRDS(args_file)

args$condition1 <- unlist(strsplit(args$condition1, split = ","))
args$condition2 <- unlist(strsplit(args$condition2, split = ","))
args$chrName <- unlist(strsplit(args$chrName, split = ","))

print(args)

# Set context-specific minProportionDifference to pass to computeDMRs or computeDMRsReplicates
if(args$context == "CpG") {
  minProportionDifference_context <- 0.4
} else if(args$context == "CHG") {
  minProportionDifference_context <- 0.2
} else if(args$context == "CHH") {
  minProportionDifference_context <- 0.1
}
print(paste0(args$context, " minProportionDifference = ", minProportionDifference_context))

# Load methylation data
condition1_Reps <- mclapply(seq_along(args$condition1), function(x) {
  readBismark(paste0(args$condition1[x], "_MappedOn_", args$refbase, "_dedup_", args$context, ".CX_report.txt.gz"))
}, mc.cores = length(args$condition1))
condition2_Reps <- mclapply(seq_along(args$condition2), function(x) {
  readBismark(paste0(args$condition2[x], "_MappedOn_", args$refbase, "_dedup_", args$context, ".CX_report.txt.gz"))
}, mc.cores = length(args$condition2))

# Calculate conversion rate and adjust coverage values accordingly
# See https://doi.org/10.1007/978-1-0716-1134-0_21
# (full link: https://link.springer.com/protocol/10.1007%2F978-1-0716-1134-0_21 )
for(x in 1:length(condition1_Reps)) {
  print(paste0(x, ">>"))
  # Get ranges corresponding to the given context
  condition1_Reps[[x]] <- condition1_Reps[[x]][condition1_Reps[[x]]$context == sub("p", "", args$context)]

  # Extract chloroplast methylation data
  condition1_Repx_ChrC <- condition1_Reps[[x]][condition1_Reps[[x]]$context == sub("p", "", args$context) &
                                               seqnames(condition1_Reps[[x]]) == "ChrC"]

  # Calculate conversion
  condition1_Repx_conv <- 1 - ( sum(mcols(condition1_Repx_ChrC)$readsM) / sum(mcols(condition1_Repx_ChrC)$readsN) )

  # Adjust methylated and total read counts (readsM and readsN)
  condition1_Reps[[x]]$readsM <- round( condition1_Reps[[x]]$readsM - condition1_Reps[[x]]$readsN * (1 - condition1_Repx_conv) )
  condition1_Reps[[x]]$readsM[condition1_Reps[[x]]$readsM < 0] <- 0 
  condition1_Reps[[x]]$readsN <- round( condition1_Reps[[x]]$readsN * condition1_Repx_conv )

  # Remove superfluous ranges
  # Get ranges corresponding to those in chrName
  condition1_Reps[[x]] <- keepSeqlevels(condition1_Reps[[x]], args$chrName, pruning.mode = "coarse")

  # Sort by seqnames, start and end
  condition1_Reps[[x]] <- sortSeqlevels(condition1_Reps[[x]])
  condition1_Reps[[x]] <- sort(condition1_Reps[[x]], ignore.strand = TRUE)

  print(condition1_Reps[[x]])
  print(paste0("<<", x))
}

for(x in 1:length(condition2_Reps)) {
  print(paste0(x, ">>"))
  # Get ranges corresponding to the given context
  condition2_Reps[[x]] <- condition2_Reps[[x]][condition2_Reps[[x]]$context == sub("p", "", args$context)]

  # Extract chloroplast methylation data
  condition2_Repx_ChrC <- condition2_Reps[[x]][condition2_Reps[[x]]$context == sub("p", "", args$context) &
                                               seqnames(condition2_Reps[[x]]) == "ChrC"]

  # Calculate conversion
  condition2_Repx_conv <- 1 - ( sum(mcols(condition2_Repx_ChrC)$readsM) / sum(mcols(condition2_Repx_ChrC)$readsN) )

  # Adjust methylated and total read counts (readsM and readsN)
  condition2_Reps[[x]]$readsM <- round( condition2_Reps[[x]]$readsM - condition2_Reps[[x]]$readsN * (1 - condition2_Repx_conv) )
  condition2_Reps[[x]]$readsM[condition2_Reps[[x]]$readsM < 0] <- 0 
  condition2_Reps[[x]]$readsN <- round( condition2_Reps[[x]]$readsN * condition2_Repx_conv )

  # Remove superfluous ranges
  # Get ranges corresponding to those in chrName
  condition2_Reps[[x]] <- keepSeqlevels(condition2_Reps[[x]], args$chrName, pruning.mode = "coarse")

  # Sort by seqnames, start and end
  condition2_Reps[[x]] <- sortSeqlevels(condition2_Reps[[x]])
  condition2_Reps[[x]] <- sort(condition2_Reps[[x]], ignore.strand = TRUE)

  print(condition2_Reps[[x]])
  print(paste0("<<", x))
}

# Define output directories
plotDir <- "mC_coverage_plots_DMRcaller/"
hypoDMRdir <- paste0("DMRs/hypoDMRs/", paste0(args$chrName, collapse = "_"), "/")
#hyperDMRdir <- paste0("DMRs/hyperDMRs/", paste0(args$chrName, collapse = "_"), "/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))
system(paste0("[ -d ", hypoDMRdir, " ] || mkdir -p ", hypoDMRdir))
#system(paste0("[ -d ", hyperDMRdir, " ] || mkdir -p ", hyperDMRdir))

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/rds/hpc-work/nanopore/", args$refbase, "/", args$refbase, ".fa.fai"), header = F)
chrs <- fai$V1[which(fai$V1 %in% args$chrName)]
chrLens <- fai$V2[which(fai$V1 %in% args$chrName)]

CEN <- read.table(paste0("/home/ajt200/rds/hpc-work/nanopore/", args$refbase, "/", args$refbase, ".fa.centromeres"), header = T)
CEN <- CEN[which(fai$V1 %in% args$chrName),]
# Not sure Dan defined pericentromeres at 10-kb resolution (appear to be at 100-kb resolution)
periCEN <- read.table(paste0("/home/ajt200/rds/hpc-work/nanopore/", args$refbase, "/", args$refbase, ".fa.Dan_pericentromeres"), header = T)
periCEN <- periCEN[which(fai$V1 %in% args$chrName),]

CENGR <- GRanges(seqnames = chrs,
                 ranges = IRanges(start = CEN$start,
                                  end = CEN$end),
                 strand = "*")
CENGR <- CENGR[which(seqnames(CENGR)@values %in% args$chrName)]

nonCENGR <- GRanges(seqnames = rep(chrs, 2),
                    ranges = IRanges(start = c(rep(1, length(chrs)),
                                               CEN$end+1),
                                     end = c(CEN$start-1,
                                             chrLens)),
                    strand = "*")
nonCENGR <- nonCENGR[which(seqnames(nonCENGR)@values %in% args$chrName)]

periCENGR <- GRanges(seqnames = chrs,
                     ranges = IRanges(start = periCEN$start,
                                      end = periCEN$end),
                     strand = "*")
periCENGR <- periCENGR[which(seqnames(periCENGR)@values %in% args$chrName)]

armGR <- GRanges(seqnames = rep(chrs, 2),
                 ranges = IRanges(start = c(rep(1, length(chrs)),
                                            periCEN$end+1),
                                  end = c(periCEN$start-1,
                                          chrLens)),
                 strand = "*")
armGR <- armGR[which(seqnames(armGR)@values %in% args$chrName)]

# Define genomic regions to be analysed (genomeRegionGR)
if(args$genomeRegion == "arm") {
  genomeRegionGR <- GRanges(seqnames = rep(chrs, 2),
                            ranges = IRanges(start = c(rep(1, length(chrs)),
                                                       periCEN$end+1),
                                             end = c(periCEN$start-1,
                                                     chrLens)),
                            strand = "*")
  genomeRegionGR <- genomeRegionGR[which(seqnames(genomeRegionGR)@values %in% args$chrName)]
} else if(args$genomeRegion == "peri") {
  genomeRegionGR <- GRanges(seqnames = chrs,
                            ranges = IRanges(start = periCEN$start,
                                             end = periCEN$end),
                            strand = "*")
  genomeRegionGR <- genomeRegionGR[which(seqnames(genomeRegionGR)@values %in% args$chrName)]
} else if(args$genomeRegion == "CEN") {
  genomeRegionGR <- GRanges(seqnames = chrs,
                            ranges = IRanges(start = CEN$start,
                                             end = CEN$end),
                            strand = "*")
  genomeRegionGR <- genomeRegionGR[which(seqnames(genomeRegionGR)@values %in% args$chrName)]
} else if(args$genomeRegion == "nonCEN") {
  genomeRegionGR <- GRanges(seqnames = rep(chrs, 2),
                            ranges = IRanges(start = c(rep(1, length(chrs)),
                                                       CEN$end+1),
                                             end = c(CEN$start-1,
                                                     chrLens)),
                            strand = "*")
  genomeRegionGR <- genomeRegionGR[which(seqnames(genomeRegionGR)@values %in% args$chrName)]
} else if(args$genomeRegion == "genomewide") {
  genomeRegionGR <- GRanges(seqnames = chrs,
                            ranges = IRanges(start = rep(1, length(chrs)),
                                             end = chrLens),
                            strand = "*")
  genomeRegionGR <- genomeRegionGR[which(seqnames(genomeRegionGR)@values %in% args$chrName)]
} else {
  stop("genomeRegion is not arm, peri, CEN, nonCEN or genomewide")
}

# Define genomic regions to be masked from analyses (genomeMaskGR)
if(args$genomeRegion == "arm") {
  genomeMaskGR <- GRanges(seqnames = chrs,
                          ranges = IRanges(start = periCEN$start,
                                           end = periCEN$end),
                          strand = "*")
  genomeMaskGR <- genomeMaskGR[which(seqnames(genomeMaskGR)@values %in% args$chrName)]
} else if(args$genomeRegion == "peri") {
  genomeMaskGR <- GRanges(seqnames = rep(chrs, 2),
                          ranges = IRanges(start = c(rep(1, length(chrs)),
                                                     periCEN$end+1),
                                           end = c(periCEN$start-1,
                                                   chrLens)),
                          strand = "*")
  genomeMaskGR <- genomeMaskGR[which(seqnames(genomeMaskGR)@values %in% args$chrName)]
} else if(args$genomeRegion == "CEN") {
  genomeMaskGR <- GRanges(seqnames = rep(chrs, 2),
                          ranges = IRanges(start = c(rep(1, length(chrs)),
                                                     CEN$end+1),
                                           end = c(CEN$start-1,
                                                   chrLens)),
                          strand = "*")
  genomeMaskGR <- genomeMaskGR[which(seqnames(genomeMaskGR)@values %in% args$chrName)]
} else if(args$genomeRegion == "nonCEN") {
  genomeMaskGR <- GRanges(seqnames = chrs,
                          ranges = IRanges(start = CEN$start,
                                           end = CEN$end),
                          strand = "*")
  genomeMaskGR <- genomeMaskGR[which(seqnames(genomeMaskGR)@values %in% args$chrName)]
} else if(args$genomeRegion == "genomewide") {
  genomeMaskGR <- GRanges()
  genomeMaskGR <- genomeMaskGR[which(seqnames(genomeMaskGR)@values %in% args$chrName)]
} else {
  stop("genomeRegion is not arm, peri, CEN, nonCEN or genomewide")
}


if(length(condition2_Reps) == 1) {

  # Plot methylation coverage for each replicate and condition
  pdf(paste0(plotDir, "plotMethylationDataCoverage_",
             paste0(args$condition1, collapse = "_"), "_",
             paste0(args$condition2, collapse = "_"), "_",
             paste0(args$chrName, collapse = "_"), "_",
             paste0(args$context, collapse = "_"),
             ".pdf"))
  plotMethylationDataCoverage(methylationData1 = condition1_Reps[[1]],
                              methylationData2 = condition1_Reps[[2]],
                              breaks = c(1, 3, 6, 9, 12, 15),
                              regions = NULL,
                              conditionsNames = c(args$condition1[1], args$condition1[2]),
                              context = sub("p", "", args$context),
                              proportion = TRUE,
                              labels=LETTERS,
                              contextPerRow = FALSE)
  plotMethylationDataCoverage(methylationData1 = condition1_Reps[[3]],
                              methylationData2 = condition2_Reps[[1]],
                              breaks = c(1, 3, 6, 9, 12, 15),
                              regions = NULL,
                              conditionsNames = c(args$condition1[3], args$condition2[1]),
                              context = sub("p", "", args$context),
                              proportion = TRUE,
                              labels=LETTERS,
                              contextPerRow = FALSE)
  dev.off()

  set.seed(83949341)
  # Compute DMRs using "bins" method
  DMRs_per_Rep_list_bins <- lapply(seq_along(condition1_Reps), function(x) {
    computeDMRs(methylationData1 = condition1_Reps[[x]],
                methylationData2 = condition2_Reps[[1]],
                regions = NULL,
                context = sub("p", "", args$context),
                method = "bins",
                binSize = 100,
                test = "fisher",
                pValueThreshold = 0.01,
                minCytosinesCount = 4,
                minProportionDifference = minProportionDifference_context,
                minGap = 0,
                minSize = 50,
                minReadsPerCytosine = 6,
                cores = detectCores())
  })

  hypoDMRs_per_Rep_list_bins <- lapply(seq_along(DMRs_per_Rep_list_bins), function(x) {
    unique(DMRs_per_Rep_list_bins[[x]][which(DMRs_per_Rep_list_bins[[x]]$regionType == "loss"),])
  })

  #hyperDMRs_per_Rep_list_bins <- lapply(seq_along(DMRs_per_Rep_list_bins), function(x) {
  #  unique(DMRs_per_Rep_list_bins[[x]][which(DMRs_per_Rep_list_bins[[x]]$regionType == "gain"),])
  #})

  hypoDMRs_allReps_bins <- hypoDMRs_per_Rep_list_bins[[1]]
  for(x in 2:length(hypoDMRs_per_Rep_list_bins)) {
    hits <- findOverlaps(query = hypoDMRs_allReps_bins,
                         subject = hypoDMRs_per_Rep_list_bins[[x]],
                         type = "equal", select = "all",
                         ignore.strand = FALSE)
    hypoDMRs_allReps_bins <- hypoDMRs_allReps_bins[unique(queryHits(hits))]
  }
  
#  set.seed(83949341)
#  hypoDMRs_allReps_bins <- mergeDMRsIteratively(DMRs = hypoDMRs_allReps_bins,
#                                                minGap = 200,
#                                                respectSigns = TRUE,
#                                                methylationData1 = condition1_Reps[[1]],
#                                                methylationData2 = condition2_Reps[[1]],
#                                                context = sub("p", "", args$context),
#                                                minProportionDifference = minProportionDifference_context,
#                                                minReadsPerCytosine = 6,
#                                                pValueThreshold = 0.01,
#                                                test = "fisher",
#                                                alternative = "two.sided",
#                                                cores = detectCores()) 
#  hypoDMRs_allReps_bins <- unique(hypoDMRs_allReps_bins)

  # Get args$genomeRegion DMRs
  hypoDMRs_allReps_bins_genomeRegion_hits <- findOverlaps(query = genomeRegionGR,
                                                          subject = hypoDMRs_allReps_bins,
                                                          type = "any", select = "all",
                                                          ignore.strand = TRUE)
  hypoDMRs_allReps_bins <- hypoDMRs_allReps_bins[unique(subjectHits(hypoDMRs_allReps_bins_genomeRegion_hits))]
  if(args$genomeRegion %in% c("arm")) {
    hypoDMRs_allReps_bins_genomeMask_hits <- findOverlaps(query = genomeMaskGR,
                                                          subject = hypoDMRs_allReps_bins,
                                                          type = "any", select = "all",
                                                          ignore.strand = TRUE)
    if(length(hypoDMRs_allReps_bins_genomeMask_hits) > 0) {
      hypoDMRs_allReps_bins <- hypoDMRs_allReps_bins[-subjectHits(hypoDMRs_allReps_bins_genomeMask_hits)]
    }
  }

  # Get CEN and nonCEN DMRs
  hypoDMRs_allReps_bins_CEN_hits <- findOverlaps(query = CENGR,
                                                 subject = hypoDMRs_allReps_bins,
                                                 type = "any", select = "all",
                                                 ignore.strand = TRUE)
  hypoDMRs_allReps_bins_CEN <- hypoDMRs_allReps_bins[unique(subjectHits(hypoDMRs_allReps_bins_CEN_hits))]
  if(length(hypoDMRs_allReps_bins_CEN_hits) > 0) {
    hypoDMRs_allReps_bins_nonCEN <- hypoDMRs_allReps_bins[-subjectHits(hypoDMRs_allReps_bins_CEN_hits)]
  } else {
    hypoDMRs_allReps_bins_nonCEN <- hypoDMRs_allReps_bins
  }

  stopifnot(length(hypoDMRs_allReps_bins) == length(hypoDMRs_allReps_bins_CEN) + length(hypoDMRs_allReps_bins_nonCEN))

  # Sort by seqnames, start and end
  hypoDMRs_allReps_bins <- sortSeqlevels(hypoDMRs_allReps_bins)
  hypoDMRs_allReps_bins <- sort(hypoDMRs_allReps_bins, ignore.strand = TRUE)
  hypoDMRs_allReps_bins_CEN <- sortSeqlevels(hypoDMRs_allReps_bins_CEN)
  hypoDMRs_allReps_bins_CEN <- sort(hypoDMRs_allReps_bins_CEN, ignore.strand = TRUE)
  hypoDMRs_allReps_bins_nonCEN <- sortSeqlevels(hypoDMRs_allReps_bins_nonCEN)
  hypoDMRs_allReps_bins_nonCEN <- sort(hypoDMRs_allReps_bins_nonCEN, ignore.strand = TRUE)


  # DMRs in args$genomeRegion
  # Create new columns containing absolute change, log2 fold change, and relative change in args$context methylation proportion
  hypoDMRs_allReps_bins$absolute_change <- as.numeric( hypoDMRs_allReps_bins$proportion1 - hypoDMRs_allReps_bins$proportion2 )
  # + 0.01 is an offset to prevent infinite log2 fold change (or equal relative change) values where proportions = 0
  hypoDMRs_allReps_bins$log2_fold_change <- as.numeric( log2( (hypoDMRs_allReps_bins$proportion1 + 0.01) /
                                                              (hypoDMRs_allReps_bins$proportion2 + 0.01) ) )
  hypoDMRs_allReps_bins$relative_change <- as.numeric( 1 - ( (hypoDMRs_allReps_bins$proportion2 + 0.01) /
                                                             (hypoDMRs_allReps_bins$proportion1 + 0.01) ) )
  
  # Create new columns containing absolute change, log2 fold change, and relative change percentiles and quantiles
  hypoDMRs_allReps_bins$ac_percentile <- as.numeric( rank(hypoDMRs_allReps_bins$absolute_change) /
                                                    length(hypoDMRs_allReps_bins$absolute_change) )
  hypoDMRs_allReps_bins$l2fc_percentile <- as.numeric( rank(hypoDMRs_allReps_bins$log2_fold_change) /
                                                      length(hypoDMRs_allReps_bins$log2_fold_change) )
  hypoDMRs_allReps_bins$rc_percentile <- as.numeric( rank(hypoDMRs_allReps_bins$relative_change) /
                                                    length(hypoDMRs_allReps_bins$relative_change) )
  hypoDMRs_allReps_bins$ac_quantile <- as.character("") 
  hypoDMRs_allReps_bins$l2fc_quantile <- as.character("") 
  hypoDMRs_allReps_bins$rc_quantile <- as.character("") 

  # Define quantiles
  ac_quantilesStats <- data.frame()
  l2fc_quantilesStats <- data.frame()
  rc_quantilesStats <- data.frame()
  for(k in 1:args$quantiles) {
    # absolute_change (ac)
    # First quantile should span 1 to greater than, e.g., 0.75 proportions of hypoDMRs_allReps_bins
    if(k < args$quantiles) {
       hypoDMRs_allReps_bins[ !is.na(hypoDMRs_allReps_bins$ac_percentile) &
                              hypoDMRs_allReps_bins$ac_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                              hypoDMRs_allReps_bins$ac_percentile > 1 - ( k / args$quantiles ) ]$ac_quantile <- paste0("Quantile ", k)
    } else { 
    # Final quantile should span 0 to, e.g., 0.25 proportions of hypoDMRs_allReps_bins
      hypoDMRs_allReps_bins[ !is.na(hypoDMRs_allReps_bins$ac_percentile) &
                             hypoDMRs_allReps_bins$ac_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                             hypoDMRs_allReps_bins$ac_percentile >= 1 - ( k / args$quantiles ) ]$ac_quantile <- paste0("Quantile ", k)
    }
    ac_stats <- data.frame(quantile = as.integer(k),
                           n = as.integer(length(hypoDMRs_allReps_bins[hypoDMRs_allReps_bins$ac_quantile == paste0("Quantile ", k)])),
                           mean_width = as.integer(round(mean(
                             width(hypoDMRs_allReps_bins[hypoDMRs_allReps_bins$ac_quantile == paste0("Quantile ", k)]), na.rm = T))),
                           total_width = as.integer(sum(
                             width(hypoDMRs_allReps_bins[hypoDMRs_allReps_bins$ac_quantile == paste0("Quantile ", k)]), na.rm = T)),
                           mean_absolute_change = as.numeric(mean(
                             hypoDMRs_allReps_bins[hypoDMRs_allReps_bins$ac_quantile == paste0("Quantile ", k)]$absolute_change, na.rm = T)),
                           mean_log2_fold_change = as.numeric(mean(
                             hypoDMRs_allReps_bins[hypoDMRs_allReps_bins$ac_quantile == paste0("Quantile ", k)]$log2_fold_change, na.rm = T)),
                           mean_relative_change = as.numeric(mean(
                             hypoDMRs_allReps_bins[hypoDMRs_allReps_bins$ac_quantile == paste0("Quantile ", k)]$relative_change, na.rm = T)),
                           stringsAsFactors = FALSE)
    ac_quantilesStats <- rbind(ac_quantilesStats, ac_stats)
    # log2_fold_change
    # First quantile should span 1 to greater than, e.g., 0.75 proportions of hypoDMRs_allReps_bins
    if(k < args$quantiles) {
       hypoDMRs_allReps_bins[ !is.na(hypoDMRs_allReps_bins$l2fc_percentile) &
                              hypoDMRs_allReps_bins$l2fc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                              hypoDMRs_allReps_bins$l2fc_percentile > 1 - ( k / args$quantiles ) ]$l2fc_quantile <- paste0("Quantile ", k)
    } else { 
    # Final quantile should span 0 to, e.g., 0.25 proportions of hypoDMRs_allReps_bins
      hypoDMRs_allReps_bins[ !is.na(hypoDMRs_allReps_bins$l2fc_percentile) &
                             hypoDMRs_allReps_bins$l2fc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                             hypoDMRs_allReps_bins$l2fc_percentile >= 1 - ( k / args$quantiles ) ]$l2fc_quantile <- paste0("Quantile ", k)
    }
    l2fc_stats <- data.frame(quantile = as.integer(k),
                             n = as.integer(length(hypoDMRs_allReps_bins[hypoDMRs_allReps_bins$l2fc_quantile == paste0("Quantile ", k)])),
                             mean_width = as.integer(round(mean(
                               width(hypoDMRs_allReps_bins[hypoDMRs_allReps_bins$l2fc_quantile == paste0("Quantile ", k)]), na.rm = T))),
                             total_width = as.integer(sum(
                               width(hypoDMRs_allReps_bins[hypoDMRs_allReps_bins$l2fc_quantile == paste0("Quantile ", k)]), na.rm = T)),
                             mean_absolute_change = as.numeric(mean(
                               hypoDMRs_allReps_bins[hypoDMRs_allReps_bins$l2fc_quantile == paste0("Quantile ", k)]$absolute_change, na.rm = T)),
                             mean_log2_fold_change = as.numeric(mean(
                               hypoDMRs_allReps_bins[hypoDMRs_allReps_bins$l2fc_quantile == paste0("Quantile ", k)]$log2_fold_change, na.rm = T)),
                             mean_relative_change = as.numeric(mean(
                               hypoDMRs_allReps_bins[hypoDMRs_allReps_bins$l2fc_quantile == paste0("Quantile ", k)]$relative_change, na.rm = T)),
                             stringsAsFactors = FALSE)
    l2fc_quantilesStats <- rbind(l2fc_quantilesStats, l2fc_stats)
    # relative_change
    # First quantile should span 1 to greater than, e.g., 0.75 proportions of hypoDMRs_allReps_bins
    if(k < args$quantiles) {
       hypoDMRs_allReps_bins[ !is.na(hypoDMRs_allReps_bins$rc_percentile) &
                              hypoDMRs_allReps_bins$rc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                              hypoDMRs_allReps_bins$rc_percentile > 1 - ( k / args$quantiles ) ]$rc_quantile <- paste0("Quantile ", k)
    } else { 
    # Final quantile should span 0 to, e.g., 0.25 proportions of hypoDMRs_allReps_bins
      hypoDMRs_allReps_bins[ !is.na(hypoDMRs_allReps_bins$rc_percentile) &
                             hypoDMRs_allReps_bins$rc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                             hypoDMRs_allReps_bins$rc_percentile >= 1 - ( k / args$quantiles ) ]$rc_quantile <- paste0("Quantile ", k)
    }
    rc_stats <- data.frame(quantile = as.integer(k),
                           n = as.integer(length(hypoDMRs_allReps_bins[hypoDMRs_allReps_bins$rc_quantile == paste0("Quantile ", k)])),
                           mean_width = as.integer(round(mean(
                             width(hypoDMRs_allReps_bins[hypoDMRs_allReps_bins$rc_quantile == paste0("Quantile ", k)]), na.rm = T))),
                           total_width = as.integer(sum(
                             width(hypoDMRs_allReps_bins[hypoDMRs_allReps_bins$rc_quantile == paste0("Quantile ", k)]), na.rm = T)),
                           mean_absolute_change = as.numeric(mean(
                             hypoDMRs_allReps_bins[hypoDMRs_allReps_bins$rc_quantile == paste0("Quantile ", k)]$absolute_change, na.rm = T)),
                           mean_log2_fold_change = as.numeric(mean(
                             hypoDMRs_allReps_bins[hypoDMRs_allReps_bins$rc_quantile == paste0("Quantile ", k)]$log2_fold_change, na.rm = T)),
                           mean_relative_change = as.numeric(mean(
                             hypoDMRs_allReps_bins[hypoDMRs_allReps_bins$rc_quantile == paste0("Quantile ", k)]$relative_change, na.rm = T)),
                           stringsAsFactors = FALSE)
    rc_quantilesStats <- rbind(rc_quantilesStats, rc_stats)
  }
  write.table(ac_quantilesStats,
              file = paste0(hypoDMRdir,
                            "summary_", args$quantiles, "quantiles",
                            "_by_absolute_change_in_",
                            paste0(args$condition2, collapse = "_"),
                            "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                            "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                            "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_", args$genomeRegion, ".tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(l2fc_quantilesStats,
              file = paste0(hypoDMRdir,
                            "summary_", args$quantiles, "quantiles",
                            "_by_log2_fold_change_in_",
                            paste0(args$condition2, collapse = "_"),
                            "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                            "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                            "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_", args$genomeRegion, ".tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(rc_quantilesStats,
              file = paste0(hypoDMRdir,
                            "summary_", args$quantiles, "quantiles",
                            "_by_relative_change_in_",
                            paste0(args$condition2, collapse = "_"),
                            "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                            "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                            "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_", args$genomeRegion, ".tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)

  # Export DMR GRanges as annotation files
  rtracklayer::export(object = hypoDMRs_allReps_bins,
                      con = paste0(hypoDMRdir,
                                   "features_", args$quantiles, "quantiles",
                                   "_by_change_in_",
                                   paste0(args$condition2, collapse = "_"),
                                   "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                                   "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                                   "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_", args$genomeRegion, ".gff3"))
  hypoDMRs_allReps_bins_bed <- data.frame(chr = as.character(seqnames(hypoDMRs_allReps_bins)),
                                          start = as.integer(start(hypoDMRs_allReps_bins)-1),
                                          end = as.integer(end(hypoDMRs_allReps_bins)),
                                          name = as.integer(1:length(hypoDMRs_allReps_bins)),
                                          score = as.numeric(hypoDMRs_allReps_bins$log2_fold_change),
                                          strand = as.character(strand(hypoDMRs_allReps_bins)),
                                          stringsAsFactors = FALSE)
  write.table(hypoDMRs_allReps_bins_bed,
              file = paste0(hypoDMRdir,
                            "features_", args$quantiles, "quantiles",
                            "_by_change_in_",
                            paste0(args$condition2, collapse = "_"),
                            "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                            "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                            "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_", args$genomeRegion, ".bed"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


  # Define random loci of the same number and width distribution,
  # and in the same per-feature genomeRegionGR as hypoDMRs_allReps_bins
  
  # Define function to select randomly positioned loci of the same
  # width distribution as hypoDMRs_allReps_bins
  ranLocStartSelect <- function(coordinates, n) {
    sample(x = coordinates,
           size = n,
           replace = FALSE)
  }
  
  # Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
  options(scipen = 100)
  
  # Apply ranLocStartSelect() on a per-chromosome basis so that
  # ranLocGR contains the same number of loci per chromosome as hypoDMRs_allReps_bins
  chrs <- args$chrName[which(args$chrName %in% seqnames(hypoDMRs_allReps_bins)@values)]
  hypoDMRs_allReps_bins_ranLocGR <- GRanges()
  for(i in 1:length(chrs)) {
    hypoDMRs_allReps_binsChrGR <- hypoDMRs_allReps_bins[seqnames(hypoDMRs_allReps_bins) == chrs[i]]
    genomeRegionChrGR <- genomeRegionGR[seqnames(genomeRegionGR) == chrs[i]]
    # Contract genomeRegionChrGR so that random loci and 2-kb flanking regions
    # do not extend beyond chromosome ends
    end(genomeRegionChrGR) <- end(genomeRegionChrGR)-max(width(hypoDMRs_allReps_binsChrGR))-2000
    start(genomeRegionChrGR) <- start(genomeRegionChrGR)+2000
    # Define seed so that random selections are reproducible
    set.seed(93750174)
    hypoDMRs_allReps_bins_ranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(genomeRegionChrGR), function(x) {
                                                                                     start(genomeRegionChrGR[x]) : end(genomeRegionChrGR[x])
                                                                                   })),
                                                              n = length(hypoDMRs_allReps_binsChrGR))
    hypoDMRs_allReps_bins_ranLocChrGR <- GRanges(seqnames = chrs[i],
                                                 ranges = IRanges(start = hypoDMRs_allReps_bins_ranLocChrStart,
                                                                  width = width(hypoDMRs_allReps_binsChrGR)),
                                                 strand = strand(hypoDMRs_allReps_binsChrGR))
    hypoDMRs_allReps_bins_ranLocGR <- append(hypoDMRs_allReps_bins_ranLocGR, hypoDMRs_allReps_bins_ranLocChrGR)
  }
  stopifnot( length( findOverlaps(query = hypoDMRs_allReps_bins_ranLocGR,
                                  subject = genomeMaskGR,
                                  type = "any", select = "all",
                                  ignore.strand = TRUE) ) == 0 )


  # Divide hypoDMRs_allReps_bins_ranLocGR into quantiles based on hypoDMRs_allReps_bins quantile indices
  hypoDMRs_allReps_bins_ranLocGR$l2fc_random <- as.character("")
  # Get row indices for each feature quantile
  l2fc_quantileIndices <- lapply(1:args$quantiles, function(k) {
    which(hypoDMRs_allReps_bins$l2fc_quantile == paste0("Quantile ", k))
  })
  for(k in 1:args$quantiles) {
    hypoDMRs_allReps_bins_ranLocGR[l2fc_quantileIndices[[k]]]$l2fc_random <- paste0("Random ", k)
  }

  hypoDMRs_allReps_bins_ranLocGR_bed <- data.frame(chr = as.character(seqnames(hypoDMRs_allReps_bins_ranLocGR)),
                                                   start = as.integer(start(hypoDMRs_allReps_bins_ranLocGR)-1),
                                                   end = as.integer(end(hypoDMRs_allReps_bins_ranLocGR)),
                                                   name = as.integer(1:length(hypoDMRs_allReps_bins_ranLocGR)),
                                                   score = as.character(hypoDMRs_allReps_bins_ranLocGR$l2fc_random),
                                                   strand = as.character(strand(hypoDMRs_allReps_bins_ranLocGR)),
                                                   stringsAsFactors = FALSE)
  write.table(hypoDMRs_allReps_bins_ranLocGR_bed,
              file = paste0(hypoDMRdir,
                            "features_", args$quantiles, "quantiles",
                            "_by_change_in_",
                            paste0(args$condition2, collapse = "_"),
                            "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                            "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                            "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_", args$genomeRegion, "_ranLoc.bed"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

  # DMRs in CEN
  if(length(hypoDMRs_allReps_bins_CEN) > 0) {
    # Create new columns containing absolute change, log2 fold change, and relative change in args$context methylation proportion
    hypoDMRs_allReps_bins_CEN$absolute_change <- as.numeric( hypoDMRs_allReps_bins_CEN$proportion1 - hypoDMRs_allReps_bins_CEN$proportion2 )
    # + 0.01 is an offset to prevent infinite log2 fold change (or equal relative change) values where proportions = 0
    hypoDMRs_allReps_bins_CEN$log2_fold_change <- as.numeric( log2( (hypoDMRs_allReps_bins_CEN$proportion1 + 0.01) /
                                                                    (hypoDMRs_allReps_bins_CEN$proportion2 + 0.01) ) )
    hypoDMRs_allReps_bins_CEN$relative_change <- as.numeric( 1 - ( (hypoDMRs_allReps_bins_CEN$proportion2 + 0.01) /
                                                                   (hypoDMRs_allReps_bins_CEN$proportion1 + 0.01) ) )
    
    # Create new columns containing absolute change, log2 fold change, and relative change percentiles and quantiles
    hypoDMRs_allReps_bins_CEN$ac_percentile <- as.numeric( rank(hypoDMRs_allReps_bins_CEN$absolute_change) /
                                                           length(hypoDMRs_allReps_bins_CEN$absolute_change) )
    hypoDMRs_allReps_bins_CEN$l2fc_percentile <- as.numeric( rank(hypoDMRs_allReps_bins_CEN$log2_fold_change) /
                                                             length(hypoDMRs_allReps_bins_CEN$log2_fold_change) )
    hypoDMRs_allReps_bins_CEN$rc_percentile <- as.numeric( rank(hypoDMRs_allReps_bins_CEN$relative_change) /
                                                           length(hypoDMRs_allReps_bins_CEN$relative_change) )
    hypoDMRs_allReps_bins_CEN$ac_quantile <- as.character("") 
    hypoDMRs_allReps_bins_CEN$l2fc_quantile <- as.character("") 
    hypoDMRs_allReps_bins_CEN$rc_quantile <- as.character("") 

    # Define quantiles
    ac_quantilesStats <- data.frame()
    l2fc_quantilesStats <- data.frame()
    rc_quantilesStats <- data.frame()
    for(k in 1:args$quantiles) {
      # absolute_change (ac)
      # First quantile should span 1 to greater than, e.g., 0.75 proportions of hypoDMRs_allReps_bins_CEN
      if(k < args$quantiles) {
         hypoDMRs_allReps_bins_CEN[ !is.na(hypoDMRs_allReps_bins_CEN$ac_percentile) &
                                    hypoDMRs_allReps_bins_CEN$ac_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                    hypoDMRs_allReps_bins_CEN$ac_percentile > 1 - ( k / args$quantiles ) ]$ac_quantile <- paste0("Quantile ", k)
      } else { 
      # Final quantile should span 0 to, e.g., 0.25 proportions of hypoDMRs_allReps_bins_CEN
        hypoDMRs_allReps_bins_CEN[ !is.na(hypoDMRs_allReps_bins_CEN$ac_percentile) &
                                   hypoDMRs_allReps_bins_CEN$ac_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                   hypoDMRs_allReps_bins_CEN$ac_percentile >= 1 - ( k / args$quantiles ) ]$ac_quantile <- paste0("Quantile ", k)
      }
      ac_stats <- data.frame(quantile = as.integer(k),
                             n = as.integer(length(hypoDMRs_allReps_bins_CEN[hypoDMRs_allReps_bins_CEN$ac_quantile == paste0("Quantile ", k)])),
                             mean_width = as.integer(round(mean(
                               width(hypoDMRs_allReps_bins_CEN[hypoDMRs_allReps_bins_CEN$ac_quantile == paste0("Quantile ", k)]), na.rm = T))),
                             total_width = as.integer(sum(
                               width(hypoDMRs_allReps_bins_CEN[hypoDMRs_allReps_bins_CEN$ac_quantile == paste0("Quantile ", k)]), na.rm = T)),
                             mean_absolute_change = as.numeric(mean(
                               hypoDMRs_allReps_bins_CEN[hypoDMRs_allReps_bins_CEN$ac_quantile == paste0("Quantile ", k)]$absolute_change, na.rm = T)),
                             mean_log2_fold_change = as.numeric(mean(
                               hypoDMRs_allReps_bins_CEN[hypoDMRs_allReps_bins_CEN$ac_quantile == paste0("Quantile ", k)]$log2_fold_change, na.rm = T)),
                             mean_relative_change = as.numeric(mean(
                               hypoDMRs_allReps_bins_CEN[hypoDMRs_allReps_bins_CEN$ac_quantile == paste0("Quantile ", k)]$relative_change, na.rm = T)),
                             stringsAsFactors = FALSE)
      ac_quantilesStats <- rbind(ac_quantilesStats, ac_stats)
      # log2_fold_change
      # First quantile should span 1 to greater than, e.g., 0.75 proportions of hypoDMRs_allReps_bins_CEN
      if(k < args$quantiles) {
         hypoDMRs_allReps_bins_CEN[ !is.na(hypoDMRs_allReps_bins_CEN$l2fc_percentile) &
                                    hypoDMRs_allReps_bins_CEN$l2fc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                    hypoDMRs_allReps_bins_CEN$l2fc_percentile > 1 - ( k / args$quantiles ) ]$l2fc_quantile <- paste0("Quantile ", k)
      } else { 
      # Final quantile should span 0 to, e.g., 0.25 proportions of hypoDMRs_allReps_bins_CEN
        hypoDMRs_allReps_bins_CEN[ !is.na(hypoDMRs_allReps_bins_CEN$l2fc_percentile) &
                                   hypoDMRs_allReps_bins_CEN$l2fc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                   hypoDMRs_allReps_bins_CEN$l2fc_percentile >= 1 - ( k / args$quantiles ) ]$l2fc_quantile <- paste0("Quantile ", k)
      }
      l2fc_stats <- data.frame(quantile = as.integer(k),
                               n = as.integer(length(hypoDMRs_allReps_bins_CEN[hypoDMRs_allReps_bins_CEN$l2fc_quantile == paste0("Quantile ", k)])),
                               mean_width = as.integer(round(mean(
                                 width(hypoDMRs_allReps_bins_CEN[hypoDMRs_allReps_bins_CEN$l2fc_quantile == paste0("Quantile ", k)]), na.rm = T))),
                               total_width = as.integer(sum(
                                 width(hypoDMRs_allReps_bins_CEN[hypoDMRs_allReps_bins_CEN$l2fc_quantile == paste0("Quantile ", k)]), na.rm = T)),
                               mean_absolute_change = as.numeric(mean(
                                 hypoDMRs_allReps_bins_CEN[hypoDMRs_allReps_bins_CEN$l2fc_quantile == paste0("Quantile ", k)]$absolute_change, na.rm = T)),
                               mean_log2_fold_change = as.numeric(mean(
                                 hypoDMRs_allReps_bins_CEN[hypoDMRs_allReps_bins_CEN$l2fc_quantile == paste0("Quantile ", k)]$log2_fold_change, na.rm = T)),
                               mean_relative_change = as.numeric(mean(
                                 hypoDMRs_allReps_bins_CEN[hypoDMRs_allReps_bins_CEN$l2fc_quantile == paste0("Quantile ", k)]$relative_change, na.rm = T)),
                               stringsAsFactors = FALSE)
      l2fc_quantilesStats <- rbind(l2fc_quantilesStats, l2fc_stats)
      # relative_change
      # First quantile should span 1 to greater than, e.g., 0.75 proportions of hypoDMRs_allReps_bins_CEN
      if(k < args$quantiles) {
         hypoDMRs_allReps_bins_CEN[ !is.na(hypoDMRs_allReps_bins_CEN$rc_percentile) &
                                    hypoDMRs_allReps_bins_CEN$rc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                    hypoDMRs_allReps_bins_CEN$rc_percentile > 1 - ( k / args$quantiles ) ]$rc_quantile <- paste0("Quantile ", k)
      } else { 
      # Final quantile should span 0 to, e.g., 0.25 proportions of hypoDMRs_allReps_bins_CEN
        hypoDMRs_allReps_bins_CEN[ !is.na(hypoDMRs_allReps_bins_CEN$rc_percentile) &
                                   hypoDMRs_allReps_bins_CEN$rc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                   hypoDMRs_allReps_bins_CEN$rc_percentile >= 1 - ( k / args$quantiles ) ]$rc_quantile <- paste0("Quantile ", k)
      }
      rc_stats <- data.frame(quantile = as.integer(k),
                             n = as.integer(length(hypoDMRs_allReps_bins_CEN[hypoDMRs_allReps_bins_CEN$rc_quantile == paste0("Quantile ", k)])),
                             mean_width = as.integer(round(mean(
                               width(hypoDMRs_allReps_bins_CEN[hypoDMRs_allReps_bins_CEN$rc_quantile == paste0("Quantile ", k)]), na.rm = T))),
                             total_width = as.integer(sum(
                               width(hypoDMRs_allReps_bins_CEN[hypoDMRs_allReps_bins_CEN$rc_quantile == paste0("Quantile ", k)]), na.rm = T)),
                             mean_absolute_change = as.numeric(mean(
                               hypoDMRs_allReps_bins_CEN[hypoDMRs_allReps_bins_CEN$rc_quantile == paste0("Quantile ", k)]$absolute_change, na.rm = T)),
                             mean_log2_fold_change = as.numeric(mean(
                               hypoDMRs_allReps_bins_CEN[hypoDMRs_allReps_bins_CEN$rc_quantile == paste0("Quantile ", k)]$log2_fold_change, na.rm = T)),
                             mean_relative_change = as.numeric(mean(
                               hypoDMRs_allReps_bins_CEN[hypoDMRs_allReps_bins_CEN$rc_quantile == paste0("Quantile ", k)]$relative_change, na.rm = T)),
                             stringsAsFactors = FALSE)
      rc_quantilesStats <- rbind(rc_quantilesStats, rc_stats)
    }
    write.table(ac_quantilesStats,
                file = paste0(hypoDMRdir,
                              "summary_", args$quantiles, "quantiles",
                              "_by_absolute_change_in_",
                              paste0(args$condition2, collapse = "_"),
                              "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                              "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                              "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_CEN.tsv"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(l2fc_quantilesStats,
                file = paste0(hypoDMRdir,
                              "summary_", args$quantiles, "quantiles",
                              "_by_log2_fold_change_in_",
                              paste0(args$condition2, collapse = "_"),
                              "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                              "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                              "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_CEN.tsv"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(rc_quantilesStats,
                file = paste0(hypoDMRdir,
                              "summary_", args$quantiles, "quantiles",
                              "_by_relative_change_in_",
                              paste0(args$condition2, collapse = "_"),
                              "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                              "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                              "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_CEN.tsv"),
                quote = FALSE, sep = "\t", row.names = FALSE)

    # Export DMR GRanges as annotation files
    rtracklayer::export(object = hypoDMRs_allReps_bins_CEN,
                        con = paste0(hypoDMRdir,
                                     "features_", args$quantiles, "quantiles",
                                     "_by_change_in_",
                                     paste0(args$condition2, collapse = "_"),
                                     "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                                     "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                                     "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_CEN.gff3"))
    hypoDMRs_allReps_bins_CEN_bed <- data.frame(chr = as.character(seqnames(hypoDMRs_allReps_bins_CEN)),
                                                start = as.integer(start(hypoDMRs_allReps_bins_CEN)-1),
                                                end = as.integer(end(hypoDMRs_allReps_bins_CEN)),
                                                name = as.integer(1:length(hypoDMRs_allReps_bins_CEN)),
                                                score = as.numeric(hypoDMRs_allReps_bins_CEN$log2_fold_change),
                                                strand = as.character(strand(hypoDMRs_allReps_bins_CEN)),
                                                stringsAsFactors = FALSE)
    write.table(hypoDMRs_allReps_bins_CEN_bed,
                file = paste0(hypoDMRdir,
                              "features_", args$quantiles, "quantiles",
                              "_by_change_in_",
                              paste0(args$condition2, collapse = "_"),
                              "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                              "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                              "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_CEN.bed"),
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


    # Define random loci of the same number and width distribution,
    # and in the same per-feature CENGR as hypoDMRs_allReps_bins_CEN
    
    # Define function to select randomly positioned loci of the same
    # width distribution as hypoDMRs_allReps_bins_CEN
    ranLocStartSelect <- function(coordinates, n) {
      sample(x = coordinates,
             size = n,
             replace = FALSE)
    }
    
    # Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
    options(scipen = 100)
    
    # Apply ranLocStartSelect() on a per-chromosome basis so that
    # ranLocGR contains the same number of loci per chromosome as hypoDMRs_allReps_bins_CEN
    chrs <- args$chrName[which(args$chrName %in% seqnames(hypoDMRs_allReps_bins_CEN)@values)]
    hypoDMRs_allReps_bins_CEN_ranLocGR <- GRanges()
    for(i in 1:length(chrs)) {
      hypoDMRs_allReps_bins_CENChrGR <- hypoDMRs_allReps_bins_CEN[seqnames(hypoDMRs_allReps_bins_CEN) == chrs[i]]
      CENChrGR <- CENGR[seqnames(CENGR) == chrs[i]]
      # Contract CENChrGR so that random loci and 2-kb flanking regions
      # do not extend beyond chromosome ends
      end(CENChrGR) <- end(CENChrGR)-max(width(hypoDMRs_allReps_bins_CENChrGR))-2000
      start(CENChrGR) <- start(CENChrGR)+2000
      # Define seed so that random selections are reproducible
      set.seed(93750174)
      hypoDMRs_allReps_bins_CEN_ranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(CENChrGR), function(x) {
                                                                                           start(CENChrGR[x]) : end(CENChrGR[x])
                                                                                         })),
                                                                    n = length(hypoDMRs_allReps_bins_CENChrGR))
      hypoDMRs_allReps_bins_CEN_ranLocChrGR <- GRanges(seqnames = chrs[i],
                                                       ranges = IRanges(start = hypoDMRs_allReps_bins_CEN_ranLocChrStart,
                                                                        width = width(hypoDMRs_allReps_bins_CENChrGR)),
                                                       strand = strand(hypoDMRs_allReps_bins_CENChrGR))
      hypoDMRs_allReps_bins_CEN_ranLocGR <- append(hypoDMRs_allReps_bins_CEN_ranLocGR, hypoDMRs_allReps_bins_CEN_ranLocChrGR)
    }
    #stopifnot( length( findOverlaps(query = hypoDMRs_allReps_bins_CEN_ranLocGR,
    #                                subject = nonCENGR,
    #                                type = "any", select = "all",
    #                                ignore.strand = TRUE) ) == 0 )


    # Divide hypoDMRs_allReps_bins_CEN_ranLocGR into quantiles based on hypoDMRs_allReps_bins_CEN quantile indices
    hypoDMRs_allReps_bins_CEN_ranLocGR$l2fc_random <- as.character("")
    # Get row indices for each feature quantile
    l2fc_quantileIndices <- lapply(1:args$quantiles, function(k) {
      which(hypoDMRs_allReps_bins_CEN$l2fc_quantile == paste0("Quantile ", k))
    })
    for(k in 1:args$quantiles) {
      hypoDMRs_allReps_bins_CEN_ranLocGR[l2fc_quantileIndices[[k]]]$l2fc_random <- paste0("Random ", k)
    }

    hypoDMRs_allReps_bins_CEN_ranLocGR_bed <- data.frame(chr = as.character(seqnames(hypoDMRs_allReps_bins_CEN_ranLocGR)),
                                                         start = as.integer(start(hypoDMRs_allReps_bins_CEN_ranLocGR)-1),
                                                         end = as.integer(end(hypoDMRs_allReps_bins_CEN_ranLocGR)),
                                                         name = as.integer(1:length(hypoDMRs_allReps_bins_CEN_ranLocGR)),
                                                         score = as.character(hypoDMRs_allReps_bins_CEN_ranLocGR$l2fc_random),
                                                         strand = as.character(strand(hypoDMRs_allReps_bins_CEN_ranLocGR)),
                                                         stringsAsFactors = FALSE)
    write.table(hypoDMRs_allReps_bins_CEN_ranLocGR_bed,
                file = paste0(hypoDMRdir,
                              "features_", args$quantiles, "quantiles",
                              "_by_change_in_",
                              paste0(args$condition2, collapse = "_"),
                              "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                              "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                              "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_CEN_ranLoc.bed"),
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }

  # DMRs in nonCEN
  if(length(hypoDMRs_allReps_bins_nonCEN) > 0) {
    # Create new columns containing absolute change, log2 fold change, and relative change in args$context methylation proportion
    hypoDMRs_allReps_bins_nonCEN$absolute_change <- as.numeric( hypoDMRs_allReps_bins_nonCEN$proportion1 - hypoDMRs_allReps_bins_nonCEN$proportion2 )
    # + 0.01 is an offset to prevent infinite log2 fold change (or equal relative change) values where proportions = 0
    hypoDMRs_allReps_bins_nonCEN$log2_fold_change <- as.numeric( log2( (hypoDMRs_allReps_bins_nonCEN$proportion1 + 0.01) /
                                                                       (hypoDMRs_allReps_bins_nonCEN$proportion2 + 0.01) ) )
    hypoDMRs_allReps_bins_nonCEN$relative_change <- as.numeric( 1 - ( (hypoDMRs_allReps_bins_nonCEN$proportion2 + 0.01) /
                                                                      (hypoDMRs_allReps_bins_nonCEN$proportion1 + 0.01) ) )
    
    # Create new columns containing absolute change, log2 fold change, and relative change percentiles and quantiles
    hypoDMRs_allReps_bins_nonCEN$ac_percentile <- as.numeric( rank(hypoDMRs_allReps_bins_nonCEN$absolute_change) /
                                                              length(hypoDMRs_allReps_bins_nonCEN$absolute_change) )
    hypoDMRs_allReps_bins_nonCEN$l2fc_percentile <- as.numeric( rank(hypoDMRs_allReps_bins_nonCEN$log2_fold_change) /
                                                                length(hypoDMRs_allReps_bins_nonCEN$log2_fold_change) )
    hypoDMRs_allReps_bins_nonCEN$rc_percentile <- as.numeric( rank(hypoDMRs_allReps_bins_nonCEN$relative_change) /
                                                              length(hypoDMRs_allReps_bins_nonCEN$relative_change) )
    hypoDMRs_allReps_bins_nonCEN$ac_quantile <- as.character("") 
    hypoDMRs_allReps_bins_nonCEN$l2fc_quantile <- as.character("") 
    hypoDMRs_allReps_bins_nonCEN$rc_quantile <- as.character("") 

    # Define quantiles
    ac_quantilesStats <- data.frame()
    l2fc_quantilesStats <- data.frame()
    rc_quantilesStats <- data.frame()
    for(k in 1:args$quantiles) {
      # absolute_change (ac)
      # First quantile should span 1 to greater than, e.g., 0.75 proportions of hypoDMRs_allReps_bins_nonCEN
      if(k < args$quantiles) {
         hypoDMRs_allReps_bins_nonCEN[ !is.na(hypoDMRs_allReps_bins_nonCEN$ac_percentile) &
                                       hypoDMRs_allReps_bins_nonCEN$ac_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                       hypoDMRs_allReps_bins_nonCEN$ac_percentile > 1 - ( k / args$quantiles ) ]$ac_quantile <- paste0("Quantile ", k)
      } else { 
      # Final quantile should span 0 to, e.g., 0.25 proportions of hypoDMRs_allReps_bins_nonCEN
        hypoDMRs_allReps_bins_nonCEN[ !is.na(hypoDMRs_allReps_bins_nonCEN$ac_percentile) &
                                      hypoDMRs_allReps_bins_nonCEN$ac_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                      hypoDMRs_allReps_bins_nonCEN$ac_percentile >= 1 - ( k / args$quantiles ) ]$ac_quantile <- paste0("Quantile ", k)
      }
      ac_stats <- data.frame(quantile = as.integer(k),
                             n = as.integer(length(hypoDMRs_allReps_bins_nonCEN[hypoDMRs_allReps_bins_nonCEN$ac_quantile == paste0("Quantile ", k)])),
                             mean_width = as.integer(round(mean(
                               width(hypoDMRs_allReps_bins_nonCEN[hypoDMRs_allReps_bins_nonCEN$ac_quantile == paste0("Quantile ", k)]), na.rm = T))),
                             total_width = as.integer(sum(
                               width(hypoDMRs_allReps_bins_nonCEN[hypoDMRs_allReps_bins_nonCEN$ac_quantile == paste0("Quantile ", k)]), na.rm = T)),
                             mean_absolute_change = as.numeric(mean(
                               hypoDMRs_allReps_bins_nonCEN[hypoDMRs_allReps_bins_nonCEN$ac_quantile == paste0("Quantile ", k)]$absolute_change, na.rm = T)),
                             mean_log2_fold_change = as.numeric(mean(
                               hypoDMRs_allReps_bins_nonCEN[hypoDMRs_allReps_bins_nonCEN$ac_quantile == paste0("Quantile ", k)]$log2_fold_change, na.rm = T)),
                             mean_relative_change = as.numeric(mean(
                               hypoDMRs_allReps_bins_nonCEN[hypoDMRs_allReps_bins_nonCEN$ac_quantile == paste0("Quantile ", k)]$relative_change, na.rm = T)),
                             stringsAsFactors = FALSE)
      ac_quantilesStats <- rbind(ac_quantilesStats, ac_stats)
      # log2_fold_change
      # First quantile should span 1 to greater than, e.g., 0.75 proportions of hypoDMRs_allReps_bins_nonCEN
      if(k < args$quantiles) {
         hypoDMRs_allReps_bins_nonCEN[ !is.na(hypoDMRs_allReps_bins_nonCEN$l2fc_percentile) &
                                       hypoDMRs_allReps_bins_nonCEN$l2fc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                       hypoDMRs_allReps_bins_nonCEN$l2fc_percentile > 1 - ( k / args$quantiles ) ]$l2fc_quantile <- paste0("Quantile ", k)
      } else { 
      # Final quantile should span 0 to, e.g., 0.25 proportions of hypoDMRs_allReps_bins_nonCEN
        hypoDMRs_allReps_bins_nonCEN[ !is.na(hypoDMRs_allReps_bins_nonCEN$l2fc_percentile) &
                                      hypoDMRs_allReps_bins_nonCEN$l2fc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                      hypoDMRs_allReps_bins_nonCEN$l2fc_percentile >= 1 - ( k / args$quantiles ) ]$l2fc_quantile <- paste0("Quantile ", k)
      }
      l2fc_stats <- data.frame(quantile = as.integer(k),
                               n = as.integer(length(hypoDMRs_allReps_bins_nonCEN[hypoDMRs_allReps_bins_nonCEN$l2fc_quantile == paste0("Quantile ", k)])),
                               mean_width = as.integer(round(mean(
                                 width(hypoDMRs_allReps_bins_nonCEN[hypoDMRs_allReps_bins_nonCEN$l2fc_quantile == paste0("Quantile ", k)]), na.rm = T))),
                               total_width = as.integer(sum(
                                 width(hypoDMRs_allReps_bins_nonCEN[hypoDMRs_allReps_bins_nonCEN$l2fc_quantile == paste0("Quantile ", k)]), na.rm = T)),
                               mean_absolute_change = as.numeric(mean(
                                 hypoDMRs_allReps_bins_nonCEN[hypoDMRs_allReps_bins_nonCEN$l2fc_quantile == paste0("Quantile ", k)]$absolute_change, na.rm = T)),
                               mean_log2_fold_change = as.numeric(mean(
                                 hypoDMRs_allReps_bins_nonCEN[hypoDMRs_allReps_bins_nonCEN$l2fc_quantile == paste0("Quantile ", k)]$log2_fold_change, na.rm = T)),
                               mean_relative_change = as.numeric(mean(
                                 hypoDMRs_allReps_bins_nonCEN[hypoDMRs_allReps_bins_nonCEN$l2fc_quantile == paste0("Quantile ", k)]$relative_change, na.rm = T)),
                               stringsAsFactors = FALSE)
      l2fc_quantilesStats <- rbind(l2fc_quantilesStats, l2fc_stats)
      # relative_change
      # First quantile should span 1 to greater than, e.g., 0.75 proportions of hypoDMRs_allReps_bins_nonCEN
      if(k < args$quantiles) {
         hypoDMRs_allReps_bins_nonCEN[ !is.na(hypoDMRs_allReps_bins_nonCEN$rc_percentile) &
                                       hypoDMRs_allReps_bins_nonCEN$rc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                       hypoDMRs_allReps_bins_nonCEN$rc_percentile > 1 - ( k / args$quantiles ) ]$rc_quantile <- paste0("Quantile ", k)
      } else { 
      # Final quantile should span 0 to, e.g., 0.25 proportions of hypoDMRs_allReps_bins_nonCEN
        hypoDMRs_allReps_bins_nonCEN[ !is.na(hypoDMRs_allReps_bins_nonCEN$rc_percentile) &
                                      hypoDMRs_allReps_bins_nonCEN$rc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                      hypoDMRs_allReps_bins_nonCEN$rc_percentile >= 1 - ( k / args$quantiles ) ]$rc_quantile <- paste0("Quantile ", k)
      }
      rc_stats <- data.frame(quantile = as.integer(k),
                             n = as.integer(length(hypoDMRs_allReps_bins_nonCEN[hypoDMRs_allReps_bins_nonCEN$rc_quantile == paste0("Quantile ", k)])),
                             mean_width = as.integer(round(mean(
                               width(hypoDMRs_allReps_bins_nonCEN[hypoDMRs_allReps_bins_nonCEN$rc_quantile == paste0("Quantile ", k)]), na.rm = T))),
                             total_width = as.integer(sum(
                               width(hypoDMRs_allReps_bins_nonCEN[hypoDMRs_allReps_bins_nonCEN$rc_quantile == paste0("Quantile ", k)]), na.rm = T)),
                             mean_absolute_change = as.numeric(mean(
                               hypoDMRs_allReps_bins_nonCEN[hypoDMRs_allReps_bins_nonCEN$rc_quantile == paste0("Quantile ", k)]$absolute_change, na.rm = T)),
                             mean_log2_fold_change = as.numeric(mean(
                               hypoDMRs_allReps_bins_nonCEN[hypoDMRs_allReps_bins_nonCEN$rc_quantile == paste0("Quantile ", k)]$log2_fold_change, na.rm = T)),
                             mean_relative_change = as.numeric(mean(
                               hypoDMRs_allReps_bins_nonCEN[hypoDMRs_allReps_bins_nonCEN$rc_quantile == paste0("Quantile ", k)]$relative_change, na.rm = T)),
                             stringsAsFactors = FALSE)
      rc_quantilesStats <- rbind(rc_quantilesStats, rc_stats)
    }
    write.table(ac_quantilesStats,
                file = paste0(hypoDMRdir,
                              "summary_", args$quantiles, "quantiles",
                              "_by_absolute_change_in_",
                              paste0(args$condition2, collapse = "_"),
                              "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                              "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                              "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_nonCEN.tsv"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(l2fc_quantilesStats,
                file = paste0(hypoDMRdir,
                              "summary_", args$quantiles, "quantiles",
                              "_by_log2_fold_change_in_",
                              paste0(args$condition2, collapse = "_"),
                              "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                              "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                              "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_nonCEN.tsv"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(rc_quantilesStats,
                file = paste0(hypoDMRdir,
                              "summary_", args$quantiles, "quantiles",
                              "_by_relative_change_in_",
                              paste0(args$condition2, collapse = "_"),
                              "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                              "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                              "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_nonCEN.tsv"),
                quote = FALSE, sep = "\t", row.names = FALSE)

    # Export DMR GRanges as annotation files
    rtracklayer::export(object = hypoDMRs_allReps_bins_nonCEN,
                        con = paste0(hypoDMRdir,
                                     "features_", args$quantiles, "quantiles",
                                     "_by_change_in_",
                                     paste0(args$condition2, collapse = "_"),
                                     "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                                     "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                                     "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_nonCEN.gff3"))
    hypoDMRs_allReps_bins_nonCEN_bed <- data.frame(chr = as.character(seqnames(hypoDMRs_allReps_bins_nonCEN)),
                                                   start = as.integer(start(hypoDMRs_allReps_bins_nonCEN)-1),
                                                   end = as.integer(end(hypoDMRs_allReps_bins_nonCEN)),
                                                   name = as.integer(1:length(hypoDMRs_allReps_bins_nonCEN)),
                                                   score = as.numeric(hypoDMRs_allReps_bins_nonCEN$log2_fold_change),
                                                   strand = as.character(strand(hypoDMRs_allReps_bins_nonCEN)),
                                                   stringsAsFactors = FALSE)
    write.table(hypoDMRs_allReps_bins_nonCEN_bed,
                file = paste0(hypoDMRdir,
                              "features_", args$quantiles, "quantiles",
                              "_by_change_in_",
                              paste0(args$condition2, collapse = "_"),
                              "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                              "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                              "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_nonCEN.bed"),
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


    # Define random loci of the same number and width distribution,
    # and in the same per-feature nonCENGR as hypoDMRs_allReps_bins_nonCEN
    
    # Define function to select randomly positioned loci of the same
    # width distribution as hypoDMRs_allReps_bins_nonCEN
    ranLocStartSelect <- function(coordinates, n) {
      sample(x = coordinates,
             size = n,
             replace = FALSE)
    }
    
    # Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
    options(scipen = 100)
    
    # Apply ranLocStartSelect() on a per-chromosome basis so that
    # ranLocGR contains the same number of loci per chromosome as hypoDMRs_allReps_bins_nonCEN
    chrs <- args$chrName[which(args$chrName %in% seqnames(hypoDMRs_allReps_bins_nonCEN)@values)]
    hypoDMRs_allReps_bins_nonCEN_ranLocGR <- GRanges()
    for(i in 1:length(chrs)) {
      hypoDMRs_allReps_bins_nonCENChrGR <- hypoDMRs_allReps_bins_nonCEN[seqnames(hypoDMRs_allReps_bins_nonCEN) == chrs[i]]
      nonCENChrGR <- nonCENGR[seqnames(nonCENGR) == chrs[i]]
      # Contract nonCENChrGR so that random loci and 2-kb flanking regions
      # do not extend beyond chromosome ends
      end(nonCENChrGR) <- end(nonCENChrGR)-max(width(hypoDMRs_allReps_bins_nonCENChrGR))-2000
      start(nonCENChrGR) <- start(nonCENChrGR)+2000
      # Define seed so that random selections are reproducible
      set.seed(93750174)
      hypoDMRs_allReps_bins_nonCEN_ranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(nonCENChrGR), function(x) {
                                                                                           start(nonCENChrGR[x]) : end(nonCENChrGR[x])
                                                                                         })),
                                                                    n = length(hypoDMRs_allReps_bins_nonCENChrGR))
      hypoDMRs_allReps_bins_nonCEN_ranLocChrGR <- GRanges(seqnames = chrs[i],
                                                       ranges = IRanges(start = hypoDMRs_allReps_bins_nonCEN_ranLocChrStart,
                                                                        width = width(hypoDMRs_allReps_bins_nonCENChrGR)),
                                                       strand = strand(hypoDMRs_allReps_bins_nonCENChrGR))
      hypoDMRs_allReps_bins_nonCEN_ranLocGR <- append(hypoDMRs_allReps_bins_nonCEN_ranLocGR, hypoDMRs_allReps_bins_nonCEN_ranLocChrGR)
    }
    #stopifnot( length( findOverlaps(query = hypoDMRs_allReps_bins_nonCEN_ranLocGR,
    #                                subject = CENGR,
    #                                type = "any", select = "all",
    #                                ignore.strand = TRUE) ) == 0 )


    # Divide hypoDMRs_allReps_bins_nonCEN_ranLocGR into quantiles based on hypoDMRs_allReps_bins_nonCEN quantile indices
    hypoDMRs_allReps_bins_nonCEN_ranLocGR$l2fc_random <- as.character("")
    # Get row indices for each feature quantile
    l2fc_quantileIndices <- lapply(1:args$quantiles, function(k) {
      which(hypoDMRs_allReps_bins_nonCEN$l2fc_quantile == paste0("Quantile ", k))
    })
    for(k in 1:args$quantiles) {
      hypoDMRs_allReps_bins_nonCEN_ranLocGR[l2fc_quantileIndices[[k]]]$l2fc_random <- paste0("Random ", k)
    }

    hypoDMRs_allReps_bins_nonCEN_ranLocGR_bed <- data.frame(chr = as.character(seqnames(hypoDMRs_allReps_bins_nonCEN_ranLocGR)),
                                                            start = as.integer(start(hypoDMRs_allReps_bins_nonCEN_ranLocGR)-1),
                                                            end = as.integer(end(hypoDMRs_allReps_bins_nonCEN_ranLocGR)),
                                                            name = as.integer(1:length(hypoDMRs_allReps_bins_nonCEN_ranLocGR)),
                                                            score = as.character(hypoDMRs_allReps_bins_nonCEN_ranLocGR$l2fc_random),
                                                            strand = as.character(strand(hypoDMRs_allReps_bins_nonCEN_ranLocGR)),
                                                            stringsAsFactors = FALSE)
    write.table(hypoDMRs_allReps_bins_nonCEN_ranLocGR_bed,
                file = paste0(hypoDMRdir,
                              "features_", args$quantiles, "quantiles",
                              "_by_change_in_",
                              paste0(args$condition2, collapse = "_"),
                              "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                              "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                              "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_nonCEN_ranLoc.bed"),
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }

}

#} else {
#  # For biological replicate analysis using computeDMRsReplicates:
#  # (Note that this requires equal numbers of replicates for each condition,
#  # or > 1 replicate for each condition. If this is not the case,
#  # computeDMRsReplicates fails with an error message:
#  # "Error in apply(readsM2, 1, sum) : dim(X) must have a positive length")
#  conditions_Reps_list <- c(condition1_Reps, condition2_Reps)
#  
#  # Combine replicates into a single GRanges object containing
#  # data for each condition and replicate
#  joined_Reps <- joinReplicates(methylationData1 = conditions_Reps_list[[1]],
#                                methylationData2 = conditions_Reps_list[[2]],
#                                usecomplete = FALSE)
#  for(x in 3:length(conditions_Reps_list)) {
#    joined_Reps <- joinReplicates(methylationData1 = joined_Reps,
#                                  methylationData2 = conditions_Reps_list[[x]]) 
#  }
#  
#  ## Get ranges corresponding to the given context
#  #joined_Reps <- joined_Reps[joined_Reps$context == sub("p", "", args$context)]
#  #
#  ## Get ranges corresponding to those in chrName
#  ##joined_Reps <- joined_Reps[seqnames(joined_Reps) %in% args$chrName]
#  #joined_Reps <- keepSeqlevels(joined_Reps, args$chrName, pruning.mode="coarse")
#  #
#  ## Sort by seqnames, start and end
#  #joined_Reps <- sortSeqlevels(joined_Reps)
#  #joined_Reps <- sort(joined_Reps, ignore.strand = TRUE)
#  #
#  #print(joined_Reps)
#  
#  
#  # Create condition vector
#  joined_Reps_conditions <- gsub("_.+", "", c(args$condition1, args$condition2))
#  
#  print("joined_Reps conditions:")
#  print(joined_Reps_conditions)
#  
#  # Compute DMRs using "bins" method
#  DMRsReplicates_bins <- computeDMRsReplicates(methylationData = joined_Reps,
#                                               condition = joined_Reps_conditions,
#                                               regions = NULL,
#                                               context = sub("p", "", args$context),
#                                               method = "bins",
#                                               binSize = 100,
#                                               test = "betareg",
#                                               pseudocountM = 1,
#                                               pseudocountN = 2,
#                                               pValueThreshold = 0.01,
#                                               minCytosinesCount = 4,
#                                               minProportionDifference = minProportionDifference_context,
#                                               minGap = 200,
#                                               minSize = 50,
#                                               minReadsPerCytosine = 4,
#                                               cores = 1)
#  
#}
#
#
## For pooled-replicate analysis
#condition1_Reps_pooled <- poolMethylationDatasets(GRangesList(condition1_Reps))
#condition2_Reps_pooled <- poolMethylationDatasets(GRangesList(condition2_Reps))
