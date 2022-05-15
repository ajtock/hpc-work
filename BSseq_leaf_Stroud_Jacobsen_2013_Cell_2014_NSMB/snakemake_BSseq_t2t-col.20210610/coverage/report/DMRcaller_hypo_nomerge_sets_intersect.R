#!/usr/bin/env Rscript

# Get intersection of DMRs for given sets

# Usage:
# source ~/.bashrc
# conda activate python_3.9.6
# ./DMRcaller_hypo_nomerge_sets_intersect.R --sets 'cmt3_BSseq_Rep1,met1_BSseq_Rep1' \
#                                           --refbase t2t-col.20210610 \
#                                           --chrName 'Chr1,Chr2,Chr3,Chr4,Chr5' \
#                                           --genomeRegion genomewide \
#                                           --quantiles 6 \
#                                           --context CHG
# conda deactivate

library(argparse)
library(DMRcaller)
library(betareg) # required by DMRcaller::computeDMRsReplicates 
library(rtracklayer)

# Create parser object
parser <- ArgumentParser()

# Specify arguments
# ArgumentParser will add a help option by default
parser$add_argument("--sets", type = "character",
                    help="Sample replicate prefixes for sets to be intersected.")
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
##system("./DMRcaller_hypo_nomerge_sets_intersect.R --sets 'cmt3_BSseq_Rep1,met1_BSseq_Rep1' --refbase 't2t-col.20210610' --chrName 'Chr1,Chr2,Chr3,Chr4,Chr5' --genomeRegion genomewide --quantiles 6 --context CHG")
#args <- readRDS(args_file)

args$sets <- unlist(strsplit(args$sets, split = ","))
args$chrName <- unlist(strsplit(args$chrName, split = ","))

print(args)


# Define output directories
hypoDMRdir <- paste0("DMRs/hypoDMRs/", paste0(args$chrName, collapse = "_"), "/")
#hyperDMRdir <- paste0("DMRs/hyperDMRs/", paste0(args$chrName, collapse = "_"), "/")
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

hypoDMRs_perSet_list_bins <- lapply(1:length(args$sets), function(x) {
  # Load DMRs and convert into GRanges
  GRanges(rtracklayer::readGFF(paste0(hypoDMRdir,
                                      "features_", args$quantiles, "quantiles",
                                      "_by_change_in_",
                                      paste0(args$sets[x], collapse = "_"),
                                      "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                                      "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                                      "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_", args$genomeRegion, ".gff3")))
})

hypoDMRs_set1Only_bins <- GRanges()
hypoDMRs_allSets_bins <- hypoDMRs_perSet_list_bins[[1]]
for(x in 2:length(hypoDMRs_perSet_list_bins)) {
  hits <- findOverlaps(query = hypoDMRs_allSets_bins,
                       subject = hypoDMRs_perSet_list_bins[[x]],
                       type = "equal", select = "all",
                       ignore.strand = FALSE)
  hypoDMRs_set1Only_bins <- hypoDMRs_allSets_bins[-queryHits(hits)]
  hypoDMRs_allSets_bins <- hypoDMRs_allSets_bins[unique(queryHits(hits))]
}
stopifnot(length(hypoDMRs_set1Only_bins) + length(hypoDMRs_allSets_bins) == length(hypoDMRs_perSet_list_bins[[1]]))

# Get args$genomeRegion DMRs
hypoDMRs_allSets_bins_genomeRegion_hits <- findOverlaps(query = genomeRegionGR,
                                                        subject = hypoDMRs_allSets_bins,
                                                        type = "any", select = "all",
                                                        ignore.strand = TRUE)
hypoDMRs_allSets_bins <- hypoDMRs_allSets_bins[unique(subjectHits(hypoDMRs_allSets_bins_genomeRegion_hits))]
if(args$genomeRegion %in% c("arm")) {
  hypoDMRs_allSets_bins_genomeMask_hits <- findOverlaps(query = genomeMaskGR,
                                                        subject = hypoDMRs_allSets_bins,
                                                        type = "any", select = "all",
                                                        ignore.strand = TRUE)
  if(length(hypoDMRs_allSets_bins_genomeMask_hits) > 0) {
    hypoDMRs_allSets_bins <- hypoDMRs_allSets_bins[-subjectHits(hypoDMRs_allSets_bins_genomeMask_hits)]
  }
}

# Sort by seqnames, start and end
hypoDMRs_allSets_bins <- sortSeqlevels(hypoDMRs_allSets_bins)
hypoDMRs_allSets_bins <- sort(hypoDMRs_allSets_bins, ignore.strand = TRUE)

# DMRs in args$genomeRegion
# Create new columns containing absolute change, log2 fold change, and relative change in args$context methylation proportion
hypoDMRs_allSets_bins$absolute_change <- as.numeric( hypoDMRs_allSets_bins$proportion1 - hypoDMRs_allSets_bins$proportion2 )
# + 0.01 is an offset to prevent infinite log2 fold change (or equal relative change) values where proportions = 0
hypoDMRs_allSets_bins$log2_fold_change <- as.numeric( log2( (hypoDMRs_allSets_bins$proportion1 + 0.01) /
                                                            (hypoDMRs_allSets_bins$proportion2 + 0.01) ) )
hypoDMRs_allSets_bins$relative_change <- as.numeric( 1 - ( (hypoDMRs_allSets_bins$proportion2 + 0.01) /
                                                           (hypoDMRs_allSets_bins$proportion1 + 0.01) ) )

# Create new columns containing absolute change, log2 fold change, and relative change percentiles and quantiles
hypoDMRs_allSets_bins$ac_percentile <- as.numeric( rank(hypoDMRs_allSets_bins$absolute_change) /
                                                  length(hypoDMRs_allSets_bins$absolute_change) )
hypoDMRs_allSets_bins$l2fc_percentile <- as.numeric( rank(hypoDMRs_allSets_bins$log2_fold_change) /
                                                    length(hypoDMRs_allSets_bins$log2_fold_change) )
hypoDMRs_allSets_bins$rc_percentile <- as.numeric( rank(hypoDMRs_allSets_bins$relative_change) /
                                                  length(hypoDMRs_allSets_bins$relative_change) )
hypoDMRs_allSets_bins$ac_quantile <- as.character("") 
hypoDMRs_allSets_bins$l2fc_quantile <- as.character("") 
hypoDMRs_allSets_bins$rc_quantile <- as.character("") 

# Define quantiles
ac_quantilesStats <- data.frame()
l2fc_quantilesStats <- data.frame()
rc_quantilesStats <- data.frame()
for(k in 1:args$quantiles) {
  # absolute_change (ac)
  # First quantile should span 1 to greater than, e.g., 0.75 proportions of hypoDMRs_allSets_bins
  if(k < args$quantiles) {
     hypoDMRs_allSets_bins[ !is.na(hypoDMRs_allSets_bins$ac_percentile) &
                            hypoDMRs_allSets_bins$ac_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                            hypoDMRs_allSets_bins$ac_percentile > 1 - ( k / args$quantiles ) ]$ac_quantile <- paste0("Quantile ", k)
  } else { 
  # Final quantile should span 0 to, e.g., 0.25 proportions of hypoDMRs_allSets_bins
    hypoDMRs_allSets_bins[ !is.na(hypoDMRs_allSets_bins$ac_percentile) &
                           hypoDMRs_allSets_bins$ac_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                           hypoDMRs_allSets_bins$ac_percentile >= 1 - ( k / args$quantiles ) ]$ac_quantile <- paste0("Quantile ", k)
  }
  ac_stats <- data.frame(quantile = as.integer(k),
                         n = as.integer(length(hypoDMRs_allSets_bins[hypoDMRs_allSets_bins$ac_quantile == paste0("Quantile ", k)])),
                         mean_width = as.integer(round(mean(
                           width(hypoDMRs_allSets_bins[hypoDMRs_allSets_bins$ac_quantile == paste0("Quantile ", k)]), na.rm = T))),
                         total_width = as.integer(sum(
                           width(hypoDMRs_allSets_bins[hypoDMRs_allSets_bins$ac_quantile == paste0("Quantile ", k)]), na.rm = T)),
                         mean_absolute_change = as.numeric(mean(
                           hypoDMRs_allSets_bins[hypoDMRs_allSets_bins$ac_quantile == paste0("Quantile ", k)]$absolute_change, na.rm = T)),
                         mean_log2_fold_change = as.numeric(mean(
                           hypoDMRs_allSets_bins[hypoDMRs_allSets_bins$ac_quantile == paste0("Quantile ", k)]$log2_fold_change, na.rm = T)),
                         mean_relative_change = as.numeric(mean(
                           hypoDMRs_allSets_bins[hypoDMRs_allSets_bins$ac_quantile == paste0("Quantile ", k)]$relative_change, na.rm = T)),
                         stringsAsFactors = FALSE)
  ac_quantilesStats <- rbind(ac_quantilesStats, ac_stats)
  # log2_fold_change
  # First quantile should span 1 to greater than, e.g., 0.75 proportions of hypoDMRs_allSets_bins
  if(k < args$quantiles) {
     hypoDMRs_allSets_bins[ !is.na(hypoDMRs_allSets_bins$l2fc_percentile) &
                            hypoDMRs_allSets_bins$l2fc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                            hypoDMRs_allSets_bins$l2fc_percentile > 1 - ( k / args$quantiles ) ]$l2fc_quantile <- paste0("Quantile ", k)
  } else { 
  # Final quantile should span 0 to, e.g., 0.25 proportions of hypoDMRs_allSets_bins
    hypoDMRs_allSets_bins[ !is.na(hypoDMRs_allSets_bins$l2fc_percentile) &
                           hypoDMRs_allSets_bins$l2fc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                           hypoDMRs_allSets_bins$l2fc_percentile >= 1 - ( k / args$quantiles ) ]$l2fc_quantile <- paste0("Quantile ", k)
  }
  l2fc_stats <- data.frame(quantile = as.integer(k),
                           n = as.integer(length(hypoDMRs_allSets_bins[hypoDMRs_allSets_bins$l2fc_quantile == paste0("Quantile ", k)])),
                           mean_width = as.integer(round(mean(
                             width(hypoDMRs_allSets_bins[hypoDMRs_allSets_bins$l2fc_quantile == paste0("Quantile ", k)]), na.rm = T))),
                           total_width = as.integer(sum(
                             width(hypoDMRs_allSets_bins[hypoDMRs_allSets_bins$l2fc_quantile == paste0("Quantile ", k)]), na.rm = T)),
                           mean_absolute_change = as.numeric(mean(
                             hypoDMRs_allSets_bins[hypoDMRs_allSets_bins$l2fc_quantile == paste0("Quantile ", k)]$absolute_change, na.rm = T)),
                           mean_log2_fold_change = as.numeric(mean(
                             hypoDMRs_allSets_bins[hypoDMRs_allSets_bins$l2fc_quantile == paste0("Quantile ", k)]$log2_fold_change, na.rm = T)),
                           mean_relative_change = as.numeric(mean(
                             hypoDMRs_allSets_bins[hypoDMRs_allSets_bins$l2fc_quantile == paste0("Quantile ", k)]$relative_change, na.rm = T)),
                           stringsAsFactors = FALSE)
  l2fc_quantilesStats <- rbind(l2fc_quantilesStats, l2fc_stats)
  # relative_change
  # First quantile should span 1 to greater than, e.g., 0.75 proportions of hypoDMRs_allSets_bins
  if(k < args$quantiles) {
     hypoDMRs_allSets_bins[ !is.na(hypoDMRs_allSets_bins$rc_percentile) &
                            hypoDMRs_allSets_bins$rc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                            hypoDMRs_allSets_bins$rc_percentile > 1 - ( k / args$quantiles ) ]$rc_quantile <- paste0("Quantile ", k)
  } else { 
  # Final quantile should span 0 to, e.g., 0.25 proportions of hypoDMRs_allSets_bins
    hypoDMRs_allSets_bins[ !is.na(hypoDMRs_allSets_bins$rc_percentile) &
                           hypoDMRs_allSets_bins$rc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                           hypoDMRs_allSets_bins$rc_percentile >= 1 - ( k / args$quantiles ) ]$rc_quantile <- paste0("Quantile ", k)
  }
  rc_stats <- data.frame(quantile = as.integer(k),
                         n = as.integer(length(hypoDMRs_allSets_bins[hypoDMRs_allSets_bins$rc_quantile == paste0("Quantile ", k)])),
                         mean_width = as.integer(round(mean(
                           width(hypoDMRs_allSets_bins[hypoDMRs_allSets_bins$rc_quantile == paste0("Quantile ", k)]), na.rm = T))),
                         total_width = as.integer(sum(
                           width(hypoDMRs_allSets_bins[hypoDMRs_allSets_bins$rc_quantile == paste0("Quantile ", k)]), na.rm = T)),
                         mean_absolute_change = as.numeric(mean(
                           hypoDMRs_allSets_bins[hypoDMRs_allSets_bins$rc_quantile == paste0("Quantile ", k)]$absolute_change, na.rm = T)),
                         mean_log2_fold_change = as.numeric(mean(
                           hypoDMRs_allSets_bins[hypoDMRs_allSets_bins$rc_quantile == paste0("Quantile ", k)]$log2_fold_change, na.rm = T)),
                         mean_relative_change = as.numeric(mean(
                           hypoDMRs_allSets_bins[hypoDMRs_allSets_bins$rc_quantile == paste0("Quantile ", k)]$relative_change, na.rm = T)),
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
rtracklayer::export(object = hypoDMRs_allSets_bins,
                    con = paste0(hypoDMRdir,
                                 "features_", args$quantiles, "quantiles",
                                 "_by_change_in_",
                                 paste0(args$condition2, collapse = "_"),
                                 "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                                 "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                                 "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_", args$genomeRegion, ".gff3"))
hypoDMRs_allSets_bins_bed <- data.frame(chr = as.character(seqnames(hypoDMRs_allSets_bins)),
                                        start = as.integer(start(hypoDMRs_allSets_bins)-1),
                                        end = as.integer(end(hypoDMRs_allSets_bins)),
                                        name = as.integer(1:length(hypoDMRs_allSets_bins)),
                                        score = as.numeric(hypoDMRs_allSets_bins$log2_fold_change),
                                        strand = as.character(strand(hypoDMRs_allSets_bins)),
                                        stringsAsFactors = FALSE)
write.table(hypoDMRs_allSets_bins_bed,
            file = paste0(hypoDMRdir,
                          "features_", args$quantiles, "quantiles",
                          "_by_change_in_",
                          paste0(args$condition2, collapse = "_"),
                          "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                          "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                          "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_", args$genomeRegion, ".bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


# Define random loci of the same number and width distribution,
# and in the same per-feature genomeRegionGR as hypoDMRs_allSets_bins

# Define function to select randomly positioned loci of the same
# width distribution as hypoDMRs_allSets_bins
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLocGR contains the same number of loci per chromosome as hypoDMRs_allSets_bins
chrs <- args$chrName[which(args$chrName %in% seqnames(hypoDMRs_allSets_bins)@values)]
hypoDMRs_allSets_bins_ranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  hypoDMRs_allSets_binsChrGR <- hypoDMRs_allSets_bins[seqnames(hypoDMRs_allSets_bins) == chrs[i]]
  genomeRegionChrGR <- genomeRegionGR[seqnames(genomeRegionGR) == chrs[i]]
  # Contract genomeRegionChrGR so that random loci and 2-kb flanking regions
  # do not extend beyond chromosome ends
  end(genomeRegionChrGR) <- end(genomeRegionChrGR)-max(width(hypoDMRs_allSets_binsChrGR))-2000
  start(genomeRegionChrGR) <- start(genomeRegionChrGR)+2000
  # Define seed so that random selections are reproducible
  set.seed(93750174)
  hypoDMRs_allSets_bins_ranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(genomeRegionChrGR), function(x) {
                                                                                   start(genomeRegionChrGR[x]) : end(genomeRegionChrGR[x])
                                                                                 })),
                                                            n = length(hypoDMRs_allSets_binsChrGR))
  hypoDMRs_allSets_bins_ranLocChrGR <- GRanges(seqnames = chrs[i],
                                               ranges = IRanges(start = hypoDMRs_allSets_bins_ranLocChrStart,
                                                                width = width(hypoDMRs_allSets_binsChrGR)),
                                               strand = strand(hypoDMRs_allSets_binsChrGR))
  hypoDMRs_allSets_bins_ranLocGR <- append(hypoDMRs_allSets_bins_ranLocGR, hypoDMRs_allSets_bins_ranLocChrGR)
}
stopifnot( length( findOverlaps(query = hypoDMRs_allSets_bins_ranLocGR,
                                subject = genomeMaskGR,
                                type = "any", select = "all",
                                ignore.strand = TRUE) ) == 0 )


# Divide hypoDMRs_allSets_bins_ranLocGR into quantiles based on hypoDMRs_allSets_bins quantile indices
hypoDMRs_allSets_bins_ranLocGR$l2fc_random <- as.character("")
# Get row indices for each feature quantile
l2fc_quantileIndices <- lapply(1:args$quantiles, function(k) {
  which(hypoDMRs_allSets_bins$l2fc_quantile == paste0("Quantile ", k))
})
for(k in 1:args$quantiles) {
  hypoDMRs_allSets_bins_ranLocGR[l2fc_quantileIndices[[k]]]$l2fc_random <- paste0("Random ", k)
}

hypoDMRs_allSets_bins_ranLocGR_bed <- data.frame(chr = as.character(seqnames(hypoDMRs_allSets_bins_ranLocGR)),
                                                 start = as.integer(start(hypoDMRs_allSets_bins_ranLocGR)-1),
                                                 end = as.integer(end(hypoDMRs_allSets_bins_ranLocGR)),
                                                 name = as.integer(1:length(hypoDMRs_allSets_bins_ranLocGR)),
                                                 score = as.character(hypoDMRs_allSets_bins_ranLocGR$l2fc_random),
                                                 strand = as.character(strand(hypoDMRs_allSets_bins_ranLocGR)),
                                                 stringsAsFactors = FALSE)
write.table(hypoDMRs_allSets_bins_ranLocGR_bed,
            file = paste0(hypoDMRdir,
                          "features_", args$quantiles, "quantiles",
                          "_by_change_in_",
                          paste0(args$condition2, collapse = "_"),
                          "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                          "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                          "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_", args$genomeRegion, "_ranLoc.bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# DMRs in CEN
if(length(hypoDMRs_allSets_bins_CEN) > 0) {
  # Create new columns containing absolute change, log2 fold change, and relative change in args$context methylation proportion
  hypoDMRs_allSets_bins_CEN$absolute_change <- as.numeric( hypoDMRs_allSets_bins_CEN$proportion1 - hypoDMRs_allSets_bins_CEN$proportion2 )
  # + 0.01 is an offset to prevent infinite log2 fold change (or equal relative change) values where proportions = 0
  hypoDMRs_allSets_bins_CEN$log2_fold_change <- as.numeric( log2( (hypoDMRs_allSets_bins_CEN$proportion1 + 0.01) /
                                                                  (hypoDMRs_allSets_bins_CEN$proportion2 + 0.01) ) )
  hypoDMRs_allSets_bins_CEN$relative_change <- as.numeric( 1 - ( (hypoDMRs_allSets_bins_CEN$proportion2 + 0.01) /
                                                                 (hypoDMRs_allSets_bins_CEN$proportion1 + 0.01) ) )
  
  # Create new columns containing absolute change, log2 fold change, and relative change percentiles and quantiles
  hypoDMRs_allSets_bins_CEN$ac_percentile <- as.numeric( rank(hypoDMRs_allSets_bins_CEN$absolute_change) /
                                                         length(hypoDMRs_allSets_bins_CEN$absolute_change) )
  hypoDMRs_allSets_bins_CEN$l2fc_percentile <- as.numeric( rank(hypoDMRs_allSets_bins_CEN$log2_fold_change) /
                                                           length(hypoDMRs_allSets_bins_CEN$log2_fold_change) )
  hypoDMRs_allSets_bins_CEN$rc_percentile <- as.numeric( rank(hypoDMRs_allSets_bins_CEN$relative_change) /
                                                         length(hypoDMRs_allSets_bins_CEN$relative_change) )
  hypoDMRs_allSets_bins_CEN$ac_quantile <- as.character("") 
  hypoDMRs_allSets_bins_CEN$l2fc_quantile <- as.character("") 
  hypoDMRs_allSets_bins_CEN$rc_quantile <- as.character("") 

  # Define quantiles
  ac_quantilesStats <- data.frame()
  l2fc_quantilesStats <- data.frame()
  rc_quantilesStats <- data.frame()
  for(k in 1:args$quantiles) {
    # absolute_change (ac)
    # First quantile should span 1 to greater than, e.g., 0.75 proportions of hypoDMRs_allSets_bins_CEN
    if(k < args$quantiles) {
       hypoDMRs_allSets_bins_CEN[ !is.na(hypoDMRs_allSets_bins_CEN$ac_percentile) &
                                  hypoDMRs_allSets_bins_CEN$ac_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                  hypoDMRs_allSets_bins_CEN$ac_percentile > 1 - ( k / args$quantiles ) ]$ac_quantile <- paste0("Quantile ", k)
    } else { 
    # Final quantile should span 0 to, e.g., 0.25 proportions of hypoDMRs_allSets_bins_CEN
      hypoDMRs_allSets_bins_CEN[ !is.na(hypoDMRs_allSets_bins_CEN$ac_percentile) &
                                 hypoDMRs_allSets_bins_CEN$ac_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                 hypoDMRs_allSets_bins_CEN$ac_percentile >= 1 - ( k / args$quantiles ) ]$ac_quantile <- paste0("Quantile ", k)
    }
    ac_stats <- data.frame(quantile = as.integer(k),
                           n = as.integer(length(hypoDMRs_allSets_bins_CEN[hypoDMRs_allSets_bins_CEN$ac_quantile == paste0("Quantile ", k)])),
                           mean_width = as.integer(round(mean(
                             width(hypoDMRs_allSets_bins_CEN[hypoDMRs_allSets_bins_CEN$ac_quantile == paste0("Quantile ", k)]), na.rm = T))),
                           total_width = as.integer(sum(
                             width(hypoDMRs_allSets_bins_CEN[hypoDMRs_allSets_bins_CEN$ac_quantile == paste0("Quantile ", k)]), na.rm = T)),
                           mean_absolute_change = as.numeric(mean(
                             hypoDMRs_allSets_bins_CEN[hypoDMRs_allSets_bins_CEN$ac_quantile == paste0("Quantile ", k)]$absolute_change, na.rm = T)),
                           mean_log2_fold_change = as.numeric(mean(
                             hypoDMRs_allSets_bins_CEN[hypoDMRs_allSets_bins_CEN$ac_quantile == paste0("Quantile ", k)]$log2_fold_change, na.rm = T)),
                           mean_relative_change = as.numeric(mean(
                             hypoDMRs_allSets_bins_CEN[hypoDMRs_allSets_bins_CEN$ac_quantile == paste0("Quantile ", k)]$relative_change, na.rm = T)),
                           stringsAsFactors = FALSE)
    ac_quantilesStats <- rbind(ac_quantilesStats, ac_stats)
    # log2_fold_change
    # First quantile should span 1 to greater than, e.g., 0.75 proportions of hypoDMRs_allSets_bins_CEN
    if(k < args$quantiles) {
       hypoDMRs_allSets_bins_CEN[ !is.na(hypoDMRs_allSets_bins_CEN$l2fc_percentile) &
                                  hypoDMRs_allSets_bins_CEN$l2fc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                  hypoDMRs_allSets_bins_CEN$l2fc_percentile > 1 - ( k / args$quantiles ) ]$l2fc_quantile <- paste0("Quantile ", k)
    } else { 
    # Final quantile should span 0 to, e.g., 0.25 proportions of hypoDMRs_allSets_bins_CEN
      hypoDMRs_allSets_bins_CEN[ !is.na(hypoDMRs_allSets_bins_CEN$l2fc_percentile) &
                                 hypoDMRs_allSets_bins_CEN$l2fc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                 hypoDMRs_allSets_bins_CEN$l2fc_percentile >= 1 - ( k / args$quantiles ) ]$l2fc_quantile <- paste0("Quantile ", k)
    }
    l2fc_stats <- data.frame(quantile = as.integer(k),
                             n = as.integer(length(hypoDMRs_allSets_bins_CEN[hypoDMRs_allSets_bins_CEN$l2fc_quantile == paste0("Quantile ", k)])),
                             mean_width = as.integer(round(mean(
                               width(hypoDMRs_allSets_bins_CEN[hypoDMRs_allSets_bins_CEN$l2fc_quantile == paste0("Quantile ", k)]), na.rm = T))),
                             total_width = as.integer(sum(
                               width(hypoDMRs_allSets_bins_CEN[hypoDMRs_allSets_bins_CEN$l2fc_quantile == paste0("Quantile ", k)]), na.rm = T)),
                             mean_absolute_change = as.numeric(mean(
                               hypoDMRs_allSets_bins_CEN[hypoDMRs_allSets_bins_CEN$l2fc_quantile == paste0("Quantile ", k)]$absolute_change, na.rm = T)),
                             mean_log2_fold_change = as.numeric(mean(
                               hypoDMRs_allSets_bins_CEN[hypoDMRs_allSets_bins_CEN$l2fc_quantile == paste0("Quantile ", k)]$log2_fold_change, na.rm = T)),
                             mean_relative_change = as.numeric(mean(
                               hypoDMRs_allSets_bins_CEN[hypoDMRs_allSets_bins_CEN$l2fc_quantile == paste0("Quantile ", k)]$relative_change, na.rm = T)),
                             stringsAsFactors = FALSE)
    l2fc_quantilesStats <- rbind(l2fc_quantilesStats, l2fc_stats)
    # relative_change
    # First quantile should span 1 to greater than, e.g., 0.75 proportions of hypoDMRs_allSets_bins_CEN
    if(k < args$quantiles) {
       hypoDMRs_allSets_bins_CEN[ !is.na(hypoDMRs_allSets_bins_CEN$rc_percentile) &
                                  hypoDMRs_allSets_bins_CEN$rc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                  hypoDMRs_allSets_bins_CEN$rc_percentile > 1 - ( k / args$quantiles ) ]$rc_quantile <- paste0("Quantile ", k)
    } else { 
    # Final quantile should span 0 to, e.g., 0.25 proportions of hypoDMRs_allSets_bins_CEN
      hypoDMRs_allSets_bins_CEN[ !is.na(hypoDMRs_allSets_bins_CEN$rc_percentile) &
                                 hypoDMRs_allSets_bins_CEN$rc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                 hypoDMRs_allSets_bins_CEN$rc_percentile >= 1 - ( k / args$quantiles ) ]$rc_quantile <- paste0("Quantile ", k)
    }
    rc_stats <- data.frame(quantile = as.integer(k),
                           n = as.integer(length(hypoDMRs_allSets_bins_CEN[hypoDMRs_allSets_bins_CEN$rc_quantile == paste0("Quantile ", k)])),
                           mean_width = as.integer(round(mean(
                             width(hypoDMRs_allSets_bins_CEN[hypoDMRs_allSets_bins_CEN$rc_quantile == paste0("Quantile ", k)]), na.rm = T))),
                           total_width = as.integer(sum(
                             width(hypoDMRs_allSets_bins_CEN[hypoDMRs_allSets_bins_CEN$rc_quantile == paste0("Quantile ", k)]), na.rm = T)),
                           mean_absolute_change = as.numeric(mean(
                             hypoDMRs_allSets_bins_CEN[hypoDMRs_allSets_bins_CEN$rc_quantile == paste0("Quantile ", k)]$absolute_change, na.rm = T)),
                           mean_log2_fold_change = as.numeric(mean(
                             hypoDMRs_allSets_bins_CEN[hypoDMRs_allSets_bins_CEN$rc_quantile == paste0("Quantile ", k)]$log2_fold_change, na.rm = T)),
                           mean_relative_change = as.numeric(mean(
                             hypoDMRs_allSets_bins_CEN[hypoDMRs_allSets_bins_CEN$rc_quantile == paste0("Quantile ", k)]$relative_change, na.rm = T)),
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
  rtracklayer::export(object = hypoDMRs_allSets_bins_CEN,
                      con = paste0(hypoDMRdir,
                                   "features_", args$quantiles, "quantiles",
                                   "_by_change_in_",
                                   paste0(args$condition2, collapse = "_"),
                                   "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                                   "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                                   "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_CEN.gff3"))
  hypoDMRs_allSets_bins_CEN_bed <- data.frame(chr = as.character(seqnames(hypoDMRs_allSets_bins_CEN)),
                                              start = as.integer(start(hypoDMRs_allSets_bins_CEN)-1),
                                              end = as.integer(end(hypoDMRs_allSets_bins_CEN)),
                                              name = as.integer(1:length(hypoDMRs_allSets_bins_CEN)),
                                              score = as.numeric(hypoDMRs_allSets_bins_CEN$log2_fold_change),
                                              strand = as.character(strand(hypoDMRs_allSets_bins_CEN)),
                                              stringsAsFactors = FALSE)
  write.table(hypoDMRs_allSets_bins_CEN_bed,
              file = paste0(hypoDMRdir,
                            "features_", args$quantiles, "quantiles",
                            "_by_change_in_",
                            paste0(args$condition2, collapse = "_"),
                            "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                            "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                            "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_CEN.bed"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


  # Define random loci of the same number and width distribution,
  # and in the same per-feature CENGR as hypoDMRs_allSets_bins_CEN
  
  # Define function to select randomly positioned loci of the same
  # width distribution as hypoDMRs_allSets_bins_CEN
  ranLocStartSelect <- function(coordinates, n) {
    sample(x = coordinates,
           size = n,
           replace = FALSE)
  }
  
  # Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
  options(scipen = 100)
  
  # Apply ranLocStartSelect() on a per-chromosome basis so that
  # ranLocGR contains the same number of loci per chromosome as hypoDMRs_allSets_bins_CEN
  chrs <- args$chrName[which(args$chrName %in% seqnames(hypoDMRs_allSets_bins_CEN)@values)]
  hypoDMRs_allSets_bins_CEN_ranLocGR <- GRanges()
  for(i in 1:length(chrs)) {
    hypoDMRs_allSets_bins_CENChrGR <- hypoDMRs_allSets_bins_CEN[seqnames(hypoDMRs_allSets_bins_CEN) == chrs[i]]
    CENChrGR <- CENGR[seqnames(CENGR) == chrs[i]]
    # Contract CENChrGR so that random loci and 2-kb flanking regions
    # do not extend beyond chromosome ends
    end(CENChrGR) <- end(CENChrGR)-max(width(hypoDMRs_allSets_bins_CENChrGR))-2000
    start(CENChrGR) <- start(CENChrGR)+2000
    # Define seed so that random selections are reproducible
    set.seed(93750174)
    hypoDMRs_allSets_bins_CEN_ranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(CENChrGR), function(x) {
                                                                                         start(CENChrGR[x]) : end(CENChrGR[x])
                                                                                       })),
                                                                  n = length(hypoDMRs_allSets_bins_CENChrGR))
    hypoDMRs_allSets_bins_CEN_ranLocChrGR <- GRanges(seqnames = chrs[i],
                                                     ranges = IRanges(start = hypoDMRs_allSets_bins_CEN_ranLocChrStart,
                                                                      width = width(hypoDMRs_allSets_bins_CENChrGR)),
                                                     strand = strand(hypoDMRs_allSets_bins_CENChrGR))
    hypoDMRs_allSets_bins_CEN_ranLocGR <- append(hypoDMRs_allSets_bins_CEN_ranLocGR, hypoDMRs_allSets_bins_CEN_ranLocChrGR)
  }
  #stopifnot( length( findOverlaps(query = hypoDMRs_allSets_bins_CEN_ranLocGR,
  #                                subject = nonCENGR,
  #                                type = "any", select = "all",
  #                                ignore.strand = TRUE) ) == 0 )


  # Divide hypoDMRs_allSets_bins_CEN_ranLocGR into quantiles based on hypoDMRs_allSets_bins_CEN quantile indices
  hypoDMRs_allSets_bins_CEN_ranLocGR$l2fc_random <- as.character("")
  # Get row indices for each feature quantile
  l2fc_quantileIndices <- lapply(1:args$quantiles, function(k) {
    which(hypoDMRs_allSets_bins_CEN$l2fc_quantile == paste0("Quantile ", k))
  })
  for(k in 1:args$quantiles) {
    hypoDMRs_allSets_bins_CEN_ranLocGR[l2fc_quantileIndices[[k]]]$l2fc_random <- paste0("Random ", k)
  }

  hypoDMRs_allSets_bins_CEN_ranLocGR_bed <- data.frame(chr = as.character(seqnames(hypoDMRs_allSets_bins_CEN_ranLocGR)),
                                                       start = as.integer(start(hypoDMRs_allSets_bins_CEN_ranLocGR)-1),
                                                       end = as.integer(end(hypoDMRs_allSets_bins_CEN_ranLocGR)),
                                                       name = as.integer(1:length(hypoDMRs_allSets_bins_CEN_ranLocGR)),
                                                       score = as.character(hypoDMRs_allSets_bins_CEN_ranLocGR$l2fc_random),
                                                       strand = as.character(strand(hypoDMRs_allSets_bins_CEN_ranLocGR)),
                                                       stringsAsFactors = FALSE)
  write.table(hypoDMRs_allSets_bins_CEN_ranLocGR_bed,
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
if(length(hypoDMRs_allSets_bins_nonCEN) > 0) {
  # Create new columns containing absolute change, log2 fold change, and relative change in args$context methylation proportion
  hypoDMRs_allSets_bins_nonCEN$absolute_change <- as.numeric( hypoDMRs_allSets_bins_nonCEN$proportion1 - hypoDMRs_allSets_bins_nonCEN$proportion2 )
  # + 0.01 is an offset to prevent infinite log2 fold change (or equal relative change) values where proportions = 0
  hypoDMRs_allSets_bins_nonCEN$log2_fold_change <- as.numeric( log2( (hypoDMRs_allSets_bins_nonCEN$proportion1 + 0.01) /
                                                                     (hypoDMRs_allSets_bins_nonCEN$proportion2 + 0.01) ) )
  hypoDMRs_allSets_bins_nonCEN$relative_change <- as.numeric( 1 - ( (hypoDMRs_allSets_bins_nonCEN$proportion2 + 0.01) /
                                                                    (hypoDMRs_allSets_bins_nonCEN$proportion1 + 0.01) ) )
  
  # Create new columns containing absolute change, log2 fold change, and relative change percentiles and quantiles
  hypoDMRs_allSets_bins_nonCEN$ac_percentile <- as.numeric( rank(hypoDMRs_allSets_bins_nonCEN$absolute_change) /
                                                            length(hypoDMRs_allSets_bins_nonCEN$absolute_change) )
  hypoDMRs_allSets_bins_nonCEN$l2fc_percentile <- as.numeric( rank(hypoDMRs_allSets_bins_nonCEN$log2_fold_change) /
                                                              length(hypoDMRs_allSets_bins_nonCEN$log2_fold_change) )
  hypoDMRs_allSets_bins_nonCEN$rc_percentile <- as.numeric( rank(hypoDMRs_allSets_bins_nonCEN$relative_change) /
                                                            length(hypoDMRs_allSets_bins_nonCEN$relative_change) )
  hypoDMRs_allSets_bins_nonCEN$ac_quantile <- as.character("") 
  hypoDMRs_allSets_bins_nonCEN$l2fc_quantile <- as.character("") 
  hypoDMRs_allSets_bins_nonCEN$rc_quantile <- as.character("") 

  # Define quantiles
  ac_quantilesStats <- data.frame()
  l2fc_quantilesStats <- data.frame()
  rc_quantilesStats <- data.frame()
  for(k in 1:args$quantiles) {
    # absolute_change (ac)
    # First quantile should span 1 to greater than, e.g., 0.75 proportions of hypoDMRs_allSets_bins_nonCEN
    if(k < args$quantiles) {
       hypoDMRs_allSets_bins_nonCEN[ !is.na(hypoDMRs_allSets_bins_nonCEN$ac_percentile) &
                                     hypoDMRs_allSets_bins_nonCEN$ac_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                     hypoDMRs_allSets_bins_nonCEN$ac_percentile > 1 - ( k / args$quantiles ) ]$ac_quantile <- paste0("Quantile ", k)
    } else { 
    # Final quantile should span 0 to, e.g., 0.25 proportions of hypoDMRs_allSets_bins_nonCEN
      hypoDMRs_allSets_bins_nonCEN[ !is.na(hypoDMRs_allSets_bins_nonCEN$ac_percentile) &
                                    hypoDMRs_allSets_bins_nonCEN$ac_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                    hypoDMRs_allSets_bins_nonCEN$ac_percentile >= 1 - ( k / args$quantiles ) ]$ac_quantile <- paste0("Quantile ", k)
    }
    ac_stats <- data.frame(quantile = as.integer(k),
                           n = as.integer(length(hypoDMRs_allSets_bins_nonCEN[hypoDMRs_allSets_bins_nonCEN$ac_quantile == paste0("Quantile ", k)])),
                           mean_width = as.integer(round(mean(
                             width(hypoDMRs_allSets_bins_nonCEN[hypoDMRs_allSets_bins_nonCEN$ac_quantile == paste0("Quantile ", k)]), na.rm = T))),
                           total_width = as.integer(sum(
                             width(hypoDMRs_allSets_bins_nonCEN[hypoDMRs_allSets_bins_nonCEN$ac_quantile == paste0("Quantile ", k)]), na.rm = T)),
                           mean_absolute_change = as.numeric(mean(
                             hypoDMRs_allSets_bins_nonCEN[hypoDMRs_allSets_bins_nonCEN$ac_quantile == paste0("Quantile ", k)]$absolute_change, na.rm = T)),
                           mean_log2_fold_change = as.numeric(mean(
                             hypoDMRs_allSets_bins_nonCEN[hypoDMRs_allSets_bins_nonCEN$ac_quantile == paste0("Quantile ", k)]$log2_fold_change, na.rm = T)),
                           mean_relative_change = as.numeric(mean(
                             hypoDMRs_allSets_bins_nonCEN[hypoDMRs_allSets_bins_nonCEN$ac_quantile == paste0("Quantile ", k)]$relative_change, na.rm = T)),
                           stringsAsFactors = FALSE)
    ac_quantilesStats <- rbind(ac_quantilesStats, ac_stats)
    # log2_fold_change
    # First quantile should span 1 to greater than, e.g., 0.75 proportions of hypoDMRs_allSets_bins_nonCEN
    if(k < args$quantiles) {
       hypoDMRs_allSets_bins_nonCEN[ !is.na(hypoDMRs_allSets_bins_nonCEN$l2fc_percentile) &
                                     hypoDMRs_allSets_bins_nonCEN$l2fc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                     hypoDMRs_allSets_bins_nonCEN$l2fc_percentile > 1 - ( k / args$quantiles ) ]$l2fc_quantile <- paste0("Quantile ", k)
    } else { 
    # Final quantile should span 0 to, e.g., 0.25 proportions of hypoDMRs_allSets_bins_nonCEN
      hypoDMRs_allSets_bins_nonCEN[ !is.na(hypoDMRs_allSets_bins_nonCEN$l2fc_percentile) &
                                    hypoDMRs_allSets_bins_nonCEN$l2fc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                    hypoDMRs_allSets_bins_nonCEN$l2fc_percentile >= 1 - ( k / args$quantiles ) ]$l2fc_quantile <- paste0("Quantile ", k)
    }
    l2fc_stats <- data.frame(quantile = as.integer(k),
                             n = as.integer(length(hypoDMRs_allSets_bins_nonCEN[hypoDMRs_allSets_bins_nonCEN$l2fc_quantile == paste0("Quantile ", k)])),
                             mean_width = as.integer(round(mean(
                               width(hypoDMRs_allSets_bins_nonCEN[hypoDMRs_allSets_bins_nonCEN$l2fc_quantile == paste0("Quantile ", k)]), na.rm = T))),
                             total_width = as.integer(sum(
                               width(hypoDMRs_allSets_bins_nonCEN[hypoDMRs_allSets_bins_nonCEN$l2fc_quantile == paste0("Quantile ", k)]), na.rm = T)),
                             mean_absolute_change = as.numeric(mean(
                               hypoDMRs_allSets_bins_nonCEN[hypoDMRs_allSets_bins_nonCEN$l2fc_quantile == paste0("Quantile ", k)]$absolute_change, na.rm = T)),
                             mean_log2_fold_change = as.numeric(mean(
                               hypoDMRs_allSets_bins_nonCEN[hypoDMRs_allSets_bins_nonCEN$l2fc_quantile == paste0("Quantile ", k)]$log2_fold_change, na.rm = T)),
                             mean_relative_change = as.numeric(mean(
                               hypoDMRs_allSets_bins_nonCEN[hypoDMRs_allSets_bins_nonCEN$l2fc_quantile == paste0("Quantile ", k)]$relative_change, na.rm = T)),
                             stringsAsFactors = FALSE)
    l2fc_quantilesStats <- rbind(l2fc_quantilesStats, l2fc_stats)
    # relative_change
    # First quantile should span 1 to greater than, e.g., 0.75 proportions of hypoDMRs_allSets_bins_nonCEN
    if(k < args$quantiles) {
       hypoDMRs_allSets_bins_nonCEN[ !is.na(hypoDMRs_allSets_bins_nonCEN$rc_percentile) &
                                     hypoDMRs_allSets_bins_nonCEN$rc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                     hypoDMRs_allSets_bins_nonCEN$rc_percentile > 1 - ( k / args$quantiles ) ]$rc_quantile <- paste0("Quantile ", k)
    } else { 
    # Final quantile should span 0 to, e.g., 0.25 proportions of hypoDMRs_allSets_bins_nonCEN
      hypoDMRs_allSets_bins_nonCEN[ !is.na(hypoDMRs_allSets_bins_nonCEN$rc_percentile) &
                                    hypoDMRs_allSets_bins_nonCEN$rc_percentile <= 1 - ( (k - 1) / args$quantiles ) &
                                    hypoDMRs_allSets_bins_nonCEN$rc_percentile >= 1 - ( k / args$quantiles ) ]$rc_quantile <- paste0("Quantile ", k)
    }
    rc_stats <- data.frame(quantile = as.integer(k),
                           n = as.integer(length(hypoDMRs_allSets_bins_nonCEN[hypoDMRs_allSets_bins_nonCEN$rc_quantile == paste0("Quantile ", k)])),
                           mean_width = as.integer(round(mean(
                             width(hypoDMRs_allSets_bins_nonCEN[hypoDMRs_allSets_bins_nonCEN$rc_quantile == paste0("Quantile ", k)]), na.rm = T))),
                           total_width = as.integer(sum(
                             width(hypoDMRs_allSets_bins_nonCEN[hypoDMRs_allSets_bins_nonCEN$rc_quantile == paste0("Quantile ", k)]), na.rm = T)),
                           mean_absolute_change = as.numeric(mean(
                             hypoDMRs_allSets_bins_nonCEN[hypoDMRs_allSets_bins_nonCEN$rc_quantile == paste0("Quantile ", k)]$absolute_change, na.rm = T)),
                           mean_log2_fold_change = as.numeric(mean(
                             hypoDMRs_allSets_bins_nonCEN[hypoDMRs_allSets_bins_nonCEN$rc_quantile == paste0("Quantile ", k)]$log2_fold_change, na.rm = T)),
                           mean_relative_change = as.numeric(mean(
                             hypoDMRs_allSets_bins_nonCEN[hypoDMRs_allSets_bins_nonCEN$rc_quantile == paste0("Quantile ", k)]$relative_change, na.rm = T)),
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
  rtracklayer::export(object = hypoDMRs_allSets_bins_nonCEN,
                      con = paste0(hypoDMRdir,
                                   "features_", args$quantiles, "quantiles",
                                   "_by_change_in_",
                                   paste0(args$condition2, collapse = "_"),
                                   "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                                   "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                                   "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_nonCEN.gff3"))
  hypoDMRs_allSets_bins_nonCEN_bed <- data.frame(chr = as.character(seqnames(hypoDMRs_allSets_bins_nonCEN)),
                                                 start = as.integer(start(hypoDMRs_allSets_bins_nonCEN)-1),
                                                 end = as.integer(end(hypoDMRs_allSets_bins_nonCEN)),
                                                 name = as.integer(1:length(hypoDMRs_allSets_bins_nonCEN)),
                                                 score = as.numeric(hypoDMRs_allSets_bins_nonCEN$log2_fold_change),
                                                 strand = as.character(strand(hypoDMRs_allSets_bins_nonCEN)),
                                                 stringsAsFactors = FALSE)
  write.table(hypoDMRs_allSets_bins_nonCEN_bed,
              file = paste0(hypoDMRdir,
                            "features_", args$quantiles, "quantiles",
                            "_by_change_in_",
                            paste0(args$condition2, collapse = "_"),
                            "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                            "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                            "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_nonCEN.bed"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


  # Define random loci of the same number and width distribution,
  # and in the same per-feature nonCENGR as hypoDMRs_allSets_bins_nonCEN
  
  # Define function to select randomly positioned loci of the same
  # width distribution as hypoDMRs_allSets_bins_nonCEN
  ranLocStartSelect <- function(coordinates, n) {
    sample(x = coordinates,
           size = n,
           replace = FALSE)
  }
  
  # Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
  options(scipen = 100)
  
  # Apply ranLocStartSelect() on a per-chromosome basis so that
  # ranLocGR contains the same number of loci per chromosome as hypoDMRs_allSets_bins_nonCEN
  chrs <- args$chrName[which(args$chrName %in% seqnames(hypoDMRs_allSets_bins_nonCEN)@values)]
  hypoDMRs_allSets_bins_nonCEN_ranLocGR <- GRanges()
  for(i in 1:length(chrs)) {
    hypoDMRs_allSets_bins_nonCENChrGR <- hypoDMRs_allSets_bins_nonCEN[seqnames(hypoDMRs_allSets_bins_nonCEN) == chrs[i]]
    nonCENChrGR <- nonCENGR[seqnames(nonCENGR) == chrs[i]]
    # Contract nonCENChrGR so that random loci and 2-kb flanking regions
    # do not extend beyond chromosome ends
    end(nonCENChrGR) <- end(nonCENChrGR)-max(width(hypoDMRs_allSets_bins_nonCENChrGR))-2000
    start(nonCENChrGR) <- start(nonCENChrGR)+2000
    # Define seed so that random selections are reproducible
    set.seed(93750174)
    hypoDMRs_allSets_bins_nonCEN_ranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(nonCENChrGR), function(x) {
                                                                                         start(nonCENChrGR[x]) : end(nonCENChrGR[x])
                                                                                       })),
                                                                  n = length(hypoDMRs_allSets_bins_nonCENChrGR))
    hypoDMRs_allSets_bins_nonCEN_ranLocChrGR <- GRanges(seqnames = chrs[i],
                                                     ranges = IRanges(start = hypoDMRs_allSets_bins_nonCEN_ranLocChrStart,
                                                                      width = width(hypoDMRs_allSets_bins_nonCENChrGR)),
                                                     strand = strand(hypoDMRs_allSets_bins_nonCENChrGR))
    hypoDMRs_allSets_bins_nonCEN_ranLocGR <- append(hypoDMRs_allSets_bins_nonCEN_ranLocGR, hypoDMRs_allSets_bins_nonCEN_ranLocChrGR)
  }
  #stopifnot( length( findOverlaps(query = hypoDMRs_allSets_bins_nonCEN_ranLocGR,
  #                                subject = CENGR,
  #                                type = "any", select = "all",
  #                                ignore.strand = TRUE) ) == 0 )


  # Divide hypoDMRs_allSets_bins_nonCEN_ranLocGR into quantiles based on hypoDMRs_allSets_bins_nonCEN quantile indices
  hypoDMRs_allSets_bins_nonCEN_ranLocGR$l2fc_random <- as.character("")
  # Get row indices for each feature quantile
  l2fc_quantileIndices <- lapply(1:args$quantiles, function(k) {
    which(hypoDMRs_allSets_bins_nonCEN$l2fc_quantile == paste0("Quantile ", k))
  })
  for(k in 1:args$quantiles) {
    hypoDMRs_allSets_bins_nonCEN_ranLocGR[l2fc_quantileIndices[[k]]]$l2fc_random <- paste0("Random ", k)
  }

  hypoDMRs_allSets_bins_nonCEN_ranLocGR_bed <- data.frame(chr = as.character(seqnames(hypoDMRs_allSets_bins_nonCEN_ranLocGR)),
                                                          start = as.integer(start(hypoDMRs_allSets_bins_nonCEN_ranLocGR)-1),
                                                          end = as.integer(end(hypoDMRs_allSets_bins_nonCEN_ranLocGR)),
                                                          name = as.integer(1:length(hypoDMRs_allSets_bins_nonCEN_ranLocGR)),
                                                          score = as.character(hypoDMRs_allSets_bins_nonCEN_ranLocGR$l2fc_random),
                                                          strand = as.character(strand(hypoDMRs_allSets_bins_nonCEN_ranLocGR)),
                                                          stringsAsFactors = FALSE)
  write.table(hypoDMRs_allSets_bins_nonCEN_ranLocGR_bed,
              file = paste0(hypoDMRdir,
                            "features_", args$quantiles, "quantiles",
                            "_by_change_in_",
                            paste0(args$condition2, collapse = "_"),
                            "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                            "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                            "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_nonCEN_ranLoc.bed"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}
