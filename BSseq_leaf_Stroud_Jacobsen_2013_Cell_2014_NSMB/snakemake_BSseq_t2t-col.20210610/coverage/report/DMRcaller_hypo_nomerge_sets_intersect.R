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

# Sort by seqnames, start and end
hypoDMRs_allSets_bins <- sortSeqlevels(hypoDMRs_allSets_bins)
hypoDMRs_allSets_bins <- sort(hypoDMRs_allSets_bins, ignore.strand = TRUE)
hypoDMRs_set1Only_bins <- sortSeqlevels(hypoDMRs_set1Only_bins)
hypoDMRs_set1Only_bins <- sort(hypoDMRs_set1Only_bins, ignore.strand = TRUE)

# Export hypoDMRs_allSets_bins DMR GRanges as annotation files
rtracklayer::export(object = hypoDMRs_allSets_bins,
                    con = paste0(hypoDMRdir,
                                 "features_", args$quantiles, "quantiles_",
                                 paste0(args$sets, collapse = "_and_"),
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
                          "features_", args$quantiles, "quantiles_",
                          paste0(args$sets, collapse = "_and_"),
                          "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                          "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                          "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_", args$genomeRegion, ".bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Export hypoDMRs_set1Only_bins DMR GRanges as annotation files
rtracklayer::export(object = hypoDMRs_set1Only_bins,
                    con = paste0(hypoDMRdir,
                                 "features_", args$quantiles, "quantiles_",
                                 paste0(args$sets, collapse = "_not_"),
                                 "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                                 "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                                 "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_", args$genomeRegion, ".gff3"))
hypoDMRs_set1Only_bins_bed <- data.frame(chr = as.character(seqnames(hypoDMRs_set1Only_bins)),
                                         start = as.integer(start(hypoDMRs_set1Only_bins)-1),
                                         end = as.integer(end(hypoDMRs_set1Only_bins)),
                                         name = as.integer(1:length(hypoDMRs_set1Only_bins)),
                                         score = as.numeric(hypoDMRs_set1Only_bins$log2_fold_change),
                                         strand = as.character(strand(hypoDMRs_set1Only_bins)),
                                         stringsAsFactors = FALSE)
write.table(hypoDMRs_set1Only_bins_bed,
            file = paste0(hypoDMRdir,
                          "features_", args$quantiles, "quantiles_",
                          paste0(args$sets, collapse = "_not_"),
                          "_hypo", sub("p", "", args$context), "_DMRs_vs3reps",
                          "_mbins_bS100_tfisher_pVT0.01_mCC4_mRPC4_mPD_0.4_0.2_0.1_mG200",
                          "_in_", args$refbase, "_", paste0(args$chrName, collapse = "_"), "_", args$genomeRegion, ".bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
