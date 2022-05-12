#!/usr/bin/env Rscript

# This R script is called by ../Snakefile
# Calculate mean DNA methylation proportions in adjacent windows
# for generating chromosome-scale plots

# Usage:
# /applications/R/R-4.0.0/bin/Rscript ./DNAmeth_proportion_perbin.R WT_BSseq_Rep1_2014 t2t-col.20210610 200000 5

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
genomeBinSize <- as.integer(args[3])
# smoothN must be an odd number
# see profile smoothing below
smoothN <- as.numeric(args[4])

if(floor(log10(genomeBinSize)) + 1 < 4) {
  genomeBinName <- paste0(genomeBinSize, "bp")
} else if(floor(log10(genomeBinSize)) + 1 >= 4 &
          floor(log10(genomeBinSize)) + 1 <= 6) {
  genomeBinName <- paste0(genomeBinSize/1e3, "kb")
} else if(floor(log10(genomeBinSize)) + 1 >= 7) {
  genomeBinName <- paste0(genomeBinSize/1e6, "Mb")
}

options(stringsAsFactors = F)
library(segmentSeq)
library(yaml)
config <- read_yaml("config.yaml")

# Load gzip-compressed DNA methylation percentage tables
# and convert percentages into proportions
CG <- read.table(gzfile(paste0("coverage/", sampleName,
                               "_MappedOn_", refbase, "_dedup_CpG.gz.bismark.cov.gz")),
                 colClasses = c(NA, NA, "NULL", NA, "NULL", "NULL"))
CG$V4 <- CG$V4
CHG <- read.table(gzfile(paste0("coverage/", sampleName,
                                "_MappedOn_", refbase, "_dedup_CHG.gz.bismark.cov.gz")),
                  colClasses = c(NA, NA, "NULL", NA, "NULL", "NULL"))
CHG$V4 <- CHG$V4
CHH <- read.table(gzfile(paste0("coverage/", sampleName,
                                "_MappedOn_", refbase, "_dedup_CHH.gz.bismark.cov.gz")),
                  colClasses = c(NA, NA, "NULL", NA, "NULL", "NULL"))
CHH$V4 <- CHH$V4

# Genomic definitions
fai <- read.table(paste0("data/index/", refbase, ".fa.fai"), header = F)
ignoreChrs <- unlist(strsplit(config$GENOMEPROFILES$ignoreChrs,
                              split = " "))
fai <- fai[!(fai$V1 %in% ignoreChrs),]
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1])
} else {
  chrs <- fai[,1]
}
chrLens <- fai[,2]

# Make chromosomal coordinates cumulative
# such that the first coordinate of Chr2 is
# equal to the last coordinate of Chr1 + 1
sumchr <- cumsum(c(0, chrLens))
print(sumchr)

# For each genomeBinSize-bp adjacent window, calculate the mean of the genomeBinSize
# per-base DNA methylation proportions within that window
print(genomeBinName)
cumMethDat <- NULL
for(i in seq_along(chrs)) {
  # Define adjacent windows
  winSeq <- seq(from = 1, to = chrLens[i], by = genomeBinSize)
  winCum <- winSeq + sumchr[i]
  winIR <- IRanges(start = winSeq,
                   width = genomeBinSize)
  winIR <- winIR[-length(winIR)]
  winIR <- append(winIR,
                  IRanges(start = winSeq[length(winSeq)],
                          end = chrLens[i]))
  winGR <- GRanges(seqnames = chrs[i],
                   ranges = winIR,
                   strand = "*")
  print(winGR)

  # Define DNA methylation coordinates as GRanges objects
  # and calculate mean methylation proportions in each window
  # CG
  chrCG <- CG[CG[,1] == chrs[i],]
  chrCG_GR <- GRanges(seqnames = chrs[i],
                      ranges = IRanges(start = chrCG[,2],
                                       width = 1),
                      strand = "*")
  CGoverlaps <- getOverlaps(coordinates = winGR,
                            segments = chrCG_GR,
                            overlapType = "overlapping",
                            whichOverlaps = TRUE,
                            ignoreStrand = TRUE)
  CGwinVals <- sapply(CGoverlaps, function(x) mean(as.numeric(chrCG[,3][x])))
  # CHG
  chrCHG <- CHG[CHG[,1] == chrs[i],]
  chrCHG_GR <- GRanges(seqnames = chrs[i],
                       ranges = IRanges(start = chrCHG[,2],
                                        width = 1),
                       strand = "*")
  CHGoverlaps <- getOverlaps(coordinates = winGR,
                             segments = chrCHG_GR,
                             overlapType = "overlapping",
                             whichOverlaps = TRUE,
                             ignoreStrand = TRUE)
  CHGwinVals <- sapply(CHGoverlaps, function(x) mean(as.numeric(chrCHG[,3][x])))
  # CHH
  chrCHH <- CHH[CHH[,1] == chrs[i],]
  chrCHH_GR <- GRanges(seqnames = chrs[i],
                       ranges = IRanges(start = chrCHH[,2],
                                        width = 1),
                       strand = "*")
  CHHoverlaps <- getOverlaps(coordinates = winGR,
                             segments = chrCHH_GR,
                             overlapType = "overlapping",
                             whichOverlaps = TRUE,
                             ignoreStrand = TRUE)
  CHHwinVals <- sapply(CHHoverlaps, function(x) mean(as.numeric(chrCHH[,3][x])))
  # Mean of all 3 contexts
  CwinVals <- sapply(seq_along(CGwinVals), function(x) {
                      mean(c(CGwinVals[x], CHGwinVals[x], CHHwinVals[x]))
                    })
  # Combine in data frame
  methDat <- data.frame(chr = as.character(chrs[i]),
                        window = as.integer(start(winGR)),
                        cumwindow = as.integer(start(winGR) + sumchr[i]),
                        mCG  = as.numeric(CGwinVals),
                        mCHG = as.numeric(CHGwinVals),
                        mCHH = as.numeric(CHHwinVals),
                        mC   = as.numeric(CwinVals))
  cumMethDat <- rbind(cumMethDat, methDat)
}
write.table(cumMethDat,
            file = paste0("coverage/tsv/",
                          "DNAmeth_", sampleName, "_MappedOn_", refbase,
                          "_dedup_binSize", genomeBinName, "_unsmoothed.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)

chrMethDat <- lapply(seq_along(chrs), function(x) {
  cumMethDat[cumMethDat$chr == chrs[x],]
})

# Calculate moving average (MA) of current window,
# (smoothN/2)-0.5 previous windows (where smoothN is odd),
# and
# (smoothN/2)-0.5 subsequent windows (where smoothN is odd)
# The higher smoothN is, the greater the smoothing
# Use modulo operation (%%) to confirm that the remainder of
# integer division smoothN %/% 2 is 1 (i.e., smoothN is odd)
stopifnot(smoothN %% 2 == 1)
flank <- smoothN %/% 2
# Define MA filter coefficients
f <- rep(x = 1/smoothN, times = smoothN)

filt_methDat <- NULL
for(i in seq_along(chrs)) {
  # mCG
  filt_chr_mCG <- stats::filter(chrMethDat[[i]]$mCG,
                                filter = f,
                                sides = 2)
  filt_chr_mCG[1:flank] <- filt_chr_mCG[flank+1]
  filt_chr_mCG[(length(filt_chr_mCG)-flank+1):length(filt_chr_mCG)] <- filt_chr_mCG[(length(filt_chr_mCG)-flank)]
  # mCHG
  filt_chr_mCHG <- stats::filter(chrMethDat[[i]]$mCHG,
                                 filter = f,
                                 sides = 2)
  filt_chr_mCHG[1:flank] <- filt_chr_mCHG[flank+1]
  filt_chr_mCHG[(length(filt_chr_mCHG)-flank+1):length(filt_chr_mCHG)] <- filt_chr_mCHG[(length(filt_chr_mCHG)-flank)]
  # mCHH
  filt_chr_mCHH <- stats::filter(chrMethDat[[i]]$mCHH,
                                 filter = f,
                                 sides = 2)
  filt_chr_mCHH[1:flank] <- filt_chr_mCHH[flank+1]
  filt_chr_mCHH[(length(filt_chr_mCHH)-flank+1):length(filt_chr_mCHH)] <- filt_chr_mCHH[(length(filt_chr_mCHH)-flank)]
  # mC
  filt_chr_mC <- stats::filter(chrMethDat[[i]]$mC,
                               filter = f,
                               sides = 2)
  filt_chr_mC[1:flank] <- filt_chr_mC[flank+1]
  filt_chr_mC[(length(filt_chr_mC)-flank+1):length(filt_chr_mC)] <- filt_chr_mC[(length(filt_chr_mC)-flank)]
  # Combine in data frame
  filt_chrMethDat <- data.frame(chr = chrMethDat[[i]]$chr,
                                window = as.integer(chrMethDat[[i]]$window),
                                cumwindow = as.integer(chrMethDat[[i]]$cumwindow),
                                filt_mCG  = as.numeric(filt_chr_mCG),
                                filt_mCHG = as.numeric(filt_chr_mCHG),
                                filt_mCHH = as.numeric(filt_chr_mCHH),
                                filt_mC   = as.numeric(filt_chr_mC))
  filt_methDat <- rbind(filt_methDat, filt_chrMethDat)
}
write.table(filt_methDat,
            file = paste0("coverage/tsv/",
                          "DNAmeth_", sampleName, "_MappedOn_", refbase,
                          "_dedup_binSize", genomeBinName, "_smoothed.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)
