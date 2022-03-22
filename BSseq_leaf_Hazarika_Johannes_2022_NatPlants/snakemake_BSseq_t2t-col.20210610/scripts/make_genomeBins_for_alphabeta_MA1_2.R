#!/home/ajt200/miniconda3/envs/R-4.1.2/bin/Rscript

# Make genomic bins for later estimation of methylation divergence with AlphaBeta,
# using output TXT files from METHimpute run on Bismark-processed
# BS-seq data from mutation accumulation (MA) lines

# Usage:
# ./make_genomeBins_for_alphabeta_MA1_2.R t2t-col.20210610 CpG 10000 1000

args <- commandArgs(trailingOnly = T)
refbase <- args[1]
context <- args[2]
genomeBinSize <- as.numeric(args[3])
genomeStepSize <- as.numeric(args[4])

#refbase <- "t2t-col.20210610"
#context <- "CpG"
#genomeBinSize <- 10000
#genomeStepSize <- 1000

options(stringsAsFactors = F)
options(scipen=999)
library(dplyr)
library(data.table)
library(yaml)
config <- read_yaml("config.yaml")

if(floor(log10(genomeBinSize)) + 1 < 4) {
  genomeBinName <- paste0(genomeBinSize, "bp")
  genomeBinNamePlot <- paste0(genomeBinSize, "-bp")
} else if(floor(log10(genomeBinSize)) + 1 >= 4 &
          floor(log10(genomeBinSize)) + 1 <= 6) {
  genomeBinName <- paste0(genomeBinSize/1e3, "kb")
  genomeBinNamePlot <- paste0(genomeBinSize/1e3, "-kb")
} else if(floor(log10(genomeBinSize)) + 1 >= 7) {
  genomeBinName <- paste0(genomeBinSize/1e6, "Mb")
  genomeBinNamePlot <- paste0(genomeBinSize/1e6, "-Mb")
}

if(floor(log10(genomeStepSize)) + 1 < 4) {
  genomeStepName <- paste0(genomeStepSize, "bp")
  genomeStepNamePlot <- paste0(genomeStepSize, "-bp")
} else if(floor(log10(genomeStepSize)) + 1 >= 4 &
          floor(log10(genomeStepSize)) + 1 <= 6) {
  genomeStepName <- paste0(genomeStepSize/1e3, "kb")
  genomeStepNamePlot <- paste0(genomeStepSize/1e3, "-kb")
} else if(floor(log10(genomeStepSize)) + 1 >= 7) {
  genomeStepName <- paste0(genomeStepSize/1e6, "Mb")
  genomeStepNamePlot <- paste0(genomeStepSize/1e6, "-Mb")
}

inDir <- paste0("coverage/report/methimpute/")
inDirBin <- paste0(inDir, "genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName, "/")
outDir <- paste0("coverage/report/alphabeta/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", inDirBin, " ] || mkdir -p ", inDirBin))
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# Genomic definitions
fai <- read.table(paste0("data/index/", refbase, ".fa.fai"), header = F)
chromosomes <- fai[,1:2]
colnames(chromosomes) <- c("chromosome", "length")
ignoreChrs <- unlist(strsplit(config$GENOMEPROFILES$ignoreChrs,
                              split = " "))
fai <- fai[!(fai$V1 %in% ignoreChrs),]
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1])
} else {
  chrs <- fai[,1]
}
chrLens <- fai[,2]

# Define paths to methylome TXT files generated with methimpute
filePathsGlobal <- paste0(inDir, config$SAMPLES, "_MappedOn_", refbase, "_dedup_", context, "_methylome.txt")

binDF_file <- paste0(outDir, 
# Define genomic windows
binDF <- data.frame()
for(i in seq_along(chrs)) {
  # Define sliding windows of width genomeBinSize bp,
  # with a step of genomeStepSize bp
  ## Note: the active code creates windows of genomeBinSize,
  ## and the last window of a chromosome may be smaller than genomeBinSize
  binStarts <- seq(from = 1,
                   to = chrLens[i]-genomeBinSize,
                   by = genomeStepSize)
  if(chrLens[i] - binStarts[length(binStarts)] >= genomeBinSize) {
    binStarts <- c(binStarts,
                   binStarts[length(binStarts)]+genomeStepSize)
  }
  binEnds <- seq(from = binStarts[1]+genomeBinSize-1,
                 to = chrLens[i],
                 by = genomeStepSize)
  binEnds <- c(binEnds,
               rep(chrLens[i], times = length(binStarts)-length(binEnds)))                                                                    
  stopifnot(binEnds[length(binEnds)] == chrLens[i])
  stopifnot(length(binStarts) == length(binEnds))

  chr_binDF <- data.frame(chr = chrs[i],
                          start = binStarts,
                          end = binEnds)
  binDF <- rbind(binDF, chr_binDF)
}
fwrite(binDF,
       file = paste0(outDir, "mD_genomeBins_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
                     "_MA1_2_MappedOn_", refbase, "_", context, ".tsv"),
       quote = F, sep = "\t", row.names = F, col.names = T)
