#!/usr/bin/env Rscript

# Use METHimpute to call context-specific methylation status at each cytosine using a HMM-based binomial test.
# "Besides improved accuracy over the classical binomial test, the HMM allows imputation of the methylation
# status of all cytosines in the genome. It achieves this by borrowing information from neighboring covered
# cytosines. The confidence in the methylation status call is reported as well."

# The output TXT file from METHimpute can be used as an input file for AlphaBeta to calculate

# Usage:
# ./methimpute_per_cytosine.R MA1_2_G3_L1_BSseq_Rep1_SRR342347 t2t-col.20210610 CpG

args <- commandsArgs(trailingOnly = T)
libName <- args[1]
refbase <- args[2]
context <- args[3]

#libName <- "MA1_2_G3_L1_BSseq_Rep1_SRR342347"
#refbase <- "t2t-col.20210610"
#context <- "CpG"

options(stringsAsFactors = F)
library(methimpute)
library(yaml)
config <- read_yaml("config.yaml")

outDir <- "coverage/report/methimpute/"
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# Genomic definitions
fai <- read.table(paste0("data/index/", refbase, ".fa.fai"), header = F)
#ignoreChrs <- unlist(strsplit(config$GENOMEPROFILES$ignoreChrs,
#                              split = " "))
#fai <- fai[!(fai$V1 %in% ignoreChrs),]
#if(!grepl("Chr", fai[,1][1])) {
#  chrs <- paste0("Chr", fai[,1])
#} else {
#  chrs <- fai[,1]
#}
#chrLens <- fai[,2]
chromosomes <- fai[,1:2]
colnames(chromosomes) <- c("chromosome", "length")


#### Step 1: Import the data

# Define path to input file for methimpute
filePath <- paste0("coverage/report/",
                   libName, "_MappedOn_", refbase, "_dedup_", context, ".CX_report.txt.gz")
# Import the data
bm <- importBismark(filePath, chrom.lengths = chromosomes)
print(bm)

# "Because most methylation extractor programs report only covered cytosines,
# we need to inflate the data to include all cytosines (including non-covered sites)"
fastaPath <- paste0("data/index/", refbase, ".fa")
cytosinePos <- extractCytosinesFromFASTA(fastaPath,
                                         contexts = c("CG", "CHG", "CHH"))
methylome <- inflateMethylome(bm, cytosinePos)
print(methylome)


#### Step 2: Obtain correlation parameters

# "The correlation of methylation levels between neighboring cytosines is an important
# parameter for the methylation status calling, so we need to get it first"
distcor <- distanceCorrelation(methylome, separate.contexts = TRUE)
fit <- estimateTransDist(distcor)
print(fit$transDist)
ggsave(file = paste0(plotDir, libName, "_MappedOn_", refbase, "_dedup_", context, "_distanceCorrelation.pdf"),
       plot = fit$plot,
       height = 2.5*3, width = 3.5*3, limitsize = FALSE)


#### Step 3: Call and impute methylation status 

model <- callMethylationSeparate(data = methylome, transDist = fit$transDist,
                                 verbosity = 0)
# "The confidence in the methylation status call is given in the column "posteriorMax".
# For further analysis one could split the results into high-confidence
# (posteriorMax >= 0.98) and low-confidence calls (posteriorMax < 0.98) for instance."
print(model)




