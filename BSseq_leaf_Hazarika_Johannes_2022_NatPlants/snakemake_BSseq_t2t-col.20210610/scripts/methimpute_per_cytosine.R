#!/usr/bin/env Rscript

# Use METHimpute to call context-specific methylation status at each cytosine using a HMM-based binomial test.
# "Besides improved accuracy over the classical binomial test, the HMM allows imputation of the methylation
# status of all cytosines in the genome. It achieves this by borrowing information from neighboring covered
# cytosines. The confidence in the methylation status call is reported as well."

# The output TXT file from METHimpute can be used as an input file for AlphaBeta to calculate

# Usage:
# ./methimpute_per_cytosine.R MA1_2_G3_L1_BSseq_Rep1_SRR342347 t2t-col.20210610 CpG

args <- commandArgs(trailingOnly = T)
libName <- args[1]
refbase <- args[2]
context <- args[3]

#libName <- "MA1_2_G3_L1_BSseq_Rep1_SRR342347"
#refbase <- "t2t-col.20210610"
#context <- "CpG"

options(stringsAsFactors = F)
library(methimpute)
library(dplyr)
library(data.table)
library(stringr)
library(yaml)
config <- read_yaml("config.yaml")

outDir <- paste0("coverage/report/methimpute/")
#plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
#system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

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
methylome <- methylome[seqnames(methylome) %in% chrs]
print(methylome)


#### Step 2: Obtain correlation parameters

## Interacting-context model

# "The interacting-context model runs a single HMM for all contexts. This takes into account
# the within-context and between-context correlations and should be more accurate than the
# separate-context model if sufficient data is available. However, we have observed that in low
# coverage settings too much information from well covered contexts is diffusing into the low
#covered contexts (e.g. CHH and CHG will look like CG with very low coverage). In this case,
# please use the separate-context model"

# "The correlation of methylation levels between neighboring cytosines is an important
# parameter for the methylation status calling, so we need to get it first"
## NOTE "separate":
distCor <- distanceCorrelation(methylome,
                               separate.contexts = FALSE)
fit <- estimateTransDist(distCor)
print(fit$transDist)
plotDir <- paste0(outDir, "plots/distanceCorrelation/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))
ggsave(file = paste0(plotDir, libName, "_MappedOn_", refbase, "_dedup_", context, "_distanceCorrelation.pdf"),
       plot = fit$plot,
       height = 2.5*3, width = 3.5*3, limitsize = FALSE)


#### Step 3: Call and impute methylation status 

## NOTE "separate":
model <- callMethylation(data = methylome,
                         transDist = fit$transDist,
                         include.intermediate = TRUE,
                         update = "constrained",
                         verbosity = 1)
# "The confidence in the methylation status call is given in the column "posteriorMax".
# For further analysis one could split the results into high-confidence
# (posteriorMax >= 0.98) and low-confidence calls (posteriorMax < 0.98) for instance."
print(model)

# "Bisulfite conversion rates can be obtained with"
1 - model$params$emissionParams$Unmethylated

# Retain cytosines with a maximum posterior probability of >= 0.99
# to avoid low-quality methylation status calls
model$data <- model$data[model$data$posteriorMax >= 0.99]

# "You can also check several properties of the fitted Hidden Markov Model, such as convergence
# or transition probabilities, and check how well the fitted distributions describe the data."
plotDir <- paste0(outDir, "plots/convergence/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))
ggsave(file = paste0(plotDir, libName, "_MappedOn_", refbase, "_dedup_", context, "_convergence.pdf"),
       plot = plotConvergence(model),
       height = 2.5, width = 3.5*3, limitsize = FALSE)
plotDir <- paste0(outDir, "plots/transitionProbs/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))
ggsave(file = paste0(plotDir, libName, "_MappedOn_", refbase, "_dedup_", context, "_transitionProbs.pdf"),
       plot = plotTransitionProbs(model),
       height = 2.5, width = 3.5*3, limitsize = FALSE)
plotDir <- paste0(outDir, "plots/fittedDists/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))
ggsave(file = paste0(plotDir, libName, "_MappedOn_", refbase, "_dedup_", context, "_fittedDists.pdf"),
       plot = plotHistogram(model, total.counts = 10),
       height = 2.5, width = 3.5*3, limitsize = FALSE)


# Adapt model data for use by alphabeta (Bioconductor package for epimutation rate calculation)
model_data <- model$data
model_df <- methods::as(model_data, "data.frame")
model_df <- model_df[,c("seqnames", "start", "strand", "context", "counts.methylated", "counts.total",
                        "posteriorMax", "posteriorMeth", "posteriorUnmeth", "status", "rc.meth.lvl")]
model_df <- model_df %>%
  dplyr::mutate_if(is.factor, as.character)
#final_dataset <- model_df

cx <- fread(filePath, skip = 0, sep = "\t",
            select = c(1, 2, 7),
            col.names = c("seqnames", "start", "context.trinucleotide"),
            stringsAsFactors = F)
# Append cx colmns
final_dataset <- model_df %>%
  dplyr::left_join(cx, by = c("seqnames", "start"))

# Drop columns not needed by alphabeta
dropCols <- c("posteriorMeth", "posteriorUnmeth")
final_dataset <- final_dataset[ , !(names(final_dataset) %in% dropCols)]

# Abbreviate methylation status names for use with alphabeta
final_dataset$status <- stringr::str_replace_all(final_dataset$status, pattern = "Unmethylated", replacement = "U")
final_dataset$status <- stringr::str_replace_all(final_dataset$status, pattern = "Intermediate", replacement = "I")
final_dataset$status <- stringr::str_replace_all(final_dataset$status, pattern = "Methylated", replacement = "M")
 
# Take 4 digits of decimal values in "posteriorMax" and "rc.meth.lvl" 
floor_dec <- function(x, level = 1) round(x - 5*10^(-level-1), level)
final_dataset$posteriorMax <- floor_dec(as.numeric(as.character(final_dataset$posteriorMax)), 4)
final_dataset$rc.meth.lvl <- floor_dec(as.numeric(as.character(final_dataset$rc.meth.lvl)), 4)
final_dataset <- final_dataset[final_dataset$context == gsub("p", "", context),]
print(paste0("Writing to file: ", outDir, "methylome_", libName, "_MappedOn_", refbase, "_dedup_", context, ".txt"))
fwrite(final_dataset,
       file = paste0(outDir, libName, "_MappedOn_", refbase, "_dedup_", context, "_methylome.txt"),
       quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
rm(bm, cytosinePos, methylome, distCor, model, model_data, model_df, cx, final_dataset); gc()
