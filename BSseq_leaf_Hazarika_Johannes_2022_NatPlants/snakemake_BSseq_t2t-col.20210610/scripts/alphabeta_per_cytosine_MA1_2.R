#!/usr/bin/env Rscript

# Estimate epimutation rates and spectra with AlphaBeta,
# using output TXT files from METHimpute run on Bismark-processed
# BS-seq data from mutation accumulation (MA) lines

# Usage:
# ./alphabeta_per_cytosine_MA1_2.R t2t-col.20210610 CpG

args <- commandArgs(trailingOnly = T)
refbase <- args[1]
context <- args[2]

#refbase <- "t2t-col.20210610"
#context <- "CpG"

options(stringsAsFactors = F)
library(AlphaBeta)
library(dplyr)
library(data.table)
#library(stringr)
library(yaml)
config <- read_yaml("config.yaml")

inDir <- paste0("coverage/report/methimpute/")
outDir <- paste0("coverage/report/alphabeta/")
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


# Define paths to methylome TXT files generated with methimpute
filePaths <- paste0(inDir, config$SAMPLES, "_MappedOn_", refbase, "_dedup_", context, "_methylome.txt")
# Extract node, generation and methylome info from filePaths
node <- gsub("coverage/report/methimpute/MA\\d+_\\d+_G", "", filePaths)
node <- gsub("_SRR.+", "", node)
node <- gsub("_BSseq", "", node)
node <- gsub("_L", "_", node)
generation <- as.integer(gsub("_.+", "", node))
methylome <- rep("Y", length(node))

# Make pedigree files for the MA lineage, to enable alphabeta to
# reconstruct the topology of the underlying pedigree

# "The nodes of the network correspond to 'individuals' whose
# methylomes have been sampled (i.e. type S* nodes), or of the
# common ancestors of these individuals, whose methylomes have
# typically not been sampled (i.e. type S nodes)"
# e.g., "nodelist_MA1_2_MappedOn_t2t-col.20210610_dedup_CpG.fn"
node_df_GM <- data.frame(filename = filePaths,
                         node = node,
                         gen = generation,
                         meth = methylome)
node_df_G0 <- data.frame(filename = "-",
                         node = "0_0",
                         gen  = as.integer(0),
                         meth = "N")
node_df_G2 <- data.frame(filename = "-",
                         node = gsub("^3_", "2_",
                                     node[ grep( "Rep1", node ) ] [ grep( "^3_", node[ grep( "Rep1", node ) ] ) ] ),
                         gen = as.integer(2),
                         meth = "N") 
node_df_G3 <- node_df_GM[node_df_GM$gen == 3,]
node_df_G30 <- data.frame(filename = "-",
                          node = gsub("^31_", "30_",
                                      node[ grep( "Rep1", node ) ] [ grep( "^31_", node[ grep( "Rep1", node ) ] ) ] ),
                          gen = as.integer(30),
                          meth = "N") 
node_df_G31 <- node_df_GM[node_df_GM$gen == 31,]

rm(node_df_GM)

# Combine node generations into one data.frame
node_df <- dplyr::bind_rows(mget(sort(ls(pattern = "node_df"))))
print(node_df)
fwrite(node_df,
       file = paste0(outDir, "nodelist_MA1_2_MappedOn_", refbase, "_", context, ".fn"),
       quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
node_df <- fread(paste0(outDir, "nodelist_MA1_2_MappedOn_", refbase, "_", context, ".fn"))

node_df$node[which(node_df$meth == "N")[1]]








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

## Retain cytosines with a maximum posterior probability of >= 0.99
## to avoid low-quality methylation status calls
#model$data <- model$data[model$data$posteriorMax >= 0.99]

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
