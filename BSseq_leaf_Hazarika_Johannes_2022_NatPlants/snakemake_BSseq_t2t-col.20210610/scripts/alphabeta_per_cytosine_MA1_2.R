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
config <- read_yaml("../config.yaml")

inDir <- paste0("../coverage/report/methimpute/")
outDir <- paste0("../coverage/report/alphabeta/")
#plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
#system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# Genomic definitions
fai <- read.table(paste0("../data/index/", refbase, ".fa.fai"), header = F)
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
node <- gsub("../coverage/report/methimpute/MA\\d+_\\d+_G", "", filePaths)
node <- gsub("_SRR.+", "", node)
node <- gsub("_BSseq", "", node)
node <- gsub("_L", "_", node)
generation <- as.integer(gsub("_.+", "", node))
methylome <- rep("Y", length(node))

# Make pedigree files for the MA lineage, to enable alphabeta to
# reconstruct the topology of the underlying pedigree

# "AlphaBeta requires two types of input files: 'nodeslist.fn'
# and 'edgelist.fn'. The structure of these files follows the
# standard file format required by the R network package igraph."

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
print(node_df)


# Create edges file, corresponding to a sparse directed acyclic graph (DAG)
# of connections between individuals (nodes) from different generations
# e.g., "edgelist_MA1_2_MappedOn_t2t-col.20210610_dedup_CpG.fn"
# "from and to: Specifies the network edges, which are any direct connections
# between type S and type S* nodes in the pedigree. Only unique pairs of nodes
# need to be supplied. These 2 columns are mandatory."
# "gendiff (optional): Specifies the number of generations that separate the
# two nodes. This column is useful only for plotting purposes and it can be
# omitted for epimutation rate estimation. However, we recommened that this
# column be supplied because it is useful for accurately scaling the edge lengths
# when plotting certain pedigrees with progenitor.endpoint and sibling design"
# "group (optional): Along with "gendiff" column, groupings supplied in this column
# will help in scaling the edge lengths when plotting the pedigree."

edge_df_G0 <- data.frame(from = node_df$node[which(node_df$meth == "N")][1],
                         to = node_df$node[which(node_df$meth == "N")][-1])

edge_df_GN <- data.frame()
for(x in node_df$node[which(node_df$meth == "N")][-1]) {
  edge_df_x <- data.frame(from = x,
                          to = node_df$node[which(node_df$meth == "Y")] [
                                 grep( 
                                      paste0( paste0(as.integer(gsub("_\\d+", "", gsub("_Rep\\d+", "", x))) + 1),
                                              gsub("^\\d+", "", gsub("_Rep\\d+", "", x)),
                                              "_" ),
                                      node_df$node[which(node_df$meth == "Y")]
                                     )
                               ]
                         )
  edge_df_GN <- rbind(edge_df_GN, edge_df_x)
}

edge_df <- rbind(edge_df_G0, edge_df_GN)
#edge_df$line <- as.numeric( paste0(
##                                   gsub("_", "",
##                                        gsub("_Rep.+", "",
##                                             edge_df$from) ),
#                                   gsub("_", ".",
#                                        gsub("_Rep", "",
#                                             edge_df$to) )
#                                  ) )
edge_df$gendiff <- as.integer(gsub("_.+", "", edge_df$to)) -
                   as.integer(gsub("_.+", "", edge_df$from))
edge_df$gendiff[edge_df$gendiff == 2] <- 1
edge_df$gento <- as.integer(gsub("_.+", "", edge_df$to))
edge_df$group <- "A"
edge_df[edge_df$gento >= 30,]$group <- "B"
edge_df <- edge_df[ with(edge_df, order(to)), ]

dropCols <- c("gento")
edge_df <- edge_df[ , !(names(edge_df) %in% dropCols)]
print(edge_df)

fwrite(edge_df,
       file = paste0(outDir, "edgelist_MA1_2_MappedOn_", refbase, "_", context, ".fn"),
       quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
edge_df <- fread(paste0(outDir, "edgelist_MA1_2_MappedOn_", refbase, "_", context, ".fn"))
print(edge_df)



