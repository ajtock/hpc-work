#!/home/ajt200/miniconda3/envs/R-4.1.2/bin/Rscript

# Estimate methylation divergence and epimutation rates and spectra with AlphaBeta,
# using output TXT files from METHimpute run on Bismark-processed
# BS-seq data from mutation accumulation (MA) lines

# Target outputs are genomic binned methylation divergence values,
# to enable analysis of relationships with binned among-read variation in
# ONT DeepSignal-derived DNA methylation (Fleiss' kappa), and
# definition of epimutation hotspots and coldspots corresponding to
# the top 10% and bottom 10% of bins with regard to methylation divergence

# Usage:
# ./alphabeta_per_cytosine_MA1_2.R t2t-col.20210610 CpG 10000 1000

args <- commandArgs(trailingOnly = T)
refbase <- args[1]
context <- args[2]
genomeBinSize <- as.numeric(args[3]
genomeStepSize <- as.numeric(args[4])

#refbase <- "t2t-col.20210610"
#context <- "CpG"
#genomeBinSize <- 10000
#genomeStepSize <- 1000

options(stringsAsFactors = F)
library(AlphaBeta)
library(dplyr)
library(data.table)
library(parallel)
library(yaml)
config <- read_yaml("../config.yaml")

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

inDir <- paste0("../coverage/report/methimpute/")
inDirBin <- paste0(inDir, "genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName, "/")
outDir <- paste0("../coverage/report/alphabeta/genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName, "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", inDirBin, " ] || mkdir -p ", inDirBin))
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

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
filePathsGlobal <- paste0(inDir, config$SAMPLES, "_MappedOn_", refbase, "_dedup_", context, "_methylome.txt")

methylomesGlobalList <- mclapply(1:length(filePathsGlobal), function(x) {
#methylomesGlobalList <- mclapply(1, function(x) {
  fread(filePathsGlobal[x])
#}, mc.cores = length(filePathsGlobal[[1]]), mc.preschedule = F)
}, mc.cores = length(filePathsGlobal), mc.preschedule = F)

# Define genomic windows
print(genomeBinName)
print(genomeStepName)
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

targetDF <- data.frame()
for(i in 1:nrow(binDF)) {
  bin_i <- binDF[i,]
  filePaths_bin_i_trunc <- gsub(pattern = "methylome.txt", replacement = paste0("methylome_", paste0(bin_i, collapse = "_"), ".txt"),
                                x = gsub(pattern = "../coverage/report/methimpute/", replacement = "",
                                         x = filePathsGlobal))
  filePaths_bin_i <- paste0(inDirBin, filePaths_bin_i_trunc)
  mclapply(1:length(methylomesGlobalList), function(x) {

    methylome_bin_i <- methylomesGlobalList[[x]] %>%
      dplyr::filter(seqnames == bin_i$chr) %>%
      dplyr::filter(start >= bin_i$start & start <= bin_i$end)

    fwrite(methylome_bin_i,
           file = filePaths_bin_i[x],
           quote = F, sep = "\t", row.names = F, col.names = T)

  }, mc.cores = length(methylomesGlobalList), mc.preschedule = F)

 



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
node_file <- paste0(outDir, "nodelist_MA1_2_MappedOn_", refbase, "_", context, ".fn")
write.table(node_df,
            file = node_file,
            quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
#fwrite(node_df,
#       file = node_file,
#       quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
#node_df <- fread(node_file)
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

edge_file <- paste0(outDir, "edgelist_MA1_2_MappedOn_", refbase, "_", context, ".fn")
write.table(edge_df,
            file = edge_file,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
#fwrite(edge_df,
#       file = edge_file,
#       quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
#edge_df <- fread(edge_file)
#print(edge_df)


# Build the pedigree of the MA lines

output_file <- paste0(outDir, "pedigree_output_MA1_2_MappedOn_", refbase, "_", context, ".RData")
output <- buildPedigree(nodelist = node_file,
                        edgelist = edge_file,
                        cytosine = sub("p", "", context),
                        posteriorMaxFilter = 0.99)
outputTmp <- output
save(output,
     file = output_file)
rm(output)
load(output_file)
stopifnot(identical(output, outputTmp))
rm(outputTmp); gc()

# Plot the pedigree of the MA lines

# plotPedigree doesn't work with sampling.design set to "progenitor.endpoint":
#plotPedigree(nodelist = node_file,
#             edgelist = edge_file,
#             sampling.design = "progenitor.endpoint",
#             output.dir = plotDir,
#             plot.width = 5, plot.height = 5, aspect.ratio = 1,
#             vertex.size = 6, vertex.label = FALSE,
#             out.pdf = paste0("pedigree_output_MA1_2_MappedOn_", refbase, "_", context))
#sampleFile <- "../nodelist_MA1_1.fn"
#edgesFile <- "../edgelist_MA1_1.fn"
#sampleFile <- system.file("extdata/vg", "nodelist.fn", package = "AlphaBeta")
#edgesFile <- system.file("extdata/vg", "edgelist.fn", package = "AlphaBeta")
#plotPedigree(nodelist = sampleFile, edgelist = edgesFile, sampling.design = "progenitor.endpoint", output.dir = plotDir, plot.width = 5, plot.height = 5, aspect.ratio = 1, vertex.size = 6, vertex.label = FALSE, out.pdf = "MA1_1")
## Error: Can't convert from <tbl_df<
##   V1          : character
##   V2          : character
##   weight.scale: character
##   grp         : character
##   relevel     : double
## >> to <tbl_df<
##   V1          : character
##   V2          : character
##   weight.scale: character
##   grp         : character
## >> due to loss of precision.
## Dropped variables: `relevel`


# Plot divergence time (delta.t) versus methylome divergence (D.value)
# Interactive plot for inspecting the divergence data and removing outlier samples (if any)

pedigree <- output$Pdata
dt <- pedigree[, 2] + pedigree[, 3] - 2 * pedigree[, 1]
dt_max <- max(dt, na.rm = T)

pdf(paste0(plotDir, "divergence_time_vs_methylome_divergence_MA1_2_MappedOn_", refbase, "_", context, ".pdf"))
plot(x = dt,
     y = pedigree[, "D.value"],
     ylab = "Methylome divergence value",
     xlab = expression(paste(Delta, italic("t"), sep = "")))
dev.off()


# Epimutation rate estimation in selfing systems

# "Models ABneutral, ABselectMM, ABselectUU and ABnull can be used to estimate
# the rate of spontaneous epimutations in selfing-derived MA lines.
# The models are currently restricted to diploids."

# Run model with no selection (ABneutral)
# "ABneutral fits a neutral epimutation model. The model assumes that epimutation
# accumulation is under no selective constraint. Returned are estimates of the
# methylation gain and loss rates and the proportion of epi-heterozygote loci
# in the pedigree founder genome."

# "NOTE: it is recommended to use at least 50 Nstarts to achieve best solutions"

print("Initial proportions of unmethylated cytosines:")
p0uu_in <- output$tmpp0
print(p0uu_in)

ABneutral_out <- ABneutral(pedigree.data = pedigree,
                           p0uu = p0uu_in, eqp = p0uu_in, eqp.weight = 1,
                           Nstarts = 50, out.dir = outDir,
                           out.name = paste0("ABneutral_global_estimates_MA1_2_MappedOn_", refbase, "_", context))
print(summary(ABneutral_out))
head(ABneutral_out$pedigree)
ABneutral_file <- paste0(outDir, "ABneutral_global_estimates_MA1_2_MappedOn_", refbase, "_", context, ".Rdata")
ABneutral_outTest <- dget(ABneutral_file)
#stopifnot(identical(ABneutral_out, ABneutral_outTest))
rm(ABneutral_outTest); gc()

# Plot estimates of ABneutral model
# "In 'ABplot' function you can set parameters to customize the pdf output."
ABplot(pedigree.names = ABneutral_file,
       output.dir = plotDir,
       out.name = paste0("ABneutral_global_estimates_MA1_2_MappedOn_", refbase, "_", context),
       plot.height = 8, plot.width = 11)


# Boostrap analysis

# "NOTE: it is recommended to use at least 50 Nboot to achieve best solutions" 

ABneutral_BOOTout <- BOOTmodel(pedigree.data = ABneutral_file,
                               Nboot = 50,
                               out.dir = outDir,
                               out.name = paste0("ABneutral_BOOT_global_estimates_MA1_2_MappedOn_", refbase, "_", context))
print(summary(ABneutral_BOOTout))
ABneutralBOOT_file <- paste0(outDir, "ABneutral_BOOT_global_estimates_MA1_2_MappedOn_", refbase, "_", context, ".Rdata")
ABneutral_BOOToutTest <- dget(ABneutralBOOT_file)
#stopifnot(identical(ABneutral_BOOTout, ABneutral_BOOToutTest))
rm(ABneutral_BOOToutTest); gc()

