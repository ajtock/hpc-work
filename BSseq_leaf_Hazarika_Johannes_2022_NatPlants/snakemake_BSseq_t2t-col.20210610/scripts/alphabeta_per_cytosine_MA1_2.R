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
# ./alphabeta_per_cytosine_MA1_2.R t2t-col.20210610 CpG 10000 1000 Chr1

args <- commandArgs(trailingOnly = T)
refbase <- args[1]
context <- args[2]
genomeBinSize <- as.numeric(args[3])
genomeStepSize <- as.numeric(args[4])
chrName <- args[5]

#refbase <- "t2t-col.20210610"
#context <- "CpG"
#genomeBinSize <- 10000
#genomeStepSize <- 1000
#chrName <- "Chr1"

options(stringsAsFactors = F)
options(scipen=999)
library(AlphaBeta)
library(dplyr)
library(data.table)
library(parallel)
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
  chrs <- paste0("Chr", fai[,1])[which(paste0("Chr", fai[,1]) %in% chrName)]
} else {
  chrs <- fai[,1][which(fai[,1] %in% chrName)]
}
chrLens <- fai[,2][which(fai[,1] %in% chrName)]

# Define paths to methylome TXT files generated with methimpute
filePathsGlobal <- paste0(inDir, config$SAMPLES, "_MappedOn_", refbase, "_dedup_", context, "_methylome.txt")

# Load methylome TXT files for subsetting cytosines (rows) by genomic bin
# NOTE: might be worth changing methimpute script to filter out cytosines
# that are below a given coverage threshold;
# e.g., >= 4; see https://www.pnas.org/doi/10.1073/pnas.1424254112
# However, coverage threshold in https://www.nature.com/articles/s41477-021-01086-7
# seems to be >= 1 (as implemented by default in MethylStar) and with a maximum posterior probability >= 0.99 
methylomesGlobalList <- mclapply(1:length(filePathsGlobal), function(x) {
  fread(filePathsGlobal[x])
}, mc.cores = length(filePathsGlobal), mc.preschedule = F)

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


targetDF_list <- lapply(1:nrow(binDF), function(i) {

  print(i)
  bin_i <- binDF[i,]
  filePaths_bin_i_trunc <- gsub(pattern = "methylome.txt", replacement = paste0("methylome_", paste0(bin_i, collapse = "_"), ".txt"),
                                x = gsub(pattern = "coverage/report/methimpute/", replacement = "",
                                         x = filePathsGlobal))
  filePaths_bin_i <- paste0(inDirBin, filePaths_bin_i_trunc)

  # Write filePaths_bin_i
  # NOTE: remove these files after use by buildPedigree() due to large file numbers (> 1M)
  mclapply(1:length(methylomesGlobalList), function(x) {
#  lapply(1:length(methylomesGlobalList), function(x) {

    methylome_bin_i <- methylomesGlobalList[[x]] %>%
      dplyr::filter(seqnames == bin_i$chr) %>%
      dplyr::filter(start >= bin_i$start & start <= bin_i$end)

    fwrite(methylome_bin_i,
           file = filePaths_bin_i[x],
           quote = F, sep = "\t", row.names = F, col.names = T)
#  })
  }, mc.cores = length(methylomesGlobalList), mc.preschedule = F)


  # Extract node, generation and methylome info from filePaths_bin_i
  node <- gsub(paste0(inDirBin, "MA\\d+_\\d+_G"), "", filePaths_bin_i)
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
  # e.g., "nodelist_MA1_2_MappedOn_t2t-col.20210610_CpG_Chr1_1_10000.fn"
  
  node_df_GM <- data.frame(filename = filePaths_bin_i,
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
  node_file <- paste0(outDir, "nodelist_MA1_2_MappedOn_", refbase, "_", context, "_",
                      paste0(bin_i, collapse = "_"), ".fn")
  # NOTE: remove this file after use by buildPedigree() due to large file numbers (> 100k)
  fwrite(node_df,
         file = node_file,
         quote = F, sep = ",", row.names = F, col.names = T)


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
  
  edge_file <- paste0(outDir, "edgelist_MA1_2_MappedOn_", refbase, "_", context, "_",
                      paste0(bin_i, collapse = "_"), ".fn")
  # NOTE: remove this file after use by buildPedigree() due to large file numbers (> 100k)
  fwrite(edge_df,
         file = edge_file,
         quote = F, sep = "\t", row.names = F, col.names = T)


  # Build the pedigree of the MA lines in the given genomic bin
  #buildPedigree_file <- paste0(outDir, "buildPedigree_MA1_2_MappedOn_", refbase, "_", context, "_",
  #                             paste0(bin_i, collapse = "_"), ".RData")
  invisible(capture.output(buildPedigree_out <- suppressMessages(buildPedigree(nodelist = node_file,
                                                                               edgelist = edge_file,
                                                                               cytosine = sub("p", "", context),
                                                                               posteriorMaxFilter = 0.99))))
  system(paste0("rm ", inDirBin, "*", "_", paste0(bin_i, collapse = "_"), ".txt"))
  system(paste0("rm ", node_file, " ", edge_file))
  ## NOTE: remove this file after use by buildPedigree() due to large file numbers (> 100k)
  #save(buildPedigree_out,
  #     file = buildPedigree_file)


  # Get the mean, minimum and maximum methylation divergence values at delta.t = 62
  # (MA1_2_mean.D, MA1_2_min.D, MA1_2_max.D)
  pedigree <- buildPedigree_out$Pdata
  delta.t <- pedigree[, 2] + pedigree[, 3] - 2 * pedigree[, 1]
  pedigree <- data.frame(pedigree, delta.t) 
  D_at_dt62 <- pedigree[pedigree$delta.t == 62,]$D.value
  MA1_2_mean.D <- mean(D_at_dt62, na.rm = T)
  MA1_2_min.D <- min(D_at_dt62, na.rm = T)
  MA1_2_max.D <- max(D_at_dt62, na.rm = T)

  data.frame(bin_i,
             MA1_2_mean.D,
             MA1_2_min.D,
             MA1_2_max.D)


  # Plot the pedigree of the MA lines
  
  # plotPedigree doesn't work with sampling.design set to "progenitor.endpoint":
  #plotPedigree(nodelist = node_file,
  #             edgelist = edge_file,
  #             sampling.design = "progenitor.endpoint",
  #             output.dir = plotDir,
  #             plot.width = 5, plot.height = 5, aspect.ratio = 1,
  #             vertex.size = 6, vertex.label = FALSE,
  #             out.pdf = paste0("pedigree_output_MA1_2_MappedOn_", refbase, "_", context, "_",
  #                              paste0(bin_i, collapse = "_")))


  ## Plot divergence time (delta.t) versus methylation divergence (D.value)
  ## Interactive plot for inspecting the divergence data and removing outlier samples (if any)
  #
  #pedigree <- buildPedigree_out$Pdata
  #dt <- pedigree[, 2] + pedigree[, 3] - 2 * pedigree[, 1]
  ##dt_max <- max(dt, na.rm = T)
  #
  #pdf(paste0(plotDir, "divergence_time_vs_methylation_divergence_MA1_2_MappedOn_", refbase, "_", context, "_",
  #           paste0(bin_i, collapse = "_"), ".pdf"))
  #plot(x = dt,
  #     y = pedigree[, "D.value"],
  #     ylab = paste0(gsub("p", "", context), " methyation divergence"),
  #     xlab = expression(paste(Delta, italic("t"), " (generations)", sep = "")))
  #dev.off()


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
  
  #ABneutral_file <- paste0("ABneutral_MA1_2_MappedOn_", refbase, "_", context, "_",
  #                         paste0(bin_i, collapse = "_"))
  #ABneutral_out <- ABneutral(pedigree.data = buildPedigree_out$Pdata,
  #                           p0uu = buildPedigree_out$tmpp0, eqp = buildPedigree_out$tmpp0, eqp.weight = 1,
  #                           Nstarts = 50, out.dir = outDir,
  #                           out.name = ABneutral_file)
  ### NOTE: not an actual .RData file but a list that can be loaded with dget
  ##ABneutral_file_full <- paste0(outDir, ABneutral_file, ".Rdata")
  ##ABneutral_outTest <- dget(ABneutral_file_full)
  ##stopifnot(all.equal(ABneutral_out, ABneutral_outTest))
  ##rm(ABneutral_outTest); gc()
  #head(ABneutral_out$pedigree)
  
  ## Plot estimates of ABneutral model
  ## "In 'ABplot' function you can set parameters to customize the pdf output."
  #ABplot(pedigree.names = paste0(outDir, ABneutral_file, ".Rdata"),
  #       output.dir = plotDir,
  #       out.name = ABneutral_file,
  #       plot.height = 8, plot.width = 11)
  #system(paste0("rm ", outDir, ABneutral_file, ".Rdata"))
  
  
  # Boostrap analysis
  
  # "NOTE: it is recommended to use at least 50 Nboot to achieve best solutions" 
  
  #ABneutral_BOOTfile <- paste0("ABneutral_BOOT_MA1_2_MappedOn_", refbase, "_", context, "_",
  #                             paste0(bin_i, collapse = "_"))
  #ABneutral_BOOTout <- BOOTmodel(pedigree.data = paste0(outDir, ABneutral_file, ".Rdata"),
  #                               Nboot = 50,
  #                               out.dir = outDir,
  #                               out.name = ABneutral_BOOTfile)
  #print(summary(ABneutral_BOOTout))
  ##ABneutral_BOOToutTest <- dget(paste0(outDir, ABneutral_BOOTfile, ".Rdata"))
  ##stopifnot(all.equal(ABneutral_BOOTout, ABneutral_BOOToutTest))
  ##rm(ABneutral_BOOToutTest); gc()

})
#}, mc.cores = detectCores(), mc.preschedule = T)


targetDF <- dplyr::bind_rows(targetDF_list)
#targetDF <- do.call(rbind, targetDF_list)
fwrite(targetDF,
       file = paste0(outDir, "mD_at_dt62_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
                     "_MA1_2_MappedOn_", refbase, "_", chrName, "_", context, ".tsv"),
       quote = F, sep = "\t", row.names = F, col.names = T)
