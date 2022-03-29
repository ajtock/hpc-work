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
# ./alphabeta_per_cytosine_MA1_2_dopar.R t2t-col.20210610 CpG 10000 10000 Chr1

args <- commandArgs(trailingOnly = T)
refbase <- args[1]
context <- args[2]
genomeBinSize <- as.numeric(args[3])
genomeStepSize <- as.numeric(args[4])
chrName <- args[5]
#cores <- as.numeric(args[6]) - 4

#refbase <- "t2t-col.20210610"
#context <- "CpG"
#genomeBinSize <- 1000000
#genomeStepSize <- 100000
#chrName <- "Chr1"
##cores <- as.numeric(48) - 4

options(stringsAsFactors = F)
options(scipen=999)
library(AlphaBeta)
library(dplyr)
library(data.table)
#library(parallel)
#library(doParallel)
#library(doFuture)
#library(snow)
#library(Rmpi)
library(doMPI)
#library(iterators)
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

## Load methylome TXT files for subsetting cytosines (rows) by genomic bin
## NOTE: might be worth changing methimpute script to filter out cytosines
## that are below a given coverage threshold;
## e.g., >= 4; see https://www.pnas.org/doi/10.1073/pnas.1424254112
## However, coverage threshold in https://www.nature.com/articles/s41477-021-01086-7
## seems to be >= 1 (as implemented by default in MethylStar) and with a maximum posterior probability >= 0.99 
#methylomesGlobalList <- lapply(1:length(filePathsGlobal), function(x) {
#  fread(filePathsGlobal[x])
#})
##}, mc.cores = length(filePathsGlobal), mc.preschedule = F)

# Define genomic windows
binDF <- data.frame()
for(i in 1:length(chrs)) {
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

nrow_binDF <- nrow(binDF)


# Define function to get binned methylation divergence values
bin_mD <- function(i, bins) {

  bin_i <- binDF[i,]

  methylomesGlobal_1 <- fread(filePathsGlobal[1])

  methylome1_bin_i <- methylomesGlobal_1 %>%
    dplyr::filter(seqnames == bin_i$chr) %>%
    dplyr::filter(start >= bin_i$start & start <= bin_i$end)

  nrow_methylome1_bin_i <- nrow(methylome1_bin_i)

  rm(methylomesGlobal_1, methylome1_bin_i); invisible(gc())

  if(nrow_methylome1_bin_i > 0) {

    filePaths_bin_i_trunc <- gsub(pattern = "methylome.txt", replacement = paste0("methylome_", paste0(bin_i, collapse = "_"), ".txt"),
                                  x = gsub(pattern = "coverage/report/methimpute/", replacement = "",
                                           x = filePathsGlobal))
    filePaths_bin_i <- paste0(inDirBin, filePaths_bin_i_trunc)

    rm(filePaths_bin_i_trunc); invisible(gc())

    # Write filePaths_bin_i
    # NOTE: remove these files after use by buildPedigree() due to large file numbers (> 1M)
    for(x in 1:length(filePathsGlobal)) {

      methylomesGlobal_x <- fread(filePathsGlobal[x])

      methylome_bin_i <- methylomesGlobal_x %>%
        dplyr::filter(seqnames == bin_i$chr) %>%
        dplyr::filter(start >= bin_i$start & start <= bin_i$end)

      rm(methylomesGlobal_x); invisible(gc())

      write.table(methylome_bin_i,
                  file = filePaths_bin_i[x],
                  quote = F, sep = "\t", row.names = F, col.names = T)

      rm(methylome_bin_i); invisible(gc())
    }


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

    rm(filePaths_bin_i, generation, methylome); invisible(gc())

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
 
    rm(node); invisible(gc())

    node_df_G31 <- node_df_GM[node_df_GM$gen == 31,]
    ##
    rm(node_df_GM); invisible(gc())
    ##
    # Combine node generations into one data.frame
    node_df <- dplyr::bind_rows(mget(sort(ls(pattern = "node_df"))))

    rm(node_df_G0, node_df_G2, node_df_G3, node_df_G30, node_df_G31); invisible(gc())

    node_file <- paste0(outDir, "nodelist_MA1_2_MappedOn_", refbase, "_", context, "_",
                        paste0(bin_i, collapse = "_"), ".fn")
    # NOTE: remove this file after use by buildPedigree() due to large file numbers (> 100k)
    write.table(node_df,
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

    rm(node_df); invisible(gc())

    edge_df <- rbind(edge_df_G0, edge_df_GN)

    rm(edge_df_G0, edge_df_GN); invisible(gc())

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

    rm(dropCols); invisible(gc())

    edge_file <- paste0(outDir, "edgelist_MA1_2_MappedOn_", refbase, "_", context, "_",
                        paste0(bin_i, collapse = "_"), ".fn")
    # NOTE: remove this file after use by buildPedigree() due to large file numbers (> 100k)
    write.table(edge_df,
                file = edge_file,
                quote = F, sep = "\t", row.names = F, col.names = T)

    rm(edge_df); invisible(gc())

    ## Build the pedigree of the MA lines in the given genomic bin
    invisible(capture.output(buildPedigree_out <- suppressMessages(buildPedigree(nodelist = node_file,
                                                                                 edgelist = edge_file,
                                                                                 cytosine = sub("p", "", context),
                                                                                 posteriorMaxFilter = 0.99))))
    #buildPedigree_out <- buildPedigree(nodelist = node_file,
    #                                   edgelist = edge_file,
    #                                   cytosine = sub("p", "", context),
    #                                   posteriorMaxFilter = 0.99)
#    system(paste0("rm ", inDirBin, "*", "_", paste0(bin_i, collapse = "_"), ".txt"))
#    system(paste0("rm ", node_file, " ", edge_file))

    unlink(paste0(inDirBin, "*_", paste0(bin_i, collapse = "_"), ".txt"))
    unlink(node_file)
    unlink(edge_file)

    rm(node_file, edge_file); invisible(gc()) 


    # Get the mean, minimum and maximum methylation divergence values at delta.t = 62
    # (MA1_2_mean.D, MA1_2_min.D, MA1_2_max.D)
    pedigree <- buildPedigree_out$Pdata
    delta.t <- pedigree[, 2] + pedigree[, 3] - 2 * pedigree[, 1]
    pedigree <- data.frame(pedigree, delta.t)

    rm(delta.t); invisible(gc())

    D_at_dt62 <- pedigree[pedigree$delta.t == 62,]$D.value

    rm(pedigree); invisible(gc())

    MA1_2_mean.D <- suppressWarnings(mean(D_at_dt62, na.rm = T))
    MA1_2_min.D <- suppressWarnings(min(D_at_dt62, na.rm = T))
    MA1_2_max.D <- suppressWarnings(max(D_at_dt62, na.rm = T))
    if( MA1_2_mean.D %in% c(NaN, Inf, -Inf) ) { MA1_2_mean.D <- NA }
    if( MA1_2_min.D %in% c(NaN, Inf, -Inf) ) { MA1_2_min.D <- NA }
    if( MA1_2_max.D %in% c(NaN, Inf, -Inf) ) { MA1_2_max.D <- NA }

    rm(D_at_dt62); invisible(gc())

    data.frame(bin_i,
               MA1_2_mean.D,
               MA1_2_min.D,
               MA1_2_max.D)
    
  } else {

    data.frame(bin_i,
               MA1_2_mean.D = NA,
               MA1_2_min.D = NA,
               MA1_2_max.D = NA)

  }

}

bin_mD_test <- function(i, bins) {

  bin_i <- binDF[i,]

  data.frame(bin_i,
             MA1_2_mean.D = NA,
             MA1_2_min.D = NA,
             MA1_2_max.D = NA)

}

#registerDoParallel(cores = cores)
#print("Currently registered parallel backend name, version and cores")
#print(getDoParName())
#print(getDoParVersion())
#print(getDoParWorkers())


#cl <- makeCluster(cores, type = "FORK", outfile = "./info_dopar.log")
#registerDoParallel(cl)
#print("Currently registered parallel backend name, version and cores")
#print(getDoParName())
#print(getDoParVersion())
#print(getDoParWorkers())
#
#stopCluster(cl)


#registerDoFuture()
#cl <- parallel::makeCluster(cores, type = "FORK", outfile = "./info_dopar.log")
#plan(cluster, workers = cl)
##plan(multicore)
#print("Currently registered parallel backend name, version and cores")
#print(getDoParName())
#print(getDoParVersion())
#print(getDoParWorkers())
##plan(future.batchtools::batchjobs_slurm)
###plan(future.batchtools::batchtools_slurm)


# Create and register an MPI cluster
cl <- startMPIcluster()
registerDoMPI(cl)
#library(doFuture)
#registerDoFuture()
#cl <- makeCluster(cores, type = "MPI")
#plan(cluster, workers = cl)
#stopCluster(cl)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())


# Define chunkSize so that each cluster worker gets a single "task chunk"
chunkSize <- ceiling(nrow_binDF / getDoParWorkers())
mpiopts <- list(chunkSize=chunkSize)

start <- proc.time()

targetDF <- foreach(i = icount(nrow_binDF), .options.mpi=mpiopts,
                    .combine = "rbind", .maxcombine = nrow_binDF+1e1,
                    .inorder = F, .errorhandling = "pass") %dopar% {
  bin_mD_test(i = i, bins = binDF)
}

dopar_loop <- proc.time()-start

print(dopar_loop)

targetDF <- targetDF[ with(targetDF,
                           order(chr, start, end)), ]

write.table(targetDF,
            file = paste0(outDir, "mD_at_dt62_genomeBinSize", genomeBinName, "_genomeStepSize", genomeStepName,
                          "_MA1_2_MappedOn_", refbase, "_", chrName, "_", context, ".tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)

# Shutdown the cluster and quit
closeCluster(cl)
#stopCluster(cl)
mpi.quit()
