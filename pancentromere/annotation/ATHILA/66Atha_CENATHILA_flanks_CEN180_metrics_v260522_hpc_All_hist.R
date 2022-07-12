#!/usr/bin/env Rscript

# Compare average CEN180 metrics (HOR membership and divergence) in regions flanking
# centromeric ATHILA and matched centromeric random loci

# Usage:
# conda activate R-4.1.2
# ./66Atha_CENATHILA_flanks_CEN180_metrics_v260522_hpc_All_hist.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 1000 1e4 Flanks
# conda deactivate

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#flankSize <- 1000
#perms <- 1e4
#regionName <- "Flanks"

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
flankSize <- as.numeric(args[2])
perms <- as.numeric(args[3])
regionName <- args[4]

#print(chrName)
#print(flankSize)
#print(perms)
#print(regionName)

# Set minimum possible P-value for permutation test result with
# perms sets of random loci
minPval <- 1 - ( (perms - 1) / perms)

if(floor(log10(flankSize)) + 1 < 4) {
  flankName <- paste0(flankSize, "bp")
} else if(floor(log10(flankSize)) + 1 >= 4 &
          floor(log10(flankSize)) + 1 <= 6) {
  flankName <- paste0(flankSize/1e3, "kb")
} else if(floor(log10(flankSize)) + 1 >= 7) {
  flankName <- paste0(flankSize/1e6, "Mb")
}
flankNamePlot <- paste0(c(strsplit(flankName, split = "")[[1]][1:(length(strsplit(flankName, split = "")[[1]])-2)],
                          "-",
                          strsplit(flankName, split = "")[[1]][(length(strsplit(flankName, split = "")[[1]])-1):(length(strsplit(flankName, split = "")[[1]]))]),
                          collapse = "")

options(stringsAsFactors = F)
source("/home/ajt200/rds/hpc-work/pancentromere/annotation/ATHILA/66Atha_CENATHILA_flanks_CEN180_metrics_v260522_permTestAllList_function_hpc.R")
source("/home/ajt200/rds/hpc-work/pancentromere/annotation/ATHILA/TTSplus.R")
suppressMessages(library(data.table, quietly = T))
suppressMessages(library(GenomicRanges, quietly = T))
suppressMessages(library(segmentSeq, quietly = T))
suppressMessages(library(parallel, quietly = T))
suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(scales, quietly = T))
suppressMessages(library(pals, quietly = T))
suppressMessages(library(RColorBrewer, quietly = T))
suppressMessages(library(ggplot2, quietly = T))
suppressMessages(library(cowplot, quietly = T))
suppressMessages(library(ggbeeswarm, quietly = T))
suppressMessages(library(viridis, quietly = T))
suppressMessages(library(doParallel, quietly = T))
suppressMessages(library(doRNG, quietly = T))
suppressMessages(library(doFuture, quietly = T))
suppressMessages(library(snow, quietly = T))
suppressMessages(library(Rmpi, quietly = T))
suppressMessages(library(doMPI, quietly = T))
suppressMessages(library(iterators, quietly = T))

#library(data.table)
#library(GenomicRanges)
#library(segmentSeq)
#library(parallel)
#library(dplyr)
#library(scales)
#library(pals)
#library(RColorBrewer)
#library(ggplot2)
#library(cowplot)
#library(ggbeeswarm)
#library(viridis)
#library(doParallel)
#library(doRNG)
#library(doFuture)
#library(snow)
#library(Rmpi)
#library(doMPI)
#library(iterators)

# Create and register a doParallel parallel backend
registerDoParallel()
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())
options(future.globals.maxSize = +Inf)

## Create and register a doFuture parallel backend
#registerDoFuture()
#plan(multicore, workers = detectCores())
#print("Currently registered parallel backend name, version and cores")
#print(getDoParName())
#print(getDoParVersion())
#print(getDoParWorkers())
#options(future.globals.maxSize = +Inf)

## Create and register an MPI cluster
#cl <- startMPIcluster(verbose = T, logdir = "logs/", bcast = F)
##cl <- startMPIcluster(verbose = F, bcast = T)
#registerDoMPI(cl)
##registerDoFuture()
##cl <- snow::makeCluster(mpi.universe.size() - 1, type = "MPI", outfile = "logs/alphabeta_per_cytosine_MA1_2_CpG_Chr2_snow_mpi.log")
##plan(cluster, workers = cl)
###cl <- snow::makeMPIcluster(count = mpi.comm.size(0) - 1, outfile = "logs/alphabeta_per_cytosine_MA1_2_CpG_Chr2_snow_mpi.log")
###cl <- parallel::makeCluster(cores, type = "FORK", outfile = "logs/alphabeta_per_cytosine_MA1_2_CpG_Chr2_parallel_fork.log")
###plan(multicore)
#print("Currently registered parallel backend name, version and cores")
#print(getDoParName())
#print(getDoParVersion())
#print(getDoParWorkers())

#options(future.globals.maxSize = +Inf)
# Define chunkSize so that each cluster worker gets a single "task chunk"
# (i.e. one task chunk covering multiple loop iterations (rows of binDF)),
# which is faster than if each cluster worker gets multiple task chunks
# (i.e., one task chunk per loop iteration (row of binDF))
chunkSize <- ceiling(perms / getDoParWorkers())
#initWorkers <- function() options(scipen = 999, stringsAsFactors = F)
mpiopts <- list(chunkSize = chunkSize)

plotDir <- paste0("ATHILA/plots/")
plotDirHORlengthsSum <- paste0(plotDir, "HORlengthsSum/")
plotDirHORcount <- paste0(plotDir, "HORcount/")
plotDirWeightedConsensusScore <- paste0(plotDir, "WeightedConsensusScore/")
plotDirEditDistance <- paste0(plotDir, "EditDistance/")
plotDirAllMetrics <- paste0(plotDir, "AllMetrics/")
plotDirAllAccessions <- paste0(plotDir, "AllAccessions/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))
system(paste0("[ -d ", plotDirHORlengthsSum, " ] || mkdir -p ", plotDirHORlengthsSum))
system(paste0("[ -d ", plotDirHORcount, " ] || mkdir -p ", plotDirHORcount))
system(paste0("[ -d ", plotDirWeightedConsensusScore, " ] || mkdir -p ", plotDirWeightedConsensusScore))
system(paste0("[ -d ", plotDirEditDistance, " ] || mkdir -p ", plotDirEditDistance))
system(paste0("[ -d ", plotDirAllMetrics, " ] || mkdir -p ", plotDirAllMetrics))
system(paste0("[ -d ", plotDirAllAccessions, " ] || mkdir -p ", plotDirAllAccessions))

# Load and define accession names
acc_full <- system("ls /home/ajt200/rds/hpc-work/pancentromere/annotation/CEN180/repeats/*.fa*", intern = T)
acc_full <- gsub("/home/ajt200/rds/hpc-work/pancentromere/annotation/CEN180/repeats/cen180.consensus.repetitiveness", "", acc_full)
acc_full <- acc_full[-grep("SUPER_", acc_full)]
acc_uniq <- unique( gsub("\\.fa\\..+", "", acc_full) )
acc_uniq_len <- NULL
for(i in acc_uniq) {
  acc_uniq_len <- c(acc_uniq_len, length( grep(i , acc_full)) )
}
acc <- acc_uniq[which(acc_uniq_len == length(chrName))]

chrs_list <- foreach(x = 1:length(acc), .inorder = T) %dopar% {
  acc_chrs <- read.table(paste0("/home/ajt200/rds/hpc-work/pancentromere/assemblies/",
                                acc[x], ".fa.fai"),
                         header = F)[,1]
  acc_chrs <- gsub("_RagTag_RagTag", "", acc_chrs)
  acc_chrs <- gsub("chr", "Chr", acc_chrs)
  acc_chrs <- gsub("SUPER_", "Chr", acc_chrs)
  acc_chrs <- acc_chrs[which(acc_chrs %in% chrName)]
  acc_chrs[sort.int(acc_chrs, index.return = T)$ix]
}

chrLens_list <- foreach(x = 1:length(acc), .inorder = T) %dopar% {
  acc_chrs <- read.table(paste0("/home/ajt200/rds/hpc-work/pancentromere/assemblies/",
                                acc[x], ".fa.fai"),
                         header = F)[,1]
  acc_chrs <- gsub("_RagTag_RagTag", "", acc_chrs)
  acc_chrs <- gsub("chr", "Chr", acc_chrs)
  acc_chrs <- gsub("SUPER_", "Chr", acc_chrs)
  acc_chrLens <- read.table(paste0("/home/ajt200/rds/hpc-work/pancentromere/assemblies/",
                                   acc[x], ".fa.fai"),
                            header = F)[,2]
  acc_chrLens <- acc_chrLens[which(acc_chrs %in% chrName)]
  acc_chrs <- acc_chrs[which(acc_chrs %in% chrName)]
  print(acc_chrs)
  print(acc_chrLens)
  acc_chrLens <- acc_chrLens[sort.int(acc_chrs, index.return = T)$ix]
  acc_chrs <- acc_chrs[sort.int(acc_chrs, index.return = T)$ix]
  print(acc_chrs)
  print(acc_chrLens)
  acc_chrLens
}

# Load CEN coordinates for each accession
CEN <- read.csv(paste0("/home/ajt200/rds/hpc-work/pancentromere/centromeric_coordinates/",
                       "centromere_manual_EDTA4_fa.csv"),
                header = T)
CEN$fasta.name <- gsub(".fa", "", CEN$fasta.name)
CEN_GR_list <- foreach(x = 1:length(acc), .inorder = T) %dopar% {
  acc_CEN <- CEN[grep(acc[x], CEN$fasta.name),]
  GRanges(seqnames = acc_CEN$chr,
          ranges = IRanges(start = acc_CEN$start,
                           end = acc_CEN$end),
          strand = "*",
          acc = acc_CEN$fasta.name,
          region = acc_CEN$region)
}

# Load CEN180 for each accession
CEN180_list <- foreach(x = 1:length(acc), .inorder = T) %dopar% {
  tab_list <- lapply(1:length(chrName), function(y) {
    read.csv(paste0("/home/ajt200/rds/hpc-work/pancentromere/annotation/CEN180/repeats/",
                    "cen180.consensus.repetitiveness", acc[x], ".fa.", chrName[y], ".csv"),
             header = T)
  })
  if(length(chrName) > 1) {
    tab <- dplyr::bind_rows(tab_list)
  } else {
    tab <- tab_list[[1]]
  }
  tab <- tab[which(tab$class == "aTha178"),]
  colnames(tab)[10] <- "chr"
  tab$chr <- gsub("_RagTag_RagTag", "", tab$chr)
  tab$chr <- gsub("chr", "Chr", tab$chr)
  tab$chr <- gsub("SUPER_", "Chr", tab$chr)
  tab$assembly <- gsub("\\.fasta", "", tab$assembly)
  tab$assembly <- gsub("\\.fa", "", tab$assembly)
  tab <- tab[
             with( tab, order(chr, start, end) ),
            ]
  tab
}

CEN180_GR_list <- foreach(x = 1:length(acc), .inorder = T) %dopar% { 
  GRanges(seqnames = as.character(CEN180_list[[x]]$chr),
          ranges = IRanges(start = as.integer(CEN180_list[[x]]$start),
                           end = as.integer(CEN180_list[[x]]$end)),
          strand = as.character(CEN180_list[[x]]$strand),
          acc = as.character(CEN180_list[[x]]$assembly),
          HORlengthsSum = as.numeric(CEN180_list[[x]]$HORlengthsSum),
          HORcount = as.numeric(CEN180_list[[x]]$HORcount),
          WeightedConsensusScore = as.numeric(CEN180_list[[x]]$weighted.consensus.score),
          EditDistance = as.numeric(CEN180_list[[x]]$edit.distance))
}

# Load CENATHILA for each accession
CENATHILA_list <- foreach(x = 1:length(acc), .inorder = T) %dopar% {
  tab <- read.table(paste0("/home/ajt200/rds/hpc-work/pancentromere/annotation/ATHILA/ATHILA/",
                           acc[x], "/CENATHILA_in_", acc[x], "_",
                           paste0(chrName, collapse = "_"), ".bed"),
                    header = F)
  tab$V2 <- tab$V2+1
  colnames(tab) <- c("chr", "start", "end", "name", "phylo", "strand")
  tab <- tab[
             with( tab, order(chr, start, end) ),
            ]
  tab
}

CENATHILA_DF_phylo <- dplyr::bind_rows(CENATHILA_list) 
phylo <- sort(unique(CENATHILA_DF_phylo$phylo))

phylo_colFun <- cols25(n = 18)[-c(7:16)][1:length(phylo)]
stopifnot(length(phylo_colFun) == length(phylo))
names(phylo_colFun) <- phylo

CENATHILA_GR_list <- foreach(x = 1:length(acc), .inorder = T) %dopar% {
  GRanges(seqnames = as.character(CENATHILA_list[[x]]$chr),
          ranges = IRanges(start = as.integer(CENATHILA_list[[x]]$start),
                           end = as.integer(CENATHILA_list[[x]]$end)),
          strand = as.character(CENATHILA_list[[x]]$strand),
          phylo = as.character(CENATHILA_list[[x]]$phylo))
}


# Define function to select randomly positioned loci of the same
# width distribution as CENgapAllAthila_bed
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Function to define, for each accession, "perms" sets of centromeric random loci (acc_CENranLoc_GR)
# of the same number and width distribution as acc_CENATHILA_GR
defineCENranLoc <- function(acc_idx, chrs_list, chrLens_list, CEN_GR_list, CENATHILA_GR_list, seed) {
  print(acc[acc_idx])

  acc_chrs <- chrs_list[[acc_idx]]
  acc_chrLens <- chrLens_list[[acc_idx]]
  acc_CEN_GR <- CEN_GR_list[[acc_idx]]
  acc_CENATHILA_GR <- CENATHILA_GR_list[[acc_idx]]

  # Apply ranLocStartSelect() on a per-chromosome basis so that
  # acc_CENranLoc_GR contains the same number of loci per chromosome as acc_CENATHILA_GR
  acc_CENranLoc_GR <- GRanges()
  for(j in 1:length(acc_chrs)) {
    print(acc_chrs[j])

    chr_acc_CEN_GR <- acc_CEN_GR[seqnames(acc_CEN_GR) == acc_chrs[j]]

    chr_acc_CENATHILA_GR <- acc_CENATHILA_GR[seqnames(acc_CENATHILA_GR) == acc_chrs[j]]
    if(length(chr_acc_CENATHILA_GR) > 0) {
      set.seed(seed + 1e6)
      chr_acc_CENranLoc_Start <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(chr_acc_CEN_GR), function(x) {
                                                                          ( start(chr_acc_CEN_GR[x]) + max(width(chr_acc_CENATHILA_GR)) ) :
                                                                          ( end(chr_acc_CEN_GR[x]) - max(width(chr_acc_CENATHILA_GR)) )
                                                                        })),
                                                   n = length(chr_acc_CENATHILA_GR))
      chr_acc_CENranLoc_GR <- GRanges(seqnames = acc_chrs[j],
                                      ranges = IRanges(start = chr_acc_CENranLoc_Start,
                                                       width = width(chr_acc_CENATHILA_GR)),
                                      strand = strand(chr_acc_CENATHILA_GR),
                                      phylo = as.character(chr_acc_CENATHILA_GR$phylo))
      acc_CENranLoc_GR <- append(acc_CENranLoc_GR, chr_acc_CENranLoc_GR)
    }
  }
  stopifnot(identical(width(acc_CENranLoc_GR), width(acc_CENATHILA_GR)))
  acc_CENranLoc_GR
}

# Define seed so that random selections are reproducible
set.seed(76492749)
CENranLoc_GR_acc_list <- foreach(acc_idx = 1:length(acc), .inorder = T) %do% {
  acc_perms_GR_list <- foreach(seed = iter(1:perms),
#                               .options.mpi = mpiopts,
                               .multicombine = T,
                               .maxcombine = perms+1e1,
                               .inorder = F) %dorng% {
    defineCENranLoc(acc_idx = acc_idx,
                    chrs_list = chrs_list,
                    chrLens_list = chrLens_list,
                    CEN_GR_list = CEN_GR_list,
                    CENATHILA_GR_list = CENATHILA_GR_list,
                    seed = seed)
  }
  acc_perms_GR_list
}

# Function to calculate mean CEN180 metrics for regions
# upstream and downstream of CENATHILA and CENranLoc bodies
CEN180metricsAtCENATHILA <- function(acc_idx, CEN180_GR, CENATHILA_GR, featureName) {
  print(acc[acc_idx])

  # Get flankSize bp upstream of start coordinates
  CENATHILA_up_GR <- promoters(CENATHILA_GR, upstream = flankSize, downstream = 0)
  CENATHILA_down_GR <- TTSplus(CENATHILA_GR, upstream = -1, downstream = flankSize)
  stopifnot(identical(CENATHILA_GR$phylo, CENATHILA_up_GR$phylo))
  stopifnot(identical(CENATHILA_GR$phylo, CENATHILA_down_GR$phylo))
  stopifnot(identical(seqnames(CENATHILA_GR), seqnames(CENATHILA_up_GR)))
  stopifnot(identical(seqnames(CENATHILA_GR), seqnames(CENATHILA_down_GR)))
  stopifnot(identical(strand(CENATHILA_GR), strand(CENATHILA_up_GR)))
  stopifnot(identical(strand(CENATHILA_GR), strand(CENATHILA_down_GR)))

  CENATHILA_CEN180_metrics <- data.frame()
  for(i in 1:length(chrName)) {
    print(chrName[i])

    CEN180_GR_chr <- CEN180_GR[seqnames(CEN180_GR) == chrName[i]]
    CENATHILA_GR_chr <- CENATHILA_GR[seqnames(CENATHILA_GR) == chrName[i]]
    CENATHILA_up_GR_chr <- CENATHILA_up_GR[seqnames(CENATHILA_up_GR) == chrName[i]]
    CENATHILA_down_GR_chr <- CENATHILA_down_GR[seqnames(CENATHILA_down_GR) == chrName[i]]
    stopifnot(identical(CENATHILA_GR_chr$phylo, CENATHILA_up_GR_chr$phylo))
    stopifnot(identical(CENATHILA_GR_chr$phylo, CENATHILA_down_GR_chr$phylo))
    stopifnot(identical(strand(CENATHILA_GR_chr), strand(CENATHILA_up_GR_chr)))
    stopifnot(identical(strand(CENATHILA_GR_chr), strand(CENATHILA_down_GR_chr)))

    if(length(CENATHILA_GR_chr) > 0) {
      CENATHILA_up_CEN180 <- getOverlaps(coordinates = CENATHILA_up_GR_chr,
                                         segments = CEN180_GR_chr,
                                         overlapType = "overlapping",
                                         whichOverlaps = TRUE,
                                         ignoreStrand = TRUE)
      CENATHILA_up_CEN180_HORlengthsSum <- sapply(1:length(CENATHILA_up_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_up_CEN180[[x]]]$HORlengthsSum, na.rm = T)
      })
      CENATHILA_up_CEN180_HORcount <- sapply(1:length(CENATHILA_up_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_up_CEN180[[x]]]$HORcount, na.rm = T)
      })
      CENATHILA_up_CEN180_WeightedConsensusScore <- sapply(1:length(CENATHILA_up_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_up_CEN180[[x]]]$WeightedConsensusScore, na.rm = T)
      })
      CENATHILA_up_CEN180_EditDistance <- sapply(1:length(CENATHILA_up_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_up_CEN180[[x]]]$EditDistance, na.rm = T)
      })

      CENATHILA_down_CEN180 <- getOverlaps(coordinates = CENATHILA_down_GR_chr,
                                           segments = CEN180_GR_chr,
                                           overlapType = "overlapping",
                                           whichOverlaps = TRUE,
                                           ignoreStrand = TRUE)
      CENATHILA_down_CEN180_HORlengthsSum <- sapply(1:length(CENATHILA_down_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_down_CEN180[[x]]]$HORlengthsSum, na.rm = T)
      })
      CENATHILA_down_CEN180_HORcount <- sapply(1:length(CENATHILA_down_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_down_CEN180[[x]]]$HORcount, na.rm = T)
      })
      CENATHILA_down_CEN180_WeightedConsensusScore <- sapply(1:length(CENATHILA_down_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_down_CEN180[[x]]]$WeightedConsensusScore, na.rm = T)
      })
      CENATHILA_down_CEN180_EditDistance <- sapply(1:length(CENATHILA_down_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_down_CEN180[[x]]]$EditDistance, na.rm = T)
      })

      # Don't concatenate CENATHILA_up_GR_chr and CENATHILA_down_GR_chr;
      # better to keep start and end coordinates and annotate with "Upstream" or "Downstream" instead
      ##CENATHILA_reg_GR_chr <- c(CENATHILA_up_GR_chr, CENATHILA_down_GR_chr)
      stopifnot(length(CENATHILA_up_GR_chr) == length(CENATHILA_down_GR_chr))
      CENATHILA_reg_GR_chr <- c(CENATHILA_GR_chr, CENATHILA_GR_chr)
      CENATHILA_chr <- data.frame(CENATHILA_reg_GR_chr,
                                  feature = rep(featureName, length(CENATHILA_reg_GR_chr)),
                                  accession = rep(acc[acc_idx], length(CENATHILA_reg_GR_chr)),
                                  region = rep(c("Upstream", "Downstream"), each = length(CENATHILA_GR_chr)),
                                  HORlengthsSum = c(CENATHILA_up_CEN180_HORlengthsSum, CENATHILA_down_CEN180_HORlengthsSum),
                                  HORcount = c(CENATHILA_up_CEN180_HORcount, CENATHILA_down_CEN180_HORcount),
                                  WeightedConsensusScore = c(CENATHILA_up_CEN180_WeightedConsensusScore, CENATHILA_down_CEN180_WeightedConsensusScore),
                                  EditDistance = c(CENATHILA_up_CEN180_EditDistance, CENATHILA_down_CEN180_EditDistance))
      colnames(CENATHILA_chr)[1] <- "chr"
      CENATHILA_chr$chr <- as.character(CENATHILA_chr$chr)
      CENATHILA_chr$strand <- as.character(CENATHILA_chr$strand)
      colnames(CENATHILA_chr)[which(colnames(CENATHILA_chr) == "phylo")] <- "Family"
      colnames(CENATHILA_chr)[which(colnames(CENATHILA_chr) == "feature")] <- "Feature"
      colnames(CENATHILA_chr)[which(colnames(CENATHILA_chr) == "region")] <- "Region"

      CENATHILA_CEN180_metrics <- rbind(CENATHILA_CEN180_metrics, CENATHILA_chr)
    } 
  }
  CENATHILA_CEN180_metrics
}

# For each acc, apply CEN180metricsAtCENATHILA() function
# CENATHILA
CENATHILA_CEN180_metrics_list <- foreach(acc_idx = 1:length(acc), .inorder = T) %dopar% {
  CEN180metricsAtCENATHILA(acc_idx = acc_idx,
                           CEN180_GR = CEN180_GR_list[[acc_idx]],
                           CENATHILA_GR = CENATHILA_GR_list[[acc_idx]],
                           featureName = "CENATHILA")
}

# CENranLoc
CENranLoc_CEN180_metrics_list <- foreach(acc_idx = 1:length(acc), .inorder = T) %do% {
  foreach(x = iter(1:perms),
#          .options.mpi = mpiopts,
          .multicombine = T,
          .maxcombine = perms+1e1,
          .inorder = F) %dopar% {
    CEN180metricsAtCENATHILA(acc_idx = acc_idx,
                             CEN180_GR = CEN180_GR_list[[acc_idx]],
                             CENATHILA_GR = CENranLoc_GR_acc_list[[acc_idx]][[x]],
                             featureName = "CENranLoc")
  }
}


# Do permutation tests across all accessions
source("/home/ajt200/rds/hpc-work/pancentromere/annotation/ATHILA/66Atha_CENATHILA_flanks_CEN180_metrics_v260522_permTestAllList_function_hpc.R")

# HORlengthsSum
permTestAllList_HORlengthsSum <- permTestAllList(
                                                 CENATHILA_CEN180_metrics_list = CENATHILA_CEN180_metrics_list,
                                                 CENranLoc_CEN180_metrics_list = CENranLoc_CEN180_metrics_list,
                                                 region_name = regionName,
                                                 metric_name = "HORlengthsSum"
                                                )
permTestAllList_HORlengthsSum_permDistDF <- data.frame(
                                                       Accession = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                         rep(permTestAllList_HORlengthsSum[[y]]@accession,
                                                             times = length(permTestAllList_HORlengthsSum[[y]]@permDist))
                                                       })),
                                                       Family = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                         rep(permTestAllList_HORlengthsSum[[y]]@fam,
                                                             times = length(permTestAllList_HORlengthsSum[[y]]@permDist))
                                                       })),
                                                       Metric = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                         rep(permTestAllList_HORlengthsSum[[y]]@metric,
                                                             times = length(permTestAllList_HORlengthsSum[[y]]@permDist))
                                                       })),
                                                       Region = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                         rep(permTestAllList_HORlengthsSum[[y]]@region,
                                                             times = length(permTestAllList_HORlengthsSum[[y]]@permDist))
                                                       })),
                                                       Permuted = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                         permTestAllList_HORlengthsSum[[y]]@permDist
                                                       }))
                                                      ) 
permTestAllList_HORlengthsSum_permDF <- data.frame(
                                                   Accession = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@accession
                                                   })),
                                                   Family = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@fam
                                                   })),
                                                   Metric = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@metric
                                                   })),
                                                   Region = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@region
                                                   })),
                                                   Observed = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@observed
                                                   })),
                                                   Expected = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@expected
                                                   })),
                                                   Alpha = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@alphaThreshold
                                                   })),
                                                   Pvalue = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@pval
                                                   })),
                                                   Alternative = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@alternative
                                                   })),
                                                   Features = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@features
                                                   })),
                                                   Log2ObservedExpected = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@log2obsexp
                                                   })),
                                                   Log2Alpha = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@log2alpha
                                                   }))
                                                  ) 

# EditDistance
permTestAllList_EditDistance <- permTestAllList(
                                                CENATHILA_CEN180_metrics_list = CENATHILA_CEN180_metrics_list,
                                                CENranLoc_CEN180_metrics_list = CENranLoc_CEN180_metrics_list,
                                                region_name = regionName,
                                                metric_name = "EditDistance"
                                               )
permTestAllList_EditDistance_permDistDF <- data.frame(
                                                      Family = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                        rep(permTestAllList_EditDistance[[y]]@fam,
                                                            times = length(permTestAllList_EditDistance[[y]]@permDist))
                                                      })),
                                                      Accession = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                        rep(permTestAllList_EditDistance[[y]]@accession,
                                                            times = length(permTestAllList_EditDistance[[y]]@permDist))
                                                      })),
                                                      Metric = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                        rep(permTestAllList_EditDistance[[y]]@metric,
                                                            times = length(permTestAllList_EditDistance[[y]]@permDist))
                                                      })),
                                                      Region = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                        rep(permTestAllList_EditDistance[[y]]@region,
                                                            times = length(permTestAllList_EditDistance[[y]]@permDist))
                                                      })),
                                                      Permuted = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                        permTestAllList_EditDistance[[y]]@permDist
                                                      }))
                                                     ) 
permTestAllList_EditDistance_permDF <- data.frame(
                                                  Family = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@fam
                                                  })),
                                                  Alternative = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@alternative
                                                  })),
                                                  Alpha = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@alphaThreshold
                                                  })),
                                                  Pvalue = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@pval
                                                  })),
                                                  Observed = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@observed
                                                  })),
                                                  Expected = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@expected
                                                  })),
                                                  Log2ObservedExpected = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@log2obsexp
                                                  })),
                                                  Log2Alpha = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@log2alpha
                                                  })),
                                                  Features = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@features
                                                  })),
                                                  Accession = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@accession
                                                  })),
                                                  Metric = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@metric
                                                  })),
                                                  Region = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@region
                                                  }))
                                                 ) 

## HORcount
#permTestAllList_HORcount <- permTestAllList(
#                                            CENATHILA_CEN180_metrics_list = CENATHILA_CEN180_metrics_list,
#                                            CENranLoc_CEN180_metrics_list = CENranLoc_CEN180_metrics_list,
#                                            region_name = regionName,
#                                            metric_name = "HORcount"
#                                           )
#permTestAllList_HORcount_permDistDF <- data.frame(
#                                                  Family = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                    rep(permTestAllList_HORcount[[y]]@fam,
#                                                        times = length(permTestAllList_HORcount[[y]]@permDist))
#                                                  })),
#                                                  Metric = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                    rep(permTestAllList_HORcount[[y]]@metric,
#                                                        times = length(permTestAllList_HORcount[[y]]@permDist))
#                                                  })),
#                                                  Permuted = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                    permTestAllList_HORcount[[y]]@permDist
#                                                  }))
#                                                 ) 
#permTestAllList_HORcount_permDistDF <- data.frame(
#                                                  Family = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                    rep(permTestAllList_HORcount[[y]]@fam,
#                                                        times = length(permTestAllList_HORcount[[y]]@permDist))
#                                                  })),
#                                                  Accession = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                    rep(permTestAllList_HORcount[[y]]@accession,
#                                                        times = length(permTestAllList_HORcount[[y]]@permDist))
#                                                  })),
#                                                  Metric = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                    rep(permTestAllList_HORcount[[y]]@metric,
#                                                        times = length(permTestAllList_HORcount[[y]]@permDist))
#                                                  })),
#                                                  Region = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                    rep(permTestAllList_HORcount[[y]]@region,
#                                                        times = length(permTestAllList_HORcount[[y]]@permDist))
#                                                  })),
#                                                  Permuted = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                    permTestAllList_HORcount[[y]]@permDist
#                                                  }))
#                                                 ) 
#permTestAllList_HORcount_permDF <- data.frame(
#                                              Family = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@fam
#                                              })),
#                                              Alternative = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@alternative
#                                              })),
#                                              Alpha = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@alphaThreshold
#                                              })),
#                                              Pvalue = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@pval
#                                              })),
#                                              Observed = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@observed
#                                              })),
#                                              Expected = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@expected
#                                              })),
#                                              Log2ObservedExpected = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@log2obsexp
#                                              })),
#                                              Log2Alpha = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@log2alpha
#                                              })),
#                                              Features = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@features
#                                              })),
#                                              Accession = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@accession
#                                              })),
#                                              Metric = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@metric
#                                              })),
#                                              Region = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@region
#                                              }))
#                                             ) 
#
## WeightedConsensusScore
#permTestAllList_WeightedConsensusScore <- permTestAllList(
#                                                          CENATHILA_CEN180_metrics_list = CENATHILA_CEN180_metrics_list,
#                                                          CENranLoc_CEN180_metrics_list = CENranLoc_CEN180_metrics_list,
#                                                          region_name = regionName,
#                                                          metric_name = "WeightedConsensusScore"
#                                                         )
#permTestAllList_WeightedConsensusScore_permDistDF <- data.frame(
#                                                                Family = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                                  rep(permTestAllList_WeightedConsensusScore[[y]]@fam,
#                                                                      times = length(permTestAllList_WeightedConsensusScore[[y]]@permDist))
#                                                                })),
#                                                                Accession = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                                  rep(permTestAllList_WeightedConsensusScore[[y]]@accession,
#                                                                      times = length(permTestAllList_WeightedConsensusScore[[y]]@permDist))
#                                                                })),
#                                                                Metric = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                                  rep(permTestAllList_WeightedConsensusScore[[y]]@metric,
#                                                                      times = length(permTestAllList_WeightedConsensusScore[[y]]@permDist))
#                                                                })),
#                                                                Region = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                                  rep(permTestAllList_WeightedConsensusScore[[y]]@region,
#                                                                      times = length(permTestAllList_WeightedConsensusScore[[y]]@permDist))
#                                                                })),
#                                                                Permuted = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                                  permTestAllList_WeightedConsensusScore[[y]]@permDist
#                                                                }))
#                                                               ) 
#permTestAllList_WeightedConsensusScore_permDF <- data.frame(
#                                                            Family = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@fam
#                                                            })),
#                                                            Alternative = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@alternative
#                                                            })),
#                                                            Alpha = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@alphaThreshold
#                                                            })),
#                                                            Pvalue = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@pval
#                                                            })),
#                                                            Observed = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@observed
#                                                            })),
#                                                            Expected = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@expected
#                                                            })),
#                                                            Log2ObservedExpected = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@log2obsexp
#                                                            })),
#                                                            Log2Alpha = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@log2alpha
#                                                            })),
#                                                            Features = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@features
#                                                            })),
#                                                            Accession = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@accession
#                                                            })),
#                                                            Metric = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@metric
#                                                            })),
#                                                            Region = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@region
#                                                            }))
#                                                           ) 


# Combine into one data.frame for plotting with ggplot2
combined_permDistDF <- rbind(
                             permTestAllList_HORlengthsSum_permDistDF,
                             permTestAllList_EditDistance_permDistDF
#                             permTestAllList_HORcount_permDistDF,
#                             permTestAllList_WeightedConsensusScore_permDistDF
                            )
combined_permDF <- rbind(
                         permTestAllList_HORlengthsSum_permDF,
                         permTestAllList_EditDistance_permDF
#                         permTestAllList_HORcount_permDF,
#                         permTestAllList_WeightedConsensusScore_permDF
                        )

write.table(combined_permDistDF,
            file = paste0(plotDirAllMetrics,
                          "CENATHILA_", flankName, "_regions_CEN180_AllMetrics_all_accessions_combined_",
                          paste0(chrName, collapse = "_"), "_",
                          perms, "perms_",
                          regionName,
                          "_permDistDF.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
combined_permDistDF <- fread(file = paste0(plotDirAllMetrics,
                                           "CENATHILA_", flankName, "_regions_CEN180_AllMetrics_all_accessions_combined_",
                                           paste0(chrName, collapse = "_"), "_",
                                           perms, "perms_",
                                           regionName,
                                           "_permDistDF.tsv"),
                             data.table = F)
write.table(combined_permDF,
            file = paste0(plotDirAllMetrics,
                          "CENATHILA_", flankName, "_regions_CEN180_AllMetrics_all_accessions_combined_",
                          paste0(chrName, collapse = "_"), "_",
                          perms, "perms_",
                          regionName,
                          "_permDF.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
combined_permDF <- fread(file = paste0(plotDirAllMetrics,
                                       "CENATHILA_", flankName, "_regions_CEN180_AllMetrics_all_accessions_combined_",
                                       paste0(chrName, collapse = "_"), "_",
                                       perms, "perms_",
                                       regionName,
                                       "_permDF.tsv"),
                         data.table = F)

combined_permDistDF$Family <- factor(combined_permDistDF$Family,
                                     levels = rev(sort(unique(combined_permDistDF$Family))))
combined_permDistDF$Accession <- factor(combined_permDistDF$Accession,
                                        levels = sort(unique(combined_permDistDF$Accession)))
combined_permDistDF$Metric <- factor(combined_permDistDF$Metric,
                                     levels =  unique(combined_permDistDF$Metric))
combined_permDistDF$Region <- factor(combined_permDistDF$Region,
                                     levels = unique(combined_permDistDF$Region))
combined_permDF$Family <- factor(combined_permDF$Family,
                                 levels = rev(sort(unique(combined_permDF$Family))))
combined_permDF$Accession <- factor(combined_permDF$Accession,
                                    levels = sort(unique(combined_permDF$Accession)))
combined_permDF$Metric <- factor(combined_permDF$Metric,
                                 levels =  unique(combined_permDF$Metric))
combined_permDF$Region <- factor(combined_permDF$Region,
                                 levels = unique(combined_permDF$Region))


vp_all <- ggplot(data = combined_permDistDF,
                 mapping = aes(x = Family,
                               y = Permuted)) +
  xlab("Family") +
  ylab(bquote(italic("CEN178") ~ "metric in" ~ .(flankNamePlot) ~ "flanking regions")) +
  geom_violin(trim = F,
              scale = "count",
              fill = "grey80") +
  geom_point(data = combined_permDF,
             mapping = aes(x = Family,
                           y = Observed),
             position = position_nudge(x = -0.04),
             shape = "|", colour = "dodgerblue2", size = 12) +
  geom_point(data = combined_permDF,
             mapping = aes(x = Family,
                           y = Expected),
             position = position_nudge(x = -0.04),
             shape = "|", colour = "black", size = 12) +
  geom_point(data = combined_permDF,
             mapping = aes(x = Family,
                           y = Alpha),
             position = position_nudge(x = -0.04),
             shape = "|", colour = "darkorange2", size = 12) +

  coord_flip() +
  theme_bw() +
  theme(
#        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks.y = element_line(size = 0.5, colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black", face = "italic", hjust = 0.5, vjust = 0.5, angle = 0),
        axis.title.y = element_text(size = 20, colour = "black"),
#        axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.ticks.x = element_line(size = 0.5, colour = "black"),
        axis.text.x = element_text(size = 20, colour = "black", hjust = 0.5, vjust = 0.5, angle = 0),
        axis.title.x = element_text(size = 20, colour = "black", margin = margin(t = 20)),
        strip.text.x = element_text(size = 20, colour = "black"),
        strip.text.y = element_text(size = 20, colour = "black"),
        strip.background = element_rect(fill = "white", colour = "white"),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        plot.title = element_text(size = 15, colour = "black", hjust = 0.5)) +
  ggtitle(bquote("Permutations =" ~
                 .(prettyNum(perms,
                             big.mark = ",",
                             trim = T)) ~
                 "sets of randomly positioned centromeric loci"))

vp_all <- vp_all +
  facet_grid(cols = vars(Metric), scales = "free_x")
ggsave(paste0(plotDirAllMetrics,
              "CENATHILA_", flankName, "_regions_CEN180_AllMetrics_all_accessions_violin_",
              paste0(chrName, collapse = "_"), "_",
              perms, "perms_",
              paste0(unique(as.vector(combined_permDF$Region)), collapse = "_"),
              ".pdf"),
       plot = vp_all,
       width = 4*length(unique(combined_permDF$Metric)), height = 8, limitsize = F)



combined$Accession <- factor(combined$Accession,
                             levels = c(paste0(length(unique(combined$Accession))-1, " accessions"),
                                        sort(unique(combined$Accession))[-grep(" accessions", sort(unique(combined$Accession)))]))
combined$Metric <- factor(combined$Metric,
                          levels = unique(combined$Metric))
combined$Region <- factor(combined$Region,
                          levels = unique(combined$Region))
combined$Family <- factor(combined$Family,
                          levels = sort(unique(combined$Family)))


# Per-accession plots (including All accessions plot)
bp_per <- ggplot(data = combined,
                 mapping = aes(x = Family,
                               y = log2obsexp,
                               fill = Metric)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
#  scale_fill_discrete(name = "") +
  scale_fill_brewer(name = "",
                    palette = "Dark2") +
  geom_point(mapping = aes(x = Family,
                           y = log2alpha),
             position = position_dodge(0.9),
             shape = "-", colour  = "grey70", size = 12) +
  xlab(bquote(italic("CENATHILA") ~ .(flankNamePlot) ~ "flanking regions")) +
  ylab(bquote("Log"[2] * "(observed/expected)" ~ italic("CEN178") ~ "metric")) +
#  scale_y_continuous(limits = c(-4.0, 4.0)) +
  scale_x_discrete(position = "bottom") +
  guides(fill = guide_legend(direction = "vertical",
                             label.position = "right",
                             label.theme = element_text(size = 16, hjust = 0, vjust = 0.5, angle = 0),
                             nrow = length(unique(combined$Family)),
                             byrow = TRUE)) +

  theme_bw() +
  theme(axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks.y = element_line(size = 0.5, colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black", hjust = 0.5, vjust = 0.5, angle = 90),
        axis.title.y = element_text(size = 20, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 20, colour = "black", hjust = 0.5, vjust = 0.5, angle = 45),
        axis.title.x = element_text(size = 20, colour = "black", margin = margin(t = 20)),
        strip.text.x = element_text(size = 20, colour = "white"),
        strip.text.y = element_text(size = 20, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"),
#        panel.grid = element_blank(),
#        panel.border = element_blank(),
#        panel.background = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        plot.title = element_text(size = 20, colour = "black", hjust = 0.5)) +
  ggtitle(bquote("Permutations =" ~
                 .(prettyNum(perms,
                             big.mark = ",",
                             trim = T)) ~ "sets of randomly positioned centromeric loci"))
bp_per <- bp_per +
  facet_grid(cols = vars(Accession), rows = vars(Region), scales = "fixed")
ggsave(paste0(plotDirAllMetrics,
              "CENATHILA_", flankName, "_regions_CEN180_AllMetrics_combined_bargraph_",
              paste0(chrName, collapse = "_"), "_",
              perms, "perms_",
              paste0(unique(as.vector(combined$Region)), collapse = "_"),
              ".pdf"),
       plot = bp_per,
       width = 10*length(unique(as.vector(combined$Accession))), height = 8, limitsize = F)

# All accessions plot (excluding per-accession plots)
all_accessions <- combined[grep(" accessions", combined$Accession),] 
bp_all <- ggplot(data = all_accessions,
                 mapping = aes(x = Family,
                               y = log2obsexp,
                               fill = Metric)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
#  scale_fill_discrete(name = "") +
  scale_fill_brewer(name = "",
                    palette = "Dark2") +
  geom_point(mapping = aes(x = Family,
                           y = log2alpha),
             position = position_dodge(0.9),
             shape = "-", colour  = "grey70", size = 12) +
  xlab(bquote(italic("CENATHILA") ~ .(flankNamePlot) ~ "flanking regions")) +
  ylab(bquote("Log"[2] * "(observed/expected)" ~ italic("CEN178") ~ "metric")) +
#  scale_y_continuous(limits = c(-4.0, 4.0)) +
  scale_x_discrete(position = "bottom") +
  guides(fill = guide_legend(direction = "vertical",
                             label.position = "right",
                             label.theme = element_text(size = 16, hjust = 0, vjust = 0.5, angle = 0),
                             nrow = length(unique(all_accessions$Family)),
                             byrow = TRUE)) +

  theme_bw() +
  theme(axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks.y = element_line(size = 0.5, colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black", hjust = 0.5, vjust = 0.5, angle = 90),
        axis.title.y = element_text(size = 20, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 20, colour = "black", hjust = 0.5, vjust = 0.5, angle = 45),
        axis.title.x = element_text(size = 20, colour = "black", margin = margin(t = 20)),
        strip.text.x = element_text(size = 20, colour = "white"),
        strip.text.y = element_text(size = 20, colour = "white"),
        strip.background = element_rect(fill = "black", colour = "black"),
#        panel.grid = element_blank(),
#        panel.border = element_blank(),
#        panel.background = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        plot.title = element_text(size = 20, colour = "black", hjust = 0.5)) +
  ggtitle(bquote("Permutations =" ~
                 .(prettyNum(perms,
                             big.mark = ",",
                             trim = T)) ~ "sets of randomly positioned centromeric loci"))
bp_all <- bp_all +
  facet_grid(cols = vars(Accession), scales = "fixed")
ggsave(paste0(plotDirAllMetrics,
              "CENATHILA_", flankName, "_regions_CEN180_AllMetrics_combined_all_accessions_bargraph_",
              paste0(chrName, collapse = "_"), "_",
              perms, "perms_",
              paste0(unique(as.vector(all_accessions$Region)), collapse = "_"),
              ".pdf"),
       plot = bp_all,
       width = 14*length(unique(as.vector(all_accessions$Accession))), height = 8, limitsize = F)

print("warnings 1")
print(warnings())

## Shutdown the cluster and quit
##stopCluster(cl) # use if cl made with makeCluster() (i.e., using doFuture package)
#closeCluster(cl) # use if cl made with startMPIcluster() (i.e., using doMPI package - faster than doFuture)
#mpi.quit()
#
#print("warnings 2")
#print(warnings())
