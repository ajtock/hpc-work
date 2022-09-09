#!/usr/bin/env Rscript

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 09/09/2022

# Plot putative crossover events detected in Col-0/Ler-0
# hybrid ONT reads


# Usage:
# ./circlize_co_nco.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 10000 Col-0.ragtag_scaffolds Ler-0_110x.ragtag_scaffolds co 090922

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#genomeBinSize <- 10000
#acc1 <- "Col-0.ragtag_scaffolds"
#acc2 <- "Ler-0_110x.ragtag_scaffolds"
#recombType <- "co"
#date <- "090922"

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
genomeBinSize <- as.integer(args[2])
acc1 <- args[3]
acc2 <- args[4]
recombType <- args[5]
date <- as.character(args[6])

if(floor(log10(genomeBinSize)) + 1 < 4) {
  genomeBinName <- paste0(genomeBinSize, "bp")
} else if(floor(log10(genomeBinSize)) + 1 >= 4 &
          floor(log10(genomeBinSize)) + 1 <= 6) {
  genomeBinName <- paste0(genomeBinSize/1e3, "kb")
} else if(floor(log10(genomeBinSize)) + 1 >= 7) {
  genomeBinName <- paste0(genomeBinSize/1e6, "Mb")
}

options(stringsAsFactors = F)
options(scipen=999)
library(parallel)
library(GenomicRanges)
library(circlize)
library(ComplexHeatmap)
library(gridBase)
library(viridis)

plotDir <- paste0("plots/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# CEN coordinates
CEN <- read.csv(paste0("/home/ajt200/rds/hpc-work/pancentromere/centromeric_coordinates/",
                       "centromere_manual_EDTA4_fa.csv"),
                header = T)
CEN$fasta.name <- gsub(".fa", "", CEN$fasta.name)

# Genomic definitions

#acc1
acc1_fai <- read.table(paste0("index/", acc1, ".fa.fai"), header = F)
acc1_chrs <- acc1_fai[which(acc1_fai$V1 %in% chrName),]$V1
acc1_chrLens <- acc1_fai[which(acc1_fai$V1 %in% chrName),]$V2

acc1_CEN <- CEN[grep(acc1, CEN$fasta.name),]
acc1_CEN <- acc1_CEN[,which(colnames(acc1_CEN) %in% c("chr", "start", "end"))]
acc1_CEN_new <- data.frame()
for(i in 1:length(acc1_chrs)) {
  acc1_CEN_chr <- acc1_CEN[which(acc1_CEN$chr == acc1_chrs[i]),]
  if(nrow(acc1_CEN_chr) > 1) {
    acc1_CEN_chr <- data.frame(chr = acc1_CEN_chr$chr[1],
                               start = acc1_CEN_chr$start[1],
                               end = acc1_CEN_chr$end[nrow(acc1_CEN_chr)])
  }
  acc1_CEN_new <- rbind(acc1_CEN_new, acc1_CEN_chr)
}
acc1_CEN <- acc1_CEN_new
acc1_CENstart <- acc1_CEN$start
acc1_CENend <- acc1_CEN$end

#acc2
acc2_fai <- read.table(paste0("index/", acc2, ".fa.fai"), header = F)
acc2_chrs <- acc2_fai[which(acc2_fai$V1 %in% chrName),]$V1
acc2_chrLens <- acc2_fai[which(acc2_fai$V1 %in% chrName),]$V2

acc2_CEN <- CEN[grep(acc2, CEN$fasta.name),]
acc2_CEN <- acc2_CEN[,which(colnames(acc2_CEN) %in% c("chr", "start", "end"))]
acc2_CEN_new <- data.frame()
for(i in 1:length(acc2_chrs)) {
  acc2_CEN_chr <- acc2_CEN[which(acc2_CEN$chr == acc2_chrs[i]),]
  if(nrow(acc2_CEN_chr) > 1) {
    acc2_CEN_chr <- data.frame(chr = acc2_CEN_chr$chr[1],
                               start = acc2_CEN_chr$start[1],
                               end = acc2_CEN_chr$end[nrow(acc2_CEN_chr)])
  }
  acc2_CEN_new <- rbind(acc2_CEN_new, acc2_CEN_chr)
}
acc2_CEN <- acc2_CEN_new
acc2_CENstart <- acc2_CEN$start
acc2_CENend <- acc2_CEN$end

# Accession names
acc1_name <- strsplit( strsplit(acc1, split = "\\.")[[1]][1],
                       split = "_")[[1]][1]
acc2_name <- strsplit( strsplit(acc2, split = "\\.")[[1]][1],
                       split = "_")[[1]][1]

# Directories containing read segment alignment files
acc1_indir <- paste0("segments/", acc1_name, "/", recombType, "/")
acc2_indir <- paste0("segments/", acc2_name, "/", recombType, "/")

# Load read segments
load_pafs <- function(acc_name=acc1_name, indir=acc1_indir suffix="_wm_ont.paf") {

system(paste0("ls -1 *", 
read.table(paste0(acc1_indir
system(

bed1 = generateRandomBed(nr = 100)
head(bed1)
bed1 = bed1[sample(nrow(bed1), 20),]
head(bed1)







# Genomic definitions
fai <- read.table("/home/ajt200/analysis/nanopore/t2t-col.20210610/t2t-col.20210610.fa.fai", header = F)                                      
chrs <- fai$V1[which(fai$V1 %in% chrName)]
chrLens <- fai$V2[which(fai$V1 %in% chrName)]
regionGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = 1,
                                     end = chrLens),
                    strand = "*")
CENstart <- c(14840750, 3724530, 13597090, 4203495, 11783990)[which(fai$V1 %in% chrName)]                                                     
CENend <- c(17558182, 5946091, 15733029, 6977107, 14551874)[which(fai$V1 %in% chrName)]                                                       
CENGR <- GRanges(seqnames = chrs,
                 ranges = IRanges(start = CENstart,
                                  end = CENend),
                 strand = "*")



