#!/usr/bin/env Rscript

# Convert CEN180 sequence coordinates identified in acc
# from CSV (1-based start coordinates) into BED format (0-based start coordinates)
# Convert coordinates corresponding to the intervening sequences between CEN180
# sequences from CSV into BED format

# Usage:
# conda activate R-4.0.0
# ./CEN180_acc_CSVtoBED.R Cvi-0.ragtag_scaffolds 'Chr1,Chr2,Chr3,Chr4,Chr5'
# conda deactivate

#acc <- "Cvi-0.ragtag_scaffolds"
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
acc <- args[1]
chrName <- unlist(strsplit(args[2],
                           split = ","))

options(stringsAsFactors = F)
library(GenomicRanges)

outDir <- paste0(acc, "/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/analysis/pancentromere/assemblies/", acc, ".fa.fai"),
                  header = F)
chrs <- fai$V1[which(fai$V1 %in% chrName)]
chrLens <- fai$V2[which(fai$V1 %in% chrName)]
regionGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = 1,
                                     end = chrLens),
                    strand = "*")

# Load table of centromeric coordinates
allaccsCEN <- read.csv(paste0("/home/ajt200/analysis/pancentromere/centromeric_coordinates/",
                              "centromeric.coordinates.300K.prune.robin.03.02.22.csv"),
                       header = T)
CEN <- allaccsCEN[grep(acc, allaccsCEN$region.acc),]
CENGR <- GRanges(seqnames = CEN$chromosome,
                 ranges = IRanges(start = CEN$start,
                                  end = CEN$end),
                 strand = "*")
CENGR <- CENGR[which(seqnames(CENGR) %in% chrName)]

# Load table of CEN180 sequence coordinates
# with weighted SNVs vs consensus

tab_list <- lapply(1:length(chrName), function(y) {
  read.csv(paste0("/home/ajt200/analysis/pancentromere/annotation/CEN180/repeats/",
                  "cen180.consensus.repetitiveness", acc, ".fa.", chrName[y], ".csv"),
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
tab$fasta.file.name <- gsub("\\.fa.+", "", tab$assembly)
tab <- tab[
           with( tab, order(chr, start, end) ),
          ]

CEN180GR <- GRanges(seqnames = tab$chr,
                    ranges = IRanges(start = tab$start,
                                     end = tab$end),
                    strand = tab$strand,
                    index = as.integer(rownames(tab)),
                    weighted.consensus.score = tab$weighted.consensus.score,
                    HORlengthsSum = tab$HORlengthsSum,
                    HORcount = tab$HORcount,
                    edit.distance = tab$edit.distance)

# Sort GRanges object (by chromosome, strand, start coordinate and, finally, end coordinate)
# Sorting data.frame by multiple columns (including that correspodning to strand)
# would be more complicated
# Necessary to include sorting by strand because intervening CEN180 sequences on the
# opposite strand would otherwise result in non-detection of tandem repeats on the same strand
CEN180GR <- sort(CEN180GR)
CEN180GR <- CEN180GR[seqnames(CEN180GR) %in% chrName]
CEN180_bed <- data.frame(chr = as.character(seqnames(CEN180GR)),
                         start = start(CEN180GR)-1,
                         end = end(CEN180GR),
                         name = CEN180GR$index,
                         score = CEN180GR$weighted.consensus.score,
                         strand = as.character(strand(CEN180GR)),
                         HORlengthsSum = CEN180GR$HORlengthsSum,
                         HORcount = CEN180GR$HORcount,
                         edit.distance = CEN180GR$edit.distance)
write.table(CEN180_bed,
            file = paste0(outDir, "CEN180_in_", acc, "_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Define function to select randomly positioned loci of the same
# width distribution as CEN180_bed
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLocGR contains the same number of loci per chromosome as CEN180GR
ranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  CEN180ChrGR <- CEN180GR[seqnames(CEN180GR) == chrs[i]]
  regionChrGR <- regionGR[seqnames(regionGR) == chrs[i]]
  # Contract regionChrGR so that random loci and 2-kb flanking regions
  # do not extend beyond chromosome ends
  end(regionChrGR) <- end(regionChrGR)-max(width(CEN180ChrGR))-2000
  start(regionChrGR) <- start(regionChrGR)+2000
  # Define seed so that random selections are reproducible
  set.seed(76492749)
  ranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(regionChrGR), function(x) {           
                                                             start(regionChrGR[x]) : end(regionChrGR[x])          
                                                           })),
                                      n = length(CEN180ChrGR))
  ranLocChrGR <- GRanges(seqnames = chrs[i],
                         ranges = IRanges(start = ranLocChrStart,
                                          width = width(CEN180ChrGR)),
                         strand = strand(CEN180ChrGR))
  ranLocGR <- append(ranLocGR, ranLocChrGR)
}
stopifnot(identical(width(ranLocGR), width(CEN180GR)))
stopifnot(identical(as.character(seqnames(ranLocGR)), as.character(seqnames(CEN180GR))))
stopifnot(identical(strand(ranLocGR), strand(CEN180GR)))
ranLoc_bed <- data.frame(chr = as.character(seqnames(ranLocGR)),
                         start = start(ranLocGR)-1,
                         end = end(ranLocGR),
                         name = 1:length(ranLocGR),
                         score = "NA",
                         strand = strand(ranLocGR))
write.table(ranLoc_bed,
            file = paste0(outDir, "CEN180_in_", acc, "_",
                          paste0(chrName, collapse = "_"), "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# CENranLocGR contains the same number of loci per chromosome as CEN180GR
CENranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  CEN180ChrGR <- CEN180GR[seqnames(CEN180GR) == chrs[i]]
  CENChrGR <- CENGR[seqnames(CENGR) == chrs[i]]
  # Contract CENChrGR so that CENrandom loci and 2-kb flanking regions
  # do not extend beyond chromosome ends
  end(CENChrGR) <- end(CENChrGR)-max(width(CEN180ChrGR))-2000
  start(CENChrGR) <- start(CENChrGR)+2000
  # Define seed so that random selections are reproducible
  set.seed(76492749)
  CENranLocChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(CENChrGR), function(x) {           
                                                                start(CENChrGR[x]) : end(CENChrGR[x])          
                                                              })),
                                         n = length(CEN180ChrGR))
  CENranLocChrGR <- GRanges(seqnames = chrs[i],
                            ranges = IRanges(start = CENranLocChrStart,
                                             width = width(CEN180ChrGR)),
                            strand = strand(CEN180ChrGR))
  CENranLocGR <- append(CENranLocGR, CENranLocChrGR)
}
stopifnot(identical(width(CENranLocGR), width(CEN180GR)))
stopifnot(identical(as.character(seqnames(CENranLocGR)), as.character(seqnames(CEN180GR))))
stopifnot(identical(strand(CENranLocGR), strand(CEN180GR)))
CENranLoc_bed <- data.frame(chr = as.character(seqnames(CENranLocGR)),
                            start = start(CENranLocGR)-1,
                            end = end(CENranLocGR),
                            name = 1:length(CENranLocGR),
                            score = "NA",
                            strand = strand(CENranLocGR))
write.table(CENranLoc_bed,
            file = paste0(outDir, "CEN180_in_", acc, "_",
                          paste0(chrName, collapse = "_"), "_CENrandomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
