#!/applications/R/R-4.0.0/bin/Rscript

# Convert centromeric and non-centromeric Athila LTR element coordinates identified in t2t-col.20210610
# from TAB (1-based start coordinates) into BED format (0-based start coordinates)
# Strand and additional information for these elements was determined by Alexandros Bousios (University of Sussex)

# Usage:
# ./58Atha_ATHILA_soloLTR_AlexBousios_TABtoBED_v110222.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))

options(stringsAsFactors = F)
library(GenomicRanges)
library(scales)
# From https://github.com/dmarcelinobr/ggdecor/blob/master/R/color.R :
rich10 <- function() {manual_pal(values = c("#000041","#0000A9","#0049FF","#00A4DE","#03E070","#5DFC21","#F6F905","#FFD701","#FF9500","#FF3300"))}
scale_colour_rich10 <- function(...) { discrete_scale("colour", "rich10", rich10(), ...) }
scale_fill_rich10 <- function(...) { discrete_scale("fill", "rich10", rich10(), ...) }
rich12 <- function() {manual_pal(values = c("#000040","#000093","#0020E9","#0076FF","#00B8C2","#04E466","#49FB25","#E7FD09","#FEEA02","#FFC200","#FF8500","#FF3300"))}
scale_colour_rich12 <- function(...) { discrete_scale("colour", "rich12", rich12(), ...) }
scale_fill_rich12 <- function(...) { discrete_scale("fill", "rich12", rich12(), ...) }

# Load table of centromeric coordinates
CEN <- read.csv(paste0("/home/ajt200/analysis/pancentromere/centromeric_coordinates/",
                       "centromeric.coordinates.300K.prune.robin.03.02.22.csv"),
                header = T)

# Load table of ATHILA sequence coordinates
tab <- read.table("58Atha_ATHILA_soloLTR_AlexBousios_v110222.tsv",
                  header = T)
print(dim(tab))
tab <- tab[tab$chromosome %in% chrName,]
print(dim(tab))
tab$family_FL[which(tab$family_FL == "ATHILA2_check")] <- "ATHILA2"
tab$acc <- gsub("Atha_", "", tab$species)
tab$acc_trunc <- tab$acc
tab$species <- gsub("_.+", "", tab$species)
tab$family_FL <- toupper(tab$family_FL)

acc_full <- system("ls /home/ajt200/analysis/pancentromere/assemblies/*.fa", intern = T)
acc_full <- gsub("/home/ajt200/analysis/pancentromere/assemblies/", "", acc_full)
acc_full <- gsub(".fa", "", acc_full)
acc_full_trunc <- gsub("\\..+", "", acc_full)
acc_full_trunc <- gsub("_TAIR10", "", acc_full_trunc)
acc_full <- acc_full[which(acc_full_trunc %in% unique(tab$acc))]
acc_full_trunc <- acc_full_trunc[which(acc_full_trunc %in% unique(tab$acc))]

for(acc_trunc in unique(tab$acc)) {
  print(acc_trunc)
  acc_full_match <- acc_full[which(acc_full_trunc == acc_trunc)]
  stopifnot(length(acc_full_match) == 1)
  tab[tab$acc == acc_trunc,]$acc <- acc_full_match
} 

# Sanity check to make sure newly added full accession names match
# Alex's truncated accession names in each row of tab
grepl_acc_trunc_acc <- 0
for(i in 1:nrow(tab)) {
  grepl_acc_trunc_acc <- grepl_acc_trunc_acc + grepl(tab$acc_trunc[i], tab$acc[i])
}
stopifnot(grepl_acc_trunc_acc == nrow(tab))

ATHILA <- tab[tab$quality == "intact",]
soloLTR <- tab[tab$quality == "solo",]
stopifnot(nrow(ATHILA) + nrow(soloLTR) == nrow(tab))

ATHILA <- data.frame(species = ATHILA$species,
                     acc = ATHILA$acc,
                     chr = ATHILA$chromosome,
                     start = ATHILA$genome_left_coord_FL,
                     end = ATHILA$genome_right_coord_FL,
                     strand = ATHILA$direction,
                     phylo = ATHILA$family_FL,
                     TE_ID = ATHILA$TE_ID)

ATHILA <- ATHILA[
                 with( ATHILA, order(species, acc, chr, start, end) ),
                ]

ATHILA_GR <- GRanges(seqnames = ATHILA$chr,
                     ranges = IRanges(start = ATHILA$start, end = ATHILA$end),
                     strand = ATHILA$strand,
                     acc = ATHILA$acc,
                     phylo = ATHILA$phylo,
                     TE_ID = ATHILA$TE_ID)
ATHILA_GR <- unique(ATHILA_GR)

ATHILA_BED <- data.frame(chr = as.character(seqnames(ATHILA_GR)),
                         start = start(ATHILA_GR)-1,
                         end = end(ATHILA_GR),
                         name = ATHILA_GR$acc,
                         score = ATHILA_GR$phylo,
                         strand = as.character(strand(ATHILA_GR)))


soloLTR <- data.frame(species = soloLTR$species,
                      acc = soloLTR$acc,
                      chr = soloLTR$chromosome,
                      start = soloLTR$genome_left_coord_FL,
                      end = soloLTR$genome_right_coord_FL,
                      strand = soloLTR$direction,
                      phylo = soloLTR$family_FL,
                      TE_ID = soloLTR$TE_ID)

soloLTR <- soloLTR[
                   with( soloLTR, order(species, acc, chr, start, end) ),
                  ]

soloLTR_GR <- GRanges(seqnames = soloLTR$chr,
                      ranges = IRanges(start = soloLTR$start, end = soloLTR$end),
                      strand = soloLTR$strand,
                      acc = soloLTR$acc,
                      phylo = soloLTR$phylo,
                      TE_ID = soloLTR$TE_ID)
soloLTR_GR <- unique(soloLTR_GR)

soloLTR_BED <- data.frame(chr = as.character(seqnames(soloLTR_GR)),
                          start = start(soloLTR_GR)-1,
                          end = end(soloLTR_GR),
                          name = soloLTR_GR$acc,
                          score = soloLTR_GR$phylo,
                          strand = as.character(strand(soloLTR_GR)))

ATHILADir <- "ATHILA/"
soloLTRDir <- "soloLTR/"
system(paste0("[ -d ", ATHILADir, " ] || mkdir -p ", ATHILADir))
system(paste0("[ -d ", soloLTRDir, " ] || mkdir -p ", soloLTRDir))
allDir_ATHILA <- paste0(ATHILADir, "58Atha/")
allDir_soloLTR <- paste0(soloLTRDir, "58Atha/")
system(paste0("[ -d ", allDir_ATHILA, " ] || mkdir -p ", allDir_ATHILA))
system(paste0("[ -d ", allDir_soloLTR, " ] || mkdir -p ", allDir_soloLTR))
acc <- unique(tab$acc)
accDir_ATHILA <- sapply(seq_along(acc), function(x) {
  paste0(ATHILADir, acc[x], "/")
})
sapply(seq_along(accDir_ATHILA), function(x) {
  system(paste0("[ -d ", accDir_ATHILA[x], " ] || mkdir -p ", accDir_ATHILA[x])) 
})
accDir_soloLTR <- sapply(seq_along(acc), function(x) {
  paste0(soloLTRDir, acc[x], "/")
})
sapply(seq_along(accDir_soloLTR), function(x) {
  system(paste0("[ -d ", accDir_soloLTR[x], " ] || mkdir -p ", accDir_soloLTR[x])) 
})

write.table(ATHILA_BED,
            file = paste0(allDir_ATHILA, "ATHILA_in_58Atha_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
write.table(soloLTR_BED,
            file = paste0(allDir_soloLTR, "soloLTR_in_58Atha_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Write BED with colour-coded family information, for use with pyGenomeTracks (BED specified in *.ini file)
ATHILA_phylo <- sort(unique(ATHILA_BED$score))
print(ATHILA_phylo)
# [1] "ATHILA0"  "ATHILA1"  "ATHILA2"  "ATHILA3"  "ATHILA4"  "ATHILA4A"
# [7] "ATHILA4C" "ATHILA5"  "ATHILA6A" "ATHILA6B" "ATHILA7A" "ATHILA8A"

ATHILA_BED_colophylo <- data.frame(ATHILA_BED,
                                   thickStart = as.integer(0),
                                   thickEnd = as.integer(0),
                                   itemRgb = ".")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA0",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[12])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA1",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[11])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA2",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[10])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA3",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[9])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA4",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[8])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA4A",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[7])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA4C",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[6])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA5",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[5])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA6A",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[4])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA6B",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[3])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA7A",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[2])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA8A",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[1])), collapse = ",")
ATHILA_BED_colophylo$score <- as.integer(0)
write.table(ATHILA_BED_colophylo,
            file = paste0(allDir_ATHILA, "ATHILA_in_58Atha_colophylo_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)


# Define function to select randomly positioned loci of the same
# width distribution as CENgapAllAthila_bed
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)
  

# Make separate files for each accession and chromosome,
# and generate random loci equivalent to acc_CENATHILA and acc_nonCENATHILA
for(i in 1:length(acc)) {
  print(acc[i])
  acc_chrs <- read.table(paste0("/home/ajt200/analysis/pancentromere/assemblies/",
                                acc[i], ".fa.fai"),
                         header = F)[,1]
  acc_chrs <- gsub("_RagTag_RagTag", "", acc_chrs)
  acc_chrs <- gsub("chr", "Chr", acc_chrs)
  acc_chrLens <- read.table(paste0("/home/ajt200/analysis/pancentromere/assemblies/",
                                   acc[i], ".fa.fai"),
                            header = F)[,2]
  acc_chrLens <- acc_chrLens[which(acc_chrs %in% chrName)]
  acc_chrs <- acc_chrs[which(acc_chrs %in% chrName)]
  print(acc_chrs)
  print(acc_chrLens)

  acc_CEN <- CEN[grep(acc[i], CEN$region.acc),]
  acc_CEN_GR <- GRanges(seqnames = acc_CEN$chromosome,
                        ranges = IRanges(start = acc_CEN$start,
                                         end = acc_CEN$end),
                        strand = "*",
                        acc = acc_CEN$region.acc)
  print(acc_CEN_GR)

  acc_ATHILA_GR <- ATHILA_GR[ATHILA_GR$acc == acc[i]]
  print(acc_ATHILA_GR)
 
  acc_CEN_acc_ATHILA_overlaps <- findOverlaps(query = acc_CEN_GR,
                                              subject = acc_ATHILA_GR,
                                              type = "any",
                                              select = "all",
                                              ignore.strand = TRUE)
  acc_CENATHILA_GR <- acc_ATHILA_GR[unique(subjectHits(acc_CEN_acc_ATHILA_overlaps))]
  acc_nonCENATHILA_GR <- acc_ATHILA_GR[-subjectHits(acc_CEN_acc_ATHILA_overlaps)]
  stopifnot(length(acc_CENATHILA_GR) + length(acc_nonCENATHILA_GR) == length(acc_ATHILA_GR))

  acc_CENATHILA_BED <- data.frame(chr = as.character(seqnames(acc_CENATHILA_GR)),
                                  start = start(acc_CENATHILA_GR)-1,
                                  end = end(acc_CENATHILA_GR),
                                  name = acc_CENATHILA_GR$TE_ID,
                                  score = acc_CENATHILA_GR$phylo,
                                  strand = strand(acc_CENATHILA_GR))
  write.table(acc_CENATHILA_BED,
              file = paste0(accDir_ATHILA[i], "CENATHILA_in_", acc[i], "_",
                            paste0(chrName, collapse = "_"), ".bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)

  acc_nonCENATHILA_BED <- data.frame(chr = as.character(seqnames(acc_nonCENATHILA_GR)),
                                     start = start(acc_nonCENATHILA_GR)-1,
                                     end = end(acc_nonCENATHILA_GR),
                                     name = acc_nonCENATHILA_GR$TE_ID,
                                     score = acc_nonCENATHILA_GR$phylo,
                                     strand = strand(acc_nonCENATHILA_GR))
  write.table(acc_nonCENATHILA_BED,
              file = paste0(accDir_ATHILA[i], "nonCENATHILA_in_", acc[i], "_",
                            paste0(chrName, collapse = "_"), ".bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)


  acc_soloLTR_GR <- soloLTR_GR[soloLTR_GR$acc == acc[i]]
  print(acc_soloLTR_GR)
 
  acc_CEN_acc_soloLTR_overlaps <- findOverlaps(query = acc_CEN_GR,
                                               subject = acc_soloLTR_GR,
                                               type = "any",
                                               select = "all",
                                               ignore.strand = TRUE)
  acc_CENsoloLTR_GR <- acc_soloLTR_GR[unique(subjectHits(acc_CEN_acc_soloLTR_overlaps))]
  acc_nonCENsoloLTR_GR <- acc_soloLTR_GR[-subjectHits(acc_CEN_acc_soloLTR_overlaps)]
  stopifnot(length(acc_CENsoloLTR_GR) + length(acc_nonCENsoloLTR_GR) == length(acc_soloLTR_GR))

  acc_CENsoloLTR_BED <- data.frame(chr = as.character(seqnames(acc_CENsoloLTR_GR)),
                                   start = start(acc_CENsoloLTR_GR)-1,
                                   end = end(acc_CENsoloLTR_GR),
                                   name = acc_CENsoloLTR_GR$TE_ID,
                                   score = acc_CENsoloLTR_GR$phylo,
                                   strand = strand(acc_CENsoloLTR_GR))
  write.table(acc_CENsoloLTR_BED,
              file = paste0(accDir_soloLTR[i], "CENsoloLTR_in_", acc[i], "_",
                            paste0(chrName, collapse = "_"), ".bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)

  acc_nonCENsoloLTR_BED <- data.frame(chr = as.character(seqnames(acc_nonCENsoloLTR_GR)),
                                      start = start(acc_nonCENsoloLTR_GR)-1,
                                      end = end(acc_nonCENsoloLTR_GR),
                                      name = acc_nonCENsoloLTR_GR$TE_ID,
                                      score = acc_nonCENsoloLTR_GR$phylo,
                                      strand = strand(acc_nonCENsoloLTR_GR))
  write.table(acc_nonCENsoloLTR_BED,
              file = paste0(accDir_soloLTR[i], "nonCENsoloLTR_in_", acc[i], "_",
                            paste0(chrName, collapse = "_"), ".bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)

  # Apply ranLocStartSelect() on a per-chromosome basis so that
  # acc_CENranLoc_GR contains the same number of loci per chromosome as acc_CENATHILA_GR
  acc_CENranLoc_GR <- GRanges()
  acc_nonCENranLoc_GR <- GRanges()
  for(j in 1:length(acc_chrs)) {
    print(acc_chrs[j])

    chr_acc_CEN_GR <- acc_CEN_GR[seqnames(acc_CEN_GR) == acc_chrs[j]]

    chr_acc_CENATHILA_GR <- acc_CENATHILA_GR[seqnames(acc_CENATHILA_GR) == acc_chrs[j]]
    if(length(chr_acc_CENATHILA_GR) > 0) {

      chr_acc_CENATHILA_BED <- data.frame(chr = as.character(seqnames(chr_acc_CENATHILA_GR)),
                                          start = start(chr_acc_CENATHILA_GR)-1,
                                          end = end(chr_acc_CENATHILA_GR),
                                          name = chr_acc_CENATHILA_GR$TE_ID,
                                          score = chr_acc_CENATHILA_GR$phylo,
                                          strand = strand(chr_acc_CENATHILA_GR))
      write.table(chr_acc_CENATHILA_BED,
                  file = paste0(accDir_ATHILA[i], "CENATHILA_in_", acc[i], "_",
                                acc_chrs[j], ".bed"),
                  quote = F, sep = "\t", row.names = F, col.names = F)

      ## Contract chr_acc_CEN_GR so that acc_CENranLoc_GR and 2-kb flanking regions
      ## do not extend beyond centromeric coordinates
      #end(chr_acc_CEN_GR) <- end(chr_acc_CEN_GR)-max(width(chr_acc_CENATHILA_GR))-2000
      #start(chr_acc_CEN_GR) <- start(chr_acc_CEN_GR)+2000

      # Define seed so that random selections are reproducible
      set.seed(76492749)
      chr_acc_CENranLoc_Start <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(chr_acc_CEN_GR), function(x) {
                                                                          start(chr_acc_CEN_GR[x]) : end(chr_acc_CEN_GR[x])
                                                                        })),
                                                   n = length(chr_acc_CENATHILA_GR))
      chr_acc_CENranLoc_GR <- GRanges(seqnames = acc_chrs[j],
                                      ranges = IRanges(start = chr_acc_CENranLoc_Start,
                                                       width = width(chr_acc_CENATHILA_GR)),
                                      strand = strand(chr_acc_CENATHILA_GR))
      acc_CENranLoc_GR <- append(acc_CENranLoc_GR, chr_acc_CENranLoc_GR)

      chr_acc_CENranLoc_BED <- data.frame(chr = as.character(seqnames(chr_acc_CENranLoc_GR)),
                                          start = start(chr_acc_CENranLoc_GR)-1,
                                          end = end(chr_acc_CENranLoc_GR),
                                          name = 1:length(chr_acc_CENranLoc_GR),
                                          score = chr_acc_CENATHILA_GR$phylo,
                                          strand = strand(chr_acc_CENranLoc_GR))
      write.table(chr_acc_CENranLoc_BED,
                  file = paste0(accDir_ATHILA[i], "CENATHILA_in_", acc[i], "_",
                                acc_chrs[j], "_CENrandomLoci.bed"),
                  quote = F, sep = "\t", row.names = F, col.names = F)
    }

    chr_acc_nonCENATHILA_GR <- acc_nonCENATHILA_GR[seqnames(acc_nonCENATHILA_GR) == acc_chrs[j]]
    if(length(chr_acc_nonCENATHILA_GR) > 0) {
      chr_acc_nonCENATHILA_BED <- data.frame(chr = as.character(seqnames(chr_acc_nonCENATHILA_GR)),
                                             start = start(chr_acc_nonCENATHILA_GR)-1,
                                             end = end(chr_acc_nonCENATHILA_GR),
                                             name = chr_acc_nonCENATHILA_GR$TE_ID,
                                             score = chr_acc_nonCENATHILA_GR$phylo,
                                             strand = strand(chr_acc_nonCENATHILA_GR))
      write.table(chr_acc_nonCENATHILA_BED,
                  file = paste0(accDir_ATHILA[i], "nonCENATHILA_in_", acc[i], "_",
                                acc_chrs[j], ".bed"),
                  quote = F, sep = "\t", row.names = F, col.names = F)

      # Define seed so that random selections are reproducible
      set.seed(76492749)
      chr_acc_nonCENranLoc_Start <- ranLocStartSelect(coordinates = unique(unlist(lapply(seq_along(chr_acc_CEN_GR), function(x) {
                                                                             c(1:acc_chrLens[j])[
                                                                               which(
                                                                                 !( 1:acc_chrLens[j] %in%
                                                                                    start(chr_acc_CEN_GR[x]) : end(chr_acc_CEN_GR[x])
                                                                                  )
                                                                               )
                                                                             ]
                                                                           }))),
                                                      n = length(chr_acc_nonCENATHILA_GR))
      chr_acc_nonCENranLoc_GR <- GRanges(seqnames = acc_chrs[j],
                                         ranges = IRanges(start = chr_acc_nonCENranLoc_Start,
                                                          width = width(chr_acc_nonCENATHILA_GR)),
                                         strand = strand(chr_acc_nonCENATHILA_GR))
      acc_nonCENranLoc_GR <- append(acc_nonCENranLoc_GR, chr_acc_nonCENranLoc_GR)

      chr_acc_nonCENranLoc_BED <- data.frame(chr = as.character(seqnames(chr_acc_nonCENranLoc_GR)),
                                             start = start(chr_acc_nonCENranLoc_GR)-1,
                                             end = end(chr_acc_nonCENranLoc_GR),
                                             name = 1:length(chr_acc_nonCENranLoc_GR),
                                             score = chr_acc_nonCENATHILA_GR$phylo,
                                             strand = strand(chr_acc_nonCENranLoc_GR))
      write.table(chr_acc_nonCENranLoc_BED,
                  file = paste0(accDir_ATHILA[i], "nonCENATHILA_in_", acc[i], "_",
                                acc_chrs[j], "_nonCENrandomLoci.bed"),
                  quote = F, sep = "\t", row.names = F, col.names = F)
    }

    chr_acc_CENsoloLTR_GR <- acc_CENsoloLTR_GR[seqnames(acc_CENsoloLTR_GR) == acc_chrs[j]]
    if(length(chr_acc_CENsoloLTR_GR) > 0) {
      chr_acc_CENsoloLTR_BED <- data.frame(chr = as.character(seqnames(chr_acc_CENsoloLTR_GR)),
                                          start = start(chr_acc_CENsoloLTR_GR)-1,
                                          end = end(chr_acc_CENsoloLTR_GR),
                                          name = chr_acc_CENsoloLTR_GR$TE_ID,
                                          score = chr_acc_CENsoloLTR_GR$phylo,
                                          strand = strand(chr_acc_CENsoloLTR_GR))
      write.table(chr_acc_CENsoloLTR_BED,
                  file = paste0(accDir_soloLTR[i], "CENsoloLTR_in_", acc[i], "_",
                                acc_chrs[j], ".bed"),
                  quote = F, sep = "\t", row.names = F, col.names = F)
    }

    chr_acc_nonCENsoloLTR_GR <- acc_nonCENsoloLTR_GR[seqnames(acc_nonCENsoloLTR_GR) == acc_chrs[j]]
    if(length(chr_acc_nonCENsoloLTR_GR) > 0) {
      chr_acc_nonCENsoloLTR_BED <- data.frame(chr = as.character(seqnames(chr_acc_nonCENsoloLTR_GR)),
                                              start = start(chr_acc_nonCENsoloLTR_GR)-1,
                                              end = end(chr_acc_nonCENsoloLTR_GR),
                                              name = chr_acc_nonCENsoloLTR_GR$TE_ID,
                                              score = chr_acc_nonCENsoloLTR_GR$phylo,
                                              strand = strand(chr_acc_nonCENsoloLTR_GR))
      write.table(chr_acc_nonCENsoloLTR_BED,
                  file = paste0(accDir_soloLTR[i], "nonCENsoloLTR_in_", acc[i], "_",
                                acc_chrs[j], ".bed"),
                  quote = F, sep = "\t", row.names = F, col.names = F)
    }

  }
  stopifnot(identical(width(acc_CENranLoc_GR), width(acc_CENATHILA_GR)))
  stopifnot(identical(as.character(seqnames(acc_CENranLoc_GR)), as.character(seqnames(acc_CENATHILA_GR))))
  stopifnot(identical(strand(acc_CENranLoc_GR), strand(acc_CENATHILA_GR)))
  acc_CENranLoc_BED <- data.frame(chr = as.character(seqnames(acc_CENranLoc_GR)),
                                  start = start(acc_CENranLoc_GR)-1,
                                  end = end(acc_CENranLoc_GR),
                                  name = 1:length(acc_CENranLoc_GR),
                                  score = acc_CENATHILA_GR$phylo,
                                  strand = strand(acc_CENranLoc_GR))
  write.table(acc_CENranLoc_BED,
              file = paste0(accDir_ATHILA[i], "CENATHILA_in_", acc[i], "_",
                            paste0(chrName, collapse = "_"), "_CENrandomLoci.bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)

  stopifnot(identical(width(acc_nonCENranLoc_GR), width(acc_nonCENATHILA_GR)))
  stopifnot(identical(as.character(seqnames(acc_nonCENranLoc_GR)), as.character(seqnames(acc_nonCENATHILA_GR))))
  stopifnot(identical(strand(acc_nonCENranLoc_GR), strand(acc_nonCENATHILA_GR)))
  acc_nonCENranLoc_BED <- data.frame(chr = as.character(seqnames(acc_nonCENranLoc_GR)),
                                     start = start(acc_nonCENranLoc_GR)-1,
                                     end = end(acc_nonCENranLoc_GR),
                                     name = 1:length(acc_nonCENranLoc_GR),
                                     score = acc_nonCENATHILA_GR$phylo,
                                     strand = strand(acc_nonCENranLoc_GR))
  write.table(acc_nonCENranLoc_BED,
              file = paste0(accDir_ATHILA[i], "nonCENATHILA_in_", acc[i], "_",
                            paste0(chrName, collapse = "_"), "_nonCENrandomLoci.bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)

}
print(warnings())

