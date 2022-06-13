#!/applications/R/R-4.0.0/bin/Rscript

# Convert centromeric and non-centromeric Athila LTR element coordinates identified in t2t-col.20210610
# from TAB (1-based start coordinates) into BED format (0-based start coordinates)
# Strand and additional information for these elements was determined by Alexandros Bousios (University of Sussex)

# Usage:
# ./66Atha_ATHILA_soloLTR_AlexBousios_TABtoBED_v250522.R 'Chr1,Chr2,Chr3,Chr4,Chr5'

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))

options(stringsAsFactors = F)
library(GenomicRanges)
library(scales)
library(doParallel)
library(doRNG)
library(doFuture)
registerDoFuture()
plan(multicore, workers = detectCores())
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

# From https://github.com/dmarcelinobr/ggdecor/blob/master/R/color.R :
rich10 <- function() {manual_pal(values = c("#000041","#0000A9","#0049FF","#00A4DE","#03E070","#5DFC21","#F6F905","#FFD701","#FF9500","#FF3300"))}
scale_colour_rich10 <- function(...) { discrete_scale("colour", "rich10", rich10(), ...) }
scale_fill_rich10 <- function(...) { discrete_scale("fill", "rich10", rich10(), ...) }
rich12 <- function() {manual_pal(values = c("#000040","#000093","#0020E9","#0076FF","#00B8C2","#04E466","#49FB25","#E7FD09","#FEEA02","#FFC200","#FF8500","#FF3300"))}
scale_colour_rich12 <- function(...) { discrete_scale("colour", "rich12", rich12(), ...) }
scale_fill_rich12 <- function(...) { discrete_scale("fill", "rich12", rich12(), ...) }

# Load table of centromeric coordinates
CEN <- read.csv(paste0("/home/ajt200/analysis/pancentromere/centromeric_coordinates/",
                       "centromere_manual_EDTA4_fa.csv"),
                header = T)
CEN$fasta.name <- gsub(".fa", "", CEN$fasta.name)

# Load table of ATHILA sequence coordinates
tab <- read.table("66Atha_fulllength+soloLTRs.matrix+phylo+hmm+cengap+solofam_new_v240522.tsv",
                  header = T)
print(dim(tab))
tab <- tab[which(tab$chromosome %in% chrName),]
print(dim(tab))
print(unique(tab$FL_fam))
tab$FL_fam[which(tab$FL_fam == "ATHILA2_check")] <- "ATHILA2"
tab$FL_fam[which(tab$FL_fam == "ATHILA1_maybe")] <- "ATHILA1"
tab$FL_fam[which(tab$FL_fam == "ATHILA5_maybe")] <- "ATHILA5"
tab$FL_fam[which(tab$FL_fam == "ATHILA6a_maybe")] <- "ATHILA6a"
print(unique(tab$FL_fam))
print(unique(tab$combo_fam))
print(unique(tab$quality))
print(unique(tab$centromere))
tab$accession <- gsub(".fa", "", tab$accession_fasta)
tab$accession_trunc <- gsub("Atha_", "", tab$species)
tab$species <- gsub("_.+", "", tab$species)
tab <- tab[which(tab$species == "Atha"),]
#tab$FL_fam <- toupper(tab$FL_fam)

acc_full <- system("ls /home/ajt200/analysis/pancentromere/assemblies/*.fa", intern = T)
acc_full <- gsub("/home/ajt200/analysis/pancentromere/assemblies/", "", acc_full)
acc_full <- gsub(".fa", "", acc_full)
acc_full_trunc <- gsub("\\..+", "", acc_full)
acc_full <- acc_full[which(acc_full %in% unique(tab$accession))]
acc_full_trunc <- acc_full_trunc[which(acc_full_trunc %in% unique(tab$accession_trunc))]

# Sanity check to ensure correct number of accessions
stopifnot(length(unique(acc_full)) == length(unique(tab$accession)))
for(acc in unique(tab$accession)) {
  print(acc)
  acc_full_match <- acc_full[which(acc_full == acc)]
  stopifnot(length(acc_full_match) == 1)
} 

# Sanity check to make sure full accession names match
# Alex's truncated accession names in each row of tab
grepl_acc_trunc_acc <- 0
for(i in 1:nrow(tab)) {
  grepl_acc_trunc_acc <- grepl_acc_trunc_acc + grepl(tab$accession_trunc[i], tab$accession[i])
}
stopifnot(grepl_acc_trunc_acc == nrow(tab))


ATHILA <- tab[which(tab$quality == "intact"),]
soloLTR <- tab[which(tab$quality == "solo"),]
stopifnot(nrow(ATHILA) + nrow(soloLTR) == nrow(tab))


ATHILA <- data.frame(cen = ATHILA$centromere,
                     acc = ATHILA$accession,
                     chr = ATHILA$chromosome,
                     start = ATHILA$genome_left_coord_FL,
                     end = ATHILA$genome_right_coord_FL,
                     strand = ATHILA$direction,
                     phylo = ATHILA$FL_fam,
                     TE_ID = ATHILA$TE_ID)

ATHILA <- ATHILA[
                 with( ATHILA, order(acc, chr, start, end) ),
                ]

ATHILA_GR <- GRanges(seqnames = ATHILA$chr,
                     ranges = IRanges(start = ATHILA$start, end = ATHILA$end),
                     strand = ATHILA$strand,
                     acc = ATHILA$acc,
                     phylo = ATHILA$phylo,
                     TE_ID = ATHILA$TE_ID,
                     cen = ATHILA$cen)
ATHILA_GR <- unique(ATHILA_GR)

ATHILA_BED <- data.frame(chr = as.character(seqnames(ATHILA_GR)),
                         start = start(ATHILA_GR)-1,
                         end = end(ATHILA_GR),
                         name = ATHILA_GR$TE_ID,
                         score = ATHILA_GR$phylo,
                         strand = as.character(strand(ATHILA_GR)),
                         cen = ATHILA_GR$cen)

CENATHILA <- ATHILA[which(ATHILA$cen == "in"),]
CENATHILA_GR <- ATHILA_GR[which(ATHILA_GR$cen == "in")]
CENATHILA_BED <- ATHILA_BED[which(ATHILA_BED$cen == "in"),]
nonCENATHILA <- ATHILA[which(ATHILA$cen == "out"),]
nonCENATHILA_GR <- ATHILA_GR[which(ATHILA_GR$cen == "out")]
nonCENATHILA_BED <- ATHILA_BED[which(ATHILA_BED$cen == "out"),]

ATHILA_BED <- ATHILA_BED[,-7]
CENATHILA_BED <- CENATHILA_BED[,-7]
nonCENATHILA_BED <- nonCENATHILA_BED[,-7]


soloLTR <- data.frame(cen = soloLTR$centromere,
                      acc = soloLTR$accession,
                      chr = soloLTR$chromosome,
                      start = soloLTR$genome_left_coord_FL,
                      end = soloLTR$genome_right_coord_FL,
                      strand = soloLTR$direction,
                      phylo = soloLTR$combo_fam,
                      TE_ID = soloLTR$TE_ID)

soloLTR <- soloLTR[
                   with( soloLTR, order(acc, chr, start, end) ),
                  ]

soloLTR_GR <- GRanges(seqnames = soloLTR$chr,
                      ranges = IRanges(start = soloLTR$start, end = soloLTR$end),
                      strand = soloLTR$strand,
                      acc = soloLTR$acc,
                      phylo = soloLTR$phylo,
                      TE_ID = soloLTR$TE_ID,
                      cen = soloLTR$cen)
soloLTR_GR <- unique(soloLTR_GR)

soloLTR_BED <- data.frame(chr = as.character(seqnames(soloLTR_GR)),
                          start = start(soloLTR_GR)-1,
                          end = end(soloLTR_GR),
                          name = soloLTR_GR$TE_ID,
                          score = soloLTR_GR$phylo,
                          strand = as.character(strand(soloLTR_GR)),
                          cen = soloLTR_GR$cen)

CENsoloLTR <- soloLTR[which(soloLTR$cen == "in"),]
CENsoloLTR_GR <- soloLTR_GR[which(soloLTR_GR$cen == "in")]
CENsoloLTR_BED <- soloLTR_BED[which(soloLTR_BED$cen == "in"),]
nonCENsoloLTR <- soloLTR[which(soloLTR$cen == "out"),]
nonCENsoloLTR_GR <- soloLTR_GR[which(soloLTR_GR$cen == "out")]
nonCENsoloLTR_BED <- soloLTR_BED[which(soloLTR_BED$cen == "out"),]

soloLTR_BED <- soloLTR_BED[,-7]
CENsoloLTR_BED <- CENsoloLTR_BED[,-7]
nonCENsoloLTR_BED <- nonCENsoloLTR_BED[,-7]


ATHILADir <- "ATHILA/"
soloLTRDir <- "soloLTR/"
system(paste0("[ -d ", ATHILADir, " ] || mkdir -p ", ATHILADir))
system(paste0("[ -d ", soloLTRDir, " ] || mkdir -p ", soloLTRDir))
allDir_ATHILA <- paste0(ATHILADir, "66Atha/")
allDir_soloLTR <- paste0(soloLTRDir, "66Atha/")
system(paste0("[ -d ", allDir_ATHILA, " ] || mkdir -p ", allDir_ATHILA))
system(paste0("[ -d ", allDir_soloLTR, " ] || mkdir -p ", allDir_soloLTR))
acc <- unique(tab$accession)
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
            file = paste0(allDir_ATHILA, "ATHILA_in_66Atha_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
write.table(soloLTR_BED,
            file = paste0(allDir_soloLTR, "soloLTR_in_66Atha_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Write BED with colour-coded family information, for use with pyGenomeTracks (BED specified in *.ini file)
ATHILA_phylo <- sort(unique(ATHILA_BED$score))
print(ATHILA_phylo)
# [1] "ATHILA_general" "ATHILA0"        "ATHILA1"        "ATHILA2"       
# [5] "ATHILA3"        "ATHILA4"        "ATHILA4a"       "ATHILA4c"      
# [9] "ATHILA5"        "ATHILA6a"       "ATHILA6b"       "ATHILA7"       
#[13] "ATHILA7a"       "ATHILA8a"       "ATHILA9"

ATHILA_BED_colophylo <- data.frame(ATHILA_BED,
                                   thickStart = as.integer(0),
                                   thickEnd = as.integer(0),
                                   itemRgb = ".")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA_general",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[12])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA0",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[11])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA1",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[10])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA2",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[9])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA3",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[8])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA4",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[7])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA4a",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[7])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA4c",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[7])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA5",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[6])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA6a",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[5])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA6b",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[5])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA7",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[4])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA7a",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[4])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA8a",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[3])), collapse = ",")
ATHILA_BED_colophylo[ATHILA_BED_colophylo$score == "ATHILA9",]$itemRgb <- paste(as.vector(col2rgb(rich12()(12)[2])), collapse = ",")
ATHILA_BED_colophylo$score <- as.integer(0)
write.table(ATHILA_BED_colophylo,
            file = paste0(allDir_ATHILA, "ATHILA_in_66Atha_colophylo_",
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
foreach(i = 1:length(acc)) %dorng% {
#for(i in 1:length(acc)) {
  print(acc[i])
  acc_chrs <- read.table(paste0("/home/ajt200/analysis/pancentromere/assemblies/",
                                acc[i], ".fa.fai"),
                         header = F)[,1]
  acc_chrs <- gsub("_RagTag_RagTag", "", acc_chrs)
  acc_chrs <- gsub("chr", "Chr", acc_chrs)
  acc_chrs <- gsub("SUPER_", "Chr", acc_chrs)
  acc_chrLens <- read.table(paste0("/home/ajt200/analysis/pancentromere/assemblies/",
                                   acc[i], ".fa.fai"),
                            header = F)[,2]
  acc_chrLens <- acc_chrLens[which(acc_chrs %in% chrName)]
  acc_chrs <- acc_chrs[which(acc_chrs %in% chrName)]
  print(acc_chrs)
  print(acc_chrLens)
  acc_chrLens <- acc_chrLens[sort.int(acc_chrs, index.return = T)$ix]
  acc_chrs <- acc_chrs[sort.int(acc_chrs, index.return = T)$ix]
  print(acc_chrs)
  print(acc_chrLens)

  acc_CEN <- CEN[grep(acc[i], CEN$fasta.name),]
  acc_CEN_GR <- GRanges(seqnames = acc_CEN$chr,
                        ranges = IRanges(start = acc_CEN$start,
                                         end = acc_CEN$end),
                        strand = "*",
                        acc = acc_CEN$fasta.name,
                        region = acc_CEN$region)
  print(acc_CEN_GR)

  acc_ATHILA <- ATHILA[which(ATHILA$acc == acc[i]),]
  acc_ATHILA_GR <- ATHILA_GR[which(ATHILA_GR$acc == acc[i])]
  acc_ATHILA_BED <- data.frame(chr = as.character(seqnames(acc_ATHILA_GR)),
                               start = start(acc_ATHILA_GR)-1,
                               end = end(acc_ATHILA_GR),
                               name = acc_ATHILA_GR$TE_ID,
                               score = acc_ATHILA_GR$phylo,
                               strand = strand(acc_ATHILA_GR))

  acc_CENATHILA <- CENATHILA[which(CENATHILA$acc == acc[i]),]
  acc_CENATHILA_GR <- CENATHILA_GR[which(CENATHILA_GR$acc == acc[i])]
  acc_CENATHILA_BED <- data.frame(chr = as.character(seqnames(acc_CENATHILA_GR)),
                                  start = start(acc_CENATHILA_GR)-1,
                                  end = end(acc_CENATHILA_GR),
                                  name = acc_CENATHILA_GR$TE_ID,
                                  score = acc_CENATHILA_GR$phylo,
                                  strand = strand(acc_CENATHILA_GR))

  acc_nonCENATHILA <- nonCENATHILA[which(nonCENATHILA$acc == acc[i]),]
  acc_nonCENATHILA_GR <- nonCENATHILA_GR[which(nonCENATHILA_GR$acc == acc[i])]
  acc_nonCENATHILA_BED <- data.frame(chr = as.character(seqnames(acc_nonCENATHILA_GR)),
                                     start = start(acc_nonCENATHILA_GR)-1,
                                     end = end(acc_nonCENATHILA_GR),
                                     name = acc_nonCENATHILA_GR$TE_ID,
                                     score = acc_nonCENATHILA_GR$phylo,
                                     strand = strand(acc_nonCENATHILA_GR))

  stopifnot(length(acc_CENATHILA_GR) + length(acc_nonCENATHILA_GR) == length(acc_ATHILA_GR))

  write.table(acc_CENATHILA_BED,
              file = paste0(accDir_ATHILA[i], "CENATHILA_in_", acc[i], "_",
                            paste0(chrName, collapse = "_"), ".bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)

  write.table(acc_nonCENATHILA_BED,
              file = paste0(accDir_ATHILA[i], "nonCENATHILA_in_", acc[i], "_",
                            paste0(chrName, collapse = "_"), ".bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)


  acc_soloLTR <- soloLTR[which(soloLTR$acc == acc[i]),]
  acc_soloLTR_GR <- soloLTR_GR[which(soloLTR_GR$acc == acc[i])]
  acc_soloLTR_BED <- data.frame(chr = as.character(seqnames(acc_soloLTR_GR)),
                                start = start(acc_soloLTR_GR)-1,
                                end = end(acc_soloLTR_GR),
                                name = acc_soloLTR_GR$TE_ID,
                                score = acc_soloLTR_GR$phylo,
                                strand = strand(acc_soloLTR_GR))

  acc_CENsoloLTR <- CENsoloLTR[which(CENsoloLTR$acc == acc[i]),]
  acc_CENsoloLTR_GR <- CENsoloLTR_GR[which(CENsoloLTR_GR$acc == acc[i])]
  acc_CENsoloLTR_BED <- data.frame(chr = as.character(seqnames(acc_CENsoloLTR_GR)),
                                   start = start(acc_CENsoloLTR_GR)-1,
                                   end = end(acc_CENsoloLTR_GR),
                                   name = acc_CENsoloLTR_GR$TE_ID,
                                   score = acc_CENsoloLTR_GR$phylo,
                                   strand = strand(acc_CENsoloLTR_GR))

  acc_nonCENsoloLTR <- nonCENsoloLTR[which(nonCENsoloLTR$acc == acc[i]),]
  acc_nonCENsoloLTR_GR <- nonCENsoloLTR_GR[which(nonCENsoloLTR_GR$acc == acc[i])]
  acc_nonCENsoloLTR_BED <- data.frame(chr = as.character(seqnames(acc_nonCENsoloLTR_GR)),
                                      start = start(acc_nonCENsoloLTR_GR)-1,
                                      end = end(acc_nonCENsoloLTR_GR),
                                      name = acc_nonCENsoloLTR_GR$TE_ID,
                                      score = acc_nonCENsoloLTR_GR$phylo,
                                      strand = strand(acc_nonCENsoloLTR_GR))

  stopifnot(length(acc_CENsoloLTR_GR) + length(acc_nonCENsoloLTR_GR) == length(acc_soloLTR_GR))

  write.table(acc_CENsoloLTR_BED,
              file = paste0(accDir_soloLTR[i], "CENsoloLTR_in_", acc[i], "_",
                            paste0(chrName, collapse = "_"), ".bed"),
              quote = F, sep = "\t", row.names = F, col.names = F)

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

      ## Contract chr_acc_CEN_GR so that acc_CENranLoc_GR
      ## do not extend beyond centromeric coordinates
      #end(chr_acc_CEN_GR) <- end(chr_acc_CEN_GR)-max(width(chr_acc_CENATHILA_GR))
      #start(chr_acc_CEN_GR) <- start(chr_acc_CEN_GR)+max(width(chr_acc_CENATHILA_GR))

      # Define seed so that random selections are reproducible
      set.seed(76492749)
      chr_acc_CENranLoc_Start <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(chr_acc_CEN_GR), function(x) {
                                                                          ( start(chr_acc_CEN_GR[x]) + max(width(chr_acc_CENATHILA_GR)) ) :
                                                                          ( end(chr_acc_CEN_GR[x]) - max(width(chr_acc_CENATHILA_GR)) )
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
                                                                             c( ( max(width(chr_acc_nonCENATHILA_GR)) + 2000 ):
                                                                                ( acc_chrLens[j] - max(width(chr_acc_nonCENATHILA_GR)) - 2000 ) )[
                                                                               which(
                                                                                 !( ( max(width(chr_acc_nonCENATHILA_GR)) + 2000 ):
                                                                                    ( acc_chrLens[j] - max(width(chr_acc_nonCENATHILA_GR)) - 2000 ) %in%
                                                                                    ( start(chr_acc_CEN_GR[x]) - max(width(chr_acc_nonCENATHILA_GR)) ) :
                                                                                    ( end(chr_acc_CEN_GR[x]) + max(width(chr_acc_nonCENATHILA_GR)) )
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

