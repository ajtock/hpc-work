#!/usr/bin/env Rscript

# Get centromere and non-centromere coordinates and output in BED format

# Usage:
# conda activate python_3.9.6
# ./get_chr_region_coordinates.R Col-0.ragtag_scaffolds 'Chr1,Chr2,Chr3,Chr4,Chr5'
# ./get_chr_region_coordinates.R Ler-0_110x.ragtag_scaffolds 'Chr1,Chr2,Chr3,Chr4,Chr5'
# conda deactivate

#accName = "Col-0.ragtag_scaffolds"
#chrName = unlist(strsplit("Chr1",
#                          split = ","))

args = commandArgs(trailingOnly = T)
accName = args[1]
chrName = unlist(strsplit(args[2],
                          split = ","))

options(stringsAsFactors = F)

outDir = "bed/"
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system("[ -d logs/segments ] || mkdir -p logs/segments/")

# Define chromosomes and chromosome lengths
acc_chrs = read.table(paste0("/rds/project/rds-O5Ty9yVfQKg/pancentromere/assemblies/",
                             accName, ".fa.fai"),
                      header = F)[,1]
acc_chrs = gsub("_RagTag_RagTag", "", acc_chrs)
acc_chrs = gsub("chr", "Chr", acc_chrs)
acc_chrs = gsub("SUPER_", "Chr", acc_chrs)
acc_chrLens = read.table(paste0("/rds/project/rds-O5Ty9yVfQKg/pancentromere/assemblies/",
                                accName, ".fa.fai"),
                         header = F)[,2]
acc_ptgs = acc_chrs[-grep("Chr", acc_chrs)]
acc_ptgLens = acc_chrLens[-grep("Chr", acc_chrs)]
acc_chrLens = acc_chrLens[which(acc_chrs %in% chrName)]
acc_chrs = acc_chrs[which(acc_chrs %in% chrName)]
print(acc_ptgs)
print(acc_ptgLens)
print(acc_chrs)
print(acc_chrLens)
acc_ptgLens = acc_ptgLens[sort.int(acc_ptgs, index.return = T)$ix]
acc_ptgs = acc_ptgs[sort.int(acc_ptgs, index.return = T)$ix]
acc_chrLens = acc_chrLens[sort.int(acc_chrs, index.return = T)$ix]
acc_chrs = acc_chrs[sort.int(acc_chrs, index.return = T)$ix]
print(acc_ptgs)
print(acc_ptgLens)
print(acc_chrs)
print(acc_chrLens)

# Load CEN coordinates
CEN = read.csv(paste0("/rds/project/rds-O5Ty9yVfQKg/pancentromere/centromeric_coordinates/",
                      "centromere_manual_EDTA4_fa.csv"),
               header = T)
CEN$fasta.name = gsub(".fa", "", CEN$fasta.name)
acc_CEN = CEN[grep(accName, CEN$fasta.name),]
acc_CEN = acc_CEN[,which(colnames(acc_CEN) %in% c("chr", "start", "end"))]
acc_CEN_new = data.frame()
for(i in 1:length(acc_chrs)) {
  acc_CEN_chr = acc_CEN[which(acc_CEN$chr == acc_chrs[i]),]
  if(nrow(acc_CEN_chr) > 1) {
    acc_CEN_chr = data.frame(chr = acc_CEN_chr$chr[1],
                             start = acc_CEN_chr$start[1],
                             end = acc_CEN_chr$end[nrow(acc_CEN_chr)])
  }
  acc_CEN_new = rbind(acc_CEN_new, acc_CEN_chr)
}
acc_CEN = acc_CEN_new 

acc_nonCEN = data.frame(chr = rep(acc_CEN$chr, 2),
                        start = c(rep(1, nrow(acc_CEN)), acc_CEN$end+1),
                        end = c(acc_CEN$start-1, acc_chrLens))

acc_other = data.frame(chr = acc_ptgs,
                       start = rep(1, length(acc_ptgs)),
                       end = acc_ptgLens)   

 
# Order first by chromosome, then start, then end
acc_CEN = acc_CEN[ with(acc_CEN, order(chr, start)), ]
acc_nonCEN = acc_nonCEN[ with(acc_nonCEN, order(chr, start)), ]
acc_other = acc_other[ with(acc_other, order(chr, start)), ]

# Substract 1 from start coordinates for output in BED format
acc_CEN$start = acc_CEN$start-1
acc_nonCEN$start = acc_nonCEN$start-1
acc_other$start = acc_other$start-1

# Output in BED format
write.table(acc_CEN,
            file = paste0(outDir, accName,
                          "_centromere_",
                          paste0(chrName, collapse="_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
write.table(acc_nonCEN,
            file = paste0(outDir, accName,
                          "_not_centromere_",
                          paste0(chrName, collapse="_"), ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
#write.table(acc_other,
#            file = paste0(outDir, accName,
#                          "_other.bed"),
#            quote = F, sep = "\t", row.names = F, col.names = F)
