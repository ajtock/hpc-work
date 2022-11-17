#!/usr/bin/env Rscript

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 24/10/2022

# Extract putative recombination events detected in Col-0/Ler-0
# hybrid ONT reads identified on the basis that they contain >=
# a threshold number of accession-specific, chromosome-specific k-mers

# Usage:
# conda activate python_3.9.6
# ./get_chr_specific_segment_pairs_alnToSame.R ColLerF1pollen_1000bp_minq90 Col-0.ragtag_scaffolds Ler-0_110x.ragtag_scaffolds not_centromere Col-0.ragtag_scaffolds_Chr 24 0.9 11 30000 0.90 co 'Chr1,Chr2,Chr3,Chr4,Chr5'
# conda deactivate

#readsPrefix = "ColLerF1pollen_1000bp_minq90"
#acc1 = "Col-0.ragtag_scaffolds"
#acc2 = "Ler-0_110x.ragtag_scaffolds"
#region = "not_centromere"
#alnTo = "Col-0.ragtag_scaffolds_Chr"
#kmerSize = 24
#overlapProp = 0.9
#minHits = 11
#maxDist = 30000
#alenTOqlen = 0.90
#recombType = "co"
#chrName = unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split=","))

args = commandArgs(trailingOnly=T)
readsPrefix = args[1]
acc1 = args[2]
acc2 = args[3]
region = args[4]
alnTo = args[5]
kmerSize = as.integer(args[6])
overlapProp = as.numeric(args[7])
minHits = as.integer(args[8])
maxDist = as.integer(args[9])
alenTOqlen = as.numeric(args[10])
recombType = args[11]
chrName = unlist(strsplit(args[12],
                          split=","))


if(floor(log10(maxDist)) + 1 < 4) {
    maxDistName = paste0(maxDist, "bp")
} else if(floor(log10(maxDist)) + 1 >= 4 &
          floor(log10(maxDist)) + 1 <= 6) {
    maxDistName = paste0(maxDist/1e3, "kb")
} else if(floor(log10(maxDist)) + 1 >= 7) {
    maxDistName = paste0(maxDist/1e6, "Mb")
}

options(stringsAsFactors=F)
options(scipen=999)
library(parallel)
library(GenomicRanges)
library(circlize)
library(ComplexHeatmap)
library(gridBase)
library(viridis)
library(colorspace)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(seqinr)


outDir = paste0(region, "/segment_pairs/", recombType, "/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))

# Accession names
acc1_name = strsplit( strsplit(acc1, split="\\.")[[1]][1],
                      split="_" )[[1]][1]
acc2_name = strsplit( strsplit(acc2, split="\\.")[[1]][1],
                      split="_" )[[1]][1]

# Directories containing read segment alignment files
acc1_indir_list = lapply(1:length(chrName), function(x) {
    paste0(region, "/", chrName[x], "/segments/", acc1_name, "/", recombType, "/")
})
acc2_indir_list = lapply(1:length(chrName), function(x) {
    paste0(region, "/", chrName[x], "/segments/", acc2_name, "/", recombType, "/")
})

# CEN coordinates
CEN = read.csv(paste0("/rds/project/rds-O5Ty9yVfQKg/pancentromere/centromeric_coordinates/",
                      "centromere_manual_EDTA4_fa.csv"),
               header = T)
CEN$fasta.name = gsub(".fa", "", CEN$fasta.name)

# Genomic definitions
#acc1
acc1_fai = read.table(paste0("index/", acc1, ".fa.fai"), header=F)
acc1_chrs = acc1_fai[which(acc1_fai$V1 %in% chrName),]$V1
acc1_chrLens = acc1_fai[which(acc1_fai$V1 %in% chrName),]$V2

acc1_CEN = CEN[grep(acc1, CEN$fasta.name),]
acc1_CEN = acc1_CEN[,which(colnames(acc1_CEN) %in% c("chr", "start", "end"))]
acc1_CEN_new = data.frame()
for(i in 1:length(acc1_chrs)) {
    acc1_CEN_chr = acc1_CEN[which(acc1_CEN$chr == acc1_chrs[i]),]
    if(nrow(acc1_CEN_chr) > 1) {
        acc1_CEN_chr = data.frame(chr=acc1_CEN_chr$chr[1],
                                  start=acc1_CEN_chr$start[1],
                                  end=acc1_CEN_chr$end[nrow(acc1_CEN_chr)])
    }
    acc1_CEN_new = rbind(acc1_CEN_new, acc1_CEN_chr)
}
acc1_CEN = acc1_CEN_new
acc1_CENstart = acc1_CEN$start
acc1_CENend = acc1_CEN$end
acc1_chrs = paste0(acc1_name, "_", acc1_chrs)

#acc2
acc2_fai = read.table(paste0("index/", acc2, ".fa.fai"), header=F)
acc2_chrs = acc2_fai[which(acc2_fai$V1 %in% chrName),]$V1
acc2_chrLens = acc2_fai[which(acc2_fai$V1 %in% chrName),]$V2

acc2_CEN = CEN[grep(acc2, CEN$fasta.name),]
acc2_CEN = acc2_CEN[,which(colnames(acc2_CEN) %in% c("chr", "start", "end"))]
acc2_CEN_new = data.frame()
for(i in 1:length(acc2_chrs)) {
    acc2_CEN_chr = acc2_CEN[which(acc2_CEN$chr == acc2_chrs[i]),]
    if(nrow(acc2_CEN_chr) > 1) {
        acc2_CEN_chr = data.frame(chr=acc2_CEN_chr$chr[1],
                                  start=acc2_CEN_chr$start[1],
                                  end=acc2_CEN_chr$end[nrow(acc2_CEN_chr)])
    }
    acc2_CEN_new = rbind(acc2_CEN_new, acc2_CEN_chr)
}
acc2_CEN = acc2_CEN_new
acc2_CENstart = acc2_CEN$start
acc2_CENend = acc2_CEN$end
acc2_chrs = paste0(acc2_name, "_", acc2_chrs)


# Load read segment alignment files as a combined data.frame
load_pafs = function(indir, acc_name, aln_acc, suffix, aligner) {
    #indir=acc1_indir_list[[1]]
    #acc_name=acc1_name
    #suffix=paste0("_alnTo_", alnTo, "_wm_ont.paf")
    #aligner="wm"
    files = system(paste0("find ", indir, " -type f -name ", "'*", acc_name, suffix, "'", " -print"), intern=T)
    if(length(files) > 0) {
        aln_DF = data.frame()
        for(h in 1:length(files)) {
            aln = fread(files[h],
                        header=F, fill=T, sep="\t", data.table=F)[,1:13]
            aln_DF = rbind(aln_DF, aln)
        }
        aln_DF$aligner = aligner
        colnames(aln_DF) = c("qname", "qlen", "qstart0", "qend0",
                             "strand", "tname", "tlen", "tstart", "tend",
                             "nmatch", "alen", "mapq", "atype", "aligner")

    return(aln_DF)
    }
}


## wm alignments
#acc1_wm_list = mclapply(1:length(acc1_indir_list), function(x) {
#    load_pafs(indir=acc1_indir_list[[x]], acc_name=acc1_name, suffix=paste0("_alnTo_", alnTo, "_wm_ont.paf"), aligner="wm")
#}, mc.preschedule=F, mc.cores=length(acc1_indir_list))
#acc2_wm_list = mclapply(1:length(acc2_indir_list), function(x) {
#    load_pafs(indir=acc2_indir_list[[x]], acc_name=acc2_name, suffix=paste0("_alnTo_", alnTo, "_wm_ont.paf"), aligner="wm")
#}, mc.preschedule=F, mc.cores=length(acc2_indir_list))

# mm alignments
acc1_mm_list = mclapply(1:length(acc1_indir_list), function(x) {
    load_pafs(indir=acc1_indir_list[[x]], acc_name=acc1_name, suffix=paste0("_alnTo_", alnTo, "_mm_ont.paf"), aligner="mm")
}, mc.preschedule=F, mc.cores=length(acc1_indir_list))
acc2_mm_list = mclapply(1:length(acc2_indir_list), function(x) {
    load_pafs(indir=acc2_indir_list[[x]], acc_name=acc2_name, suffix=paste0("_alnTo_", alnTo, "_mm_ont.paf"), aligner="mm")
}, mc.preschedule=F, mc.cores=length(acc2_indir_list))

# sr alignments
acc1_sr_list = mclapply(1:length(acc1_indir_list), function(x) {
    load_pafs(indir=acc1_indir_list[[x]], acc_name=acc1_name, suffix=paste0("_alnTo_", alnTo, "_mm_sr.paf"), aligner="sr")
}, mc.preschedule=F, mc.cores=length(acc1_indir_list))
acc2_sr_list = mclapply(1:length(acc2_indir_list), function(x) {
    load_pafs(indir=acc2_indir_list[[x]], acc_name=acc2_name, suffix=paste0("_alnTo_", alnTo, "_mm_sr.paf"), aligner="sr")
}, mc.preschedule=F, mc.cores=length(acc2_indir_list))


# acc alignments - chromosome list of aligner lists
acc1_aln_chr_list_of_lists = lapply(1:length(acc1_mm_list), function(x) {
    list(acc1_mm_list[[x]], acc1_sr_list[[x]]) ## acc1_wm_list[[x]], removed as wm alignment not produced 
})
for(x in 1:length(acc1_aln_chr_list_of_lists)) {
    names(acc1_aln_chr_list_of_lists[[x]]) = c("mm", "sr") ## "wm"
}

acc2_aln_chr_list_of_lists = lapply(1:length(acc2_mm_list), function(x) {
    list(acc2_mm_list[[x]], acc2_sr_list[[x]]) ## acc2_wm_list[[x]], removed as wm alignment not produced
})
for(x in 1:length(acc2_aln_chr_list_of_lists)) {
    names(acc2_aln_chr_list_of_lists[[x]]) = c("mm", "sr") ## "wm" removed
}


# Get best pair of acc1 and acc2 read segment alignments, based on:
# 1. The alignment length (alen)
# 2. The alignment number of matching bases (nmatch)
# 3. The alignment strand
aln_best_pair = function(acc1_aln_DF_list, acc2_aln_DF_list) {
    #acc1_aln_DF_list=acc1_aln_list
    #acc2_aln_DF_list=acc2_aln_list
    # Each of the 3 list elements in acc1_aln_DF_list and acc2_aln_DF_list is a
    # a data.frame of alignments done by wm_ont, mm_ont or mm_sr
    acc1_aln_DF_bind_rows = dplyr::bind_rows(acc1_aln_DF_list)
    acc2_aln_DF_bind_rows = dplyr::bind_rows(acc2_aln_DF_list)

    # Get the best alignment from acc1_aln_DF and corresponding row from acc2_aln_DF
    acc1_aln_DF_best = data.frame()
    acc2_aln_DF_best = data.frame()
    for(read_id in unique(acc1_aln_DF_bind_rows$qname)) {
        acc1_aln_DF_read_id = acc1_aln_DF_bind_rows %>%
            dplyr::filter(qname == read_id)
        acc1_aln_DF_read_id = acc1_aln_DF_read_id[ with(acc1_aln_DF_read_id,
                                                        order(alen, nmatch, decreasing=T)), ]
        acc1_aln_DF_read_id_select = acc1_aln_DF_read_id[1, ]
        acc1_aln_DF_best = rbind(acc1_aln_DF_best, acc1_aln_DF_read_id_select)

        acc2_aln_DF_read_id = acc2_aln_DF_bind_rows %>%
            dplyr::filter(qname == read_id)
        acc2_aln_DF_read_id = acc2_aln_DF_read_id[ with(acc2_aln_DF_read_id,
                                                        order(alen, nmatch, decreasing=T)), ]
        acc2_aln_DF_read_id_strand = acc2_aln_DF_read_id[ which(acc2_aln_DF_read_id$strand == acc1_aln_DF_read_id_select$strand), ]
        if(nrow(acc2_aln_DF_read_id_strand) > 0) {
            if(acc2_aln_DF_read_id_strand[1, ]$alen == acc2_aln_DF_read_id[1, ]$alen &
               acc2_aln_DF_read_id_strand[1, ]$nmatch == acc2_aln_DF_read_id[1, ]$nmatch) {
                acc2_aln_DF_read_id_select = acc2_aln_DF_read_id_strand[1, ]
            } else {
                acc2_aln_DF_read_id_select = acc2_aln_DF_read_id[1, ]
            }
        } else {
            acc2_aln_DF_read_id_select = acc2_aln_DF_read_id[1, ]
        }
        acc2_aln_DF_best = rbind(acc2_aln_DF_best, acc2_aln_DF_read_id_select)
    }

    colnames(acc1_aln_DF_best) = paste0("acc1_", colnames(acc1_aln_DF_best))
    colnames(acc2_aln_DF_best) = paste0("acc2_", colnames(acc2_aln_DF_best))
    stopifnot(identical(acc1_aln_DF_best$acc1_qname, acc2_aln_DF_best$acc2_qname))

    #aln_best_pair_DF = cbind(acc1_aln_DF_best, acc2_aln_DF_best)
    #aln_best_pair_DF = aln_best_pair_DF[ , -which(colnames(aln_best_pair_DF) == "acc2_qname") ]

    aln_best_pair_DF = base::merge(x = acc1_aln_DF_best,
                                   y = acc2_aln_DF_best,
                                   by.x = "acc1_qname",
                                   by.y = "acc2_qname")

    colnames(aln_best_pair_DF)[which(colnames(aln_best_pair_DF) == "acc1_qname")] = "qname"
    aln_best_pair_DF = aln_best_pair_DF[ with(aln_best_pair_DF,
                                              order(acc1_tname, acc1_tstart, acc1_tend, decreasing=F)), ]

    return(aln_best_pair_DF)
}


# Get best pair of aligned read segments for each read
aln_best_pair_DF = dplyr::bind_rows(
    mclapply(1:length(acc1_aln_chr_list_of_lists), function(x) {
        aln_best_pair(acc1_aln_DF_list=acc1_aln_chr_list_of_lists[[x]], acc2_aln_DF_list=acc2_aln_chr_list_of_lists[[x]])
    }, mc.preschedule=F, mc.cores=length(acc1_aln_chr_list_of_lists))
)

#aln_best_pair_DF = aln_best_pair_DF %>%
#    dplyr::filter(acc1_alen >= 200) %>%
#    dplyr::filter(acc2_alen >= 200)

print(paste0(nrow(aln_best_pair_DF), " putative ", recombType, " events"))

# Filter to retain putative recombination events between homologous chromosomes only
aln_best_pair_hom_DF = aln_best_pair_DF[ which(aln_best_pair_DF$acc1_tname == aln_best_pair_DF$acc2_tname), ]

print(paste0(nrow(aln_best_pair_hom_DF), " putative ", recombType, " events are between homologous chromosomes"))

print(paste0( round( ( nrow(aln_best_pair_hom_DF) / nrow(aln_best_pair_DF) ), 4 ) * 100, "% of putative ", recombType, " events are between homologous chromosomes"))

# Filter to retain putative recombination events between homologous chromosomes where
# the per-accession read segments align to within maxDist bp of each other in the same reference assembly
aln_dist_acc1_tstart_acc2_tstart = abs(aln_best_pair_hom_DF$acc1_tstart - aln_best_pair_hom_DF$acc2_tstart) + 1
aln_dist_acc1_tstart_acc2_tend = abs(aln_best_pair_hom_DF$acc1_tstart - aln_best_pair_hom_DF$acc2_tend) + 1
aln_dist_acc1_tend_acc2_tstart = abs(aln_best_pair_hom_DF$acc1_tend - aln_best_pair_hom_DF$acc2_tstart) + 1
aln_dist_acc1_tend_acc2_tend = abs(aln_best_pair_hom_DF$acc1_tend - aln_best_pair_hom_DF$acc2_tend) + 1

aln_dist_min = pmin(aln_dist_acc1_tstart_acc2_tstart, aln_dist_acc1_tstart_acc2_tend,
                    aln_dist_acc1_tend_acc2_tstart, aln_dist_acc1_tend_acc2_tend, na.rm = T)
aln_dist_max = pmax(aln_dist_acc1_tstart_acc2_tstart, aln_dist_acc1_tstart_acc2_tend,
                    aln_dist_acc1_tend_acc2_tstart, aln_dist_acc1_tend_acc2_tend, na.rm = T)

event_start = pmin(aln_best_pair_hom_DF$acc1_tstart, aln_best_pair_hom_DF$acc1_tend,
                   aln_best_pair_hom_DF$acc2_tstart, aln_best_pair_hom_DF$acc2_tend, na.rm = T)
event_end = pmax(aln_best_pair_hom_DF$acc1_tstart, aln_best_pair_hom_DF$acc1_tend,
                 aln_best_pair_hom_DF$acc2_tstart, aln_best_pair_hom_DF$acc2_tend, na.rm = T)
event_midpoint = round((event_start + event_end) / 2)

aln_best_pair_hom_DF$aln_dist_min = aln_dist_min
aln_best_pair_hom_DF$aln_dist_max = aln_dist_max
aln_best_pair_hom_DF$event_start = event_start
aln_best_pair_hom_DF$event_end = event_end
aln_best_pair_hom_DF$event_midpoint = event_midpoint

# Get read lengths for hybrid read IDs in aln_best_pair_hom_DF$qname
# to be used for retaining alignment pairs where the Col and Ler
# read segments align to within a given distance of each other
# (e.g., the given hybrid read length) in the same assembly
hybrid_read_lengths = distinct(dplyr::bind_rows(
    mclapply(1:length(chrName), function(x) {
        tmp_list = read.fasta(paste0("fasta/", readsPrefix,
                                     "_match_", acc1, "_", region, "_", chrName[x],
                                     "_specific_k", kmerSize, "_downsampled_op", overlapProp, "_hits", minHits,
                                     "_match_", acc2, "_", region, "_", chrName[x],
                                     "_specific_k", kmerSize, "_downsampled_op", overlapProp, "_hits", minHits,
                                     ".fa"),
                              forceDNAtolower=F)
        read_length_DF = data.frame()
        for(i in 1:length(tmp_list)) {
            if( attr(tmp_list[[i]], which="name", exact=T) %in% aln_best_pair_hom_DF$qname ) {
                read_length_DF_i = data.frame(read_id = attr(tmp_list[[i]], which="name", exact=T),
                                              read_len = length(getSequence(tmp_list[[i]])))
                read_length_DF = rbind(read_length_DF, read_length_DF_i)
            } 
        }                                                              
        read_length_DF
    }, mc.preschedule=F, mc.cores=length(chrName))
))

aln_best_pair_hom_DF = base::merge(x = aln_best_pair_hom_DF,
                                   y = hybrid_read_lengths,
                                   by.x = "qname",
                                   by.y = "read_id",
                                   sort = F)

aln_best_pair_hom_maxDist_DF = data.frame()
for(x in 1:nrow(aln_best_pair_hom_DF)) {
    if(aln_best_pair_hom_DF[x, ]$aln_dist_min <= aln_best_pair_hom_DF[x, ]$read_len * 2) {
        aln_best_pair_hom_maxDist_DF = rbind(aln_best_pair_hom_maxDist_DF, aln_best_pair_hom_DF[x, ])
    }
}

print(paste0(nrow(aln_best_pair_hom_maxDist_DF), " putative ", recombType, " events are between homologous chromosomes where the per-accession read segments align to within hybrid-read-length bp of each other in the same reference assembly"))

print(paste0( round( ( nrow(aln_best_pair_hom_maxDist_DF) / nrow(aln_best_pair_DF) ), 4 ) * 100, "% of putative ", recombType, " events are between homologous chromosomes where the per-accession read segments align to within hybrid-read-length bp of each other in the same reference assembly"))

print(paste0( round( ( nrow(aln_best_pair_hom_maxDist_DF) / nrow(aln_best_pair_hom_DF) ), 4 ) * 100, "% of putative ", recombType, " events that are between homologous chromosomes where the per-accession read segments align to within hybrid-read-length bp of each other in the same reference assembly"))



# Filter to retain putative recombination events between homologous chromosomes where
# the per-accession read segments align to within maxDist bp of each other in the same reference assembly, AND
# where the per-accession alignment length is >= alenTOqlen of the segment length
aln_best_pair_hom_maxDist_alenTOqlen_DF = aln_best_pair_hom_maxDist_DF[ which(aln_best_pair_hom_maxDist_DF$acc1_alen /
                                                                              aln_best_pair_hom_maxDist_DF$acc1_qlen >= alenTOqlen &
                                                                              aln_best_pair_hom_maxDist_DF$acc2_alen /
                                                                              aln_best_pair_hom_maxDist_DF$acc2_qlen >= alenTOqlen), ]

print(paste0(nrow(aln_best_pair_hom_maxDist_alenTOqlen_DF), " putative ", recombType, " events are between homologous chromosomes where the per-accession read segments align to within hybrid-read-length bp of each other in the same reference assembly, and where the per-accession alignment length is >= ", alenTOqlen * 100, "% of the segment length"))

print(paste0( round( ( nrow(aln_best_pair_hom_maxDist_alenTOqlen_DF) / nrow(aln_best_pair_DF) ), 4 ) * 100, "% of putative ", recombType, " events are between homologous chromosomes where the per-accession read segments align to within hybrid-read-length bp of each other in the same reference assembly, and where the per-accession alignment length is >= ", alenTOqlen * 100, "% of the segment length"))

print(paste0( round( ( nrow(aln_best_pair_hom_maxDist_alenTOqlen_DF) / nrow(aln_best_pair_hom_maxDist_DF) ), 4 ) * 100, "% of putative ", recombType, " events that are between homologous chromosomes and where the per-accession read segments align to within hybrid-read-length bp of each other in the same reference assembly, are those where the per-accession alignment length is >= ", alenTOqlen * 100, "% of the segment length"))



#write.table(aln_best_pair_DF,
#            paste0(outDir, readsPrefix,
#                   "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
#                   "_", recombType,
#                   "_alnTo_", alnTo, "_",
#                   paste0(chrName, collapse="_"), ".tsv"),
#            quote=F, sep="\t", col.names=T, row.names=F)
#write.table(aln_best_pair_hom_DF,
#            paste0(outDir, readsPrefix,
#                   "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
#                   "_hom_", recombType,
#                   "_alnTo_", alnTo, "_",
#                   paste0(chrName, collapse="_"), ".tsv"),
#            quote=F, sep="\t", col.names=T, row.names=F)
write.table(aln_best_pair_hom_maxDist_DF,
            paste0(outDir, readsPrefix,
                   "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                   "_hom_maxDist_", recombType,
                   "_alnTo_", alnTo, "_",
                   paste0(chrName, collapse="_"), ".tsv"),
            quote=F, sep="\t", col.names=T, row.names=F)
write.table(aln_best_pair_hom_maxDist_alenTOqlen_DF,
            paste0(outDir, readsPrefix,
                   "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                   "_hom_maxDist_aTOq", alenTOqlen, "_", recombType,
                   "_alnTo_", alnTo, "_",
                   paste0(chrName, collapse="_"), ".tsv"),
            quote=F, sep="\t", col.names=T, row.names=F)


#aln_best_pair_DF = read.table(paste0(outDir, readsPrefix,
#                                     "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
#                                     "_", recombType,
#                                     "_alnTo_", alnTo, "_",
#                                     paste0(chrName, collapse="_"), ".tsv"),
#                              header=T)
#aln_best_pair_hom_DF = read.table(paste0(outDir, readsPrefix,
#                                         "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
#                                         "_hom_", recombType,
#                                         "_alnTo_", alnTo, "_",
#                                         paste0(chrName, collapse="_"), ".tsv"),
#                                  header=T)
aln_best_pair_hom_maxDist_DF = read.table(paste0(outDir, readsPrefix,
                                                 "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                                                 "_hom_maxDist_", recombType,
                                                 "_alnTo_", alnTo, "_",
                                                 paste0(chrName, collapse="_"), ".tsv"),
                                          header=T)
aln_best_pair_hom_maxDist_alenTOqlen_DF = read.table(paste0(outDir, readsPrefix,
                                                            "_", acc1, "_", acc2, "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                                                            "_hom_maxDist_aTOq", alenTOqlen, "_", recombType,
                                                            "_alnTo_", alnTo, "_",
                                                            paste0(chrName, collapse="_"), ".tsv"),
                                                     header=T)
