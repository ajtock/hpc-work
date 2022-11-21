#!/usr/bin/env Rscript

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 17/11/2022

# Do permutation tests to evaluate overlap between candidate recombinant
# read segment pair alignment intervals and features of interest

# Usage:
# conda activate python_3.9.6
# ./recomb_overlap_permTest.R ColLerF1pollen_1000bp_minq90 Col-0.ragtag_scaffolds Ler-0_110x.ragtag_scaffolds 24 0.9 11 Col-0.ragtag_scaffolds_Chr 0.90 not_centromere 'Chr1,Chr2,Chr3,Chr4,Chr5' 10000 co_and_nco
# conda deactivate

args = commandArgs(trailingOnly=T)
readsPrefix = args[1]
acc1 = args[2]
acc2 = args[3]
kmerSize = as.integer(args[4])
overlapProp = as.numeric(args[5])
minHits = as.integer(args[6])
alnTo = args[7]
alenTOqlen = args[8]
region = args[9]
chrom = unlist(strsplit(args[10], split=","))
perms = as.integer(args[11])
recombType = args[12] 

#readsPrefix = "ColLerF1pollen_1000bp_minq90"
#acc1 = "Col-0.ragtag_scaffolds"
#acc2 = "Ler-0_110x.ragtag_scaffolds"
#kmerSize = as.integer(24)
#overlapProp = 0.9
#minHits = as.integer(11)
#alnTo = "Col-0.ragtag_scaffolds_Chr"
#alenTOqlen = 0.90
#region = "not_centromere"
#chrom = unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                        split=","))
#perms = as.integer(10000)
#recombType = "co_and_nco"

# Set minimum possible P-value for permutation test result with
# perms sets of random loci
min_pval = 1 - ( (perms - 1) / perms)


library(regioneR)
library(rtracklayer)
library(dplyr)


indir = paste0(region, "/segment_pairs/")
outdir = paste0(indir, "perm_tests/")
plotdir = paste0(outdir, "plots/")
system(paste0("[ -d ", plotdir, " ] || mkdir -p ", plotdir))


# Accession names
acc1_name = strsplit( strsplit(acc1, split="\\.")[[1]][1],
                      split="_" )[[1]][1]
acc2_name = strsplit( strsplit(acc2, split="\\.")[[1]][1],
                      split="_" )[[1]][1]
alnTo_name =  strsplit( strsplit(alnTo, split="\\.")[[1]][1],
                        split="_" )[[1]][1] 
acc1_name = acc1.split(".")[0].split("_")[0]
acc2_name = acc2.split(".")[0].split("_")[0]
alnTo_name = alnTo.split(".")[0].split("_")[0]


# Filtered read segment alignment pairs TSV (candidate recombinant reads)
# COs
COs_DF_list = lapply(1:length(chrom), function(x) {
    COs_tsv = paste0(indir, "co/", readsPrefix,
                     "_", acc1, "_", acc2,
                     "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                     "_hom_maxDist_aTOq", alenTOqlen, "_co",
                     "_alnTo_", alnTo, "_",
                     chrom[x], ".tsv")
    # File exists sanity check
    COs_DF_x = read.table(COs_tsv, header=T)
    COs_DF_x_filt = COs_DF_x[which(COs_DF_x$acc1_tname == chrom[x]),]
    COs_DF_x_filt
})
if(length(COs_DF_list) > 1) {
    COs_DF = dplyr::bind_rows(COs_DF_list)
} else {
    COs_DF = COs_DF_list[[1]]
}

COs_GR = GRanges(seqnames=COs_DF$acc1_tname,
                 ranges=IRanges(start=COs_DF$event_start,
                                end=COs_DF$event_end),
                 strand="*")

# NCOs
NCOs_DF_list = lapply(1:length(chrom), function(x) {
    NCOs_tsv = paste0(indir, "nco/", readsPrefix,
                      "_", acc1, "_", acc2,
                      "_k", kmerSize, "_op", overlapProp, "_h", minHits,
                      "_hom_maxDist_aTOq", alenTOqlen, "_nco",
                      "_alnTo_", alnTo, "_",
                      chrom[x], ".tsv")
    # File exists sanity check
    NCOs_DF_x = read.table(NCOs_tsv, header=T)
    NCOs_DF_x_filt = NCOs_DF_x[which(NCOs_DF_x$acc1_tname == chrom[x]),]
    NCOs_DF_x_filt
})
if( length(NCOs_DF_list) > 1 ) {
    NCOs_DF = dplyr::bind_rows(NCOs_DF_list)
} else {
    NCOs_DF = NCOs_DF_list[[1]]
}

NCOs_GR = GRanges(seqnames=NCOs_DF$acc1_tname,
                  ranges=IRanges(start=NCOs_DF$event_start,
                                 end=NCOs_DF$event_end),
                  strand="*")


if(recombType == "co_and_nco") {
   recomb_GR = c(COs_GR, NCOs_GR) 
} else if(recombType == "co") {
   recomb_GR = COs_GR
} else if(recombtype == "nco") {
   recomb_GR = NCOs_GR
} else {
   stop("recombtype is not 'co_and_nco', 'co' or 'nco'")
}


# Genomic definitions

# CEN coordinates
CEN = read.csv(paste0("/rds/project/rds-O5Ty9yVfQKg/pancentromere/centromeric_coordinates/",
                      "centromere_manual_EDTA4_fa.csv"),
               header=T)
CEN$fasta.name = gsub(".fa", "", CEN$fasta.name)

# alnTo accession
alnTo_fai = read.table(paste0("index/", alnTo, ".fa.fai"),
                       header=F)
alnTo_fai =  alnTo_fai[which(alnTo_fai$V1 %in% chrom),]
alnTo_chrs = alnTo_fai$V1
alnTo_chrLens = alnTo_fai$V2


alnTo_CEN = CEN[which(CEN$fasta.name == gsub("_Chr", "", alnTo)),]
alnTo_CEN = alnTo_CEN[,which(colnames(alnTo_CEN) %in% c("chr", "start", "end"))]
alnTo_CEN_new = data.frame()
for(i in 1:length(alnTo_chrs)) {
  alnTo_CEN_chr = alnTo_CEN[which(alnTo_CEN$chr == alnTo_chrs[i]),]
  if(nrow(alnTo_CEN_chr) > 1) {
    alnTo_CEN_chr = data.frame(chr=alnTo_CEN_chr$chr[1],
                               start=alnTo_CEN_chr$start[1],
                               end=alnTo_CEN_chr$end[nrow(alnTo_CEN_chr)])
  }
  alnTo_CEN_new = rbind(alnTo_CEN_new, alnTo_CEN_chr)
}
alnTo_CEN = alnTo_CEN_new
alnTo_CEN = alnTo_CEN_new

alnTo_CEN_GR = GRanges(seqnames=alnTo_CEN$chr,
                       ranges=IRanges(start=alnTo_CEN$start,
                                      end=alnTo_CEN$end),
                       strand="*")

alnTo_nonCEN = data.frame(chr=rep(alnTo_chrs, 2),
                          start=c( rep(1, length(alnTo_chrs)), (alnTo_CEN$end + 1) ),
                          end=c( (alnTo_CEN$start - 1), alnTo_chrLens))

alnTo_nonCEN_GR = GRanges(seqnames=alnTo_nonCEN$chr,
                          ranges=IRanges(start=alnTo_nonCEN$start,
                                         end=alnTo_nonCEN$end),
                          strand="*")

# Define region to be analysed
if(region == "not_centromere") {
    region_GR = alnTo_nonCEN_GR
    mask_GR = alnTo_CEN_GR
} else if(region == "centromere") {
    region_GR = alnTo_CEN_GR
    mask_GR = alnTo_nonCEN_GR
} else {
    stop("region is not 'not_centromere' or 'centromere'")
}

genome_GR = GRanges(seqnames=alnTo_chrs,
                    ranges=IRanges(start=rep(1, length(alnTo_chrs)),
                                   end=alnTo_chrLens),
                    strand="*")


## Subset to include only those not overlapping masked region (e.g., centromere)                                                          
#mask_recomb_overlap = findOverlaps(query=mask_GR,
#                                subject=recomb_GR,
#                                type = "any",
#                                select = "all",
#                                ignore.strand = TRUE)
#recomb_GR = recomb_GR[-subjectHits(mask_recomb_overlap)]
print("Candidate recombination events:")
print(recomb_GR)


# Load genes
genes = readGFF(paste0("/rds/project/rds-O5Ty9yVfQKg/pancentromere/annotation/genes/",
                        alnTo, "/", alnTo, ".genes.gff3"))
genes = genes[which(genes$seqid %in% alnTo_chrs),]
genes = genes[which(genes$type == "mRNA"),]
print(dim(genes))

genes_GR = GRanges(seqnames=genes$seqid,
                   ranges=IRanges(start=genes$start,
                                  end=genes$end),
                   strand=genes$strand)

# Subset to include only those not overlapping masked region
mask_genes_overlap <- findOverlaps(query = mask_GR,
                                   subject = genes_GR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
genes_GR <- genes_GR[-subjectHits(mask_genes_overlap)]

genes_GR = unique(genes_GR)

# NOTE: Retain strand information until after obtaining promoters, etc.
# Overlap analysis should be strand-unaware

# Obtain 1000-bp gene promoters
promoters_GR = promoters(genes_GR, upstream=1000, downstream=0)
promoters_GR = unique(promoters_GR)
strand(promoters_GR) = "*"
print(promoters_GR)

# Obtain regions immediately downstream of gene TSSs (gene 5' ends: TSS to TSS+499 bp)
g5ends_GR = promoters(genes_GR, upstream=0, downstream=500)
g5ends_GR = unique(g5ends_GR)
strand(g5ends_GR) = "*"
print(g5ends_GR)

# Obtain regions relative to TTS using TTSplus()
TTSplus = function(x, upstream=100, downstream=1000, ...) {
    if(!isSingleNumber(upstream))
        stop("'upstream' must be a single integer")
    if(!is.integer(upstream))
        upstream = as.numeric(upstream)
    if(!isSingleNumber(downstream))
        stop("'downstream' must be a single integer")
    if(!is.integer(downstream))
        downstream = as.numeric(downstream)
    if(downstream < 0)
        stop("'downstream' must be an integer >= 0")
    if(any(strand(x) == "*"))
        warning("'*' ranges were treated as '+'")
    on_plus = which(strand(x) == "+" | strand(x) == "*")
    on_plus_TTS = end(x)[on_plus]
    start(x)[on_plus] = on_plus_TTS - upstream
    end(x)[on_plus] = on_plus_TTS + downstream
    on_minus = which(strand(x) == "-")
    on_minus_TTS = start(x)[on_minus]
    end(x)[on_minus] = on_minus_TTS + upstream
    start(x)[on_minus] = on_minus_TTS - downstream
    return(x)
}

# Obtain regions immediately upstream of gene TTSs (gene 3' ends: TTS to TTS-499 bp)
g3ends_GR = TTSplus(genes_GR, upstream=499, downstream=0)
g3ends_GR = unique(g3ends_GR)
strand(g3ends_GR) = "*"
print(g3ends_GR)

# Obtain 1000-bp gene terminators
terminators_GR = TTSplus(genes_GR, upstream=-1, downstream=1000)
terminators_GR = unique(terminators_GR)
strand(terminators_GR) = "*"
print(terminators_GR)

# Remove strand information from genes_GR
strand(genes_GR) = "*"
print(genes_GR)

# Feature names
features_names = c(
                   "genes",
                   "1kb_upstream_TSS",                  
                   "500bp_5prime_ends",
                   "500bp_3prime_ends",
                   "1kb_downstream_TTS"
                  )
# GRanges list
features_GR_list = c(
                     "genes"=genes_GR,
                     "1kb_upstream_TSS"=promoters_GR,
                     "500bp_5prime_ends"=g5ends_GR,
                     "500bp_3prime_ends"=g3ends_GR,
                     "1kb_downstream_TTS"=terminators_GR
                    )


# Perform permutation tests with randomized regions generated on a per chromosome basis
set.seed(47393573)
pt_recomb_vs_features = lapply(1:length(features_GR_list), function(x) {
    permTest(A=recomb_GR,
             B=features_GR_list[[x]],
             alternative="auto",
             ntimes=perms,
             randomize.function=randomizeRegions,
             genome=genome_GR,
             mask=mask_GR,
             allow.overlaps=TRUE,
             per.chromosome=TRUE,
             evaluate.function=numOverlaps,
             count.once=TRUE,
             mc.set.seed=FALSE,
             mc.cores=detectCores())
})
                                   
for(i in 1:length(ptPeaksOtherPerChrom)) {
  assign(paste0(otherNames[i]), ptPeaksOtherPerChrom[[i]])
}
save(ptPeaksOtherPerChrom,
     file = paste0("permTest_", as.character(perms), "perms_",
                   libNameChIP, "_peaks_vs_others_in_",
                   genomeName, "genome_", region,
                   ".RData"))

pt_dist_DF = data.frame(





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


# Define perms sets of random loci in region
# of the same number and width distribution as the candidate recombination events 
def ranLoc_start_select(coordinates, n):
    rng = default_rng(238435)    
    rng.choice(a=coordinates,
               size=n,
               replace=False)


def define_ranLoc(region_PR, features_PR):
#regions_PR = alnTo_nonCEN_PR 
#features_PR = COs_PR

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






# Concatenate all alignment files for the given chromosome,
# accession and aligner
# NOTE: for some reason the "find ..." approach below doesn't work
# from within python, so need to create and run an equivalent bash script:
#cat_cmd = ["find"] + \
#          ["/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + indir + "/"] + \
#          ["-mindepth", "1"] + \
#          ["-maxdepth", "1"] + \
#          ["-type", "f"] + \
#          ["-name", "*" + acc_name + suffix] + \
#          ["-exec cat {} + >> cat.paf"]
#subprocess.run(cat_cmd)
# For details on -exec, see https://stackoverflow.com/questions/2961673/find-missing-argument-to-exec
def cat_pafs(indir, acc_name, suffix):
    ##indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc1_indir_list[1]
    #indir=acc1_indir_list[1]
    #acc_name=acc1_name
    #suffix="_alnTo_" + alnTo + "_mm_ont.paf"
    """
    Concatenate all alignment files for the given chromosome,
    accession and aligner.
    """
    outdir_paf = indir + "/cat_paf"
    if not os.path.exists(outdir_paf):
        os.makedirs(outdir_paf)
    #
    out_paf = outdir_paf + "/all_segments__" + acc_name + suffix
    if os.path.exists(out_paf):
        print("Concatenated alignment file " + out_paf + " already exists!")
        return
    else:
        cat_pafs_script = indir + "/find_cat_pafs.sh"
        with open(cat_pafs_script, "w") as cat_pafs_script_handle:
            cat_pafs_script_handle.write("#!/bin/bash\n\n" + \
                                         "find " + indir + "/ \\\n" + \
                                         "  -mindepth 1 \\\n" +\
                                         "  -maxdepth 1 \\\n" +\
                                         "  -type f \\\n" + \
                                         "  -name '*" + acc_name + suffix + "' \\\n" + \
                                         "  -exec cat {} + >> " + out_paf)
        #
        subprocess.run(["bash", cat_pafs_script])
        return


for x in range(0, len(acc1_indir_list)):
    print(acc1_indir_list[x])
    cat_pafs(indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc1_indir_list[x],
             acc_name=acc1_name,
             suffix="_alnTo_" + alnTo + "_mm_ont.paf")
    cat_pafs(indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc1_indir_list[x],
             acc_name=acc1_name,
             suffix="_alnTo_" + alnTo + "_mm_sr.paf")

for x in range(0, len(acc2_indir_list)):
    print(acc2_indir_list[x])
    cat_pafs(indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc2_indir_list[x],
             acc_name=acc2_name,
             suffix="_alnTo_" + alnTo + "_mm_ont.paf")
    cat_pafs(indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc2_indir_list[x],
             acc_name=acc2_name,
             suffix="_alnTo_" + alnTo + "_mm_sr.paf")


# Load concatenated read segment alignment file as a DataFrame
def load_cat_paf(indir, acc_name, suffix, aligner):
    ##indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc1_indir_list[1]
    #indir=acc1_indir_list[1]
    #acc_name=acc1_name
    #suffix="_alnTo_" + alnTo + "_mm_ont.paf"
    #aligner="mm"
    """
    Load concatenated read segment alignment file as a DataFrame.
    """
    #cat_pafs(indir=indir, acc_name=acc_name, suffix=suffix)
    indir_paf = indir + "/cat_paf"
    in_paf = indir_paf + "/all_segments__" + acc_name + suffix
    aln_DF = pd.read_csv(in_paf,
                         sep="\t", header=None, usecols=list(range(0, 13)))
    aln_DF["aligner"] = aligner
    aln_DF.columns = ["qname", "qlen", "qstart0", "qend0",
                      "strand", "tname", "tlen", "tstart", "tend",
                      "nmatch", "alen", "mapq", "atype", "aligner"]
    #
    return aln_DF


# Load read segment alignment files as a combined DataFrame
def load_pafs_slowly(indir, acc_name, suffix, aligner):
    #indir=acc1_indir_list[0]
    #acc_name=acc1_name
    #suffix="_alnTo_" + alnTo + "_mm_ont.paf"
    #aligner="mm"
    find_cmd = ["find"] + \
               ["/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + indir + "/"] + \
               ["-type", "f"] + \
               ["-name", "*" + acc_name + suffix] + \
               ["-print"]
    find_output = subprocess.run(find_cmd, capture_output=True)
    files_bytes = find_output.stdout.splitlines()
    files = [x.decode("utf-8") for x in files_bytes]
    #
    if len(files) > 0:
        aln_DF = pd.DataFrame()
        for h in range(0, len(files)):
            aln = pd.read_csv(files[h],
                              sep="\t", header=None, usecols=list(range(0, 13)))
            aln_DF = pd.concat(objs=[aln_DF, aln],
                               axis=0,
                               ignore_index=True)
        #aln_DF = aln_DF.iloc[:, :13]
        aln_DF["aligner"] = aligner
        aln_DF.columns = ["qname", "qlen", "qstart0", "qend0",
                          "strand", "tname", "tlen", "tstart", "tend",
                          "nmatch", "alen", "mapq", "atype", "aligner"]
        #
        return aln_DF


# mm alignments
acc1_mm_list = []
for x in range(0, len(acc1_indir_list)):
    acc1_mm_Chr = load_cat_paf(indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc1_indir_list[x],
                               acc_name=acc1_name,
                               suffix="_alnTo_" + alnTo + "_mm_ont.paf",
                               aligner="mm")
    acc1_mm_list.append(acc1_mm_Chr)
    del acc1_mm_Chr
    gc.collect()

acc2_mm_list = []
for x in range(0, len(acc2_indir_list)):
    acc2_mm_Chr = load_cat_paf(indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc2_indir_list[x],
                               acc_name=acc2_name,
                               suffix="_alnTo_" + alnTo + "_mm_ont.paf",
                               aligner="mm")
    acc2_mm_list.append(acc2_mm_Chr)
    del acc2_mm_Chr
    gc.collect()

# sr alignments
acc1_sr_list = []
for x in range(0, len(acc1_indir_list)):
    acc1_sr_Chr = load_cat_paf(indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc1_indir_list[x],
                               acc_name=acc1_name,
                               suffix="_alnTo_" + alnTo + "_mm_sr.paf",
                               aligner="sr")
    acc1_sr_list.append(acc1_sr_Chr)
    del acc1_sr_Chr
    gc.collect()

acc2_sr_list = []
for x in range(0, len(acc2_indir_list)):
    acc2_sr_Chr = load_cat_paf(indir="/rds/project/rds-O5Ty9yVfQKg/Col_Ler_F1_pollen_data/nanopore_recomb/chr_specific/" + acc2_indir_list[x],
                               acc_name=acc2_name,
                               suffix="_alnTo_" + alnTo + "_mm_sr.paf",
                               aligner="sr")
    acc2_sr_list.append(acc2_sr_Chr)
    del acc2_sr_Chr
    gc.collect()


# acc alignments - chromosome list of aligner lists
acc1_aln_chr_nested_list = []
for x in range(0, len(acc1_mm_list)):
    acc1_aln_chr_nested_list.append([acc1_mm_list[x], acc1_sr_list[x]])

acc2_aln_chr_nested_list = []
for x in range(0, len(acc2_mm_list)):
    acc2_aln_chr_nested_list.append([acc2_mm_list[x], acc2_sr_list[x]])


# Get best pair of acc1 and acc2 read segment alignments, based on:
# 1. The alignment strand
# 2. The alignment length (alen)
# 3. The alignment number of matching bases (nmatch)
def aln_best_pair(acc1_aln_DF_list, acc2_aln_DF_list):
    #acc1_aln_DF_list=acc1_aln_chr_nested_list[0]
    #acc2_aln_DF_list=acc2_aln_chr_nested_list[0]
    """
    Get the best pair of acc1 and acc2 read segment alignments, based on:
    1. The alignment length (alen)
    2. The alignment number of matching bases (nmatch)
    3. The alignment strand
    """
    # Each of the 2 list elements in acc1_aln_DF_list and acc2_aln_DF_list is
    # a DataFrame of alignments done by mm_ont or mm_sr
    acc1_aln_DF_concat = pd.concat(objs=acc1_aln_DF_list, axis=0, ignore_index=True)
    acc2_aln_DF_concat = pd.concat(objs=acc2_aln_DF_list, axis=0, ignore_index=True)
    # 
    # For each read ID, get the best alignment from each of acc1_aln_DF_concat and
    # and acc2_aln_DF_concat
    # acc1
    acc1_aln_DF_concat_sort = acc1_aln_DF_concat.sort_values(by=["qname", "alen", "nmatch"],
                                                             axis=0,
                                                             ascending=[True, False, False],
                                                             kind="quicksort",
                                                             ignore_index=True)
    acc1_aln_DF_concat_sort_list = list(acc1_aln_DF_concat_sort.groupby("qname"))
    acc1_aln_DF_best = pd.DataFrame()
    for read_tuple in acc1_aln_DF_concat_sort_list:
        read_tuple_aln_DF_best = read_tuple[1].iloc[[0]]
        acc1_aln_DF_best = pd.concat(objs=[acc1_aln_DF_best, read_tuple_aln_DF_best],
                                     axis=0,
                                     ignore_index=True)
    del acc1_aln_DF_concat, acc1_aln_DF_concat_sort, acc1_aln_DF_concat_sort_list, read_tuple_aln_DF_best
    gc.collect()
    # acc2
    acc2_aln_DF_best = pd.DataFrame()
    for read_id in list(acc1_aln_DF_best["qname"]):
        #print(read_id)
        acc1_aln_DF_read_id = acc1_aln_DF_best[acc1_aln_DF_best["qname"] == read_id] 
        acc2_aln_DF_read_id = acc2_aln_DF_concat[acc2_aln_DF_concat["qname"] == read_id] 
        acc2_aln_DF_read_id_sort = acc2_aln_DF_read_id.sort_values(by=["alen", "nmatch"],
                                                                   axis=0,
                                                                   ascending=[False, False],
                                                                   kind="quicksort",
                                                                   ignore_index=True)
        acc2_aln_DF_read_id_sort_strand = acc2_aln_DF_read_id_sort[acc2_aln_DF_read_id_sort["strand"] == acc1_aln_DF_read_id["strand"].iloc[0]]
        if acc2_aln_DF_read_id_sort_strand.shape[0] > 0:
            acc2_aln_DF_read_id_sort_select = acc2_aln_DF_read_id_sort_strand.iloc[[0]]
            ## If "alen" and "nmatch" are to be given priority over finding pairs where both alignments are to the same strand:
            #if acc2_aln_DF_read_id_sort_strand.iloc[0]["alen"] == acc2_aln_DF_read_id_sort.iloc[0]["alen"] and \
            #   acc2_aln_DF_read_id_sort_strand.iloc[0]["nmatch"] == acc2_aln_DF_read_id_sort.iloc[0]["nmatch"]:
            #    acc2_aln_DF_read_id_sort_select = acc2_aln_DF_read_id_sort_strand.iloc[[0]]
            #else:
            #    acc2_aln_DF_read_id_sort_select = acc2_aln_DF_read_id_sort.iloc[[0]]
        else:
            acc2_aln_DF_read_id_sort_select = acc2_aln_DF_read_id_sort.iloc[[0]]
        acc2_aln_DF_best = pd.concat(objs=[acc2_aln_DF_best, acc2_aln_DF_read_id_sort_select],
                                     axis=0,
                                     ignore_index=True)
    del acc2_aln_DF_concat, acc1_aln_DF_read_id, acc2_aln_DF_read_id, acc2_aln_DF_read_id_sort, acc2_aln_DF_read_id_sort_strand, acc2_aln_DF_read_id_sort_select
    gc.collect()
    #
    acc1_aln_DF_best.columns = "acc1_" + acc1_aln_DF_best.columns 
    acc2_aln_DF_best.columns = "acc2_" + acc2_aln_DF_best.columns 
    #
    # Stop if read ID order differs between acc1_aln_DF_best and acc2_aln_DF_best
    #list(acc1_aln_DF_best["acc1_qname"]) == list(acc2_aln_DF_best["acc2_qname"])
    if not acc1_aln_DF_best["acc1_qname"].equals(acc2_aln_DF_best["acc2_qname"]):
        print("Stopping because read ID order differs between acc1_aln_DF_best and acc2_aln_DF_best")
        return
    #
    aln_best_pair_DF = pd.merge(left=acc1_aln_DF_best, right=acc2_aln_DF_best,
                                how="inner", left_on="acc1_qname", right_on="acc2_qname")
    aln_best_pair_DF = aln_best_pair_DF.rename(columns = {"acc1_qname":"qname"})
    aln_best_pair_DF = aln_best_pair_DF.drop(columns="acc2_qname")
    #
    aln_best_pair_DF_sort = aln_best_pair_DF.sort_values(by=["acc1_tname", "acc1_tstart", "acc1_tend"],
                                                         axis=0,
                                                         ascending=[True, True, True],
                                                         kind="quicksort",
                                                         ignore_index=True)
    del acc1_aln_DF_best, acc2_aln_DF_best, aln_best_pair_DF
    gc.collect()
    #
    return aln_best_pair_DF_sort


# Get best pair of aligned read segments for each read
aln_best_pair_DF_list = []
for x in range(0, len(acc1_aln_chr_nested_list)):
    aln_best_pair_DF_x = aln_best_pair(acc1_aln_DF_list=acc1_aln_chr_nested_list[x],
                                       acc2_aln_DF_list=acc2_aln_chr_nested_list[x])
    aln_best_pair_DF_list.append(aln_best_pair_DF_x)

# Concatenate list elements (per-chromosome pd.DataFrames) into one pd.DataFrame
if len(aln_best_pair_DF_list) > 1:
    aln_best_pair_DF = pd.concat(objs=aln_best_pair_DF_list, axis=0, ignore_index=True)
else:
    aln_best_pair_DF = aln_best_pair_DF_list[0]


# Report total
str(aln_best_pair_DF.shape[0]) + " validly aligning hybrid read segment pairs where unaligned segments are of '" + recombType + " type', \n" + \
    "based on the sequence of accession-specific, chromosome-specific k-mers"

# Filter to retain hybrid read segments pairs where each segment aligns to the same chromosome
aln_best_pair_hom_DF = aln_best_pair_DF.loc[aln_best_pair_DF["acc1_tname"] == aln_best_pair_DF["acc2_tname"]]
str(aln_best_pair_hom_DF.shape[0]) + " '" + recombType + "-type' hybrid read segments align to the same chromosome"
str( round( aln_best_pair_hom_DF.shape[0] / aln_best_pair_DF.shape[0], 4 ) * 100 ) + "% of '" + recombType + "-type' hybrid read segment pairs align to the same chromosome"

# Filter to retain hybrid read segments pairs where the per-accession read segments align to within
# 2 * the given read length of each other in the same reference assembly
aln_dist_acc1_tstart_acc2_tstart = list(abs(aln_best_pair_hom_DF["acc1_tstart"] - aln_best_pair_hom_DF["acc2_tstart"]) + 1)
aln_dist_acc1_tstart_acc2_tend = list(abs(aln_best_pair_hom_DF["acc1_tstart"] - aln_best_pair_hom_DF["acc2_tend"]) + 1)
aln_dist_acc1_tend_acc2_tstart = list(abs(aln_best_pair_hom_DF["acc1_tend"] - aln_best_pair_hom_DF["acc2_tstart"]) + 1)
aln_dist_acc1_tend_acc2_tend = list(abs(aln_best_pair_hom_DF["acc1_tend"] - aln_best_pair_hom_DF["acc2_tend"]) + 1)

aln_dist_nparray = np.array([aln_dist_acc1_tstart_acc2_tstart,
                             aln_dist_acc1_tstart_acc2_tend,
                             aln_dist_acc1_tend_acc2_tstart,
                             aln_dist_acc1_tend_acc2_tend])
aln_dist_min = list(aln_dist_nparray.min(axis=0))
aln_dist_max = list(aln_dist_nparray.max(axis=0))
del aln_dist_nparray
gc.collect()

## Skin and overcook a cat
#aln_dist_tuple_list = list(zip(aln_dist_acc1_tstart_acc2_tstart,
#                               aln_dist_acc1_tstart_acc2_tend,
#                               aln_dist_acc1_tend_acc2_tstart,
#                               aln_dist_acc1_tend_acc2_tend))
#aln_dist_min = list(map(min, aln_dist_tuple_list))
#aln_dist_max = list(map(max, aln_dist_tuple_list))

aln_coords_nparray = np.array([list(aln_best_pair_hom_DF["acc1_tstart"]),
                               list(aln_best_pair_hom_DF["acc1_tend"]),
                               list(aln_best_pair_hom_DF["acc2_tstart"]),
                               list(aln_best_pair_hom_DF["acc2_tend"])])
#aln_coords_nparray.sort(axis=0)
sidx = aln_coords_nparray.argsort(axis=0)
aln_coords_nparray_sort = aln_coords_nparray[sidx, np.arange(sidx.shape[1])]

# Define recombination interval as the inner boundaries of segment alignment coordinates
event_start = list(aln_coords_nparray_sort[1])
event_end = list(aln_coords_nparray_sort[2])
event_start_end_nparray = np.array([event_start, event_end])
event_midpoint = (event_start_end_nparray[0] + event_start_end_nparray[1]) / 2
event_midpoint = list(event_midpoint.round().astype(int))
del aln_coords_nparray, sidx, aln_coords_nparray_sort, event_start_end_nparray
gc.collect()

aln_best_pair_hom_DF_cp = aln_best_pair_hom_DF.copy()
del aln_best_pair_hom_DF
gc.collect()
aln_best_pair_hom_DF = aln_best_pair_hom_DF_cp
#aln_best_pair_hom_DF.reset_index(drop=True, inplace=True)
aln_best_pair_hom_DF["aln_dist_min"] = aln_dist_min
aln_best_pair_hom_DF["aln_dist_max"] = aln_dist_max
aln_best_pair_hom_DF["event_start"] = event_start
aln_best_pair_hom_DF["event_end"] = event_end
aln_best_pair_hom_DF["event_midpoint"] = event_midpoint


# Get read lengths for hybrid read IDs in aln_best_pair_hom_DF$qname
# to be used for retaining alignment pairs where the Col and Ler
# read segments align to within a given distance of each other
# (e.g., the given hybrid read length) in the same assembly
hybrid_read_lengths_DF_list = []
for x in range(0, len(chrom)):
    reads_fa = "fasta/" + readsPrefix + \
        "_match_" + acc1 + "_" + region + "_" + chrom[x] + \
        "_specific_k" + str(kmerSize) + "_downsampled_op" + str(overlapProp) + "_hits" + str(minHits) + \
        "_match_" + acc2 + "_" + region + "_" + chrom[x] + \
        "_specific_k" + str(kmerSize) + "_downsampled_op" + str(overlapProp) + "_hits" + str(minHits) + \
        ".fa"
    reads_dict = SeqIO.index(reads_fa, "fasta")
    reads = [v for i, v in enumerate(reads_dict.values()) if v.id in list(aln_best_pair_hom_DF["qname"])]
    #reads_iter = SeqIO.parse(reads_fa, "fasta")
    #reads2 = [x for x in reads_iter if x.id in list(aln_best_pair_hom_DF["qname"])]
    read_ids = [x.id for x in reads]
    read_lens = [len(x) for x in reads]
    hybrid_read_lengths_DF = pd.DataFrame({ "read_id": read_ids,
                                            "read_len": read_lens })
    hybrid_read_lengths_DF_list.append(hybrid_read_lengths_DF) 

del hybrid_read_lengths_DF
gc.collect()

# Concatenate list elements (per-chromosome pd.DataFrames) into one pd.DataFrame
if len(hybrid_read_lengths_DF_list) > 1:
    hybrid_read_lengths_DF = pd.concat(objs=hybrid_read_lengths_DF_list, axis=0, ignore_index=True)
else:
    hybrid_read_lengths_DF = hybrid_read_lengths_DF_list[0]


aln_best_pair_hom_DF_read_lens = pd.merge(left=aln_best_pair_hom_DF, right=hybrid_read_lengths_DF,
                                          how="inner", left_on="qname", right_on="read_id")
del aln_best_pair_hom_DF
gc.collect()
aln_best_pair_hom_DF = aln_best_pair_hom_DF_read_lens.drop(columns="read_id")

aln_best_pair_hom_maxDist_DF = aln_best_pair_hom_DF.loc[aln_best_pair_hom_DF["aln_dist_min"] <= aln_best_pair_hom_DF["read_len"] * 2]


# Filter to retain putative recombination events between homologous chromosomes where
# the per-accession read segments align to within maxDist bp of each other in the same reference assembly, AND
# where the per-accession alignment length is >= alenTOqlen of the segment length
aln_best_pair_hom_maxDist_alenTOqlen_DF = aln_best_pair_hom_maxDist_DF.loc[ ( aln_best_pair_hom_maxDist_DF["acc1_alen"] / \
                                                                              aln_best_pair_hom_maxDist_DF["acc1_qlen"] >= alenTOqlen ) & \
                                                                            ( aln_best_pair_hom_maxDist_DF["acc2_alen"] / \
                                                                              aln_best_pair_hom_maxDist_DF["acc2_qlen"] >= alenTOqlen ) ] 


# Write to TSV
aln_best_pair_hom_maxDist_DF_filename = outdir + "/" + readsPrefix + \
    "_" + acc1 + "_" + acc2 + \
    "_k" + str(kmerSize) + "_op" + str(overlapProp) + "_h" + str(minHits) + \
    "_hom_maxDist_" + recombType + \
    "_alnTo_" + alnTo + "_" + \
    re.sub(",", "_", chrom) + ".tsv"
aln_best_pair_hom_maxDist_DF.to_csv(aln_best_pair_hom_maxDist_DF_filename, sep="\t", header=True, index=False)

aln_best_pair_hom_maxDist_alenTOqlen_DF_filename = outdir + "/" + readsPrefix + \
    "_" + acc1 + "_" + acc2 + \
    "_k" + str(kmerSize) + "_op" + str(overlapProp) + "_h" + str(minHits) + \
    "_hom_maxDist_aTOq" + str(alenTOqlen) + "_" + recombType + \
    "_alnTo_" + alnTo + "_" + \
    re.sub(",", "_", chrom) + ".tsv"
aln_best_pair_hom_maxDist_alenTOqlen_DF.to_csv(aln_best_pair_hom_maxDist_alenTOqlen_DF_filename, sep="\t", header=True, index=False)


# Write filtered hybrid reads to FASTA
for x in range(0, len(chrom)):
    # Input FASTA path to all hybrid reads for the given chromosome
    reads_fa = "fasta/" + readsPrefix + \
        "_match_" + acc1 + "_" + region + "_" + chrom[x] + \
        "_specific_k" + str(kmerSize) + "_downsampled_op" + str(overlapProp) + "_hits" + str(minHits) + \
        "_match_" + acc2 + "_" + region + "_" + chrom[x] + \
        "_specific_k" + str(kmerSize) + "_downsampled_op" + str(overlapProp) + "_hits" + str(minHits) + \
        ".fa"
    reads_dict = SeqIO.index(reads_fa, "fasta")
    # Get the reads where the segments align to within the given distance of each other
    reads_maxDist = [v for i, v in enumerate(reads_dict.values()) if v.id in list(aln_best_pair_hom_maxDist_DF["qname"])]
    # Get the reads where the segments align to within the given distance of each other, and
    # where the alen / qlen ratio >= alenTOqlen
    reads_maxDist_alenTOqlen = [v for i, v in enumerate(reads_dict.values()) if v.id in list(aln_best_pair_hom_maxDist_alenTOqlen_DF["qname"])]
    # Output FASTA path to maxDist-filtered hybrid reads for the given chromosome
    reads_maxDist_fa = outdir + "/" + readsPrefix + \
        "_" + acc1 + "_" + acc2 + \
        "_k" + str(kmerSize) + "_op" + str(overlapProp) + "_h" + str(minHits) + \
        "_hom_maxDist_" + recombType + \
        "_alnTo_" + alnTo + "_" + \
        chrom[x] + "_hybrid_reads.fa"
    # Output FASTA path to maxDist_alenTOqlen-filtered hybrid reads for the given chromosome
    reads_maxDist_alenTOqlen_fa = outdir + "/" + readsPrefix + \
        "_" + acc1 + "_" + acc2 + \
        "_k" + str(kmerSize) + "_op" + str(overlapProp) + "_h" + str(minHits) + \
        "_hom_maxDist_aTOq" + str(alenTOqlen) + "_" + recombType + \
        "_alnTo_" + alnTo + "_" + \
        chrom[x] + "_hybrid_reads.fa"
    # Write outputs
    with open(reads_maxDist_fa, "w") as reads_maxDist_fa_handle:
        SeqIO.write(reads_maxDist, reads_maxDist_fa_handle, "fasta")
    with open(reads_maxDist_alenTOqlen_fa, "w") as reads_maxDist_alenTOqlen_fa_handle:
        SeqIO.write(reads_maxDist_alenTOqlen, reads_maxDist_alenTOqlen_fa_handle, "fasta")

