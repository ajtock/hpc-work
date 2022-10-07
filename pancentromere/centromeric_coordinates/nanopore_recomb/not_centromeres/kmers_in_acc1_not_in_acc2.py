#!/usr/bin/env python3

# Author: Andy Tock (ajt200@cam.ac.uk)
# Date: 08/08/2022

# Usage:
# ./kmers_in_acc1_not_in_acc2.py -a1c Col-0.ragtag_scaffolds_centromeres -a2c Ler-0_110x.ragtag_scaffolds_centromeres -a1nc Col-0.ragtag_scaffolds_not_centromeres -a2nc Ler-0_110x.ragtag_scaffolds_not_centromeres -k 24 

# Write accession-specific centromeric k-mers
# to FASTA to supply as input file to bbduk.sh ("ref" parameter)
# containing query sequences to be matched against the FASTQ of ONT reads 

# ==== Import libraries
import sys
import os
import argparse
import screed
import pickle
import re
import subprocess

from matplotlib import pyplot as plt
from matplotlib.pyplot import subplots
from venn import venn
from time import time, sleep
import timeit


# ==== Capture user input as command-line arguments
# https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
def create_parser():
    parser = argparse.ArgumentParser(description="Fasta filename variables.")
    #### Define command-line arguments
    parser.add_argument("-a1c", "--acc1c", type=str, default="Col-0.ragtag_scaffolds_centromeres",
                        help="The prefix of the first accession's centromeric sequences. Default: Col-0.ragtag_scaffolds_centromeres")
    parser.add_argument("-a2c", "--acc2c", type=str, default="Ler-0_110x.ragtag_scaffolds_centromeres",
                        help="The prefix of the second accession's centromeric sequences. Default: Ler-0_110x.ragtag_scaffolds_centromeres")
    parser.add_argument("-a1nc", "--acc1nc", type=str, default="Col-0.ragtag_scaffolds_not_centromeres",
                        help="The prefix of the first accession's non-centromeric sequences. Default: Col-0.ragtag_scaffolds_not_centromeres")
    parser.add_argument("-a2nc", "--acc2nc", type=str, default="Ler-0_110x.ragtag_scaffolds_not_centromeres",
                        help="The prefix of the second accession's non-centromeric sequences. Default: Ler-0_110x.ragtag_scaffolds_not_centromeres")
    parser.add_argument("-k", "--kmerSize", type=int, default="24",
                        help="The size of the k-mers to be found and counted in the FASTA file.")
    #### Create parser
    return parser

parser = create_parser().parse_args()
print(parser)


outDir = "fasta"
plotDir = outDir + "/plots"

if not os.path.exists(outDir):
    os.makedirs(outDir)

if not os.path.exists(plotDir):
    os.makedirs(plotDir)


# Define a list containing the union of
# elements in an arbitrary number of lists 
def union_lists(*lists):
    """
    Get the union of elements in an arbitrary number of lists
    """
    return list(set.union(*map(set, lists)))


# Write dictionary of accession-specific centromeric k-mers
# to FASTA to supply to bbduk.sh as input k-mer database file
def write_fasta(kmer_dict, acc_name, outfile):
    """
    Write dictionary of k-mers to FASTA
    """
    with open(outfile, "w") as fa_object:
        for s in kmer_dict.keys():
            fa_object.write(">" + str(s) + "_" + acc_name + "\n")
            fa_object.write(kmer_dict[s] + "\n")


# Convert a k-mer into a hash
# See https://sourmash.readthedocs.io/en/latest/kmers-and-minhash.html
# But using hash() instead of mmh3.hash64() because since version 3.4,
# Python hash() uses SipHash, which is more secure and less vulnerable
# to hash collision attacks
def hash_kmer(kmer):
    """
    Convert a k-mer into a hash.
    """
    # Get the reverse complement
    rc_kmer = screed.rc(kmer)
    #
    # Get the lesser of kmer and rc_kmer, based on lexicographical order
    if kmer < rc_kmer:
        canonical_kmer = kmer
    else:
        canonical_kmer = rc_kmer
    #
    # Calculate hash
    kmer_hash = hash(canonical_kmer)
    #
    return kmer_hash

# Convert a list of k-mers into a list of hashes
# See https://sourmash.readthedocs.io/en/latest/kmers-and-minhash.html
def hash_kmers(kmers):
    """
    Convert a list of k-mers into a list of hashes.
    """
    kmer_hashes = []
    for kmer in kmers:
        kmer_hashes.append(hash_kmer(kmer))
    #
    return kmer_hashes


# Load k-mer count dictionaries (saved as pickle files)
with open("kmers_in_fasta_dict/" + parser.acc1c + ".fa_k" + str(parser.kmerSize) + ".pickle", "rb") as handle:
  acc1cen_dict = pickle.load(handle)

with open("kmers_in_fasta_dict/" + parser.acc2c + ".fa_k" + str(parser.kmerSize) + ".pickle", "rb") as handle:
  acc2cen_dict = pickle.load(handle)

with open("kmers_in_fasta_dict/" + parser.acc1nc + ".fa_k" + str(parser.kmerSize) + ".pickle", "rb") as handle:
  acc1notcen_dict = pickle.load(handle)

with open("kmers_in_fasta_dict/" + parser.acc2nc + ".fa_k" + str(parser.kmerSize) + ".pickle", "rb") as handle:
  acc2notcen_dict = pickle.load(handle)

with open("kmers_in_fasta_dict/t2t-col.20210610_ChrM_ChrC.fa_k" + str(parser.kmerSize) + ".pickle", "rb") as handle:
  mitochloro_dict = pickle.load(handle)


# Make keys (k-mers) in each dictionary a list
acc1cen = list(acc1cen_dict.keys())
acc2cen = list(acc2cen_dict.keys())
acc1notcen = list(acc1notcen_dict.keys())
acc2notcen = list(acc2notcen_dict.keys())
mitochloro = list(mitochloro_dict.keys())

print(len(acc1cen))
# 660531
print(len(acc2cen))
# 659434 
print(len(acc1notcen))
# 109988424
print(len(acc2notcen))
# 110021023
print(len(mitochloro))
# 475021

## Convert each list of k-mers into a list of hashes
#acc1cen_hashes = hash_kmers(acc1cen)
#acc2cen_hashes = hash_kmers(acc2cen)
#acc1notcen_hashes = hash_kmers(acc1notcen)
#acc2notcen_hashes = hash_kmers(acc2notcen)
#mitochloro_hashes = hash_kmers(mitochloro)


# Make venn to show k-mer overlap
cmap = "plasma"

#hash_dict = {
#    "Col-0 cen": set(acc1cen_hashes),
#    "Ler-0 cen": set(acc2cen_hashes),
#    "Col-0 arm": set(acc1notcen_hashes),
#    "Ler-0 arm": set(acc2notcen_hashes)
#} 

dataset_dict = {
    "Col-0 cen": set(acc1cen),
    "Ler-0 cen": set(acc2cen),
    "Col-0 arm": set(acc1notcen),
    "Ler-0 arm": set(acc2notcen)
} 


venn_plot = venn(dataset_dict, cmap="plasma")
venn_fig = venn_plot.get_figure()
venn_fig.savefig(plotDir + "/venn.png")

#hash_venn_plot = venn(hash_dict, cmap="plasma")
#hash_venn_fig = venn_plot.get_figure()
#hash_venn_fig.savefig(plotDir + "/venn_hash.png")


## acc1c_kmers (acc1 centromere-specific k-mers)
# Define union of acc1notcen, acc2cen, acc2notcen and mitochloro
acc1notcen_acc2cen_acc2notcen_mitochloro_union = union_lists(acc1notcen, acc2cen, acc2notcen, mitochloro)

# Get k-mers unique to acc1cen by computing difference between
# acc1cen and the above-defined union
# (i.e., acc1cen \ acc1notcen_acc2cen_acc2notcen_mitochloro_union): acc1c_kmers
acc1c_kmers = list(set(acc1cen).difference(set(acc1notcen_acc2cen_acc2notcen_mitochloro_union)))
acc1c_kmers = sorted(acc1c_kmers)
print(len(acc1c_kmers))
# 309849
print(len(acc1cen) - len(acc1c_kmers))
# 350682

# Convert into dictionary
acc1c_kmers_1tolen = list(range(1, len(acc1c_kmers)+1)) 
acc1c_kmers_dict = dict(zip(acc1c_kmers_1tolen, acc1c_kmers))


# acc2c_kmers (acc2 centromere-specific k-mers)
# Define union of acc2notcen, acc1cen, acc1notcen and mitochloro
acc2notcen_acc1cen_acc1notcen_mitochloro_union = union_lists(acc2notcen, acc1cen, acc1notcen, mitochloro)

# Get k-mers unique to acc2cen by computing difference between
# acc2cen and the second above-defined union
# (i.e., acc2cen \ acc2notcen_acc1cen_acc1notcen_mitochloro_union): acc2c_kmers
acc2c_kmers = list(set(acc2cen).difference(set(acc2notcen_acc1cen_acc1notcen_mitochloro_union)))
acc2c_kmers = sorted(acc2c_kmers)
print(len(acc2c_kmers))
# 353341 
print(len(acc2cen) - len(acc2c_kmers))
# 306093

# Convert into dictionary
acc2c_kmers_1tolen = list(range(1, len(acc2c_kmers)+1)) 
acc2c_kmers_dict = dict(zip(acc2c_kmers_1tolen, acc2c_kmers))


# acc1nc_kmers (acc1 non-centromere-specific k-mers)
# Define union of acc1cen, acc2cen, acc2notcen and mitochloro
acc1cen_acc2cen_acc2notcen_mitochloro_union = union_lists(acc1cen, acc2cen, acc2notcen, mitochloro)

# Get k-mers unique to acc1notcen by computing difference between
# acc1notcen and the above-defined union
# (i.e., acc1notcen \ acc1cen_acc2cen_acc2notcen_mitochloro_union): acc1nc_kmers
acc1nc_kmers = list(set(acc1notcen).difference(set(acc1cen_acc2cen_acc2notcen_mitochloro_union)))
acc1nc_kmers = sorted(acc1nc_kmers)
print(len(acc1nc_kmers))
# 19183610 
print(len(acc1notcen) - len(acc1nc_kmers))
# 90804814 

# Convert into dictionary
acc1nc_kmers_1tolen = list(range(1, len(acc1nc_kmers)+1)) 
acc1nc_kmers_dict = dict(zip(acc1nc_kmers_1tolen, acc1nc_kmers))


# acc2nc_kmers (acc2 non-centromere-specific k-mers)
# Define union of acc1cen, acc2cen, acc1notcen and mitochloro
acc1cen_acc2cen_acc1notcen_mitochloro_union = union_lists(acc1cen, acc2cen, acc1notcen, mitochloro)

# Get k-mers unique to acc2notcen by computing difference between
# acc2notcen and the above-defined union
# (i.e., acc2notcen \ acc1cen_acc2cen_acc1notcen_mitochloro_union): acc2nc_kmers
acc2nc_kmers = list(set(acc2notcen).difference(set(acc1cen_acc2cen_acc1notcen_mitochloro_union)))
acc2nc_kmers = sorted(acc2nc_kmers)
print(len(acc2nc_kmers))
# 19452914 
print(len(acc2notcen) - len(acc2nc_kmers))
# 90568109

# Convert into dictionary
acc2nc_kmers_1tolen = list(range(1, len(acc2nc_kmers)+1)) 
acc2nc_kmers_dict = dict(zip(acc2nc_kmers_1tolen, acc2nc_kmers))


# Write to FASTA
write_fasta(kmer_dict=acc1c_kmers_dict,
            acc_name=parser.acc1c[0:5],
            outfile=outDir + "/" + parser.acc1c + "_specific_k" + str(parser.kmerSize) + ".fa")

write_fasta(kmer_dict=acc2c_kmers_dict,
            acc_name=parser.acc2c[0:5],
            outfile=outDir + "/" + parser.acc2c + "_specific_k" + str(parser.kmerSize) + ".fa")

write_fasta(kmer_dict=acc1nc_kmers_dict,
            acc_name=parser.acc1nc[0:5],
            outfile=outDir + "/" + parser.acc1nc + "_specific_k" + str(parser.kmerSize) + ".fa")

write_fasta(kmer_dict=acc2nc_kmers_dict,
            acc_name=parser.acc2nc[0:5],
            outfile=outDir + "/" + parser.acc2nc + "_specific_k" + str(parser.kmerSize) + ".fa")


# Align accession-specific k-mers to respesctive genome
def align_kmers_bowtie(kmers_fasta, genome):
    #kmers_fasta = outDir + "/" + parser.acc1nc + "_specific_k" + str(parser.kmerSize) + ".fa"
    #genome = re.sub(r"(_scaffolds)_.+", r"\1", parser.acc1nc)
    """
    Align accession-specific kmers in FASTA format to respective genome using bowtie
    and convert SAM into sorted BAM.
    """
    out_sam_err = re.sub(".fa", "_bowtie.err", kmers_fasta)
    out_bam = re.sub(".fa", "_bowtie_sorted.bam", kmers_fasta)
    out_bam_err = re.sub(".fa", "_bowtie_sorted.err", kmers_fasta)
    aln_cmd = ["bowtie"] + \
              ["-x", "index/" + genome] + \
              ["-f", kmers_fasta] + \
              ["-v", "0"] + \
              ["-a", "--best", "--strata"] + \
              ["--sam"] + \
              ["--no-unal"] + \
              ["--threads", "32"]
    bam_cmd = ["samtools", "sort"] + \ 
              ["-@", "32"] + \
              ["-m", "3G"] + \
              ["-o", out_bam]
    with open(out_sam_err, "w") as out_sam_err_handle, open(out_bam_err, "w") as out_bam_err_handle:
        sam = subprocess.Popen(aln_cmd, stdout=subprocess.PIPE, stderr=out_sam_err_handle)
        subprocess.check_output(bam_cmd, stdin=sam.stdout, stderr=out_bam_err_handle)
        sam.wait()


# Convert BAM into BED and coordinate-sort BED
def bam_to_bed(in_bam):
in_bam = outDir + "/" + parser.acc1nc + "_specific_k" + str(parser.kmerSize) + "_bowtie_sorted.bam"
"""
Convert BAM into BED and coordinate-sort BED.
"""
out_bed = re.sub(".bam", ".bed", in_bam)
out_bed_err = re.sub(".bam", "_bam2bed.err", in_bam)
out_grep_err = re.sub(".bam", "_grepbed.err", in_bam)
out_sort_err = re.sub(".bam", "_sortbed.err", in_bam)
bed_cmd = ["bedtools", "bamtobed"] + \
          ["-i", in_bam] + \
          ["-tag", "NM"]
grep_cmd = ["grep", "^Chr"]
sort_cmd = ["sort", "-k1,1", "-k2,2n"]
#awk_cmd =  ["awk", "'BEGIN{FS="\t";OFS="\t"}", "{print", "$1,", "$2,", "$3}'"]
with open(out_bed, "w") as out_bed_handle, \
    open(out_bed_err, "w") as out_bed_err_handle, \
    open(out_grep_err, "w") as out_grep_err_handle, \
    open(out_sort_err, "w") as out_sort_err_handle:
    bed = subprocess.Popen(bed_cmd, stdout=subprocess.PIPE, stderr=out_bed_err_handle)
    grep = subprocess.Popen(grep_cmd, stdin=bed.stdout, stdout=subprocess.PIPE, stderr=out_grep_err_handle)
    subprocess.call(sort_cmd, stdin=grep.stdout, stdout=out_bed_handle, stderr=out_sort_err_handle)
    bed.wait()
    grep.wait()

 
with open(out_bed_err, "w") as out_bed_err_handle, open(out_sort_err, "w") as out_sort_err_handle:
    bed = subprocess.Popen(bed_cmd, stdout=subprocess.PIPE, stderr=out_bed_err_handle)
    subprocess.check_output(sort_cmd, stdin=bed.stdout, stderr=out_sort_err_handle)


    subprocess.run(bam_cmd, stderr=out_bam_err_handle)
    subprocess.run(bed_cmd, stdout=out_bed_handle, stderr=out_bed_err_handle)
    #subprocess.run(["rm", out_bam])

        subprocess.Popen(bed.

#    # Delete file(s) if unmapped
#    sam_view_cmd = ["samtools"] + \
#                   ["view", outsam]
#    sam_record = subprocess.Popen(sam_view_cmd, stdout=subprocess.PIPE)
#    sam_flag = int( subprocess.check_output(["cut", "-f2"], stdin=sam_record.stdout) )
#    sam_rm_cmd = ["rm"] + \
#                 [outsam, outerr]
#    if sam_flag == 4:
#        subprocess.run(sam_rm_cmd)




# Align accession-specific k-mers to respective genome
# Convert BAM into BED
# acc1c_kmers
align_kmers_bowtie(kmers_fasta=outDir + "/" + parser.acc1c + "_specific_k" + str(parser.kmerSize) + ".fa",
                   genome=re.sub(r"(_scaffolds)_.+", r"\1", parser.acc1c))
bam_to_bed(in_bam=outDir + "/" + parser.acc1c + "_specific_k" + str(parser.kmerSize) + "_bowtie_sorted.bam")
# acc2c_kmers
align_kmers_bowtie(kmers_fasta=outDir + "/" + parser.acc2c + "_specific_k" + str(parser.kmerSize) + ".fa",
                   genome=re.sub(r"(_scaffolds)_.+", r"\1", parser.acc2c))
bam_to_bed(in_bam=outDir + "/" + parser.acc2c + "_specific_k" + str(parser.kmerSize) + "_bowtie_sorted.bam")

# acc1nc_kmers
align_kmers_bowtie(kmers_fasta=outDir + "/" + parser.acc1nc + "_specific_k" + str(parser.kmerSize) + ".fa",
                   genome=re.sub(r"(_scaffolds)_.+", r"\1", parser.acc1nc))
bam_to_bed(in_bam=outDir + "/" + parser.acc1nc + "_specific_k" + str(parser.kmerSize) + "_bowtie_sorted.bam")
# acc2nc_kmers
align_kmers_bowtie(kmers_fasta=outDir + "/" + parser.acc2nc + "_specific_k" + str(parser.kmerSize) + ".fa",
                   genome=re.sub(r"(_scaffolds)_.+", r"\1", parser.acc2nc))
bam_to_bed(in_bam=outDir + "/" + parser.acc2nc + "_specific_k" + str(parser.kmerSize) + "_bowtie_sorted.bam")


