#!/bin/bash

# From an ONT FASTQ file, get reads that contain subsequences matching
# query k-mers given in a FASTA file
# (e.g., accession-specific centromeric k-mers, the output of kmers_in_acc1_not_in_acc2.py)
 
# Usage:
# ./acc_specific_kmers_bbduk.sh 256g 76 Col_ler_f1_pollen_500bp_minq99 Col-0.ragtag_scaffolds_centromeres_specific_k 24 3
# ./acc_specific_kmers_bbduk.sh 256g 76 Col_ler_f1_pollen_500bp_minq99 Ler-0_110x.ragtag_scaffolds_centromeres_specific_k 24 3
# ./acc_specific_kmers_bbduk.sh 256g 76 Col_ler_f1_pollen_500bp_minq99_match_Col-0.ragtag_scaffolds_centromeres_specific_k24_hits3 Ler-0_110x.ragtag_scaffolds_centromeres_specific_k 24 3
# ./acc_specific_kmers_bbduk.sh 256g 76 Col_ler_f1_pollen_500bp_minq99_match_Ler-0_110x.ragtag_scaffolds_centromeres_specific_k24_hits3 Col-0.ragtag_scaffolds_centromeres_specific_k 24 3

#MEMORY=256g
#THREADS=76
#FQ_PREFIX=Col_ler_f1_pollen_500bp_minq99
#FA_PREFIX=Col-0.ragtag_scaffolds_centromeres_specific_k
#K=24
#HITS=3

MEMORY=$1
THREADS=$2
FQ_PREFIX=$3
FA_PREFIX=$4
K=$5
HITS=$6

echo ${MEMORY}
echo ${THREADS}
echo ${FQ_PREFIX}
echo ${FA_PREFIX}
echo ${K}
echo ${FA_PREFIX}${K} 
echo ${HITS}

echo $(which bbduk.sh)

bbduk.sh -Xmx${MEMORY} \
         in=fastq/${FQ_PREFIX}.fq \
         outmatch=fastq/${FQ_PREFIX}_match_${FA_PREFIX}${K}_hits${HITS}.fq \
         outnonmatch=fastq/${FQ_PREFIX}_nonmatch_${FA_PREFIX}${K}_hits${HITS}.fq \
         k=${K} \
         ref=fasta/${FA_PREFIX}${K}.fa \
         maskmiddle=f \
         minkmerhits=${HITS} \
         threads=${THREADS} &> ${FQ_PREFIX}_match_${FA_PREFIX}${K}_hits${HITS}_bbduk.log

reformat.sh in=fastq/${FQ_PREFIX}_match_${FA_PREFIX}${K}_hits${HITS}.fq \
            out=fasta/${FQ_PREFIX}_match_${FA_PREFIX}${K}_hits${HITS}.fa &> ${FQ_PREFIX}_match_${FA_PREFIX}${K}_hits${HITS}_reformat.log
