#!/bin/bash


# Usage:
# ./acc_specific_kmers_bbduk.sh 512g 76 Col_ler_f1_pollen_500bp_minq99 Col-0.ragtag_scaffolds_centromere_specific_k 178

#MEMORY=512g
#THREADS=76
#FQ_PREFIX=Col_ler_f1_pollen_500bp_minq99
#FA_PREFIX=Col-0.ragtag_scaffolds_centromere_specific_k
#K=178

MEMORY=$1
THREADS=$2
FQ_PREFIX=$3
FA_PREFIX=$4
K=$5

echo ${MEMORY}
echo ${THREADS}
echo ${FQ_PREFIX}
echo ${FA_PREFIX}
echo ${K}
echo ${FA_PREFIX}${K} 

echo $(which bbduk.sh)

bbduk.sh -Xmx${MEMORY} \
         in=fastq/${FQ_PREFIX}.fq \
         outmatch=fastq/${FQ_PREFIX}_match_${FA_PREFIX}${K}.fq \
         outnonmatch=fastq/${FQ_PREFIX}_nonmatch_${FA_PREFIX}${K}.fq \
         k=${K} \
         ref=${FA_PREFIX}${K}.fa \
         threads=${THREADS} &> ${FQ_PREFIX}_match_${FA_PREFIX}${K}.log
