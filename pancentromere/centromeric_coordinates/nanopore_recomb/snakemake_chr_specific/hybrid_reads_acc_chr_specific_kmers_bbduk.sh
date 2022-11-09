#!/bin/bash

# From an ONT FASTQ file, get reads that contain subsequences matching
# query k-mers given in a FASTA file
# (e.g., Col-0-Chr1-specific centromeric k-mers, the output of acc_chr_specific_kmers.py)
 
# Usage (requires a lot of RAM and CPUs so
# sbatch hybrid_reads_acc_chr_specific_kmers_bbduk_icelake_slurm
# sbatch hybrid_reads_acc1_AND_acc2_chr_specific_kmers_bbduk_icelake_slurm ):
# conda activate python_3.9.6
# ./hybrid_reads_acc_chr_specific_kmers_bbduk.sh 48g 16 ColLerF1pollen_1000bp_minq90 Col-0.ragtag_scaffolds_not_centromere_Chr1_specific_k 24 0.9 11 not_centromere Chr1
# ./hybrid_reads_acc_chr_specific_kmers_bbduk.sh 24g 8 ColLerF1pollen_1000bp_minq90_match_Col-0.ragtag_scaffolds_not_centromere_Chr1_specific_k24_downsampled_op0.9_hits11 Ler-0_110x.ragtag_scaffolds_not_centromere_Chr1_specific_k 24 0.9 11 not_centromere Chr1
# conda deactivate

#MEMORY=48g
#THREADS=16
#FQ_PREFIX=ColLerF1pollen_1000bp_minq90
#FA_PREFIX=Col-0.ragtag_scaffolds_not_centromere_Chr1_specific_k
#K=24
#OP=0.9
#HITS=11
#REGION=not_centromere
#CHROM=Chr1

MEMORY=$1
THREADS=$2
FQ_PREFIX=$3
FA_PREFIX=$4
K=$5
OP=$6
HITS=$7
REGION=$8
CHROM=$9

echo ${MEMORY}
echo ${THREADS}
echo ${FQ_PREFIX}
echo ${FA_PREFIX}
echo ${K}
echo ${FA_PREFIX}${K} 
echo ${OP}
echo ${HITS}
echo ${REGION}
echo ${CHROM}

echo $(which bbduk.sh)

[ -d ${REGION}/${CHROM}/fasta/ ] || mkdir -p ${REGION}/${CHROM}/fasta/
[ -d ${REGION}/${CHROM}/kmer_loc_tsv/ ] || mkdir -p ${REGION}/${CHROM}/kmer_loc_tsv/
[ -d ${REGION}/${CHROM}/segments/Col-0/co/ ] || mkdir -p ${REGION}/${CHROM}/segments/Col-0/co/
[ -d ${REGION}/${CHROM}/segments/Col-0/nco/ ] || mkdir -p ${REGION}/${CHROM}/segments/Col-0/nco/
[ -d ${REGION}/${CHROM}/segments/Ler-0/co/ ] || mkdir -p ${REGION}/${CHROM}/segments/Ler-0/co/
[ -d ${REGION}/${CHROM}/segments/Ler-0/nco/ ] || mkdir -p ${REGION}/${CHROM}/segments/Ler-0/nco/
[ -d logs/segments/${REGION}/${CHROM}/ ] || mkdir -p logs/segments/${REGION}/${CHROM}/

bbduk.sh -Xmx${MEMORY} \
         in=fastq/${FQ_PREFIX}.fq \
         outmatch=fastq/${FQ_PREFIX}_match_${FA_PREFIX}${K}_downsampled_op${OP}_hits${HITS}.fq \
         outnonmatch=fastq/${FQ_PREFIX}_nonmatch_${FA_PREFIX}${K}_downsampled_op${OP}_hits${HITS}.fq \
         k=${K} \
         ref=fasta/${FA_PREFIX}${K}_downsampled_op${OP}.fa \
         maskmiddle=f \
         minkmerhits=${HITS} \
         threads=${THREADS} &> logs/${FQ_PREFIX}_match_${FA_PREFIX}${K}_downsampled_op${OP}_hits${HITS}_bbduk.log

reformat.sh in=fastq/${FQ_PREFIX}_match_${FA_PREFIX}${K}_downsampled_op${OP}_hits${HITS}.fq \
            out=fasta/${FQ_PREFIX}_match_${FA_PREFIX}${K}_downsampled_op${OP}_hits${HITS}.fa &> logs/${FQ_PREFIX}_match_${FA_PREFIX}${K}_downsampled_op${OP}_hits${HITS}_reformat.log

if [[ fasta/${FQ_PREFIX}_match_${FA_PREFIX}${K}_downsampled_op${OP}_hits${HITS}.fa == *"scaffolds_${REGION}_${CHROM}_"* ]]; then
    cd ${REGION}/${CHROM}/fasta/
    ln -s ../../../fasta/${FQ_PREFIX}_match_${FA_PREFIX}${K}_downsampled_op${OP}_hits${HITS}.fa .
    cd  ../../../
fi
