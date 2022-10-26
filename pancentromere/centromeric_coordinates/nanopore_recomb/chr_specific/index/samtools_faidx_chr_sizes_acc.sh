#!/bin/bash

# Usage:
# ./samtools_faidx_chr_sizes_acc.sh Col-0.ragtag_scaffolds_Chr
# ./samtools_faidx_chr_sizes_acc.sh Ler-0_110x.ragtag_scaffolds_Chr

source ~/.bashrc
conda activate python_3.9.6

GENOME=$1

echo $(which samtools)

samtools faidx ${GENOME}.fa
cut -f1,2 ${GENOME}.fa.fai > ${GENOME}.fa.chrom.sizes
ln -s ${GENOME}.fa.chrom.sizes ${GENOME}.fa.sizes 
for CHROM in `cut -f 1 ${GENOME}.fa.sizes`;
do
  echo ${CHROM}
  grep "^${CHROM}" ${GENOME}.fa.sizes > ${GENOME}_${CHROM}.fa.sizes
done

conda deactivate
