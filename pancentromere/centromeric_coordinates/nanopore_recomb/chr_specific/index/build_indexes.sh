#!/bin/bash

# Pre-compute high-frequency k-mers (e.g., top 0.02% most frequent)
# in a reference genome using the meryl v1.3 k-mer counting tool,
# to enable whole-genome alignment of reads using winnowmap v2.03.
# Then build an index for the genome using winnowmap v2.03,
# and separately minimap2 v2.24 and bowtie v1.3.1

# Usage:
# (./build_indexes.sh Col-0.ragtag_scaffolds_Chr 15 32) &> Col-0.ragtag_scaffolds_Chr_indexes.log
# (./build_indexes.sh Ler-0_110x.ragtag_scaffolds_Chr 15 32) &> Ler-0_110x.ragtag_scaffolds_Chr_indexes.log

source ~/.bashrc
conda activate python_3.9.6

GENOME=$1
K=$2
THREADS=$3

echo $(which meryl)
echo $(which winnowmap)
echo $(which minimap2)
echo $(which bowtie-build)

/home/ajt200/miniconda3/envs/python_3.9.6/bin/meryl count \
  k=${K} \
  threads=${THREADS} \
  output ${GENOME}_merylDB \
  ${GENOME}.fa

/home/ajt200/miniconda3/envs/python_3.9.6/bin/meryl print greater-than distinct=0.9998 \
  ${GENOME}_merylDB > ${GENOME}_repetitive_k${K}.txt

# winnowmap map-ont index
/home/ajt200/miniconda3/envs/python_3.9.6/bin/winnowmap \
  -W ${GENOME}_repetitive_k${K}.txt \
  -x map-ont \
  -k ${K} \
  -d ${GENOME}__ont.wmi \
  -t ${THREADS} \
  ${GENOME}.fa

# minimap2 map-ont index
/home/ajt200/miniconda3/envs/python_3.9.6/bin/minimap2 \
  -x map-ont \
  -k ${K} \
  -d ${GENOME}__ont.mmi \
  -t ${THREADS} \
  ${GENOME}.fa

# minimap2 sr (short-read) index
/home/ajt200/miniconda3/envs/python_3.9.6/bin/minimap2 \
  -x sr \
  -k ${K} \
  -d ${GENOME}__sr.mmi \
  -t ${THREADS} \
  ${GENOME}.fa

# bowtie index
/home/ajt200/miniconda3/envs/python_3.9.6/bin/bowtie-build \
  --threads ${THREADS} \
  ${GENOME}.fa \
  ${GENOME}

conda deactivate
