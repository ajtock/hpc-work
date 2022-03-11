#!/bin/bash

source ~/.bashrc

conda activate sRNAseq_mapping
echo $(which bowtie)
echo $(which R)
snakemake -p --cores 32
conda deactivate
