#!/bin/bash

source ~/.bashrc

conda activate sRNAseq_mapping
echo $(which snakemake)
echo $(which bowtie)
echo $(which R)
snakemake -p --profile profile/
conda deactivate
