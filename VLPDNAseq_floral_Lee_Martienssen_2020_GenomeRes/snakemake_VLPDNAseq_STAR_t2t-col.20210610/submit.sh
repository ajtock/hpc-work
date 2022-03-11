#!/bin/bash

source ~/.bashrc

conda activate RNAseq_mapping
echo $(which snakemake)
echo $(which STAR)
echo $(which R)
snakemake -p --profile profile/
conda deactivate
