#!/bin/bash

source ~/.bashrc

conda activate RNAseq_mapping
echo $(which deeptools)
snakemake -p --profile profile/
conda deactivate
