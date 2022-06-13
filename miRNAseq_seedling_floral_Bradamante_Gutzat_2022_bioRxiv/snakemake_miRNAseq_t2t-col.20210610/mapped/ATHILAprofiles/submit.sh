#!/bin/bash

source ~/.bashrc

conda activate sRNAseq_mapping
echo $(which deeptools)
snakemake -p --profile profile/
conda deactivate
