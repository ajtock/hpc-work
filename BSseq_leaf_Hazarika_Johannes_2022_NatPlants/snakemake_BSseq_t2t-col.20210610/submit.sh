#!/bin/bash

source ~/.bashrc

conda activate BSseq_mapping
echo $(which snakemake)
echo $(which bismark)
echo $(which R)
snakemake -p --profile profile/
conda deactivate
