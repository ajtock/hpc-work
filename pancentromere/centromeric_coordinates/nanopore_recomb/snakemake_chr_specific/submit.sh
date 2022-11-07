#!/bin/bash

source ~/.bashrc

conda activate python_3.9.6
echo $(which snakemake)
echo $(which python)
echo $(which minimap2)
snakemake --printshellcmds --profile profile/ --dryrun
conda deactivate
