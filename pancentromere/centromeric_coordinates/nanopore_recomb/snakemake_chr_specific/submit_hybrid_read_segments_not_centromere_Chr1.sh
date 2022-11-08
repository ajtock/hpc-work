#!/bin/bash

source ~/.bashrc

conda activate python_3.9.6
echo $(which snakemake)
echo $(which python)
echo $(which minimap2)
snakemake --configfile config_not_centromere_Chr1.yaml --profile profile/ --quiet rules
conda deactivate
