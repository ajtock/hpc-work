#!/bin/bash

source ~/.bashrc

conda activate sratools

# VLP DNA-seq Col_0 Rep1
fasterq-dump SRR8792541 -O ./
#mv SRR8792541_1.fastq.gz Col_0_VLPDNAseqs_Rep1_SRR8792541_R1.fastq.gz

conda deactivate
