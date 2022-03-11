#!/bin/bash

for i in *.fa; do
  /home/ajt200/miniconda3/envs/ChIPseq_mapping/bin/samtools faidx ${i}
#cut -f1,2 ${i}.fa.fai > ${i}.fa.chrom.sizes
#ln -s ${i}.fa.chrom.sizes ${i}.fa.sizes 
done
