#!/bin/bash

source ~/.bashrc

conda activate BSseq_mapping

rename -v 's/7d_/D7_sRNA_Rep/' *.fastq.gz
rename -v 's/35d_/D35_sRNA_Rep/' *.fastq.gz
rename -v 's/.fastq/_R1.fastq/' *.fastq.gz

conda deactivate

