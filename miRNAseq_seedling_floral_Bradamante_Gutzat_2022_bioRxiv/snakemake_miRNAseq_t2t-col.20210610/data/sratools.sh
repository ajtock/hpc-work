#!/bin/bash

source ~/.bashrc

conda activate sratools

# D7 (seedling shoot apicical meristem)

fasterq-dump SRX13507587 -O ./
mv SRX13507587.fastq.gz Col_0_AGO5_D7_sRNA_Rep1_SRX13507587_R1.fastq.gz

fasterq-dump SRX13507588 -O ./
mv SRX13507588.fastq.gz Col_0_AGO5_D7_sRNA_Rep2_SRX13507588_R1.fastq.gz

fasterq-dump SRX13507589 -O ./
mv SRX13507589.fastq.gz Col_0_AGO9_D7_sRNA_Rep1_SRX13507589_R1.fastq.gz

fasterq-dump SRX13507590 -O ./
mv SRX13507590.fastq.gz Col_0_AGO9_D7_sRNA_Rep2_SRX13507590_R1.fastq.gz

fasterq-dump SRX13507591 -O ./
mv SRX13507591.fastq.gz Col_0_input_D7_sRNA_Rep1_SRX13507591_R1.fastq.gz

fasterq-dump SRX13507592 -O ./
mv SRX13507592.fastq.gz Col_0_input_D7_sRNA_Rep2_SRX13507592_R1.fastq.gz


# D35 (floral meristem) 

fasterq-dump SRX13507595 -O ./
mv SRX13507595.fastq.gz Col_0_AGO5_D35_sRNA_Rep1_SRX13507595_R1.fastq.gz

fasterq-dump SRX13507596 -O ./
mv SRX13507596.fastq.gz Col_0_AGO5_D35_sRNA_Rep2_SRX13507596_R1.fastq.gz

fasterq-dump SRX13507597 -O ./
mv SRX13507597.fastq.gz Col_0_AGO9_D35_sRNA_Rep1_SRX13507597_R1.fastq.gz

fasterq-dump SRX13507598 -O ./
mv SRX13507598.fastq.gz Col_0_AGO9_D35_sRNA_Rep2_SRX13507598_R1.fastq.gz

fasterq-dump SRX13507599 -O ./
mv SRX13507599.fastq.gz Col_0_input_D35_sRNA_Rep1_SRX13507599_R1.fastq.gz

fasterq-dump SRX13507600 -O ./
mv SRX13507600.fastq.gz Col_0_input_D35_sRNA_Rep2_SRX13507600_R1.fastq.gz


conda deactivate
