#!/bin/bash

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/001/SRR1042171/SRR1042171.fastq.gz
mv SRR1042171.fastq.gz Col_0_sRNA_SRR1042171_R1.fastq.gz 

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/002/SRR1042172/SRR1042172.fastq.gz
mv SRR1042172.fastq.gz ddm1_2_sRNA_SRR1042172_R1.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/003/SRR1042173/SRR1042173.fastq.gz
mv SRR1042173.fastq.gz rdr6_15_sRNA_SRR1042173_R1.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/004/SRR1042174/SRR1042174.fastq.gz
mv SRR1042174.fastq.gz ddm1_2_rdr6_15_sRNA_SRR1042174_R1.fastq.gz 
