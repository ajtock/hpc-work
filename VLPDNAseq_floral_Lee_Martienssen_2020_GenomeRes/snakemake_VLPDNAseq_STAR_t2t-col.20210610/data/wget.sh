#!/bin/bash

# Download raw short-read VLP DNA-seq FASTQ files from ENA

# Col_0 VLP DNA-seq Rep1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/001/SRR8792541/SRR8792541_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/001/SRR8792541/SRR8792541_2.fastq.gz
rename -v SRR8792541_ Col_0_VLPDNAseq_Rep1_SRR8792541_R SRR8792541_*.fastq.gz

# Col_0 VLP DNA-seq Rep2
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/002/SRR8792542/SRR8792542_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/002/SRR8792542/SRR8792542_2.fastq.gz
rename -v SRR8792542_ Col_0_VLPDNAseq_Rep2_SRR8792542_R SRR8792542_*.fastq.gz

# Col_0 VLP DNA-seq Rep3
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/003/SRR8792543/SRR8792543_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/003/SRR8792543/SRR8792543_2.fastq.gz
rename -v SRR8792543_ Col_0_VLPDNAseq_Rep3_SRR8792543_R SRR8792543_*.fastq.gz


# ddm1_2 VLP DNA-seq Rep1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/004/SRR8792544/SRR8792544_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/004/SRR8792544/SRR8792544_2.fastq.gz
rename -v SRR8792544_ ddm1_2_VLPDNAseq_Rep1_SRR8792544_R SRR8792544_*.fastq.gz

# ddm1_2 VLP DNA-seq Rep2
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/005/SRR8792545/SRR8792545_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/005/SRR8792545/SRR8792545_2.fastq.gz
rename -v SRR8792545_ ddm1_2_VLPDNAseq_Rep2_SRR8792545_R SRR8792545_*.fastq.gz


# ddm1_2_rdr6_15 VLP DNA-seq Rep1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/006/SRR8792546/SRR8792546_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/006/SRR8792546/SRR8792546_2.fastq.gz
rename -v SRR8792546_ ddm1_2_rdr6_15_VLPDNAseq_Rep1_SRR8792546_R SRR8792546_*.fastq.gz

# ddm1_2_rdr6_15 VLP DNA-seq Rep2
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/007/SRR8792547/SRR8792547_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/007/SRR8792547/SRR8792547_2.fastq.gz
rename -v SRR8792547_ ddm1_2_rdr6_15_VLPDNAseq_Rep2_SRR8792547_R SRR8792547_*.fastq.gz

# ddm1_2_rdr6_15 VLP DNA-seq Rep3
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/008/SRR8792548/SRR8792548_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/008/SRR8792548/SRR8792548_2.fastq.gz
rename -v SRR8792548_ ddm1_2_rdr6_15_VLPDNAseq_Rep3_SRR8792548_R SRR8792548_*.fastq.gz

