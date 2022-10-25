#!/bin/bash

# Extract nucleotide sequence at given coordinates (specified in a BED file)
# in the genome in fasta format

# Note that BED files use 0-based half-open coordinates;
# start coordinates are 0-based and end coordinates are 1-based,
# such that a feature START coordinate that cooresponds to the first base in
# a chromosome is numbered 0, while a feature END coordinate that corresponds
# to the second base in a chromosome is numbered 2 in the BED file:
# see https://genome.ucsc.edu/FAQ/FAQformat.html#format1

./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Col-0.ragtag_scaffolds.fa Col-0.ragtag_scaffolds_centromere_Chr1 
./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Ler-0_110x.ragtag_scaffolds.fa Ler-0_110x.ragtag_scaffolds_centromere_Chr1
./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Col-0.ragtag_scaffolds.fa Col-0.ragtag_scaffolds_not_centromere_Chr1 
./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Ler-0_110x.ragtag_scaffolds.fa Ler-0_110x.ragtag_scaffolds_not_centromere_Chr1

./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Col-0.ragtag_scaffolds.fa Col-0.ragtag_scaffolds_centromere_Chr2 
./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Ler-0_110x.ragtag_scaffolds.fa Ler-0_110x.ragtag_scaffolds_centromere_Chr2
./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Col-0.ragtag_scaffolds.fa Col-0.ragtag_scaffolds_not_centromere_Chr2 
./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Ler-0_110x.ragtag_scaffolds.fa Ler-0_110x.ragtag_scaffolds_not_centromere_Chr2

./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Col-0.ragtag_scaffolds.fa Col-0.ragtag_scaffolds_centromere_Chr3 
./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Ler-0_110x.ragtag_scaffolds.fa Ler-0_110x.ragtag_scaffolds_centromere_Chr3
./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Col-0.ragtag_scaffolds.fa Col-0.ragtag_scaffolds_not_centromere_Chr3 
./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Ler-0_110x.ragtag_scaffolds.fa Ler-0_110x.ragtag_scaffolds_not_centromere_Chr3

./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Col-0.ragtag_scaffolds.fa Col-0.ragtag_scaffolds_centromere_Chr4 
./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Ler-0_110x.ragtag_scaffolds.fa Ler-0_110x.ragtag_scaffolds_centromere_Chr4
./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Col-0.ragtag_scaffolds.fa Col-0.ragtag_scaffolds_not_centromere_Chr4 
./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Ler-0_110x.ragtag_scaffolds.fa Ler-0_110x.ragtag_scaffolds_not_centromere_Chr4

./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Col-0.ragtag_scaffolds.fa Col-0.ragtag_scaffolds_centromere_Chr5 
./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Ler-0_110x.ragtag_scaffolds.fa Ler-0_110x.ragtag_scaffolds_centromere_Chr5
./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Col-0.ragtag_scaffolds.fa Col-0.ragtag_scaffolds_not_centromere_Chr5 
./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/Ler-0_110x.ragtag_scaffolds.fa Ler-0_110x.ragtag_scaffolds_not_centromere_Chr5

./bedtools_getfasta.sh /home/ajt200/rds/hpc-work/pancentromere/assemblies/t2t-col.20210610.fasta t2t-col.20210610_ChrM_ChrC 
