#!/bin/bash

source ~/.bashrc
conda activate python_3.9.6

./kmers_in_fasta.py -f Col-0.ragtag_scaffolds_centromere_Chr1.fa -k 24
./kmers_in_fasta.py -f Ler-0_110x.ragtag_scaffolds_centromere_Chr1.fa -k 24
./kmers_in_fasta.py -f Col-0.ragtag_scaffolds_not_centromere_Chr1.fa -k 24
./kmers_in_fasta.py -f Ler-0_110x.ragtag_scaffolds_not_centromere_Chr1.fa -k 24

./kmers_in_fasta.py -f Col-0.ragtag_scaffolds_centromere_Chr2.fa -k 24
./kmers_in_fasta.py -f Ler-0_110x.ragtag_scaffolds_centromere_Chr2.fa -k 24
./kmers_in_fasta.py -f Col-0.ragtag_scaffolds_not_centromere_Chr2.fa -k 24
./kmers_in_fasta.py -f Ler-0_110x.ragtag_scaffolds_not_centromere_Chr2.fa -k 24

./kmers_in_fasta.py -f Col-0.ragtag_scaffolds_centromere_Chr3.fa -k 24
./kmers_in_fasta.py -f Ler-0_110x.ragtag_scaffolds_centromere_Chr3.fa -k 24
./kmers_in_fasta.py -f Col-0.ragtag_scaffolds_not_centromere_Chr3.fa -k 24
./kmers_in_fasta.py -f Ler-0_110x.ragtag_scaffolds_not_centromere_Chr3.fa -k 24

./kmers_in_fasta.py -f Col-0.ragtag_scaffolds_centromere_Chr4.fa -k 24
./kmers_in_fasta.py -f Ler-0_110x.ragtag_scaffolds_centromere_Chr4.fa -k 24
./kmers_in_fasta.py -f Col-0.ragtag_scaffolds_not_centromere_Chr4.fa -k 24
./kmers_in_fasta.py -f Ler-0_110x.ragtag_scaffolds_not_centromere_Chr4.fa -k 24

./kmers_in_fasta.py -f Col-0.ragtag_scaffolds_centromere_Chr5.fa -k 24
./kmers_in_fasta.py -f Ler-0_110x.ragtag_scaffolds_centromere_Chr5.fa -k 24
./kmers_in_fasta.py -f Col-0.ragtag_scaffolds_not_centromere_Chr5.fa -k 24
./kmers_in_fasta.py -f Ler-0_110x.ragtag_scaffolds_not_centromere_Chr5.fa -k 24

./kmers_in_fasta.py -f t2t-col.20210610_ChrM_ChrC.fa -k 24

conda deactivate
