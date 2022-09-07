#!/bin/bash

./read_seqments_acc1_and_acc2_specific_kmers.py \
 -r Col_ler_f1_pollen_500bp_minq99 \
 -a1 Col-0.ragtag_scaffolds_centromeres \
 -a2 Ler-0_110x.ragtag_scaffolds_centromeres \
 -k 24 \
 -mh 3
