#!/bin/bash

source activate R-4.1.2

./alphabeta_per_cytosine_MA1_2_dopar.R t2t-col.20210610 CpG 10000 10000 Chr1 48
./alphabeta_per_cytosine_MA1_2_dopar.R t2t-col.20210610 CpG 10000 10000 Chr2 48
./alphabeta_per_cytosine_MA1_2_dopar.R t2t-col.20210610 CpG 10000 10000 Chr3 48
./alphabeta_per_cytosine_MA1_2_dopar.R t2t-col.20210610 CpG 10000 10000 Chr4 48
./alphabeta_per_cytosine_MA1_2_dopar.R t2t-col.20210610 CpG 10000 10000 Chr5 48

conda deactivate
