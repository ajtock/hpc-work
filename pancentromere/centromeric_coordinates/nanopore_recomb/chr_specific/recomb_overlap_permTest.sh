#!/bin/bash

source ~/.bashrc

conda activate python_3.9.6
./recomb_overlap_permTest.R ColLerF1pollen_1000bp_minq90 Col-0.ragtag_scaffolds Ler-0_110x.ragtag_scaffolds 24 0.9 11 Col-0.ragtag_scaffolds_Chr 0.9 not_centromere 'Chr1,Chr2,Chr3,Chr4,Chr5' 10000 COs
echo "Done COs"
./recomb_overlap_permTest.R ColLerF1pollen_1000bp_minq90 Col-0.ragtag_scaffolds Ler-0_110x.ragtag_scaffolds 24 0.9 11 Col-0.ragtag_scaffolds_Chr 0.9 not_centromere 'Chr1,Chr2,Chr3,Chr4,Chr5' 10000 NCOs
echo "Done NCOs"
./recomb_overlap_permTest.R ColLerF1pollen_1000bp_minq90 Col-0.ragtag_scaffolds Ler-0_110x.ragtag_scaffolds 24 0.9 11 Col-0.ragtag_scaffolds_Chr 0.9 not_centromere 'Chr1,Chr2,Chr3,Chr4,Chr5' 10000 COs_and_NCOs
echo "Done COs_and_NCOs"
conda deactivate
