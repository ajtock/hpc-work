#!/bin/bash

source ~/.bashrc

conda activate python_3.9.6

./get_chr_specific_segment_pairs_alnToSame.py -r ColLerF1pollen_1000bp_minq90 -a1 Col-0.ragtag_scaffolds -a2 Ler-0_110x.ragtag_scaffolds -k 24 -op 0.9 -mh 11 -at Col-0.ragtag_scaffolds_Chr -aq 0.90 -rt co -reg not_centromere -c 'Chr1'
./get_chr_specific_segment_pairs_alnToSame.py -r ColLerF1pollen_1000bp_minq90 -a1 Col-0.ragtag_scaffolds -a2 Ler-0_110x.ragtag_scaffolds -k 24 -op 0.9 -mh 11 -at Col-0.ragtag_scaffolds_Chr -aq 0.90 -rt nco -reg not_centromere -c 'Chr1'

conda deactivate
