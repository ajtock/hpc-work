#!/bin/bash

source ~/.bashrc

conda activate python_3.9.6

# COs
./get_chr_specific_segment_pairs_alnToDifferent.py -r ColLerF1pollen_1000bp_minq90 -a1 Col-0.ragtag_scaffolds -a2 Ler-0_110x.ragtag_scaffolds -k 24 -op 0.9 -mh 11 -mmq 2 -aq 0.90 -nq 0.90 -rt co -reg centromere -c 'Chr1'
./get_chr_specific_segment_pairs_alnToDifferent.py -r ColLerF1pollen_1000bp_minq90 -a1 Col-0.ragtag_scaffolds -a2 Ler-0_110x.ragtag_scaffolds -k 24 -op 0.9 -mh 11 -mmq 2 -aq 0.90 -nq 0.90 -rt co -reg centromere -c 'Chr2'
./get_chr_specific_segment_pairs_alnToDifferent.py -r ColLerF1pollen_1000bp_minq90 -a1 Col-0.ragtag_scaffolds -a2 Ler-0_110x.ragtag_scaffolds -k 24 -op 0.9 -mh 11 -mmq 2 -aq 0.90 -nq 0.90 -rt co -reg centromere -c 'Chr3'
./get_chr_specific_segment_pairs_alnToDifferent.py -r ColLerF1pollen_1000bp_minq90 -a1 Col-0.ragtag_scaffolds -a2 Ler-0_110x.ragtag_scaffolds -k 24 -op 0.9 -mh 11 -mmq 2 -aq 0.90 -nq 0.90 -rt co -reg centromere -c 'Chr4'
./get_chr_specific_segment_pairs_alnToDifferent.py -r ColLerF1pollen_1000bp_minq90 -a1 Col-0.ragtag_scaffolds -a2 Ler-0_110x.ragtag_scaffolds -k 24 -op 0.9 -mh 11 -mmq 2 -aq 0.90 -nq 0.90 -rt co -reg centromere -c 'Chr5'

# NCOs
./get_chr_specific_segment_pairs_alnToDifferent.py -r ColLerF1pollen_1000bp_minq90 -a1 Col-0.ragtag_scaffolds -a2 Ler-0_110x.ragtag_scaffolds -k 24 -op 0.9 -mh 11 -mmq 2 -aq 0.90 -nq 0.90 -rt nco -reg centromere -c 'Chr1'
./get_chr_specific_segment_pairs_alnToDifferent.py -r ColLerF1pollen_1000bp_minq90 -a1 Col-0.ragtag_scaffolds -a2 Ler-0_110x.ragtag_scaffolds -k 24 -op 0.9 -mh 11 -mmq 2 -aq 0.90 -nq 0.90 -rt nco -reg centromere -c 'Chr2'
./get_chr_specific_segment_pairs_alnToDifferent.py -r ColLerF1pollen_1000bp_minq90 -a1 Col-0.ragtag_scaffolds -a2 Ler-0_110x.ragtag_scaffolds -k 24 -op 0.9 -mh 11 -mmq 2 -aq 0.90 -nq 0.90 -rt nco -reg centromere -c 'Chr3'
./get_chr_specific_segment_pairs_alnToDifferent.py -r ColLerF1pollen_1000bp_minq90 -a1 Col-0.ragtag_scaffolds -a2 Ler-0_110x.ragtag_scaffolds -k 24 -op 0.9 -mh 11 -mmq 2 -aq 0.90 -nq 0.90 -rt nco -reg centromere -c 'Chr4'
./get_chr_specific_segment_pairs_alnToDifferent.py -r ColLerF1pollen_1000bp_minq90 -a1 Col-0.ragtag_scaffolds -a2 Ler-0_110x.ragtag_scaffolds -k 24 -op 0.9 -mh 11 -mmq 2 -aq 0.90 -nq 0.90 -rt nco -reg centromere -c 'Chr5'

conda deactivate
