#!/bin/bash

sbatch centromere_Chr1_hybrid_read_segments_acc_chr_specific_kmers_py_icelake_array01of04_slurm && sleep 2;
sbatch centromere_Chr1_hybrid_read_segments_acc_chr_specific_kmers_py_icelake_array02of04_slurm && sleep 2;
sbatch centromere_Chr1_hybrid_read_segments_acc_chr_specific_kmers_py_icelake_array03of04_slurm && sleep 2;
sbatch centromere_Chr1_hybrid_read_segments_acc_chr_specific_kmers_py_icelake_array04of04_slurm && sleep 2;