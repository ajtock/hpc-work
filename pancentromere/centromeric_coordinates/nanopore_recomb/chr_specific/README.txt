# ONT F1 hybrid reads recombination detection pipeline

1. Run get_chr_region_coordinates.R for each assembly and each chromosome

2. Run bedtools_getfasta.sh to get the fasta sequence based on BEDs of centromere and not_centromere coordinates (per accesion, per chromosome)

3. Run index/build_indexes.sh to generate genome indexes for various aligners

4. Run kmers_in_fasta.py on the per-accession, per-chromosome (and ChrM_ChrC), per-region FASTA files

5. Run acc_chr_specific_kmers.py on the pickle output files of kmers_in_fasta.py to obtain accession-specific, chromosome-specific, region-specific, downsampled k-mers

6.1. Via slurm array job, run hybrid_reads_acc_chr_specific_kmers_bbduk.sh to obtain reads containing acc1-specific k-mers

6.2. Via slurm array job, run hybrid_reads_acc_chr_specific_kmers_bbduk.sh to obtain reads containing acc1-specific and acc2-specific

7. Via multiple slurm array jobs (number of tasks informed by number of hybrid reads found), run hybrid_read_segments_acc_chr_specific_kmers.py to get accession-specific read segments and align them either to one genome assembly (for not_centromere) or respective genome assembly (for centromere). Note customisable variables maxAlns, minMAPQ and bt2maxLen. For this step, slurm submission scripts should be moved into main "chr_specific" directory from centromere_hybrid_read_segments_slurm_submission_scripts/ and not_centromere_hybrid_read_segments_slurm_submission_scripts/

8. Run get_chr_specific_segment_pairs_alnToSame.py for not_centromere and run get_chr_specific_segment_pairs_alnToDifferent.py for centromere (each with either "co" or "nco" specified for "--recombType" [-rt] option). These run for each chromosome and recombination type using get_chr_specific_segment_pairs_alnToSame.sh and get_chr_specific_segment_pairs_alnToDifferent.sh. As matching bases is more of a concern for centromere, it may be best to use the filtered file that's based on the nmatchTOqlen ratio rather than the alenTOqlen ratio, the latter of which is used for not_centromere due to alignment to one assembly.

9. Remove machine.files, and per-read-segment PAF and FASTA files (e.g., using centromere/find_and_remove_paf_fa.sh)

