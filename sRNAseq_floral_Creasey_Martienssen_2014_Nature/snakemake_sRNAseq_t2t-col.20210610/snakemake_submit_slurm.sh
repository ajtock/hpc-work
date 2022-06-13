#!/bin/bash

#! Working directory
#!SBATCH -D /home/ajt200/analysis/sRNAseq_floral_Creasey_Martienssen_2014_Nature/snakemake_sRNAseq_t2t-col.20210610

#! Account to use (users can be part of one or more accounts).
#! Currently we have only two: 'bioinf' and 'it-support'.
#SBATCH -A bioinf

#! Partition to run on
#SBATCH -p production

#! Email notification for job conditions (e.g., START, END,FAIL, HOLD or none)
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ajt200@cam.ac.uk # Doesn't seem to send these to my email inbox

#! Do not re-queue job after system fail
#SBATCH --no-requeue

#! Per-array-task-ID log file (NOTE: enclosing directory must exist)
#SBATCH -o logs/snakemake_submit.out

#! Per-array-task-ID error file (NOTE: enclosing directory must exist)
#SBATCH -e logs/snakemake_submit.err

#! Number of nodes to allocate
#SBATCH --nodes=1

#! Number of CPUs per task. Default: 1
#SBATCH --cpus-per-task=32

#! Minimum RAM needed for all tasks
#! NOTE: Doesn't work with > 1M ; e.g., with 10M, 100M, 1G, get these errors:
#! sbatch: error: Memory specification can not be satisfied
#! sbatch: error: Batch job submission failed: Requested node configuration is not available
#SBATCH --mem=120G

#! Time in HH:MM:SS. Default: 14 days (currently)
#SBATCH -t 99:00:00

#! Output some informative messages
echo "Number of CPUs used: $SLURM_CPUS_PER_TASK"
echo "This job is running on:"
hostname

#! Execute
./condor_submit.sh
