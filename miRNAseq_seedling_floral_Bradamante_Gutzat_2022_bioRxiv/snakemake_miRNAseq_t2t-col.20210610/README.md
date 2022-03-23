# Automated workflow for small RNA-seq data processing and alignment

This is a Snakemake workflow for automated processing and alignment of small RNA-seq (sRNA-seq) data.

### Requirements

- Installation of [Snakemake](https://snakemake.readthedocs.io/en/stable/) and optionally [conda](https://conda.io/docs/)
- Demultiplexed sRNA-seq reads in gzipped FASTQ format located in the `data/` directory. These should be named according to the following naming convention: `{sample}_R1.fastq.gz`
- A samtools-indexed reference genome in FASTA format and a chromosome sizes file (e.g., `T2T_Col.fa`, `T2T_Col.fa.fai`, and `T2T_Col.fa.sizes`, the latter two of which generated with `samtools faidx T2T_Col.fa; cut -f1,2 T2T_Col.fa.fai > T2T_Col.fa.sizes`), each located in `data/index/`
- A list of potential contaminant ribosomal RNA sequences to be removed, located in `contaminants/`, described in `contaminants/ribokmers_README.txt`, which includes a link to download `ribokmers.fa.gz`
- `Snakefile` in this repository. This contains "rules" that each execute a step in the workflow
- `config.yaml` in this repository. This contains customizable parameters including `reference`, which should be the relative path to the reference genome file name without the `.fa` extension (e.g., `data/index/T2T_Col`)
- R scripts `scripts/bin_bamTPM_sRNAsizes.R` and `scripts/genomeBin_bamTPM_sRNAsizes.R` in this repository; these need to be made executable first (e.g., `chmod +x bin_bamTPM_sRNAsizes.R`)
- Optional: `environment.yaml` in this repository, used to create the software environment if conda is used
- If conda is not used, the tools listed in environment.yaml must be specified in the PATH variable

These files can be downloaded together by cloning the repository:

```
git clone https://github.com/ajtock/sRNAseq_leaf_Rigal_Mathieu_2016_PNAS.git
```

### Creating the conda environment

```
conda env create --file environment.yaml --name srna_mapping
```

### Usage

In a Unix shell, navigate to the base directory containing `Snakefile`, `config.yaml`, `environment.yaml`, and the `data/`, `scripts/` and `contaminants/` subdirectories, which should have a directory tree structure like this:

```
.
├── config.yaml
├── contaminants
│   ├── ribokmers.fa.gz
│   └── ribokmers_README.txt
├── data
│   ├── Col_0_sRNA_ERR966148_R1.fastq.gz
│   ├── index
│   │   ├── T2T_Col.1.ebwt
│   │   ├── T2T_Col.2.ebwt
│   │   ├── T2T_Col.3.ebwt
│   │   ├── T2T_Col.4.ebwt
│   │   ├── T2T_Col.fa
│   │   ├── T2T_Col.fa.fai
│   │   ├── T2T_Col.fa.sizes
│   │   ├── T2T_Col.rev.1.ebwt
│   │   └── T2T_Col.rev.2.ebwt
│   └── met1_3_sRNA_ERR966149_R1.fastq.gz
├── environment.yaml
├── LICENSE
├── README.md
├── scripts
│   ├── bin_bamTPM_sRNAsizes.R
│   └── genomeBin_bamTPM_sRNAsizes.R
└── Snakefile
```

Then run the following commands in the base directory (`--cores` should match the `THREADS` parameter in `config.yaml`):

```
conda activate srna_mapping
snakemake -p --cores 48
conda deactivate
```

### Useful Snakemake parameters

- `--cores` specifies the maximum number of threads
- `-n` performs a dry run
- `-p` prints commands
- `--use-conda`
- `--conda-prefix ~/.myconda`
- `--forcerun genomeBin_bedgraphTPM` forces rerun of a given rule (e.g., `genomeBin_bedgraphTPM`)

### Updating the conda environment

```
conda env update --file environment.yaml --name srna_mapping
```
