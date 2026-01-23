# ONT Metagenomics Pipeline

## Overview
Automated Nextflow pipeline for analyzing Oxford Nanopore Technologies (ONT) metagenomic sequencing data. The workflow performs quality control, host decontamination, taxonomic classification, per-taxa assembly/consensus generation, and creporting.

## Workflow Steps

1. **Host Depletion** (Optional): Removes host reads using Kraken2 against a specified host database (kraken2).
2. **Taxonomic Classification**: Classifies reads using Kraken2.
3. **Abundance Estimation**: Estimates species abundance using Bracken.
4. **Visualization**: Generates interactive Krona plots.
5. **Binning & Extraction**: Splits reads into separate files based on their TaxID using Bracken results.
6. ** Generation**: Assemblies are generated per TaxID using megathit.
1. **Validation**: Consensus sequences are BLASTed against a reference database.
2. **Reporting**: Generates an interactive HTML report containing:
   - Summary statistics.
   - Per-sample abundance tables.
   - Split BLAST results by organism (filtered by E-value < 1e-25, Identity > 90%, Coverage > 50%).

## Requirements

- Nextflow (>=22.12.0)
- Docker
- Reference Databases:
  - Kraken2 Host Database
  - Kraken2 Microbial Database
  - BLAST Database (e.g., custom microbial DB)

## Usage

```bash
nextflow run dhineshp565/ont_metagenomics \
    --input /path/to/fastqs/ \
    --out_dir /path/to/output/ \
    --kraken_db /path/to/kraken_db \
    --host_db /path/to/host_db \
    --blastdb_path /path/to/blast_db_dir \
    --blastdb_name "blast_db_prefix"
```

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Directory containing input FASTQ files | `filepath` |
| `--out_dir` | Output directory for results | `output` |
| `--dehost` | Enable host removal (`true` or `false`) | `true` |
| `--host_db` | Path to Kraken2 host database | `/data/referenceDB/kraken/k2_bos_taurus_db_12012026/` |
| `--kraken_db` | Path to Kraken2 microbial database | `/data/referenceDB/kraken/k2_pluspf_08gb_20230605` |
| `--blastdb_path` | Directory containing BLAST database files | `/data/referenceDB/blast/microbe_db` |
| `--blastdb_name` | Prefix name of the BLAST database | `microbe_db` |
| `--qscore` | Minimum Q-score for filtering (if applicable) | `10` |

## Output Structure

- **binned_reads/**: FASTQ files split by TaxID.
- **blast/**: Raw and processed BLAST results.
- **kraken2/**: Kraken2 classification outputs.
- **bracken/**: Bracken abundance estimation files.
- **krona_kraken/**: Interactive Krona charts.
- **megahit_metagenomes/**: Polished consensus sequences (FASTA).
- **execution/**: Nextflow execution reports (timeline, resource usage).
- **[SampleName]_report.html**: Final interactive summary report.

## Software Used

The pipeline leverages the following open-source bioinformatics tools:

| Tool | Purpose | Source |
|------|---------|--------|
| **Nextflow** | Workflow management and orchestration | [Link](https://www.nextflow.io/) |
| **Kraken2** | Taxonomic classification of reads | [Link](https://ccb.jhu.edu/software/kraken2/) |
| **Bracken** | Bayesian abundance estimation | [Link](https://ccb.jhu.edu/software/bracken/) |
| **KrakenTools** | Extracting reads by TaxID | [Link](https://github.com/jenniferlu717/KrakenTools) |
| **MEGAHIT** | Ultra-fast metagenomic assembler | [Link](https://github.com/voutcn/megahit) |
| **Krona** | Hierarchical data visualization | [Link](https://github.com/marbl/Krona) |
| **BLAST+** | Consensus verification against reference | [Link](https://blast.ncbi.nlm.nih.gov/Blast.cgi) |
| **R** | Report generation (rmarkdown, knitr, kableExtra) | [Link](https://www.r-project.org/) |


