# RNA-seq Analysis Pipeline - Standard Operating Procedure (SOP)

## 1. Pipeline Overview

### Purpose
This Nextflow pipeline performs RNA-seq data analysis for differential gene expression. The pipeline is optimized for speed using Kallisto pseudoalignment instead of traditional genome alignment, reducing runtime from hours to under 1 hour while maintaining accuracy for transcript quantification. It uses data from [Lupus nephritis serum induces changes in gene expression in human glomerular endothelial cells, which is modulated by L-sepiapterin: implications for redox-mediated endothelial dysfunction](https://pmc.ncbi.nlm.nih.gov/articles/PMC12161315/). This paper looked at differential expression for LN flare patients who used L-sepiapterin versus those who did not. They also had data for non-flare patients, and this pipeline focuses on comparing flare vs no-flare conditions. 

### Pipeline Workflow

![DAG Image](DAG.png)

The pipeline consists of two parallel processing chains:

#### Main Processing Chain (Sequential)
1. **Read Trimming (fastp)** - Adapter trimming and quality filtering
2. **Transcriptome Indexing (Kallisto)** - Build transcriptome index (if not provided)
3. **Create tx2gene Mapping** - Generate transcript-to-gene mapping from transcriptome FASTA
4. **Pseudoalignment (Kallisto)** - Ultra-fast transcript quantification
5. **Merge Counts** - Aggregate transcript counts to gene level using custom R script
6. **Differential Expression (DESeq2)** - Identify differentially expressed genes based on Flare status

#### Quality Control Chain (Parallel, Independent)
7. **Quality Control (FastQC)** - Quality assessment of raw reads
8. **MultiQC Report** - Comprehensive quality report aggregating all QC metrics

The QC steps run independently in parallel with the main processing chain for maximum efficiency.

## 2. Quick Start Guide
to download and prepare data, as well as the reference transcriptome and annotation files, run the `setup_and_download.sh` script included in the repository. I have it set up to downsample each FASTQ file to 20% of the original reads for faster testing. **If this is still too slow to download, go to the `input.csv` file and remove rows there, but make sure to at least have one flare and one non-flare sample.** after running the setup script, run the `run_pipeline.sh` script to execute the analysis pipeline.

## 3. Software Requirements and Versions

### Core Workflow Engine
I use these versions for reproducibility:
- **Nextflow**: version 21.04.0
- **Container Engine**: Singularity version 1.3.6-1

### Bioinformatics Tools (Containerized)

All analysis tools run in Singularity containers - no manual installation of bioinformatics software needed!

| Tool | Version | Container | Purpose |
|------|---------|-----------|---------|
| FastQC | 0.11.9 | `biocontainers/fastqc:v0.11.9_cv8` | Quality control of raw reads |
| fastp | 0.20.1 | `biocontainers/fastp:v0.20.1_cv1` | Read trimming and filtering |
| Kallisto | 0.46.2 | `zlskidmore/kallisto:0.46.2` | Pseudoalignment and quantification |
| MultiQC | 1.19 | `ewels/multiqc:v1.19` | Aggregated QC report generation |
| R/Bioconductor | 4.3 (R 4.3.1) | `rnakato/rumball` | R environment for DESeq2 analysis |

### R Packages (in rnakato/rumball container)
- **DESeq2**: Differential expression analysis
- **Base R**: Data manipulation and visualization

## 4. Environment Setup

### Prerequisites

This pipeline is designed to run on an HPC cluster with SLURM workload manager and Singularity container support.

## 5. Usage Instructions

### Complete Analysis Workflow

#### Step 1: Prepare Data and Environment

Run the automated setup script (only needed once):

```bash
./setup_and_download.sh
```

#### Step 2: Run Pipeline with Downsampled Data (Testing)

First run recommended with downsampled data for faster testing:

```bash
./run_pipeline.sh -i input_downsampled.csv
```

#### Step 3: Run Pipeline with Full Data (Production)

Once testing is successful, run with full dataset:

```bash
./run_pipeline.sh -i input_raw.csv
```

## 6. Expected Results

### Output Directory Structure

```
results/
├── reference/                 # Reference files (if auto-generated)
│   └── tx2gene.tsv
├── final_counts/              # Gene-level count matrices
│   ├── kallisto_count_matrix.txt    # Raw counts for DESeq2
│   └── kallisto_tpm_matrix.txt      # TPM normalized values
├── differential_expression/   # Differential expression results
│   ├── differential_expression_results.csv   # Full DE results
│   ├── normalized_counts.csv                 # DESeq2 normalized counts
│   └── volcano_plot.png                      # Volcano plot visualization
├── multiqc/                   # Quality control summary
│   ├── multiqc_report.html
│   └── multiqc_data/
└── pipeline_info/             # Execution reports
    ├── execution_report.html  # Resource usage and timing
    ├── timeline.html          # Process execution timeline
    ├── trace.txt              # Detailed trace of all processes
    └── dag.html               # Directed acyclic graph of workflow
```
