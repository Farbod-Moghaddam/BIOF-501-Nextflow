#!/bin/bash

set -e  # Exit on error

SUBSAMPLE_FRACTION=0.2
SEED=42
PARALLEL_JOBS=9  # Number of parallel downloads

echo "=========================================="
echo "Setup and Data Download Script"
echo "=========================================="

# ============================================================================
# 1. Install required tools
# ============================================================================

echo ""
echo "Step 1: Installing required tools..."
echo "--------------------------------------"

# Install SRA Toolkit
if [ ! -d "$HOME/sratoolkit.3.0.0-ubuntu64" ]; then
    echo "Installing SRA Toolkit..."
    cd ~
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-ubuntu64.tar.gz
    tar -xzf sratoolkit.3.0.0-ubuntu64.tar.gz
    rm sratoolkit.3.0.0-ubuntu64.tar.gz
    echo "SRA Toolkit installed"
else
    echo "SRA Toolkit already installed"
fi

# Install seqtk
if [ ! -d "$HOME/seqtk" ]; then
    echo "Installing seqtk..."
    cd ~
    git clone https://github.com/lh3/seqtk.git
    cd seqtk
    make
    echo "seqtk installed"
else
    echo "seqtk already installed"
fi

export PATH=$PATH:~/sratoolkit.3.0.0-ubuntu64/bin:~/seqtk

# ============================================================================
# 2. Download reference files
# ============================================================================

echo ""
echo "Step 2: Downloading reference files..."
echo "--------------------------------------"

cd ${SLURM_SUBMIT_DIR}
mkdir -p reference

# Download transcriptome
if [ ! -f "reference/Homo_sapiens.GRCh38.cdna.all.fa.gz" ]; then
    echo "Downloading human transcriptome (cDNA)..."
    wget -P reference ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
    echo "Transcriptome downloaded"
else
    echo "Transcriptome already exists"
fi

# Download GTF annotation
if [ ! -f "reference/Homo_sapiens.GRCh38.110.gtf" ]; then
    echo "Downloading GTF annotation..."
    wget -P reference ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
    gunzip reference/Homo_sapiens.GRCh38.110.gtf.gz
    echo "GTF downloaded"
else
    echo "GTF already exists"
fi

# ============================================================================
# 3. Download and process SRA data
# ============================================================================

echo ""
echo "Step 3: Downloading and processing SRA files..."
echo "--------------------------------------"

# Create directories
mkdir -p input/raw
mkdir -p input/downsampled

# Create temporary directory
TEMP_DIR="./tmp_sra_download"
mkdir -p "$TEMP_DIR"

# Extract unique SRR IDs from input.csv
tail -n +2 input.csv | cut -d',' -f1 | sort -u > /tmp/srr_list.txt

# Function to download, convert, and downsample a single SRR
process_srr() {
    local srr=$1
    
    echo "[$srr] Prefetching..."
    prefetch "$srr" --output-directory "$TEMP_DIR"
    
    echo "[$srr] Converting to FASTQ..."
    fasterq-dump --split-files --threads 3 --outdir "$TEMP_DIR" --temp "$TEMP_DIR" "${TEMP_DIR}/${srr}/${srr}.sra"
    
    echo "[$srr] Compressing raw files..."
    pigz -p 3 "${TEMP_DIR}/${srr}".sra_*.fastq
    
    echo "[$srr] Copying raw files..."
    cp "${TEMP_DIR}/${srr}".sra_*.fastq.gz "./input/raw/"
    
    echo "[$srr] Downsampling to ${SUBSAMPLE_FRACTION}..."
    for fastq_gz in "${TEMP_DIR}/${srr}".sra_*.fastq.gz; do
        basename_file=$(basename "$fastq_gz")
        seqtk sample -s${SEED} "$fastq_gz" ${SUBSAMPLE_FRACTION} | gzip > "./input/downsampled/${basename_file}"
    done
    
    echo "[$srr] Cleaning up temp space..."
    rm -rf "${TEMP_DIR}/${srr}"
    rm -f "${TEMP_DIR}/${srr}"_*.fastq.gz
    
    echo "[$srr] Complete!"
}

export -f process_srr
export TEMP_DIR SUBSAMPLE_FRACTION SEED

# Process all SRR files in parallel
cat /tmp/srr_list.txt | xargs -P ${PARALLEL_JOBS} -I {} bash -c 'process_srr "$@"' _ {}

# Clean up
rm -rf "$TEMP_DIR"
rm /tmp/srr_list.txt

# ============================================================================
# 4. Update input CSV files
# ============================================================================

echo ""
echo "Step 4: Creating input CSV files..."
echo "--------------------------------------"

# Create input_raw.csv for full data
echo "Sample,R1,R2,Flare" > input_raw.csv
tail -n +2 input.csv | while IFS=',' read -r sample r1 r2 flare; do
    # Extract SRR ID from sample name
    srr=$(echo "$sample" | grep -oE 'SRR[0-9]+')
    echo "${sample},input/raw/${srr}_1.fastq.gz,input/raw/${srr}_2.fastq.gz,${flare}"
done >> input_raw.csv

# Create input_downsampled.csv for downsampled data
echo "Sample,R1,R2,Flare" > input_downsampled.csv
tail -n +2 input.csv | while IFS=',' read -r sample r1 r2 flare; do
    srr=$(echo "$sample" | grep -oE 'SRR[0-9]+')
    echo "${sample},input/downsampled/${srr}_1.fastq.gz,input/downsampled/${srr}_2.fastq.gz,${flare}"
done >> input_downsampled.csv

echo ""
echo "=========================================="
echo "Setup Complete!"
echo "=========================================="
echo ""
echo "Tools installed:"
echo "  - SRA Toolkit: ~/sratoolkit.3.0.0-ubuntu64"
echo "  - seqtk: ~/seqtk"
echo ""
echo "Reference files downloaded:"
echo "  - Transcriptome: reference/Homo_sapiens.GRCh38.cdna.all.fa.gz"
echo "  - GTF: reference/Homo_sapiens.GRCh38.110.gtf"
echo ""
echo "Data downloaded and processed:"
echo "  - Raw FASTQ files: input/raw/"
echo "  - Downsampled FASTQ files (${SUBSAMPLE_FRACTION}): input/downsampled/"
echo ""
echo "CSV files created:"
echo "  - input_raw.csv - for full dataset analysis"
echo "  - input_downsampled.csv - for quick testing"
echo ""
echo "To run the pipeline with full data:"
echo "  ./run_pipeline.sh --input input_raw.csv"
echo ""
echo "To run the pipeline with downsampled data:"
echo "  ./run_pipeline.sh --input input_downsampled.csv"
echo ""
echo "=========================================="
