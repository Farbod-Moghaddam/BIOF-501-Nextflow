#!/bin/bash

set -e  # Exit on error

echo "================================"
echo "RNA-seq Analysis Pipeline"
echo "(Singularity-based for maximum speed)"
echo "================================"

# Default parameters
RESUME=""
INPUT_CSV="input.csv"
HELP=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -r|--resume)
            RESUME="-resume"
            shift
            ;;
        -i|--input)
            INPUT_CSV="$2"
            shift 2
            ;;
        -h|--help)
            HELP=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            HELP=true
            shift
            ;;
    esac
done

# Display help
if [ "$HELP" = true ]; then
    echo ""
    echo "Usage: ./run_pipeline.sh [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -i, --input     Input CSV file (default: input_downsampled.csv)"
    echo "  -r, --resume    Resume previous run"
    echo "  -h, --help      Display this help message"
    echo ""
    echo "Examples:"
    echo "  ./run_pipeline.sh"
    echo "  ./run_pipeline.sh -r"
    echo "  ./run_pipeline.sh -i input_downsampled.csv"
    echo "  ./run_pipeline.sh -i input_raw.csv -r"
    echo ""
    exit 0
fi

# Check if Nextflow is installed
if ! command -v nextflow &> /dev/null; then
    echo "ERROR: Nextflow is not installed or not in PATH"
    echo "Please install Nextflow: https://www.nextflow.io/docs/latest/getstarted.html"
    exit 1
fi

# Check if input CSV exists
if [ ! -f "$INPUT_CSV" ]; then
    echo "ERROR: $INPUT_CSV not found"
    echo "Please create an input CSV file with sample information"
    exit 1
fi


# Check if Singularity is installed
if ! command -v singularity &> /dev/null; then
    echo "ERROR: Singularity is not installed or not in PATH"
    echo "Please install Singularity: https://sylabs.io/guides/latest/user-guide/"
    exit 1
fi

echo ""
echo "Starting pipeline with Singularity"
echo "Input CSV: $INPUT_CSV"
echo "Resume option: ${RESUME:-disabled}"
echo ""

nextflow run main.nf \
    --input "$INPUT_CSV" \
    $RESUME \
    -with-report results/pipeline_info/execution_report.html \
    -with-timeline results/pipeline_info/timeline.html \
    -with-trace results/pipeline_info/trace.txt \
    -with-dag results/pipeline_info/dag.html

echo ""
echo "================================"
echo "Pipeline execution completed!"
echo "Check results in: results/"
echo "================================"
