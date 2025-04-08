#!/bin/bash
set -e

# Print welcome message
echo "======================================================"
echo "  LPA Prediction Validation Pipeline                  "
echo "======================================================"
echo ""

# Detect if running in Singularity/Apptainer
RUNNING_IN_SINGULARITY=0
if [ -e "/.singularity.d" ] || [ -e "/.apptainer.d" ] || [ -d "/singularity" ]; then
  RUNNING_IN_SINGULARITY=1
  echo "Detected Singularity/Apptainer environment"
fi

# Setup directories for mounted data
if [ -d "/data" ]; then
  echo "Found mounted /data directory, setting up environment..."

  # Set up input files
  if [ ! -d "/data/input" ]; then
    mkdir -p /data/input
    echo "Created /data/input directory"
  fi

  # Create Results directory if it doesn't exist
  if [ ! -d "/data/Results" ]; then
    mkdir -p /data/Results
    echo "Created /data/Results directory"
  fi

  # For Singularity, we'll work directly with mounted directories
  # For Docker, we'll use symlinks
  if [ $RUNNING_IN_SINGULARITY -eq 1 ]; then
    echo "Using direct file access for Singularity"
    INPUT_DIR="/data/input"
    RESULTS_DIR="/data/Results"
  else
    echo "Creating symlinks for Docker environment"
    # Safe symlink creation - remove if exists first
    rm -f /app/input /app/Results
    ln -sf /data/input /app/input
    ln -sf /data/Results /app/Results

    INPUT_DIR="/app/input"
    RESULTS_DIR="/app/Results"
  fi
else
  echo "No data volume mounted. Using container's internal directories."
  # Ensure directories exist
  mkdir -p /app/input /app/Results

  INPUT_DIR="/app/input"
  RESULTS_DIR="/app/Results"
fi

# Check for required input files
if [ ! -f "${INPUT_DIR}/genotypes.vcf.gz" ]; then
  echo "ERROR: Required input file not found: ${INPUT_DIR}/genotypes.vcf.gz"
  echo "Please make sure to place your VCF file at input/genotypes.vcf.gz"
  exit 1
fi

if [ ! -f "${INPUT_DIR}/measured.csv" ]; then
  echo "ERROR: Required input file not found: ${INPUT_DIR}/measured.csv"
  echo "Please make sure to place your measured LPA values at input/measured.csv"
  exit 1
fi

echo "Found required input files:"
echo "- ${INPUT_DIR}/genotypes.vcf.gz"
echo "- ${INPUT_DIR}/measured.csv"
echo ""

# Run the pipeline
echo "Starting LPA prediction validation pipeline..."
echo "This may take some time depending on your dataset size."
echo ""

cd /app

# Create a temporary targets script with correct paths
if [ $RUNNING_IN_SINGULARITY -eq 1 ]; then
  # For Singularity, modify the targets paths without modifying the original file
  TEMP_TARGETS=$(mktemp)
  cat _targets.R > $TEMP_TARGETS

  # Replace input file paths
  sed -i "s|\"input/genotypes.vcf.gz\"|\"${INPUT_DIR}/genotypes.vcf.gz\"|g" $TEMP_TARGETS
  sed -i "s|\"input/measured.csv\"|\"${INPUT_DIR}/measured.csv\"|g" $TEMP_TARGETS

  # Replace output paths
  sed -i "s|\"Results/|\"${RESULTS_DIR}/|g" $TEMP_TARGETS

  # Run the pipeline with the modified file
  Rscript -e "renv::activate(); targets::tar_source('${TEMP_TARGETS}'); targets::tar_make()"

  # Clean up
  rm $TEMP_TARGETS
else
  # For Docker, we can use the original file with symlinks
  Rscript -e "renv::activate(); targets::tar_make()"
fi

echo ""
echo "Pipeline execution complete."

# Check if report was generated
if [ -f "${RESULTS_DIR}/lpa_validation_report.html" ]; then
  echo "Validation report generated successfully."
  echo "You can find all results in the Results directory."
else
  echo "WARNING: Pipeline completed but report was not generated."
  echo "Check the Results directory for any output files and error messages."
fi

# Create the ZIP archive if it doesn't exist yet
if [ ! -f "/app/lpa_validation_results.zip" ]; then
  echo "Creating results archive..."
  cd /app

  # For Singularity, zip from the mounted directory
  # For Docker, zip from the symlinked directory
  if [ $RUNNING_IN_SINGULARITY -eq 1 ]; then
    zip -r lpa_validation_results.zip ${RESULTS_DIR}/
  else
    zip -r lpa_validation_results.zip Results/
  fi

  # Copy to mounted directory if available
  if [ -d "/data" ]; then
    cp /app/lpa_validation_results.zip /data/
  fi

  echo "Results archive created: lpa_validation_results.zip"
fi

echo ""
echo "======================================================"
echo "  Pipeline execution finished                         "
echo "======================================================"
