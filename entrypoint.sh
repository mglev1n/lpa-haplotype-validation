#!/bin/bash
set -e

# Print welcome message
echo "======================================================"
echo "  LPA Prediction Validation Pipeline                  "
echo "======================================================"
echo ""

# Setup directories for mounted data
if [ -d "/data" ]; then
  echo "Found mounted /data directory, setting up environment..."

  # Set up input files
  if [ ! -d "/data/input" ]; then
    mkdir -p /data/input
    echo "Created /data/input directory"
  fi

  # Create symbolic links to input and Results directories
  ln -sf /data/input /app/input

  # Check if Results directory exists in mounted volume
  if [ ! -d "/data/Results" ]; then
    mkdir -p /data/Results
    echo "Created /data/Results directory"
  fi

  # Link Results directory
  ln -sf /data/Results /app/Results

  echo "Environment setup complete."
else
  echo "No data volume mounted. Using container's internal directories."
  # Ensure directories exist
  mkdir -p /app/input /app/Results
fi

# Check for required input files
if [ ! -f "/app/input/genotypes.vcf.gz" ]; then
  echo "ERROR: Required input file not found: /app/input/genotypes.vcf.gz"
  echo "Please make sure to place your VCF file at input/genotypes.vcf.gz"
  exit 1
fi

if [ ! -f "/app/input/measured.csv" ]; then
  echo "ERROR: Required input file not found: /app/input/measured.csv"
  echo "Please make sure to place your measured LPA values at input/measured.csv"
  exit 1
fi

echo "Found required input files:"
echo "- input/genotypes.vcf.gz"
echo "- input/measured.csv"
echo ""

# Run the pipeline
echo "Starting LPA prediction validation pipeline..."
echo "This may take some time depending on your dataset size."
echo ""

cd /app
Rscript -e 'renv::activate(); targets::tar_make()'

echo ""
echo "Pipeline execution complete."

# Check if report was generated
if [ -f "/app/Results/lpa_validation_report.html" ]; then
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
  zip -r lpa_validation_results.zip Results/

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
