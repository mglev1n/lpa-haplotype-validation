#!/bin/bash
set -e

# Print welcome message
echo "======================================================"
echo "  LPA Prediction Validation Pipeline                  "
echo "======================================================"
echo ""

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
echo "- /app/input/genotypes.vcf.gz"
echo "- /app/input/measured.csv"
echo ""

# Run the pipeline
echo "Starting LPA prediction validation pipeline..."
echo "This may take some time depending on your dataset size."
echo ""

cd /app
R --vanilla -s -e 'targets::tar_make()'

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

# Create the ZIP archive
echo "Creating results archive..."
cd /app
zip -r lpa_validation_results.zip Results/

echo ""
echo "======================================================"
echo "  Pipeline execution finished                         "
echo "======================================================"
