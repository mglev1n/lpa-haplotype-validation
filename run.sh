#!/bin/bash
set -e

# Print welcome message
echo "======================================================"
echo "  LPA Prediction Validation Pipeline                  "
echo "======================================================"
echo ""

# Check for renv and activate it
echo "Checking renv setup..."
if [ ! -f "renv.lock" ]; then
  echo "ERROR: renv.lock file not found"
  echo "Please make sure you're in the correct directory with the renv.lock file"
  exit 1
fi

# Ensure renv is installed and activate it
Rscript -e '
if(!requireNamespace("renv", quietly = TRUE)) {
  cat("Installing renv package...\n")
  install.packages("renv")
}
if(!file.exists("renv/activate.R")) {
  cat("Initializing renv...\n")
  renv::init()
}
cat("Restoring packages from lockfile...\n")
renv::restore()
cat("renv setup complete.\n")
'
if [ $? -ne 0 ]; then
  echo "ERROR: Failed to set up renv environment. Please check the error messages above."
  exit 1
fi

# Create required directories
mkdir -p input Results rmarkdown

# Check for required input files
if [ ! -f "input/genotypes.vcf.gz" ]; then
  echo "ERROR: Required input file not found: input/genotypes.vcf.gz"
  echo "Please make sure to place your VCF file at input/genotypes.vcf.gz"
  exit 1
fi

if [ ! -f "input/measured.csv" ]; then
  echo "ERROR: Required input file not found: input/measured.csv"
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

Rscript -e 'renv::activate(); targets::tar_make()'

echo ""
echo "Pipeline execution complete."

# Check if report was generated
if [ -f "Results/lpa_validation_report.html" ]; then
  echo "Validation report generated successfully."
  echo "You can find all results in the Results directory."
else
  echo "WARNING: Pipeline completed but report was not generated."
  echo "Check the Results directory for any output files and error messages."
fi

# Check if ZIP archive was created
if [ -f "lpa_validation_results.zip" ]; then
  echo "Results archive created: lpa_validation_results.zip"
else
  echo "Creating results archive..."
  zip -r lpa_validation_results.zip Results/
  echo "Results archive created: lpa_validation_results.zip"
fi

echo ""
echo "======================================================"
echo "  Pipeline execution finished                         "
echo "======================================================"
