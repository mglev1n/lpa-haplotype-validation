FROM bioconductor/tidyverse:3.17

LABEL maintainer="Michael Levin <michael.levin@pennmedicine.upenn.edu>"
LABEL description="LPA Prediction Validation Pipeline"
LABEL version="0.1.0"

# Accept GitHub PAT as build argument
ARG GITHUB_PAT
ENV GITHUB_PAT=$GITHUB_PAT

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgsl-dev \
    zip \
    unzip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install bcftools from source + shapeit5 from github
RUN git clone --recurse-submodules https://github.com/samtools/htslib.git /tmp/htslib \
    && git clone https://github.com/samtools/bcftools.git /tmp/bcftools \
    && cd /tmp/bcftools \
    && make \
    && make install \
    && cd / \
    && rm -rf /tmp/htslib /tmp/bcftools \
    # Download and install phase_common_static
    && wget -O /usr/local/bin/phase_common_static https://github.com/odelaneau/shapeit5/releases/download/v5.1.1/phase_common_static \
    && chmod +x /usr/local/bin/phase_common_static

# Install renv and required packages
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

# Copy renv lockfile to a temporary location
COPY renv.lock /tmp/renv.lock

# Pre-install all packages during build time
RUN mkdir -p /opt/lpa-pipeline && \
    cd /opt/lpa-pipeline && \
    cp /tmp/renv.lock renv.lock && \
    R -e "Sys.setenv(GITHUB_PAT=Sys.getenv('GITHUB_PAT')); \
         renv::restore(library='/opt/R-packages')"

# Set the R library path to use our pre-installed packages
ENV R_LIBS_USER=/opt/R-packages

# Copy project files to the package directory
COPY _targets.R /opt/lpa-pipeline/
COPY rmarkdown/*.Rmd /opt/lpa-pipeline/rmarkdown/
COPY rmarkdown/*.yaml /opt/lpa-pipeline/rmarkdown/
COPY Scripts/ /opt/lpa-pipeline/Scripts/
COPY Resources/ /opt/lpa-pipeline/Resources/

# Create an entrypoint script that runs in the current directory
RUN echo '#!/bin/bash \n\
set -e \n\
\n\
# Print welcome message \n\
echo "======================================================" \n\
echo "  LPA Prediction Validation Pipeline                  " \n\
echo "======================================================" \n\
echo "" \n\
\n\
# Detect if running in Singularity/Apptainer \n\
RUNNING_IN_SINGULARITY=0 \n\
if [ -e "/.singularity.d" ] || [ -e "/.apptainer.d" ] || [ -d "/singularity" ]; then \n\
  RUNNING_IN_SINGULARITY=1 \n\
  echo "Detected Singularity/Apptainer environment" \n\
fi \n\
\n\
# Check for required directories and files \n\
if [ ! -d "input" ]; then \n\
  echo "ERROR: Required directory not found: ./input" \n\
  echo "Please create the input directory in your current working directory" \n\
  exit 1 \n\
fi \n\
\n\
if [ ! -f "input/genotypes.vcf.gz" ]; then \n\
  echo "ERROR: Required input file not found: ./input/genotypes.vcf.gz" \n\
  echo "Please place your VCF file in the input directory" \n\
  exit 1 \n\
fi \n\
\n\
if [ ! -f "input/measured.csv" ]; then \n\
  echo "ERROR: Required input file not found: ./input/measured.csv" \n\
  echo "Please place your CSV file in the input directory" \n\
  exit 1 \n\
fi \n\
\n\
# Create required directories if they don'\''t exist \n\
mkdir -p Results _targets \n\
\n\
echo "Found required input files:" \n\
echo "- ./input/genotypes.vcf.gz" \n\
echo "- ./input/measured.csv" \n\
echo "" \n\
\n\
# Copy the targets script to current directory for execution \n\
cp /opt/lpa-pipeline/_targets.R ./ \n\
mkdir -p rmarkdown Scripts Resources \n\
cp /opt/lpa-pipeline/rmarkdown/* ./rmarkdown/ \n\
cp -r /opt/lpa-pipeline/Scripts/* ./Scripts/ \n\
cp -r /opt/lpa-pipeline/Resources/* ./Resources/ \n\
\n\
# Run the pipeline \n\
echo "Starting LPA prediction validation pipeline..." \n\
echo "This may take some time depending on your dataset size." \n\
echo "" \n\
\n\
  R --vanilla -e "targets::tar_make(callr_arguments = list(cmdargs = c(\"--slave\", \"--no-save\", \"--no-restore\", \"--vanilla\")))" \n\
\n\
# Check if report was generated \n\
if [ -f "Results/lpa_validation_report.html" ]; then \n\
  echo "Validation report generated successfully." \n\
  echo "You can find all results in the Results directory." \n\
else \n\
  echo "WARNING: Pipeline completed but report was not generated." \n\
  echo "Check the Results directory for any output files and error messages." \n\
fi \n\
\n\
echo "" \n\
echo "======================================================" \n\
echo "  Pipeline execution finished                         " \n\
echo "======================================================" \n\
' > /usr/local/bin/run-lpa-pipeline.sh && \
chmod +x /usr/local/bin/run-lpa-pipeline.sh

# Set the working directory to /work which will be bound to the user's current directory
WORKDIR /work

# Set the entrypoint
ENTRYPOINT ["/usr/local/bin/run-lpa-pipeline.sh"]
