FROM bioconductor/tidyverse:3.17

LABEL maintainer="Michael Levin <michael.levin@pennmedicine.upenn.edu>"
LABEL description="LPA Prediction Validation Pipeline"
LABEL version="0.2.0"

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

# Create version file
RUN echo "0.2.0" > /opt/lpa-pipeline/VERSION && \
    echo "Built: $(date)" >> /opt/lpa-pipeline/VERSION && \
    echo "Git commit: ${GITHUB_SHA:-unknown}" >> /opt/lpa-pipeline/VERSION

# Copy the enhanced entrypoint script
COPY entrypoint.sh /usr/local/bin/run-lpa-pipeline.sh
RUN chmod +x /usr/local/bin/run-lpa-pipeline.sh

# Set the working directory to /work which will be bound to the user's current directory
WORKDIR /work

# Set the entrypoint
ENTRYPOINT ["/usr/local/bin/run-lpa-pipeline.sh"]
