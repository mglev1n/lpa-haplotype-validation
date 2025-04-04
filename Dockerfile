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

# Install renv
RUN R -e 'install.packages("renv")'

# Copy renv lockfile
COPY renv.lock /app/

# Initialize renv and restore packages from lockfile
# Use the provided GitHub PAT to access private repositories
RUN cd /app && \
    R -e 'Sys.setenv(GITHUB_PAT=Sys.getenv("GITHUB_PAT")); renv::init(); renv::restore()'

# Create directories
RUN mkdir -p /app/input /app/Results /app/rmarkdown

# Copy project files
COPY _targets.R /app/
COPY rmarkdown/*.Rmd /app/rmarkdown/
COPY rmarkdown/*.yaml /app/rmarkdown/
COPY entrypoint.sh /app/

# Make the entrypoint script executable
RUN chmod +x /app/entrypoint.sh

# Set working directory
WORKDIR /app

# Set the entrypoint
ENTRYPOINT ["/app/entrypoint.sh"]
