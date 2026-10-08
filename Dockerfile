# syntax=docker/dockerfile:1
FROM bioconductor/tidyverse:3.17

LABEL maintainer="Michael Levin <michael.levin@pennmedicine.upenn.edu>"
LABEL description="LPA Prediction Validation Pipeline"

# The GitHub token needed to install private R packages (lpapredictr) is
# supplied as a BuildKit secret, not a build argument or environment variable,
# so it is never stored in an image layer or the image configuration.
# Local build:
#   docker build --secret id=github_pat,env=GITHUB_PAT -t lpa-validation .

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

# Install bcftools from source
RUN git clone --recurse-submodules https://github.com/samtools/htslib.git /tmp/htslib \
    && git clone https://github.com/samtools/bcftools.git /tmp/bcftools \
    && cd /tmp/bcftools \
    && make \
    && make install \
    && cd / \
    && rm -rf /tmp/htslib /tmp/bcftools

# Copy ShapeIt5 and impute5 static binaries from the Applications folder.
# ShapeIt5 was previously downloaded from the odelaneau/shapeit5 GitHub
# releases; GitHub disabled that repository, and the download URL then
# returned an HTML page that was installed in place of the binary.
COPY Applications/shapeit5_v5.1.1/phase_common_static /usr/local/bin/phase_common_static
COPY Applications/impute5_v1.2.0/impute5_v1.2.0_static /usr/local/bin/impute5

# Fail the build unless both files are Linux executables (ELF), so a failed
# download or a Git LFS pointer can never be installed as a binary again
RUN chmod +x /usr/local/bin/phase_common_static /usr/local/bin/impute5 \
    && for bin in /usr/local/bin/phase_common_static /usr/local/bin/impute5; do \
        if [ "$(head -c 4 "$bin" | tail -c 3)" != "ELF" ]; then \
            echo "ERROR: $bin is not an ELF executable" >&2; head -c 200 "$bin" >&2; exit 1; \
        fi; \
    done

# Install renv and required packages
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

# Copy renv lockfile to a temporary location
COPY renv.lock /tmp/renv.lock

# Pre-install all packages during build time. The token is read from the
# secret mount and exported only for this command.
RUN --mount=type=secret,id=github_pat \
    mkdir -p /opt/lpa-pipeline && \
    cd /opt/lpa-pipeline && \
    cp /tmp/renv.lock renv.lock && \
    if [ ! -s /run/secrets/github_pat ]; then \
        echo "WARNING: github_pat build secret not provided; private GitHub packages will fail to install" >&2; \
    fi && \
    GITHUB_PAT="$(cat /run/secrets/github_pat 2>/dev/null || true)" \
    R -e "renv::restore(library='/opt/R-packages')"

# Set the R library path to use our pre-installed packages
ENV R_LIBS_USER=/opt/R-packages

# Copy project files to the package directory
COPY _targets.R /opt/lpa-pipeline/
COPY rmarkdown/*.Rmd /opt/lpa-pipeline/rmarkdown/
COPY rmarkdown/*.yaml /opt/lpa-pipeline/rmarkdown/
COPY Scripts/ /opt/lpa-pipeline/Scripts/
COPY Resources/ /opt/lpa-pipeline/Resources/

# Version information, passed in by the Docker build workflow.
# Declared here (not at the top) so a new version or commit only invalidates
# the layers below, and the slow renv restore above stays cached.
ARG PIPELINE_VERSION=dev
ARG GIT_SHA=unknown
LABEL version="${PIPELINE_VERSION}"
LABEL org.opencontainers.image.revision="${GIT_SHA}"

# Create version file (shown by `--version`)
RUN echo "Version: ${PIPELINE_VERSION}" > /opt/lpa-pipeline/VERSION && \
    echo "Built: $(date -u +'%Y-%m-%d %H:%M:%S UTC')" >> /opt/lpa-pipeline/VERSION && \
    echo "Git commit: ${GIT_SHA}" >> /opt/lpa-pipeline/VERSION

# Copy the enhanced entrypoint script
COPY entrypoint.sh /usr/local/bin/run-lpa-pipeline.sh
RUN chmod +x /usr/local/bin/run-lpa-pipeline.sh

# Set the working directory to /work which will be bound to the user's current directory
WORKDIR /work

# Set the entrypoint
ENTRYPOINT ["/usr/local/bin/run-lpa-pipeline.sh"]
