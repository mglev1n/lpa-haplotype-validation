#!/bin/bash
set -e

# Enhanced entrypoint script with command line argument support

# Default values
OVERWRITE_FILES=true
DEBUG_MODE=false
VERBOSE=false
SHOW_HELP=false
SHOW_VERSION=false
RUN_PIPELINE=true

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to display usage
usage() {
    cat << EOF
LPA Prediction Validation Pipeline Container

Usage: singularity run [singularity-options] container.sif [OPTIONS]

Options:
    --no-overwrite      Don't overwrite existing pipeline files in working directory
    --debug             Enable debug mode (verbose output, keep temp files)
    --verbose           Enable verbose output
    --version           Show container version information
    --help              Show this help message
    --shell             Start an interactive shell instead of running pipeline

File Management:
    By default, the container copies pipeline files (_targets.R, Scripts/, Resources/,
    rmarkdown/) to your working directory. Use --no-overwrite to skip this if you
    have local modifications you want to preserve.

Examples:
    # Standard run (overwrites pipeline files)
    singularity run container.sif

    # Preserve local modifications to pipeline files
    singularity run container.sif --no-overwrite

    # Debug mode with verbose output
    singularity run container.sif --debug --verbose

    # Start interactive shell for debugging
    singularity run container.sif --shell

    # Check container version
    singularity run container.sif --version

Requirements:
    Your working directory must contain:
    - input/genotypes.vcf.gz
    - input/measured.csv

    Results will be written to Results/ directory.
EOF
}

# Function to show version information
show_version() {
    echo "LPA Prediction Validation Pipeline Container"
    if [[ -f "/opt/lpa-pipeline/VERSION" ]]; then
        cat /opt/lpa-pipeline/VERSION
    else
        echo "Version: unknown"
        echo "Built: unknown"
    fi
    echo ""
    echo "Container files last modified:"
    find /opt/lpa-pipeline -name "*.R" -o -name "*.sh" -o -name "*.Rmd" | head -5 | while read file; do
        echo "  $(basename $file): $(stat -c %y "$file" 2>/dev/null || stat -f %Sm "$file" 2>/dev/null || echo 'unknown')"
    done
}

# Function to log messages
log() {
    echo -e "${GREEN}[$(date +'%H:%M:%S')]${NC} $1"
}

debug_log() {
    if [[ "$DEBUG_MODE" == "true" ]] || [[ "$VERBOSE" == "true" ]]; then
        echo -e "${BLUE}[DEBUG]${NC} $1"
    fi
}

warn() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1"
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --no-overwrite)
            OVERWRITE_FILES=false
            shift
            ;;
        --debug)
            DEBUG_MODE=true
            VERBOSE=true
            shift
            ;;
        --verbose)
            VERBOSE=true
            shift
            ;;
        --help)
            SHOW_HELP=true
            shift
            ;;
        --version)
            SHOW_VERSION=true
            shift
            ;;
        --shell)
            RUN_PIPELINE=false
            shift
            ;;
        *)
            error "Unknown option: $1. Use --help for usage information."
            ;;
    esac
done

# Handle special flags
if [[ "$SHOW_HELP" == "true" ]]; then
    usage
    exit 0
fi

if [[ "$SHOW_VERSION" == "true" ]]; then
    show_version
    exit 0
fi

# Print welcome message
echo "======================================================"
echo "  LPA Prediction Validation Pipeline                  "
echo "======================================================"
echo ""

# If shell mode, start interactive shell
if [[ "$RUN_PIPELINE" == "false" ]]; then
    log "Starting interactive shell..."
    log "Pipeline files are available in /opt/lpa-pipeline/"
    log "Current working directory: $(pwd)"
    exec /bin/bash
fi

# Detect if running in Singularity/Apptainer
RUNNING_IN_SINGULARITY=0
if [ -e "/.singularity.d" ] || [ -e "/.apptainer.d" ] || [ -d "/singularity" ]; then
    RUNNING_IN_SINGULARITY=1
    debug_log "Detected Singularity/Apptainer environment"
fi

# Check for required directories and files
if [ ! -d "input" ]; then
    error "Required directory not found: ./input"
fi

if [ ! -f "input/genotypes.vcf.gz" ]; then
    error "Required input file not found: ./input/genotypes.vcf.gz"
fi

if [ ! -f "input/measured.csv" ]; then
    error "Required input file not found: ./input/measured.csv"
fi

# Create required directories if they don't exist
debug_log "Creating required directories..."
mkdir -p Results _targets

# Handle file copying based on --no-overwrite flag
copy_files() {
    local source_dir="$1"
    local dest_dir="$2"
    local file_type="$3"

    if [[ "$OVERWRITE_FILES" == "true" ]]; then
        debug_log "Copying $file_type files from $source_dir to $dest_dir"
        mkdir -p "$dest_dir"
        cp -r "$source_dir"/* "$dest_dir"/

        # Make scripts executable
        if [[ "$file_type" == "Scripts" ]]; then
            chmod +x "$dest_dir"/*.sh 2>/dev/null || true
        fi
    else
        debug_log "Skipping $file_type copy (--no-overwrite specified)"

        # Check if essential files exist
        case "$file_type" in
            "pipeline")
                if [[ ! -f "_targets.R" ]]; then
                    warn "No _targets.R found in working directory and --no-overwrite specified"
                    warn "Copying essential pipeline file..."
                    cp "$source_dir"/_targets.R ./
                fi
                ;;
            "Scripts")
                if [[ ! -d "Scripts" ]] || [[ ! -f "Scripts/lpa_processing.sh" ]]; then
                    warn "No Scripts directory found and --no-overwrite specified"
                    warn "Copying essential scripts..."
                    mkdir -p Scripts
                    cp -r "$source_dir"/* Scripts/
                    chmod +x Scripts/*.sh 2>/dev/null || true
                fi
                ;;
            "Resources")
                if [[ ! -d "Resources" ]]; then
                    warn "No Resources directory found and --no-overwrite specified"
                    warn "Copying essential resources..."
                    mkdir -p Resources
                    cp -r "$source_dir"/* Resources/
                fi
                ;;
            "rmarkdown")
                if [[ ! -d "rmarkdown" ]]; then
                    warn "No rmarkdown directory found and --no-overwrite specified"
                    warn "Copying essential rmarkdown files..."
                    mkdir -p rmarkdown
                    cp -r "$source_dir"/* rmarkdown/
                fi
                ;;
        esac
    fi
}

log "Found required input files:"
log "- ./input/genotypes.vcf.gz"
log "- ./input/measured.csv"
echo ""

# Copy pipeline files with respect to --no-overwrite flag
if [[ "$OVERWRITE_FILES" == "true" ]]; then
    log "Copying pipeline files to working directory..."
else
    log "Using existing pipeline files (--no-overwrite specified)..."
fi

copy_files "/opt/lpa-pipeline" "." "pipeline"
copy_files "/opt/lpa-pipeline/rmarkdown" "rmarkdown" "rmarkdown"
copy_files "/opt/lpa-pipeline/Scripts" "Scripts" "Scripts"
copy_files "/opt/lpa-pipeline/Resources" "Resources" "Resources"

# Show what files are being used
if [[ "$VERBOSE" == "true" ]]; then
    log "Pipeline file status:"
    echo "  _targets.R: $(if [[ -f '_targets.R' ]]; then echo 'found'; else echo 'missing'; fi)"
    echo "  Scripts/: $(if [[ -d 'Scripts' ]]; then echo "$(ls Scripts/*.sh 2>/dev/null | wc -l) scripts found"; else echo 'missing'; fi)"
    echo "  Resources/: $(if [[ -d 'Resources' ]]; then echo "$(ls Resources/*.gz Resources/*.vcf* 2>/dev/null | wc -l) files found"; else echo 'missing'; fi)"
    echo "  rmarkdown/: $(if [[ -d 'rmarkdown' ]]; then echo "$(ls rmarkdown/*.Rmd 2>/dev/null | wc -l) files found"; else echo 'missing'; fi)"
    echo ""
fi

# Run the pipeline
log "Starting LPA prediction validation pipeline..."
if [[ "$DEBUG_MODE" == "true" ]]; then
    log "Debug mode enabled - using verbose R output"
fi
echo "This may take some time depending on your dataset size."
echo ""

# Set R command arguments based on debug mode
if [[ "$DEBUG_MODE" == "true" ]]; then
    R_ARGS="--vanilla"
    R_CMD="targets::tar_make()"
else
    R_ARGS="--slave --no-save --no-restore --vanilla"
    R_CMD="targets::tar_make(callr_arguments = list(cmdargs = c('--slave', '--no-save', '--no-restore', '--vanilla')))"
fi

# Run R with appropriate arguments
if [[ "$DEBUG_MODE" == "true" ]]; then
    debug_log "Running R command: R $R_ARGS -e \"$R_CMD\""
fi

R $R_ARGS -e "$R_CMD"

# Check if report was generated
if [ -f "Results/lpa_validation_report.html" ]; then
    log "Validation report generated successfully."
    log "You can find all results in the Results directory."
else
    warn "Pipeline completed but report was not generated."
    warn "Check the Results directory for any output files and error messages."
fi

# Clean up debug information
if [[ "$DEBUG_MODE" == "true" ]]; then
    log "Debug mode: Temporary files may be preserved for inspection"
    log "Check .temp directory for intermediate files"

    # Show some debug information
    if [[ -d ".temp" ]]; then
        debug_log "Temporary directory contents:"
        find .temp -type f -name "*.txt" | head -5 | while read file; do
            debug_log "  $file"
        done
    fi
fi

# Zip results folder

if [[ -d "Results" ]]; then
    log "Zipping results directory..."
    zip -r Results.zip Results
    log "Results zipped to Results.zip"
else
    warn "Results directory not found, skipping zipping."
fi

echo ""
echo "======================================================"
echo "  Pipeline execution finished                         "
echo "======================================================"
