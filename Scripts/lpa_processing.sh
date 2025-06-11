#!/bin/bash
set -euo pipefail

# Script: preprocess_lpa_genotypes.sh
# Purpose: Preprocess genotype data for LPA prediction
# Requirements: bcftools, shapeit5

# Function to detect available CPU cores
get_available_threads() {
    local threads=1

    # Try nproc first (most common on Linux)
    if command -v nproc &> /dev/null; then
        threads=$(nproc)
    # Try sysctl for macOS
    elif command -v sysctl &> /dev/null; then
        threads=$(sysctl -n hw.ncpu 2>/dev/null || echo 1)
    # Fallback to /proc/cpuinfo
    elif [[ -r /proc/cpuinfo ]]; then
        threads=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || echo 1)
    fi

    echo "$threads"
}

# Default values - auto-detect available threads
THREADS=$(get_available_threads)
REGION="chr6:159500000-161700000"
EXTRACTION_REGION="chr6:160400000-160800000"
TEMP_BASE_DIR=""  # New variable for custom temp directory

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to display usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Required arguments:
    -i, --input          Input VCF/BCF file (must be phased)
    -o, --output-dir     Output directory
    -x, --sites          Model sites VCF file (e.g., Lpa_hap_model.flank100.sites.vcf.gz)
    -r, --reference      Reference panel VCF for imputation (e.g., 1000G)
    -m, --map            Genetic map file for chromosome 6

Optional arguments:
    -t, --threads        Number of threads (default: auto-detect, currently $THREADS)
    -s, --shapeit5       Path to ShapeIt5 phase_common_static binary
    --region             Target region (default: chr6:159500000-161700000)
    --extract-region     Extraction region (default: chr6:160400000-160800000)
    --temp-dir           Base directory for temporary files (default: current directory)
    -h, --help           Display this help message

Example:
    $0 -i input.vcf.gz -o /path/to/output -x model.sites.vcf.gz \\
       -r /path/to/reference.vcf.gz -m /path/to/chr6.gmap.gz
EOF
    exit 1
}

# Function to log messages
log() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1" >&2
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
    exit 1
}

warn() {
    echo -e "${YELLOW}[WARNING]${NC} $1" >&2
}

# Function to check if file exists
check_file() {
    local file=$1
    local desc=$2
    if [[ ! -f "$file" ]]; then
        error "$desc file not found: $file"
    fi
}

# Function to check if command exists
check_command() {
    local cmd=$1
    if ! command -v "$cmd" &> /dev/null; then
        error "Required command not found: $cmd. Please ensure it's installed and in PATH."
    fi
}

# Function to cleanup temporary directory
cleanup_temp() {
    if [[ -n "${TEMP_DIR:-}" ]] && [[ -d "$TEMP_DIR" ]]; then
        log "Cleaning up temporary directory..."
        rm -rf "$TEMP_DIR"
    fi
}

# Parse command line arguments
PARAMS=""
while (( "$#" )); do
    case "$1" in
        -i|--input)
            INPUT_VCF=$2
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR=$2
            shift 2
            ;;
        -x|--sites)
            SITES_VCF=$2
            shift 2
            ;;
        -r|--reference)
            REFERENCE_PANEL=$2
            shift 2
            ;;
        -m|--map)
            GENETIC_MAP=$2
            shift 2
            ;;
        -t|--threads)
            THREADS=$2
            shift 2
            ;;
        -s|--shapeit5)
            SHAPEIT5_BIN=$2
            shift 2
            ;;
        --region)
            REGION=$2
            shift 2
            ;;
        --extract-region)
            EXTRACTION_REGION=$2
            shift 2
            ;;
        --temp-dir)
            TEMP_BASE_DIR=$2
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        -*|--*=)
            error "Unsupported flag $1"
            ;;
        *)
            PARAMS="$PARAMS $1"
            shift
            ;;
    esac
done

# Validate that threads is a positive integer
if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    error "Invalid number of threads: $THREADS. Must be a positive integer."
fi

# Validate required arguments
[[ -z "${INPUT_VCF:-}" ]] && error "Input VCF/BCF file is required (-i)"
[[ -z "${OUTPUT_DIR:-}" ]] && error "Output directory is required (-o)"
[[ -z "${SITES_VCF:-}" ]] && error "Model sites VCF file is required (-x)"
[[ -z "${REFERENCE_PANEL:-}" ]] && error "Reference panel is required (-r)"
[[ -z "${GENETIC_MAP:-}" ]] && error "Genetic map is required (-m)"

# Check for required commands
log "Checking required tools..."
check_command "bcftools"

# Set default ShapeIt5 path if not provided
if [[ -z "${SHAPEIT5_BIN:-}" ]]; then
    if command -v phase_common_static &> /dev/null; then
        SHAPEIT5_BIN="phase_common_static"
    else
        error "ShapeIt5 binary not found. Please specify with -s flag or ensure phase_common_static is in PATH"
    fi
fi

# Validate input files
log "Validating input files..."
check_file "$INPUT_VCF" "Input VCF/BCF"
check_file "$SITES_VCF" "Model sites VCF"
check_file "$REFERENCE_PANEL" "Reference panel"
check_file "$GENETIC_MAP" "Genetic map"
check_file "$SHAPEIT5_BIN" "ShapeIt5 binary"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Create temporary directory for intermediate files
# Use current directory if no temp base directory specified
if [[ -z "$TEMP_BASE_DIR" ]]; then
    TEMP_BASE_DIR="$(pwd)"
fi

# Ensure temp base directory exists and is writable
if [[ ! -d "$TEMP_BASE_DIR" ]]; then
    error "Temporary directory base does not exist: $TEMP_BASE_DIR"
fi

if [[ ! -w "$TEMP_BASE_DIR" ]]; then
    error "Temporary directory base is not writable: $TEMP_BASE_DIR"
fi

# Create unique temporary directory in the specified base directory
TEMP_DIR=$(mktemp -d "$TEMP_BASE_DIR/lpa_preprocess.XXXXXX")
if [[ ! -d "$TEMP_DIR" ]]; then
    error "Failed to create temporary directory in $TEMP_BASE_DIR"
fi

trap cleanup_temp EXIT

log "Using temporary directory: $TEMP_DIR"
log "Using $THREADS threads for parallel processing"

# Set up file paths
BASENAME=$(basename "$INPUT_VCF" | sed 's/\.[^.]*$//')
CHR_FIXED_VCF="$TEMP_DIR/${BASENAME}.chr6.bcf"
EXTRACTED_BCF="$TEMP_DIR/${BASENAME}.extracted.bcf"
IMPUTED_BCF="$TEMP_DIR/${BASENAME}.imputed.bcf"
FINAL_BCF="$OUTPUT_DIR/${BASENAME}.processed.bcf"
CHR_RENAME_FILE="$TEMP_DIR/chr_rename.txt"
LOG_FILE="$OUTPUT_DIR/preprocessing.log"

# Redirect stdout and stderr to log file while still displaying on console
exec > >(tee -a "$LOG_FILE")
exec 2>&1

log "Starting LPA genotype preprocessing pipeline"
log "Input: $INPUT_VCF"
log "Output directory: $OUTPUT_DIR"
log "Model sites: $SITES_VCF"

# Step 1: Index input VCF
log "Creating index for input VCF/BCF file..."
bcftools index -f "$INPUT_VCF"

# Step 2: Determine chromosome naming convention
log "Determining chromosome naming convention..."
CHR=$(bcftools index -s "$INPUT_VCF" | awk '$1==6 || $1=="chr6" {print $1; exit}')

if [[ -z "$CHR" ]]; then
    error "VCF/BCF file does not contain chromosome '6' or 'chr6'"
fi

log "Detected chromosome naming: $CHR"

# Step 3: Standardize chromosome naming to chr6 if needed
if [[ "$CHR" == "6" ]]; then
    log "Converting chromosome naming from 6 to chr6..."
    echo -e "6\tchr6" > "$CHR_RENAME_FILE"
    bcftools annotate --rename-chrs "$CHR_RENAME_FILE" \
        "$INPUT_VCF" -Ob -o "$CHR_FIXED_VCF"
    bcftools index -f "$CHR_FIXED_VCF"
    WORKING_VCF="$CHR_FIXED_VCF"
else
    log "Chromosome naming already uses chr6 format"
    WORKING_VCF="$INPUT_VCF"
fi

# Step 4: Extract variants at model sites
log "Extracting variants at model sites..."
bcftools isec -c none -r "$EXTRACTION_REGION" -n=2 -w1 -Ou \
    "$WORKING_VCF" \
    "$SITES_VCF" \
    | bcftools annotate -x INFO,^FMT/GT -Ob -o "$EXTRACTED_BCF"

bcftools index -f "$EXTRACTED_BCF"

# Step 5: Check coverage and missing data
log "Checking coverage of model sites and missing genotypes..."
N_MODEL_SITES=$(bcftools view -H "$SITES_VCF" | wc -l)
N_EXTRACTED=$(bcftools view -H "$EXTRACTED_BCF" | wc -l)
COVERAGE_PCT=$(awk "BEGIN {printf \"%.1f\", ($N_EXTRACTED/$N_MODEL_SITES)*100}")

log "Model sites: $N_MODEL_SITES"
log "Extracted sites: $N_EXTRACTED"
log "Coverage: $COVERAGE_PCT%"

# Check for missing genotypes per sample
bcftools stats -s- "$EXTRACTED_BCF" > "$TEMP_DIR/extracted_stats.txt"
N_SAMPLES_WITH_MISSING=$(grep "^PSC" "$TEMP_DIR/extracted_stats.txt" | awk '$14>0' | wc -l)
TOTAL_SAMPLES=$(bcftools query -l "$EXTRACTED_BCF" | wc -l)

if [[ "$N_SAMPLES_WITH_MISSING" -gt 0 ]]; then
    warn "Found $N_SAMPLES_WITH_MISSING out of $TOTAL_SAMPLES samples with missing genotypes"
    # Extract detailed missing data stats (sample name and missing count)
    grep "^PSC" "$TEMP_DIR/extracted_stats.txt" | awk '$14>0 {print $3"\t"$14}' > "$OUTPUT_DIR/samples_with_missing.txt"
fi

# Step 6: Determine if imputation is needed
NEED_IMPUTATION=false

if [[ "$N_EXTRACTED" -lt "$N_MODEL_SITES" ]]; then
    warn "Missing $(($N_MODEL_SITES - $N_EXTRACTED)) model sites. Imputation required."
    NEED_IMPUTATION=true
elif [[ "$N_SAMPLES_WITH_MISSING" -gt 0 ]]; then
    warn "All sites present but $N_SAMPLES_WITH_MISSING samples have missing genotypes. Imputation required."
    NEED_IMPUTATION=true
else
    log "All model sites present and no missing genotypes. No imputation needed."
fi

# Step 7: Perform imputation if needed
if [[ "$NEED_IMPUTATION" == "true" ]]; then
    # Check if reference panel index exists
    if [[ ! -f "${REFERENCE_PANEL}.tbi" && ! -f "${REFERENCE_PANEL}.csi" ]]; then
        log "Creating reference panel index..."
        bcftools index -f "$REFERENCE_PANEL"
    fi

    log "Running ShapeIt5 imputation with $THREADS threads..."
    log "Region: $REGION"

    "$SHAPEIT5_BIN" \
        --input "$WORKING_VCF" \
        --region "$REGION" \
        --map "$GENETIC_MAP" \
        --reference "$REFERENCE_PANEL" \
        --output "$IMPUTED_BCF" \
        --thread "$THREADS" 2>&1 | tee -a "$LOG_FILE"

    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        error "ShapeIt5 imputation failed"
    fi

    # Re-extract after imputation
    log "Re-extracting variants after imputation..."
    bcftools isec -c none -r "$EXTRACTION_REGION" -n=2 -w1 -Ou \
        "$IMPUTED_BCF" \
        "$SITES_VCF" \
        | bcftools annotate -x INFO,^FMT/GT -Ob -o "$FINAL_BCF"

    bcftools index -f "$FINAL_BCF"

    # Verify results after imputation
    N_FINAL=$(bcftools view -H "$FINAL_BCF" | wc -l)
    FINAL_COVERAGE_PCT=$(awk "BEGIN {printf \"%.1f\", ($N_FINAL/$N_MODEL_SITES)*100}")

    log "Final sites: $N_FINAL"
    log "Final coverage: $FINAL_COVERAGE_PCT%"

    if [[ "$N_FINAL" -lt "$N_MODEL_SITES" ]]; then
        warn "Still missing $(($N_MODEL_SITES - $N_FINAL)) sites after imputation"
    fi

    # Check missing genotypes after imputation
    bcftools stats -s- "$FINAL_BCF" > "$OUTPUT_DIR/final_stats.txt"
    N_SAMPLES_WITH_MISSING_FINAL=$(grep "^PSC" "$OUTPUT_DIR/final_stats.txt" | awk '$14>0' | wc -l)

    if [[ "$N_SAMPLES_WITH_MISSING_FINAL" -gt 0 ]]; then
        warn "After imputation, $N_SAMPLES_WITH_MISSING_FINAL samples still have missing genotypes"
        grep "^PSC" "$OUTPUT_DIR/final_stats.txt" | awk '$16>0 {print $3"\t"$14}' > "$OUTPUT_DIR/samples_with_missing_after_imputation.txt"
    else
        log "All missing genotypes successfully imputed"
    fi
else
    log "No imputation performed. Using extracted data as final."
    cp "$EXTRACTED_BCF" "$FINAL_BCF"
    cp "${EXTRACTED_BCF}.csi" "${FINAL_BCF}.csi"
    cp "$TEMP_DIR/extracted_stats.txt" "$OUTPUT_DIR/final_stats.txt"
fi

# Step 8: Final validation
log "Performing final validation..."

# Check for multiallelic sites
N_MULTIALLELIC=$(bcftools view -H "$FINAL_BCF" | awk '$5 ~ /,/' | wc -l)
if [[ "$N_MULTIALLELIC" -gt 0 ]]; then
    warn "Found $N_MULTIALLELIC multiallelic sites. Consider splitting with bcftools norm."
fi

# Get final sample count with missing data
FINAL_MISSING=$(grep "^PSC" "$OUTPUT_DIR/final_stats.txt" | awk '$14>0' | wc -l 2>/dev/null || echo "0")

# Create summary report
log "Creating summary report..."
cat > "$OUTPUT_DIR/preprocessing_summary.txt" << EOF
LPA Genotype Preprocessing Summary
==================================
Date: $(date)
Input: $INPUT_VCF
Output: $FINAL_BCF
Threads used: $THREADS
Temporary directory: $TEMP_DIR

Original chromosome naming: $CHR
Standardized to: chr6
Model sites required: $N_MODEL_SITES
Initial extracted sites: $N_EXTRACTED (${COVERAGE_PCT}%)
Samples with missing data (before imputation): $N_SAMPLES_WITH_MISSING / $TOTAL_SAMPLES
Imputation performed: $NEED_IMPUTATION
Final sites: $(bcftools view -H "$FINAL_BCF" | wc -l)
Total samples: $(bcftools query -l "$FINAL_BCF" | wc -l)
Multiallelic sites: $N_MULTIALLELIC
Samples with missing data (final): $FINAL_MISSING
EOF

cat "$OUTPUT_DIR/preprocessing_summary.txt"

log "Preprocessing completed successfully!"
log "Ready-to-use BCF file: $FINAL_BCF"

echo $FINAL_BCF
