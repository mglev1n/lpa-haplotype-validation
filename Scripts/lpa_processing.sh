#!/bin/bash
set -euo pipefail

# Script: preprocess_lpa_genotypes.sh
# Purpose: Preprocess genotype data for LPA prediction using ShapeIt5 for phasing and impute5 for formal imputation
# Requirements: bcftools, shapeit5, impute5

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
OUTPUT_FILENAME=""  # New variable for custom output filename
IMPUTE5_BIN=""  # Path to impute5 binary

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
    -f, --output-file    Output filename (default: {input_basename}.processed.bcf)
    -t, --threads        Number of threads (default: auto-detect, currently $THREADS)
    -s, --shapeit5       Path to ShapeIt5 phase_common_static binary
    -p, --impute5        Path to impute5 binary (for formal imputation if needed)
    --region             Target region (default: chr6:159500000-161700000)
    --extract-region     Extraction region (default: chr6:160400000-160800000)
    --temp-dir           Base directory for temporary files (default: current directory)
                         Note: Temporary files will be created in a .temp subdirectory
    -h, --help           Display this help message

Example:
    $0 -i input.vcf.gz -o /path/to/output -x model.sites.vcf.gz \\
       -r /path/to/reference.vcf.gz -m /path/to/chr6.gmap.gz

    $0 -i input.vcf.gz -o /path/to/output -f custom_output.bcf \\
       -x model.sites.vcf.gz -r /path/to/reference.vcf.gz -m /path/to/chr6.gmap.gz \\
       -p /usr/local/bin/impute5 -s /usr/local/bin/phase_common_static
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

    # Also clean up the .temp directory if it's empty
    if [[ -n "${TEMP_BASE_DIR:-}" ]] && [[ -d "$TEMP_BASE_DIR/.temp" ]]; then
        # Only remove if empty (will fail silently if not empty, which is fine)
        rmdir "$TEMP_BASE_DIR/.temp" 2>/dev/null || true
    fi
}

# Function to fill AN/AC tags
fill_tags() {
    local input_file=$1
    local output_file=$2
    local description=$3

    log "Filling AN/AC tags for $description..."
    bcftools +fill-tags "$input_file" -Ob -o "$output_file" -- -t AN,AC
    bcftools index -f "$output_file"
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
        -f|--output-file)
            OUTPUT_FILENAME=$2
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
        -p|--impute5)
            IMPUTE5_BIN=$2
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

# Validate output filename if provided
if [[ -n "$OUTPUT_FILENAME" ]]; then
    # Check if filename has appropriate extension
    if [[ ! "$OUTPUT_FILENAME" =~ \.(bcf|vcf|vcf\.gz)$ ]]; then
        warn "Output filename '$OUTPUT_FILENAME' does not have a standard extension (.bcf, .vcf, or .vcf.gz)"
        warn "Proceeding anyway, but you may want to use a standard extension"
    fi

    # Check if file already exists
    if [[ -f "$OUTPUT_DIR/$OUTPUT_FILENAME" ]]; then
        warn "Output file already exists: $OUTPUT_DIR/$OUTPUT_FILENAME"
        warn "It will be overwritten"
    fi
fi

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

# Set default impute5 path if not provided
if [[ -z "${IMPUTE5_BIN:-}" ]]; then
    if command -v impute5 &> /dev/null; then
        IMPUTE5_BIN="impute5"
        log "Found impute5 in PATH: $(which impute5)"
    else
        warn "impute5 binary not found in PATH. Formal imputation will not be available if ShapeIt5 cannot fill all missing sites."
        warn "To enable full imputation, specify path with -p flag or ensure impute5 is in PATH"
    fi
else
    # Validate the provided impute5 path
    if [[ ! -f "$IMPUTE5_BIN" ]]; then
        error "Specified impute5 binary not found: $IMPUTE5_BIN"
    fi
    if [[ ! -x "$IMPUTE5_BIN" ]]; then
        error "Specified impute5 binary is not executable: $IMPUTE5_BIN"
    fi
    log "Using specified impute5 binary: $IMPUTE5_BIN"
fi

# Validate input files
log "Validating input files..."
check_file "$INPUT_VCF" "Input VCF/BCF"
check_file "$SITES_VCF" "Model sites VCF"
check_file "$REFERENCE_PANEL" "Reference panel"
check_file "$GENETIC_MAP" "Genetic map"
check_file "$SHAPEIT5_BIN" "ShapeIt5 binary"

# Validate impute5 binary if provided
if [[ -n "${IMPUTE5_BIN:-}" ]]; then
    check_file "$IMPUTE5_BIN" "impute5 binary"
fi

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

# Create .temp subdirectory within the base directory
TEMP_PARENT="$TEMP_BASE_DIR/.temp"
mkdir -p "$TEMP_PARENT"

if [[ ! -d "$TEMP_PARENT" ]]; then
    error "Failed to create .temp directory in $TEMP_BASE_DIR"
fi

# Create unique temporary directory within the .temp directory
TEMP_DIR=$(mktemp -d "$TEMP_PARENT/lpa_preprocess.XXXXXX")
if [[ ! -d "$TEMP_DIR" ]]; then
    error "Failed to create temporary directory in $TEMP_PARENT"
fi

# Set up comprehensive cleanup traps for various signals
# This helps ensure cleanup even if the process is killed
trap cleanup_temp EXIT
trap cleanup_temp SIGTERM
trap cleanup_temp SIGINT
trap cleanup_temp SIGQUIT

log "Using temporary directory: $TEMP_DIR"
log "Using $THREADS threads for parallel processing"

# Note: Create a file to mark this run in case cleanup is needed later
echo "LPA processing started at $(date)" > "$TEMP_DIR/process_info.txt"
echo "PID: $$" >> "$TEMP_DIR/process_info.txt"

# Set up file paths
BASENAME=$(basename "$INPUT_VCF" | sed 's/\.[^.]*$//')
CHR_FIXED_VCF="$TEMP_DIR/${BASENAME}.chr6.bcf"
CHR_FIXED_TAGGED_VCF="$TEMP_DIR/${BASENAME}.chr6.tagged.bcf"
EXTRACTED_BCF="$TEMP_DIR/${BASENAME}.extracted.bcf"
IMPUTED_BCF="$TEMP_DIR/${BASENAME}.imputed.bcf"
IMPUTED_TAGGED_BCF="$TEMP_DIR/${BASENAME}.imputed.tagged.bcf"
TEMP_FINAL_BCF="$TEMP_DIR/${BASENAME}.final.bcf"

# Set output filename - use custom name if provided, otherwise use default
if [[ -n "$OUTPUT_FILENAME" ]]; then
    FINAL_BCF="$OUTPUT_DIR/$OUTPUT_FILENAME"
else
    FINAL_BCF="$OUTPUT_DIR/${BASENAME}.processed.bcf"
fi

CHR_RENAME_FILE="$TEMP_DIR/chr_rename.txt"
LOG_FILE="$OUTPUT_DIR/preprocessing.log"

# Redirect stdout and stderr to log file while still displaying on console
exec > >(tee -a "$LOG_FILE")
exec 2>&1

log "Starting LPA genotype preprocessing pipeline"
log "Input: $INPUT_VCF"
log "Output directory: $OUTPUT_DIR"
log "Output file: $FINAL_BCF"
log "Model sites: $SITES_VCF"

# Step 1: Check if input VCF is indexed
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

    # Fill AN/AC tags after chromosome renaming
    fill_tags "$CHR_FIXED_VCF" "$CHR_FIXED_TAGGED_VCF" "chromosome-renamed data"
    WORKING_VCF="$CHR_FIXED_TAGGED_VCF"
else
    log "Chromosome naming already uses chr6 format"
    # Still need to fill AN/AC tags for the original file
    fill_tags "$INPUT_VCF" "$CHR_FIXED_TAGGED_VCF" "input data"
    WORKING_VCF="$CHR_FIXED_TAGGED_VCF"
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

    # Fill AN/AC tags after imputation
    fill_tags "$IMPUTED_BCF" "$IMPUTED_TAGGED_BCF" "imputed data"

    # Re-extract after imputation
    log "Re-extracting variants after imputation..."
    bcftools isec -c none -r "$EXTRACTION_REGION" -n=2 -w1 -Ou \
        "$IMPUTED_TAGGED_BCF" \
        "$SITES_VCF" \
        | bcftools annotate -x INFO,^FMT/GT -Ob -o "$TEMP_FINAL_BCF"

    # Check coverage after shapeit5
    N_SHAPEIT5_FINAL=$(bcftools view -H "$TEMP_FINAL_BCF" | wc -l)
    SHAPEIT5_COVERAGE_PCT=$(awk "BEGIN {printf \"%.1f\", ($N_SHAPEIT5_FINAL/$N_MODEL_SITES)*100}")

    log "Sites after ShapeIt5: $N_SHAPEIT5_FINAL"
    log "Coverage after ShapeIt5: $SHAPEIT5_COVERAGE_PCT%"

    # Step 7b: If sites are still missing after shapeit5, use impute5 for formal imputation
    if [[ "$N_SHAPEIT5_FINAL" -lt "$N_MODEL_SITES" ]]; then
        warn "Still missing $(($N_MODEL_SITES - $N_SHAPEIT5_FINAL)) sites after ShapeIt5 phasing."
        log "Running formal imputation with impute5..."

        # Check if impute5 binary is available
        if [[ -z "${IMPUTE5_BIN:-}" ]]; then
            error "impute5 binary not found. Please specify with -p flag or ensure impute5 is in PATH"
        fi

        # Define impute5 output files
        IMPUTE5_BCF="$TEMP_DIR/${BASENAME}.impute5.bcf"
        IMPUTE5_TAGGED_BCF="$TEMP_DIR/${BASENAME}.impute5.tagged.bcf"

        # Run impute5 with the phased data from shapeit5
        log "Running impute5 with $THREADS threads..."
        "$IMPUTE5_BIN" \
            --h "$REFERENCE_PANEL" \
            --m "$GENETIC_MAP" \
            --g "$IMPUTED_TAGGED_BCF" \
            --r "$REGION" \
            --o "$IMPUTE5_BCF" \
            --threads "$THREADS" 2>&1 | tee -a "$LOG_FILE"

        if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
            error "impute5 imputation failed"
        fi

        # Fill AN/AC tags after impute5
        fill_tags "$IMPUTE5_BCF" "$IMPUTE5_TAGGED_BCF" "impute5 imputed data"

        # Re-extract after impute5
        log "Re-extracting variants after impute5..."
        bcftools isec -c none -r "$EXTRACTION_REGION" -n=2 -w1 -Ou \
            "$IMPUTE5_TAGGED_BCF" \
            "$SITES_VCF" \
            | bcftools annotate -x INFO,^FMT/GT -Ob -o "$TEMP_FINAL_BCF"

        # Update the imputed BCF reference for downstream steps
        IMPUTED_TAGGED_BCF="$IMPUTE5_TAGGED_BCF"

        log "impute5 formal imputation completed"
    fi

    # Fill AN/AC tags for final extracted data
    fill_tags "$TEMP_FINAL_BCF" "$FINAL_BCF" "final extracted data"

    # Verify results after all imputation steps
    N_FINAL=$(bcftools view -H "$FINAL_BCF" | wc -l)
    FINAL_COVERAGE_PCT=$(awk "BEGIN {printf \"%.1f\", ($N_FINAL/$N_MODEL_SITES)*100}")

    log "Final sites: $N_FINAL"
    log "Final coverage: $FINAL_COVERAGE_PCT%"

    if [[ "$N_FINAL" -lt "$N_MODEL_SITES" ]]; then
        warn "Still missing $(($N_MODEL_SITES - $N_FINAL)) sites after all imputation attempts"
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
    log "No imputation performed. Filling AN/AC tags and using extracted data as final."
    # Fill AN/AC tags for the extracted data before making it final
    fill_tags "$EXTRACTED_BCF" "$FINAL_BCF" "final extracted data"
    cp "$TEMP_DIR/extracted_stats.txt" "$OUTPUT_DIR/final_stats.txt"
    # Set N_SHAPEIT5_FINAL for summary report
    N_SHAPEIT5_FINAL="N/A (no imputation needed)"
fi

# Step 8: Final validation
log "Performing final validation..."

# Verify AN/AC tags are present
log "Verifying AN/AC tags are present in final output..."
AN_COUNT=$(bcftools view -h "$FINAL_BCF" | grep -c "##INFO=.*ID=AN" || echo "0")
AC_COUNT=$(bcftools view -h "$FINAL_BCF" | grep -c "##INFO=.*ID=AC" || echo "0")

if [[ "$AN_COUNT" -eq 0 || "$AC_COUNT" -eq 0 ]]; then
    error "AN/AC tags not found in final output. This is required for LPA prediction."
else
    log "AN/AC tags successfully added to final output"
fi

# Check for multiallelic sites
N_MULTIALLELIC=$(bcftools view -H "$FINAL_BCF" | awk '$5 ~ /,/' | wc -l)
if [[ "$N_MULTIALLELIC" -gt 0 ]]; then
    warn "Found $N_MULTIALLELIC multiallelic sites. Consider splitting with bcftools norm."
fi

# Get final sample count with missing data
FINAL_MISSING=$(grep "^PSC" "$OUTPUT_DIR/final_stats.txt" | awk '$14>0' | wc -l 2>/dev/null || echo "0")

# Create summary report
log "Creating summary report..."

# Check if impute5 was used
IMPUTE5_USED="No"
if [[ -f "$TEMP_DIR/${BASENAME}.impute5.bcf" ]]; then
    IMPUTE5_USED="Yes"
fi

cat > "$OUTPUT_DIR/preprocessing_summary.txt" << EOF
LPA Genotype Preprocessing Summary
==================================
Date: $(date)
Input: $INPUT_VCF
Output: $FINAL_BCF
Threads used: $THREADS
Temporary directory: $TEMP_DIR

Tools used:
- bcftools: $(which bcftools 2>/dev/null || echo "path not found")
- ShapeIt5: ${SHAPEIT5_BIN:-not found}
- impute5: ${IMPUTE5_BIN:-not found}

Original chromosome naming: $CHR
Standardized to: chr6
Model sites required: $N_MODEL_SITES
Initial extracted sites: $N_EXTRACTED (${COVERAGE_PCT}%)
Samples with missing data (before imputation): $N_SAMPLES_WITH_MISSING / $TOTAL_SAMPLES
ShapeIt5 phasing/imputation performed: $NEED_IMPUTATION
Sites after ShapeIt5: ${N_SHAPEIT5_FINAL:-N/A}
impute5 formal imputation performed: $IMPUTE5_USED
Final sites: $(bcftools view -H "$FINAL_BCF" | wc -l)
Total samples: $(bcftools query -l "$FINAL_BCF" | wc -l)
Multiallelic sites: $N_MULTIALLELIC
Samples with missing data (final): $FINAL_MISSING
AN/AC tags present: $(if [[ "$AN_COUNT" -gt 0 && "$AC_COUNT" -gt 0 ]]; then echo "Yes"; else echo "No"; fi)
EOF

cat "$OUTPUT_DIR/preprocessing_summary.txt"

log "Preprocessing completed successfully!"
log "Ready-to-use BCF file with AN/AC tags: $FINAL_BCF"

echo $FINAL_BCF
