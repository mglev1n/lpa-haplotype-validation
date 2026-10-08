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

# Function to extract model sites from a VCF/BCF.
# Writes two files:
#   - a genotype BCF with INFO removed (GT only), used for LPA prediction
#   - a sites-only BCF that keeps the full INFO field, which carries the
#     imputation quality measures (e.g. R2 from upstream imputation, INFO from impute5)
extract_model_sites() {
    local input_file=$1
    local out_gt_bcf=$2
    local out_sites_bcf=$3
    local isec_bcf="${out_gt_bcf%.bcf}.isec.bcf"

    bcftools isec -c none -r "$EXTRACTION_REGION" -n=2 -w1 -Ob -o "$isec_bcf" \
        "$input_file" \
        "$SITES_VCF"
    bcftools index -f "$isec_bcf"

    bcftools view -G -Ob -o "$out_sites_bcf" "$isec_bcf"
    bcftools index -f "$out_sites_bcf"

    bcftools annotate -x INFO,^FMT/GT -Ob -o "$out_gt_bcf" "$isec_bcf"
    bcftools index -f "$out_gt_bcf"
}

# Function to add a prefix to every INFO tag in a sites-only file, so that
# quality measures from different sources (input data, impute5) can be stored
# in one file without name clashes. AN/AC are dropped; they describe allele
# counts, not quality, and are recomputed in the genotype output.
prefix_info_tags() {
    local input_file=$1
    local output_file=$2
    local prefix=$3
    local rename_file="${output_file}.rename.txt"
    local info_tags drop_tags

    info_tags=$(bcftools view -h "$input_file" | sed -n 's/^##INFO=<ID=\([^,]*\),.*/\1/p')

    # Drop AN/AC only if present (avoids bcftools warnings for undefined tags)
    drop_tags=$(echo "$info_tags" | awk '$1 == "AN" || $1 == "AC" {printf "%sINFO/%s", sep, $1; sep=","}')

    echo "$info_tags" \
        | awk -v p="$prefix" 'NF && $1 != "AN" && $1 != "AC" {print "INFO/"$1"\t"p$1}' \
        > "$rename_file"

    if [[ -n "$drop_tags" ]]; then
        bcftools annotate -x "$drop_tags" -Ob -o "${output_file}.tmp.bcf" "$input_file"
    else
        bcftools view -Ob -o "${output_file}.tmp.bcf" "$input_file"
    fi

    if [[ -s "$rename_file" ]]; then
        bcftools annotate --rename-annots "$rename_file" -Ob -o "$output_file" "${output_file}.tmp.bcf"
    else
        mv "${output_file}.tmp.bcf" "$output_file"
    fi
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
TEMP_FINAL_BCF="$TEMP_DIR/${BASENAME}.final.bcf"
INPUT_SITES_BCF="$TEMP_DIR/${BASENAME}.input.sites.bcf"
SHAPEIT5_SITES_BCF="$TEMP_DIR/${BASENAME}.shapeit5.sites.bcf"
IMPUTE5_SITES_BCF="$TEMP_DIR/${BASENAME}.impute5.sites.bcf"

# Set output filename - use custom name if provided, otherwise use default
if [[ -n "$OUTPUT_FILENAME" ]]; then
    FINAL_BCF="$OUTPUT_DIR/$OUTPUT_FILENAME"
else
    FINAL_BCF="$OUTPUT_DIR/${BASENAME}.processed.bcf"
fi

# Imputation quality outputs, named after the final genotype file
QUALITY_VCF="${FINAL_BCF%.*}.imputation_quality.vcf.gz"
QUALITY_TSV="${FINAL_BCF%.*}.imputation_quality.tsv"

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
# INPUT_SITES_BCF keeps the INFO field of the input data at the model sites
extract_model_sites "$WORKING_VCF" "$EXTRACTED_BCF" "$INPUT_SITES_BCF"
FINAL_STAGE_SITES_BCF="$INPUT_SITES_BCF"

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

    bcftools index -f "$IMPUTED_BCF"

    # Re-extract after imputation (before filling tags)
    log "Re-extracting variants after ShapeIt5 imputation..."
    extract_model_sites "$IMPUTED_BCF" "$TEMP_FINAL_BCF" "$SHAPEIT5_SITES_BCF"
    FINAL_STAGE_SITES_BCF="$SHAPEIT5_SITES_BCF"

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

        # Run impute5 with the phased data from shapeit5
        log "Running impute5 with $THREADS threads..."
        "$IMPUTE5_BIN" \
            --h "$REFERENCE_PANEL" \
            --m "$GENETIC_MAP" \
            --g "$IMPUTED_BCF" \
            --r "$REGION" \
            --buffer-region "$REGION" \
            --o "$IMPUTE5_BCF" \
            --threads "$THREADS" 2>&1 | tee -a "$LOG_FILE"

        if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
            error "impute5 imputation failed"
        fi

        # bcftools index -f "$IMPUTE5_BCF"

        # Re-extract after impute5 (before filling tags)
        log "Re-extracting variants after impute5..."
        # IMPUTE5_SITES_BCF keeps the impute5 INFO field (imputation quality)
        extract_model_sites "$IMPUTE5_BCF" "$TEMP_FINAL_BCF" "$IMPUTE5_SITES_BCF"
        FINAL_STAGE_SITES_BCF="$IMPUTE5_SITES_BCF"

        log "impute5 imputation completed"
    fi

    # Fill AN/AC tags for final extracted data only
    log "Filling AN/AC tags for final extracted sites..."
    fill_tags "$TEMP_FINAL_BCF" "$FINAL_BCF" "final extracted data"

    # Verify results after all imputation steps
    N_FINAL=$(bcftools view -H "$FINAL_BCF" | wc -l)
    FINAL_COVERAGE_PCT=$(awk "BEGIN {printf \"%.1f\", ($N_FINAL/$N_MODEL_SITES)*100}")

    log "Final sites: $N_FINAL"
    log "Final coverage: $FINAL_COVERAGE_PCT%"

    if [[ "$N_FINAL" -lt "$N_MODEL_SITES" ]]; then
        warn "Still missing $(($N_MODEL_SITES - $N_FINAL)) sites after all imputation attempts"
    else
        log "All required model sites successfully recovered"
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
    log "No imputation performed. Filling AN/AC tags for extracted data and using as final."
    # Fill AN/AC tags for the extracted data
    fill_tags "$EXTRACTED_BCF" "$FINAL_BCF" "final extracted data"
    cp "$TEMP_DIR/extracted_stats.txt" "$OUTPUT_DIR/final_stats.txt"
    # Set N_SHAPEIT5_FINAL for summary report
    N_SHAPEIT5_FINAL="N/A (no imputation needed)"
fi

# Step 8: Write imputation quality measures for the final model sites
# The sites-only VCF contains one record per site in the final genotype output:
#   - INPUT_*  : INFO tags from the input data (e.g. INPUT_R2 from upstream imputation)
#   - IMPUTE5_*: INFO tags from impute5 (only if impute5 ran)
#   - IN_INPUT : flag set when the site was present in the input data; sites
#                without the flag were added by impute5
log "Writing imputation quality measures..."
QUALITY_BASE_BCF="$TEMP_DIR/quality.base.bcf"
INPUT_SITES_PREFIXED_BCF="$TEMP_DIR/input.sites.prefixed.bcf"
QUALITY_HEADER="$TEMP_DIR/quality.header.txt"

prefix_info_tags "$INPUT_SITES_BCF" "$INPUT_SITES_PREFIXED_BCF" "INPUT_"

if [[ "$FINAL_STAGE_SITES_BCF" == "$IMPUTE5_SITES_BCF" ]]; then
    prefix_info_tags "$IMPUTE5_SITES_BCF" "$QUALITY_BASE_BCF" "IMPUTE5_"
else
    # ShapeIt5 phasing does not produce quality measures; keep only the site list
    bcftools annotate -x INFO -Ob -o "$QUALITY_BASE_BCF" "$FINAL_STAGE_SITES_BCF"
    bcftools index -f "$QUALITY_BASE_BCF"
fi

echo '##INFO=<ID=IN_INPUT,Number=0,Type=Flag,Description="Site present in the input genotypes; sites without this flag were imputed by impute5">' \
    > "$QUALITY_HEADER"

bcftools annotate -a "$INPUT_SITES_PREFIXED_BCF" -c INFO -m +IN_INPUT -h "$QUALITY_HEADER" \
    -Oz -o "$QUALITY_VCF" "$QUALITY_BASE_BCF"
bcftools index -t -f "$QUALITY_VCF"

# Flat table of the same measures, one column per INFO tag. Flag tags are
# written as 1 (set) or 0 (not set); missing values are written as NA.
QUALITY_TAGS=$(bcftools view -h "$QUALITY_VCF" | sed -n 's/^##INFO=<ID=\([^,]*\),.*/\1/p' | grep -v -x "IN_INPUT" || true)
QUALITY_FLAG_TAGS=$(bcftools view -h "$QUALITY_VCF" | sed -n 's/^##INFO=<ID=\([^,]*\),.*Type=Flag.*/\1/p' | tr '\n' ' ')
QUALITY_FORMAT='%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/IN_INPUT'
QUALITY_COLUMNS='CHROM\tPOS\tID\tREF\tALT\tIN_INPUT'
for tag in $QUALITY_TAGS; do
    QUALITY_FORMAT="${QUALITY_FORMAT}\t%INFO/${tag}"
    QUALITY_COLUMNS="${QUALITY_COLUMNS}\t${tag}"
done
{
    echo -e "$QUALITY_COLUMNS"
    bcftools query -f "${QUALITY_FORMAT}\n" "$QUALITY_VCF"
} | awk -F '\t' -v OFS='\t' -v flags="$QUALITY_FLAG_TAGS" '
    BEGIN { n = split(flags, f, " "); for (i = 1; i <= n; i++) is_flag[f[i]] = 1 }
    NR == 1 { for (i = 1; i <= NF; i++) col_is_flag[i] = ($i in is_flag); print; next }
    {
        for (i = 6; i <= NF; i++) {
            if (col_is_flag[i]) { $i = ($i == "1") ? 1 : 0 }
            else if ($i == ".") { $i = "NA" }
        }
        print
    }' > "$QUALITY_TSV"

N_IMPUTED_BY_IMPUTE5=$(bcftools view -H -i 'INFO/IN_INPUT=0' "$QUALITY_VCF" | wc -l)
log "Imputation quality VCF: $QUALITY_VCF"
log "Imputation quality table: $QUALITY_TSV"
log "Sites added by impute5: $N_IMPUTED_BY_IMPUTE5"

# Step 9: Final validation
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
Sites added by impute5: $N_IMPUTED_BY_IMPUTE5
Imputation quality (sites-only VCF): $QUALITY_VCF
Imputation quality (table): $QUALITY_TSV
EOF

cat "$OUTPUT_DIR/preprocessing_summary.txt"

log "Preprocessing completed successfully!"
log "Ready-to-use BCF file with AN/AC tags: $FINAL_BCF"

echo $FINAL_BCF
