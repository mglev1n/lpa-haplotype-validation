# LPA Prediction Validation Pipeline

This pipeline compares predicted Lipoprotein(a) values from genotype data against measured values across various genetic ancestry groups.

## Overview

The pipeline:
1. Processes phased genotype data from VCF files; imputation will automatically be performed to fill-in missing LPA variants using SHAPEIT5 (if necessary)
2. Uses the `lpapredictr` R package to predict LPA levels from the phased genotype data
3. Compares predictions to measured LPA values, and compares diagnostic potential vs. usual care among the full cohort
4. Generates validation statistics and visualizations

**Expected Runtime:** 1 to 3 hours depending on dataset size and compute resources.

## Prerequisites

- **For containerized runs (recommended):** Singularity/Apptainer or Docker
- **For local runs:** R 4.3.0+, bcftools, shapeit5, and reference panels
- **System requirements:** 4+ CPU cores, 8+ GB RAM recommended
- **Storage:** ~5-10 GB free space for intermediate files

## Directory Structure

### Initial Setup

```
lpa-validation-project/
├── input/
│   ├── genotypes.vcf.gz        # Your phased genotype data
│   ├── genotypes.vcf.gz.tbi    # Index file (auto-generated if missing)
│   └── measured.csv            # Measured Lp(a) values + demographics
└── (pipeline files will be copied here when run)
```

## Input File Requirements

The pipeline requires two input files in the `input/` directory:

### 1. Genotype VCF File: `input/genotypes.vcf.gz`

> [!IMPORTANT]  
> This file should include the **full genotyped biobank cohort** (both with and without measured Lp(a) values)

**File specifications:**
- **Format:** Phased VCF, bgzip compressed (.vcf.gz)
- **Region:** Chromosome 6, positions ~159.5-161.7 Mb (GRCh38)
- **Phasing:** Must be phased (if using imputed data, this is typically already phased)
- **Size expectation:** 100 MB - 10 GB depending on cohort size

**Preparation example using PLINK2:**
```bash
# Extract LPA region and export to VCF
plink2 --pfile your_dataset \
  --chr 6 \
  --from-bp 159500000 \
  --to-bp 161700000 \
  --export vcf bgz id-paste=iid vcf-dosage='HDS' \
  --output-chr chr26 \
  --out input/genotypes

# Verify the file was created correctly
bcftools index input/genotypes.vcf.gz
bcftools view -H input/genotypes.vcf.gz | head -5
```

**Validation checklist:**
- [ ] File is phased (check for `|` separators in GT field, not `/`)
- [ ] Sample IDs match those in `measured.csv`
- [ ] Contains chromosome 6 LPA region variants
- [ ] File is properly indexed (`.tbi` or `.csi` file present)

### 2. Measured LPA Values + Demographics: `input/measured.csv`

> [!IMPORTANT]  
> Include the **full genotyped biobank cohort** - the `lpa_nmol_max` column can be `NA` for samples without measured Lp(a) values. There should be 1 row per participant. For individuals with more than 1 available measurement, select the maximum (convert to nmol/L first).

> [!CAUTION]
> **Unit conversion required:** Lp(a) values measured in mg/dL must be converted to nmol/L by multiplying by 2.15

**LPA columns:**
- `ID`: Sample identifiers (must exactly match VCF sample IDs)
- `lpa_nmol_max`: Measured LPA values in nmol/L (convert from mg/dL × 2.15). Use `NA` for unmeasured samples. 
  - If multiple measurements exist for a participant, use the maximum value.
  - If no measurements are available, leave as `NA`.
- `GIA`: Genetic ancestry group - must be one of: `"EUR"`, `"AFR"`, `"EAS"`, `"CSA"`, `"AMR"`, `"MID"`, `"Other"`

**Demographic columns**:
- `age`: Age at enrollment (numeric, years)
- `sex`: Sex (`"Male"`, `"Female"`, or `"Other"`)
- `IHD`: Ischemic heart disease (TRUE/FALSE, ≥2 instances of phecode 411)
- `HTN`: Hypertension (TRUE/FALSE, ≥2 instances of phecode 401)
- `HLD`: Hyperlipidemia (TRUE/FALSE, ≥2 instances of phecode 272.1)
- `HF`: Heart failure (TRUE/FALSE, ≥2 instances of phecode 428)
- `CKD`: Chronic kidney disease (TRUE/FALSE, ≥2 instances of phecode 585.3)
- `DM`: Diabetes mellitus (TRUE/FALSE, ≥2 instances of phecode 250)
- `AF`: Atrial fibrillation (TRUE/FALSE, ≥2 instances of phecode 427.2)
- `PVD`: Peripheral vascular disease (TRUE/FALSE, ≥2 instances of phecode 443)
- `IS`: Ischemic stroke (TRUE/FALSE, ≥2 instances of phecode 433.21)

**Example file structure:**
```csv
ID,lpa_nmol_max,GIA,age,sex,IHD,HTN
SAMPLE001,120.5,EUR,65,Male,TRUE,FALSE
SAMPLE002,NA,AFR,45,Female,FALSE,TRUE
SAMPLE003,250.3,EAS,55,Male,FALSE,FALSE
```

**Validation checklist:**
- [ ] Sample IDs in CSV exist in the VCF file
- [ ] `lpa_nmol_max` values are in nmol/L (not mg/dL)
- [ ] `GIA` values use only valid ancestry codes

## Running the Pipeline

### Option 1: Singularity/Apptainer (Recommended for HPC)

**Step 1:** Prepare your workspace
```bash
# Create project directory and navigate to it
mkdir lpa-validation-project
cd lpa-validation-project

# Create input directory
mkdir -p input

# Copy your prepared files
cp /path/to/your/genotypes.vcf.gz input/
cp /path/to/your/measured.csv input/

# Verify files are present and readable
ls -la input/
```

**Step 2:** Download and run the container
```bash
# Download the container (one-time setup)
singularity pull oras://ghcr.io/mglev1n/lpa-validation-singularity:latest

# Run the pipeline
singularity run --pwd /work -B $(pwd):/work lpa-validation-singularity_latest.sif
```

**For HPC clusters with job schedulers:**
```bash
# Example LSF submission (4 cores, 16 GB memory)
bsub -n 4 -M 16000 -R "rusage[mem=16000]" \
  "module load singularity; \
   singularity run --pwd /work -B $(pwd):/work lpa-validation-singularity_latest.sif"
```

**Troubleshooting tips:**
- If you get temp directory errors, add: `-B $(pwd)/tmp:/tmp`
- If you get permission errors, ensure your current directory is writable

### Option 2: Docker

```bash
# Navigate to your project directory (with input/ subdirectory)
cd /path/to/lpa-validation-project

# Run the pipeline
docker run --rm -v $(pwd):/work ghcr.io/mglev1n/lpa-validation:latest
```

### Option 3: Local Installation (Advanced Users)

> [!CAUTION]
> Local installation requires additional system dependencies (bcftools, SHAPEIT5, reference panels) that are pre-installed in the containers. Only recommended if you're familiar with these tools.

**Step 1:** Clone and setup
```bash
git clone https://github.com/mglev1n/lpa-validation.git
cd lpa-validation

# Install R dependencies
Rscript -e 'install.packages("renv")'
Rscript -e 'renv::restore()'
```

**Step 2:** Prepare input files
```bash
mkdir -p input
cp your_genotypes.vcf.gz input/genotypes.vcf.gz
cp your_measured.csv input/measured.csv
```

**Step 3:** Run the pipeline
```bash
# Option A: Using R directly
Rscript -e 'targets::tar_make()'

# Option B: Using R interactively
R
> library(targets)
> tar_make()
> quit()
```

## Output Files

All results are saved to the `Results/` directory:

### Main Report
- `lpa_validation_report.html` - **Primary output:** Comprehensive validation report

### Archive
- `lpa_validation_results.zip` - All results (plots/tables) in a single downloadable archive


## License

MIT
