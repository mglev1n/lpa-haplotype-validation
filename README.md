# LPA Prediction Validation Pipeline

This pipeline validates LPA (Lipoprotein(a)) prediction models by comparing predicted values from genotype data against measured values across various genetic ancestry groups.

## Overview

The pipeline:
1. Processes phased genotype data from VCF files
2. Uses the `lpapredictr` package to predict LPA levels
3. Compares predictions to measured LPA values
4. Generates validation statistics and visualizations

## Input File Requirements

The pipeline requires two input files:
1. `input/genotypes.vcf.gz` containing phased genotypes for the LPA gene region
2. `input/measured.csv` containing measured LPA values and covariates

### 1. Genotype VCF File

> [!IMPORTANT]  
> This file should include the full genotyped biobank cohort

Place a phased VCF file at `input/genotypes.vcf.gz` containing phased genotypes for the LPA gene region. 

**Preparation Example using PLINK2:**
```bash
plink2 --pfile your_dataset --chr 6 --from-bp 159500000 --to-bp 161700000 \
  --export vcf bgz id-paste=iid vcf-dosage='HDS' \
  --out input/genotypes
```

**Requirements:**
- File must be phased and contain the LPA gene region; if you're using imputed genotype data, this should already be phased
- Chromosome 6, approximately position 159.5-161.7 Mb (GRCh38)
- File must be gzipped (.vcf.gz)
- Sample IDs should match those in the `input/measured.csv` file containing Lp(a) measurements and covariates.

### 2. Measured LPA Values + Comorbidities

> [!IMPORTANT]  
> This file should include the full genotyped biobank cohort - the `lpa_nmol_max` column can be `NA` for samples without values

> [!CAUTION]
> Lp(a) levels may be measured in mg/dL or nmol/L. For this study, all measurements should be converted to nmol/L. For values measured in mg/dL, multiply by 2.15 to convert to nmol/L.

Place a CSV file at `input/measured.csv` containing measured LPA values in nmol/L and covariates; convert from mg/dL by multiplying by 2.15. Prevalent clinical comorbidities can be derived from a standard dataframe of phecodes, or alternative local phenotype definitions if more readily available.

**Required columns:**
- `ID`: Sample identifiers matching the VCF file
- `lpa_nmol_max`: Measured LPA values in nmol/L; this should be the maximum value if multiple measurements have been performed. Convert from mg/dL by multiplying by 2.15. Cells in this column will be empty for samples without values.
- `GIA`: Genetic ancestry group (must be one of "EUR", "AFR", "EAS", "CSA", "AMR", "MID", "Other")

**Covariate/demographic columns:**
- `age` (numeric): Age at enrollment in years
- `sex` (Male/Female/Other): Sex of participant
- `IHD` (TRUE/FALSE): 2 or more instances of phecode `411` (Ischemic Heart Disease)
- `HTN` (TRUE/FALSE): 2 or more instances of phecode `401` (Hypertension)
- `HLD` (TRUE/FALSE): 2 or more instances of phecode `272.1` (Hyperlipidemia)
- `HF` (TRUE/FALSE): 2 or more instances of phecode `428` (Heart Failure)
- `CKD` (TRUE/FALSE): 2 or more instances of phecode `585.3` (Chronic Kidney Disease)
- `DM` (TRUE/FALSE): 2 or more instances of phecode `250` (Diabetes Mellitus)
- `AF` (TRUE/FALSE): 2 or more instances of phecode `427.2` (Atrial Fibrillation)
- `PVD` (TRUE/FALSE): 2 or more instances of phecode `443` (Peripheral Vascular Disease)
- `IS` (TRUE/FALSE): 2 or more instances of phecode `433.21` (Ischemic Stroke)

**Example:**
```
ID,lpa_nmol_max,GIA,age,sex
Sample1,125.6,EUR,54,F
Sample2,34.2,AFR,63,M
...
```

## Running the Pipeline

There are three options for running the pipeline:

### Option 1: Running with Singularity/Apptainer (Recommended for HPC)

1. Prepare your directory structure:
   ```bash
   mkdir -p input
   cp /path/to/your/genotypes.vcf.gz input/
   cp /path/to/your/measured.csv input/
   ```

2. Run the container in the current directory:
   ```bash
   singularity run --pwd /work -B $(pwd):/work lpa-validation_latest.sif
   ```

   All results will be written to the `Results` directory in your current working directory.

3. For HPC environments:
   ```bash
   # Go to your project directory
   cd /project/path/lpa-validation
   
   # Create input directory and add files
   mkdir -p input
   cp /path/to/your/genotypes.vcf.gz input/
   cp /path/to/your/measured.csv input/
   
   # Run with Singularity/Apptainer
   singularity run --pwd /work -B $(pwd):/work /path/to/lpa-validation_latest.sif
   ```

   If you encounter temp directory issues, add: `-B $(pwd)/tmp:/tmp`

### Option 2: Using Docker

1. Prepare your directory structure as above.

2. Run the Docker container:
   ```bash
   docker run --rm -v $(pwd):/work ghcr.io/mglev1n/lpa-validation:latest
   ```

### Option 3: Running Locally from GitHub

1. Clone the repository:
   ```bash
   git clone https://github.com/mglev1n/lpa-validation.git
   cd lpa-validation
   ```

2. Install R (version 4.3.0 or higher): [https://cran.r-project.org/](https://cran.r-project.org/)

3. Install renv:
   ```r
   install.packages("renv")
   ```

4. Restore dependencies from the lockfile:
   ```bash
   Rscript -e 'renv::restore()'
   ```

5. Place your input files:
   ```bash
   mkdir -p input
   cp your_genotypes.vcf.gz input/genotypes.vcf.gz
   cp your_measured.csv input/measured.csv
   ```

6. Run the pipeline:
   ```bash
   ./run.sh
   ```

## Output Files

The pipeline generates various outputs in the `Results/` directory:

- `lpa_validation_report.html` - Comprehensive HTML report
- `demographics_overall.html` - Demographics table
- `demographics_by_ancestry.html` - Demographics by ancestry group
- `validation_plot.pdf` - Plot of predicted vs measured values
- `distribution_plot.pdf` - Distribution of values by ancestry
- `performance_plot.pdf` - Performance metrics visualization
- `nnt_plot.pdf` - Number Needed to Test visualization
- CSV files with detailed results and metrics
- `lpa_validation_results.zip` - All results in a single ZIP archive

## Troubleshooting

- **Input validation error**: Check that your input files match the required format
- **Missing dependencies**: If running locally, ensure `renv::restore()` completed successfully
- **Memory issues**: For large datasets, ensure adequate memory is available
- **HPC issues**: For Singularity, make sure your current directory is writable
- **Temp directory issues**: Add `-B $(pwd)/tmp:/tmp` to your Singularity command

## License

MIT
