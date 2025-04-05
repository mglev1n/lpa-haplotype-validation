# LPA Prediction Validation Pipeline

This pipeline validates LPA (Lipoprotein(a)) prediction models by comparing predicted values from genotype data against measured values across various genetic ancestry groups.

## Overview

The pipeline:
1. Processes phased genotype data from VCF files
2. Uses the `lpapredictr` package to predict LPA levels
3. Compares predictions to measured LPA values
4. Generates comprehensive validation statistics and visualizations
5. Creates a detailed HTML report of results

## Input File Requirements

### 1. Genotype VCF File

Place a phased VCF file at `input/genotypes.vcf.gz` containing genotypes for the LPA gene region.

**Preparation Example using PLINK2:**
```bash
plink2 --pfile your_dataset --chr 6 --from-bp 159500000 --to-bp 161700000 \
  --export vcf bgz id-paste=iid vcf-dosage='HDS' \
  --out input/genotypes
```

**Requirements:**
- File must be phased and contain the LPA gene region
- Chromosome 6, approximately position 159.5-161.7 Mb (GRCh38)
- File must be gzipped (.vcf.gz)
- Sample IDs should match those in the measured file

### 2. Measured LPA Values

Place a CSV file at `input/measured.csv` containing measured LPA values.

**Required columns:**
- `ID`: Sample identifiers matching the VCF file
- `lpa_nmol_max`: Measured LPA values in nmol/L
- `GIA`: Genetic ancestry group (e.g., "EUR", "AFR", "EAS", "CSA", "AMR", "MID", "Other")

**Optional columns (if available):**
- `age`: Age in years
- `sex`: Sex of participant
- Comorbidities (binary 0/1): `IHD`, `HTN`, `HLD`, `HF`, `CKD`, `DM`, `AF`, `PVD`, `IS`

**Example:**
```
ID,lpa_nmol_max,GIA,age,sex
Sample1,125.6,EUR,54,F
Sample2,34.2,AFR,63,M
...
```

## Running the Pipeline

There are two options for running the pipeline:

### Option 1: Using Docker (Recommended)

1. Install Docker on your system: [https://docs.docker.com/get-docker/](https://docs.docker.com/get-docker/)

2. Pull the Docker image:
   ```bash
   docker pull ghcr.io/mglev1n/lpa-validation:latest
   ```

3. Create a directory for your input and output:
   ```bash
   mkdir -p ~/lpa-validation/input ~/lpa-validation/Results
   ```

4. Place your input files in the input directory:
   ```bash
   cp your_genotypes.vcf.gz ~/lpa-validation/input/genotypes.vcf.gz
   cp your_measured.csv ~/lpa-validation/input/measured.csv
   ```

5. Run the container:
   ```bash
   docker run --rm -v ~/lpa-validation:/data ghcr.io/mglev1n/lpa-validation:latest
   ```

6. Check the Results directory for outputs:
   ```bash
   ls -la ~/lpa-validation/Results/
   ```

### Option 2: Running Locally from GitHub

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

   The `renv.lock` file ensures that exact package versions are used, maintaining reproducibility across different environments.

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
- **Memory issues**: Docker may need additional memory allocation for large datasets

## License

MIT
