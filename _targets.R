# _targets.R ----------------------------------------------------------------

# Load Packages ------------------------------------------------------------
library(targets)
library(tarchetypes)
library(tidyverse)
library(crew)

# Controller Setup ---------------------------------------------------------
controller_local <- crew_controller_local(
  name = "LPA-haplotype_local",
  workers = parallelly::availableCores(),
  seconds_idle = 10
)

# Global Options -----------------------------------------------------------
tar_option_set(
  packages = c(
    "tidyverse", "targets", "tarchetypes", "lpapredictr", "qs", "gt",
    "crew", "tidymodels", "patchwork", "vroom", "ggsci"
  ),
  controller = crew::crew_controller_group(controller_local),
  resources = tar_resources(
    crew = tar_resources_crew(controller = "LPA-haplotype_local")
  ),
  format = "qs",
  storage = "worker",
  retrieval = "worker"
)

# Pipeline Configuration ---------------------------------------------------
tar_option_set(
  error = "continue" # Continue pipeline despite errors in individual targets
)

# Create Directories -------------------------------------------------------
if (!dir.exists("Results")) {
  dir.create("Results")
}

if (!dir.exists("rmarkdown")) {
  dir.create("rmarkdown")
}

# Define the pipeline ------------------------------------------------------
list(
  # Input File Paths -------------------------------------------------------
  tar_target(
    vcf_file,
    "input/genotypes.vcf.gz",
    format = "file",
    description = "Path to VCF file containing phased genotypes for LPA region"
  ),
  tar_target(
    measured_file,
    "input/measured.csv",
    format = "file",
    description = "Path to CSV file containing measured LPA values and demographics"
  ),

  # Validation Checks ------------------------------------------------------
  tar_target(
    validate_inputs,
    {
      # Check if VCF file exists
      if (!file.exists(vcf_file)) {
        stop(
          "VCF file not found: ", vcf_file,
          "\nPlease place your VCF file at input/genotypes.vcf.gz"
        )
      }

      # Check if measured file exists and has required columns
      if (!file.exists(measured_file)) {
        stop(
          "Measured LPA file not found: ", measured_file,
          "\nPlease place your measured LPA file at input/measured.csv"
        )
      }

      df <- vroom::vroom(measured_file, show_col_types = FALSE)
      required_cols <- c("ID", "lpa_nmol_max", "GIA")
      missing_cols <- setdiff(required_cols, names(df))
      if (length(missing_cols) > 0) {
        stop(
          "Missing required columns in measured file: ",
          paste(missing_cols, collapse = ", ")
        )
      }

      # Return TRUE if validation passes
      TRUE
    },
    description = "Validate that input files exist and have the required format"
  ),

  # VCF File Processing ----------------------------------------------------

  tar_file(genotypes_clean,
   {
     # Run lpa_processing.sh to impute missing variants + extract model-relevant SNPs
     processx::run("Scripts/lpa_processing.sh",
                   args = c(
                     "-i", vcf_file,
                     "-o", "input",
                     "-x", "Resources/Lpa_hap_model.flank100.sites.vcf.gz",
                     "-r", "Resources/ALL.chr6.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz",
                     "-m", "Resources/chr6.b38.gmap.gz",
                     "-s", "/usr/local/bin/phase_common_static"
                   ),
                   echo = TRUE )
   },
   description = "Process VCF file to impute missing variants and extract relevant SNPs"
  ),

  # LPA Prediction ---------------------------------------------------------
  tar_target(
    lpa_predictions,
    {
      # Use lpapredictr package to predict LPA levels
      predictions <- lpapredictr::predict_lpa(vcf_file = genotypes_clean)
      predictions
    },
    description = "Run LPA prediction on phased genotypes from VCF file"
  ),


  ## Read Measured LPA Values -----------------------------------------------
  tar_target(
    lpa_measured,
    {
      # Read the CSV file with measured values
      vroom::vroom(measured_file, show_col_types = FALSE)
    },
    description = "Read measured LPA values and demographics from CSV file"
  ),
  tar_target(
    lpa_pred_summary,
    {
      lpa_predictions %>%
        group_by(ID) %>%
        summarize(
          lpa_pred_nm = sum(Lpa_pred_nM),
          lpa_pred_se = sqrt(sum(min_uncertainty_nM^2)),
          .groups = "drop"
        ) %>%
        mutate(
          lpa_pred_lci95 = lpa_pred_nm - 1.96 * lpa_pred_se,
          lpa_pred_uci95 = lpa_pred_nm + 1.96 * lpa_pred_se
        )
    },
    description = "Summarize LPA predictions by ID with confidence intervals"
  ),

  ## Combine Predicted and Measured Values ----------------------------------
  tar_target(
    lpa_pred_measured_df,
    {
      # Join predicted and measured values
      lpa_pred_summary %>%
        left_join(lpa_measured, by = "ID") %>%
        select(-ID)
    },
    description = "Combine predicted and measured LPA values for validation analysis"
  ),
  tar_target(
    lpa_validation_df,
    {
      lpa_pred_measured_df %>%
        filter(!is.na(lpa_pred_nm) & !is.na(lpa_nmol_max)) %>%
        add_count(GIA) %>%
        filter(n > 20)
    },
    description = "Filter combined dataset to individuals with both measured and predicted values"
  ),


  # Demographics Tables ----------------------------------------------------
  tar_target(
    lpa_demographics_overall,
    {
      # Get the available columns from the validation dataframe
      available_cols <- names(lpa_validation_df)

      # Define columns we'd ideally like to include (if they exist)
      ideal_demo_cols <- c(
        "age", "sex", "GIA", "lpa_nmol_max", "IHD", "HTN", "HLD", "HF",
        "CKD", "DM", "AF", "PVD", "IS"
      )

      # Find which columns actually exist in the data
      demo_cols <- intersect(ideal_demo_cols, available_cols)

      # Ensure we have at least lpa_nmol_max (our primary measure)
      if (!"lpa_nmol_max" %in% demo_cols) {
        warning("lpa_nmol_max column missing from validation data, demographics will be limited")
      }

      # Create minimal column list if few demographic columns exist
      if (length(demo_cols) < 3) {
        message("Few demographic columns available. Creating minimal summary table.")
        # Ensure we include at least GIA, even if not in ideal_demo_cols
        if ("GIA" %in% available_cols && !("GIA" %in% demo_cols)) {
          demo_cols <- c(demo_cols, "GIA")
        }
      }

      # Check if we have any columns to summarize
      if (length(demo_cols) == 0) {
        # Create a dummy table if no useful columns exist
        warning("No demographic columns found. Creating placeholder table.")
        df_raw <- tibble(
          variable = "cohort_size",
          label = "Cohort Size",
          Value = nrow(lpa_validation_df)
        )

        # Create a simple GT table
        gt_tbl <- df_raw %>%
          gt::gt() %>%
          gt::tab_header("Study Cohort Summary (Limited Data Available)")

        return(list(
          gt_table = gt_tbl,
          raw_data = df_raw,
          cohort = "Overall"
        ))
      }

      # Build label mappings for available columns
      label_list <- list(
        lpa_nmol_max = "Lp(a) (nmol/L)",
        GIA = "Genetic ancestry"
      )

      # Add conditional labels for other possible columns
      if ("age" %in% demo_cols) label_list$age <- "Age (years)"
      if ("sex" %in% demo_cols) label_list$sex <- "Sex"
      if ("IHD" %in% demo_cols) label_list$IHD <- "Ischemic heart disease"
      if ("HTN" %in% demo_cols) label_list$HTN <- "Hypertension"
      if ("HLD" %in% demo_cols) label_list$HLD <- "Hyperlipidemia"
      if ("HF" %in% demo_cols) label_list$HF <- "Heart failure"
      if ("CKD" %in% demo_cols) label_list$CKD <- "Chronic kidney disease"
      if ("DM" %in% demo_cols) label_list$DM <- "Diabetes"
      if ("AF" %in% demo_cols) label_list$AF <- "Atrial fibrillation"
      if ("PVD" %in% demo_cols) label_list$PVD <- "Peripheral vascular disease"
      if ("IS" %in% demo_cols) label_list$IS <- "Ischemic stroke"

      # Convert list to format expected by tbl_summary
      label_expr <- lapply(names(label_list), function(n) {
        as.formula(paste0(n, " ~ '", label_list[[n]], "'"))
      })

      # Create summary table, wrapped in tryCatch for extra safety
      tbl <- tryCatch(
        {
          lpa_validation_df %>%
            select(all_of(demo_cols)) %>%
            gtsummary::tbl_summary(
              missing = "no",
              label = label_expr
            )
        },
        error = function(e) {
          # Fallback if gtsummary fails
          warning(
            "Error in creating gtsummary table: ", e$message,
            "\nCreating simple summary instead."
          )

          # Create a simplified table manually
          basic_stats <- tibble(
            variable = "cohort_size",
            label = "Cohort Size",
            Value = nrow(lpa_validation_df)
          )

          # Add mean Lp(a) if available
          if ("lpa_nmol_max" %in% demo_cols) {
            basic_stats <- bind_rows(basic_stats, tibble(
              variable = "lpa_nmol_max",
              label = "Mean Lp(a) (nmol/L)",
              Value = sprintf("%.1f", mean(lpa_validation_df$lpa_nmol_max, na.rm = TRUE))
            ))
          }

          return(basic_stats)
        }
      )

      # Extract raw data frame depending on what type tbl is
      if (inherits(tbl, "tbl_summary")) {
        df_raw <- tbl$table_body %>%
          select(variable, label, stat_0 = `stat_0`) %>%
          rename(Value = stat_0)

        gt_table <- tbl %>% gtsummary::as_gt()
      } else {
        # tbl is already a tibble from the fallback
        df_raw <- tbl
        gt_table <- gt::gt(df_raw)
      }

      # Return both the formatted table and raw data
      list(
        gt_table = gt_table,
        raw_data = df_raw,
        cohort = "Overall"
      )
    },
    description = "Generate demographics table for the overall study population with robust error handling"
  ),

  ## Save Demographics Tables -----------------------------------------------
  tar_file(
    lpa_demographics_overall_html,
    {
      file_path <- "Results/demographics_overall.html"
      gt::gtsave(lpa_demographics_overall$gt_table, file_path)
      file_path
    },
    description = "Save overall demographics table to HTML file"
  ),
  tar_file(
    lpa_demographics_overall_csv,
    {
      file_path <- "Results/demographics_overall.csv"
      write.csv(lpa_demographics_overall$raw_data, file_path, row.names = FALSE)
      file_path
    },
    description = "Save overall demographics table to CSV file"
  ),


  ## Demographics By Ancestry -----------------------------------------------
  tar_target(
    lpa_demographics_by_ancestry,
    {
      # Get the available columns from the validation dataframe
      available_cols <- names(lpa_validation_df)

      # Verify GIA column exists
      if (!"GIA" %in% available_cols) {
        stop("GIA column is required for ancestry-stratified demographics but is missing")
      }

      # Define columns we'd ideally like to include (if they exist)
      ideal_demo_cols <- c(
        "age", "sex", "lpa_nmol_max", "IHD", "HTN", "HLD", "HF",
        "CKD", "DM", "AF", "PVD", "IS"
      )

      # Find which columns actually exist in the data
      demo_cols <- intersect(ideal_demo_cols, available_cols)

      # Ensure we have at least lpa_nmol_max (our primary measure)
      if (!"lpa_nmol_max" %in% demo_cols) {
        warning("lpa_nmol_max column missing from validation data, demographics will be limited")
      }

      # Create minimal column list if few demographic columns exist
      if (length(demo_cols) < 2) {
        message("Few demographic columns available. Creating minimal summary table.")
        # Add any available columns that might be useful
        extra_cols <- setdiff(available_cols, c(
          "GIA", demo_cols, "ID", "n", "lpa_pred_nm",
          "lpa_pred_se", "lpa_pred_lci95", "lpa_pred_uci95"
        ))
        if (length(extra_cols) > 0) {
          # Add up to 3 extra columns
          demo_cols <- c(demo_cols, head(extra_cols, 3))
        }
      }

      # Check if we have any columns to summarize
      if (length(demo_cols) == 0) {
        # Create a dummy table if no useful columns exist
        warning("No demographic columns found. Creating minimal table with counts by ancestry.")

        # Just count by ancestry
        ancestry_counts <- lpa_validation_df %>%
          count(GIA) %>%
          rename(Count = n)

        # Create a simple GT table
        gt_tbl <- ancestry_counts %>%
          gt::gt() %>%
          gt::tab_header("Study Cohort Counts by Ancestry (Limited Data Available)")

        # Create dummy data for the expected output format
        dummy_variable <- tibble(
          variable = "count",
          label = "Count",
          GIA = ancestry_counts$GIA,
          value = ancestry_counts$Count
        )

        return(list(
          gt_table = gt_tbl,
          raw_data = ancestry_counts,
          tidy_data = dummy_variable,
          cohort = "By Ancestry"
        ))
      }

      # Build label mappings for available columns
      label_list <- list()

      # Add conditional labels for possible columns
      if ("lpa_nmol_max" %in% demo_cols) label_list$lpa_nmol_max <- "Lp(a) (nmol/L)"
      if ("age" %in% demo_cols) label_list$age <- "Age (years)"
      if ("sex" %in% demo_cols) label_list$sex <- "Sex"
      if ("IHD" %in% demo_cols) label_list$IHD <- "Ischemic heart disease"
      if ("HTN" %in% demo_cols) label_list$HTN <- "Hypertension"
      if ("HLD" %in% demo_cols) label_list$HLD <- "Hyperlipidemia"
      if ("HF" %in% demo_cols) label_list$HF <- "Heart failure"
      if ("CKD" %in% demo_cols) label_list$CKD <- "Chronic kidney disease"
      if ("DM" %in% demo_cols) label_list$DM <- "Diabetes"
      if ("AF" %in% demo_cols) label_list$AF <- "Atrial fibrillation"
      if ("PVD" %in% demo_cols) label_list$PVD <- "Peripheral vascular disease"
      if ("IS" %in% demo_cols) label_list$IS <- "Ischemic stroke"

      # Convert list to format expected by tbl_summary
      label_expr <- lapply(names(label_list), function(n) {
        as.formula(paste0(n, " ~ '", label_list[[n]], "'"))
      })

      # Create stratified table, wrapped in tryCatch for extra safety
      tbl <- tryCatch(
        {
          lpa_validation_df %>%
            select(GIA, all_of(demo_cols)) %>%
            gtsummary::tbl_summary(
              by = GIA,
              missing = "no",
              label = label_expr
            )
        },
        error = function(e) {
          # Fallback if gtsummary fails
          warning(
            "Error in creating gtsummary table: ", e$message,
            "\nCreating simple summary by ancestry instead."
          )

          # Create a simplified table manually - counts and mean Lp(a) by ancestry
          ancestry_counts <- lpa_validation_df %>%
            group_by(GIA) %>%
            summarize(Count = n(), .groups = "drop")

          # Add mean Lp(a) if available
          if ("lpa_nmol_max" %in% demo_cols) {
            lpa_by_ancestry <- lpa_validation_df %>%
              group_by(GIA) %>%
              summarize(
                Mean_Lpa = sprintf("%.1f", mean(lpa_nmol_max, na.rm = TRUE)),
                .groups = "drop"
              )

            ancestry_counts <- left_join(ancestry_counts, lpa_by_ancestry, by = "GIA")
          }

          return(ancestry_counts)
        }
      )

      # Extract raw data frame and create tidy data depending on what type tbl is
      if (inherits(tbl, "tbl_summary")) {
        # Extract raw data frame for easier export
        df_raw <- tbl$table_body

        # Create tidy version by reshaping
        ancestry_vals <- unique(lpa_validation_df$GIA)

        # Transform to tidy format - convert from wide to long
        tidy_data <- df_raw %>%
          select(-label) %>%
          tidyr::pivot_longer(
            cols = starts_with("stat_"),
            names_to = "ancestry_code",
            values_to = "value"
          ) %>%
          # Extract ancestry code from column name (stat_0, stat_1, etc.)
          mutate(
            ancestry_index = as.integer(gsub("stat_", "", ancestry_code)),
            ancestry = ancestry_vals[ancestry_index + 1] # +1 because R is 1-indexed
          ) %>%
          select(variable, ancestry, value) %>%
          tidyr::drop_na()

        gt_table <- tbl %>% gtsummary::as_gt()
      } else {
        # tbl is already a tibble from the fallback
        df_raw <- tbl

        # Convert to tidy format
        tidy_data <- df_raw %>%
          tidyr::pivot_longer(
            cols = -GIA,
            names_to = "variable",
            values_to = "value"
          ) %>%
          rename(ancestry = GIA)

        gt_table <- gt::gt(df_raw)
      }

      # Return both the formatted table and raw data
      list(
        gt_table = gt_table,
        raw_data = df_raw,
        tidy_data = tidy_data,
        cohort = "By Ancestry"
      )
    },
    description = "Generate demographics table stratified by genetic ancestry with robust error handling"
  ),

  ## Save Ancestry Demographics ---------------------------------------------
  tar_file(
    lpa_demographics_by_ancestry_html,
    {
      file_path <- "Results/demographics_by_ancestry.html"
      gt::gtsave(lpa_demographics_by_ancestry$gt_table, file_path)
      file_path
    },
    description = "Save ancestry-stratified demographics table to HTML file"
  ),
  tar_file(
    lpa_demographics_by_ancestry_csv,
    {
      file_path <- "Results/demographics_by_ancestry.csv"
      write.csv(lpa_demographics_by_ancestry$raw_data, file_path, row.names = FALSE)
      file_path
    },
    description = "Save ancestry-stratified demographics table to CSV file (wide format)"
  ),
  tar_file(
    lpa_demographics_by_ancestry_tidy_csv,
    {
      file_path <- "Results/demographics_by_ancestry_tidy.csv"
      write.csv(lpa_demographics_by_ancestry$tidy_data, file_path, row.names = FALSE)
      file_path
    },
    description = "Save ancestry-stratified demographics table to CSV file (tidy/long format)"
  ),


  # Visualization Setup ----------------------------------------------------
  tar_target(
    ancestry_order,
    {
      # Get counts by ancestry
      ancestry_counts <- lpa_validation_df %>%
        count(GIA) %>%
        arrange(desc(n))

      # Move "Other" to the end if it exists
      if ("Other" %in% ancestry_counts$GIA) {
        other_row <- ancestry_counts %>% filter(GIA == "Other")
        ancestry_counts <- ancestry_counts %>%
          filter(GIA != "Other") %>%
          bind_rows(other_row)
      }

      # Create ordered factor
      ancestry_order <- ancestry_counts$GIA
      ancestry_counts
    },
    description = "Determine ordering of genetic ancestry groups for plots (by sample size, with 'Other' last)"
  ),
  tar_target(
    standard_ancestry_palette,
    {
      # Define valid super-population codes from 1000G+HGDP
      valid_ancestries <- c("AFR", "AMR", "CSA", "EAS", "EUR", "MID", "Other")

      # Check which ancestries are present in our data
      present_ancestries <- unique(lpa_validation_df$GIA)

      # Validate that all ancestries in the data are valid
      invalid_ancestries <- setdiff(present_ancestries, valid_ancestries)
      if (length(invalid_ancestries) > 0) {
        stop(
          "Invalid ancestry codes detected in GIA column: ",
          paste(invalid_ancestries, collapse = ", "),
          "\nValid codes are: ", paste(valid_ancestries, collapse = ", ")
        )
      }

      # Create consistent color palette using JAMA colors
      # This ensures the same ancestry always gets the same color across cohorts
      jama_palette <- ggsci::pal_jama()(7) # Get 7 colors from JAMA palette

      # Create named vector mapping ancestries to colors
      ancestry_colors <- setNames(
        jama_palette,
        valid_ancestries
      )

      # Return the complete palette (will subset as needed for each plot)
      ancestry_colors
    },
    description = "Create standardized color palette for genetic ancestry groups"
  ),

  # Validation Plots -------------------------------------------------------
  tar_target(
    plot_validation,
    {
      # Get the subset of colors for the ancestries in our data
      data_ancestries <- c("Overall", ancestry_order$GIA)
      ancestry_colors <- standard_ancestry_palette[names(standard_ancestry_palette) %in% ancestry_order$GIA]

      # Add "Overall" color (black)
      plot_colors <- c("Overall" = "black", ancestry_colors)

      # Function to create validation plot
      plot_val <- function(df, colors, ancestry_order) {
        # Ensure GIA is an ordered factor with correct ordering
        df <- df %>%
          mutate(GIA = factor(GIA, levels = ancestry_order))

        # Create labels and plot
        df %>%
          mutate(label = glue::glue("{GIA} (n = {n})")) %>%
          # Order facets by GIA factor levels
          mutate(label = factor(label, levels = paste0(levels(GIA), " (n = ", df$n[match(levels(GIA), df$GIA)], ")"))) %>%
          ggplot(aes(lpa_pred_nm, lpa_nmol_max)) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
          geom_jitter(alpha = 0.5) +
          geom_smooth(method = "lm", aes(color = GIA)) +
          facet_grid(~label) +
          scale_color_manual(values = colors) +
          guides(color = "none") +
          labs(y = "Measured Lp(a) (nmol/L)", x = "Predicted Lp(a) (nmol/L)") +
          theme_bw(base_size = 12) +
          theme(plot.margin = margin(1, 1, 1, 1, "cm"))
      }

      # Create overall plot
      p1 <- lpa_validation_df %>%
        mutate(GIA = "Overall") %>%
        add_count(GIA, name = "n") %>%
        plot_val(colors = plot_colors, ancestry_order = "Overall")

      # Create plot by ancestry
      p2 <- lpa_validation_df %>%
        plot_val(colors = plot_colors, ancestry_order = ancestry_order$GIA)

      # Return both components and the combined plot
      list(
        overall = p1,
        by_ancestry = p2,
        combined = (p1 / p2) + patchwork::plot_layout(heights = c(3, 1))
      )
    },
    description = "Generate validation plot comparing predicted vs measured LPA values"
  ),

  ## Save Validation Plots --------------------------------------------------
  tar_file(
    plot_validation_overall_qs,
    {
      file_path <- "Results/validation_plot_overall.qs"
      qs::qsave(plot_validation$overall, file_path)
      file_path
    },
    description = "Save overall validation plot as qs file for programmatic manipulation"
  ),
  tar_file(
    plot_validation_by_ancestry_qs,
    {
      file_path <- "Results/validation_plot_by_ancestry.qs"
      qs::qsave(plot_validation$by_ancestry, file_path)
      file_path
    },
    description = "Save by-ancestry validation plot as qs file for programmatic manipulation"
  ),
  tar_file(
    plot_validation_pdf,
    {
      file_path <- "Results/validation_plot.pdf"
      ggsave(file_path, plot_validation$combined, width = 12, height = 8)
      file_path
    },
    description = "Save combined validation plot as PDF"
  ),
  tar_file(
    plot_validation_qs,
    {
      file_path <- "Results/validation_plot.qs"
      qs::qsave(plot_validation$combined, file_path)
      file_path
    },
    description = "Save combined validation plot as qs file for programmatic manipulation"
  ),

  ## Distribution Plots -----------------------------------------------------
  tar_target(
    plot_distribution,
    {
      # Get the subset of colors for the ancestries in our data
      data_ancestries <- c("Overall", ancestry_order$GIA)
      ancestry_colors <- standard_ancestry_palette[names(standard_ancestry_palette) %in% ancestry_order$GIA]

      # Add "Overall" color (black)
      plot_colors <- c("Overall" = "black", ancestry_colors)

      # Function to create density plot
      plot_density <- function(df, colors, ancestry_order) {
        # Ensure GIA is an ordered factor with correct ordering
        df <- df %>%
          mutate(GIA = factor(GIA, levels = ancestry_order))

        # Create labels and plot
        df %>%
          mutate(label = glue::glue("{GIA} (n = {n})")) %>%
          # Order facets by GIA factor levels
          mutate(label = factor(label, levels = paste0(levels(GIA), " (n = ", df$n[match(levels(GIA), df$GIA)], ")"))) %>%
          select(GIA, label, lpa_pred_nm, lpa_nmol_max) %>%
          pivot_longer(
            cols = c(lpa_pred_nm, lpa_nmol_max),
            names_to = "measurement_type",
            values_to = "concentration"
          ) %>%
          mutate(
            measurement_type = case_when(
              str_detect(measurement_type, "pred") ~ "Predicted",
              TRUE ~ "Measured"
            )
          ) %>%
          ggplot(aes(x = concentration, fill = measurement_type)) +
          geom_density(alpha = 0.5, adjust = 2) +
          facet_grid(~label) +
          ggsci::scale_fill_jama(name = "") +
          labs(x = "Lp(a) (nmol/L)", y = "Density") +
          theme_bw(base_size = 12) +
          theme(legend.position = "top")
      }

      # Create overall plot
      p1 <- lpa_validation_df %>%
        mutate(GIA = "Overall") %>%
        add_count(GIA, name = "n") %>%
        plot_density(colors = plot_colors, ancestry_order = "Overall")

      # Create plot by ancestry
      p2 <- lpa_validation_df %>%
        plot_density(colors = plot_colors, ancestry_order = ancestry_order$GIA)

      # Return both components and the combined plot
      list(
        overall = p1,
        by_ancestry = p2,
        combined = (p1 / p2) +
          patchwork::plot_layout(heights = c(2, 1), guides = "collect") &
          theme(legend.position = "bottom")
      )
    },
    description = "Generate distribution plot of predicted and measured LPA values"
  ),

  ## Save Distribution Plots ------------------------------------------------
  tar_file(
    plot_distribution_overall_qs,
    {
      file_path <- "Results/distribution_plot_overall.qs"
      qs::qsave(plot_distribution$overall, file_path)
      file_path
    },
    description = "Save overall distribution plot as qs file for programmatic manipulation"
  ),
  tar_file(
    plot_distribution_by_ancestry_qs,
    {
      file_path <- "Results/distribution_plot_by_ancestry.qs"
      qs::qsave(plot_distribution$by_ancestry, file_path)
      file_path
    },
    description = "Save by-ancestry distribution plot as qs file for programmatic manipulation"
  ),
  tar_file(
    plot_distribution_pdf,
    {
      file_path <- "Results/distribution_plot.pdf"
      ggsave(file_path, plot_distribution$combined, width = 12, height = 8)
      file_path
    },
    description = "Save distribution plot as PDF"
  ),
  tar_file(
    plot_distribution_qs,
    {
      file_path <- "Results/distribution_plot.qs"
      qs::qsave(plot_distribution$combined, file_path)
      file_path
    },
    description = "Save distribution plot as qs file for programmatic manipulation"
  ),

  ## Comparison Plot --------------------------------------------------------
  tar_target(
    plot_comparison,
    {
      # Get the subset of colors for the ancestries in our data
      ancestry_colors <- standard_ancestry_palette[names(standard_ancestry_palette) %in% ancestry_order$GIA]

      # Create the plot
      lpa_validation_df %>%
        # Ensure GIA is an ordered factor with correct ordering
        mutate(GIA = factor(GIA, levels = ancestry_order$GIA)) %>%
        mutate(label = glue::glue("{GIA} (n = {n})")) %>%
        ggplot(aes(GIA, lpa_nmol_max)) +
        geom_violin(aes(fill = GIA)) +
        geom_boxplot(width = 0.1) +
        scale_y_log10() +
        scale_fill_manual(values = ancestry_colors) +
        guides(fill = "none") +
        labs(y = "Measured Lp(a) (nmol/L)", x = "") +
        theme_bw(base_size = 12)
    },
    description = "Generate comparison plot showing LPA values by genetic ancestry"
  ),

  ## Save Comparison Plot ---------------------------------------------------
  tar_file(
    plot_comparison_pdf,
    {
      file_path <- "Results/comparison_plot.pdf"
      ggsave(file_path, plot_comparison, width = 10, height = 6)
      file_path
    },
    description = "Save comparison plot as PDF"
  ),
  tar_file(
    plot_comparison_qs,
    {
      file_path <- "Results/comparison_plot.qs"
      qs::qsave(plot_comparison, file_path)
      file_path
    },
    description = "Save comparison plot as qs file for programmatic manipulation"
  ),

  # Bootstrap Configuration ------------------------------------------------
  tar_target(
    bootstrap_config,
    {
      list(
        do_bootstrap = TRUE, # Set to FALSE for quicker, less accurate analysis
        resamples = 1000, # Number of bootstrap resamples
        thresholds = c(125, 150, 175, 200) # LPA thresholds in nmol/L
      )
    },
    description = "Configure bootstrap parameters for performance evaluation"
  ),

  ## Metrics Setup ----------------------------------------------------------
  tar_target(
    numeric_metrics,
    yardstick::metric_set(yardstick::rsq, yardstick::mape, yardstick::rmse),
    description = "Define numeric metrics for performance evaluation (R², MAPE, RMSE)"
  ),
  tar_target(
    class_metrics,
    yardstick::metric_set(yardstick::sens, yardstick::spec, yardstick::ppv, yardstick::npv, yardstick::accuracy, yardstick::f_meas),
    description = "Define classification metrics for performance evaluation (sensitivity, specificity, etc.)"
  ),

  ## Numeric Metrics --------------------------------------------------------
  tar_target(
    lpa_numeric_boot,
    {
      if (bootstrap_config$do_bootstrap) {
        rsample::bootstraps(lpa_validation_df, strata = "GIA", times = bootstrap_config$resamples)
      } else {
        # Placeholder for simplified analysis
        NULL
      }
    },
    description = "Create bootstrap samples for numeric metrics evaluation"
  ),
  tar_target(
    lpa_numeric_metrics_grouped_ci,
    {
      if (bootstrap_config$do_bootstrap) {
        # Full bootstrap analysis
        lpa_numeric_boot %>%
          mutate(metrics = map(splits, ~ {
            analysis(.x) %>%
              group_by(GIA) %>%
              numeric_metrics(truth = lpa_nmol_max, estimate = lpa_pred_nm)
          })) %>%
          unnest(metrics) %>%
          group_by(GIA, .metric) %>%
          summarize(
            estimate = mean(.estimate),
            sd = sd(.estimate),
            lower_ci = quantile(.estimate, 0.025, na.rm = TRUE),
            upper_ci = quantile(.estimate, 0.975, na.rm = TRUE),
            .groups = "drop"
          )
      } else {
        # Simple analysis without bootstrapping
        lpa_validation_df %>%
          group_by(GIA) %>%
          summarize(
            rsq = cor(lpa_nmol_max, lpa_pred_nm)^2,
            rmse = sqrt(mean((lpa_nmol_max - lpa_pred_nm)^2)),
            mape = mean(abs((lpa_nmol_max - lpa_pred_nm) / lpa_nmol_max)) * 100,
            .groups = "drop"
          ) %>%
          pivot_longer(
            cols = c(rsq, rmse, mape),
            names_to = ".metric",
            values_to = "estimate"
          ) %>%
          mutate(
            lower_ci = NA_real_,
            upper_ci = NA_real_
          )
      }
    },
    description = "Calculate numeric metrics by ancestry group with bootstrap confidence intervals"
  ),
  tar_file(
    lpa_numeric_metrics_grouped_ci_csv,
    {
      file_path <- "Results/numeric_metrics_by_group.csv"
      write.csv(lpa_numeric_metrics_grouped_ci, file_path, row.names = FALSE)
      file_path
    },
    description = "Save numeric metrics by group to CSV file"
  ),
  tar_target(
    lpa_numeric_metrics_overall_ci,
    {
      if (bootstrap_config$do_bootstrap) {
        # Full bootstrap analysis
        metrics <- lpa_numeric_boot %>%
          mutate(metrics = map(splits, ~ {
            analysis(.x) %>%
              numeric_metrics(truth = lpa_nmol_max, estimate = lpa_pred_nm)
          })) %>%
          unnest(metrics) %>%
          group_by(.metric) %>%
          summarize(
            estimate = mean(.estimate),
            sd = sd(.estimate),
            lower_ci = quantile(.estimate, 0.025, na.rm = TRUE),
            upper_ci = quantile(.estimate, 0.975, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          mutate(GIA = "Overall")
      } else {
        # Simple analysis without bootstrapping
        metrics <- lpa_validation_df %>%
          summarize(
            rsq = cor(lpa_nmol_max, lpa_pred_nm)^2,
            rmse = sqrt(mean((lpa_nmol_max - lpa_pred_nm)^2)),
            mape = mean(abs((lpa_nmol_max - lpa_pred_nm) / lpa_nmol_max)) * 100
          ) %>%
          pivot_longer(
            cols = everything(),
            names_to = ".metric",
            values_to = "estimate"
          ) %>%
          mutate(
            lower_ci = NA_real_,
            upper_ci = NA_real_,
            GIA = "Overall"
          )
      }

      metrics
    },
    description = "Calculate overall numeric metrics with bootstrap confidence intervals"
  ),
  tar_file(
    lpa_numeric_metrics_overall_ci_csv,
    {
      file_path <- "Results/numeric_metrics_overall.csv"
      write.csv(lpa_numeric_metrics_overall_ci, file_path, row.names = FALSE)
      file_path
    },
    description = "Save overall numeric metrics to CSV file"
  ),


  ## Classification Metrics -------------------------------------------------
  tar_target(
    lpa_class_df,
    {
      tidyr::expand_grid(
        lpa_validation_df %>%
          add_count(GIA, name = "pop_total"),
        tibble(threshold = bootstrap_config$thresholds)
      ) %>%
        mutate(
          true_class = factor(lpa_nmol_max > threshold, levels = c(TRUE, FALSE)),
          pred_class = factor(lpa_pred_nm > threshold, levels = c(TRUE, FALSE))
        ) %>%
        mutate(strata = glue::glue("{GIA}_{true_class}_{threshold}"))
    },
    description = "Create classification dataset for multiple LPA thresholds"
  ),
  tar_target(
    lpa_class_boot,
    {
      if (bootstrap_config$do_bootstrap) {
        rsample::bootstraps(lpa_class_df, strata = "strata", times = bootstrap_config$resamples)
      } else {
        NULL
      }
    },
    description = "Create bootstrap samples for classification metrics evaluation"
  ),
  tar_target(
    lpa_class_metrics_grouped_ci,
    {
      if (bootstrap_config$do_bootstrap) {
        # Full bootstrap analysis
        lpa_class_boot %>%
          mutate(metrics = map(splits, ~ {
            analysis(.x) %>%
              group_by(GIA, threshold) %>%
              class_metrics(truth = true_class, estimate = pred_class)
          })) %>%
          unnest(metrics) %>%
          group_by(GIA, threshold, .metric) %>%
          summarize(
            estimate = mean(.estimate),
            sd = sd(.estimate),
            lower_ci = quantile(.estimate, 0.025, na.rm = TRUE),
            upper_ci = quantile(.estimate, 0.975, na.rm = TRUE),
            .groups = "drop"
          )
      } else {
        # Simple analysis without bootstrapping
        lpa_class_df %>%
          group_by(GIA, threshold) %>%
          summarize(
            metrics = list(class_metrics(truth = true_class, estimate = pred_class)),
            .groups = "drop"
          ) %>%
          unnest(metrics) %>%
          mutate(
            lower_ci = NA_real_,
            upper_ci = NA_real_
          )
      }
    },
    description = "Calculate classification metrics by ancestry group with bootstrap confidence intervals"
  ),
  tar_file(
    lpa_class_metrics_grouped_ci_csv,
    {
      file_path <- "Results/class_metrics_by_group.csv"
      write.csv(lpa_class_metrics_grouped_ci, file_path, row.names = FALSE)
      file_path
    },
    description = "Save classification metrics by group to CSV file"
  ),
  tar_target(
    lpa_class_metrics_overall_ci,
    {
      if (bootstrap_config$do_bootstrap) {
        # Full bootstrap analysis
        lpa_class_boot %>%
          mutate(metrics = map(splits, ~ {
            analysis(.x) %>%
              group_by(threshold) %>%
              class_metrics(truth = true_class, estimate = pred_class)
          })) %>%
          unnest(metrics) %>%
          group_by(threshold, .metric) %>%
          summarize(
            estimate = mean(.estimate),
            sd = sd(.estimate),
            lower_ci = quantile(.estimate, 0.025, na.rm = TRUE),
            upper_ci = quantile(.estimate, 0.975, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          mutate(GIA = "Overall")
      } else {
        # Simple analysis without bootstrapping
        lpa_class_df %>%
          group_by(threshold) %>%
          summarize(
            metrics = list(class_metrics(truth = true_class, estimate = pred_class)),
            .groups = "drop"
          ) %>%
          unnest(metrics) %>%
          mutate(
            lower_ci = NA_real_,
            upper_ci = NA_real_,
            GIA = "Overall"
          )
      }
    },
    description = "Calculate overall classification metrics with bootstrap confidence intervals"
  ),
  tar_file(
    lpa_class_metrics_overall_ci_csv,
    {
      file_path <- "Results/class_metrics_overall.csv"
      write.csv(lpa_class_metrics_overall_ci, file_path, row.names = FALSE)
      file_path
    },
    description = "Save overall classification metrics to CSV file"
  ),


  ## Combined Metrics -------------------------------------------------------
  tar_target(
    lpa_boot_all_metrics,
    {
      bind_rows(
        bind_rows(
          lpa_numeric_metrics_grouped_ci,
          lpa_numeric_metrics_overall_ci
        ) %>% mutate(metric_type = "Numeric"),
        bind_rows(
          lpa_class_metrics_grouped_ci,
          lpa_class_metrics_overall_ci
        ) %>% mutate(metric_type = "Classification")
      )
    },
    description = "Combine all performance metrics for reporting and visualization"
  ),
  tar_file(
    lpa_boot_all_metrics_csv,
    {
      file_path <- "Results/all_metrics.csv"
      write.csv(lpa_boot_all_metrics, file_path, row.names = FALSE)
      file_path
    },
    description = "Save combined metrics to CSV file"
  ),


  ## Group Size Calculations ------------------------------------------------
  tar_target(
    lpa_group_sizes,
    {
      bind_rows(
        lpa_validation_df %>%
          group_by(GIA) %>%
          summarize(pop_size = n()),
        tibble(
          GIA = "Overall",
          pop_size = nrow(lpa_validation_df)
        )
      )
    },
    description = "Calculate sample sizes for each genetic ancestry group"
  ),

  ## Performance Plots ------------------------------------------------------
  tar_target(
    lpa_performance_plot,
    {
      # Get the subset of colors for the ancestries in our data
      data_ancestries <- c("Overall", ancestry_order$GIA)
      ancestry_colors <- standard_ancestry_palette[names(standard_ancestry_palette) %in% ancestry_order$GIA]

      # Add "Overall" color (black)
      plot_colors <- c("Overall" = "black", ancestry_colors)

      # Adjust metrics to add overall and ensure proper ordering
      metrics_adjusted <- lpa_boot_all_metrics %>%
        mutate(
          GIA = factor(GIA, levels = c("Overall", ancestry_order$GIA))
        )

      # Numeric metrics plot
      p1 <- metrics_adjusted %>%
        filter(metric_type == "Numeric") %>%
        left_join(lpa_group_sizes) %>%
        ggplot(aes(
          x = estimate, xmin = lower_ci, xmax = upper_ci,
          y = GIA, color = GIA
        )) +
        geom_pointrange(aes(size = pop_size), shape = 15) +
        facet_wrap(~.metric, scales = "free") +
        scale_x_continuous(
          breaks = scales::extended_breaks(),
          labels = scales::label_number()
        ) +
        scale_size_continuous(range = c(0.25, 1.25)) +
        scale_color_manual(values = plot_colors) +
        guides(color = "none", size = "none") +
        labs(x = "Estimate (95% CI)", y = "") +
        theme_bw(base_size = 12)

      # Classification metrics plot (for threshold = 150)
      p2 <- metrics_adjusted %>%
        filter(
          metric_type == "Classification",
          threshold == 150
        ) %>%
        left_join(lpa_group_sizes) %>%
        ggplot(aes(
          x = estimate, xmin = lower_ci, xmax = upper_ci,
          y = GIA, color = GIA
        )) +
        geom_pointrange(aes(size = pop_size), shape = 15) +
        facet_wrap(~.metric, scales = "free_x") +
        scale_x_continuous(
          breaks = scales::extended_breaks(),
          labels = scales::label_number()
        ) +
        scale_size_continuous(range = c(0.25, 1.25)) +
        scale_color_manual(values = plot_colors) +
        guides(color = "none", size = "none") +
        labs(x = "Estimate (95% CI)", y = "") +
        theme_bw(base_size = 12)

      # Return both components and the combined plot
      list(
        numeric = p1,
        classification = p2,
        combined = (p1 / p2) +
          patchwork::plot_layout(heights = c(1, 3)) +
          patchwork::plot_annotation(tag_levels = "A") &
          theme(plot.tag = element_text(size = 20, face = "bold"))
      )
    },
    description = "Generate performance plots showing numeric and classification metrics"
  ),

  ## Save Performance Plots -------------------------------------------------
  tar_file(
    lpa_performance_plot_numeric_qs,
    {
      file_path <- "Results/performance_plot_numeric.qs"
      qs::qsave(lpa_performance_plot$numeric, file_path)
      file_path
    },
    description = "Save numeric performance plot as qs file for programmatic manipulation"
  ),
  tar_file(
    lpa_performance_plot_classification_qs,
    {
      file_path <- "Results/performance_plot_classification.qs"
      qs::qsave(lpa_performance_plot$classification, file_path)
      file_path
    },
    description = "Save classification performance plot as qs file for programmatic manipulation"
  ),
  tar_file(
    lpa_performance_plot_pdf,
    {
      file_path <- "Results/performance_plot.pdf"
      ggsave(file_path, lpa_performance_plot$combined, width = 12, height = 10)
      file_path
    },
    description = "Save performance plot as PDF"
  ),
  tar_file(
    lpa_performance_plot_qs,
    {
      file_path <- "Results/performance_plot.qs"
      qs::qsave(lpa_performance_plot$combined, file_path)
      file_path
    },
    description = "Save performance plot as qs file for programmatic manipulation"
  ),

  ## NNT Calculations -------------------------------------------------------
  tar_target(
    lpa_nnt,
    {
      lpa_boot_all_metrics %>%
        filter(threshold == 150) %>%
        filter(.metric == "ppv") %>%
        mutate(
          nnt_estimate = 1 / estimate,
          nnt_lower = 1 / upper_ci, # CI bounds flip when taking reciprocal
          nnt_upper = 1 / lower_ci
        ) %>%
        select(GIA, nnt_estimate, nnt_lower, nnt_upper) %>%
        mutate(
          .metric = "nnt",
          estimate = nnt_estimate,
          lower_ci = nnt_lower,
          upper_ci = nnt_upper,
          metric_type = "Classification"
        ) %>%
        select(GIA, .metric, estimate, lower_ci, upper_ci, metric_type)
    },
    description = "Calculate Number Needed to Test (NNT) from positive predictive value"
  ),
  tar_file(
    lpa_nnt_csv,
    {
      file_path <- "Results/nnt_metrics.csv"
      write.csv(lpa_nnt, file_path, row.names = FALSE)
      file_path
    },
    description = "Save NNT metrics to CSV file"
  ),


  ## NNT Plot ---------------------------------------------------------------
  tar_target(
    lpa_nnt_plot,
    {
      # Get the subset of colors for the ancestries in our data
      data_ancestries <- c("Overall", ancestry_order$GIA)
      ancestry_colors <- standard_ancestry_palette[names(standard_ancestry_palette) %in% ancestry_order$GIA]

      # Add "Overall" color (black)
      plot_colors <- c("Overall" = "black", ancestry_colors)

      # Make sure overall comes first and others follow in order by sample size
      ordered_nnt <- lpa_nnt %>%
        mutate(
          GIA = factor(GIA, levels = c("Overall", ancestry_order$GIA))
        )

      ordered_nnt %>%
        left_join(lpa_group_sizes) %>%
        ggplot(aes(
          x = estimate, xmin = lower_ci, xmax = upper_ci,
          y = GIA, color = GIA
        )) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_pointrange(aes(size = pop_size), shape = 15) +
        scale_color_manual(values = plot_colors) +
        scale_x_continuous(
          breaks = scales::extended_breaks(),
          labels = scales::label_number()
        ) +
        scale_size_continuous(range = c(0.25, 1.25)) +
        guides(color = "none", size = "none") +
        labs(
          x = "Number Needed to Test (95% CI)",
          y = "",
          title = "Number Needed to Test to Find One High Lp(a) Individual (>150 nmol/L)"
        ) +
        theme_bw(base_size = 12)
    },
    description = "Generate NNT plot showing Number Needed to Test across genetic ancestry groups"
  ),

  ## Save NNT Plot ----------------------------------------------------------
  tar_file(
    lpa_nnt_plot_pdf,
    {
      file_path <- "Results/nnt_plot.pdf"
      ggsave(file_path, lpa_nnt_plot, width = 10, height = 6)
      file_path
    },
    description = "Save NNT plot as PDF"
  ),
  tar_file(
    lpa_nnt_plot_qs,
    {
      file_path <- "Results/nnt_plot.qs"
      qs::qsave(lpa_nnt_plot, file_path)
      file_path
    },
    description = "Save NNT plot as qs file for programmatic manipulation"
  ),

  # Prevalence Analysis ----------------------------------------------------
  tar_target(
    lpa_prevalence_table,
    {
      # First prepare the data for each threshold
      thresholds <- bootstrap_config$thresholds

      # Use purrr::map_dfr to iterate over thresholds and bind results
      purrr::map_dfr(thresholds, function(current_threshold) {
        # Overall rates first - using ALL genotyped individuals as denominator
        overall_df <- lpa_pred_measured_df %>%
          mutate(total_count = n()) %>%
          summarize(
            total_count = unique(total_count),
            # Measured LPA (clinical) - using total genotyped as denominator
            measured_available_count = sum(!is.na(lpa_nmol_max)),
            measured_elevated_count = sum(lpa_nmol_max > current_threshold, na.rm = TRUE),
            measured_elevated_pct = measured_elevated_count / total_count * 100,

            # Predicted LPA (genetic) - all genotyped individuals
            predicted_available_count = sum(!is.na(lpa_pred_nm)),
            predicted_elevated_count = sum(lpa_pred_nm > current_threshold, na.rm = TRUE),
            predicted_elevated_pct = predicted_elevated_count / total_count * 100
          ) %>%
          # Join with PPV from validation statistics
          left_join(
            lpa_class_metrics_overall_ci %>%
              filter(.metric == "ppv", threshold == current_threshold, GIA == "Overall") %>%
              select(ppv = estimate),
            by = character()
          ) %>%
          # Adjust prediction by PPV when available
          mutate(
            # Calculate both unadjusted and adjusted values
            predicted_elevated_unadjusted = predicted_elevated_count,
            predicted_elevated_unadjusted_pct = predicted_elevated_count / total_count * 100,

            # Adjust by PPV only if PPV is available
            predicted_elevated_adjusted = if_else(!is.na(ppv),
                                                  predicted_elevated_count * ppv,
                                                  NA_real_
            ),
            predicted_elevated_adjusted_pct = if_else(!is.na(ppv),
                                                      predicted_elevated_adjusted / total_count * 100,
                                                      NA_real_
            ),
            GIA = "Overall",
            threshold = current_threshold,
            has_ppv_adjustment = !is.na(ppv)
          )

        # Then by ancestry group
        ancestry_df <- lpa_pred_measured_df %>%
          # Include only individuals with known genetic ancestry
          filter(!is.na(GIA)) %>%
          add_count(GIA, name = "total_count") %>%
          group_by(GIA) %>%
          summarize(
            # Total genotyped individuals in each ancestry group
            total_count = unique(total_count),

            # Measured LPA (clinical) - using total genotyped in group as denominator
            measured_available_count = sum(!is.na(lpa_nmol_max)),
            measured_elevated_count = sum(lpa_nmol_max > current_threshold, na.rm = TRUE),
            measured_elevated_pct = measured_elevated_count / total_count * 100,

            # Predicted LPA (genetic)
            predicted_available_count = sum(!is.na(lpa_pred_nm)),
            predicted_elevated_count = sum(lpa_pred_nm > current_threshold, na.rm = TRUE),
            predicted_elevated_pct = predicted_elevated_count / total_count * 100
          ) %>%
          # Join with PPV values for each ancestry group
          left_join(
            lpa_class_metrics_grouped_ci %>%
              filter(.metric == "ppv", threshold == current_threshold) %>%
              select(GIA, ppv = estimate),
            by = "GIA"
          ) %>%
          # Adjust prediction by PPV when available
          mutate(
            # Calculate both unadjusted and adjusted values
            predicted_elevated_unadjusted = predicted_elevated_count,
            predicted_elevated_unadjusted_pct = predicted_elevated_count / total_count * 100,

            # Adjust by PPV only if PPV is available
            predicted_elevated_adjusted = if_else(!is.na(ppv),
                                                  predicted_elevated_count * ppv,
                                                  NA_real_
            ),
            predicted_elevated_adjusted_pct = if_else(!is.na(ppv),
                                                      predicted_elevated_adjusted / total_count * 100,
                                                      NA_real_
            ),
            threshold = current_threshold,
            has_ppv_adjustment = !is.na(ppv)
          )

        # Combine overall and ancestry-specific results
        bind_rows(overall_df, ancestry_df)
      })
    },
    description = "Calculate LPA prevalence using all genotyped individuals as denominator, with PPV adjustment when available"
  ),

  ## Save Prevalence Table --------------------------------------------------
  tar_file(
    lpa_prevalence_csv,
    {
      file_path <- "Results/lpa_prevalence_by_threshold.csv"
      write.csv(lpa_prevalence_table, file_path, row.names = FALSE)
      file_path
    },
    description = "Save LPA prevalence table to CSV file"
  ),

  # Detection Improvement Analysis -----------------------------------------
  tar_target(
    lpa_detection_improvement,
    lpa_prevalence_table %>%
      # Calculate improvement ratio with confidence intervals
      mutate(
        # Clinical detection rate (per 1000)
        clinical_rate_per_1000 = measured_elevated_count / total_count * 1000,

        # Genetic detection rates (per 1000)
        genetic_rate_unadjusted_per_1000 = predicted_elevated_unadjusted / total_count * 1000,
        genetic_rate_per_1000 = if_else(!is.na(predicted_elevated_adjusted),
                                        predicted_elevated_adjusted / total_count * 1000,
                                        NA_real_
        ),

        # Improvement factors
        improvement_factor_unadjusted = genetic_rate_unadjusted_per_1000 / clinical_rate_per_1000,
        improvement_factor = if_else(!is.na(genetic_rate_per_1000),
                                     genetic_rate_per_1000 / clinical_rate_per_1000,
                                     NA_real_
        ),

        # Fisher's exact test for significance - use appropriate count based on PPV availability
        p_value = purrr::pmap_dbl(
          list(
            predicted_count = if_else(has_ppv_adjustment,
                                      as.numeric(predicted_elevated_adjusted),
                                      as.numeric(predicted_elevated_unadjusted)
            ),
            measured_count = measured_elevated_count,
            n_total = total_count
          ),
          function(predicted_count, measured_count, n_total) {
            # Skip if either count is zero or NA
            if (is.na(predicted_count) || is.na(measured_count) ||
                predicted_count == 0 || measured_count == 0) {
              return(NA_real_)
            }

            # Create contingency table (rounded counts for Fisher's test)
            predicted_count <- round(predicted_count)
            matrix <- matrix(
              c(
                predicted_count, n_total - predicted_count,
                measured_count, n_total - measured_count
              ),
              nrow = 2, byrow = TRUE
            )

            # Run Fisher's exact test
            test <- fisher.test(matrix)
            return(test$p.value)
          }
        ),

        # Add a note about which rate was used
        rate_note = if_else(has_ppv_adjustment,
                            "PPV adjusted",
                            "Unadjusted (no PPV available)"
        )
      ),
    description = "Calculate detection improvement comparing genetic prediction to clinical measurements"
  ),

  ## Comparison Table -------------------------------------------------------
  tar_target(
    lpa_detection_comparison_150,
    lpa_detection_improvement %>%
      filter(threshold == 150) %>%
      # Choose appropriate columns based on PPV availability
      mutate(
        # Use adjusted values when available, otherwise use unadjusted
        genetic_count = if_else(!is.na(predicted_elevated_adjusted),
                                predicted_elevated_adjusted,
                                predicted_elevated_unadjusted
        ),
        genetic_pct = if_else(!is.na(predicted_elevated_adjusted_pct),
                              predicted_elevated_adjusted_pct,
                              predicted_elevated_unadjusted_pct
        ),
        improvement = if_else(!is.na(improvement_factor),
                              improvement_factor,
                              improvement_factor_unadjusted
        )
      ) %>%
      select(
        GIA,
        total_count,
        clinical_count = measured_elevated_count,
        clinical_pct = measured_elevated_pct,
        genetic_count,
        genetic_pct,
        improvement,
        p_value,
        rate_note
      ) %>%
      # Order rows with Overall first, then by sample size
      mutate(
        GIA = factor(GIA, levels = c("Overall", ancestry_order$GIA))
      ) %>%
      arrange(GIA) %>%
      # Format percentages and p-values nicely
      mutate(
        clinical_pct = sprintf("%.1f%%", clinical_pct),
        genetic_pct = sprintf("%.1f%%", genetic_pct),
        improvement = case_when(
          is.na(improvement) ~ "NA",
          is.infinite(improvement) ~ "∞",
          TRUE ~ sprintf("%.2fx", improvement)
        ),
        p_value = case_when(
          is.na(p_value) ~ "NA",
          p_value < 0.001 ~ "p < 0.001",
          p_value < 0.01 ~ "p < 0.01",
          p_value < 0.05 ~ "p < 0.05",
          TRUE ~ "p ≥ 0.05"
        )
      )
  ),
  tar_file(
    lpa_detection_comparison_150_csv,
    {
      file_path <- "Results/lpa_detection_comparison_150.csv"
      write.csv(lpa_detection_comparison_150, file_path, row.names = FALSE)
      file_path
    },
    description = "Save LPA detection comparison table at 150 nmol/L to CSV file"
  ),

  # Report Generation ------------------------------------------------------
  tar_render(
    name = lpa_validation_report,
    path = "rmarkdown/lpa_validation_report.Rmd"
  ),

  # Create ZIP Archive -----------------------------------------------------
  tar_file(
    zip_results,
    {
      # Zip the Results directory and report
      zip_file <- "lpa_validation_results.zip"

      # Ensure the report is included in the zip
      report_file <- lpa_validation_report
      file.copy(report_file[1], "Results/lpa_validation_report.html", overwrite = TRUE)

      # Check if zip is available; if not, use alternative method
      if (requireNamespace("zip", quietly = TRUE)) {
        zip::zip(
          zipfile = zip_file,
          files = list.files("Results", full.names = TRUE),
          recurse = TRUE,
          compression_level = 9
        )
      } else {
        utils::zip(
          zipfile = zip_file,
          files = list.files("Results", full.names = TRUE)
        )
      }

      zip_file
    },
    description = "Create ZIP archive of all results"
  )
)
