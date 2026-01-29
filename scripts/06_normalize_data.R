#!/usr/bin/env Rscript
# ==============================================================================
# 06_normalize_data.R - Proteomics Data Normalisation
# ==============================================================================
#
# Purpose:
#   Normalises proteomics data to remove technical variability and ensure sample
#   comparability within a single batch. This step performs EXCLUSIVELY within-batch
#   median normalisation, regardless of single-batch or multi-batch mode.
#   - **Within-Batch Normalisation**: Median normalisation (standard intra-batch step)
#   - **Rationale**: Median normalisation ensures samples are comparable within a batch
#     before performing statistical tests. This is a standard preprocessing step for
#     proteomics data.
#   - **Note**: Cross-batch normalisation (bridge, ComBat) is handled in step 07
#     (bridge_normalization.R). This step does NOT perform cross-batch normalisation.
#   - **Expected Performance**: Typically achieves ~9.7% SD reduction
#
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Date: December 2025
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(OlinkAnalyze)

  library(yaml)
  library(logger)
  library(ggplot2)
})

# Source path utilities for batch-aware paths
# Get script directory safely (handles both direct execution and sourcing)
script_dir <- tryCatch({
  env_script <- Sys.getenv("SCRIPT_NAME", "")
  if (env_script != "" && file.exists(env_script)) {
    dirname(normalizePath(env_script))
  } else {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      script_path <- sub("^--file=", "", file_arg)
      dirname(normalizePath(script_path))
    } else {
      getwd()
    }
  }
}, error = function(e) getwd())
source(file.path(script_dir, "path_utils.R"))

# Get batch context
batch_id <- Sys.getenv("PIPELINE_BATCH_ID", "batch_01")
step_num <- get_step_number()

# Load configuration
config_file <- Sys.getenv("PIPELINE_CONFIG", "")
if (config_file == "" || !file.exists(config_file)) {
  stop("PIPELINE_CONFIG environment variable not set or config file not found. Please provide path to config file.")
}
config <- read_yaml(config_file)

# Set up logging with batch-aware path
log_path <- get_log_path(step_num, batch_id, config = config)
log_appender(appender_file(log_path))
log_info("Starting data normalisation for batch: {batch_id}")

# Function for median normalisation (within-batch)
normalize_median <- function(npx_matrix, by_plate = FALSE, plate_info = NULL) {

  if(by_plate && !is.null(plate_info)) {
    log_info("=== BY-PLATE MEDIAN NORMALIZATION ===")
    log_info("Number of plates: {length(unique(plate_info$PlateID))}")

    # Group samples by plate
    normalized_list <- list()
    scaling_factors <- NULL  # Initialize for return

    for(plate in unique(plate_info$PlateID)) {
      plate_samples <- plate_info[PlateID == plate]$SampleID
      plate_samples <- intersect(plate_samples, rownames(npx_matrix))

      if(length(plate_samples) > 0) {
        log_info("  Plate {plate}: {length(plate_samples)} samples")
        plate_matrix <- npx_matrix[plate_samples, , drop = FALSE]

        # Median normalise within plate
        col_medians <- apply(plate_matrix, 2, median, na.rm = TRUE)
        global_median <- median(col_medians, na.rm = TRUE)
        plate_scaling_factors <- global_median / col_medians
        plate_scaling_factors[is.na(plate_scaling_factors) | is.infinite(plate_scaling_factors)] <- 1

        # Cap extreme values
        plate_scaling_factors[plate_scaling_factors < 0.1] <- 0.1
        plate_scaling_factors[plate_scaling_factors > 10] <- 10

        normalized_list[[plate]] <- sweep(plate_matrix, 2, plate_scaling_factors, "*")
      }
    }

    # Combine normalised plates
    normalized_matrix <- do.call(rbind, normalized_list)
    log_info("Combined normalized matrix: {nrow(normalized_matrix)} × {ncol(normalized_matrix)}")

  } else {
    log_info("=== GLOBAL MEDIAN NORMALIZATION ===")

    # Calculate column medians
    col_medians <- apply(npx_matrix, 2, median, na.rm = TRUE)
    log_info("Protein median statistics:")
    log_info("  Min: {round(min(col_medians, na.rm=TRUE), 3)}")
    log_info("  Max: {round(max(col_medians, na.rm=TRUE), 3)}")
    log_info("  Median: {round(median(col_medians, na.rm=TRUE), 3)}")

    # Calculate global median
    global_median <- median(col_medians, na.rm = TRUE)
    log_info("Global median (target): {round(global_median, 3)}")

    # Calculate scaling factors
    scaling_factors <- global_median / col_medians
    scaling_factors[is.na(scaling_factors) | is.infinite(scaling_factors)] <- 1

    # Cap extreme scaling factors to prevent outliers (keep within 0.1 to 10 range)
    n_capped_low <- sum(scaling_factors < 0.1)
    n_capped_high <- sum(scaling_factors > 10)
    scaling_factors[scaling_factors < 0.1] <- 0.1
    scaling_factors[scaling_factors > 10] <- 10

    log_info("Scaling factor statistics:")
    log_info("  Min: {round(min(scaling_factors), 3)}")
    log_info("  Max: {round(max(scaling_factors), 3)}")
    log_info("  Median: {round(median(scaling_factors), 3)}")
    if (n_capped_low > 0 || n_capped_high > 0) {
      log_info("  Capped: {n_capped_low} low (<0.1), {n_capped_high} high (>10)")
    }

    # Apply scaling
    log_info("Applying scaling factors to {ncol(npx_matrix)} proteins...")
    normalized_matrix <- sweep(npx_matrix, 2, scaling_factors, "*")
    log_info("Output matrix: {nrow(normalized_matrix)} × {ncol(normalized_matrix)}")
  }

  # Note: NPX data is already log2-transformed by Olink, so we do NOT apply additional log transformation
  log_info("Note: No additional log transformation applied (NPX is already log2-transformed)")

  return(list(
    normalized_matrix = normalized_matrix,
    scaling_factors = if(exists("scaling_factors")) scaling_factors else NULL,
    method = ifelse(by_plate, "median_by_plate", "median_global")
  ))
}

# Master normalization function (within-batch median normalization only)
normalize_data <- function(npx_matrix, method = "median", plate_info = NULL) {

  log_info("╔══════════════════════════════════════════════════════════════════╗")
  log_info("║ WITHIN-BATCH NORMALIZATION (Step 06)                             ║")
  log_info("╚══════════════════════════════════════════════════════════════════╝")
  log_info("Method: {method}")
  log_info("Input matrix: {nrow(npx_matrix)} samples × {ncol(npx_matrix)} proteins")

  result <- switch(method,
    "median" = {
      normalize_median(npx_matrix, by_plate = FALSE)
    },
    "median_plate" = {
      if(is.null(plate_info)) {
        log_warn("Plate info not provided, using global median normalization")
        normalize_median(npx_matrix, by_plate = FALSE)
      } else {
        normalize_median(npx_matrix, by_plate = TRUE, plate_info = plate_info)
      }
    },
    {
      log_error("Unknown normalization method: {method}")
      log_error("Step 06 only supports 'median' and 'median_plate' methods")
      log_error("For cross-batch methods (bridge, ComBat), use step 07")
      stop("Invalid normalization method for step 06")
    }
  )

  log_info("════════════════════════════════════════════════════════════════════")
  return(result)
}

# Function to evaluate normalisation
evaluate_normalization <- function(raw_matrix, normalized_matrix, metadata) {
  log_info("Evaluating normalisation performance")

  # Validate inputs
  if (is.null(normalized_matrix)) {
    log_error("normalized_matrix is NULL - cannot evaluate normalization")
    stop("normalized_matrix is NULL in evaluate_normalization()")
  }

  if (!is.matrix(normalized_matrix) && !is.data.frame(normalized_matrix)) {
    log_error("normalized_matrix is not a matrix or data.frame")
    stop("normalized_matrix must be a matrix or data.frame")
  }

  if (nrow(normalized_matrix) == 0 || ncol(normalized_matrix) == 0) {
    log_error("normalized_matrix has zero dimensions: {nrow(normalized_matrix)} x {ncol(normalized_matrix)}")
    stop("normalized_matrix has invalid dimensions")
  }

  if (nrow(normalized_matrix) != nrow(raw_matrix) || ncol(normalized_matrix) != ncol(raw_matrix)) {
    log_error("Dimension mismatch: raw_matrix {nrow(raw_matrix)}x{ncol(raw_matrix)} vs normalized_matrix {nrow(normalized_matrix)}x{ncol(normalized_matrix)}")
    stop("raw_matrix and normalized_matrix must have the same dimensions")
  }

  # For log-transformed data like NPX, use SD and MAD instead of CV
  # (CV is not meaningful when data is centered around 0)

  # Calculate per-protein SD before and after
  sd_before <- apply(raw_matrix, 2, sd, na.rm = TRUE)
  sd_after <- apply(normalized_matrix, 2, sd, na.rm = TRUE)

  # Calculate per-protein MAD (Median Absolute Deviation) before and after
  mad_before <- apply(raw_matrix, 2, mad, na.rm = TRUE)
  mad_after <- apply(normalized_matrix, 2, mad, na.rm = TRUE)

  # Calculate per-protein IQR before and after
  iqr_before <- apply(raw_matrix, 2, IQR, na.rm = TRUE)
  iqr_after <- apply(normalized_matrix, 2, IQR, na.rm = TRUE)

  # Calculate summary statistics
  mean_sd_before <- mean(sd_before, na.rm = TRUE)
  median_sd_before <- median(sd_before, na.rm = TRUE)
  mean_sd_after <- mean(sd_after, na.rm = TRUE)
  median_sd_after <- median(sd_after, na.rm = TRUE)
  sd_reduction <- mean(sd_before - sd_after, na.rm = TRUE)

  mean_mad_before <- mean(mad_before, na.rm = TRUE)
  median_mad_before <- median(mad_before, na.rm = TRUE)
  mean_mad_after <- mean(mad_after, na.rm = TRUE)
  median_mad_after <- median(mad_after, na.rm = TRUE)
  mad_reduction <- mean(mad_before - mad_after, na.rm = TRUE)

  mean_iqr_before <- mean(iqr_before, na.rm = TRUE)
  median_iqr_before <- median(iqr_before, na.rm = TRUE)
  mean_iqr_after <- mean(iqr_after, na.rm = TRUE)
  median_iqr_after <- median(iqr_after, na.rm = TRUE)
  iqr_reduction <- mean(iqr_before - iqr_after, na.rm = TRUE)

  # Calculate percentage reductions
  # For absolute metrics: (before - after) / before * 100
  # For reduction metrics: reduction / before_value * 100
  pct_reduction_mean_sd <- ifelse(mean_sd_before > 0, (mean_sd_before - mean_sd_after) / mean_sd_before * 100, NA_real_)
  pct_reduction_median_sd <- ifelse(median_sd_before > 0, (median_sd_before - median_sd_after) / median_sd_before * 100, NA_real_)
  pct_reduction_sd <- ifelse(mean_sd_before > 0, sd_reduction / mean_sd_before * 100, NA_real_)

  pct_reduction_mean_mad <- ifelse(mean_mad_before > 0, (mean_mad_before - mean_mad_after) / mean_mad_before * 100, NA_real_)
  pct_reduction_median_mad <- ifelse(median_mad_before > 0, (median_mad_before - median_mad_after) / median_mad_before * 100, NA_real_)
  pct_reduction_mad <- ifelse(mean_mad_before > 0, mad_reduction / mean_mad_before * 100, NA_real_)

  pct_reduction_mean_iqr <- ifelse(mean_iqr_before > 0, (mean_iqr_before - mean_iqr_after) / mean_iqr_before * 100, NA_real_)
  pct_reduction_median_iqr <- ifelse(median_iqr_before > 0, (median_iqr_before - median_iqr_after) / median_iqr_before * 100, NA_real_)
  pct_reduction_iqr <- ifelse(mean_iqr_before > 0, iqr_reduction / mean_iqr_before * 100, NA_real_)

  # Summary statistics
  # Note: For "reduction" metrics (SD reduction, MAD reduction, IQR reduction),
  # the "before" column is set to NA because reduction is a difference metric
  # (before - after), not an absolute value. The "after" column contains the
  # actual reduction value.
  eval_stats <- data.table(
    metric = c("Mean SD", "Median SD", "SD reduction",
               "Mean MAD", "Median MAD", "MAD reduction",
               "Mean IQR", "Median IQR", "IQR reduction"),
    before = c(mean_sd_before,
              median_sd_before,
              NA_real_,  # Reduction is a difference, not an absolute value
              mean_mad_before,
              median_mad_before,
              NA_real_,  # Reduction is a difference, not an absolute value
              mean_iqr_before,
              median_iqr_before,
              NA_real_), # Reduction is a difference, not an absolute value
    after = c(mean_sd_after,
             median_sd_after,
             sd_reduction,  # SD reduction = before - after
             mean_mad_after,
             median_mad_after,
             mad_reduction,  # MAD reduction = before - after
             mean_iqr_after,
             median_iqr_after,
             iqr_reduction),  # IQR reduction = before - after
    pct_reduction = c(pct_reduction_mean_sd,
                     pct_reduction_median_sd,
                     pct_reduction_sd,
                     pct_reduction_mean_mad,
                     pct_reduction_median_mad,
                     pct_reduction_mad,
                     pct_reduction_mean_iqr,
                     pct_reduction_median_iqr,
                     pct_reduction_iqr)
  )

  log_info("Mean SD - Before: {round(eval_stats$before[1], 3)}, After: {round(eval_stats$after[1], 3)}")
  log_info("SD reduction: {round(eval_stats$after[3], 3)}")

  return(eval_stats)
}

# Function to create normalization plots
create_normalization_plots <- function(raw_matrix, normalized_matrix, title_suffix = "", n_samples = 100) {
  log_info("Creating normalisation visualisation plots")

  # Sample distributions before and after
  set.seed(123)
  n_samples_actual <- min(n_samples, nrow(raw_matrix))
  sample_idx <- sample(nrow(raw_matrix), n_samples_actual)

  # Before normalisation - convert to data.table first
  df_before <- as.data.table(raw_matrix[sample_idx, ], keep.rownames = TRUE)
  setnames(df_before, "rn", "Sample")
  df_before <- melt(df_before, id.vars = "Sample", variable.name = "Protein", value.name = "Expression")
  df_before$Stage <- "Before"

  # After normalisation - convert to data.table first
  df_after <- as.data.table(normalized_matrix[sample_idx, ], keep.rownames = TRUE)
  setnames(df_after, "rn", "Sample")
  df_after <- melt(df_after, id.vars = "Sample", variable.name = "Protein", value.name = "Expression")
  df_after$Stage <- "After"

  # Combine
  df_combined <- rbind(df_before, df_after)
  df_combined$Stage <- factor(df_combined$Stage, levels = c("Before", "After"))

  # Calculate summary statistics for before and after
  stats_before <- df_before[, .(
    Mean = round(mean(Expression, na.rm = TRUE), 2),
    Median = round(median(Expression, na.rm = TRUE), 2),
    SD = round(sd(Expression, na.rm = TRUE), 2),
    IQR = round(IQR(Expression, na.rm = TRUE), 2)
  )]

  stats_after <- df_after[, .(
    Mean = round(mean(Expression, na.rm = TRUE), 2),
    Median = round(median(Expression, na.rm = TRUE), 2),
    SD = round(sd(Expression, na.rm = TRUE), 2),
    IQR = round(IQR(Expression, na.rm = TRUE), 2)
  )]

  # Create subtitle with sample info
  total_samples <- nrow(raw_matrix)
  total_proteins <- ncol(raw_matrix)
  subtitle <- sprintf(
    "Showing %d randomly selected samples (out of %d total) × %d proteins\nBefore: Mean=%.2f, SD=%.2f, IQR=%.2f | After: Mean=%.2f, SD=%.2f, IQR=%.2f",
    n_samples_actual, total_samples, total_proteins,
    stats_before$Mean, stats_before$SD, stats_before$IQR,
    stats_after$Mean, stats_after$SD, stats_after$IQR
  )

  # Create boxplot
  p <- ggplot(df_combined, aes(x = as.factor(Sample), y = Expression, fill = Stage)) +
    geom_boxplot(alpha = 0.6, outlier.size = 0.5) +
    facet_wrap(~ Stage, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("Before" = "#FF6B6B", "After" = "#3A5F8A")) +
    labs(title = paste("Normalization Effect", title_suffix),
         subtitle = subtitle,
         x = "Sample (randomly selected, ordered by index)",
         y = "NPX Expression Level (log2 scale)") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.subtitle = element_text(size = 9, color = "gray40")
    )

  return(p)
}

# Main execution
main <- function() {

  # Load data
  log_info("Loading data from previous steps")
  # Use final cleaned matrix (all outliers removed) from step 05
  # Fallback to step 00 QC matrix if final cleaned matrix doesn't exist yet
  final_cleaned_path <- get_output_path("05d", "05d_npx_matrix_all_qc_passed", batch_id, "phenotypes", "rds", config = config)
  if (file.exists(final_cleaned_path)) {
    log_info("Using final cleaned matrix (all outliers removed) from step 05d")
    npx_matrix <- readRDS(final_cleaned_path)
  } else {
    log_warn("Final cleaned matrix not found: {final_cleaned_path}")
    log_warn("Falling back to step 00 analysis-ready matrix (outliers may not be removed)")
    npx_matrix_path <- get_output_path("00", "npx_matrix_analysis_ready", batch_id, "qc", config = config)
    if (!file.exists(npx_matrix_path)) {
      stop("NPX matrix file not found: {npx_matrix_path}")
    }
    npx_matrix <- readRDS(npx_matrix_path)
  }
  # Load metadata (required for evaluation)
  metadata_path <- get_output_path("00", "metadata", batch_id, "qc", config = config)
  if (!file.exists(metadata_path)) {
    stop("Metadata file not found: {metadata_path}")
  }
  metadata <- readRDS(metadata_path)

  # CRITICAL: Step 6 performs EXCLUSIVELY within-batch median normalisation
  # Cross-batch normalisation (bridge, ComBat) is handled in step 07 (bridge_normalization.R)
  # This ensures clear separation: step 6 = within-batch, step 7 = cross-batch
  #
  # Note: bridge_result and samples_data are NOT required for median normalisation
  # They are only needed for bridge/ComBat methods which are handled in step 07

  log_info("=== WITHIN-BATCH NORMALISATION (Step 06) ===")
  log_info("This step performs EXCLUSIVELY within-batch median normalisation")
  log_info("Cross-batch normalisation (bridge, ComBat) is handled in step 07")

  # Apply EXCLUSIVELY within-batch median normalisation
  # This step does NOT perform cross-batch normalisation (bridge, ComBat)
  # Cross-batch normalisation is handled in step 07 (bridge_normalization.R)

  log_info("Applying within-batch median normalisation (standard intra-batch step)")
  log_info("Median normalisation ensures samples are comparable within the batch")

  # Primary normalisation: Median (within-batch only)
  # Note: normalize_data() now only accepts method and plate_info parameters
  # Cross-batch methods (bridge, combat) are handled in step 07
  normalized_result <- normalize_data(npx_matrix, method = "median", plate_info = NULL)

  if (is.null(normalized_result)) {
    log_error("Median normalization failed")
    stop("Normalization failed - cannot proceed")
  }

  # Set method used (always median for step 6)
  normalization_method_used <- "median"

  log_info("Within-batch median normalisation completed")

  # Evaluate within-batch median normalisation
  log_info("Evaluating within-batch median normalisation")
  eval_median <- evaluate_normalization(npx_matrix, normalized_result$normalized_matrix, metadata)

  # Combine evaluations (step 6 only evaluates median normalisation)
  all_evaluations <- cbind(method = "median", eval_median)

  # Create plots (step 6 only plots median normalisation)
  # Note: Bridge and ComBat plots are generated in step 07 for cross-batch normalisation
  plot_median <- create_normalization_plots(npx_matrix, normalized_result$normalized_matrix, "- Median (Within-Batch)")

  # Save outputs with batch-aware paths
  log_info("Saving normalisation results")

  # Save primary normalisation with 06_ prefix
  normalized_matrix_path <- get_output_path(step_num, "npx_matrix_normalized", batch_id, "normalized", config = config)
  normalization_result_path <- get_output_path(step_num, "normalization_result", batch_id, "normalized", config = config)
  ensure_output_dir(normalized_matrix_path)
  ensure_output_dir(normalization_result_path)

  saveRDS(normalized_result$normalized_matrix, normalized_matrix_path)
  saveRDS(normalized_result, normalization_result_path)

  log_info("Saved primary normalisation ({normalization_method_used}): {normalized_matrix_path}")

  # Step 6 does NOT save alternative normalisations (bridge, ComBat)
  # These are cross-batch methods handled in step 07
  log_info("Step 6: Only median normalisation saved (within-batch only)")
  log_info("Cross-batch normalisations (bridge, ComBat) are handled in step 07")

  # Save evaluations with 06_ prefix
  evaluations_path <- get_output_path(step_num, "normalization_evaluations", batch_id, "normalized", "tsv", config = config)
  ensure_output_dir(evaluations_path)
  fwrite(all_evaluations, evaluations_path, sep = "\t")
  log_info("Saved normalisation evaluations: {evaluations_path}")

  # Save plots (step 6 only saves median normalisation plot)
  # Bridge and ComBat plots are generated in step 07 for cross-batch normalisation
  median_plot_path <- get_output_path(step_num, "normalization_effect_median", batch_id, "normalized", "pdf", config = config)
  ensure_output_dir(median_plot_path)

  if (!is.null(plot_median)) {
    ggsave(median_plot_path, plot_median, width = 14, height = 9)
    log_info("Saved median normalisation plot: {median_plot_path}")
  } else {
    log_warn("Skipping median normalization plot (normalization failed)")
  }

  # Print summary
  cat("\n=== WITHIN-BATCH NORMALIZATION SUMMARY (Step 06) ===\n")
  cat("Method used: Median (within-batch only)\n")
  cat("Note: Cross-batch normalisation (bridge, ComBat) is handled in step 07\n")
  cat("Matrix dimensions:", nrow(normalized_result$normalized_matrix), "x",
      ncol(normalized_result$normalized_matrix), "\n")

  cat("\n=== NORMALIZATION EVALUATION ===\n")
  print(all_evaluations)

  cat("\nNormalized data saved to: ../output/normalized/\n")

  log_info("Normalisation completed")

  # Step 6 does NOT return alternative normalizations (bridge, ComBat)
  # These are cross-batch methods handled in step 07
  alternative_normalizations <- list()

  return(list(
    normalized_matrix = normalized_result$normalized_matrix,
    normalization_result = normalized_result,
    evaluations = all_evaluations,
    alternative_normalizations = alternative_normalizations,
    mode = "within_batch"
  ))
}

# Run if executed directly
if (!interactive()) {
  result <- main()
}
