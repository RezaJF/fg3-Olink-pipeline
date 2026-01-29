#!/usr/bin/env Rscript
# ==============================================================================
# 07_bridge_normalization.R - Enhanced Bridge Sample Normalisation
# ==============================================================================
#
# Purpose:
#   Enhanced bridge sample normalisation for multi-batch harmonisation. Provides
#   median-based and quantile-based normalisation methods using bridge samples
#   from multiple batches. Optional step that is only used in multi-batch mode
#   for cross-batch data integration.
#
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Date: December 2025
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(OlinkAnalyze)
  library(ggplot2)
  library(ggpubr)
  library(yaml)
  library(logger)
  library(sva)  # For ComBat normalization
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
# Use step_num (which is "07") instead of hardcoded "09_bridge_normalization"
log_path <- get_log_path(step_num, batch_id, config = config)
log_appender(appender_file(log_path))
log_info("Starting enhanced bridge normalisation for batch: {batch_id}")

# Set theme for plots
theme_set(theme_bw())

# Function to validate bridge samples for normalisation
validate_bridge_samples <- function(npx_matrix, bridge_samples) {
  log_info("Validating bridge samples for normalisation")

  # Check bridge sample availability
  available_bridges <- intersect(rownames(npx_matrix), bridge_samples)
  missing_bridges <- setdiff(bridge_samples, rownames(npx_matrix))

  log_info("Available bridge samples: {length(available_bridges)}")
  if(length(missing_bridges) > 0) {
    log_warn("Missing bridge samples: {length(missing_bridges)}")
  }

  # Extract bridge sample data
  bridge_data <- npx_matrix[available_bridges, ]

  # Check quality of bridge samples
  bridge_quality <- data.table(
    SampleID = available_bridges,
    missing_rate = rowSums(is.na(bridge_data)) / ncol(bridge_data),
    mean_expr = rowMeans(bridge_data, na.rm = TRUE),
    sd_expr = apply(bridge_data, 1, sd, na.rm = TRUE)
  )

  # Flag low-quality bridge samples
  bridge_quality[, high_quality := missing_rate < 0.2]

  high_quality_bridges <- bridge_quality[high_quality == TRUE]$SampleID

  log_info("High-quality bridge samples: {length(high_quality_bridges)}")

  return(list(
    available_bridges = available_bridges,
    high_quality_bridges = high_quality_bridges,
    bridge_quality = bridge_quality
  ))
}

# Function for cross-batch bridge normalisation (loads bridge samples from both batches)
cross_batch_bridge_normalization <- function(batch1_matrix, batch2_matrix, batch1_bridge_finngenids, batch2_bridge_finngenids,
                                             batch1_sample_mapping, batch2_sample_mapping, method = "median") {
  log_info("╔══════════════════════════════════════════════════════════════════╗")
  log_info("║ CROSS-BATCH BRIDGE NORMALIZATION                                 ║")
  log_info("╚══════════════════════════════════════════════════════════════════╝")
  log_info("Method: {method}")
  log_info("")
  log_info("=== INPUT MATRICES ===")
  log_info("  Batch 1: {nrow(batch1_matrix)} samples × {ncol(batch1_matrix)} proteins")
  log_info("  Batch 2: {nrow(batch2_matrix)} samples × {ncol(batch2_matrix)} proteins")

  # Match bridge samples by FINNGENID across batches
  common_bridge_finngenids <- intersect(batch1_bridge_finngenids, batch2_bridge_finngenids)
  log_info("")
  log_info("=== BRIDGE SAMPLES ===")
  log_info("  Batch 1 bridge FINNGENIDs available: {length(batch1_bridge_finngenids)}")
  log_info("  Batch 2 bridge FINNGENIDs available: {length(batch2_bridge_finngenids)}")
  log_info("  Common bridge FINNGENIDs: {length(common_bridge_finngenids)}")

  if(length(common_bridge_finngenids) < 10) {
    log_error("Insufficient common bridge samples ({length(common_bridge_finngenids)} < 10)")
    return(NULL)
  }

  # Get SampleIDs for bridge samples in each batch
  batch1_bridge_samples <- batch1_sample_mapping[FINNGENID %in% common_bridge_finngenids]$SampleID
  batch2_bridge_samples <- batch2_sample_mapping[FINNGENID %in% common_bridge_finngenids]$SampleID

  # Filter to samples that exist in matrices
  batch1_bridge_in_matrix <- intersect(batch1_bridge_samples, rownames(batch1_matrix))
  batch2_bridge_in_matrix <- intersect(batch2_bridge_samples, rownames(batch2_matrix))

  log_info("  Batch 1 bridge samples in matrix: {length(batch1_bridge_in_matrix)}")
  log_info("  Batch 2 bridge samples in matrix: {length(batch2_bridge_in_matrix)}")

  if(length(batch1_bridge_in_matrix) < 10 || length(batch2_bridge_in_matrix) < 10) {
    log_error("Insufficient bridge samples in matrices (batch1: {length(batch1_bridge_in_matrix)}, batch2: {length(batch2_bridge_in_matrix)})")
    return(NULL)
  }

  # Find common proteins
  common_proteins <- intersect(colnames(batch1_matrix), colnames(batch2_matrix))
  batch1_only_proteins <- setdiff(colnames(batch1_matrix), common_proteins)
  batch2_only_proteins <- setdiff(colnames(batch2_matrix), common_proteins)

  log_info("")
  log_info("=== PROTEINS ===")
  log_info("  Common proteins: {length(common_proteins)}")
  log_info("  Batch 1 only: {length(batch1_only_proteins)}")
  log_info("  Batch 2 only: {length(batch2_only_proteins)}")

  if(length(common_proteins) < 100) {
    log_error("Too few common proteins ({length(common_proteins)} < 100)")
    return(NULL)
  }

  # Extract bridge data from both batches (common proteins only)
  batch1_bridge_data <- batch1_matrix[batch1_bridge_in_matrix, common_proteins, drop = FALSE]
  batch2_bridge_data <- batch2_matrix[batch2_bridge_in_matrix, common_proteins, drop = FALSE]

  if(method == "median") {
    log_info("")
    log_info("=== CALCULATING BATCH OFFSETS (ADDITIVE ADJUSTMENT) ===")
    log_info("NOTE: NPX values are already log2-transformed. Using ADDITIVE adjustment (not multiplicative)")

    # Calculate per-protein medians from bridge samples in each batch
    batch1_medians <- apply(batch1_bridge_data, 2, median, na.rm = TRUE)
    batch2_medians <- apply(batch2_bridge_data, 2, median, na.rm = TRUE)

    log_info("Bridge sample median statistics (per-protein, then aggregated):")
    log_info("  Batch 1: min={round(min(batch1_medians, na.rm=TRUE), 3)}, max={round(max(batch1_medians, na.rm=TRUE), 3)}, median={round(median(batch1_medians, na.rm=TRUE), 3)}")
    log_info("  Batch 2: min={round(min(batch2_medians, na.rm=TRUE), 3)}, max={round(max(batch2_medians, na.rm=TRUE), 3)}, median={round(median(batch2_medians, na.rm=TRUE), 3)}")

    # Combined reference median (average of batch medians per protein)
    reference_medians <- (batch1_medians + batch2_medians) / 2
    log_info("Reference (combined) median range: [{round(min(reference_medians, na.rm=TRUE), 3)}, {round(max(reference_medians, na.rm=TRUE), 3)}]")

    # Calculate ADDITIVE offsets for each batch (what to subtract to align with reference)
    # offset = batch_median - reference_median
    # NPX_adjusted = NPX - offset  (shifts distribution to align with reference)
    batch1_offsets <- batch1_medians - reference_medians
    batch2_offsets <- batch2_medians - reference_medians

    # Handle invalid offsets (set to 0 = no adjustment)
    n_invalid_b1 <- sum(is.na(batch1_offsets))
    n_invalid_b2 <- sum(is.na(batch2_offsets))
    batch1_offsets[is.na(batch1_offsets)] <- 0
    batch2_offsets[is.na(batch2_offsets)] <- 0

    # Cap extreme offsets (more than ±3 log2 units is suspicious)
    offset_cap <- 3.0
    n_capped_b1 <- sum(abs(batch1_offsets) > offset_cap)
    n_capped_b2 <- sum(abs(batch2_offsets) > offset_cap)
    batch1_offsets[batch1_offsets > offset_cap] <- offset_cap
    batch1_offsets[batch1_offsets < -offset_cap] <- -offset_cap
    batch2_offsets[batch2_offsets > offset_cap] <- offset_cap
    batch2_offsets[batch2_offsets < -offset_cap] <- -offset_cap

    log_info("")
    log_info("=== OFFSET STATISTICS (log2 units) ===")
    log_info("Batch 1 offsets:")
    log_info("  Range: [{round(min(batch1_offsets), 4)}, {round(max(batch1_offsets), 4)}]")
    log_info("  Mean: {round(mean(batch1_offsets), 4)}, Median: {round(median(batch1_offsets), 4)}")
    log_info("  Invalid (set to 0): {n_invalid_b1}, Capped (|offset| > {offset_cap}): {n_capped_b1}")
    log_info("Batch 2 offsets:")
    log_info("  Range: [{round(min(batch2_offsets), 4)}, {round(max(batch2_offsets), 4)}]")
    log_info("  Mean: {round(mean(batch2_offsets), 4)}, Median: {round(median(batch2_offsets), 4)}")
    log_info("  Invalid (set to 0): {n_invalid_b2}, Capped (|offset| > {offset_cap}): {n_capped_b2}")

    # Apply ADDITIVE adjustment to each batch (common proteins only)
    # NPX_adjusted = NPX - offset (subtracting shifts toward reference)
    log_info("")
    log_info("=== APPLYING ADDITIVE NORMALISATION ===")
    log_info("Formula: NPX_adjusted = NPX - offset (shifts distributions to align with reference)")
    log_info("Adjusting {length(common_proteins)} common proteins...")

    batch1_normalized_common <- sweep(batch1_matrix[, common_proteins, drop = FALSE], 2, batch1_offsets, "-")
    batch2_normalized_common <- sweep(batch2_matrix[, common_proteins, drop = FALSE], 2, batch2_offsets, "-")

    # Reconstruct full matrices (non-common proteins unchanged)
    batch1_normalized <- batch1_matrix
    batch1_normalized[, common_proteins] <- batch1_normalized_common
    batch2_normalized <- batch2_matrix
    batch2_normalized[, common_proteins] <- batch2_normalized_common

    log_info("Non-common proteins: unchanged (no cross-batch reference available)")

    # NO additional log transformation - NPX is already log2-transformed!
    log_info("NOTE: No log2 transformation applied - NPX values are already on log2 scale")

    log_info("")
    log_info("=== OUTPUT MATRICES ===")
    log_info("  Batch 1 normalized: {nrow(batch1_normalized)} × {ncol(batch1_normalized)}")
    log_info("  Batch 2 normalized: {nrow(batch2_normalized)} × {ncol(batch2_normalized)}")
    log_info("════════════════════════════════════════════════════════════════════")

    return(list(
      batch1_normalized = batch1_normalized,       # Already on log2 scale
      batch2_normalized = batch2_normalized,       # Already on log2 scale
      batch1_normalized_raw = batch1_normalized,   # Same as above (no transform needed)
      batch2_normalized_raw = batch2_normalized,   # Same as above (no transform needed)
      batch1_offsets = batch1_offsets,             # Renamed from scaling_factors
      batch2_offsets = batch2_offsets,             # Renamed from scaling_factors
      batch1_bridge_samples_used = batch1_bridge_in_matrix,
      batch2_bridge_samples_used = batch2_bridge_in_matrix,
      common_proteins = common_proteins,
      method = "bridge"
    ))
  } else if(method == "quantile") {
    log_info("Performing cross-batch quantile normalisation")
    # Quantile normalization using combined bridge samples
    combined_bridge_data <- rbind(batch1_bridge_data, batch2_bridge_data)

    # Get reference distribution from combined bridge samples
    bridge_quantiles <- apply(combined_bridge_data, 2, function(x) {
      quantile(x, probs = seq(0, 1, 0.01), na.rm = TRUE)
    })
    reference_quantiles <- rowMeans(bridge_quantiles, na.rm = TRUE)

    # Apply quantile normalisation to both batches
    normalize_quantile_batch <- function(matrix_subset, reference_quantiles) {
      apply(matrix_subset, 2, function(x) {
        x_ranks <- rank(x, na.last = "keep", ties.method = "average")
        x_probs <- x_ranks / (sum(!is.na(x)) + 1)
        x_normalized <- approx(seq(0, 1, 0.01), reference_quantiles, xout = x_probs)$y
        return(x_normalized)
      })
    }

    batch1_normalized_common <- normalize_quantile_batch(batch1_matrix[, common_proteins, drop = FALSE], reference_quantiles)
    batch2_normalized_common <- normalize_quantile_batch(batch2_matrix[, common_proteins, drop = FALSE], reference_quantiles)

    # Reconstruct full matrices
    batch1_normalized <- batch1_matrix
    batch1_normalized[, common_proteins] <- batch1_normalized_common
    batch2_normalized <- batch2_matrix
    batch2_normalized[, common_proteins] <- batch2_normalized_common

    # Apply log2 transformation: log2(NPX_normalized + 1)
    batch1_normalized_log2 <- log2(batch1_normalized + 1)
    batch2_normalized_log2 <- log2(batch2_normalized + 1)

    return(list(
      batch1_normalized = batch1_normalized_log2,
      batch2_normalized = batch2_normalized_log2,
      batch1_normalized_raw = batch1_normalized,  # Keep raw normalized for evaluation
      batch2_normalized_raw = batch2_normalized,  # Keep raw normalized for evaluation
      batch1_bridge_samples_used = batch1_bridge_in_matrix,
      batch2_bridge_samples_used = batch2_bridge_in_matrix,
      common_proteins = common_proteins,
      reference_quantiles = reference_quantiles,
      method = "quantile"
    ))
  } else if(method == "combat") {
    log_info("Performing cross-batch ComBat normalization")

    # Combine matrices with batch labels
    combined_matrix <- rbind(
      batch1_matrix[, common_proteins, drop = FALSE],
      batch2_matrix[, common_proteins, drop = FALSE]
    )
    batch_vector <- c(rep("batch_01", nrow(batch1_matrix)),
                      rep("batch_02", nrow(batch2_matrix)))

    # Transpose for ComBat (proteins as rows, samples as columns)
    expr_matrix <- t(combined_matrix)

    # Remove rows with all NAs
    complete_proteins <- rowSums(is.na(expr_matrix)) < ncol(expr_matrix)
    expr_matrix <- expr_matrix[complete_proteins, ]

    # Impute remaining NAs with row means
    for(i in 1:nrow(expr_matrix)) {
      expr_matrix[i, is.na(expr_matrix[i, ])] <- mean(expr_matrix[i, ], na.rm = TRUE)
    }

    # Apply ComBat
    tryCatch({
      combat_expr <- ComBat(dat = expr_matrix,
                           batch = batch_vector,
                           mod = NULL,  # No covariates for now
                           par.prior = TRUE,
                           prior.plots = FALSE)

      # Transpose back
      combined_normalized <- t(combat_expr)

      # Split back into batches
      batch1_normalized_common <- combined_normalized[1:nrow(batch1_matrix), ]
      batch2_normalized_common <- combined_normalized[(nrow(batch1_matrix)+1):nrow(combined_normalized), ]

      # Reconstruct full matrices (non-common proteins unchanged)
      batch1_normalized <- batch1_matrix
      batch1_normalized[, common_proteins] <- batch1_normalized_common
      batch2_normalized <- batch2_matrix
      batch2_normalized[, common_proteins] <- batch2_normalized_common

      # Apply log2 transformation: log2(NPX_normalized + 1)
      batch1_normalized_log2 <- log2(batch1_normalized + 1)
      batch2_normalized_log2 <- log2(batch2_normalized + 1)

      log_info("ComBat normalization successful")

      return(list(
        batch1_normalized = batch1_normalized_log2,
        batch2_normalized = batch2_normalized_log2,
        batch1_normalized_raw = batch1_normalized,  # Keep raw normalized for evaluation
        batch2_normalized_raw = batch2_normalized,  # Keep raw normalized for evaluation
        batch1_bridge_samples_used = batch1_bridge_in_matrix,
        batch2_bridge_samples_used = batch2_bridge_in_matrix,
        common_proteins = common_proteins,
        method = "combat"
      ))
    }, error = function(e) {
      log_error("ComBat normalization failed: {e$message}")
      return(NULL)
    })
  }

  return(NULL)
}

# Function for enhanced bridge normalisation (single-batch fallback, deprecated in cross-batch mode)
enhanced_bridge_normalization <- function(npx_matrix, bridge_samples, method = "median") {
  log_info("Performing enhanced bridge normalisation with method: {method}")

  # Validate bridge samples
  validation <- validate_bridge_samples(npx_matrix, bridge_samples)

  if(length(validation$high_quality_bridges) < 10) {
    log_error("Insufficient high-quality bridge samples for normalisation")
    return(NULL)
  }

  # Use high-quality bridge samples
  bridge_data <- npx_matrix[validation$high_quality_bridges, ]

  if(method == "median") {
    # Median-based normalisation
    bridge_medians <- apply(bridge_data, 2, median, na.rm = TRUE)
    global_median <- median(bridge_medians, na.rm = TRUE)
    scaling_factors <- global_median / bridge_medians

  } else if(method == "mean") {
    # Mean-based normalisation
    bridge_means <- apply(bridge_data, 2, mean, na.rm = TRUE)
    global_mean <- mean(bridge_means, na.rm = TRUE)
    scaling_factors <- global_mean / bridge_means

  } else if(method == "quantile") {
    # Quantile normalisation using bridge samples
    log_info("Performing quantile normalisation")

    # Get reference distribution from bridge samples
    bridge_quantiles <- apply(bridge_data, 2, function(x) {
      quantile(x, probs = seq(0, 1, 0.01), na.rm = TRUE)
    })
    reference_quantiles <- rowMeans(bridge_quantiles, na.rm = TRUE)

    # Apply quantile normalisation to full matrix
    normalized_matrix <- apply(npx_matrix, 2, function(x) {
      x_ranks <- rank(x, na.last = "keep", ties.method = "average")
      x_probs <- x_ranks / (sum(!is.na(x)) + 1)
      x_normalized <- approx(seq(0, 1, 0.01), reference_quantiles, xout = x_probs)$y
      return(x_normalized)
    })

    return(list(
      normalized_matrix = normalized_matrix,
      method = "quantile",
      bridge_samples_used = validation$high_quality_bridges,
      reference_quantiles = reference_quantiles
    ))
  }

  # Apply scaling factors for median/mean methods
  if(method %in% c("median", "mean")) {
    scaling_factors[is.na(scaling_factors) | is.infinite(scaling_factors)] <- 1

    # Cap extreme scaling factors to prevent outliers (keep within 0.1 to 10 range)
    scaling_factors[scaling_factors < 0.1] <- 0.1
    scaling_factors[scaling_factors > 10] <- 10

    normalized_matrix <- sweep(npx_matrix, 2, scaling_factors, "*")

    return(list(
      normalized_matrix = normalized_matrix,
      scaling_factors = scaling_factors,
      method = method,
      bridge_samples_used = validation$high_quality_bridges
    ))
  }
}

# Function to harmonise with reference batch
harmonize_with_reference <- function(current_batch, reference_batch = NULL, bridge_samples) {
  log_info("Harmonising with reference batch")

  if(is.null(reference_batch)) {
    log_info("No reference batch provided, using internal normalisation")
    return(enhanced_bridge_normalization(current_batch, bridge_samples, method = "median"))
  }

  # This function would be implemented when batch 1 data is available
  # For now, return placeholder
  log_info("Reference batch harmonization to be implemented")

  return(list(
    harmonized_matrix = current_batch,
    message = "Reference batch harmonization pending implementation"
  ))
}

# Function to evaluate normalisation quality
evaluate_bridge_normalization <- function(raw_matrix, normalized_matrix, bridge_samples) {
  log_info("Evaluating bridge normalisation quality")

  # Get bridge sample indices
  bridge_idx <- which(rownames(raw_matrix) %in% bridge_samples)
  non_bridge_idx <- which(!rownames(raw_matrix) %in% bridge_samples)

  # Calculate CV before and after for bridge samples
  cv_bridge_before <- apply(raw_matrix[bridge_idx, ], 2, function(x) {
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  })

  cv_bridge_after <- apply(normalized_matrix[bridge_idx, ], 2, function(x) {
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  })

  # Calculate CV for non-bridge samples
  cv_non_bridge_before <- apply(raw_matrix[non_bridge_idx, ], 2, function(x) {
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  })

  cv_non_bridge_after <- apply(normalized_matrix[non_bridge_idx, ], 2, function(x) {
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  })

  # Calculate summary statistics
  mean_cv_bridge_before <- mean(cv_bridge_before, na.rm = TRUE)
  mean_cv_bridge_after <- mean(cv_bridge_after, na.rm = TRUE)
  mean_cv_non_bridge_before <- mean(cv_non_bridge_before, na.rm = TRUE)
  mean_cv_non_bridge_after <- mean(cv_non_bridge_after, na.rm = TRUE)

  median_cv_bridge_before <- median(cv_bridge_before, na.rm = TRUE)
  median_cv_bridge_after <- median(cv_bridge_after, na.rm = TRUE)
  median_cv_non_bridge_before <- median(cv_non_bridge_before, na.rm = TRUE)
  median_cv_non_bridge_after <- median(cv_non_bridge_after, na.rm = TRUE)

  # Calculate percentage reductions
  pct_reduction_mean_cv_bridge <- ifelse(mean_cv_bridge_before > 0,
                                        (mean_cv_bridge_before - mean_cv_bridge_after) / mean_cv_bridge_before * 100,
                                        NA_real_)
  pct_reduction_mean_cv_non_bridge <- ifelse(mean_cv_non_bridge_before > 0,
                                             (mean_cv_non_bridge_before - mean_cv_non_bridge_after) / mean_cv_non_bridge_before * 100,
                                             NA_real_)
  pct_reduction_median_cv_bridge <- ifelse(median_cv_bridge_before > 0,
                                          (median_cv_bridge_before - median_cv_bridge_after) / median_cv_bridge_before * 100,
                                          NA_real_)
  pct_reduction_median_cv_non_bridge <- ifelse(median_cv_non_bridge_before > 0,
                                               (median_cv_non_bridge_before - median_cv_non_bridge_after) / median_cv_non_bridge_before * 100,
                                               NA_real_)

  # Summary statistics
  eval_stats <- data.table(
    sample_type = c("Bridge", "Bridge", "Non-Bridge", "Non-Bridge"),
    stage = c("Before", "After", "Before", "After"),
    mean_cv = c(mean_cv_bridge_before,
                mean_cv_bridge_after,
                mean_cv_non_bridge_before,
                mean_cv_non_bridge_after),
    median_cv = c(median_cv_bridge_before,
                  median_cv_bridge_after,
                  median_cv_non_bridge_before,
                  median_cv_non_bridge_after),
    pct_reduction_mean_cv = c(NA_real_,  # Before stage has no reduction
                              pct_reduction_mean_cv_bridge,
                              NA_real_,  # Before stage has no reduction
                              pct_reduction_mean_cv_non_bridge),
    pct_reduction_median_cv = c(NA_real_,  # Before stage has no reduction
                               pct_reduction_median_cv_bridge,
                               NA_real_,  # Before stage has no reduction
                               pct_reduction_median_cv_non_bridge)
  )

  # Calculate improvement
  bridge_improvement <- mean(cv_bridge_before - cv_bridge_after, na.rm = TRUE)
  overall_improvement <- mean(c(cv_bridge_before, cv_non_bridge_before) -
                             c(cv_bridge_after, cv_non_bridge_after), na.rm = TRUE)

  log_info("Bridge sample CV improvement: {round(bridge_improvement, 4)}")
  log_info("Overall CV improvement: {round(overall_improvement, 4)}")

  return(list(
    eval_stats = eval_stats,
    bridge_improvement = bridge_improvement,
    overall_improvement = overall_improvement
  ))
}

# Function to load all batch matrices for cross-batch comparison
load_all_batch_matrices <- function(current_batch_id, config, input_choice = "qc_passed") {
  log_info("Loading all batch matrices for cross-batch comparison")

  # Get all batch IDs from config
  # Check if config$batches exists and is not NULL
  if (is.null(config$batches) || !is.list(config$batches)) {
    log_warn("config$batches is not available or not a list - cannot load all batches")
    return(list())
  }

  all_batches <- tryCatch(names(config$batches), error = function(e) NULL)
  if (is.null(all_batches) || length(all_batches) < 2) {
    log_warn("Multi-batch mode not properly configured - cannot load all batches")
    log_warn("  Found {length(all_batches) %||% 0} batches in config")
    return(list())
  }

  batch_matrices <- list()

  for (batch_id in all_batches) {
    if (input_choice == "qc_passed") {
      matrix_path <- get_output_path("05d", "05d_npx_matrix_all_qc_passed", batch_id, "phenotypes", "rds", config = config)
      if (!file.exists(matrix_path)) {
        matrix_path <- get_output_path("00", "npx_matrix_analysis_ready", batch_id, "qc", config = config)
      }
    } else {
      matrix_path <- get_output_path("00", "npx_matrix_analysis_ready", batch_id, "qc", config = config)
    }

    if (file.exists(matrix_path)) {
      batch_matrices[[batch_id]] <- readRDS(matrix_path)
      log_info("Loaded {batch_id}: {nrow(batch_matrices[[batch_id]])} samples × {ncol(batch_matrices[[batch_id]])} proteins")
    } else {
      log_warn("Matrix not found for {batch_id}: {matrix_path}")
    }
  }

  return(batch_matrices)
}

# Function to create cross-batch comparison box plots using ggpubr
create_cross_batch_comparison_plots <- function(batch_matrices_raw, batch_matrices_normalized, method_name = "Bridge") {
  log_info("Creating cross-batch comparison plots using ggpubr")

  if (length(batch_matrices_raw) < 2) {
    log_warn("Need at least 2 batches for cross-batch comparison plots")
    return(NULL)
  }

  # Ensure batch_matrices_normalized has the same batches as batch_matrices_raw
  common_batches <- intersect(names(batch_matrices_raw), names(batch_matrices_normalized))
  if (length(common_batches) < 2) {
    log_warn("Insufficient common batches between raw and normalized matrices ({length(common_batches)} < 2)")
    return(NULL)
  }

  # Filter to only common batches
  batch_matrices_raw <- batch_matrices_raw[common_batches]
  batch_matrices_normalized <- batch_matrices_normalized[common_batches]

  # Find common proteins across all batches
  all_proteins <- Reduce(intersect, lapply(batch_matrices_raw, colnames))
  if (length(all_proteins) < 100) {
    log_warn("Too few common proteins ({length(all_proteins)}) for meaningful comparison")
    return(NULL)
  }

  log_info("Using {length(all_proteins)} common proteins for cross-batch comparison across {length(common_batches)} batches")

  # Randomly sample proteins for plotting
  set.seed(123)
  sample_proteins <- sample(all_proteins, min(50, length(all_proteins)))

  # Prepare data for each batch
  plot_data_list <- list()

  for (batch_id in names(batch_matrices_raw)) {
    # Defensive check: ensure batch_id exists in both lists
    if (!batch_id %in% names(batch_matrices_normalized)) {
      log_warn("Skipping {batch_id} - not found in normalized matrices")
      next
    }

    raw_matrix <- batch_matrices_raw[[batch_id]][, sample_proteins, drop = FALSE]
    norm_matrix <- batch_matrices_normalized[[batch_id]][, sample_proteins, drop = FALSE]

    # Calculate mean expression per sample (across sampled proteins)
    raw_means <- rowMeans(raw_matrix, na.rm = TRUE)
    norm_means <- rowMeans(norm_matrix, na.rm = TRUE)

    plot_data_list[[batch_id]] <- data.table(
      Batch = batch_id,
      Before = raw_means,
      After = norm_means
    )
  }

  # Combine all batches
  if (length(plot_data_list) == 0) {
    log_warn("No plot data generated - all batches were skipped")
    return(NULL)
  }

  plot_data <- rbindlist(plot_data_list)
  if (nrow(plot_data) == 0) {
    log_warn("Empty plot data after combining batches")
    return(NULL)
  }

  plot_data <- melt(plot_data, id.vars = "Batch",
                     measure.vars = c("Before", "After"),
                     variable.name = "Stage", value.name = "MeanExpression")

  # Create comparison plot using ggpubr
  p1 <- ggboxplot(plot_data, x = "Batch", y = "MeanExpression",
                  fill = "Stage", palette = c("#FF6B6B", "#3A5F8A"),
                  add = "jitter", add.params = list(alpha = 0.3, size = 0.5),
                  outlier.shape = NA) +
    stat_compare_means(aes(group = Stage), method = "t.test",
                       label = "p.format", label.y = max(plot_data$MeanExpression, na.rm = TRUE) * 1.1) +
    labs(title = paste("Cross-Batch Normalisation Comparison:", method_name),
         subtitle = "Mean expression per sample (across 50 randomly selected proteins)",
         x = "Batch", y = "Mean NPX Expression",
         fill = "Stage") +
    theme_bw() +
    theme(legend.position = "right",
          plot.subtitle = element_text(size = 9, color = "gray40"))

  # Create before/after comparison plot (boxplot faceted by batch)
  p2 <- ggboxplot(plot_data, x = "Stage", y = "MeanExpression",
                  fill = "Stage", palette = c("#FF6B6B", "#3A5F8A"),
                  add = "jitter", add.params = list(alpha = 0.3, size = 0.5),
                  facet.by = "Batch", ncol = min(3, length(unique(plot_data$Batch))),
                  outlier.shape = NA) +
    stat_compare_means(aes(group = Stage), method = "t.test",
                       label = "p.format", label.y = max(plot_data$MeanExpression, na.rm = TRUE) * 1.15,
                       comparisons = list(c("Before", "After"))) +
    labs(title = paste("Before vs After Normalisation by Batch:", method_name),
         subtitle = "Comparison of mean expression per sample across batches",
         x = "Stage", y = "Mean NPX Expression",
         fill = "Stage") +
    theme_bw() +
    theme(legend.position = "right",
          plot.subtitle = element_text(size = 9, color = "gray40"),
          strip.background = element_rect(fill = "gray95"),
          strip.text = element_text(face = "bold"))

  return(list(
    cross_batch_boxplot = p1,
    before_after_paired = p2
  ))
}

# Function to create normalization plots with alternating batches (similar to step 6)
create_alternating_batch_plots <- function(batch_matrices_raw, batch_matrices_normalized,
                                          title_suffix = "", n_samples_per_batch = 50) {
  log_info("Creating normalisation plots with alternating batches")

  if (length(batch_matrices_raw) < 2) {
    log_warn("Need at least 2 batches for alternating batch plots")
    return(NULL)
  }

  # Ensure batch_matrices_normalized has the same batches as batch_matrices_raw
  common_batches <- intersect(names(batch_matrices_raw), names(batch_matrices_normalized))
  if (length(common_batches) < 2) {
    log_warn("Insufficient common batches between raw and normalized matrices ({length(common_batches)} < 2)")
    return(NULL)
  }

  # Filter to only common batches
  batch_matrices_raw <- batch_matrices_raw[common_batches]
  batch_matrices_normalized <- batch_matrices_normalized[common_batches]

  # Find common proteins
  all_proteins <- Reduce(intersect, lapply(batch_matrices_raw, colnames))
  if (length(all_proteins) < 100) {
    log_warn("Too few common proteins for plotting")
    return(NULL)
  }

  # Randomly sample proteins
  set.seed(123)
  sample_proteins <- sample(all_proteins, min(100, length(all_proteins)))

  # Sample samples from each batch (alternating pattern)
  set.seed(456)
  sampled_data_list <- list()
  sample_counter <- 1

  for (batch_id in names(batch_matrices_raw)) {
    # Defensive check: ensure batch_id exists in both lists
    if (!batch_id %in% names(batch_matrices_normalized)) {
      log_warn("Skipping {batch_id} - not found in normalized matrices")
      next
    }

    raw_matrix <- batch_matrices_raw[[batch_id]][, sample_proteins, drop = FALSE]
    norm_matrix <- batch_matrices_normalized[[batch_id]][, sample_proteins, drop = FALSE]

    n_samples_actual <- min(n_samples_per_batch, nrow(raw_matrix))
    sample_idx <- sample(nrow(raw_matrix), n_samples_actual)

    # Before normalisation
    df_before <- as.data.table(raw_matrix[sample_idx, ], keep.rownames = TRUE)
    setnames(df_before, "rn", "Sample")
    df_before <- melt(df_before, id.vars = "Sample", variable.name = "Protein", value.name = "Expression")
    df_before$Stage <- "Before"
    df_before$Batch <- batch_id
    # Repeat SampleOrder for each protein (after melting, each sample has n_proteins rows)
    sample_order_vec <- sample_counter:(sample_counter + n_samples_actual - 1)
    df_before$SampleOrder <- rep(sample_order_vec, each = length(sample_proteins))
    sample_counter <- sample_counter + n_samples_actual

    # After normalisation
    df_after <- as.data.table(norm_matrix[sample_idx, ], keep.rownames = TRUE)
    setnames(df_after, "rn", "Sample")
    df_after <- melt(df_after, id.vars = "Sample", variable.name = "Protein", value.name = "Expression")
    df_after$Stage <- "After"
    df_after$Batch <- batch_id
    # Repeat SampleOrder for each protein (same as df_before)
    df_after$SampleOrder <- rep(sample_order_vec, each = length(sample_proteins))

    sampled_data_list[[batch_id]] <- rbind(df_before, df_after)
  }

  # Combine all batches
  if (length(sampled_data_list) == 0) {
    log_warn("No sampled data generated - all batches were skipped")
    return(NULL)
  }

  df_combined <- rbindlist(sampled_data_list)
  if (nrow(df_combined) == 0) {
    log_warn("Empty combined data after merging batches")
    return(NULL)
  }

  df_combined$Stage <- factor(df_combined$Stage, levels = c("Before", "After"))
  df_combined$Batch <- factor(df_combined$Batch, levels = names(batch_matrices_raw))

  # Calculate summary statistics
  stats_before <- df_combined[Stage == "Before", .(
    Mean = round(mean(Expression, na.rm = TRUE), 2),
    Median = round(median(Expression, na.rm = TRUE), 2),
    SD = round(sd(Expression, na.rm = TRUE), 2),
    IQR = round(IQR(Expression, na.rm = TRUE), 2)
  )]

  stats_after <- df_combined[Stage == "After", .(
    Mean = round(mean(Expression, na.rm = TRUE), 2),
    Median = round(median(Expression, na.rm = TRUE), 2),
    SD = round(sd(Expression, na.rm = TRUE), 2),
    IQR = round(IQR(Expression, na.rm = TRUE), 2)
  )]

  # Create subtitle
  total_samples <- sum(sapply(batch_matrices_raw, nrow))
  total_proteins <- length(sample_proteins)
  subtitle <- sprintf(
    "Showing %d randomly selected samples per batch (alternating) × %d proteins\nBefore: Mean=%.2f, SD=%.2f, IQR=%.2f | After: Mean=%.2f, SD=%.2f, IQR=%.2f",
    n_samples_per_batch, total_proteins,
    stats_before$Mean, stats_before$SD, stats_before$IQR,
    stats_after$Mean, stats_after$SD, stats_after$IQR
  )

  # Create boxplot with alternating batches
  p <- ggplot(df_combined, aes(x = as.factor(SampleOrder), y = Expression, fill = Batch)) +
    geom_boxplot(alpha = 0.6, outlier.size = 0.5) +
    facet_wrap(~ Stage, scales = "free_y", ncol = 1) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = paste("Normalization Effect with Alternating Batches", title_suffix),
         subtitle = subtitle,
         x = "Sample (alternating between batches, ordered by index)",
         y = "NPX Expression Level (log2 scale)",
         fill = "Batch") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      legend.position = "right"
    )

  return(p)
}

# Function to create bridge normalisation plots (original function, kept for backward compatibility)
create_bridge_plots <- function(raw_matrix, normalized_matrix, bridge_samples) {
  log_info("Creating bridge normalisation visualisation")

  # Get bridge and non-bridge samples
  bridge_idx <- rownames(raw_matrix) %in% bridge_samples

  # Sample data for plotting
  set.seed(123)
  sample_proteins <- sample(ncol(raw_matrix), min(100, ncol(raw_matrix)))

  # Prepare data
  plot_data <- rbind(
    data.table(
      value = as.vector(raw_matrix[, sample_proteins]),
      stage = "Before",
      is_bridge = rep(bridge_idx, length(sample_proteins))
    ),
    data.table(
      value = as.vector(normalized_matrix[, sample_proteins]),
      stage = "After",
      is_bridge = rep(bridge_idx, length(sample_proteins))
    )
  )

  plot_data[, sample_type := ifelse(is_bridge, "Bridge", "Non-Bridge")]

  # Create boxplot
  p1 <- ggplot(plot_data[!is.na(value)], aes(x = sample_type, y = value, fill = stage)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    scale_fill_manual(values = c("Before" = "#FF6B6B", "After" = "#3A5F8A")) +
    labs(title = "Bridge Normalization Effect",
         x = "Sample Type", y = "Expression Value",
         fill = "Stage") +
    theme_bw() +
    facet_wrap(~ stage, scales = "free_y")

  # Density plot
  p2 <- ggplot(plot_data[!is.na(value)], aes(x = value, fill = sample_type)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("Bridge" = "#5C9EAD", "Non-Bridge" = "#87CEEB")) +
    labs(title = "Expression Distribution",
         x = "Expression Value", y = "Density",
         fill = "Sample Type") +
    theme_bw() +
    facet_wrap(~ stage, scales = "free")

  return(list(
    boxplot = p1,
    density = p2
  ))
}

# Function for median normalization (global variant) - FIXED LOGIC (ADDITIVE)
cross_batch_median_normalization <- function(batch1_matrix, batch2_matrix, common_proteins, reference_batch = "batch_02") {
  log_info("Performing cross-batch median normalization (global variant, ADDITIVE)")
  log_info("NOTE: Using ADDITIVE adjustment for log2-transformed NPX data")

  # Calculate per-protein median across ALL samples (not just bridge)
  batch1_all_medians <- apply(batch1_matrix[, common_proteins, drop = FALSE], 2, median, na.rm = TRUE)
  batch2_all_medians <- apply(batch2_matrix[, common_proteins, drop = FALSE], 2, median, na.rm = TRUE)

  # Use reference batch medians as target
  if (reference_batch == "batch_02") {
    reference_medians <- batch2_all_medians
    log_info("Reference batch: batch_02")
  } else {
    reference_medians <- batch1_all_medians
    log_info("Reference batch: batch_01")
  }

  # Calculate ADDITIVE offsets relative to reference
  # offset = batch_median - reference_median
  # NPX_adjusted = NPX - offset
  batch1_offsets <- batch1_all_medians - reference_medians
  batch2_offsets <- batch2_all_medians - reference_medians

  # Handle invalid offsets
  batch1_offsets[is.na(batch1_offsets)] <- 0
  batch2_offsets[is.na(batch2_offsets)] <- 0

  # Cap extreme offsets (±3 log2 units)
  offset_cap <- 3.0
  batch1_offsets[abs(batch1_offsets) > offset_cap] <- sign(batch1_offsets[abs(batch1_offsets) > offset_cap]) * offset_cap
  batch2_offsets[abs(batch2_offsets) > offset_cap] <- sign(batch2_offsets[abs(batch2_offsets) > offset_cap]) * offset_cap

  log_info("Batch 1 offset range: [{round(min(batch1_offsets), 3)}, {round(max(batch1_offsets), 3)}]")
  log_info("Batch 2 offset range: [{round(min(batch2_offsets), 3)}, {round(max(batch2_offsets), 3)}]")

  # Apply ADDITIVE adjustment
  batch1_normalized_common <- sweep(batch1_matrix[, common_proteins, drop = FALSE], 2, batch1_offsets, "-")
  batch2_normalized_common <- sweep(batch2_matrix[, common_proteins, drop = FALSE], 2, batch2_offsets, "-")

  # Reconstruct full matrices
  batch1_normalized <- batch1_matrix
  batch1_normalized[, common_proteins] <- batch1_normalized_common
  batch2_normalized <- batch2_matrix
  batch2_normalized[, common_proteins] <- batch2_normalized_common

  # NO log2 transformation - NPX is already on log2 scale!

  return(list(
    batch1_normalized = batch1_normalized,
    batch2_normalized = batch2_normalized,
    batch1_normalized_raw = batch1_normalized,
    batch2_normalized_raw = batch2_normalized,
    common_proteins = common_proteins,
    method = "median"
  ))
}

# Function to perform PCA on combined matrices and return results
perform_pca_analysis <- function(batch1_matrix, batch2_matrix, common_proteins, batch_labels = c("batch_01", "batch_02"),
                                 batch1_bridge_samples = NULL, batch2_bridge_samples = NULL) {
  # Combine matrices with batch labels
  batch1_sample_ids <- rownames(batch1_matrix)
  batch2_sample_ids <- rownames(batch2_matrix)

  combined_matrix <- rbind(
    batch1_matrix[, common_proteins, drop = FALSE],
    batch2_matrix[, common_proteins, drop = FALSE]
  )
  batch_vector <- c(rep(batch_labels[1], nrow(batch1_matrix)),
                    rep(batch_labels[2], nrow(batch2_matrix)))

  # Mark bridge samples
  all_sample_ids <- c(batch1_sample_ids, batch2_sample_ids)
  is_bridge <- rep(FALSE, length(all_sample_ids))
  if(!is.null(batch1_bridge_samples)) {
    is_bridge[1:length(batch1_sample_ids)][batch1_sample_ids %in% batch1_bridge_samples] <- TRUE
  }
  if(!is.null(batch2_bridge_samples)) {
    is_bridge[(length(batch1_sample_ids)+1):length(all_sample_ids)][batch2_sample_ids %in% batch2_bridge_samples] <- TRUE
  }

  # Remove samples with too many NAs
  sample_completeness <- rowSums(!is.na(combined_matrix)) / ncol(combined_matrix)
  keep_samples <- sample_completeness > 0.5
  combined_matrix <- combined_matrix[keep_samples, ]
  batch_vector <- batch_vector[keep_samples]
  is_bridge <- is_bridge[keep_samples]

  # Impute remaining NAs with column medians
  for(j in 1:ncol(combined_matrix)) {
    combined_matrix[is.na(combined_matrix[, j]), j] <- median(combined_matrix[, j], na.rm = TRUE)
  }

  # Perform PCA
  pca_result <- prcomp(combined_matrix, center = TRUE, scale. = TRUE)
  pca_scores <- pca_result$x[, 1:min(4, ncol(pca_result$x)), drop = FALSE]

  # Calculate variance explained
  var_explained <- summary(pca_result)$importance[2, 1:min(4, ncol(pca_result$x))]

  # Prepare plot data
  plot_data <- data.table(
    PC1 = pca_scores[, 1],
    PC2 = pca_scores[, 2],
    Batch = batch_vector,
    IsBridge = is_bridge
  )

  if(ncol(pca_scores) >= 4) {
    plot_data[, PC3 := pca_scores[, 3]]
    plot_data[, PC4 := pca_scores[, 4]]
  }

  return(list(
    plot_data = plot_data,
    var_explained = var_explained
  ))
}

# Function to create PCA plots (before/after normalization)
create_pca_plots <- function(batch1_matrix, batch2_matrix, common_proteins, batch_labels = c("batch_01", "batch_02"),
                             title_suffix = "", stage = "before", batch1_bridge_samples = NULL, batch2_bridge_samples = NULL) {
  log_info("Creating PCA plots ({stage} normalization)")

  # Perform PCA analysis with bridge sample information
  pca_result <- perform_pca_analysis(batch1_matrix, batch2_matrix, common_proteins, batch_labels,
                                     batch1_bridge_samples, batch2_bridge_samples)
  plot_data <- pca_result$plot_data
  var_explained <- pca_result$var_explained

  # Colourblind-friendly palette matching sex comparison plots
  batch_colors <- c("batch_01" = "#E64B35", "batch_02" = "#4DBBD5")  # Match 08_covariate_sex_comparison colors

  # Create PC1 vs PC2 plot with increased alpha and bridge samples as black circles
  p1 <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Batch)) +
    # Regular samples (colored by batch)
    geom_point(data = plot_data[IsBridge == FALSE], alpha = 0.7, size = 1.5, shape = 16) +
    # Bridge samples (black circles, on top)
    geom_point(data = plot_data[IsBridge == TRUE], color = "black", alpha = 0.8, size = 2.5, shape = 21, fill = "black", stroke = 0.5) +
    scale_color_manual(values = batch_colors,
                       labels = c("batch_01" = "Batch 1", "batch_02" = "Batch 2")) +
    labs(title = paste0("PCA: PC1 vs PC2 (", stage, " normalisation)"),
         subtitle = paste0("Variance explained: PC1=", round(var_explained[1]*100, 1), "%, PC2=", round(var_explained[2]*100, 1), "%"),
         x = paste0("PC1 (", round(var_explained[1]*100, 1), "%)"),
         y = paste0("PC2 (", round(var_explained[2]*100, 1), "%)"),
         color = "Batch") +
    theme_bw(base_size = 11) +
    theme(legend.position = "right",
          panel.grid.minor = element_blank())

  # Create PC3 vs PC4 plot if available
  p2 <- NULL
  if("PC3" %in% names(plot_data) && "PC4" %in% names(plot_data)) {
    p2 <- ggplot(plot_data, aes(x = PC3, y = PC4, color = Batch)) +
      # Regular samples (colored by batch)
      geom_point(data = plot_data[IsBridge == FALSE], alpha = 0.7, size = 1.5, shape = 16) +
      # Bridge samples (black circles, on top)
      geom_point(data = plot_data[IsBridge == TRUE], color = "black", alpha = 0.8, size = 2.5, shape = 21, fill = "black", stroke = 0.5) +
      scale_color_manual(values = batch_colors,
                         labels = c("batch_01" = "Batch 1", "batch_02" = "Batch 2")) +
      labs(title = paste0("PCA: PC3 vs PC4 (", stage, " normalisation)"),
           subtitle = paste0("Variance explained: PC3=", round(var_explained[3]*100, 1), "%, PC4=", round(var_explained[4]*100, 1), "%"),
           x = paste0("PC3 (", round(var_explained[3]*100, 1), "%)"),
           y = paste0("PC4 (", round(var_explained[4]*100, 1), "%)"),
           color = "Batch") +
      theme_bw(base_size = 11) +
      theme(legend.position = "right",
            panel.grid.minor = element_blank())
  }

  # Combine plots
  if(!is.null(p2)) {
    combined_plot <- ggarrange(p1, p2, ncol = 2, common.legend = TRUE, legend = "right")
  } else {
    combined_plot <- p1
  }

  return(list(
    pc1_pc2 = p1,
    pc3_pc4 = p2,
    combined = combined_plot,
    plot_data = plot_data,
    var_explained = var_explained
  ))
}

# Function to create side-by-side PCA comparison (before vs after normalisation)
create_pca_comparison_plots <- function(pca_before, pca_after, method_name = "Bridge") {
  log_info("Creating side-by-side PCA comparison plots for {method_name} method")

  # Colourblind-friendly palette matching 08_covariate_sex_comparison
  batch_colors <- c("batch_01" = "#E64B35", "batch_02" = "#4DBBD5")

  # ===========================================================================
  # PC1 vs PC2: Before and After side-by-side
  # ===========================================================================
  p_pc12_before <- ggplot(pca_before$plot_data, aes(x = PC1, y = PC2, color = Batch)) +
    # Regular samples (colored by batch)
    geom_point(data = pca_before$plot_data[IsBridge == FALSE], alpha = 0.7, size = 1.5, shape = 16) +
    # Bridge samples (black circles, on top)
    geom_point(data = pca_before$plot_data[IsBridge == TRUE], color = "black", alpha = 0.8, size = 2.5, shape = 21, fill = "black", stroke = 0.5) +
    scale_color_manual(values = batch_colors,
                       labels = c("batch_01" = "Batch 1", "batch_02" = "Batch 2")) +
    labs(title = "PCA: PC1 vs PC2 (before normalisation)",
         subtitle = paste0("Variance explained: PC1=", round(pca_before$var_explained[1]*100, 1), "%, PC2=", round(pca_before$var_explained[2]*100, 1), "%"),
         x = paste0("PC1 (", round(pca_before$var_explained[1]*100, 1), "%)"),
         y = paste0("PC2 (", round(pca_before$var_explained[2]*100, 1), "%)"),
         color = "Batch") +
    theme_bw(base_size = 11) +
    theme(legend.position = "right", panel.grid.minor = element_blank())

  p_pc12_after <- ggplot(pca_after$plot_data, aes(x = PC1, y = PC2, color = Batch)) +
    # Regular samples (colored by batch)
    geom_point(data = pca_after$plot_data[IsBridge == FALSE], alpha = 0.7, size = 1.5, shape = 16) +
    # Bridge samples (black circles, on top)
    geom_point(data = pca_after$plot_data[IsBridge == TRUE], color = "black", alpha = 0.8, size = 2.5, shape = 21, fill = "black", stroke = 0.5) +
    scale_color_manual(values = batch_colors,
                       labels = c("batch_01" = "Batch 1", "batch_02" = "Batch 2")) +
    labs(title = paste0("PCA: PC1 vs PC2 (after ", method_name, " normalisation)"),
         subtitle = paste0("Variance explained: PC1=", round(pca_after$var_explained[1]*100, 1), "%, PC2=", round(pca_after$var_explained[2]*100, 1), "%"),
         x = paste0("PC1 (", round(pca_after$var_explained[1]*100, 1), "%)"),
         y = paste0("PC2 (", round(pca_after$var_explained[2]*100, 1), "%)"),
         color = "Batch") +
    theme_bw(base_size = 11) +
    theme(legend.position = "right", panel.grid.minor = element_blank())

  # Combine PC1 vs PC2 (side-by-side with legends on right)
  pc12_combined <- ggarrange(p_pc12_before, p_pc12_after, ncol = 2, widths = c(1, 1))

  # ===========================================================================
  # PC3 vs PC4: Before and After side-by-side
  # ===========================================================================
  pc34_combined <- NULL
  if("PC3" %in% names(pca_before$plot_data) && "PC4" %in% names(pca_before$plot_data)) {
    p_pc34_before <- ggplot(pca_before$plot_data, aes(x = PC3, y = PC4, color = Batch)) +
      # Regular samples (colored by batch)
      geom_point(data = pca_before$plot_data[IsBridge == FALSE], alpha = 0.7, size = 1.5, shape = 16) +
      # Bridge samples (black circles, on top)
      geom_point(data = pca_before$plot_data[IsBridge == TRUE], color = "black", alpha = 0.8, size = 2.5, shape = 21, fill = "black", stroke = 0.5) +
      scale_color_manual(values = batch_colors,
                         labels = c("batch_01" = "Batch 1", "batch_02" = "Batch 2")) +
      labs(title = "PCA: PC3 vs PC4 (before normalisation)",
           subtitle = paste0("Variance explained: PC3=", round(pca_before$var_explained[3]*100, 1), "%, PC4=", round(pca_before$var_explained[4]*100, 1), "%"),
           x = paste0("PC3 (", round(pca_before$var_explained[3]*100, 1), "%)"),
           y = paste0("PC4 (", round(pca_before$var_explained[4]*100, 1), "%)"),
           color = "Batch") +
      theme_bw(base_size = 11) +
      theme(legend.position = "right", panel.grid.minor = element_blank())

    p_pc34_after <- ggplot(pca_after$plot_data, aes(x = PC3, y = PC4, color = Batch)) +
      # Regular samples (colored by batch)
      geom_point(data = pca_after$plot_data[IsBridge == FALSE], alpha = 0.7, size = 1.5, shape = 16) +
      # Bridge samples (black circles, on top)
      geom_point(data = pca_after$plot_data[IsBridge == TRUE], color = "black", alpha = 0.8, size = 2.5, shape = 21, fill = "black", stroke = 0.5) +
      scale_color_manual(values = batch_colors,
                         labels = c("batch_01" = "Batch 1", "batch_02" = "Batch 2")) +
      labs(title = paste0("PCA: PC3 vs PC4 (after ", method_name, " normalisation)"),
           subtitle = paste0("Variance explained: PC3=", round(pca_after$var_explained[3]*100, 1), "%, PC4=", round(pca_after$var_explained[4]*100, 1), "%"),
           x = paste0("PC3 (", round(pca_after$var_explained[3]*100, 1), "%)"),
           y = paste0("PC4 (", round(pca_after$var_explained[4]*100, 1), "%)"),
           color = "Batch") +
      theme_bw(base_size = 11) +
      theme(legend.position = "right", panel.grid.minor = element_blank())

    # Combine PC3 vs PC4 (side-by-side with legends on right)
    pc34_combined <- ggarrange(p_pc34_before, p_pc34_after, ncol = 2, widths = c(1, 1))
  }

  # ===========================================================================
  # Full combined plot (PC1-2 on top, PC3-4 on bottom)
  # ===========================================================================
  if(!is.null(pc34_combined)) {
    full_combined <- ggarrange(pc12_combined, pc34_combined, nrow = 2)
  } else {
    full_combined <- pc12_combined
  }

  return(list(
    pc12_comparison = pc12_combined,
    pc34_comparison = pc34_combined,
    full_combined = full_combined
  ))
}

# Function to create bridge sample scatter plots
create_bridge_scatter_plots <- function(batch1_matrix, batch2_matrix, batch1_bridge_samples, batch2_bridge_samples,
                                       batch1_sample_mapping, batch2_sample_mapping, common_proteins,
                                       batch1_normalized = NULL, batch2_normalized = NULL) {
  log_info("Creating bridge sample scatter plots")

  # Match bridge samples by FINNGENID
  batch1_bridge_mapping <- batch1_sample_mapping[SampleID %in% batch1_bridge_samples]
  batch2_bridge_mapping <- batch2_sample_mapping[SampleID %in% batch2_bridge_samples]

  bridge_pairs <- merge(
    batch1_bridge_mapping[, .(SampleID_batch1 = SampleID, FINNGENID)],
    batch2_bridge_mapping[, .(SampleID_batch2 = SampleID, FINNGENID)],
    by = "FINNGENID"
  )

  if(nrow(bridge_pairs) == 0) {
    log_warn("No matching bridge samples found for scatter plots")
    return(NULL)
  }

  # Calculate mean NPX per bridge sample (across common proteins) - BEFORE
  batch1_bridge_means_before <- rowMeans(batch1_matrix[bridge_pairs$SampleID_batch1, common_proteins, drop = FALSE], na.rm = TRUE)
  batch2_bridge_means_before <- rowMeans(batch2_matrix[bridge_pairs$SampleID_batch2, common_proteins, drop = FALSE], na.rm = TRUE)

  # Calculate correlation before
  cor_before <- cor(batch1_bridge_means_before, batch2_bridge_means_before, use = "complete.obs")

  # Create scatter plot BEFORE normalization
  scatter_before <- ggplot(data.table(
    Batch1 = batch1_bridge_means_before,
    Batch2 = batch2_bridge_means_before
  ), aes(x = Batch1, y = Batch2)) +
    geom_point(alpha = 0.7, color = "#E64B35") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    geom_smooth(method = "lm", se = TRUE, color = "blue", alpha = 0.3) +
    labs(title = "Bridge Samples: Batch 1 vs Batch 2 (Before Normalization)",
         subtitle = paste0("Correlation: R = ", round(cor_before, 3), " | n = ", length(batch1_bridge_means_before), " bridge samples"),
         x = "Batch 1 Mean NPX", y = "Batch 2 Mean NPX") +
    theme_bw()

  # Calculate mean NPX per bridge sample - AFTER (if normalized matrices provided)
  scatter_after <- NULL
  if(!is.null(batch1_normalized) && !is.null(batch2_normalized)) {
    batch1_bridge_means_after <- rowMeans(batch1_normalized[bridge_pairs$SampleID_batch1, common_proteins, drop = FALSE], na.rm = TRUE)
    batch2_bridge_means_after <- rowMeans(batch2_normalized[bridge_pairs$SampleID_batch2, common_proteins, drop = FALSE], na.rm = TRUE)

    # Calculate correlation after
    cor_after <- cor(batch1_bridge_means_after, batch2_bridge_means_after, use = "complete.obs")

    # Create scatter plot AFTER normalization
    scatter_after <- ggplot(data.table(
      Batch1 = batch1_bridge_means_after,
      Batch2 = batch2_bridge_means_after
    ), aes(x = Batch1, y = Batch2)) +
      geom_point(alpha = 0.7, color = "#4DBBD5") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
      geom_smooth(method = "lm", se = TRUE, color = "blue", alpha = 0.3) +
      labs(title = "Bridge Samples: Batch 1 vs Batch 2 (After Normalization)",
           subtitle = paste0("Correlation: R = ", round(cor_after, 3), " | n = ", length(batch1_bridge_means_after), " bridge samples"),
           x = "Batch 1 Mean NPX", y = "Batch 2 Mean NPX") +
      theme_bw()

    # Combine plots
    combined_scatter <- ggarrange(scatter_before, scatter_after, ncol = 2)
  } else {
    combined_scatter <- scatter_before
  }

  return(list(
    before = scatter_before,
    after = scatter_after,
    combined = combined_scatter,
    cor_before = cor_before,
    cor_after = if(!is.null(batch1_normalized)) cor_after else NA
  ))
}

# Function to create paired boxplot of bridging samples before and after normalisation
# Style similar to 06_normalization_effect_median plot
create_bridge_paired_boxplot <- function(batch1_matrix, batch2_matrix,
                                          batch1_normalized, batch2_normalized,
                                          batch1_bridge_samples, batch2_bridge_samples,
                                          batch1_sample_mapping, batch2_sample_mapping,
                                          common_proteins, n_samples_show = 100) {
  log_info("Creating paired boxplot for bridging samples")

  # Match bridge samples by FINNGENID
  batch1_bridge_mapping <- batch1_sample_mapping[SampleID %in% batch1_bridge_samples]
  batch2_bridge_mapping <- batch2_sample_mapping[SampleID %in% batch2_bridge_samples]

  bridge_pairs <- merge(
    batch1_bridge_mapping[, .(SampleID_batch1 = SampleID, FINNGENID)],
    batch2_bridge_mapping[, .(SampleID_batch2 = SampleID, FINNGENID)],
    by = "FINNGENID"
  )

  if(nrow(bridge_pairs) == 0) {
    log_warn("No matching bridge samples found for paired boxplot")
    return(NULL)
  }

  n_bridges <- nrow(bridge_pairs)
  log_info("Found {n_bridges} paired bridge samples for boxplot")

  # Sample proteins for plotting
  set.seed(456)
  sample_proteins <- sample(common_proteins, min(100, length(common_proteins)))

  # Collect data for each bridge sample pair (mean across sampled proteins)
  # BEFORE normalisation
  batch1_bridge_means_before <- rowMeans(batch1_matrix[bridge_pairs$SampleID_batch1, sample_proteins, drop = FALSE], na.rm = TRUE)
  batch2_bridge_means_before <- rowMeans(batch2_matrix[bridge_pairs$SampleID_batch2, sample_proteins, drop = FALSE], na.rm = TRUE)

  # AFTER normalisation
  batch1_bridge_means_after <- rowMeans(batch1_normalized[bridge_pairs$SampleID_batch1, sample_proteins, drop = FALSE], na.rm = TRUE)
  batch2_bridge_means_after <- rowMeans(batch2_normalized[bridge_pairs$SampleID_batch2, sample_proteins, drop = FALSE], na.rm = TRUE)

  # Create data for boxplot - similar to 06_normalization_effect style
  # Each bridge pair becomes a "sample" that we show before and after
  plot_data <- rbind(
    data.table(
      SampleIndex = 1:n_bridges,
      FINNGENID = bridge_pairs$FINNGENID,
      Batch1_NPX = batch1_bridge_means_before,
      Batch2_NPX = batch2_bridge_means_before,
      Stage = "Before"
    ),
    data.table(
      SampleIndex = 1:n_bridges,
      FINNGENID = bridge_pairs$FINNGENID,
      Batch1_NPX = batch1_bridge_means_after,
      Batch2_NPX = batch2_bridge_means_after,
      Stage = "After"
    )
  )

  # Calculate difference between batches (for each bridge sample)
  plot_data[, Batch_Diff := abs(Batch1_NPX - Batch2_NPX)]

  # Set factor levels
  plot_data[, Stage := factor(Stage, levels = c("Before", "After"))]

  # Melt for paired plot
  plot_data_long <- melt(plot_data,
                         id.vars = c("SampleIndex", "FINNGENID", "Stage", "Batch_Diff"),
                         measure.vars = c("Batch1_NPX", "Batch2_NPX"),
                         variable.name = "Batch", value.name = "MeanNPX")
  plot_data_long[, Batch := fifelse(Batch == "Batch1_NPX", "Batch 1", "Batch 2")]

  # Calculate summary statistics
  stats_before <- plot_data[Stage == "Before", .(
    mean_diff = mean(Batch_Diff, na.rm = TRUE),
    median_diff = median(Batch_Diff, na.rm = TRUE),
    sd_diff = sd(Batch_Diff, na.rm = TRUE)
  )]
  stats_after <- plot_data[Stage == "After", .(
    mean_diff = mean(Batch_Diff, na.rm = TRUE),
    median_diff = median(Batch_Diff, na.rm = TRUE),
    sd_diff = sd(Batch_Diff, na.rm = TRUE)
  )]

  log_info("Bridge sample batch difference - Before: mean={round(stats_before$mean_diff, 3)}, median={round(stats_before$median_diff, 3)}")
  log_info("Bridge sample batch difference - After: mean={round(stats_after$mean_diff, 3)}, median={round(stats_after$median_diff, 3)}")

  # Calculate percent reduction in batch difference
  diff_reduction_pct <- (stats_before$mean_diff - stats_after$mean_diff) / stats_before$mean_diff * 100
  log_info("Batch difference reduction: {round(diff_reduction_pct, 1)}%")

  # Colour palette
  batch_colors <- c("Batch 1" = "#56B4E9", "Batch 2" = "#E69F00")

  # ===========================================================================
  # PANEL 1: Paired boxplot showing NPX values by batch, before and after
  # ===========================================================================
  p_boxplot <- ggboxplot(plot_data_long, x = "Batch", y = "MeanNPX", fill = "Batch",
                         palette = batch_colors, add = "jitter",
                         add.params = list(alpha = 0.4, size = 1.5)) +
    facet_wrap(~ Stage, ncol = 2) +
    stat_compare_means(method = "t.test", paired = FALSE,
                       label = "p.format", label.x.npc = "center", label.y.npc = 0.95) +
    labs(title = "Bridge Sample NPX by Batch: Before vs After Normalisation",
         subtitle = paste0("n = ", n_bridges, " bridge samples (paired by FINNGENID) | ",
                          "Mean diff reduction: ", round(diff_reduction_pct, 1), "%"),
         x = "Batch", y = "Mean NPX (across 100 random proteins)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "gray90"),
          strip.text = element_text(face = "bold", size = 11))

  # ===========================================================================
  # PANEL 2: Paired line plot showing each bridge sample's values in both batches
  # ===========================================================================
  # Select subset of samples to show (avoid overcrowding)
  show_samples <- min(n_samples_show, n_bridges)
  set.seed(789)
  sample_indices <- sample(1:n_bridges, show_samples)

  plot_data_subset <- plot_data_long[SampleIndex %in% sample_indices]

  p_paired <- ggplot(plot_data_subset, aes(x = Batch, y = MeanNPX, group = SampleIndex)) +
    geom_line(alpha = 0.3, color = "gray50") +
    geom_point(aes(color = Batch), size = 2, alpha = 0.7) +
    scale_color_manual(values = batch_colors) +
    facet_wrap(~ Stage, ncol = 2) +
    labs(title = paste0("Paired Bridge Samples: Batch 1 vs Batch 2 (showing ", show_samples, " of ", n_bridges, " samples)"),
         subtitle = "Lines connect the same FINNGENID across batches",
         x = "Batch", y = "Mean NPX (across 100 random proteins)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "gray90"),
          strip.text = element_text(face = "bold", size = 11))

  # ===========================================================================
  # PANEL 3: Distribution of batch differences (smaller = better harmonisation)
  # ===========================================================================
  p_diff <- ggplot(plot_data, aes(x = Stage, y = Batch_Diff, fill = Stage)) +
    geom_violin(alpha = 0.6, width = 0.8) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
    scale_fill_manual(values = c("Before" = "#F8766D", "After" = "#00BA38")) +
    stat_compare_means(method = "t.test", paired = TRUE,
                       label = "p.format", label.x = 1.5, label.y.npc = 0.95) +
    labs(title = "Absolute Batch Difference per Bridge Sample",
         subtitle = paste0("Before: mean=", round(stats_before$mean_diff, 3),
                          " | After: mean=", round(stats_after$mean_diff, 3),
                          " | Reduction: ", round(diff_reduction_pct, 1), "%"),
         x = "Normalisation Stage", y = "|Batch 1 NPX - Batch 2 NPX|") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")

  # ===========================================================================
  # Combine into multi-panel figure
  # ===========================================================================
  combined_plot <- ggarrange(
    p_boxplot, p_paired, p_diff,
    nrow = 3, heights = c(1, 1, 0.8),
    labels = c("A", "B", "C")
  )

  combined_plot <- annotate_figure(
    combined_plot,
    top = text_grob("Bridge Sample Harmonisation Assessment", face = "bold", size = 14)
  )

  return(list(
    boxplot = p_boxplot,
    paired = p_paired,
    difference = p_diff,
    combined = combined_plot,
    stats = list(
      before = stats_before,
      after = stats_after,
      reduction_pct = diff_reduction_pct
    )
  ))
}

# Function to create faceted NPX distribution plots with histogram+density and violin plots
create_distribution_plots <- function(batch1_matrix, batch2_matrix, batch1_normalized, batch2_normalized,
                                     common_proteins, method_name = "Bridge") {
  log_info("Creating enhanced NPX distribution plots for {method_name} method")

  # Sample proteins for plotting (100 random proteins)
  set.seed(123)
  sample_proteins <- sample(common_proteins, min(100, length(common_proteins)))

  # Prepare data: combine batch1 and batch2, before and after
  plot_data <- rbind(
    data.table(
      NPX = as.vector(batch1_matrix[, sample_proteins, drop = FALSE]),
      Batch = "Batch 1",
      Stage = "Before"
    ),
    data.table(
      NPX = as.vector(batch1_normalized[, sample_proteins, drop = FALSE]),
      Batch = "Batch 1",
      Stage = "After"
    ),
    data.table(
      NPX = as.vector(batch2_matrix[, sample_proteins, drop = FALSE]),
      Batch = "Batch 2",
      Stage = "Before"
    ),
    data.table(
      NPX = as.vector(batch2_normalized[, sample_proteins, drop = FALSE]),
      Batch = "Batch 2",
      Stage = "After"
    )
  )

  # Remove NAs
  plot_data <- plot_data[!is.na(NPX)]

  # Set factor levels: Before on left, After on right
  plot_data[, Stage := factor(Stage, levels = c("Before", "After"))]

  # Define colour palette (similar to example: cyan and gold)
  batch_colors <- c("Batch 1" = "#56B4E9", "Batch 2" = "#E69F00")  # Colourblind-friendly

  # Calculate medians for vertical lines
  medians <- plot_data[, .(median_npx = median(NPX, na.rm = TRUE)), by = .(Batch, Stage)]

  # ===========================================================================
  # PANEL 1: Histogram + Density plots (Before on left, After on right)
  # ===========================================================================
  p_hist <- ggplot(plot_data, aes(x = NPX, fill = Batch, color = Batch)) +
    geom_histogram(aes(y = after_stat(density)), alpha = 0.5, position = "identity", bins = 40) +
    geom_density(alpha = 0.3, linewidth = 1) +
    geom_vline(data = medians, aes(xintercept = median_npx, color = Batch),
               linetype = "dashed", linewidth = 0.8) +
    facet_wrap(~ Stage, ncol = 2) +
    scale_fill_manual(values = batch_colors) +
    scale_color_manual(values = batch_colors) +
    labs(title = paste0("NPX Distribution: ", method_name, " Normalisation"),
         subtitle = "Histogram with density overlay. Dashed lines = batch medians",
         x = "NPX Value", y = "Density") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank()
    )

  # ===========================================================================
  # PANEL 2: Violin plots with statistical test
  # ===========================================================================
  n_batches <- length(unique(plot_data$Batch))

  # Use t-test for 2 batches, Kruskal-Wallis for >2 batches
  stat_method <- ifelse(n_batches == 2, "t.test", "kruskal.test")
  stat_label <- ifelse(n_batches == 2, "T-test", "Kruskal-Wallis")

  p_violin <- ggviolin(plot_data, x = "Batch", y = "NPX", fill = "Batch",
                       palette = batch_colors,
                       add = "boxplot", add.params = list(fill = "white", width = 0.15)) +
    facet_wrap(~ Stage, ncol = 2) +
    stat_compare_means(method = stat_method, label = "p.format",
                       label.x.npc = "center", label.y.npc = 0.95,
                       size = 3.5) +
    labs(title = paste0("Batch Comparison: ", method_name, " Normalisation"),
         subtitle = paste0("Violin plots with ", stat_label, " p-value"),
         x = "Batch", y = "NPX Value") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank()
    )

  # ===========================================================================
  # Combine into a single figure (histogram/density on top, violin on bottom)
  # ===========================================================================
  combined_plot <- ggarrange(
    p_hist, p_violin,
    nrow = 2, heights = c(1.2, 1),
    labels = c("A", "B"),
    common.legend = FALSE
  )

  # Add overall title
  combined_plot <- annotate_figure(
    combined_plot,
    top = text_grob(
      paste0("NPX Distribution Before and After ", method_name, " Normalisation"),
      face = "bold", size = 14
    )
  )

  return(combined_plot)
}

# Updated evaluation function for cross-batch normalization
# NOTE: Uses SD/MAD/IQR instead of CV because CV is unchanged by multiplicative normalization
# CV = SD/mean, and both scale proportionally with multiplicative normalization
evaluate_cross_batch_normalization <- function(batch1_raw, batch2_raw, batch1_normalized, batch2_normalized,
                                               common_proteins, method_name) {
  log_info("╔══════════════════════════════════════════════════════════════════╗")
  log_info("║ EVALUATING {toupper(method_name)} NORMALIZATION                            ║")
  log_info("╚══════════════════════════════════════════════════════════════════╝")

  # Calculate SD before normalization (common proteins only)
  sd_batch1_before <- apply(batch1_raw[, common_proteins, drop = FALSE], 2, sd, na.rm = TRUE)
  sd_batch2_before <- apply(batch2_raw[, common_proteins, drop = FALSE], 2, sd, na.rm = TRUE)

  # Calculate SD after normalization
  sd_batch1_after <- apply(batch1_normalized[, common_proteins, drop = FALSE], 2, sd, na.rm = TRUE)
  sd_batch2_after <- apply(batch2_normalized[, common_proteins, drop = FALSE], 2, sd, na.rm = TRUE)

  # Calculate MAD before/after
  mad_batch1_before <- apply(batch1_raw[, common_proteins, drop = FALSE], 2, mad, na.rm = TRUE)
  mad_batch2_before <- apply(batch2_raw[, common_proteins, drop = FALSE], 2, mad, na.rm = TRUE)
  mad_batch1_after <- apply(batch1_normalized[, common_proteins, drop = FALSE], 2, mad, na.rm = TRUE)
  mad_batch2_after <- apply(batch2_normalized[, common_proteins, drop = FALSE], 2, mad, na.rm = TRUE)

  # Calculate IQR before/after
  iqr_batch1_before <- apply(batch1_raw[, common_proteins, drop = FALSE], 2, IQR, na.rm = TRUE)
  iqr_batch2_before <- apply(batch2_raw[, common_proteins, drop = FALSE], 2, IQR, na.rm = TRUE)
  iqr_batch1_after <- apply(batch1_normalized[, common_proteins, drop = FALSE], 2, IQR, na.rm = TRUE)
  iqr_batch2_after <- apply(batch2_normalized[, common_proteins, drop = FALSE], 2, IQR, na.rm = TRUE)

  # Calculate percentage reductions for SD
  mean_sd_batch1_before <- mean(sd_batch1_before, na.rm = TRUE)
  mean_sd_batch1_after <- mean(sd_batch1_after, na.rm = TRUE)
  mean_sd_batch2_before <- mean(sd_batch2_before, na.rm = TRUE)
  mean_sd_batch2_after <- mean(sd_batch2_after, na.rm = TRUE)

  sd_reduction_batch1 <- ifelse(mean_sd_batch1_before > 0,
    (mean_sd_batch1_before - mean_sd_batch1_after) / mean_sd_batch1_before * 100, NA_real_)
  sd_reduction_batch2 <- ifelse(mean_sd_batch2_before > 0,
    (mean_sd_batch2_before - mean_sd_batch2_after) / mean_sd_batch2_before * 100, NA_real_)

  # Calculate percentage reductions for MAD
  mean_mad_batch1_before <- mean(mad_batch1_before, na.rm = TRUE)
  mean_mad_batch1_after <- mean(mad_batch1_after, na.rm = TRUE)
  mean_mad_batch2_before <- mean(mad_batch2_before, na.rm = TRUE)
  mean_mad_batch2_after <- mean(mad_batch2_after, na.rm = TRUE)

  mad_reduction_batch1 <- ifelse(mean_mad_batch1_before > 0,
    (mean_mad_batch1_before - mean_mad_batch1_after) / mean_mad_batch1_before * 100, NA_real_)
  mad_reduction_batch2 <- ifelse(mean_mad_batch2_before > 0,
    (mean_mad_batch2_before - mean_mad_batch2_after) / mean_mad_batch2_before * 100, NA_real_)

  # Calculate percentage reductions for IQR
  mean_iqr_batch1_before <- mean(iqr_batch1_before, na.rm = TRUE)
  mean_iqr_batch1_after <- mean(iqr_batch1_after, na.rm = TRUE)
  mean_iqr_batch2_before <- mean(iqr_batch2_before, na.rm = TRUE)
  mean_iqr_batch2_after <- mean(iqr_batch2_after, na.rm = TRUE)

  iqr_reduction_batch1 <- ifelse(mean_iqr_batch1_before > 0,
    (mean_iqr_batch1_before - mean_iqr_batch1_after) / mean_iqr_batch1_before * 100, NA_real_)
  iqr_reduction_batch2 <- ifelse(mean_iqr_batch2_before > 0,
    (mean_iqr_batch2_before - mean_iqr_batch2_after) / mean_iqr_batch2_before * 100, NA_real_)

  # Summary statistics table (SD/MAD/IQR metrics, not CV)
  eval_stats <- data.table(
    method = method_name,
    batch = c("batch_01", "batch_02"),
    mean_sd_before = c(mean_sd_batch1_before, mean_sd_batch2_before),
    mean_sd_after = c(mean_sd_batch1_after, mean_sd_batch2_after),
    sd_reduction_pct = c(sd_reduction_batch1, sd_reduction_batch2),
    mean_mad_before = c(mean_mad_batch1_before, mean_mad_batch2_before),
    mean_mad_after = c(mean_mad_batch1_after, mean_mad_batch2_after),
    mad_reduction_pct = c(mad_reduction_batch1, mad_reduction_batch2),
    mean_iqr_before = c(mean_iqr_batch1_before, mean_iqr_batch2_before),
    mean_iqr_after = c(mean_iqr_batch1_after, mean_iqr_batch2_after),
    iqr_reduction_pct = c(iqr_reduction_batch1, iqr_reduction_batch2)
  )

  # Calculate overall improvement using SD (primary metric)
  overall_sd_improvement <- mean(c(sd_reduction_batch1, sd_reduction_batch2), na.rm = TRUE)

  log_info("Batch 1 results:")
  log_info("  SD: {round(mean_sd_batch1_before, 3)} → {round(mean_sd_batch1_after, 3)} ({round(sd_reduction_batch1, 2)}% reduction)")
  log_info("  MAD: {round(mean_mad_batch1_before, 3)} → {round(mean_mad_batch1_after, 3)} ({round(mad_reduction_batch1, 2)}% reduction)")
  log_info("  IQR: {round(mean_iqr_batch1_before, 3)} → {round(mean_iqr_batch1_after, 3)} ({round(iqr_reduction_batch1, 2)}% reduction)")
  log_info("Batch 2 results:")
  log_info("  SD: {round(mean_sd_batch2_before, 3)} → {round(mean_sd_batch2_after, 3)} ({round(sd_reduction_batch2, 2)}% reduction)")
  log_info("  MAD: {round(mean_mad_batch2_before, 3)} → {round(mean_mad_batch2_after, 3)} ({round(mad_reduction_batch2, 2)}% reduction)")
  log_info("  IQR: {round(mean_iqr_batch2_before, 3)} → {round(mean_iqr_batch2_after, 3)} ({round(iqr_reduction_batch2, 2)}% reduction)")
  log_info("Overall SD improvement: {round(overall_sd_improvement, 2)}%")
  log_info("════════════════════════════════════════════════════════════════════")

  return(list(
    eval_stats = eval_stats,
    overall_improvement = overall_sd_improvement  # Use SD reduction as primary metric
  ))
}

# Main execution
main <- function() {

  # Check if multi-batch mode is enabled
  multi_batch_mode <- tryCatch(
    isTRUE(config$parameters$normalization$multi_batch_mode),
    error = function(e) FALSE
  )

  # Check if enhanced bridge normalisation should be run
  run_enhanced_bridge <- tryCatch(
    isTRUE(config$parameters$normalization$run_enhanced_bridge),
    error = function(e) FALSE
  )

  # CRITICAL: Step 07 is EXCLUSIVELY for cross-batch normalization
  # It should be skipped in single-batch mode (step 06 handles within-batch normalization)
  if (!multi_batch_mode) {
    # Set environment variable to indicate step was skipped
    Sys.setenv(PIPELINE_STEP_SKIPPED = "TRUE")
    log_info("Single-batch mode detected: Skipping step 07 (cross-batch bridge normalization)")
    log_info("Within-batch normalization is handled in step 06 (06_normalize_data.R)")
    log_info("Step 07 is exclusively for cross-batch data integration")
    cat("\n=== ENHANCED BRIDGE NORMALIZATION SKIPPED ===\n")
    cat("Reason: Single-batch analysis (multi_batch_mode = false)\n")
    cat("Within-batch normalization: Step 06 (06_normalize_data.R)\n")
    cat("Cross-batch normalization: Step 07 (07_bridge_normalization.R) - requires multi_batch_mode = true\n\n")
    return(NULL)
  }

  # Auto-enable in multi-batch mode if not explicitly disabled
  if (multi_batch_mode && !run_enhanced_bridge) {
    log_info("Multi-batch mode detected: Auto-enabling enhanced bridge normalization")
    run_enhanced_bridge <- TRUE
  }

  if(!run_enhanced_bridge) {
    # Set environment variable to indicate step was skipped
    Sys.setenv(PIPELINE_STEP_SKIPPED = "TRUE")
    log_info("Enhanced bridge normalization is disabled in config (run_enhanced_bridge: false)")
    log_info("Skipping step 07 - Enhanced bridge normalization")
    log_info("This step is only needed for multi-batch integration (FG2+FG3 or Batch1+Batch2)")
    cat("\n=== ENHANCED BRIDGE NORMALIZATION SKIPPED ===\n")
    cat("Reason: Single-batch analysis (config: run_enhanced_bridge = false)\n")
    cat("To enable: Set run_enhanced_bridge: true in config.yaml\n")
    cat("Use case: Multi-batch data integration only\n\n")
    return(NULL)
  }

  # CRITICAL: Get all batch IDs for cross-batch normalization
  if (is.null(config$batches) || !is.list(config$batches)) {
    log_error("config$batches is not available or not a list - cannot perform cross-batch normalization")
    stop("Multi-batch mode requires batches configuration")
  }

  all_batch_ids <- tryCatch(names(config$batches), error = function(e) NULL)
  if (is.null(all_batch_ids) || length(all_batch_ids) < 2) {
    log_error("Multi-batch mode requires at least 2 batches, found: {length(all_batch_ids) %||% 0}")
    stop("Insufficient batches for cross-batch normalization")
  }

  log_info("Cross-batch normalization: Processing {length(all_batch_ids)} batches")
  log_info("  Batches: {paste(all_batch_ids, collapse=', ')}")

  # Load data from ALL batches for cross-batch normalization
  log_info("Loading data from all batches for cross-batch normalization")

  # Determine input based on configuration
  input_choice <- tryCatch(config$parameters$normalization$enhanced_bridge_input, error = function(e) "qc_passed")

  # Load matrices, sample mappings, and bridge samples for all batches
  batch_matrices <- list()
  batch_sample_mappings <- list()
  batch_bridge_finngenids <- list()

  for (current_batch_id in all_batch_ids) {
    log_info("Loading data for batch: {current_batch_id}")

    # Load matrix
    if(input_choice == "qc_passed") {
      matrix_path <- get_output_path("05d", "05d_npx_matrix_all_qc_passed", current_batch_id, "phenotypes", "rds", config = config)
      if (!file.exists(matrix_path)) {
        log_warn("QC-passed matrix not found for {current_batch_id}: {matrix_path}")
        log_warn("Falling back to analysis-ready matrix from step 00")
        matrix_path <- get_output_path("00", "npx_matrix_analysis_ready", current_batch_id, "qc", config = config)
      }
    } else if(input_choice == "raw") {
      matrix_path <- get_output_path("00", "npx_matrix_analysis_ready", current_batch_id, "qc", config = config)
    } else {
      log_error("Invalid enhanced_bridge_input: {input_choice}")
      stop("Invalid configuration: enhanced_bridge_input must be 'qc_passed' or 'raw'")
    }

    if (!file.exists(matrix_path)) {
      log_error("NPX matrix file not found for {current_batch_id}: {matrix_path}")
      stop(paste0("NPX matrix file not found for ", current_batch_id, ": ", matrix_path))
    }
    batch_matrices[[current_batch_id]] <- readRDS(matrix_path)
    log_info("  Loaded matrix: {nrow(batch_matrices[[current_batch_id]])} samples × {ncol(batch_matrices[[current_batch_id]])} proteins")

    # Load sample mapping (needed for FINNGENID matching)
    sample_mapping_path <- get_output_path("00", "sample_mapping", current_batch_id, "qc", config = config)
    if (!file.exists(sample_mapping_path)) {
      log_error("Sample mapping file not found for {current_batch_id}: {sample_mapping_path}")
      stop(paste0("Sample mapping file not found for ", current_batch_id, ": ", sample_mapping_path))
    }
    batch_sample_mappings[[current_batch_id]] <- readRDS(sample_mapping_path)

    # Load bridge samples and extract FINNGENIDs
    bridge_result_path <- get_output_path("00", "bridge_samples_identified", current_batch_id, "normalized", config = config)
    if (!file.exists(bridge_result_path)) {
      log_error("Bridge result file not found for {current_batch_id}: {bridge_result_path}")
      stop(paste0("Bridge result file not found for ", current_batch_id, ": ", bridge_result_path))
    }
    bridge_result <- readRDS(bridge_result_path)

    # Extract FINNGENIDs from bridge_summary (not just SampleIDs)
    if (!is.null(bridge_result$bridge_summary) && "FINNGENID" %in% names(bridge_result$bridge_summary)) {
      bridge_finngenids <- unique(bridge_result$bridge_summary[!is.na(FINNGENID)]$FINNGENID)
      batch_bridge_finngenids[[current_batch_id]] <- bridge_finngenids
      log_info("  Bridge FINNGENIDs: {length(bridge_finngenids)}")
    } else {
      log_warn("  No bridge FINNGENIDs found for {current_batch_id} - bridge_summary may be missing FINNGENID column")
      batch_bridge_finngenids[[current_batch_id]] <- character(0)
    }
  }

  # CRITICAL: For cross-batch normalization, we need to identify the 24 direct Batch_01↔Batch_02 samples
  # (excluding EA5 samples that are in Batch_00)
  # Load bridging metadata to get the correct filter
  bridging_samples_file <- tryCatch(
    get_batch_input_path("bridging_samples_file", all_batch_ids[1], config),
    error = function(e) NULL
  )

  if (!is.null(bridging_samples_file) && file.exists(bridging_samples_file)) {
    log_info("Loading unified bridging metadata: {bridging_samples_file}")
    bridging_metadata <- fread(bridging_samples_file, sep = "\t")

    # Identify 24 direct Batch_01↔Batch_02 samples (excluding EA5)
    if ("in_Batch_01" %in% names(bridging_metadata) && "in_Batch_02" %in% names(bridging_metadata) && "in_Batch_00" %in% names(bridging_metadata)) {
      batch_01_true <- bridging_metadata$in_Batch_01 == TRUE | bridging_metadata$in_Batch_01 == "TRUE" | bridging_metadata$in_Batch_01 == "T"
      batch_02_true <- bridging_metadata$in_Batch_02 == TRUE | bridging_metadata$in_Batch_02 == "TRUE" | bridging_metadata$in_Batch_02 == "T"
      batch_00_false <- bridging_metadata$in_Batch_00 == FALSE | bridging_metadata$in_Batch_00 == "FALSE" | bridging_metadata$in_Batch_00 == "F"
      direct_b12_finngenids <- bridging_metadata[batch_01_true & batch_02_true & batch_00_false]$FINNGENID
      log_info("Identified {length(direct_b12_finngenids)} direct Batch_01↔Batch_02 bridging FINNGENIDs (excluding EA5)")

      # Use these for cross-batch normalization
      batch_bridge_finngenids[["batch_01"]] <- direct_b12_finngenids
      batch_bridge_finngenids[["batch_02"]] <- direct_b12_finngenids
    }
  }

  # CRITICAL: For cross-batch normalization, we need at least 2 batches
  # Currently, we'll process batch_01 and batch_02 (the two main batches)
  # If more batches exist, we'll process them pairwise
  if (length(all_batch_ids) < 2) {
    log_error("Cross-batch normalization requires at least 2 batches, found: {length(all_batch_ids)}")
    stop("Insufficient batches for cross-batch normalization")
  }

  # For now, process the first two batches (batch_01 and batch_02)
  # TODO: Extend to handle more batches if needed
  batch1_id <- all_batch_ids[1]
  batch2_id <- all_batch_ids[2]

  log_info("Performing cross-batch bridge normalization between {batch1_id} and {batch2_id}")

  # Get bridge FINNGENIDs for both batches
  batch1_bridge_finngenids <- batch_bridge_finngenids[[batch1_id]]
  batch2_bridge_finngenids <- batch_bridge_finngenids[[batch2_id]]

  if (length(batch1_bridge_finngenids) == 0 || length(batch2_bridge_finngenids) == 0) {
    log_error("No bridge FINNGENIDs found for cross-batch normalization")
    log_error("  {batch1_id}: {length(batch1_bridge_finngenids)} bridge FINNGENIDs")
    log_error("  {batch2_id}: {length(batch2_bridge_finngenids)} bridge FINNGENIDs")
    stop("Bridge samples must be identified in both batches for cross-batch normalization")
  }

  # Find common proteins
  common_proteins <- intersect(colnames(batch_matrices[[batch1_id]]), colnames(batch_matrices[[batch2_id]]))
  log_info("Common proteins across batches: {length(common_proteins)}")

  if(length(common_proteins) < 100) {
    log_error("Too few common proteins ({length(common_proteins)}) for cross-batch normalization")
    stop("Insufficient common proteins")
  }

  # Get bridge SampleIDs for visualization (before normalization)
  batch1_bridge_samples <- batch_sample_mappings[[batch1_id]][FINNGENID %in% batch1_bridge_finngenids]$SampleID
  batch2_bridge_samples <- batch_sample_mappings[[batch2_id]][FINNGENID %in% batch2_bridge_finngenids]$SampleID
  batch1_bridge_samples <- intersect(batch1_bridge_samples, rownames(batch_matrices[[batch1_id]]))
  batch2_bridge_samples <- intersect(batch2_bridge_samples, rownames(batch_matrices[[batch2_id]]))

  # Create PCA plots BEFORE normalization
  log_info("Creating PCA plots before normalization")
  pca_before <- create_pca_plots(
    batch_matrices[[batch1_id]],
    batch_matrices[[batch2_id]],
    common_proteins,
    batch_labels = c(batch1_id, batch2_id),
    stage = "before",
    batch1_bridge_samples = batch1_bridge_samples,
    batch2_bridge_samples = batch2_bridge_samples
  )

  # Perform cross-batch bridge normalization (PRIMARY METHOD)
  log_info("Performing cross-batch bridge normalization (bridge method - PRIMARY)")
  cross_batch_result_bridge <- cross_batch_bridge_normalization(
    batch_matrices[[batch1_id]],
    batch_matrices[[batch2_id]],
    batch1_bridge_finngenids,
    batch2_bridge_finngenids,
    batch_sample_mappings[[batch1_id]],
    batch_sample_mappings[[batch2_id]],
    method = "median"  # Bridge normalization uses median method
  )

  if (is.null(cross_batch_result_bridge)) {
    log_error("Cross-batch bridge normalization (bridge) failed")
    stop("Cross-batch bridge normalization failed - cannot proceed")
  }

  # Perform median normalization (COMPARISON METHOD)
  log_info("Performing cross-batch median normalization (global variant - COMPARISON)")
  median_result <- cross_batch_median_normalization(
    batch_matrices[[batch1_id]],
    batch_matrices[[batch2_id]],
    common_proteins,
    reference_batch = "batch_02"
  )

  # Perform ComBat normalization (COMPARISON METHOD)
  log_info("Performing cross-batch ComBat normalization (COMPARISON)")
  combat_result <- cross_batch_bridge_normalization(
    batch_matrices[[batch1_id]],
    batch_matrices[[batch2_id]],
    batch1_bridge_finngenids,
    batch2_bridge_finngenids,
    batch_sample_mappings[[batch1_id]],
    batch_sample_mappings[[batch2_id]],
    method = "combat"
  )

  if (is.null(combat_result)) {
    log_warn("ComBat normalization failed - will continue without ComBat results")
  }

  # Extract normalized matrices for each batch
  batch1_normalized_bridge <- cross_batch_result_bridge$batch1_normalized
  batch2_normalized_bridge <- cross_batch_result_bridge$batch2_normalized
  batch1_normalized_median <- median_result$batch1_normalized
  batch2_normalized_median <- median_result$batch2_normalized

  batch1_normalized_combat <- NULL
  batch2_normalized_combat <- NULL
  if(!is.null(combat_result)) {
    batch1_normalized_combat <- combat_result$batch1_normalized
    batch2_normalized_combat <- combat_result$batch2_normalized
  }

  log_info("Bridge samples available: Batch1={length(batch1_bridge_samples)}, Batch2={length(batch2_bridge_samples)}")

  # Evaluate all normalization methods
  log_info("Evaluating all normalization methods")
  eval_bridge <- evaluate_cross_batch_normalization(
    batch_matrices[[batch1_id]], batch_matrices[[batch2_id]],
    cross_batch_result_bridge$batch1_normalized_raw, cross_batch_result_bridge$batch2_normalized_raw,
    common_proteins, "bridge"
  )

  eval_median <- evaluate_cross_batch_normalization(
    batch_matrices[[batch1_id]], batch_matrices[[batch2_id]],
    median_result$batch1_normalized_raw, median_result$batch2_normalized_raw,
    common_proteins, "median"
  )

  eval_combat <- NULL
  if(!is.null(combat_result)) {
    eval_combat <- evaluate_cross_batch_normalization(
      batch_matrices[[batch1_id]], batch_matrices[[batch2_id]],
      combat_result$batch1_normalized_raw, combat_result$batch2_normalized_raw,
      common_proteins, "combat"
    )
  }

  # Combine all evaluations
  all_evaluations <- rbind(
    eval_bridge$eval_stats,
    eval_median$eval_stats
  )
  if(!is.null(eval_combat)) {
    all_evaluations <- rbind(all_evaluations, eval_combat$eval_stats)
  }

  # Select best method based on SD reduction (average across both batches)
  # NOTE: Using SD reduction instead of CV because CV is unchanged by multiplicative normalization
  bridge_sd_reduction <- eval_bridge$overall_improvement
  median_sd_reduction <- eval_median$overall_improvement
  combat_sd_reduction <- if(!is.null(eval_combat)) eval_combat$overall_improvement else -Inf

  best_method <- "bridge"
  best_sd_reduction <- bridge_sd_reduction

  if (median_sd_reduction > best_sd_reduction) {
    best_method <- "median"
    best_sd_reduction <- median_sd_reduction
  }
  if (combat_sd_reduction > best_sd_reduction) {
    best_method <- "combat"
    best_sd_reduction <- combat_sd_reduction
  }

  log_info("")
  log_info("╔══════════════════════════════════════════════════════════════════╗")
  log_info("║ METHOD COMPARISON (SD REDUCTION %)                               ║")
  log_info("╚══════════════════════════════════════════════════════════════════╝")
  log_info("  Bridge method: SD improvement = {round(bridge_sd_reduction, 2)}%")
  log_info("  Median method: SD improvement = {round(median_sd_reduction, 2)}%")
  if(!is.null(eval_combat)) {
    log_info("  ComBat method: SD improvement = {round(combat_sd_reduction, 2)}%")
  }
  log_info("")
  log_info("  ★ Best method: {toupper(best_method)} (SD improvement: {round(best_sd_reduction, 2)}%)")
  log_info("════════════════════════════════════════════════════════════════════")

  # Create PCA plots AFTER normalization (using bridge method as primary)
  log_info("Creating PCA plots after normalization")
  pca_after_bridge <- create_pca_plots(
    batch1_normalized_bridge, batch2_normalized_bridge,
    common_proteins,
    batch_labels = c(batch1_id, batch2_id),
    stage = "after",
    batch1_bridge_samples = batch1_bridge_samples,
    batch2_bridge_samples = batch2_bridge_samples
  )

  # Create side-by-side PCA comparison plots (before vs after)
  log_info("Creating side-by-side PCA comparison plots")
  pca_comparison <- create_pca_comparison_plots(pca_before, pca_after_bridge, "Bridge")

  # Create bridge sample scatter plots
  log_info("Creating bridge sample scatter plots")
  bridge_scatter <- create_bridge_scatter_plots(
    batch_matrices[[batch1_id]], batch_matrices[[batch2_id]],
    batch1_bridge_samples, batch2_bridge_samples,
    batch_sample_mappings[[batch1_id]], batch_sample_mappings[[batch2_id]],
    common_proteins,
    batch1_normalized_bridge, batch2_normalized_bridge
  )

  # Create paired boxplot for bridge samples (shows harmonisation effectiveness)
  log_info("Creating paired boxplot for bridge samples")
  bridge_paired_boxplot <- create_bridge_paired_boxplot(
    batch_matrices[[batch1_id]], batch_matrices[[batch2_id]],
    batch1_normalized_bridge, batch2_normalized_bridge,
    batch1_bridge_samples, batch2_bridge_samples,
    batch_sample_mappings[[batch1_id]], batch_sample_mappings[[batch2_id]],
    common_proteins
  )

  # Create distribution plots for each method
  log_info("Creating distribution plots")
  dist_plot_bridge <- create_distribution_plots(
    batch_matrices[[batch1_id]], batch_matrices[[batch2_id]],
    batch1_normalized_bridge, batch2_normalized_bridge,
    common_proteins, "Bridge"
  )

  dist_plot_median <- create_distribution_plots(
    batch_matrices[[batch1_id]], batch_matrices[[batch2_id]],
    batch1_normalized_median, batch2_normalized_median,
    common_proteins, "Median"
  )

  dist_plot_combat <- NULL
  if(!is.null(batch1_normalized_combat) && !is.null(batch2_normalized_combat)) {
    dist_plot_combat <- create_distribution_plots(
      batch_matrices[[batch1_id]], batch_matrices[[batch2_id]],
      batch1_normalized_combat, batch2_normalized_combat,
      common_proteins, "ComBat"
    )
  }

  # Save outputs with 07_ prefix using batch-aware paths
  log_info("Saving cross-batch normalization results")

  # Use first batch ID for cross-batch outputs (as per pipeline runner logic)
  output_batch_id <- batch1_id

  # ============================================================================
  # CROSS-BATCH BRIDGE NORMALISATION OUTPUTS
  # ============================================================================
  # New naming convention for clarity:
  #   - 07_npx_matrix_cross_batch_bridge_{batch}.rds: Per-batch normalized matrix
  #   - 07_cross_batch_normalization_result.rds: Combined result with metadata
  # ============================================================================

  log_info("Saving cross-batch bridge normalized matrices (per-batch)")

  # Save individual normalized matrices for EACH batch
  batch1_bridge_path <- get_output_path(step_num, "npx_matrix_cross_batch_bridge", batch1_id, "normalized", config = config)
  batch2_bridge_path <- get_output_path(step_num, "npx_matrix_cross_batch_bridge", batch2_id, "normalized", config = config)
  combined_result_path <- get_output_path(step_num, "cross_batch_normalization_result", output_batch_id, "normalized", config = config)

  ensure_output_dir(batch1_bridge_path)
  ensure_output_dir(batch2_bridge_path)
  ensure_output_dir(combined_result_path)

  # Save batch 1 cross-batch normalized matrix
  saveRDS(batch1_normalized_bridge, batch1_bridge_path)
  log_info("  Saved batch 1 cross-batch bridge normalized: {batch1_bridge_path}")
  log_info("    Dimensions: {nrow(batch1_normalized_bridge)} samples × {ncol(batch1_normalized_bridge)} proteins")

  # Save batch 2 cross-batch normalized matrix
  saveRDS(batch2_normalized_bridge, batch2_bridge_path)
  log_info("  Saved batch 2 cross-batch bridge normalized: {batch2_bridge_path}")
  log_info("    Dimensions: {nrow(batch2_normalized_bridge)} samples × {ncol(batch2_normalized_bridge)} proteins")

  # Save combined normalization result (contains both matrices + metadata)
  normalization_result <- list(
    batch1_normalized = batch1_normalized_bridge,
    batch2_normalized = batch2_normalized_bridge,
    batch1_offsets = cross_batch_result_bridge$batch1_offsets,
    batch2_offsets = cross_batch_result_bridge$batch2_offsets,
    batch1_bridge_samples_used = cross_batch_result_bridge$batch1_bridge_samples_used,
    batch2_bridge_samples_used = cross_batch_result_bridge$batch2_bridge_samples_used,
    common_proteins = common_proteins,
    method = "bridge",
    batch1_id = batch1_id,
    batch2_id = batch2_id
  )

  saveRDS(normalization_result, combined_result_path)
  log_info("  Saved combined result object: {combined_result_path}")

  # Save alternative normalizations (for comparison) - per-batch
  log_info("Saving alternative normalizations (median, combat) for comparison")

  # Median normalization - both batches
  batch1_median_path <- get_output_path(step_num, "npx_matrix_cross_batch_median", batch1_id, "normalized", config = config)
  batch2_median_path <- get_output_path(step_num, "npx_matrix_cross_batch_median", batch2_id, "normalized", config = config)
  ensure_output_dir(batch1_median_path)
  ensure_output_dir(batch2_median_path)

  saveRDS(batch1_normalized_median, batch1_median_path)
  saveRDS(batch2_normalized_median, batch2_median_path)
  log_info("  Saved median normalized: batch_01 and batch_02")

  # ComBat normalization - both batches (if available)
  if(!is.null(batch1_normalized_combat) && !is.null(batch2_normalized_combat)) {
    batch1_combat_path <- get_output_path(step_num, "npx_matrix_cross_batch_combat", batch1_id, "normalized", config = config)
    batch2_combat_path <- get_output_path(step_num, "npx_matrix_cross_batch_combat", batch2_id, "normalized", config = config)
    ensure_output_dir(batch1_combat_path)
    ensure_output_dir(batch2_combat_path)

    saveRDS(batch1_normalized_combat, batch1_combat_path)
    saveRDS(batch2_normalized_combat, batch2_combat_path)
    log_info("  Saved ComBat normalized: batch_01 and batch_02")
  }

  # Save evaluation results
  log_info("Saving evaluation results")
  evaluations_path <- get_output_path(step_num, "normalization_evaluations", output_batch_id, "normalized", "tsv", config = config)
  ensure_output_dir(evaluations_path)
  fwrite(all_evaluations, evaluations_path, sep = "\t")
  log_info("  Saved evaluation results: {evaluations_path}")

  # Save visualizations
  log_info("Saving visualizations")

  # ===========================================================================
  # PCA Comparison Plots (side-by-side before vs after)
  # ===========================================================================
  # Save combined before/after comparison as the primary PCA output
  if(!is.null(pca_comparison)) {
    # PC1 vs PC2 comparison
    pca_pc12_path <- get_output_path(step_num, "pca_pc1_pc2_before_after_comparison", output_batch_id, "normalized", "pdf", config = config)
    ensure_output_dir(pca_pc12_path)
    ggsave(pca_pc12_path, pca_comparison$pc12_comparison, width = 14, height = 7)
    log_info("  Saved PCA PC1 vs PC2 comparison: {pca_pc12_path}")

    # PC3 vs PC4 comparison (if available)
    if(!is.null(pca_comparison$pc34_comparison)) {
      pca_pc34_path <- get_output_path(step_num, "pca_pc3_pc4_before_after_comparison", output_batch_id, "normalized", "pdf", config = config)
      ensure_output_dir(pca_pc34_path)
      ggsave(pca_pc34_path, pca_comparison$pc34_comparison, width = 14, height = 7)
      log_info("  Saved PCA PC3 vs PC4 comparison: {pca_pc34_path}")
    }

    # Full combined (PC1-2 and PC3-4 in one figure)
    pca_full_path <- get_output_path(step_num, "pca_batch_comparison_full", output_batch_id, "normalized", "pdf", config = config)
    ensure_output_dir(pca_full_path)
    ggsave(pca_full_path, pca_comparison$full_combined, width = 14, height = 12)
    log_info("  Saved PCA full comparison: {pca_full_path}")
  }

  # Also save individual before/after plots for backward compatibility
  if(!is.null(pca_before)) {
    pca_before_path <- get_output_path(step_num, "pca_before_normalization_batch_comparison", output_batch_id, "normalized", "pdf", config = config)
    ensure_output_dir(pca_before_path)
    ggsave(pca_before_path, pca_before$combined, width = 14, height = 7)
  }

  if(!is.null(pca_after_bridge)) {
    pca_after_path <- get_output_path(step_num, "pca_after_normalization_batch_comparison", output_batch_id, "normalized", "pdf", config = config)
    ensure_output_dir(pca_after_path)
    ggsave(pca_after_path, pca_after_bridge$combined, width = 14, height = 7)
  }

  # ===========================================================================
  # Bridge Sample Scatter Plots
  # ===========================================================================
  if(!is.null(bridge_scatter)) {
    scatter_path <- get_output_path(step_num, "bridge_samples_scatter_before_after", output_batch_id, "normalized", "pdf", config = config)
    ensure_output_dir(scatter_path)
    ggsave(scatter_path, bridge_scatter$combined, width = 14, height = 7)
    log_info("  Saved bridge sample scatter: {scatter_path}")
    log_info("    Correlation before: {round(bridge_scatter$cor_before, 3)}, after: {round(bridge_scatter$cor_after, 3)}")
  }

  # ===========================================================================
  # Bridge Sample Paired Boxplot (harmonisation effectiveness)
  # ===========================================================================
  if(!is.null(bridge_paired_boxplot)) {
    boxplot_path <- get_output_path(step_num, "bridge_samples_paired_boxplot", output_batch_id, "normalized", "pdf", config = config)
    ensure_output_dir(boxplot_path)
    ggsave(boxplot_path, bridge_paired_boxplot$combined, width = 14, height = 14)
    log_info("  Saved bridge paired boxplot: {boxplot_path}")
    log_info("    Batch difference reduction: {round(bridge_paired_boxplot$stats$reduction_pct, 1)}%")
  }

  # ===========================================================================
  # Distribution Plots (histogram + density + violin with stats)
  # ===========================================================================
  # NOTE: New format includes histogram+density panel on top, violin with t-test on bottom
  if(!is.null(dist_plot_bridge)) {
    dist_bridge_path <- get_output_path(step_num, "normalization_effect_bridge", output_batch_id, "normalized", "pdf", config = config)
    ensure_output_dir(dist_bridge_path)
    ggsave(dist_bridge_path, dist_plot_bridge, width = 14, height = 10)
    log_info("  Saved bridge distribution plot: {dist_bridge_path}")
  }

  if(!is.null(dist_plot_median)) {
    dist_median_path <- get_output_path(step_num, "normalization_effect_median", output_batch_id, "normalized", "pdf", config = config)
    ensure_output_dir(dist_median_path)
    ggsave(dist_median_path, dist_plot_median, width = 14, height = 10)
    log_info("  Saved median distribution plot: {dist_median_path}")
  }

  if(!is.null(dist_plot_combat)) {
    dist_combat_path <- get_output_path(step_num, "normalization_effect_combat", output_batch_id, "normalized", "pdf", config = config)
    ensure_output_dir(dist_combat_path)
    ggsave(dist_combat_path, dist_plot_combat, width = 14, height = 10)
    log_info("  Saved ComBat distribution plot: {dist_combat_path}")
  }


  # Print summary with comprehensive formatting
  common_bridge_count <- length(intersect(batch1_bridge_finngenids, batch2_bridge_finngenids))

  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════════════════════╗\n")
  cat("║              CROSS-BATCH BRIDGE NORMALIZATION SUMMARY                        ║\n")
  cat("╚══════════════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  cat("┌─────────────────────────────────────────────────────────────────────────────┐\n")
  cat("│ DATA SUMMARY                                                                │\n")
  cat("├─────────────────────────────────────────────────────────────────────────────┤\n")
  cat(sprintf("│ Batches processed:        %-50s│\n", paste(all_batch_ids, collapse = ", ")))
  cat(sprintf("│ Common bridge FINNGENIDs: %-50d│\n", common_bridge_count))
  cat(sprintf("│ Common proteins:          %-50d│\n", length(common_proteins)))
  cat("│                                                                             │\n")
  cat("│ Bridge samples used:                                                        │\n")
  cat(sprintf("│   - Batch 1 (%s): %d samples                                         │\n", batch1_id, length(cross_batch_result_bridge$batch1_bridge_samples_used)))
  cat(sprintf("│   - Batch 2 (%s): %d samples                                         │\n", batch2_id, length(cross_batch_result_bridge$batch2_bridge_samples_used)))
  cat("└─────────────────────────────────────────────────────────────────────────────┘\n")
  cat("\n")
  cat("┌─────────────────────────────────────────────────────────────────────────────┐\n")
  cat("│ METHOD COMPARISON (SD REDUCTION %)                                          │\n")
  cat("├─────────────────────────────────────────────────────────────────────────────┤\n")
  cat(sprintf("│ Bridge:  Batch1=%6.2f%%, Batch2=%6.2f%%, Average=%6.2f%%              │\n",
      eval_bridge$eval_stats[batch == "batch_01"]$sd_reduction_pct,
      eval_bridge$eval_stats[batch == "batch_02"]$sd_reduction_pct,
      bridge_sd_reduction))
  cat(sprintf("│ Median:  Batch1=%6.2f%%, Batch2=%6.2f%%, Average=%6.2f%%              │\n",
      eval_median$eval_stats[batch == "batch_01"]$sd_reduction_pct,
      eval_median$eval_stats[batch == "batch_02"]$sd_reduction_pct,
      median_sd_reduction))
  if(!is.null(eval_combat)) {
    cat(sprintf("│ ComBat:  Batch1=%6.2f%%, Batch2=%6.2f%%, Average=%6.2f%%              │\n",
        eval_combat$eval_stats[batch == "batch_01"]$sd_reduction_pct,
        eval_combat$eval_stats[batch == "batch_02"]$sd_reduction_pct,
        combat_sd_reduction))
  }
  cat("└─────────────────────────────────────────────────────────────────────────────┘\n")
  cat("\n")
  cat("┌─────────────────────────────────────────────────────────────────────────────┐\n")
  cat("│ ★ BEST METHOD SELECTION                                                     │\n")
  cat("├─────────────────────────────────────────────────────────────────────────────┤\n")
  cat(sprintf("│ Best method:     %-58s│\n", toupper(best_method)))
  cat(sprintf("│ SD improvement:  %.2f%%                                                     │\n", best_sd_reduction))
  cat("│                                                                             │\n")
  cat("│ OUTPUT FILES:                                                               │\n")
  cat("│   Per-batch: 07_npx_matrix_cross_batch_bridge_{batch}.rds                   │\n")
  cat("│   Combined:  07_cross_batch_normalization_result_{batch}.rds                │\n")
  cat("│                                                                             │\n")
  cat("│ Recommendation: Use cross-batch bridge normalized matrices for downstream   │\n")
  cat("└─────────────────────────────────────────────────────────────────────────────┘\n")
  cat("\n")
  cat(sprintf("Results saved to: output/normalized/%s/\n", batch1_id))
  cat("══════════════════════════════════════════════════════════════════════════════\n")

  log_info("Cross-batch bridge normalization completed successfully")

  return(list(
    batch1_normalized_bridge = batch1_normalized_bridge,
    batch2_normalized_bridge = batch2_normalized_bridge,
    batch1_normalized_median = batch1_normalized_median,
    batch2_normalized_median = batch2_normalized_median,
    batch1_normalized_combat = batch1_normalized_combat,
    batch2_normalized_combat = batch2_normalized_combat,
    evaluation = all_evaluations,
    best_method = best_method
  ))
}

# Run if executed directly
if (!interactive()) {
  result <- main()
}
