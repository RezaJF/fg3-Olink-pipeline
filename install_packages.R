#!/usr/bin/env Rscript

# Install required packages for FinnGen 3 Olink pipeline
cat("Installing required R packages...\n")

# CRAN packages
cran_packages <- c(
  "data.table",
  "tidyverse",
  "arrow",
  "yaml",
  "logger",
  "ggplot2",
  "ggrepel",
  "ggpubr",  # For paired plots and statistical comparisons (Step 07b KPI)
  "optparse",
  "moments",
  "gridExtra",
  "paletteer",
  "ggthemes",
  "dichromat",
  "glmnet",
  "pheatmap",
  "pROC",
  "PRROC",
  "e1071",
  "rpart",
  # Statistical and clustering packages (Step 07b KPI)
  "cluster",  # For silhouette scores
  "MASS",  # For robust covariance estimation
  "irr",  # For intraclass correlation coefficient (ICC)
  "DescTools",  # For concordance correlation coefficient (CCC)
  "FNN",  # For fast k-nearest neighbour search (vectorised rank computation)
  "lisi",  # For local inverse Simpson index (batch mixing metric)
  "nationalparkcolors",  # For Acadia colour palette (bridge protein correlation heatmaps)
  # Enhanced modeling dependencies
  "randomForest",
  "xgboost",
  "caret",
  # Prefer keras3 when available; also install keras for compatibility
  "keras3",
  "keras",
  "tensorflow"
)

# Check and install CRAN packages
for(pkg in cran_packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cloud.r-project.org/", quiet = TRUE)
  } else {
    cat(pkg, "already installed\n")
  }
}

# BioConductor packages
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

# Install BioConductor packages
bioc_packages <- c("sva", "kBET")  # sva for ComBat, kBET for batch effect testing (Step 07b KPI)

for(pkg in bioc_packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "from BioConductor...\n")
    BiocManager::install(pkg, quiet = TRUE, update = FALSE)
  } else {
    cat(pkg, "already installed\n")
  }
}

# OlinkAnalyze might need special handling
if(!requireNamespace("OlinkAnalyze", quietly = TRUE)) {
  cat("Installing OlinkAnalyze...\n")
  tryCatch({
    install.packages("OlinkAnalyze", repos = "https://cloud.r-project.org/")
  }, error = function(e) {
    cat("Note: OlinkAnalyze installation failed. May need manual installation.\n")
  })
}

cat("\nPackage installation complete!\n")
cat("Testing package loading...\n")

# Test loading key packages
suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(yaml)
  library(logger)
  library(ggplot2)
  library(ggpubr)
  library(cluster)
  library(MASS)
  library(FNN)
})

cat("Core packages loaded successfully!\n")

# Test BioConductor packages
if (requireNamespace("kBET", quietly = TRUE)) {
  cat("kBET package loaded successfully\n")
} else {
  cat("Warning: kBET package not available\n")
}

if (requireNamespace("lisi", quietly = TRUE)) {
  cat("lisi package loaded successfully\n")
} else {
  cat("Warning: lisi package not available\n")
}

# Keras/TensorFlow backend setup (optional)
if (requireNamespace("keras3", quietly = TRUE) || requireNamespace("keras", quietly = TRUE)) {
  cat("Checking TensorFlow backend...\n")
  ok <- FALSE
  if (requireNamespace("reticulate", quietly = TRUE)) {
    ok <- tryCatch(reticulate::py_module_available("tensorflow"), error = function(e) FALSE)
  }
  if (!ok) {
    cat("TensorFlow not detected. You may install it via:\n")
    cat("  library(keras3); keras3::install_tensorflow()\n")
    cat("or: library(keras); keras::install_keras()\n")
  } else {
    cat("TensorFlow backend detected.\n")
  }
}



