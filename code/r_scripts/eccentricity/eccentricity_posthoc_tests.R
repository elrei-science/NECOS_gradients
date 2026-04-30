#!/usr/bin/env Rscript

# Simple reconstruction of the corrected parcel summary tables.
# This script:
# 1. loads the prepared parcel datasets for parcels 100 and 377,
# 2. fits one mixed model per gradient,
# 3. extracts the stress:scan interaction p-value,
# 4. applies a Bonferroni correction within each parcel,
# 5. writes the corrected summary CSV files.

suppressPackageStartupMessages({
  library(lmerTest)
})

helper_candidates <- c(
  file.path(getwd(), "code", "r_scripts", "shared", "project_paths.R"),
  file.path(getwd(), "..", "shared", "project_paths.R"),
  file.path(getwd(), "shared", "project_paths.R")
)
helper_path <- helper_candidates[file.exists(helper_candidates)][1]
if (is.na(helper_path)) {
  stop("Could not locate `code/r_scripts/shared/project_paths.R`.")
}
source(helper_path)

paths <- pc_analysis_paths("eccentricity")

# ---- paths ----

processed_data_dir <- paths$tables_dir
output_dir <- paths$models_dir
parcel_ids <- c(100, 377)

pc_ensure_dirs(output_dir)


# ---- helper functions ----

load_parcel_data <- function(parcel_id) {
  file_path <- file.path(processed_data_dir, sprintf("eccentricity_parcel_%d.csv", parcel_id))

  if (!file.exists(file_path)) {
    stop(sprintf("Missing parcel data file: %s", file_path))
  }

  parcel_data <- read.csv(file_path, stringsAsFactors = FALSE)
  parcel_data$scan <- factor(parcel_data$scan, levels = c("rest2", "rest3"))
  parcel_data
}

fit_models <- function(parcel_data) {
  lapply(1:3, function(i) {
    model_formula <- as.formula(
      paste0("gradient", i, " ~ stress * scan + headmotion + (1 | ID)")
    )
    lmer(model_formula, data = parcel_data)
  })
}

extract_interaction_pvalues <- function(models) {
  sapply(models, function(model) {
    coefficient_table <- coef(summary(model))
    interaction_row <- grep("^stress:scan", rownames(coefficient_table), value = TRUE)

    if (length(interaction_row) != 1) {
      stop("Could not identify exactly one stress:scan interaction term.")
    }

    coefficient_table[interaction_row, "Pr(>|t|)"]
  })
}

build_summary_table <- function(parcel_id) {
  parcel_data <- load_parcel_data(parcel_id)
  models <- fit_models(parcel_data)
  interaction_pvalues <- extract_interaction_pvalues(models)

  data.frame(
    Parcel = parcel_id,
    Gradient = 1:3,
    Interaction_PValue = interaction_pvalues,
    P_Adjusted = p.adjust(interaction_pvalues, method = "bonferroni")
  )
}

write_summary_table <- function(parcel_id, summary_table) {
  output_path <- file.path(
    output_dir,
    sprintf("eccentricity_parcel_%d_model_summaries_corrected.csv", parcel_id)
  )

  write.csv(summary_table, output_path)
  output_path
}


# ---- run analysis ----

output_paths <- vapply(parcel_ids, function(parcel_id) {
  summary_table <- build_summary_table(parcel_id)
  write_summary_table(parcel_id, summary_table)
}, character(1))

cat("Wrote files:\n")
for (path in output_paths) {
  cat(" -", path, "\n")
}
