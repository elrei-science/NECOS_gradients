#!/usr/bin/env Rscript

# NECOS rest2/rest3 post-hoc connectivity: all three mixed-model families.
#
# CONTENTS
#  1. Editable paths and analysis settings
#  2. Resolve paths and optional command-line arguments
#  3. Repeated table-reading and mixed-model helpers
#  4. Read the participant and motion covariates
#  5. Fit 399 parcelwise models per seed and save within-seed FDR results
#  6. Fit seven Yeo-network models per seed and save within-seed FDR results
#  7. Fit the ten edge-averaged network models and apply their shared FDR
#  8. Finish
#
# Read/run the sections from top to bottom. RStudio recognizes the ---- markers.
# Paths/settings can be edited in section 1, or overridden on the command line.
# Analysis loops run directly in the workspace; only repeated operations are
# functions. There is no main() wrapper.
#
# Run after NECOS_poshoc_seedbased_meannetworkconnectivity_rs2_r3.py and before
# code/r_scripts/figures/supplementary_figures.py. Example, from the repository root:
#
# Rscript code/r_scripts/functional_connectivity/NECOS_poshoc_seedbased_meannetworkconnectivity_rs2_r3_lmm.R
#
# Paths default to data/ and results/ within the repository; no workstation
# paths are embedded. Use --project-root /path/to/project if running elsewhere,
# or override --seedbased-root, --network-root, --motion-file and --label-file.
# The script can also be sourced in R (including through an R environment wrapper).
#
# Every outcome uses lmerTest's REML fit and Satterthwaite coefficient p-values:
# connectivity ~ stress * scan + headmotion + (1 | ID)
# stress: 0=control, 1=stress; reference scan: rest2.
# Thus stress:scanrest3 is the stress-minus-control difference in scan change.
#
# Primary BH-FDR families are deliberately kept separate:
#   399 parcel targets PER SEED; seven target networks PER SEED;
#   ten edge-averaged outcomes (seven within + three planned between networks).
# Additional Bonferroni/pooled network corrections are also written to the output tables;
# they are not the within-seed FDR values reported in the manuscript.
# Dependency: lmerTest (and its dependencies).



# 1. Editable paths and analysis settings ----
# NULL selects the repository-relative default. Example:
# PROJECT_ROOT <- "/path/to/project"
PROJECT_ROOT <- NULL
SEEDBASED_ROOT <- NULL
NETWORK_ROOT <- NULL
MOTION_FILE <- NULL
LABEL_FILE <- NULL
RUN_STAGE <- "all"  # Alternatives: "parcel", "seed-network", "mean-network".


# 2. Resolve paths and optional command-line arguments ----
find_project_root <- function() {
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  sourced_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
  starts <- c(dirname(sub("^--file=", "", file_argument)), if (length(sourced_files)) dirname(sourced_files), getwd())
  for (start in unique(starts)) {
    current <- normalizePath(start, winslash = "/", mustWork = FALSE)
    repeat {
      if (dir.exists(file.path(current, "code")) && dir.exists(file.path(current, "data"))) {
        return(current)
      }
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  stop("Cannot locate the repository. Supply --project-root /path/to/project.")
}

parse_arguments <- function() {
  arguments <- if (interactive()) character() else commandArgs(trailingOnly = TRUE)
  options <- list(project_root = PROJECT_ROOT, seedbased_root = SEEDBASED_ROOT, network_root = NETWORK_ROOT,
                  motion_file = MOTION_FILE, label_file = LABEL_FILE, stage = RUN_STAGE)
  while (length(arguments) > 0) {
    option <- arguments[[1]]
    arguments <- arguments[-1]
    if (option %in% c("-h", "--help")) {
      cat("NECOS post-hoc connectivity mixed models\n\n",
          "Run Python preparation first, then this script, then Python plotting.\n",
          "Options (each takes a value):\n",
          "  --project-root PATH   Repository containing code/ and data/.\n",
          "  --seedbased-root PATH Default: <project>/results/posthoc_seedbased\n",
          "  --network-root PATH   Default: <project>/results/posthoc_withinbetween_fc\n",
          "  --motion-file PATH    Default: <project>/data/mean_FD.csv\n",
          "  --label-file PATH     Default: <project>/data/yeo_networks_17_reordered_400_complete_labels.csv\n",
          "  --stage NAME          all, parcel, seed-network, or mean-network\n",
          "                        seed-network requires parcel analytic datasets.\n", sep = "")
      return(NULL)
    }
    name <- gsub("-", "_", sub("^--", "", option))
    if (!startsWith(option, "--") || !name %in% names(options) || length(arguments) == 0) {
      stop("Unknown option or missing value: ", option)
    }
    options[[name]] <- arguments[[1]]
    arguments <- arguments[-1]
  }
  if (!options$stage %in% c("all", "parcel", "seed-network", "mean-network")) {
    stop("--stage must be all, parcel, seed-network, or mean-network")
  }
  if (is.null(options$project_root)) options$project_root <- find_project_root()
  root <- options$project_root
  if (is.null(options$seedbased_root)) options$seedbased_root <- file.path(root, "results", "posthoc_seedbased")
  if (is.null(options$network_root)) options$network_root <- file.path(root, "results", "posthoc_withinbetween_fc")
  if (is.null(options$motion_file)) options$motion_file <- file.path(root, "data", "mean_FD.csv")
  if (is.null(options$label_file)) options$label_file <- file.path(root, "data", "yeo_networks_17_reordered_400_complete_labels.csv")
  options
}


options <- parse_arguments()


# 3. Repeated table-reading and mixed-model helpers ----
# The model-family loops are visible in sections 5-7.
write_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(data, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, ". Run the preceding pipeline stage.")
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

read_covariates <- function(path) {
  motion <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("ID", "stress", "rest2", "rest3") %in% names(motion))) {
    stop("Motion table requires ID, stress, rest2 and rest3.")
  }
  if (anyNA(motion$ID) || anyDuplicated(motion$ID)) stop("Motion IDs must be unique and present.")
  if (!setequal(motion$stress, c(0, 1)) || anyNA(motion$stress)) stop("stress must be 0=control or 1=stress, with both groups present.")
  if (!all(is.finite(as.matrix(motion[c("stress", "rest2", "rest3")])))) {
    stop("Motion covariates must be finite for both scans.")
  }
  rbind(
    data.frame(ID = motion$ID, scan = "rest2", headmotion = motion$rest2, stress = motion$stress),
    data.frame(ID = motion$ID, scan = "rest3", headmotion = motion$rest3, stress = motion$stress)
  )
}

read_labels <- function(path) {
  labels <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("7networks_order", "7networks") %in% names(labels)) ||
      !identical(sort(as.integer(labels$`7networks_order`)), 1:400)) {
    stop("Labels must contain each Schaefer parcel 1..400 exactly once.")
  }
  if (!setequal(labels$`7networks`, c("Visual", "Somatomotor", "Dorsal Attn",
                                    "Salience / Ventral Attn", "Limbic", "Control", "Default"))) {
    stop("Labels must contain the expected seven Yeo networks.")
  }
  labels[order(labels$`7networks_order`), , drop = FALSE]
}

join_covariates <- function(data, covariates) {
  if (anyDuplicated(data[c("ID", "scan")])) stop("Duplicate participant/scan rows in seed data.")
  joined <- merge(data, covariates, by = c("ID", "scan"), all = FALSE, sort = TRUE)
  if (nrow(joined) != nrow(data)) stop("Some scan rows lack motion covariates.")
  joined$scan <- factor(joined$scan, levels = c("rest2", "rest3"))
  joined
}

fit_outcome <- function(outcome, data, metadata) {
  # One shared implementation ensures identical model coding for all families.
  model_data <- data.frame(connectivity = as.numeric(outcome),
                           stress = as.numeric(data$stress),
                           scan = factor(data$scan, levels = c("rest2", "rest3")),
                           headmotion = as.numeric(data$headmotion), ID = factor(data$ID))
  if (!all(complete.cases(model_data)) || !all(is.finite(model_data$connectivity))) {
    stop("Nonfinite values in a modeled outcome or its covariates.")
  }
  model <- lmerTest::lmer(connectivity ~ stress * scan + headmotion + (1 | ID),
                          data = model_data, REML = TRUE)
  coefficients <- coef(summary(model))
  if (!"stress:scanrest3" %in% rownames(coefficients)) stop("Missing stress-by-scan coefficient.")
  effects <- data.frame(term = rownames(coefficients), estimate = coefficients[, "Estimate"],
                        std_error = coefficients[, "Std. Error"], df = coefficients[, "df"],
                        t_value = coefficients[, "t value"], p_value = coefficients[, "Pr(>|t|)"],
                        n_observations = nrow(model_data), n_subjects = length(unique(model_data$ID)),
                        row.names = NULL)
  cbind(metadata[rep(1, nrow(effects)), , drop = FALSE], effects, row.names = NULL)
}

interaction_only <- function(effects) {
  effects[effects$term == "stress:scanrest3", , drop = FALSE]
}



# 4. Read the participant and motion covariates ----
if (!is.null(options)) {
  if (!requireNamespace("lmerTest", quietly = TRUE)) {
    stop("Install the R package lmerTest before running models.")
  }
  covariates <- read_covariates(options$motion_file)
}


# 5. Fit 399 parcelwise models per seed and save within-seed FDR results ----
if (!is.null(options) && options$stage %in% c("all", "parcel")) {
  out <- file.path(options$seedbased_root, "LMM")
  summaries <- list()
  for (seed in c(100L, 377L)) {
    scans <- lapply(c(2L, 3L), function(rest) {
      data <- read_tsv(file.path(options$seedbased_root, paste0("rs_", rest),
                                sprintf("seed_%03d", seed), "subject_seed_maps_fisher_z.tsv"))
      data$ID <- data$subject
      data$scan <- paste0("rest", rest)
      data$seed_parcel <- seed
      data
    })
    data <- join_covariates(do.call(rbind, scans), covariates)
    targets <- setdiff(1:400, seed)
    if (!all(sprintf("parcel_%03d", targets) %in% names(data))) stop("Missing seed-target columns.")
    write_tsv(data, file.path(out, sprintf("seed_%03d_analytic_dataset.tsv", seed)))
    effects <- list()
    failures <- data.frame(seed_parcel = integer(), target_parcel = integer(), error = character())
    for (target in targets) {
      result <- tryCatch(
        fit_outcome(data[[sprintf("parcel_%03d", target)]], data,
                    data.frame(seed_parcel = seed, target_parcel = target)),
        error = function(error) error
      )
      if (inherits(result, "error")) {
        failures <- rbind(failures, data.frame(seed_parcel = seed, target_parcel = target,
                                              error = conditionMessage(result)))
      } else {
        effects[[length(effects) + 1]] <- result
      }
    }
    write_tsv(failures, file.path(out, sprintf("seed_%03d_model_failures.tsv", seed)))
    # Do not silently shrink the 399-test family if an individual model fails.
    if (nrow(failures) > 0) stop("Parcel models failed; inspect the model_failures table for seed ", seed)
    effects <- do.call(rbind, effects)
    interactions <- interaction_only(effects)
    stopifnot(nrow(interactions) == 399L)
    interactions$p_adjusted_bonferroni <- p.adjust(interactions$p_value, method = "bonferroni")
    interactions$p_adjusted_fdr_bh <- p.adjust(interactions$p_value, method = "BH")
    write_tsv(effects, file.path(out, sprintf("seed_%03d_fixed_effects.tsv", seed)))
    write_tsv(interactions, file.path(out, sprintf("seed_%03d_interaction_results.tsv", seed)))
    summaries[[length(summaries) + 1]] <- data.frame(seed_parcel = seed,
      n_rows_joined = nrow(data), n_subjects_joined = length(unique(data$ID)),
      n_models_completed = nrow(interactions), n_models_failed = nrow(failures))
    cat("Completed 399 parcel models for seed", seed, "\n")
  }
  write_tsv(do.call(rbind, summaries), file.path(out, "run_summary.tsv"))
}


# 6. Fit seven Yeo-network models per seed and save within-seed FDR results ----
if (!is.null(options) && options$stage %in% c("all", "seed-network")) {
  labels <- read_labels(options$label_file)
  root <- file.path(options$seedbased_root, "LMM")
  out <- file.path(root, "Yeo7")
  all_interactions <- list()
  for (seed in c(100L, 377L)) {
    data <- read_tsv(file.path(root, sprintf("seed_%03d_analytic_dataset.tsv", seed)))
    means <- data[c("ID", "scan", "stress", "headmotion")]
    effects <- list()
    for (network in unique(labels$`7networks`)) {
      # Average already Fisher-transformed seed edges, excluding self-connectivity.
      targets <- labels$`7networks_order`[labels$`7networks` == network]
      targets <- targets[targets != seed]
      outcome <- paste0("mean_fisher_z_", gsub("[^A-Za-z0-9]+", "_", network))
      means[[outcome]] <- rowMeans(data[sprintf("parcel_%03d", targets)], na.rm = FALSE)
      effects[[length(effects) + 1]] <- fit_outcome(means[[outcome]], data,
        data.frame(seed_parcel = seed, network = network, n_target_parcels = length(targets)))
    }
    effects <- do.call(rbind, effects)
    interactions <- interaction_only(effects)
    stopifnot(nrow(interactions) == 7L)
    interactions$p_adjusted_bonferroni_within_seed <- p.adjust(interactions$p_value, method = "bonferroni")
    interactions$p_adjusted_fdr_bh_within_seed <- p.adjust(interactions$p_value, method = "BH")
    write_tsv(means, file.path(out, sprintf("seed_%03d_yeo7_network_analytic_dataset.tsv", seed)))
    write_tsv(effects, file.path(out, sprintf("seed_%03d_yeo7_network_fixed_effects.tsv", seed)))
    write_tsv(interactions, file.path(out, sprintf("seed_%03d_yeo7_network_interactions.tsv", seed)))
    all_interactions[[length(all_interactions) + 1]] <- interactions
  }
  all_interactions <- do.call(rbind, all_interactions)
  all_interactions$p_adjusted_bonferroni_across_all_seed_networks <- p.adjust(all_interactions$p_value, method = "bonferroni")
  all_interactions$p_adjusted_fdr_bh_across_all_seed_networks <- p.adjust(all_interactions$p_value, method = "BH")
  write_tsv(all_interactions, file.path(out, "all_seed_yeo7_network_interactions.tsv"))
  cat("Completed seven Yeo-network models per seed.\n")
}


# 7. Fit the ten edge-averaged network models and apply their shared FDR ----
if (!is.null(options) && options$stage %in% c("all", "mean-network")) {
  out <- file.path(options$network_root, "LMM")
  connectivity <- read_tsv(file.path(options$network_root, "tables", "subject_network_connectivity_fisher_z.tsv"))
  if (anyDuplicated(connectivity[c("subject", "scan", "measure")])) stop("Duplicate network-outcome rows.")
  data <- merge(connectivity, covariates, by.x = c("subject", "scan"), by.y = c("ID", "scan"),
                all = FALSE, sort = TRUE)
  if (nrow(data) != nrow(connectivity)) stop("Network rows lack motion covariates.")
  data$scan <- factor(data$scan, levels = c("rest2", "rest3"))
  data$ID <- factor(data$subject)
  expected <- c("within_Visual", "within_Somatomotor", "within_Dorsal_Attn",
                "within_Salience_Ventral_Attn", "within_Limbic", "within_Control",
                "within_Default", "SN_DMN", "SN_FPN", "DMN_FPN")
  if (!setequal(data$measure, expected)) stop("Expected exactly seven within- and three between-network measures.")
  effects <- list()
  failures <- data.frame(measure = character(), error = character())
  for (measure in unique(data$measure)) {
    subset <- data[data$measure == measure, , drop = FALSE]
    metadata <- data.frame(measure = measure, measure_type = subset$measure_type[1],
                           network_or_pair = subset$network_or_pair[1], n_parcel_pairs = subset$n_parcel_pairs[1])
    result <- tryCatch(fit_outcome(subset$mean_fisher_z, subset, metadata), error = function(error) error)
    if (inherits(result, "error")) {
      failures <- rbind(failures, data.frame(measure = measure, error = conditionMessage(result)))
    } else {
      effects[[length(effects) + 1]] <- result
    }
  }
  write_tsv(failures, file.path(out, "model_failures.tsv"))
  if (nrow(failures) > 0) stop("Network models failed; inspect the model_failures table.")
  effects <- do.call(rbind, effects)
  interactions <- interaction_only(effects)
  interactions$p_adjusted_bonferroni_all_10 <- p.adjust(interactions$p_value, method = "bonferroni")
  interactions$p_adjusted_fdr_bh_all_10 <- p.adjust(interactions$p_value, method = "BH")
  write_tsv(data, file.path(out, "analytic_dataset.tsv"))
  write_tsv(effects, file.path(out, "fixed_effects.tsv"))
  write_tsv(interactions, file.path(out, "all_interactions_corrected.tsv"))
  write_tsv(interactions[interactions$measure_type == "within", ], file.path(out, "within_network_interactions.tsv"))
  write_tsv(interactions[interactions$measure_type == "between", ], file.path(out, "between_network_interactions.tsv"))
  cat("Completed ten edge-averaged network models.\n")
}


# 8. Finish ----
if (!is.null(options)) {
  cat("Finished. Next: run code/r_scripts/figures/supplementary_figures.py.\n")
}
