suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
})

# Walk upward from a starting directory until a folder containing the expected
# published project structure (`code/`, `data/`, `results/`) is found.
# This lets the same scripts work when called from the project root, from an
# analysis subfolder, or from an R Markdown render context.
pc_locate_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = FALSE)

  repeat {
    has_code <- dir.exists(file.path(current, "code"))
    has_data <- dir.exists(file.path(current, "data"))
    has_results <- dir.exists(file.path(current, "results"))

    if (has_code && has_data && has_results) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the `published_code` project root from: ", start)
    }
    current <- parent
  }
}

# Try multiple context-aware starting points to find the project root.
# Priority is:
# 1. current working directory,
# 2. active R Markdown input directory during render,
# 3. directory inferred from command-line invocation.
# This reduces brittle assumptions about where a script is launched from.
pc_source_context_root <- function() {
  candidates <- c(
    getwd(),
    tryCatch(dirname(knitr::current_input(dir = TRUE)), error = function(...) NA_character_),
    tryCatch(dirname(normalizePath(commandArgs(trailingOnly = FALSE)[1], mustWork = FALSE)), error = function(...) NA_character_)
  )

  candidates <- unique(stats::na.omit(candidates))
  for (candidate in candidates) {
    try_root <- tryCatch(pc_locate_root(candidate), error = function(...) NULL)
    if (!is.null(try_root)) {
      return(try_root)
    }
  }

  pc_locate_root(getwd())
}

# Build the standard output layout for one analysis family.
# Every analysis gets a stable folder contract under `results/<analysis>/`:
# - `tables_dir` for processed tables and intermediate datasets,
# - `models_dir` for model summaries and final numeric outputs,
# - `figures_dir` for analysis-specific plots,
# - `diagnostics_dir` for assumption checks and QC artifacts.
# Scripts create only the subdirectories they actually write to, rather than
# eagerly materializing every possible output folder.
pc_analysis_paths <- function(analysis_name) {
  project_root <- pc_source_context_root()
  results_root <- file.path(project_root, "results")
  analysis_root <- file.path(results_root, analysis_name)

  paths <- list(
    project_root = project_root,
    code_root = file.path(project_root, "code"),
    data_root = file.path(project_root, "data"),
    results_root = results_root,
    analysis_root = analysis_root,
    tables_dir = file.path(analysis_root, "tables"),
    models_dir = file.path(analysis_root, "models"),
    figures_dir = file.path(analysis_root, "figures"),
    diagnostics_dir = file.path(analysis_root, "diagnostics")
  )

  paths
}

# Create one or more directories only when a script actually needs them.
pc_ensure_dirs <- function(...) {
  dirs <- unlist(list(...), use.names = FALSE)
  dirs <- dirs[!is.na(dirs) & nzchar(dirs)]
  dirs <- unique(dirs)

  purrr::walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  invisible(dirs)
}

# Load head-motion covariates from the published input file and reshape them to
# one row per participant per scan. This is the canonical helper used by the
# eccentricity, dispersion, and figure workflows when scan-level motion is
# needed as a covariate or metadata join.
pc_mean_fd <- function(data_root = file.path(pc_source_context_root(), "data"), scans = c("rest2", "rest3")) {
  readr::read_csv(file.path(data_root, "mean_FD.csv"), show_col_types = FALSE) %>%
    rename(NECOS_ID = ID, Group = stress) %>%
    pivot_longer(
      cols = starts_with("rest"),
      names_to = "scan",
      values_to = "headmotion"
    ) %>%
    mutate(scan = as.character(scan)) %>%
    filter(scan %in% scans)
}

# Load the 400-parcel label table and normalize column naming differences across
# historical CSV variants. Some source files use `X7networks_*` while others use
# `7networks_*`; this helper maps both onto a consistent schema for downstream
# code so the analysis scripts can stay agnostic to header quirks.
pc_parcel_labels <- function(data_root = file.path(pc_source_context_root(), "data")) {
  parcel_labels <- readr::read_csv(
    file.path(data_root, "yeo_networks_17_reordered_400_complete_labels.csv"),
    show_col_types = FALSE
  )

  order_col <- c("X7networks_order", "7networks_order")[c("X7networks_order", "7networks_order") %in% names(parcel_labels)][1]
  network_col <- c("X7networks", "7networks")[c("X7networks", "7networks") %in% names(parcel_labels)][1]

  if (is.na(order_col) || is.na(network_col)) {
    stop("Could not find 7-network parcel label columns in `yeo_networks_17_reordered_400_complete_labels.csv`.")
  }

  parcel_labels %>%
    mutate(
      X7networks_order = as.integer(.data[[order_col]]),
      X7networks = as.character(.data[[network_col]]),
      parcel = dplyr::coalesce(.data$X7networks_order, row_number())
    )
}

# Index all published gradient TSV files for the requested sessions.
# The helper parses participant IDs and session labels directly from the BIDS
# paths under `data/necos_gradients/derivatives/gradients/`, then adds a
# human-readable scan label (`rest1`-`rest6`) used throughout the R workflows.
pc_list_gradient_files <- function(sessions = c("rs2", "rs3"), data_root = file.path(pc_source_context_root(), "data")) {
  gradient_root <- file.path(data_root, "necos_gradients", "derivatives", "gradients")
  files <- list.files(
    gradient_root,
    pattern = "_desc-gradientpositions\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )

  tibble(file_path = files) %>%
    mutate(
      file_path = normalizePath(file_path, winslash = "/", mustWork = FALSE),
      NECOS_ID = str_match(file_path, "sub-(NECOS\\d+)")[, 2],
      session = str_match(file_path, "ses-(rs\\d+)")[, 2],
      scan = dplyr::recode(
        session,
        rs1 = "rest1",
        rs2 = "rest2",
        rs3 = "rest3",
        rs4 = "rest4",
        rs5 = "rest5",
        rs6 = "rest6",
        .default = session
      )
    ) %>%
    filter(!is.na(NECOS_ID), !is.na(session), session %in% sessions) %>%
    arrange(NECOS_ID, session)
}

# Read one BIDS-style gradient TSV and standardize it to the common internal
# schema: `parcel`, `V1`, `V2`, `V3`, plus participant/session metadata. This
# shields downstream scripts from the on-disk column naming used by the Python
# gradient pipeline.
pc_read_gradient_file <- function(file_path, NECOS_ID = NULL, session = NULL, scan = NULL) {
  gradient_df <- readr::read_tsv(file_path, show_col_types = FALSE) %>%
    transmute(
      parcel = as.integer(parcel_index),
      V1 = as.numeric(gradient_1),
      V2 = as.numeric(gradient_2),
      V3 = as.numeric(gradient_3)
    )

  gradient_df$NECOS_ID <- NECOS_ID %||% stringr::str_match(file_path, "sub-(NECOS\\d+)")[, 2]
  gradient_df$session <- session %||% stringr::str_match(file_path, "ses-(rs\\d+)")[, 2]
  gradient_df$scan <- scan %||% dplyr::recode(
    gradient_df$session[1],
    rs1 = "rest1",
    rs2 = "rest2",
    rs3 = "rest3",
    rs4 = "rest4",
    rs5 = "rest5",
    rs6 = "rest6",
    .default = gradient_df$session[1]
  )

  gradient_df
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)[1]) y else x
}

pc_read_gradients <- function(sessions = c("rs2", "rs3"), data_root = file.path(pc_source_context_root(), "data")) {
  gradient_index <- pc_list_gradient_files(sessions = sessions, data_root = data_root)

  purrr::pmap_dfr(
    gradient_index,
    function(file_path, NECOS_ID, session, scan) {
      pc_read_gradient_file(file_path, NECOS_ID = NECOS_ID, session = session, scan = scan)
    }
  )
}

# Read one session as a list of participant-by-parcel gradient matrices.
# This is mainly used by the figure workflow, which needs matrix-like access to
# the three aligned gradients per participant rather than a single long table.
pc_read_session_gradient_matrices <- function(session = "rs2", data_root = file.path(pc_source_context_root(), "data")) {
  gradient_index <- pc_list_gradient_files(sessions = session, data_root = data_root)

  matrices <- purrr::map(gradient_index$file_path, function(file_path) {
    gradient_df <- pc_read_gradient_file(file_path)
    as.matrix(gradient_df[, c("V1", "V2", "V3")])
  })

  list(
    files = gradient_index$file_path,
    participants = gradient_index$NECOS_ID,
    data = matrices,
    session = session,
    scan = unique(gradient_index$scan)
  )
}

# Convenience wrapper to read multiple sessions at once and return a named list
# keyed by session label (`rs1` ... `rs6`).
pc_read_all_session_gradient_matrices <- function(sessions = paste0("rs", 1:6), data_root = file.path(pc_source_context_root(), "data")) {
  stats::setNames(
    purrr::map(sessions, pc_read_session_gradient_matrices, data_root = data_root),
    sessions
  )
}

# Small wrapper around `write_csv()` that guarantees parent directories exist.
# This keeps output writes terse and consistent across scripts.
pc_write_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(df, path)
  invisible(path)
}
