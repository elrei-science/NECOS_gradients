#!/usr/bin/env Rscript

# AUCi of the six peripheral stress measures and their association with the brain outcomes.
# This script:
# 1. computes AUCi (area under the curve with respect to increase) for HR, HRV, saliva
#    cortisol, ACTH, STAI and VAS from the pre-TSST anchor (-15 min) to the last sample,
# 2. transforms ACTH AUCi with asinh(x / MAD),
# 3. builds the vPFC seed-based connectivity change (rest3 - rest2),
# 4. tests Group x AUCi for the five brain outcomes (baseline-adjusted models of change),
# 5. applies Benjamini-Hochberg correction across the 30 tests.

helper_candidates <- c(
  file.path(getwd(), "code", "r_scripts", "shared", "project_io.R"),
  file.path(getwd(), "..", "shared", "project_io.R"),
  file.path(getwd(), "shared", "project_io.R")
)
helper_path <- helper_candidates[file.exists(helper_candidates)][1]
if (is.na(helper_path)) {
  stop("Could not locate `code/r_scripts/shared/project_io.R`.")
}
source(helper_path)

paths <- pc_analysis_paths("peripheral_auci")

# ---- paths ----

timeline_file <- file.path(paths$data_root, "peripheral_timeline_long.csv")
eccentricity_file <- file.path(paths$results_root, "peripheral", "tables", "peripheral_eccentricity_delta_base.csv")
dispersion_file <- file.path(paths$results_root, "peripheral", "tables", "peripheral_dispersion_delta_base.csv")
seed_dir <- "path/to/posthoc_seedbased"  # seed-based FC results (LMM/ and rs_2/, rs_3/ subfolders)
output_dir <- paths$models_dir
table_dir <- paths$tables_dir

# ---- settings ----

anchor_min <- -15        # last pre-TSST sample, AUCi baseline
min_coverage <- 0.8      # share of the expected samples required per participant
exclude_all <- "NECOS049"
exclude_endocrine <- c("NECOS049", "NECOS054")  # cortisol and ACTH only

# sample times (minutes relative to TSST onset)
endocrine_times <- tibble(
  Timepoint = as.character(3:14),
  tsst_min = c(-45, -30, -15, 5, 15, 25, 45, 60, 80, 95, 110, 130)
)
autonomic_times <- tibble(
  Timepoint = c("01_rest1a", "02_rest1b", "03_rest2a", "04_rest2b", "05_TSST_ant", "06_TSST_int",
                "07_TSST_mat", "08_rest3a", "09_rest3b", "10_rest4a", "11_rest4b", "12_rest5a",
                "13_rest5b", "14_rest6a", "15_rest6b"),
  tsst_min = c(-35, -30, -20, -15, 0, 10, 15, 40, 45, 60, 65, 100, 105, 120, 125)
)

# last value is carried forward to the final sample time for HR and HRV if it is at most
# this many minutes short of it
measures <- tribble(
  ~measure,            ~label,   ~times,          ~carry_forward_min,
  "Speichel_Cortisol", "SaCort", "endocrine",     0,
  "Plasma_ACTH",       "ACTH",   "endocrine",     0,
  "STAI",              "STAI",   "endocrine",     0,
  "VAS",               "VAS",    "endocrine",     0,
  "HR",                "HR",     "autonomic",     5,
  "HRV",               "HRV",    "autonomic",     5
)

# ---- AUCi ----

timeline <- read_csv(timeline_file, show_col_types = FALSE) %>%
  mutate(NECOS_ID = as.character(NECOS_ID), Group = as.character(Group), Timepoint = as.character(Timepoint)) %>%
  filter(!NECOS_ID %in% exclude_all)

compute_auci <- function(measure, label, times, carry_forward_min) {
  dat <- timeline %>% filter(measure == !!measure)
  if (measure %in% c("Speichel_Cortisol", "Plasma_ACTH")) dat <- dat %>% filter(!NECOS_ID %in% exclude_endocrine)

  if (times == "endocrine") {
    dat <- dat %>% mutate(Timepoint = gsub("^T", "", Timepoint)) %>% left_join(endocrine_times, by = "Timepoint")
  } else {
    dat <- dat %>% left_join(autonomic_times, by = "Timepoint")
  }

  window <- sort(unique(dat$tsst_min[!is.na(dat$tsst_min) & dat$tsst_min >= anchor_min]))
  window_end <- max(window)

  dat %>%
    filter(!is.na(tsst_min), !is.na(value), tsst_min >= anchor_min) %>%
    group_by(NECOS_ID, Group) %>%
    group_modify(function(d, key) {
      d <- arrange(d, tsst_min)
      baseline <- d$value[match(anchor_min, d$tsst_min)]
      n_observed <- nrow(d)

      gap <- window_end - max(d$tsst_min)
      if (gap > 0 && gap <= carry_forward_min) {
        d <- bind_rows(d, tibble(tsst_min = window_end, value = last(d$value)))
      }

      complete <- n_observed / length(window) >= min_coverage && !is.na(baseline) && window_end %in% d$tsst_min
      auc_ground <- sum(diff(d$tsst_min) * (head(d$value, -1) + tail(d$value, -1)) / 2)
      tibble(AUCi = if (complete) auc_ground - baseline * (window_end - anchor_min) else NA_real_)
    }) %>%
    ungroup() %>%
    mutate(label = label)
}

auci_wide <- pmap_dfr(measures, compute_auci) %>%
  pivot_wider(id_cols = c(NECOS_ID, Group), names_from = label, values_from = AUCi, names_glue = "{label}_AUCi") %>%
  mutate(ACTH_AUCi_asinh = asinh(ACTH_AUCi / mad(ACTH_AUCi, na.rm = TRUE)))

pc_write_csv(auci_wide, file.path(table_dir, "peripheral_auci_values.csv"))

# ---- vPFC seed-based connectivity change ----

# target parcels of the right vPFC seed (parcel 377) with an FDR-significant group x scan interaction
vpfc_targets <- read_tsv(file.path(seed_dir, "LMM", "seed_377_interaction_results.tsv"), show_col_types = FALSE) %>%
  filter(p_adjusted_fdr_bh < 0.05) %>%
  pull(target_parcel) %>%
  sort()

read_cluster_mean <- function(scan, rs) {
  read_tsv(file.path(seed_dir, paste0("rs_", rs), "seed_377", "subject_seed_maps_fisher_z.tsv"), show_col_types = FALSE) %>%
    transmute(
      NECOS_ID = subject,
      scan = scan,
      cluster_mean_fz = rowMeans(across(all_of(sprintf("parcel_%03d", vpfc_targets))), na.rm = TRUE)
    )
}

fc_vpfc <- bind_rows(read_cluster_mean("rest2", 2), read_cluster_mean("rest3", 3)) %>%
  inner_join(pc_mean_fd() %>% mutate(NECOS_ID = as.character(NECOS_ID), Group = as.character(Group)), by = c("NECOS_ID", "scan")) %>%
  filter(!NECOS_ID %in% exclude_all) %>%
  pivot_wider(id_cols = c(NECOS_ID, Group), names_from = scan, values_from = c(cluster_mean_fz, headmotion)) %>%
  mutate(
    delta_FC_vpfc_cluster = cluster_mean_fz_rest3 - cluster_mean_fz_rest2,
    log_mean_MFD_scaled = as.numeric(scale(log(rowMeans(cbind(headmotion_rest2, headmotion_rest3), na.rm = TRUE))))
  ) %>%
  rename(FC_vpfc_cluster_rest2 = cluster_mean_fz_rest2)

pc_write_csv(fc_vpfc, file.path(table_dir, "peripheral_auci_fc_vpfc_delta.csv"))

# ---- Group x AUCi tests ----

read_brain <- function(file) read_csv(file, show_col_types = FALSE) %>% mutate(NECOS_ID = as.character(NECOS_ID), Group = as.character(Group))
eccentricity <- read_brain(eccentricity_file)
dispersion <- read_brain(dispersion_file)

outcomes <- list(
  list(data = eccentricity, delta = "delta_parcel_100", baseline = "parcel_100_rest2"),
  list(data = eccentricity, delta = "delta_parcel_377", baseline = "parcel_377_rest2"),
  list(data = dispersion, delta = "delta_distance_default_salience", baseline = "distance_default_salience_rest2"),
  list(data = dispersion, delta = "delta_distance_control_salience", baseline = "distance_control_salience_rest2"),
  list(data = fc_vpfc, delta = "delta_FC_vpfc_cluster", baseline = "FC_vpfc_cluster_rest2")
)
predictors <- c("HR_AUCi", "HRV_AUCi", "SaCort_AUCi", "ACTH_AUCi_asinh", "STAI_AUCi", "VAS_AUCi")

test_interaction <- function(outcome, predictor) {
  dat <- outcome$data %>%
    select(NECOS_ID, Group, log_mean_MFD_scaled, all_of(c(outcome$delta, outcome$baseline))) %>%
    inner_join(auci_wide %>% select(NECOS_ID, all_of(predictor)), by = "NECOS_ID") %>%
    drop_na() %>%
    mutate(Group = factor(Group))

  fit <- lm(reformulate(c(paste("Group *", predictor), outcome$baseline, "log_mean_MFD_scaled"), response = outcome$delta), data = dat)
  term <- paste0("Group:", predictor)
  anova_table <- car::Anova(fit, type = 3)

  tibble(
    outcome = outcome$delta, predictor = predictor, n = nobs(fit),
    f_value = anova_table[term, "F value"], df = anova_table[term, "Df"], df_residual = df.residual(fit),
    p = anova_table[term, "Pr(>F)"]
  )
}

results <- map_dfr(outcomes, function(o) map_dfr(predictors, function(p) test_interaction(o, p))) %>%
  mutate(p_fdr = p.adjust(p, method = "BH")) %>%
  arrange(p)

pc_write_csv(results, file.path(output_dir, "peripheral_auci_brain_associations.csv"))
print(results, n = Inf)

# ---- correlation between the two vPFC outcomes ----

vpfc_outcomes <- fc_vpfc %>%
  select(NECOS_ID, delta_FC_vpfc_cluster) %>%
  inner_join(eccentricity %>% select(NECOS_ID, delta_parcel_377), by = "NECOS_ID") %>%
  drop_na()
correlation <- cor.test(vpfc_outcomes$delta_FC_vpfc_cluster, vpfc_outcomes$delta_parcel_377)

pc_write_csv(
  tibble(n = nrow(vpfc_outcomes), r = unname(correlation$estimate), p = correlation$p.value),
  file.path(table_dir, "peripheral_auci_vpfc_fc_vs_eccentricity.csv")
)
