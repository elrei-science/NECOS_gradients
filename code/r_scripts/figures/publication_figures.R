# =====================================
# NECOS Gradients Publication Plots
# =====================================
#
# Generates publication-ready figures for the NECOS gradients analyses.
#
# This script creates:
#
#   Figure 2A:
#     Six-resting-state eccentricity plot
#
#   Figure 2B:
#     Parcel-level eccentricity raincloud plots
#
#   Figure 2C:
#     Six-resting-state between-network dispersion plot
#
#   Figure 2D:
#     Network centroid ridge plots
#
# Group-based figures use the same Control and Stress colours throughout
# for visual consistency.
#

library(tidyverse)
library(ggplot2)
library(ggridges)
library(patchwork)
library(grid)


# ==============================
# Project paths
# ==============================

# The script needs access to shared path/helper functions.
# Because the script may be run from slightly different working directories,
# we define several possible locations where `project_io.R` might be found.
helper_candidates <- c(
  file.path(getwd(), "code", "r_scripts", "shared", "project_io.R"),
  file.path(getwd(), "..", "shared", "project_io.R"),
  file.path(getwd(), "shared", "project_io.R")
)

# Select the first existing helper file from the possible paths above.
helper_path <- helper_candidates[file.exists(helper_candidates)][1]

# Stop the script early if the helper file cannot be found.
# This prevents later errors caused by missing path functions.
if (is.na(helper_path)) {
  stop("Could not locate `code/r_scripts/shared/project_io.R`.")
}

# Load the shared project path functions, such as:
#   pc_source_context_root()
#   pc_analysis_paths()
source(helper_path)

# Identify the root folder of the project.
# This makes the script independent of the exact folder it is run from.
project_root <- pc_source_context_root()

# Store the main project folders in one list.
# This makes later file paths easier to read and less error-prone.
figure_paths <- list(
  project_root  = project_root,
  data_root     = file.path(project_root, "data"),
  results_root  = file.path(project_root, "results"),
  analysis_root = file.path(project_root, "results", "figures")
)

# Get standard analysis-specific paths for each analysis family.
# These objects usually contain paths such as tables_dir, models_dir, etc.
eccentricity_paths <- pc_analysis_paths("eccentricity")
dispersion_between_paths <- pc_analysis_paths(file.path("dispersion", "between_networks"))
peripheral_paths   <- pc_analysis_paths("peripheral")

# Main folder where publication figures will be saved.
publication_plot_dir <- figure_paths$analysis_root

# Figure-specific output folders.
# At the moment, all figure types are saved into the same publication folder.
# Keeping separate object names is useful because it makes the saving
# functions easier to read and easier to modify later.
parcel_plot_dir                 <- publication_plot_dir
six_rs_plot_dir                 <- publication_plot_dir
between_network_six_rs_plot_dir <- publication_plot_dir
physio_plot_dir                 <- publication_plot_dir
ridge_plot_dir                  <- publication_plot_dir

# Create the figure output directory if it does not already exist.
dir.create(publication_plot_dir, recursive = TRUE, showWarnings = FALSE)


# ==============================
# Shared colours and labels
# ==============================

# Consistent group colours used across the manuscript figures.
# These colours are used for Control and Stress groups in most plots.
group_colors <- c(
  Control = "#a1c9b6",
  Stress  = "#9c7cc1"
)

# Some timeline data uses lowercase group labels:
# "control" and "stress" instead of "Control" and "Stress".
timeline_group_colors <- c(
  control = group_colors[["Control"]],
  stress  = group_colors[["Stress"]]
)

# Format p-values for figure annotations.
#
# Examples:
#   NA       -> "NA"
#   0.0004   -> "<.001"
#   0.03456  -> "0.035"
fmt_p <- function(x) {
  case_when(
    is.na(x) ~ "NA",
    x < 0.001 ~ "<.001",
    TRUE ~ sprintf("%.3f", x)
  )
}

# Convert different possible group encodings into clean group labels.
#
# Input may be encoded as:
#   0 / 1
#   control / stress
#   Control / Stress
#
# Output is standardised to:
#   Control / Stress
pretty_group <- function(x) {
  x_chr <- as.character(x)
  
  case_when(
    x_chr %in% c("1", "stress", "Stress") ~ "Stress",
    x_chr %in% c("0", "control", "Control") ~ "Control",
    TRUE ~ x_chr
  )
}


# =========================
# Results figures: Eccentricity
# =========================
#
# This section generates two eccentricity figure sets:
#
#   Figure 2A:
#     Six-resting-state eccentricity plot for parcels 100 and 377
#
#   Figure 2B:
#     Parcel-level raincloud plots for parcels 100 and 377
#


# -------------------------------------------------------------------------
# Figure 2A:
# Six-resting-state eccentricity plot
# -------------------------------------------------------------------------

make_six_rs_eccentricity_plot <- function(data) {
  physical_distances <- c(
    rest1 = -35,
    rest2 = -20,
    rest3 = 35,
    rest4 = 55,
    rest5 = 90,
    rest6 = 105
  )
  
  plot_data <- data %>%
    mutate(
      scan = factor(scan, levels = names(physical_distances)),
      position = unname(physical_distances[as.character(scan)])
    )
  
  summary_data <- plot_data %>%
    group_by(parcel_label, stress, scan, position) %>%
    summarise(
      mean_value = mean(eccentricity, na.rm = TRUE),
      se_value = mean_se(eccentricity)$ymax - mean_se(eccentricity)$y,
      .groups = "drop"
    )
  
  ggplot(summary_data, aes(x = position, y = mean_value, color = stress, group = stress)) +
    geom_line(
      linewidth = 1.2,
      alpha = 0.85,
      linetype = "dotted"
    ) +
    geom_point(
      size = 3.2
    ) +
    geom_errorbar(
      aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
      width = 4,
      linewidth = 0.6
    ) +
    facet_wrap(~ parcel_label, ncol = 1) +
    scale_color_manual(values = group_colors, name = "Group") +
    scale_x_continuous(
      breaks = unname(physical_distances),
      labels = paste("Scan", 1:6)
    ) +
    labs(
      title = "Eccentricity Across Six Resting States",
      x = "Resting-state scan",
      y = "Eccentricity score"
    ) +
    theme_classic() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 12),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 11)
    )
}

save_six_rs_eccentricity_plot <- function(output_dir = six_rs_plot_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  six_rs_data <- build_six_rs_eccentricity_data()
  plot_object <- make_six_rs_eccentricity_plot(six_rs_data)
  
  save_plot_png_svg(
    plot_object,
    file_stem = "eccentricity_six_resting_states_parcel_100_377",
    output_dir = output_dir,
    width = 10,
    height = 8
  )
  
  plot_object
}


# -------------------------------------------------------------------------
# Figure 2B:
# Parcel-level eccentricity raincloud plots
# -------------------------------------------------------------------------

build_parcel_data <- function(parcel_number) {
  rs_2_bundle <- pc_read_session_gradient_matrices("rs2", data_root = figure_paths$data_root)
  rs_3_bundle <- pc_read_session_gradient_matrices("rs3", data_root = figure_paths$data_root)
  
  my_data_rs_2 <- rs_2_bundle$data
  my_data_rs_3 <- rs_3_bundle$data
  ids <- rs_2_bundle$participants
  n_participants <- length(ids)
  
  better_together <- pc_mean_fd(figure_paths$data_root, scans = c("rest2", "rest3")) %>%
    transmute(id = NECOS_ID, stress = Group, scan)
  
  parcel_metadata <- bind_rows(
    better_together %>%
      filter(scan == "rest2") %>%
      mutate(id = factor(id, levels = ids)) %>%
      arrange(id),
    better_together %>%
      filter(scan == "rest3") %>%
      mutate(id = factor(id, levels = ids)) %>%
      arrange(id)
  ) %>%
    mutate(id = as.character(id))
  
  parcel_matrix <- matrix(ncol = 3, nrow = n_participants * 2)
  
  for (index in seq_len(n_participants)) {
    parcel_matrix[index, 1:3] <- c(
      my_data_rs_2[[index]][parcel_number, 1],
      my_data_rs_2[[index]][parcel_number, 2],
      my_data_rs_2[[index]][parcel_number, 3]
    )
    
    parcel_matrix[n_participants + index, 1:3] <- c(
      my_data_rs_3[[index]][parcel_number, 1],
      my_data_rs_3[[index]][parcel_number, 2],
      my_data_rs_3[[index]][parcel_number, 3]
    )
  }
  
  parcel_data <- as.data.frame(cbind(rep(ids, 2), parcel_matrix))
  parcel_data <- cbind(parcel_data, parcel_metadata$stress, parcel_metadata$scan)
  colnames(parcel_data) <- c("id", "g1", "g2", "g3", "stress", "scan")
  
  parcel_data %>%
    mutate(
      g1 = as.numeric(g1),
      g2 = as.numeric(g2),
      g3 = as.numeric(g3),
      stress = as.character(stress),
      scan = as.character(scan),
      parcel = parcel_number,
      time = factor(if_else(scan == "rest2", "Pre", "Post"), levels = c("Pre", "Post")),
      group = factor(if_else(stress == "1", "Stress", "Control"), levels = c("Control", "Stress"))
    )
}

# Build half-violin polygon coordinates for the raincloud-style parcel plots.
build_half_violin_df <- function(data, value_col, width = 0.2, offset = 0.24) {
  split(data, interaction(data$time, data$group, drop = TRUE)) %>%
    purrr::imap_dfr(function(chunk, key) {
      values <- chunk[[value_col]]
      values <- values[is.finite(values)]
      
      if (length(values) < 2) {
        return(tibble())
      }
      
      density_values <- density(values, n = 512, na.rm = TRUE)
      scaled_width <- density_values$y / max(density_values$y) * width
      time_num <- unique(chunk$time_num)
      side <- unique(chunk$violin_side)
      x_center <- if (side == "left") time_num - offset else time_num + offset
      
      half_x <- if (side == "left") {
        x_center - scaled_width
      } else {
        x_center + scaled_width
      }
      
      tibble(
        x = c(rep(x_center, length(density_values$x)), rev(half_x)),
        y = c(density_values$x, rev(density_values$x)),
        time = unique(chunk$time),
        group = unique(chunk$group),
        polygon_id = key
      )
    })
}

# Create a parcel raincloud plot for one gradient at one parcel.
make_parcel_plot <- function(data, gradient_col, title_text) {
  plot_data <- data %>%
    transmute(
      id,
      time,
      group,
      value = .data[[gradient_col]]
    ) %>%
    filter(is.finite(value)) %>%
    mutate(
      time_num = if_else(time == "Pre", 1, 2),
      group_sign = if_else(group == "Control", -1, 1),
      violin_side = if_else(time == "Pre", "left", "right"),
      point_x = time_num + group_sign * 0.035,
      box_x = if_else(
        time == "Pre",
        time_num - 0.11 + if_else(group == "Control", -0.02, 0.02),
        time_num + 0.11 + if_else(group == "Control", -0.02, 0.02)
      ),
      summary_x = time_num + group_sign * 0.045
    ) %>%
    arrange(id, time)
  
  violin_df <- build_half_violin_df(plot_data, "value")
  
  summary_df <- plot_data %>%
    group_by(time, time_num, group, summary_x) %>%
    summarise(summary_value = median(value, na.rm = TRUE), .groups = "drop") %>%
    arrange(group, time_num)
  
  ggplot() +
    geom_polygon(
      data = violin_df,
      aes(x = x, y = y, fill = group, group = polygon_id),
      alpha = 0.65,
      color = NA
    ) +
    geom_boxplot(
      data = plot_data,
      aes(x = box_x, y = value, group = interaction(time, group), fill = group, color = group),
      width = 0.028,
      outlier.shape = NA,
      alpha = 0.55,
      linewidth = 0.5
    ) +
    geom_line(
      data = plot_data,
      aes(x = point_x, y = value, group = id, color = group),
      alpha = 0.3,
      linewidth = 0.3
    ) +
    geom_point(
      data = plot_data,
      aes(x = point_x, y = value, color = group),
      alpha = 0.4,
      size = 1.1,
      position = position_jitter(width = 0.008, height = 0)
    ) +
    geom_line(
      data = summary_df,
      aes(x = summary_x, y = summary_value, group = group, color = group),
      linewidth = 1.25,
      alpha = 0.95
    ) +
    scale_x_continuous(
      breaks = c(1, 2),
      labels = c("Pre", "Post"),
      limits = c(0.5, 2.5)
    ) +
    scale_fill_manual(values = group_colors, name = "Group") +
    scale_color_manual(values = group_colors, guide = "none") +
    labs(
      title = title_text,
      x = "Time",
      y = "Gradient score"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11)
    ) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.65, color = NA)))
}

# Save raincloud parcel plots for parcels 100 and 377 across gradients 1-3.
save_parcel_plots <- function(parcels = c(100, 377), output_dir = parcel_plot_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  gradient_specs <- tribble(
    ~gradient_col, ~gradient_label,
    "g1", "Gradient 1",
    "g2", "Gradient 2",
    "g3", "Gradient 3"
  )
  
  plots <- purrr::map(parcels, function(parcel_number) {
    parcel_data <- build_parcel_data(parcel_number)
    
    parcel_plots <- purrr::pmap(
      gradient_specs,
      function(gradient_col, gradient_label) {
        title_text <- paste(gradient_label, "Parcel", parcel_number)
        file_stem <- paste0(
          "eccentricity_parcel_",
          parcel_number,
          "_",
          tolower(gsub(" ", "_", gradient_label)),
          "_raincloud"
        )
        
        plot_object <- make_parcel_plot(parcel_data, gradient_col, title_text)
        
        save_plot_all_formats(
          plot_object = plot_object,
          file_stem = file_stem,
          output_dir = output_dir,
          width = 6.5,
          height = 6.5,
          dpi = 300
        )
        
        plot_object
      }
    )
    
    names(parcel_plots) <- gradient_specs$gradient_col
    parcel_plots
  })
  
  names(plots) <- paste0("parcel_", parcels)
  plots
}

read_process_gradients <- function(folder_path, participants = NULL) {
  session <- stringr::str_replace(folder_path, ".*rs_?", "rs")
  pc_read_session_gradient_matrices(session = session, data_root = figure_paths$data_root)
}

load_all_resting_state_gradients <- function(
    num_scans = 6,
    common_folder_path = NULL
) {
  pc_read_all_session_gradient_matrices(
    sessions = paste0("rs", seq_len(num_scans)),
    data_root = figure_paths$data_root
  )
}

process_eccentricity_data <- function(data, scan_label, parcels = c(100, 377)) {
  eccentricity_df <- tibble(ID = data$participants)
  
  for (parcel_number in parcels) {
    parcel_scores <- sapply(data$data, function(participant_data) {
      if (!is.null(participant_data) && ncol(participant_data) >= 3 && nrow(participant_data) >= parcel_number) {
        sqrt(
          participant_data[parcel_number, 1]^2 +
            participant_data[parcel_number, 2]^2 +
            participant_data[parcel_number, 3]^2
        )
      } else {
        NA_real_
      }
    })
    
    eccentricity_df[[paste0("parcel_", parcel_number)]] <- parcel_scores
  }
  
  eccentricity_df$scan <- scan_label
  eccentricity_df
}

build_six_rs_eccentricity_data <- function(
    parcels = c(100, 377),
    iv_path = file.path(figure_paths$data_root, "mean_FD.csv")
) {
  rs_data_list <- load_all_resting_state_gradients()
  scan_labels <- paste0("rest", seq_along(rs_data_list))
  
  eccentricity_data <- purrr::map2_dfr(
    rs_data_list,
    scan_labels,
    ~ process_eccentricity_data(.x, .y, parcels = parcels)
  )
  
  iv_data <- read.csv(iv_path) %>%
    select(ID, rest1, rest2, rest3, rest4, rest5, rest6, stress) %>%
    pivot_longer(
      cols = c(rest1, rest2, rest3, rest4, rest5, rest6),
      names_to = "scan",
      values_to = "headmotion"
    ) %>%
    mutate(scan = factor(scan, levels = paste0("rest", 1:6)))
  
  eccentricity_data %>%
    inner_join(iv_data, by = c("ID", "scan")) %>%
    pivot_longer(
      cols = starts_with("parcel_"),
      names_to = "parcel",
      values_to = "eccentricity"
    ) %>%
    mutate(
      stress = factor(as.character(stress), levels = c("0", "1"), labels = c("Control", "Stress")),
      parcel = factor(parcel, levels = paste0("parcel_", parcels)),
      parcel_label = recode(
        parcel,
        parcel_100 = "Parcel 100",
        parcel_377 = "Parcel 377"
      )
    )
}


# =========================
# Results figures: Dispersion
# =========================
#
# This section generates two dispersion figure sets:
#
#   Figure 2C:
#     Six-resting-state between-network distance plot
#
#   Figure 2D:
#     Network centroid ridge plots
#


# -------------------------------------------------------------------------
# Figure 2C:
# Six-resting-state between-network dispersion plot
# -------------------------------------------------------------------------

build_all_rs_gradient_dataframe <- function(
    iv_path = file.path(figure_paths$data_root, "mean_FD.csv"),
    parcel_labels_path = file.path(figure_paths$data_root, "yeo_networks_17_reordered_400_complete_labels.csv")
) {
  rs_data_list <- load_all_resting_state_gradients()
  
  flattened_data <- purrr::imap(rs_data_list, function(scan_data, scan_index) {
    purrr::imap_dfr(scan_data$data, function(participant_gradient_data, participant_index) {
      if (is.null(participant_gradient_data) || nrow(participant_gradient_data) == 0 || ncol(participant_gradient_data) < 3) {
        return(tibble())
      }
      
      tibble(
        NECOS_ID = scan_data$participants[[participant_index]],
        scan = stringr::str_replace(as.character(scan_index), "^rs", ""),
        parcel = seq_len(nrow(participant_gradient_data)),
        V1 = participant_gradient_data[, 1],
        V2 = participant_gradient_data[, 2],
        V3 = participant_gradient_data[, 3]
      )
    })
  })
  
  final_data <- bind_rows(flattened_data)
  
  iv_long <- read.csv(iv_path) %>%
    rename(NECOS_ID = ID) %>%
    pivot_longer(
      cols = starts_with("rest"),
      names_to = "scan",
      names_prefix = "rest",
      values_to = "headmotion"
    )
  
  parcel_labels <- read.csv(parcel_labels_path) %>%
    select(X7networks_order, X7networks) %>%
    rename(parcel = X7networks_order, networks = X7networks)
  
  final_data %>%
    inner_join(
      iv_long %>% select(NECOS_ID, scan, stress, headmotion),
      by = c("NECOS_ID", "scan")
    ) %>%
    inner_join(parcel_labels, by = "parcel")
}

build_between_network_six_rs_data <- function() {
  my_networks <- c("Default", "Control", "Salience")
  
  mean_centroids_all_rs <- build_all_rs_gradient_dataframe() %>%
    group_by(NECOS_ID, networks, scan, stress, headmotion) %>%
    summarise(
      mean_V1 = mean(V1, na.rm = TRUE),
      mean_V2 = mean(V2, na.rm = TRUE),
      mean_V3 = mean(V3, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(networks = if_else(networks == "Salience / Ventral Attn", "Salience", as.character(networks))) %>%
    filter(networks %in% my_networks)
  
  wide_data <- mean_centroids_all_rs %>%
    pivot_wider(names_from = networks, values_from = c(mean_V1, mean_V2, mean_V3))
  
  wide_data %>%
    rowwise() %>%
    mutate(
      distance_default_salience = sqrt(sum(
        (c(mean_V1_Default, mean_V2_Default, mean_V3_Default) -
           c(mean_V1_Salience, mean_V2_Salience, mean_V3_Salience))^2
      )),
      distance_control_salience = sqrt(sum(
        (c(mean_V1_Control, mean_V2_Control, mean_V3_Control) -
           c(mean_V1_Salience, mean_V2_Salience, mean_V3_Salience))^2
      ))
    ) %>%
    ungroup() %>%
    transmute(
      NECOS_ID,
      scan = paste0("rest", scan),
      stress = factor(as.character(stress), levels = c("0", "1"), labels = c("Control", "Stress")),
      headmotion,
      distance_default_salience,
      distance_control_salience
    ) %>%
    pivot_longer(
      cols = c(distance_default_salience, distance_control_salience),
      names_to = "comparison",
      values_to = "distance"
    ) %>%
    mutate(
      comparison = recode(
        comparison,
        distance_default_salience = "Default-Salience",
        distance_control_salience = "Control-Salience"
      )
    )
}

# Create the six-resting-state between-network distance figure.
make_between_network_six_rs_plot <- function(data) {
  physical_distances <- c(
    rest1 = -35,
    rest2 = -20,
    rest3 = 35,
    rest4 = 55,
    rest5 = 90,
    rest6 = 105
  )
  
  plot_data <- data %>%
    mutate(
      scan = factor(scan, levels = names(physical_distances)),
      position = unname(physical_distances[as.character(scan)])
    )
  
  summary_data <- plot_data %>%
    group_by(comparison, stress, scan, position) %>%
    summarise(
      mean_value = mean(distance, na.rm = TRUE),
      se_value = mean_se(distance)$ymax - mean_se(distance)$y,
      .groups = "drop"
    )
  
  ggplot(summary_data, aes(x = position, y = mean_value, color = stress, group = stress)) +
    geom_line(
      linewidth = 1.2,
      alpha = 0.85,
      linetype = "dotted"
    ) +
    geom_point(
      size = 3.2
    ) +
    geom_errorbar(
      aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
      width = 4,
      linewidth = 0.6
    ) +
    facet_wrap(~ comparison, ncol = 1, scales = "free_y") +
    scale_color_manual(values = group_colors, name = "Group") +
    scale_x_continuous(
      breaks = unname(physical_distances),
      labels = paste("Scan", 1:6)
    ) +
    labs(
      title = "Between-Network Distances Across Six Resting States",
      x = "Resting-state scan",
      y = "Mean Euclidean distance"
    ) +
    theme_classic() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 12),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 11)
    )
}

save_between_network_six_rs_plot <- function(output_dir = between_network_six_rs_plot_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  between_network_data <- build_between_network_six_rs_data()
  plot_object <- make_between_network_six_rs_plot(between_network_data)
  
  save_plot_png_svg(
    plot_object,
    file_stem = "dispersion_six_resting_states_between_network_distances",
    output_dir = output_dir,
    width = 10,
    height = 8
  )
  
  plot_object
}


# -------------------------------------------------------------------------
# Figure 2D:
# Network centroid ridge plots
# -------------------------------------------------------------------------

# Ridge plots of network centroid distributions for Gradient 1, Gradient 2,
# and Gradient 3.
build_network_centroids <- function(
    raw_data_path = file.path(dispersion_between_paths$tables_dir, "dispersion_gradients_parcel_labels_iv.csv")
) {
  raw_data <- read.csv(raw_data_path)
  
  raw_data %>%
    group_by(NECOS_ID, X7networks, scan, headmotion) %>%
    summarise(
      centroid_V1 = mean(V1, na.rm = TRUE),
      centroid_V2 = mean(V2, na.rm = TRUE),
      centroid_V3 = mean(V3, na.rm = TRUE),
      stress = first(stress),
      .groups = "drop"
    ) %>%
    filter(X7networks %in% c("Control", "Default", "Salience / Ventral Attn")) %>%
    mutate(
      group = factor(
        if_else(as.character(stress) %in% c("1", "Stress"), "Stress", "Control"),
        levels = c("Control", "Stress")
      ),
      time = factor(scan, levels = c("rest2", "rest3"), labels = c("Pre", "Post")),
      network = factor(
        X7networks,
        levels = c("Control", "Default", "Salience / Ventral Attn")
      )
    )
}

save_plot_png_svg <- function(plot_object, file_stem, output_dir, width = 8, height = 6) {
  ggsave(
    filename = file.path(output_dir, paste0(file_stem, ".png")),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(output_dir, paste0(file_stem, ".svg")),
    plot = plot_object,
    width = width,
    height = height,
    device = "svg"
  )
}

save_plot_all_formats <- function(plot_object, file_stem, output_dir, width = 8, height = 6, dpi = 300) {
  ggsave(
    filename = file.path(output_dir, paste0(file_stem, ".png")),
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi
  )
  
  ggsave(
    filename = file.path(output_dir, paste0(file_stem, ".pdf")),
    plot = plot_object,
    width = width,
    height = height
  )
  
  ggsave(
    filename = file.path(output_dir, paste0(file_stem, ".svg")),
    plot = plot_object,
    width = width,
    height = height,
    device = "svg"
  )
}

make_network_ridge_plot <- function(data, gradient_col, title_text) {
  network_colors <- c(
    "Control" = "#E69422",
    "Default" = "#CD3E4E",
    "Salience / Ventral Attn" = "#C43AFA"
  )
  
  ggplot(data, aes(x = .data[[gradient_col]], y = group, fill = network, alpha = time)) +
    geom_density_ridges(aes(height = after_stat(density)), color = NA, scale = 0.75) +
    scale_fill_manual(values = network_colors, name = "Network") +
    scale_alpha_manual(values = c(Pre = 0.8, Post = 0.5), name = "Time") +
    labs(
      title = title_text,
      x = "Gradient score",
      y = "Group"
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_classic() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 12),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9)
    )
}

save_network_ridge_plots <- function(output_dir = ridge_plot_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  centroids <- build_network_centroids()
  
  gradient_specs <- tribble(
    ~gradient_col, ~gradient_label, ~file_stem,
    "centroid_V1", "Gradient 1", "ridge_plot_g1",
    "centroid_V2", "Gradient 2", "ridge_plot_g2",
    "centroid_V3", "Gradient 3", "ridge_plot_g3"
  )
  
  plots <- purrr::pmap(
    gradient_specs,
    function(gradient_col, gradient_label, file_stem) {
      plot_object <- make_network_ridge_plot(centroids, gradient_col, gradient_label)
      save_plot_png_svg(plot_object, file_stem, output_dir)
      plot_object
    }
  )
  
  names(plots) <- c("g1", "g2", "g3")
  plots
}
# =============================
# Peripheral Measures Figures
# =============================
#
# This section generates three peripheral-measures figure sets:
#
#   Peripheral Figure 1:
#     Best brain outcome per peripheral marker
#     - grid plot
#     - single-panel versions of the same best-outcome plots
#
#   Peripheral Figure 2:
#     Stress-group-only Spearman correlation plots
#
#   Peripheral Figure 3:
#     Epoch-level physiological timeline plots
#


# -------------------------------------------------------------------------
# Shared data preparation for Peripheral Figures 1 and 2
# -------------------------------------------------------------------------

build_physio_brain_dataframe <- function(
    processed_dir = peripheral_paths$tables_dir,
    final_dir = peripheral_paths$models_dir,
    data_dir = figure_paths$data_root
) {
  eccentricity_delta <- readr::read_csv(file.path(processed_dir, "peripheral_eccentricity_delta_base.csv"), show_col_types = FALSE)
  dispersion_delta <- readr::read_csv(file.path(processed_dir, "peripheral_dispersion_delta_base.csv"), show_col_types = FALSE)
  physio_hr_hrv_delta <- readr::read_csv(file.path(data_dir, "physio_hr_hrv_delta.csv"), show_col_types = FALSE)
  physio_baselines <- readr::read_csv(file.path(data_dir, "physio_baselines.csv"), show_col_types = FALSE)
  endo_delta <- readr::read_csv(file.path(data_dir, "endo_delta.csv"), show_col_types = FALSE)
  endo_baselines <- readr::read_csv(file.path(data_dir, "endo_baselines.csv"), show_col_types = FALSE)
  stai_delta <- readr::read_csv(file.path(data_dir, "stai_delta.csv"), show_col_types = FALSE)
  stai_baselines <- readr::read_csv(file.path(data_dir, "stai_baselines.csv"), show_col_types = FALSE)
  vas_delta <- readr::read_csv(file.path(data_dir, "vas_delta.csv"), show_col_types = FALSE)
  vas_baselines <- readr::read_csv(file.path(data_dir, "vas_baselines.csv"), show_col_types = FALSE)
  results_tbl <- readr::read_csv(file.path(final_dir, "peripheral_model_results.csv"), show_col_types = FALSE)
  
  brain_df <- bind_rows(
    eccentricity_delta %>% mutate(brain_family = "Eccentricity"),
    dispersion_delta %>% mutate(brain_family = "Dispersion")
  ) %>%
    left_join(select(physio_hr_hrv_delta, NECOS_ID, HR_delta, HRV_delta), by = "NECOS_ID") %>%
    left_join(select(physio_baselines, NECOS_ID, log_base_HR, sqrt_base_HRV), by = "NECOS_ID") %>%
    left_join(select(endo_delta, NECOS_ID, SaCort_delta, asinh_ACTH_delta), by = "NECOS_ID") %>%
    left_join(select(endo_baselines, NECOS_ID, log_base_saliva_cortisol, log_base_ACTH), by = "NECOS_ID") %>%
    left_join(select(stai_delta, NECOS_ID, STAI_delta), by = "NECOS_ID") %>%
    left_join(select(stai_baselines, NECOS_ID, base_STAI), by = "NECOS_ID") %>%
    left_join(select(vas_delta, NECOS_ID, VAS_delta), by = "NECOS_ID") %>%
    left_join(select(vas_baselines, NECOS_ID, sqrt_base_VAS), by = "NECOS_ID") %>%
    mutate(Group = factor(pretty_group(Group), levels = c("Control", "Stress")))
  
  outcome_map <- c(
    delta_parcel_100 = "Insula eccentricity",
    delta_parcel_377 = "vPFC eccentricity",
    delta_distance_control_salience = "SN-FPN dispersion",
    delta_distance_default_salience = "SN-DMN dispersion"
  )
  
  outcome_baseline_map <- c(
    delta_parcel_100 = "parcel_100_rest2",
    delta_parcel_377 = "parcel_377_rest2",
    delta_distance_control_salience = "distance_control_salience_rest2",
    delta_distance_default_salience = "distance_default_salience_rest2"
  )
  
  predictor_map <- c(
    VAS_delta = "Delta VAS",
    STAI_delta = "Delta STAI",
    HR_delta = "Delta heart rate",
    HRV_delta = "Delta HRV/RMSSD",
    SaCort_delta = "Delta saliva cortisol",
    asinh_ACTH_delta = "Delta ACTH"
  )
  
  predictor_baseline_map <- c(
    VAS_delta = "sqrt_base_VAS",
    STAI_delta = "base_STAI",
    HR_delta = "log_base_HR",
    HRV_delta = "sqrt_base_HRV",
    SaCort_delta = "log_base_saliva_cortisol",
    asinh_ACTH_delta = "log_base_ACTH"
  )
  
  outcome_map <- outcome_map[names(outcome_map) %in% names(brain_df)]
  outcome_baseline_map <- outcome_baseline_map[names(outcome_map)]
  predictor_map <- predictor_map[names(predictor_map) %in% names(brain_df)]
  predictor_baseline_map <- predictor_baseline_map[names(predictor_map)]
  
  list(
    brain_df = brain_df,
    results_tbl = results_tbl,
    outcome_map = outcome_map,
    outcome_baseline_map = outcome_baseline_map,
    predictor_map = predictor_map,
    predictor_baseline_map = predictor_baseline_map
  )
}

build_physio_panel_data <- function(
    data_bundle,
    outcome_col,
    predictor_col,
    plot_mode = "raw_delta"
) {
  dat <- data_bundle$brain_df
  motion_col <- "log_mean_MFD_scaled"
  y0_col <- data_bundle$outcome_baseline_map[[outcome_col]]
  x0_col <- data_bundle$predictor_baseline_map[[predictor_col]]
  outcome_label <- unname(data_bundle$outcome_map[[outcome_col]])
  predictor_label <- unname(data_bundle$predictor_map[[predictor_col]])
  
  needed_cols <- c("NECOS_ID", "Group", "brain_family", motion_col, outcome_col, predictor_col)
  if (plot_mode == "adjusted") {
    needed_cols <- c(needed_cols, y0_col, x0_col)
  }
  
  panel <- dat %>%
    select(any_of(needed_cols)) %>%
    drop_na()
  
  if (nrow(panel) < 5) {
    return(NULL)
  }
  
  if (plot_mode == "raw_delta") {
    panel <- panel %>% mutate(x_plot = .data[[predictor_col]], y_plot = .data[[outcome_col]])
  } else {
    y_fit <- lm(as.formula(paste0(outcome_col, " ~ ", y0_col, " + ", x0_col, " + ", motion_col)), data = panel)
    x_fit <- lm(as.formula(paste0(predictor_col, " ~ ", x0_col, " + ", motion_col)), data = panel)
    panel <- panel %>% mutate(x_plot = resid(x_fit), y_plot = resid(y_fit))
  }
  
  panel %>%
    transmute(
      NECOS_ID,
      Group,
      brain_family,
      outcome_col,
      predictor_col,
      Outcome = outcome_label,
      Predictor = predictor_label,
      x_plot,
      y_plot,
      n_panel = dplyr::n()
    )
}

build_all_physio_panels <- function(data_bundle, plot_mode = "raw_delta") {
  tidyr::crossing(
    outcome_col = names(data_bundle$outcome_map),
    predictor_col = names(data_bundle$predictor_map)
  ) %>%
    mutate(data = purrr::map2(outcome_col, predictor_col, ~ build_physio_panel_data(data_bundle, .x, .y, plot_mode = plot_mode))) %>%
    pull(data) %>%
    purrr::compact() %>%
    bind_rows() %>%
    mutate(
      Outcome = factor(Outcome, levels = unname(data_bundle$outcome_map)),
      Predictor = factor(Predictor, levels = unname(data_bundle$predictor_map)),
      Group = factor(Group, levels = c("Control", "Stress"))
    )
}

build_physio_best_pairs <- function(data_bundle, selection_rule = "lowest_p") {
  best_pairs <- data_bundle$results_tbl %>%
    filter(
      predictor %in% names(data_bundle$predictor_map),
      DV %in% names(data_bundle$outcome_map),
      !is.na(p)
    )
  
  if (selection_rule == "highest_p") {
    best_pairs <- best_pairs %>% group_by(predictor) %>% slice_max(order_by = p, n = 1, with_ties = FALSE) %>% ungroup()
  } else {
    best_pairs <- best_pairs %>% group_by(predictor) %>% slice_min(order_by = p, n = 1, with_ties = FALSE) %>% ungroup()
  }
  
  best_pairs %>%
    mutate(
      Outcome = unname(data_bundle$outcome_map[DV]),
      Predictor = unname(data_bundle$predictor_map[predictor]),
      predictor = factor(predictor, levels = names(data_bundle$predictor_map))
    ) %>%
    arrange(predictor)
}


# -------------------------------------------------------------------------
# Peripheral Figure 1A:
# Best brain outcome per peripheral marker — grid version
# -------------------------------------------------------------------------

make_physio_best_outcome_grid <- function(data_bundle, plot_mode = "raw_delta", selection_rule = "lowest_p") {
  best_pairs <- build_physio_best_pairs(data_bundle, selection_rule = selection_rule)
  
  panel_df_best <- purrr::map2_dfr(best_pairs$DV, as.character(best_pairs$predictor), ~ build_physio_panel_data(data_bundle, .x, .y, plot_mode = plot_mode)) %>%
    mutate(Group = factor(Group, levels = c("Control", "Stress")))
  
  ann_df_best <- best_pairs %>%
    transmute(
      outcome_col = DV,
      predictor_col = as.character(predictor),
      Predictor,
      label = paste0(
        "n = ", n, "\n",
        "p = ", fmt_p(p), "\n",
        "p corrected = ", fmt_p(p_fdr_within_family)
      )
    )
  
  plot_object <- ggplot(panel_df_best, aes(x = x_plot, y = y_plot, color = Group)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
    geom_point(alpha = 0.8, size = 2.6) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
    geom_text(
      data = ann_df_best,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.03,
      vjust = 1.15,
      size = 3.8,
      color = "black",
      lineheight = 0.95
    ) +
    facet_wrap(~ Predictor, nrow = 1, scales = "free") +
    scale_color_manual(values = group_colors, name = "Group") +
    labs(
      title = paste0("Best brain outcome per peripheral marker (", selection_rule, ", ", plot_mode, ")"),
      x = ifelse(plot_mode == "raw_delta", "Delta peripheral marker", "Covariate-adjusted delta peripheral marker"),
      y = ifelse(plot_mode == "raw_delta", "Delta brain measure", "Covariate-adjusted delta brain measure")
    ) +
    theme_bw(base_size = 14) +
    theme(
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "top",
      legend.text = element_text(size = 11),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 10),
      plot.title = element_text(face = "bold", size = 16),
      aspect.ratio = 1
    )
  
  list(plot = plot_object, best_pairs = best_pairs, panel_df_best = panel_df_best)
}

save_physio_best_outcome_grid <- function(output_dir = physio_plot_dir, plot_mode = "raw_delta", selection_rule = "lowest_p") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  data_bundle <- build_physio_brain_dataframe()
  grid_bundle <- make_physio_best_outcome_grid(data_bundle, plot_mode = plot_mode, selection_rule = selection_rule)
  
  save_plot_all_formats(
    plot_object = grid_bundle$plot,
    file_stem = paste0("peripheral_best_brain_outcome_", selection_rule, "_", plot_mode),
    output_dir = output_dir,
    width = 28,
    height = 8
  )
  
  grid_bundle
}


# -------------------------------------------------------------------------
# Peripheral Figure 1B:
# Best brain outcome per peripheral marker — single-panel versions
# -------------------------------------------------------------------------

make_physio_single_best_plot <- function(df_sub, label_text, predictor_title, outcome_title, plot_mode = "raw_delta") {
  ggplot(df_sub, aes(x = x_plot, y = y_plot, color = Group)) +
    geom_hline(yintercept = 0, linewidth = 0.7, color = "grey78") +
    geom_vline(xintercept = 0, linewidth = 0.7, color = "grey78") +
    geom_point(alpha = 0.85, size = 3.6) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.8) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = label_text,
      hjust = -0.05,
      vjust = 1.1,
      size = 6,
      color = "black",
      lineheight = 1
    ) +
    scale_color_manual(values = group_colors, name = "Group") +
    labs(
      title = predictor_title,
      x = ifelse(plot_mode == "raw_delta", predictor_title, paste0("Covariate-adjusted ", predictor_title)),
      y = outcome_title
    ) +
    theme_bw(base_size = 20) +
    theme(
      plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
      legend.position = "top",
      legend.text = element_text(size = 18),
      axis.title = element_text(face = "bold", size = 20),
      axis.text = element_text(size = 17),
      aspect.ratio = 1,
      panel.border = element_rect(linewidth = 1.1),
      axis.line = element_line(linewidth = 0.9)
    )
}

save_physio_single_best_plots <- function(output_dir = physio_plot_dir, plot_mode = "raw_delta", selection_rule = "lowest_p") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  data_bundle <- build_physio_brain_dataframe()
  best_pairs <- build_physio_best_pairs(data_bundle, selection_rule = selection_rule)
  
  plot_list <- purrr::pmap(
    list(best_pairs$DV, as.character(best_pairs$predictor), best_pairs$Predictor, best_pairs$Outcome, best_pairs$n, best_pairs$p, best_pairs$p_fdr_within_family),
    function(outcome_col, predictor_col, predictor_title, outcome_title, n, p, p_fdr_within_family) {
      df_sub <- build_physio_panel_data(data_bundle, outcome_col, predictor_col, plot_mode = plot_mode)
      if (is.null(df_sub)) {
        return(NULL)
      }
      
      label_text <- paste0(
        "n = ", n, "\n",
        "p = ", fmt_p(p), "\n",
        "p corrected = ", fmt_p(p_fdr_within_family)
      )
      
      plot_object <- make_physio_single_best_plot(df_sub, label_text, predictor_title, outcome_title, plot_mode = plot_mode)
      file_stub <- paste0("best_outcome_", gsub("[^A-Za-z0-9]+", "_", predictor_col), "_", plot_mode)
      save_plot_all_formats(plot_object, file_stub, output_dir, width = 8, height = 8)
      plot_object
    }
  )
  
  purrr::compact(plot_list)
}


# Single panel for a chosen outcome and peripheral marker (same layout as the best-outcome panels)
save_physio_selected_pair_plot <- function(outcome_col, predictor_col, output_dir = physio_plot_dir, plot_mode = "raw_delta") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  data_bundle <- build_physio_brain_dataframe()
  pair <- data_bundle$results_tbl %>% filter(DV == outcome_col, predictor == predictor_col)
  df_sub <- build_physio_panel_data(data_bundle, outcome_col, predictor_col, plot_mode = plot_mode)
  
  label_text <- paste0(
    "n = ", pair$n, "\n",
    "p = ", fmt_p(pair$p), "\n",
    "p corrected = ", fmt_p(pair$p_fdr_within_family)
  )
  plot_object <- make_physio_single_best_plot(
    df_sub, label_text,
    predictor_title = unname(data_bundle$predictor_map[predictor_col]),
    outcome_title = unname(data_bundle$outcome_map[outcome_col]),
    plot_mode = plot_mode
  )
  save_plot_all_formats(plot_object, paste0("selected_pair_", predictor_col, "_", outcome_col, "_", plot_mode), output_dir, width = 8, height = 8)
  plot_object
}


# -------------------------------------------------------------------------
# Peripheral Figure 2:
# Stress-group-only Spearman correlation plots
# -------------------------------------------------------------------------

build_stress_only_spearman_stats <- function(panel_df, data_bundle) {
  stress_only_df <- panel_df %>% filter(Group == "Stress") %>% droplevels()
  
  spearman_stats <- stress_only_df %>%
    group_by(outcome_col, predictor_col, Outcome, Predictor, brain_family) %>%
    summarise(
      n = sum(complete.cases(x_plot, y_plot)),
      rho = if (n >= 3) suppressWarnings(cor(x_plot, y_plot, method = "spearman", use = "complete.obs")) else NA_real_,
      p = if (n >= 3) suppressWarnings(tryCatch(cor.test(x_plot, y_plot, method = "spearman", exact = FALSE)$p.value, error = function(e) NA_real_)) else NA_real_,
      .groups = "drop"
    ) %>%
    group_by(brain_family) %>%
    mutate(p_fdr_within_family = p.adjust(p, method = "fdr")) %>%
    ungroup() %>%
    mutate(
      rho_label = if_else(is.na(rho), "NA", sprintf("%.2f", rho)),
      label = paste0(
        "n=", n, "\n",
        "rho=", rho_label, "\n",
        "p=", fmt_p(p), "\n",
        "p corrected=", fmt_p(p_fdr_within_family)
      )
    )
  
  list(stress_only_df = stress_only_df, spearman_stats = spearman_stats)
}

make_physio_stress_only_grid <- function(df, ann, title_text, plot_mode = "raw_delta") {
  ggplot(df, aes(x = x_plot, y = y_plot)) +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey85") +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "grey85") +
    geom_point(alpha = 0.8, size = 1.8, color = group_colors[["Stress"]]) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = group_colors[["Stress"]]) +
    geom_text(
      data = ann,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.05,
      vjust = 1.1,
      size = 2.9,
      color = "black"
    ) +
    facet_grid(rows = vars(Outcome), cols = vars(Predictor), scales = "free") +
    labs(
      title = title_text,
      x = ifelse(plot_mode == "raw_delta", "Delta peripheral marker", "Covariate-adjusted delta peripheral marker"),
      y = ifelse(plot_mode == "raw_delta", "Delta brain measure", "Covariate-adjusted delta brain measure")
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(0.7, "lines"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 14),
      aspect.ratio = 1
    )
}

save_physio_stress_only_spearman_plots <- function(output_dir = physio_plot_dir, plot_mode = "raw_delta") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  data_bundle <- build_physio_brain_dataframe()
  panel_df <- build_all_physio_panels(data_bundle, plot_mode = plot_mode)
  spearman_bundle <- build_stress_only_spearman_stats(panel_df, data_bundle)
  
  readr::write_csv(
    spearman_bundle$spearman_stats,
    file.path(peripheral_paths$models_dir, paste0("stress_group_only_spearman_results_", plot_mode, ".csv"))
  )
  
  full_plot <- make_physio_stress_only_grid(
    spearman_bundle$stress_only_df,
    spearman_bundle$spearman_stats,
    ifelse(
      plot_mode == "raw_delta",
      "Stress group only: delta brain outcomes vs delta peripheral markers",
      "Stress group only: covariate-adjusted brain outcomes vs peripheral markers"
    ),
    plot_mode = plot_mode
  )
  save_plot_all_formats(full_plot, paste0("stress_group_only_brain_peripheral_grid_spearman_", plot_mode), output_dir, width = 18, height = 15)
  
  ecc_plot <- NULL
  disp_plot <- NULL
  stress_ecc <- spearman_bundle$stress_only_df %>% filter(brain_family == "Eccentricity")
  ann_ecc <- spearman_bundle$spearman_stats %>% filter(brain_family == "Eccentricity")
  if (nrow(stress_ecc) > 0) {
    ecc_plot <- make_physio_stress_only_grid(stress_ecc, ann_ecc, paste0("Stress group only: Eccentricity panels (", plot_mode, ")"), plot_mode = plot_mode)
    save_plot_all_formats(ecc_plot, paste0("stress_group_only_eccentricity_grid_spearman_", plot_mode), output_dir, width = 18, height = 15)
  }
  
  stress_disp <- spearman_bundle$stress_only_df %>% filter(brain_family == "Dispersion")
  ann_disp <- spearman_bundle$spearman_stats %>% filter(brain_family == "Dispersion")
  if (nrow(stress_disp) > 0) {
    disp_plot <- make_physio_stress_only_grid(stress_disp, ann_disp, paste0("Stress group only: Dispersion panels (", plot_mode, ")"), plot_mode = plot_mode)
    save_plot_all_formats(disp_plot, paste0("stress_group_only_dispersion_grid_spearman_", plot_mode), output_dir, width = 18, height = 15)
  }
  
  list(full_plot = full_plot, eccentricity_plot = ecc_plot, dispersion_plot = disp_plot)
}


# -------------------------------------------------------------------------
# Peripheral Figure 3:
# Epoch-level physiological timeline plots
# -------------------------------------------------------------------------

standardize_timeline_stress <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("0", "control", "Control", "CG") ~ "control",
    x %in% c("1", "stress", "Stress", "SG") ~ "stress",
    TRUE ~ tolower(x)
  )
}

standard_error <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) {
    return(NA_real_)
  }
  stats::sd(x) / sqrt(length(x))
}

build_physio_epoch_timeline_data <- function(path_wrapper = NULL) {
  timeline_path <- file.path(figure_paths$data_root, "peripheral_timeline_long.csv")
  
  if (!file.exists(timeline_path)) {
    message("Skipping epoch-level peripheral timeline figure because `peripheral_timeline_long.csv` is missing.")
    return(NULL)
  }
  
  timeline_df <- readr::read_csv(timeline_path, show_col_types = FALSE)
  
  endo_time_map <- tibble(
    Timepoint = as.character(c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)),
    label = c(
      "T3\n-45", "T4\n-30", "T5\n-15", "T6\n+5", "T7\n+15", "T8\n+25",
      "T9\n+45", "T10\n+60", "T11\n+80", "T12\n+95", "T13\n+110", "T14\n+130"
    ),
    tsst_min = c(-45, -30, -15, 5, 15, 25, 45, 60, 80, 95, 110, 130)
  )
  
  auto_time_map <- tibble(
    Timepoint = c(
      "01_rest1a", "02_rest1b", "03_rest2a", "04_rest2b",
      "05_TSST_ant", "06_TSST_int", "07_TSST_mat",
      "08_rest3a", "09_rest3b", "10_rest4a", "11_rest4b",
      "12_rest5a", "13_rest5b", "14_rest6a", "15_rest6b"
    ),
    label = c(
      "rest1a\n-35", "rest1b\n-30", "rest2a\n-20", "rest2b\n-15",
      "TSST1\n0", "TSST2\n+10", "TSST3\n+15",
      "rest3a\n+40", "rest3b\n+45", "rest4a\n+60", "rest4b\n+65",
      "rest5a\n+100", "rest5b\n+105", "rest6a\n+120", "rest6b\n+125"
    ),
    tsst_min = c(-35, -30, -20, -15, 0, 10, 15, 40, 45, 60, 65, 100, 105, 120, 125)
  )
  
  endo_epoch_map <- tibble(
    Timepoint = c("3", "4", "4", "5", "6", "7", "8", "9", "9", "10", "11", "12", "12", "13"),
    epoch = c("rest1", "rest1", "rest2", "rest2", "TSST1", "TSST3", "rest3", "rest3", "rest4", "rest4", "rest5", "rest5", "rest6", "rest6"),
    epoch_tsst_min = c(-35, -35, -20, -20, 0, 15, 35, 35, 55, 55, 90, 90, 115, 115),
    epoch_label = c(
      "rest1\n-35", "rest1\n-35", "rest2\n-20", "rest2\n-20", "TSST1\n+5", "TSST3\n+15",
      "rest3\n+35", "rest3\n+35", "rest4\n+55", "rest4\n+55", "rest5\n+90", "rest5\n+90", "rest6\n+115", "rest6\n+115"
    )
  )
  
  auto_epoch_map <- tibble(
    Timepoint = c(
      "01_rest1a", "02_rest1b", "03_rest2a", "04_rest2b",
      "05_TSST_ant", "06_TSST_int", "07_TSST_mat",
      "08_rest3a", "09_rest3b", "10_rest4a", "11_rest4b",
      "12_rest5a", "13_rest5b", "14_rest6a", "15_rest6b"
    ),
    epoch = c(
      "rest1", "rest1", "rest2", "rest2", "TSST1", "TSST2", "TSST3",
      "rest3", "rest3", "rest4", "rest4", "rest5", "rest5", "rest6", "rest6"
    ),
    epoch_tsst_min = c(-35, -35, -20, -20, 0, 10, 15, 35, 35, 55, 55, 90, 90, 122.5, 122.5),
    epoch_label = c(
      "rest1\n-35", "rest1\n-35", "rest2\n-20", "rest2\n-20", "TSST1\n0", "TSST2\n+10", "TSST3\n+15",
      "rest3\n+35", "rest3\n+35", "rest4\n+55", "rest4\n+55", "rest5\n+90", "rest5\n+90", "rest6\n+123", "rest6\n+123"
    )
  )
  
  timeline_base <- timeline_df %>%
    mutate(
      ID = as.character(NECOS_ID),
      stress = standardize_timeline_stress(Group),
      measure = as.character(measure),
      value = as.numeric(value),
      modality = as.character(modality),
      Timepoint = as.character(Timepoint)
    )
  
  endo_long <- timeline_base %>%
    filter(measure %in% c("Speichel_Cortisol", "Plasma_ACTH")) %>%
    mutate(
      Timepoint = gsub("^T", "", as.character(Timepoint))
    ) %>%
    left_join(endo_time_map, by = "Timepoint") %>%
    filter(!is.na(tsst_min)) %>%
    filter(!is.na(value)) %>%
    mutate(outcome = dplyr::recode(measure, Speichel_Cortisol = "cortisol\nnmol/L", Plasma_ACTH = "ACTH\npmol/L"))
  
  stai_long <- timeline_base %>%
    filter(measure == "STAI") %>%
    mutate(
      Timepoint = gsub("^T", "", as.character(Timepoint))
    ) %>%
    left_join(endo_time_map, by = "Timepoint") %>%
    filter(!is.na(tsst_min), !is.na(value)) %>%
    transmute(ID, stress, Timepoint, label, tsst_min, outcome = "STAI\n(20-80)", value)
  
  vas_long <- timeline_base %>%
    filter(measure == "VAS") %>%
    mutate(
      Timepoint = gsub("^T", "", as.character(Timepoint))
    ) %>%
    left_join(endo_time_map, by = "Timepoint") %>%
    filter(!is.na(tsst_min), !is.na(value)) %>%
    transmute(ID, stress, Timepoint, label, tsst_min, outcome = "perceived stress VAS\n(0-100)", value)
  
  auto_long <- timeline_base %>%
    filter(measure %in% c("HR", "HRV")) %>%
    mutate(
      outcome = dplyr::recode(measure, HR = "heart rate\nbpm", HRV = "HRV RMSSD\nms")
    ) %>%
    transmute(ID, stress, Timepoint, outcome, value) %>%
    left_join(auto_time_map, by = "Timepoint") %>%
    filter(!is.na(tsst_min), !is.na(value)) %>%
    mutate(outcome = factor(outcome, levels = c("heart rate\nbpm", "HRV RMSSD\nms")))
  
  make_endo_epoch <- function(long_df) {
    long_df %>%
      select(ID, stress, outcome, Timepoint, value) %>%
      inner_join(endo_epoch_map, by = "Timepoint", relationship = "many-to-many") %>%
      group_by(ID, stress, outcome, epoch, epoch_label, epoch_tsst_min) %>%
      summarise(value = if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE), .groups = "drop")
  }
  
  epoch_summary <- function(epoch_long) {
    epoch_long %>%
      group_by(stress, outcome, epoch, epoch_label, epoch_tsst_min) %>%
      summarise(mean = mean(value, na.rm = TRUE), se = standard_error(value), n = sum(!is.na(value)), .groups = "drop") %>%
      rename(label = epoch_label, tsst_min = epoch_tsst_min)
  }
  
  endo_all_epoch_sum <- bind_rows(
    epoch_summary(make_endo_epoch(endo_long)),
    epoch_summary(make_endo_epoch(stai_long)),
    epoch_summary(make_endo_epoch(vas_long))
  ) %>%
    mutate(outcome = factor(outcome, levels = c("STAI\n(20-80)", "perceived stress VAS\n(0-100)", "cortisol\nnmol/L", "ACTH\npmol/L")))
  
  auto_epoch_long <- auto_long %>%
    select(ID, stress, outcome, Timepoint, value) %>%
    inner_join(auto_epoch_map, by = "Timepoint") %>%
    group_by(ID, stress, outcome, epoch, epoch_label, epoch_tsst_min) %>%
    summarise(value = if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE), .groups = "drop")
  
  auto_epoch_sum <- epoch_summary(auto_epoch_long) %>%
    mutate(outcome = factor(outcome, levels = c("heart rate\nbpm", "HRV RMSSD\nms")))
  
  list(endo_all_epoch_sum = endo_all_epoch_sum, auto_epoch_sum = auto_epoch_sum)
}

make_physio_stacked_plot <- function(df, xlims = NULL, base_size = 11, point_size = 2.2, line_width = 0.8, error_width = 2, x_text_size = 10) {
  ggplot(df, aes(x = tsst_min, y = mean, color = stress, group = stress)) +
    annotate("rect", xmin = 0, xmax = 20, ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "grey50") +
    annotate("rect", xmin = 20, xmax = 60, ymin = -Inf, ymax = Inf, alpha = 0.06, fill = "grey50") +
    geom_line(linewidth = line_width, linetype = "dotted") +
    geom_point(size = point_size) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = error_width, linewidth = 0.4) +
    facet_grid(outcome ~ ., scales = "free_y", switch = "y") +
    scale_x_continuous(
      breaks = df %>% distinct(tsst_min, label) %>% arrange(tsst_min) %>% pull(tsst_min),
      labels = df %>% distinct(tsst_min, label) %>% arrange(tsst_min) %>% pull(label),
      limits = xlims
    ) +
    scale_color_manual(values = timeline_group_colors, breaks = c("control", "stress")) +
    labs(x = "time point", y = NULL, color = NULL) +
    theme_classic(base_size = base_size) +
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0),
      axis.text.x = element_text(size = x_text_size),
      legend.position = "right",
      panel.spacing = unit(0.6, "lines")
    )
}

make_physio_epoch_combined_plot <- function(data_bundle) {
  shared_xlim <- c(-50, 128)
  
  p_endo_epoch <- make_physio_stacked_plot(
    data_bundle$endo_all_epoch_sum,
    xlims = shared_xlim,
    base_size = 12,
    point_size = 2.5,
    line_width = 1,
    error_width = 2,
    x_text_size = 12
  ) +
    ggtitle("A") +
    theme(plot.title = element_text(face = "bold", hjust = 0))
  
  p_auto_epoch <- make_physio_stacked_plot(
    data_bundle$auto_epoch_sum,
    xlims = shared_xlim,
    base_size = 12,
    point_size = 2.5,
    line_width = 1,
    error_width = 2,
    x_text_size = 12
  ) +
    ggtitle("B") +
    theme(plot.title = element_text(face = "bold", hjust = 0)) +
    guides(x = guide_axis(n.dodge = 2))
  
  (p_endo_epoch / p_auto_epoch) + plot_layout(heights = c(1.4, 1))
}

save_physio_epoch_timeline_plot <- function(output_dir = physio_plot_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  data_bundle <- build_physio_epoch_timeline_data()
  if (is.null(data_bundle)) {
    return(NULL)
  }
  plot_object <- make_physio_epoch_combined_plot(data_bundle)
  
  save_plot_all_formats(
    plot_object,
    file_stem = "physio_epoch_stacked_plot",
    output_dir = output_dir,
    width = 10,
    height = 9
  )
  
  plot_object
}


# -------------------------------------------------------------------------
# Save all publication figures
# -------------------------------------------------------------------------

save_publication_figures <- function() {
  list(
    eccentricity = list(
      parcel_plots = save_parcel_plots(),
      six_rs_eccentricity_plot = save_six_rs_eccentricity_plot()
    ),
    dispersion = list(
      between_network_six_rs_plot = save_between_network_six_rs_plot(),
      ridge_plots = save_network_ridge_plots()
    ),
    peripheral_measures = list(
      physio_best_outcome_grid = save_physio_best_outcome_grid(),
      physio_single_best_plots = save_physio_single_best_plots(),
      physio_sn_dmn_cortisol_plot = save_physio_selected_pair_plot("delta_distance_default_salience", "SaCort_delta"),
      physio_stress_only_plots = save_physio_stress_only_spearman_plots(),
      physio_epoch_timeline_plot = save_physio_epoch_timeline_plot()
    )
  )
}

generated_plots <- save_publication_figures()
print(generated_plots$eccentricity$parcel_plots[["parcel_100"]][["g1"]])