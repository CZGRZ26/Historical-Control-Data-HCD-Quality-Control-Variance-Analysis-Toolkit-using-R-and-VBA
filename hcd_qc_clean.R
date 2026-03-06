#!/usr/bin/env Rscript
# ============================================================
# HCD QC Control Charts - For Pre-Cleaned Data
# ============================================================
# Expects clean CSV from VBA with columns:
#   Study / Study No, Study_No, Sampling day, Animal No, Mut Freq x 10-6, Log10_MF
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

CONFIG <- list(
  exclude_study1_from_limits = TRUE,
  exclude_study1_from_plot   = FALSE,
  d2_range_size_2 = 1.128,
  plot_width  = 9,
  plot_height = 4.5,
  plot_dpi    = 300
)

# ============================================================
# READ CLEAN DATA
# ============================================================

read_clean_data <- function(csv_path) {
  message("Reading clean QC data from: ", csv_path)
  
  # Read CSV with standard encoding
  data <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  message("Columns found: ", paste(names(data), collapse = ", "))
  
  # Identify columns (flexible matching)
  # Clean column names first - remove line breaks
  clean_names <- gsub("\n", " ", names(data))
  clean_names <- gsub("\r", " ", clean_names)
  clean_names <- gsub("\\s+", " ", clean_names)
  
  # Look for Study_No column (exact match preferred)
  study_no_col <- names(data)[which(clean_names == "Study_No")]
  if (length(study_no_col) == 0) {
    # Fallback: look for column with "study" and "no" but NOT "/"
    for (i in seq_along(clean_names)) {
      cn <- tolower(clean_names[i])
      if (grepl("study", cn) && grepl("no", cn) && !grepl("/", cn)) {
        study_no_col <- names(data)[i]
        break
      }
    }
  }
  
  animal_col <- names(data)[grep("animal", tolower(clean_names))]
  
  mutfreq_col <- names(data)[grep("mut.*freq.*10|mutfreq.*x.*10", tolower(clean_names))]
  
  log_col <- names(data)[grep("log10", tolower(clean_names))]
  
  # Validate
  if (length(study_no_col) == 0) stop("Could not find Study_No column")
  if (length(mutfreq_col) == 0) stop("Could not find Mut Freq x 10-6 column")
  
  message("Key columns identified:")
  message("  Study_No: ", study_no_col[1])
  message("  Animal: ", animal_col[1])
  message("  MUT_FREQ: ", mutfreq_col[1])
  if (length(log_col) > 0) message("  Log10: ", log_col[1])
  
  # Create clean data frame
  data_clean <- data.frame(
    Study_No = as.numeric(data[[study_no_col[1]]]),
    Animal_ID = data[[animal_col[1]]],
    MUT_FREQ = as.numeric(data[[mutfreq_col[1]]]),
    stringsAsFactors = FALSE
  )
  
  # Use pre-calculated log10 if available, otherwise calculate
  if (length(log_col) > 0) {
    data_clean$Log_MF <- as.numeric(data[[log_col[1]]])
  } else {
    data_clean$Log_MF <- log10(data_clean$MUT_FREQ)
  }
  
  # Remove invalid rows
  data_clean <- data_clean %>%
    filter(!is.na(Study_No), Study_No > 0) %>%
    filter(!is.na(MUT_FREQ), is.finite(MUT_FREQ)) %>%
    filter(!is.na(Log_MF), is.finite(Log_MF)) %>%
    arrange(Study_No)
  
  message("Data ready: ", nrow(data_clean), " animals from ", n_distinct(data_clean$Study_No), " studies")
  
  data_clean
}

# ============================================================
# QC CONTROL CHART FUNCTIONS
# ============================================================

calculate_control_limits <- function(data) {
  # Calculate study means
  study_means <- data %>%
    group_by(Study_No) %>%
    summarise(mean_log_mf = mean(Log_MF, na.rm = TRUE), .groups = "drop") %>%
    arrange(Study_No)
  
  x_raw <- study_means$mean_log_mf
  mean_x <- mean(x_raw, na.rm = TRUE)
  
  # Moving Range
  mr <- abs(diff(x_raw))
  mr_bar <- mean(mr, na.rm = TRUE)
  sigma <- mr_bar / CONFIG$d2_range_size_2
  
  # Stability index
  sd_all <- sd(x_raw, na.rm = TRUE)
  stability_index <- sd_all / sigma
  
  list(
    mean = mean_x,
    mr_bar = mr_bar,
    sigma = sigma,
    sd_all = sd_all,
    stability_idx = stability_index,
    lcl_2sd = mean_x - 2 * sigma,
    ucl_2sd = mean_x + 2 * sigma,
    lcl_3sd = mean_x - 3 * sigma,
    ucl_3sd = mean_x + 3 * sigma
  )
}

normalize_to_overall_mean <- function(data, overall_mean) {
  data %>%
    group_by(Study_No) %>%
    mutate(study_mean = mean(Log_MF, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(Log_MF_norm = Log_MF - study_mean + overall_mean) %>%
    select(-study_mean)
}

flag_observations <- function(data, limits) {
  data %>%
    mutate(
      Flag = case_when(
        Log_MF_norm < limits$lcl_3sd | Log_MF_norm > limits$ucl_3sd ~ "Action (outside 3σ)",
        Log_MF_norm < limits$lcl_2sd | Log_MF_norm > limits$ucl_2sd ~ "Warning (outside 2σ)",
        TRUE ~ "In control"
      )
    )
}

create_control_chart <- function(data, limits, sampling_day) {
  subtitle <- sprintf(
    "Normalised points; limits from study means (Study 1 excluded). Stability = %.2f",
    limits$stability_idx
  )
  
  ggplot(data, aes(x = Study_No, y = Log_MF_norm)) +
    geom_point(aes(shape = Flag), size = 2) +
    geom_hline(aes(yintercept = limits$mean, colour = "Mean"), linewidth = 0.6) +
    geom_hline(aes(yintercept = limits$ucl_2sd, colour = "±2 SD"), linetype = "dashed") +
    geom_hline(aes(yintercept = limits$lcl_2sd, colour = "±2 SD"), linetype = "dashed") +
    geom_hline(aes(yintercept = limits$ucl_3sd, colour = "±3 SD"), linetype = "dotted") +
    geom_hline(aes(yintercept = limits$lcl_3sd, colour = "±3 SD"), linetype = "dotted") +
    scale_colour_manual(
      name = "Control limits",
      values = c("Mean" = "grey40", "±2 SD" = "darkgoldenrod1", "±3 SD" = "red")
    ) +
    scale_x_continuous(breaks = seq(min(data$Study_No), max(data$Study_No), by = 2)) +
    labs(
      title = sprintf("Control Chart for Individual +ve Control Liver MF: Sampling Day = %d", sampling_day),
      subtitle = subtitle,
      x = "Study",
      y = "Log10 MF (normalised)"
    ) +
    theme_classic()
}

# ============================================================
# MAIN QC PIPELINE
# ============================================================

run_qc_analysis <- function(in_csv, out_sum, out_data, out_png) {
  
  message("=" %>% rep(60) %>% paste(collapse = ""))
  message("HCD QC CONTROL CHART ANALYSIS")
  message("=" %>% rep(60) %>% paste(collapse = ""))
  message()
  
  # Read clean data
  data_all <- read_clean_data(in_csv)
  
  # Detect sampling day from filename or use default
  sampling_day <- 31
  if (grepl("56|day56", tolower(in_csv))) {
    sampling_day <- 56
  }
  
  # Apply exclusions for limits
  data_limits <- data_all
  if (CONFIG$exclude_study1_from_limits) {
    data_limits <- data_limits %>% filter(Study_No != 1)
    message("Excluding Study 1 from control limit calculations")
  }
  
  # Calculate control limits
  message("Calculating control limits...")
  limits <- calculate_control_limits(data_limits)
  
  # Apply exclusions for plot
  data_plot <- data_all
  if (CONFIG$exclude_study1_from_plot) {
    data_plot <- data_plot %>% filter(Study_No != 1)
  }
  
  # Normalize and flag
  message("Normalizing individual values...")
  data_plot <- data_plot %>%
    normalize_to_overall_mean(limits$mean) %>%
    flag_observations(limits)
  
  # Create outputs
  message("Generating outputs...")
  summary_tbl <- data.frame(
    Metric = c("Mean_log10", "MR_bar", "Sigma_MR_d2", "SD_all", "Stability_index",
               "LCL_2SD", "UCL_2SD", "LCL_3SD", "UCL_3SD"),
    Value = c(limits$mean, limits$mr_bar, limits$sigma, limits$sd_all, limits$stability_idx,
              limits$lcl_2sd, limits$ucl_2sd, limits$lcl_3sd, limits$ucl_3sd)
  )
  
  plot <- create_control_chart(data_plot, limits, sampling_day)
  
  # Write files
  write.csv(summary_tbl, out_sum, row.names = FALSE)
  write.csv(data_plot, out_data, row.names = FALSE)
  ggsave(out_png, plot, width = CONFIG$plot_width, height = CONFIG$plot_height, dpi = CONFIG$plot_dpi)
  
  message()
  message("QC analysis complete!")
  message("  Summary: ", out_sum)
  message("  Data:    ", out_data)
  message("  Plot:    ", out_png)
  
  invisible(list(summary = summary_tbl, data = data_plot, limits = limits))
}

# ============================================================
# COMMAND LINE INTERFACE
# ============================================================

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 4) {
    stop("Usage: Rscript hcd_qc_clean.R <input.csv> <summary.csv> <data.csv> <plot.png>")
  }
  
  tryCatch({
    run_qc_analysis(args[1], args[2], args[3], args[4])
  }, error = function(e) {
    message("ERROR: ", e$message)
    quit(status = 1)
  })
}