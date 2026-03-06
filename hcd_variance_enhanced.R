#!/usr/bin/env Rscript
# ============================================================
# HCD Variance Components Analysis - ENHANCED VERSION
# ============================================================
# Enhancements over previous version:
#   - Outlier/jackpot mutation detection
#   - Confidence intervals for variance components
#   - CV% calculation for between-study variation
#   - Optional random slopes model
#   - Enhanced residual diagnostics
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(tidyr)
})

# ============================================================
# CONFIGURATION
# ============================================================

CONFIG <- list(
  # Analysis settings
  routine_group_size = 5,
  exclude_study1 = TRUE,
  
  # Outlier detection
  outlier_threshold_sd = 3.0,  # Z-score threshold for jackpot mutations
  
  # Model options
  fit_random_slopes = FALSE,  # Set TRUE to allow sampling day effects to vary by study
  
  # Output settings
  plot_width = 10,
  plot_height = 6,
  plot_dpi = 300
)

# ============================================================
# DATA CLEANING FUNCTIONS (Same as before)
# ============================================================

detect_header_row <- function(raw_data) {
  for (i in 1:min(10, nrow(raw_data))) {
    row_text <- paste(tolower(raw_data[i, ]), collapse = " ")
    if (grepl("study", row_text) && 
        (grepl("animal", row_text) || grepl("sampling", row_text))) {
      message("Detected header row: ", i)
      return(i)
    }
  }
  return(1)
}

fill_down_column <- function(data, col_name) {
  if (!(col_name %in% names(data))) return(data)
  
  current_value <- NA
  for (i in seq_len(nrow(data))) {
    val <- data[[col_name]][i]
    if (!is.na(val) && val != "" && val != " ") {
      current_value <- val
    } else {
      data[[col_name]][i] <- current_value
    }
  }
  data
}

extract_study_number <- function(study_col) {
  study_no <- rep(0, length(study_col))
  current_code <- ""
  current_no <- 0
  
  for (i in seq_along(study_col)) {
    val <- as.character(study_col[i])
    
    # Check if study code (starts with letter)
    if (grepl("^[A-Za-z]", val)) {
      current_code <- val
      # Look ahead for number
      if (i < length(study_col)) {
        next_val <- as.character(study_col[i + 1])
        if (!is.na(next_val) && grepl("^[0-9]+$", next_val)) {
          current_no <- as.numeric(next_val)
        }
      }
    } else if (grepl("^[0-9]+$", val)) {
      current_no <- as.numeric(val)
    }
    
    study_no[i] <- current_no
  }
  study_no
}

read_messy_hcd_data <- function(csv_path) {
  message("Reading HCD data: ", csv_path)
  
  # Try UTF-8 first, fallback to native encoding
  raw_data <- tryCatch({
    read.csv(csv_path, stringsAsFactors = FALSE, 
             check.names = FALSE, header = FALSE,
             fileEncoding = "UTF-8")
  }, error = function(e) {
    message("UTF-8 failed, trying native encoding...")
    read.csv(csv_path, stringsAsFactors = FALSE, 
             check.names = FALSE, header = FALSE,
             fileEncoding = "")
  })
  
  header_row <- detect_header_row(raw_data)
  
  if (header_row > 1) {
    data <- raw_data[-(1:(header_row - 1)), ]
    names(data) <- as.character(data[1, ])
    data <- data[-1, ]
  } else {
    names(data) <- as.character(raw_data[1, ])
    data <- raw_data[-1, ]
  }
  
  rownames(data) <- NULL
  
  # Identify columns
  study_col <- names(data)[grep("study", tolower(names(data)))[1]]
  sampling_col <- names(data)[grep("sampling", tolower(names(data)))[1]]
  animal_col <- names(data)[grep("animal", tolower(names(data)))[1]]
  
  mutfreq_patterns <- c("mut.*freq.*10", "mutfreq.*x.*10", "mut.*freq", "frequency")
  mutfreq_col <- NULL
  for (pattern in mutfreq_patterns) {
    matches <- grep(pattern, tolower(names(data)))
    if (length(matches) > 0) {
      mutfreq_col <- names(data)[matches[1]]
      break
    }
  }
  
  if (is.null(mutfreq_col)) {
    stop("Could not find MUT_FREQ column")
  }
  
  message("Key columns:")
  message("  Study: ", study_col)
  message("  Sampling: ", sampling_col)
  message("  MUT_FREQ: ", mutfreq_col)
  
  # Fill down and extract
  data <- fill_down_column(data, study_col)
  data <- fill_down_column(data, sampling_col)
  data$Study_No <- extract_study_number(data[[study_col]])
  data$Study_Code <- data[[study_col]]
  data$Sampling_Day <- as.numeric(data[[sampling_col]])
  data$Animal_ID <- data[[animal_col]]
  data$MUT_FREQ <- as.numeric(data[[mutfreq_col]])
  
  data_clean <- data %>%
    select(Study_No, Study_Code, Sampling_Day, Animal_ID, MUT_FREQ) %>%
    filter(!is.na(Study_No), Study_No > 0) %>%
    filter(!is.na(MUT_FREQ), is.finite(MUT_FREQ), MUT_FREQ > 0) %>%
    filter(!is.na(Sampling_Day)) %>%
    mutate(Log_MF = log10(MUT_FREQ)) %>%
    arrange(Study_No, Sampling_Day)
  
  if (CONFIG$exclude_study1) {
    message("Excluding Study 1 (large initial group)")
    data_clean <- data_clean %>% filter(Study_No != 1)
  }
  
  message("Data ready: ", nrow(data_clean), " animals from ", 
          n_distinct(data_clean$Study_No), " studies")
  
  data_clean
}

# ============================================================
# NEW: OUTLIER DETECTION
# ============================================================

detect_outliers <- function(data) {
  message("\n=== OUTLIER DETECTION ===")
  message("Checking for jackpot mutations (>", CONFIG$outlier_threshold_sd, " SD from study mean)...")
  
  data_flagged <- data %>%
    group_by(Study_No, Sampling_Day) %>%
    mutate(
      study_mean = mean(Log_MF, na.rm = TRUE),
      study_sd = sd(Log_MF, na.rm = TRUE),
      z_score = (Log_MF - study_mean) / study_sd,
      outlier_flag = abs(z_score) > CONFIG$outlier_threshold_sd
    ) %>%
    ungroup()
  
  outliers <- data_flagged %>% filter(outlier_flag)
  
  if (nrow(outliers) > 0) {
    message("WARNING: ", nrow(outliers), " potential jackpot mutations detected:")
    for (i in 1:min(5, nrow(outliers))) {
      message(sprintf("  Study %d, Day %d, Animal %s: MF=%.1f, Z=%.2f",
                      outliers$Study_No[i], outliers$Sampling_Day[i], 
                      outliers$Animal_ID[i], outliers$MUT_FREQ[i], outliers$z_score[i]))
    }
    if (nrow(outliers) > 5) message("  ... and ", nrow(outliers) - 5, " more")
    
    message("\nRecommendation: Review these animals. Consider removing if confirmed jackpot.")
    message("To exclude, set their MUT_FREQ to NA in source data.\n")
  } else {
    message("No extreme outliers detected.\n")
  }
  
  data_flagged
}

# ============================================================
# VARIANCE COMPONENTS ANALYSIS - ENHANCED
# ============================================================

fit_variance_model <- function(data) {
  message("\n=== VARIANCE COMPONENTS MODEL ===")
  
  # Prepare data
  data$Study_No <- factor(data$Study_No)
  data$Sampling_Day <- factor(data$Sampling_Day)
  
  # Fit base model (random intercept only)
  message("Fitting REML model: Log_MF ~ Sampling_Day + (1|Study_No)")
  model_base <- lmer(Log_MF ~ Sampling_Day + (1|Study_No), data = data, REML = TRUE)
  
  # Optional: Fit random slopes model
  model_slopes <- NULL
  if (CONFIG$fit_random_slopes) {
    message("Fitting random slopes model: Log_MF ~ Sampling_Day + (Sampling_Day|Study_No)")
    tryCatch({
      model_slopes <- lmer(Log_MF ~ Sampling_Day + (Sampling_Day|Study_No), 
                          data = data, REML = TRUE)
      
      # Compare models
      anova_result <- anova(model_base, model_slopes)
      if (anova_result$`Pr(>Chisq)`[2] < 0.05) {
        message("Random slopes model significantly better (p < 0.05)")
        model_base <- model_slopes
      } else {
        message("Random intercept model adequate (p >= 0.05)")
      }
    }, error = function(e) {
      message("Random slopes model failed to converge, using random intercept only")
    })
  }
  
  model_base
}

calculate_variance_components <- function(model, data) {
  message("\n=== VARIANCE DECOMPOSITION ===")
  
  # Extract variance components
  vc <- as.data.frame(VarCorr(model))
  sigma_between <- vc$vcov[vc$grp == "Study_No"]
  sigma_residual <- vc$vcov[vc$grp == "Residual"]
  sigma_total <- sigma_between + sigma_residual
  
  # Percentage contributions
  pct_between <- (sigma_between / sigma_total) * 100
  pct_within <- (sigma_residual / sigma_total) * 100
  
  # NEW: Coefficient of Variation (CV%)
  overall_mean <- mean(data$Log_MF, na.rm = TRUE)
  cv_percent <- (sqrt(sigma_between) / overall_mean) * 100
  
  # NEW: Confidence intervals for variance components
  message("Calculating confidence intervals (may take a moment)...")
  ci <- NULL
  tryCatch({
    ci <- confint(model, oldNames = FALSE, quiet = TRUE)
  }, error = function(e) {
    message("Could not calculate CIs: ", e$message)
  })
  
  results <- list(
    sigma_between = sigma_between,
    sigma_within = sigma_residual,
    sigma_total = sigma_total,
    pct_between = pct_between,
    pct_within = pct_within,
    cv_percent = cv_percent,
    confidence_intervals = ci
  )
  
  message(sprintf("Between-study variance: %.4f (%.1f%%)", sigma_between, pct_between))
  message(sprintf("Within-study variance:  %.4f (%.1f%%)", sigma_residual, pct_within))
  message(sprintf("Between-study CV%%:      %.1f%%", cv_percent))
  
  if (!is.null(ci)) {
    between_ci <- ci[grep("Study_No", rownames(ci)), ]
    if (length(between_ci) >= 2) {
      message(sprintf("95%% CI for between-study SD: [%.4f, %.4f]", 
                      sqrt(between_ci[1]), sqrt(between_ci[2])))
    }
  }
  
  results
}

test_sampling_day_effect <- function(model, data) {
  message("\n=== SAMPLING DAY EFFECT ===")
  
  # ANOVA F-test
  anova_result <- anova(model)
  p_value <- anova_result$`Pr(>F)`[1]
  
  # Effect size (difference between days)
  day_means <- data %>%
    group_by(Sampling_Day) %>%
    summarise(mean_log = mean(Log_MF, na.rm = TRUE), .groups = "drop")
  
  if (nrow(day_means) == 2) {
    diff_log <- abs(diff(day_means$mean_log))
    diff_orig <- 10^day_means$mean_log[2] - 10^day_means$mean_log[1]
    
    message(sprintf("Day 31 mean (log10): %.3f", day_means$mean_log[day_means$Sampling_Day == 31]))
    message(sprintf("Day 56 mean (log10): %.3f", day_means$mean_log[day_means$Sampling_Day == 56]))
    message(sprintf("Difference (log10):  %.3f", diff_log))
    message(sprintf("p-value:             %.4f %s", p_value, 
                    ifelse(p_value < 0.001, "***", 
                           ifelse(p_value < 0.01, "**", 
                                  ifelse(p_value < 0.05, "*", "ns")))))
  }
  
  list(
    p_value = p_value,
    anova_table = anova_result,
    day_means = day_means
  )
}

calculate_prediction_intervals <- function(model, data, variance_components) {
  message("\n=== PREDICTION INTERVALS ===")
  
  # Overall mean
  overall_mean_log <- mean(data$Log_MF, na.rm = TRUE)
  overall_mean_original <- 10^overall_mean_log
  
  # Prediction interval incorporating both sources of variation
  sigma_total <- sqrt(variance_components$sigma_total)
  
  # 95% prediction interval
  lower_log <- overall_mean_log - 1.96 * sigma_total
  upper_log <- overall_mean_log + 1.96 * sigma_total
  
  lower_original <- 10^lower_log
  upper_original <- 10^upper_log
  
  message(sprintf("Overall geometric mean: %.1f x 10^-6", overall_mean_original))
  message(sprintf("95%% prediction interval (log10): [%.3f, %.3f]", lower_log, upper_log))
  message(sprintf("95%% prediction interval (original): [%.1f, %.1f] x 10^-6", 
                  lower_original, upper_original))
  
  list(
    mean_log = overall_mean_log,
    mean_original = overall_mean_original,
    lower_log = lower_log,
    upper_log = upper_log,
    lower_original = lower_original,
    upper_original = upper_original
  )
}

# ============================================================
# TABLE 3: VARIANCE BY SAMPLING DAY (SEPARATE ANALYSES)
# ============================================================

calculate_variance_by_day <- function(data) {
  message("\n=== TABLE 3: VARIANCE COMPONENTS BY SAMPLING DAY ===")
  
  results_list <- list()
  
  # Analyze each sampling day separately
  for (day in unique(data$Sampling_Day)) {
    message(sprintf("\nAnalyzing Sampling Day %d...", day))
    
    day_data <- data[data$Sampling_Day == day, ]
    
    # Check if we have enough studies for this day
    n_studies <- length(unique(day_data$Study_No))
    
    if (n_studies < 2) {
      message(sprintf("  Skipping Day %d: only %d study (need at least 2 for variance analysis)", 
                      day, n_studies))
      next
    }
    
    # Simple model: just between-study and within-study for this day
    tryCatch({
      model_day <- lmer(Log_MF ~ (1|Study_No), data = day_data)
      
      # Extract variance components
      vc_day <- as.data.frame(VarCorr(model_day))
      
      sigma_between_day <- vc_day$vcov[1]  # Between study
      sigma_within_day <- vc_day$vcov[2]   # Residual (within study)
      sigma_total_day <- sigma_between_day + sigma_within_day
      
      pct_between_day <- (sigma_between_day / sigma_total_day) * 100
      pct_within_day <- (sigma_within_day / sigma_total_day) * 100
      
      message(sprintf("  Between-study: %.4f (%.1f%%)", sigma_between_day, pct_between_day))
      message(sprintf("  Within-study:  %.4f (%.1f%%)", sigma_within_day, pct_within_day))
      message(sprintf("  Total:         %.4f", sigma_total_day))
      
      results_list[[paste0("Day_", day)]] <- data.frame(
        Sampling_Day = day,
        Source = c("Between Study", "Within Study", "Total"),
        Variance_Component = c(sigma_between_day, sigma_within_day, sigma_total_day),
        Percent_of_Total = c(pct_between_day, pct_within_day, 100),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      message(sprintf("  Error analyzing Day %d: %s", day, e$message))
    })
  }
  
  # Combine all days
  if (length(results_list) > 0) {
    do.call(rbind, results_list)
  } else {
    data.frame(
      Sampling_Day = numeric(0),
      Source = character(0),
      Variance_Component = numeric(0),
      Percent_of_Total = numeric(0)
    )
  }
}

# ============================================================
# VISUALIZATION - STATISTICIAN STYLE
# ============================================================

create_statistician_boxplot <- function(data, outfile) {
  # Create plot matching statistician's Figure 1 format
  # Box plot by Study Number and Sampling Day
  
  p <- ggplot(data, aes(x = factor(Study_No), y = Log_MF, fill = factor(Sampling_Day))) +
    geom_boxplot(alpha = 0.7, outlier.shape = 1) +
    scale_fill_manual(values = c("31" = "#4CAF50", "36" = "#FFEB3B", "56" = "#2196F3"),
                     labels = c("Day 31", "Day 36", "Day 56"),
                     name = "Sampling Day") +
    labs(
      title = "Positive Control Liver Mutation Frequency by Study Number and Sampling Day",
      subtitle = paste0("n = ", nrow(data), " animals from ", n_distinct(data$Study_No), " studies"),
      x = "Study Number",
      y = "Log10 Mutant Frequency (x10^-6)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom",
      panel.grid.major.x = element_blank()
    ) +
    # Add reference line at overall mean
    geom_hline(yintercept = mean(data$Log_MF, na.rm = TRUE), 
               linetype = "dashed", color = "red", alpha = 0.5)
  
  ggsave(outfile, p, width = 12, height = 6, dpi = 300)
  message("Statistician-style boxplot saved: ", outfile)
}

create_enhanced_boxplot <- function(data, outfile) {
  # Simpler plot by sampling day only
  p <- ggplot(data, aes(x = factor(Sampling_Day), y = Log_MF)) +
    geom_boxplot(aes(fill = factor(Sampling_Day)), alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_manual(values = c("31" = "#4CAF50", "36" = "#FFEB3B", "56" = "#2196F3"),
                     name = "Sampling Day") +
    labs(
      title = "Mutation Frequency by Sampling Day (Summary)",
      x = "Sampling Day",
      y = "Log10 Mutant Frequency",
      fill = "Day"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  
  ggsave(outfile, p, width = 8, height = 6, dpi = 300)
  message("Enhanced boxplot saved: ", outfile)
}

# ============================================================
# MAIN PIPELINE
# ============================================================

run_variance_analysis <- function(in_csv, out_summary, out_plots) {
  
  message("============================================================")
  message("HCD VARIANCE COMPONENTS ANALYSIS - ENHANCED")
  message("============================================================\n")
  
  # Read data
  data <- read_messy_hcd_data(in_csv)
  
  # NEW: Detect outliers
  data <- detect_outliers(data)
  
  # Fit model
  model <- fit_variance_model(data)
  
  # Calculate variance components (with CIs and CV%)
  vc <- calculate_variance_components(model, data)
  
  # Test sampling day effect
  day_test <- test_sampling_day_effect(model, data)
  
  # Prediction intervals
  pred_int <- calculate_prediction_intervals(model, data, vc)
  
  # TABLE 3: Variance by sampling day (separate analyses)
  variance_by_day <- calculate_variance_by_day(data)
  
  # Create STATISTICIAN-STYLE summary (Table 2 format)
  stat_summary <- data.frame(
    Source = c("Sampling Day", "Between Study (within day)", "Within Study", "Total"),
    Variance_Component = c(
      NA,  # Sampling day variance not directly in this model structure
      vc$sigma_between,
      vc$sigma_within,
      vc$sigma_total
    ),
    Percent_of_Total = c(
      NA,
      vc$pct_between,
      vc$pct_within,
      100
    ),
    stringsAsFactors = FALSE
  )
  
  # Add F-test results
  ftest_table <- data.frame(
    Test = "Sampling Day Effect",
    F_statistic = sprintf("%.2f", day_test$anova_table$`F value`[1]),
    df_num = day_test$anova_table$NumDF[1],
    df_denom = day_test$anova_table$DenDF[1],
    p_value = sprintf("%.4f", day_test$p_value),
    Significant = ifelse(day_test$p_value < 0.001, "***",
                        ifelse(day_test$p_value < 0.01, "**",
                              ifelse(day_test$p_value < 0.05, "*", "ns"))),
    stringsAsFactors = FALSE
  )
  
  # Add prediction interval
  pred_table <- data.frame(
    Metric = c("Geometric Mean (x10^-6)", 
               "Geometric SD",
               "95% Prediction Interval Lower (x10^-6)",
               "95% Prediction Interval Upper (x10^-6)",
               "Between Study %",
               "Total SD (log scale)"),
    Value = c(
      sprintf("%.1f", pred_int$mean_original),
      sprintf("%.2f", 10^sqrt(vc$sigma_total)),
      sprintf("%.1f", pred_int$lower_original),
      sprintf("%.1f", pred_int$upper_original),
      sprintf("%.1f%%", vc$pct_between),
      sprintf("%.3f", sqrt(vc$sigma_total))
    ),
    stringsAsFactors = FALSE
  )
  
  # Create ENHANCED summary table (original format)
  enhanced_summary <- data.frame(
    Metric = c(
      "Between_Study_Variance",
      "Within_Study_Variance",
      "Total_Variance",
      "Between_Study_Percent",
      "Within_Study_Percent",
      "Between_Study_CV_Percent",
      "Sampling_Day_p_value",
      "Geometric_Mean_x10_6",
      "Prediction_Lower_x10_6",
      "Prediction_Upper_x10_6"
    ),
    Value = c(
      vc$sigma_between,
      vc$sigma_within,
      vc$sigma_total,
      vc$pct_between,
      vc$pct_within,
      vc$cv_percent,
      day_test$p_value,
      pred_int$mean_original,
      pred_int$lower_original,
      pred_int$upper_original
    ),
    stringsAsFactors = FALSE
  )
  
  # Save each table separately to CSV
  # Main summary file will have all sections
  sink(out_summary)
  
  cat("TABLE 2: Components of Variance for Liver (Combined Analysis)\n")
  cat("==============================================================\n")
  write.table(stat_summary, sep = ",", row.names = FALSE, quote = FALSE)
  cat("\n")
  
  cat("F-Test for Sampling Day Effect\n")
  cat("================================\n")
  write.table(ftest_table, sep = ",", row.names = FALSE, quote = FALSE)
  cat("\n")
  
  cat("Prediction Intervals & Statistics\n")
  cat("==================================\n")
  write.table(pred_table, sep = ",", row.names = FALSE, quote = FALSE)
  cat("\n")
  
  cat("TABLE 3: Components of Variance for Liver by Sampling Day (Separate Analyses)\n")
  cat("==============================================================================\n")
  write.table(variance_by_day, sep = ",", row.names = FALSE, quote = FALSE)
  cat("\n")
  
  cat("Enhanced Metrics (with CV% and CIs)\n")
  cat("====================================\n")
  write.table(enhanced_summary, sep = ",", row.names = FALSE, quote = FALSE)
  
  sink()
  message("\nSummary saved: ", out_summary)
  
  # Create BOTH plot styles
  # Statistician style (Figure 1 format)
  stat_plot <- sub("\\.png$", "_figure1.png", out_plots)
  create_statistician_boxplot(data, stat_plot)
  
  # Enhanced style (summary by day)
  enh_plot <- sub("\\.png$", "_summary.png", out_plots)
  create_enhanced_boxplot(data, enh_plot)
  
  message("\n============================================================")
  message("ANALYSIS COMPLETE")
  message("============================================================")
  
  invisible(list(
    model = model,
    variance_components = vc,
    day_test = day_test,
    prediction_intervals = pred_int,
    data = data,
    stat_plot = stat_plot,
    enh_plot = enh_plot
  ))
}

# ============================================================
# COMMAND LINE INTERFACE
# ============================================================

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 3) {
    stop("Usage: Rscript hcd_variance_enhanced.R <input.csv> <summary.csv> <plots.png>")
  }
  
  tryCatch({
    run_variance_analysis(args[1], args[2], args[3])
  }, error = function(e) {
    message("ERROR: ", e$message)
    quit(status = 1)
  })
}
