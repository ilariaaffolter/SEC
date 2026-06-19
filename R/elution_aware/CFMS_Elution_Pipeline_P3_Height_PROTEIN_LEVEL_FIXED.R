# ═══════════════════════════════════════════════════════════════════════════════
# PROTEIN-LEVEL SEC-MS PIPELINE (Simplified - No Proteoform Clustering)
# ═══════════════════════════════════════════════════════════════════════════════
#
# This is a SIMPLIFIED version that:
# 1. Peptide-to-protein aggregation per fraction per sample
#    - maxlfq: MaxLFQ across samples (most accurate)
#    - median: Median of all peptides (robust)
#    - sum: Sum of all peptides
#    - top3, top5, top10: Sum of top N peptides by intensity
#    - top10median: Median of top 10 peptides (min 4 required, robust)
# 2. NO proteoform clustering - works at protein level only
# 3. Peak detection on protein chromatograms from CONTROL
# 4. Adaptive (FWHM-based) or fixed peak boundaries for quantification
# 5. Quantify all samples at reference peak positions
# 6. Limma differential analysis with IBMT
# 7. Flexible p-value adjustment for correlated CF-MS data:
#    - BH: Benjamini-Hochberg (standard, assumes independence)
#    - qvalue: Storey's q-value (estimates pi0, less conservative)
#    - none: Raw p-values (exploratory)
#    - Permutation FDR: Empirical FDR (gold standard for correlations)
# 8. Shows individual peptide traces in plots
#
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(limma)
  library(ggplot2)
  library(ggrepel)
})

# ── Resolve namespace conflicts ──────────────────────────────────────────────
# Bioconductor packages (loaded by limma or its dependencies, e.g.
# AnnotationDbi, GenomicRanges) mask dplyr::select and dplyr::filter with
# S4 generics. Force dplyr's versions back on top.
select <- dplyr::select
filter <- dplyr::filter

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1: HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

`%||%` <- function(a, b) if (is.null(a)) b else a

is_finite_num <- function(x) !is.na(x) & is.finite(x)

clamp <- function(x, lo, hi) pmax(lo, pmin(hi, x))

#' Create volcano plot with consistent styling
#' @param data Data frame with limma_logFC and p-value column
#' @param pval_col Name of p-value column to use
#' @param pval_label Label for y-axis (e.g., "adj. P-value", "q-value", "perm FDR")
#' @param title_suffix Suffix for title
#' @param treat_condition Treatment condition name
#' @param ref_condition Reference condition name
#' @param agg_label Aggregation method label
#' @return ggplot object
create_volcano_plot <- function(data, pval_col, pval_label, title_suffix,
                                 treat_condition, ref_condition, agg_label) {
  
  volcano_data <- data %>%
    filter(!is.na(.data[[pval_col]]) & is.finite(.data[[pval_col]]) & .data[[pval_col]] > 0) %>%
    mutate(
      neg_log_p = -log10(.data[[pval_col]]),
      significant = .data[[pval_col]] < 0.05 & abs(limma_logFC) > 1,
      high_fc = abs(limma_logFC) > 2
    )
  
  if (nrow(volcano_data) == 0) return(NULL)
  
  n_significant <- sum(volcano_data$significant, na.rm = TRUE)
  
  if (n_significant > 0) {
    volcano_data <- volcano_data %>%
      mutate(
        should_label = significant | high_fc,
        label_priority = case_when(
          significant & high_fc ~ 1,
          significant ~ 2,
          high_fc ~ 3,
          TRUE ~ 4
        )
      ) %>%
      arrange(label_priority, desc(abs(limma_logFC)))
    
    labels_to_show <- volcano_data %>%
      filter(should_label) %>%
      head(50) %>%
      pull(feature_id)
    
    volcano_data$label <- ifelse(volcano_data$feature_id %in% labels_to_show,
                                  volcano_data$protein_id, "")
  } else {
    volcano_data <- volcano_data %>%
      arrange(desc(abs(limma_logFC)))
    
    labels_to_show <- volcano_data %>%
      filter(high_fc) %>%
      head(50) %>%
      pull(feature_id)
    
    volcano_data$label <- ifelse(volcano_data$feature_id %in% labels_to_show,
                                  volcano_data$protein_id, "")
  }
  
  volcano_data$point_color <- case_when(
    volcano_data$significant & volcano_data$high_fc ~ "Significant & |FC|>2",
    volcano_data$significant ~ "Significant",
    volcano_data$high_fc ~ "|FC|>2 (not sig)",
    TRUE ~ "Not significant"
  )
  
  p <- ggplot(volcano_data, aes(x = limma_logFC, y = neg_log_p)) +
    geom_point(aes(color = point_color), alpha = 0.6, size = 2) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = "orange") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_text_repel(
      data = subset(volcano_data, label != ""),
      aes(label = label),
      size = 3, max.overlaps = 50,
      segment.color = "gray50", segment.size = 0.3,
      box.padding = 0.3, point.padding = 0.2
    ) +
    scale_color_manual(values = c(
      "Significant & |FC|>2" = "darkred",
      "Significant" = "red",
      "|FC|>2 (not sig)" = "orange",
      "Not significant" = "gray50"
    )) +
    labs(x = "log2 Fold Change (Limma)", 
         y = paste0("-log10(", pval_label, ")"),
         title = paste0("Volcano Plot (", title_suffix, "): ", treat_condition, " vs ", ref_condition),
         subtitle = paste0("Labeled: significant (<0.05) + |log2FC|>2 (up to 50) | Aggregation: ", agg_label),
         color = "Status") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 11)
    )
  
  attr(p, "n_significant") <- n_significant
  attr(p, "n_high_fc") <- sum(volcano_data$high_fc, na.rm = TRUE)
  attr(p, "n_labels") <- sum(volcano_data$label != "")
  
  return(p)
}

#' Parse sample name to extract condition and replicate
parse_sample_name <- function(sample_name) {
  parts <- strsplit(sample_name, "_")[[1]]
  if (length(parts) >= 2) {
    replicate <- parts[length(parts)]
    condition <- paste(parts[-length(parts)], collapse = "_")
    return(list(condition = condition, replicate = replicate))
  } else {
    return(list(condition = sample_name, replicate = "1"))
  }
}

#' Extract trace matrix from trace object
extract_trace_matrix <- function(trace_obj) {
  if ("traces" %in% names(trace_obj)) {
    trace_dt <- trace_obj$traces
  } else if (is.data.frame(trace_obj) || is.data.table(trace_obj)) {
    trace_dt <- trace_obj
  } else {
    stop("Unexpected trace object structure")
  }
  
  fraction_cols <- grep("^[0-9]+$", colnames(trace_dt), value = TRUE)
  fraction_cols <- fraction_cols[order(as.numeric(fraction_cols))]
  
  if (is.data.table(trace_dt)) {
    trace_mat <- as.matrix(trace_dt[, ..fraction_cols])
    peptide_ids <- if ("id" %in% colnames(trace_dt)) trace_dt$id else trace_dt[[1]]
  } else {
    trace_mat <- as.matrix(trace_dt[, fraction_cols, drop = FALSE])
    peptide_ids <- if ("id" %in% colnames(trace_dt)) trace_dt$id else trace_dt[[1]]
  }
  
  rownames(trace_mat) <- peptide_ids
  colnames(trace_mat) <- fraction_cols
  
  list(matrix = trace_mat, peptide_ids = peptide_ids, fractions = fraction_cols)
}

#' Get protein mapping from traces (checks trace_annotation first)
get_protein_mapping <- function(traces_list, protein_col = NULL) {
  protein_col_names <- c("protein_id", "Protein", "ProteinId", "protein", 
                          "uniprotID", "UniprotID", "uniprot_id", "Uniprot",
                          "ProteinAccession", "Accession", "protein_accession")
  
  if (!is.null(protein_col)) {
    protein_col_names <- c(protein_col, protein_col_names)
  }
  
  first_obj <- traces_list[[1]]
  
  # Check for trace_annotation (CCprofiler/SECAT style)
  if ("trace_annotation" %in% names(first_obj)) {
    annotation_df <- first_obj$trace_annotation
    message("  Found trace_annotation in traces object")
    
    found_col <- NULL
    for (col in protein_col_names) {
      if (col %in% colnames(annotation_df)) {
        found_col <- col
        break
      }
    }
    
    if (is.null(found_col)) {
      message("  WARNING: No protein column found in trace_annotation!")
      return(data.frame(peptide_id = character(), protein_id = character()))
    }
    
    message("  Using protein column: ", found_col)
    id_col <- if ("id" %in% colnames(annotation_df)) "id" else colnames(annotation_df)[1]
    
    mappings <- data.frame(
      peptide_id = annotation_df[[id_col]],
      protein_id = annotation_df[[found_col]],
      stringsAsFactors = FALSE
    ) %>% distinct()
    
    return(mappings)
  }
  
  # Fallback: check traces directly
  message("  WARNING: No trace_annotation found, checking traces directly")
  return(data.frame(peptide_id = character(), protein_id = character()))
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2: MAXLFQ FUNCTION
# ═══════════════════════════════════════════════════════════════════════════════

#' MaxLFQ algorithm - compares peptide intensities ACROSS samples
#' 
#' For each fraction, this takes a matrix of peptide intensities (peptides x samples)
#' and returns one protein abundance per sample, normalized using pairwise ratios.
#' 
#' @param peptide_intensities_matrix Matrix with rows = peptides, cols = samples
#' @return Vector of protein abundances (one per sample)
maxlfq_across_samples <- function(peptide_intensities_matrix) {
  n_samples <- ncol(peptide_intensities_matrix)
  n_peptides <- nrow(peptide_intensities_matrix)
  
  if (n_peptides == 0 || n_samples == 0) return(rep(NA, n_samples))
  if (n_peptides == 1) return(as.numeric(peptide_intensities_matrix[1, ]))
  if (n_samples == 1) return(median(peptide_intensities_matrix[, 1], na.rm = TRUE))
  
  log_intensities <- log2(peptide_intensities_matrix + 1)
  
  # Calculate median pairwise ratios between samples
  median_ratios <- matrix(NA, nrow = n_samples, ncol = n_samples)
  
  for (i in 1:(n_samples - 1)) {
    for (j in (i + 1):n_samples) {
      ratios <- log_intensities[, j] - log_intensities[, i]
      valid_ratios <- ratios[!is.na(ratios) & is.finite(ratios)]
      if (length(valid_ratios) > 0) {
        median_ratios[i, j] <- median(valid_ratios)
        median_ratios[j, i] <- -median_ratios[i, j]
      }
    }
  }
  
  # Iterative estimation of protein abundances
  abundances <- colMeans(log_intensities, na.rm = TRUE)
  
  for (iter in 1:10) {
    new_abundances <- rep(NA, n_samples)
    for (i in 1:n_samples) {
      estimates <- c()
      for (j in 1:n_samples) {
        if (i != j && !is.na(median_ratios[i, j]) && !is.na(abundances[j])) {
          estimates <- c(estimates, abundances[j] + median_ratios[i, j])
        }
      }
      new_abundances[i] <- if (length(estimates) > 0) median(estimates) else abundances[i]
    }
    
    new_abundances <- new_abundances - mean(new_abundances, na.rm = TRUE) + 
      mean(abundances, na.rm = TRUE)
    
    if (max(abs(new_abundances - abundances), na.rm = TRUE) < 0.001) break
    abundances <- new_abundances
  }
  
  return(2^abundances)
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3: PEAK DETECTION (from working P3_Height)
# ═══════════════════════════════════════════════════════════════════════════════

#' Estimate peak half-width using FWHM (Full Width at Half Maximum)
#' This gives adaptive peak boundaries without expensive Gaussian fitting
#' 
#' @param chrom Chromatogram vector
#' @param peak_idx Index of peak center
#' @param min_hw Minimum half-width to return (default 2)
#' @param max_hw Maximum half-width to return (default 10)
#' @return Estimated half-width in fractions
estimate_halfwidth_fwhm <- function(chrom, peak_idx, min_hw = 2L, max_hw = 10L) {
  n <- length(chrom)
  if (peak_idx < 1 || peak_idx > n) return(min_hw)
  
  peak_height <- chrom[peak_idx]
  if (!is.finite(peak_height) || peak_height <= 0) return(min_hw)
  
  half_max <- peak_height / 2
  
 # Scan left to find half-max crossing
  left_hw <- 0L
  for (i in (peak_idx - 1):1) {
    if (i < 1) break
    if (!is.finite(chrom[i])) next
    if (chrom[i] <= half_max) {
      left_hw <- as.integer(peak_idx - i)
      break
    }
    left_hw <- as.integer(peak_idx - i)
  }
  
  # Scan right to find half-max crossing
  right_hw <- 0L
  for (i in (peak_idx + 1):n) {
    if (i > n) break
    if (!is.finite(chrom[i])) next
    if (chrom[i] <= half_max) {
      right_hw <- as.integer(i - peak_idx)
      break
    }
    right_hw <- as.integer(i - peak_idx)
  }
  
  # Use average of left and right, or whichever is valid
  if (left_hw > 0 && right_hw > 0) {
    hw <- as.integer(round((left_hw + right_hw) / 2))
  } else if (left_hw > 0) {
    hw <- left_hw
  } else if (right_hw > 0) {
    hw <- right_hw
  } else {
    hw <- min_hw
  }
  
  # Clamp to bounds
  return(max(min_hw, min(max_hw, hw)))
}

#' Find local maxima in chromatogram
find_local_maxima_na <- function(chrom) {
  n <- length(chrom)
  if (n < 3 || all(is.na(chrom))) return(integer(0))
  
  peaks <- c()
  for (i in 2:(n - 1)) {
    if (is.na(chrom[i]) || !is.finite(chrom[i])) next
    
    left_val <- chrom[i - 1]
    right_val <- chrom[i + 1]
    
    left_ok <- is.finite(left_val) && chrom[i] > left_val
    right_ok <- is.finite(right_val) && chrom[i] >= right_val
    
    if (left_ok && right_ok) peaks <- c(peaks, i)
  }
  
  # Check edges
  if (n >= 2) {
    if (is.finite(chrom[1]) && is.finite(chrom[2]) && chrom[1] > chrom[2]) {
      peaks <- c(peaks, 1)
    }
    if (is.finite(chrom[n]) && is.finite(chrom[n - 1]) && chrom[n] > chrom[n - 1]) {
      peaks <- c(peaks, n)
    }
  }
  
  return(sort(unique(peaks)))
}

#' Calculate prominence of a peak
prominence_na <- function(chrom, peak_idx) {
  n <- length(chrom)
  if (n == 0 || peak_idx < 1 || peak_idx > n) return(NA_real_)
  if (!is.finite(chrom[peak_idx])) return(NA_real_)
  
  peak_h <- chrom[peak_idx]
  
  left_min <- peak_h
  if (peak_idx > 1) {
    for (i in (peak_idx - 1):1) {
      if (!is.finite(chrom[i])) next
      if (chrom[i] > peak_h) break
      left_min <- min(left_min, chrom[i])
    }
  }
  
  right_min <- peak_h
  if (peak_idx < n) {
    for (i in (peak_idx + 1):n) {
      if (!is.finite(chrom[i])) next
      if (chrom[i] > peak_h) break
      right_min <- min(right_min, chrom[i])
    }
  }
  
  baseline <- max(left_min, right_min)
  return(peak_h - baseline)
}

#' Detect peaks in chromatogram with FWHM-based width estimation
detect_peaks_local_maxima <- function(chromatogram, peak_params) {
  chrom <- as.numeric(chromatogram)
  n <- length(chrom)
  
  if (n == 0 || all(is.na(chrom))) {
    return(data.frame(peak_id = integer(), center = integer(), 
                      height = numeric(), prominence = numeric(),
                      half_width = integer()))
  }
  
  min_height_fraction <- peak_params$min_height_fraction %||% 0.10
  min_prominence <- peak_params$min_prominence %||% 0.15
  min_distance <- peak_params$min_distance %||% 3
  
  max_intensity <- max(chrom, na.rm = TRUE)
  if (!is.finite(max_intensity) || max_intensity <= 0) {
    return(data.frame(peak_id = integer(), center = integer(), 
                      height = numeric(), prominence = numeric(),
                      half_width = integer()))
  }
  
  peaks <- find_local_maxima_na(chrom)
  if (length(peaks) == 0) {
    peaks <- which.max(ifelse(is.finite(chrom), chrom, -Inf))
  }
  
  # Filter by height
  height_threshold <- max_intensity * min_height_fraction
  peaks <- peaks[!is.na(chrom[peaks]) & chrom[peaks] >= height_threshold]
  
  if (length(peaks) == 0) {
    return(data.frame(peak_id = integer(), center = integer(), 
                      height = numeric(), prominence = numeric(),
                      half_width = integer()))
  }
  
  # Calculate prominence
  prom <- vapply(peaks, function(p) prominence_na(chrom, p), numeric(1))
  prom_rel <- prom / max_intensity
  
  # Filter by prominence
  valid <- !is.na(prom_rel) & prom_rel >= min_prominence
  peaks <- peaks[valid]
  prom <- prom[valid]
  
  if (length(peaks) == 0) {
    return(data.frame(peak_id = integer(), center = integer(), 
                      height = numeric(), prominence = numeric(),
                      half_width = integer()))
  }
  
  # Enforce minimum distance
  if (length(peaks) > 1) {
    order_by_height <- order(chrom[peaks], decreasing = TRUE)
    kept <- logical(length(peaks))
    
    for (idx in order_by_height) {
      pos <- peaks[idx]
      if (!any(kept & abs(peaks - pos) < min_distance)) {
        kept[idx] <- TRUE
      }
    }
    peaks <- peaks[kept]
    prom <- prom[kept]
  }
  
  # Estimate half-width for each peak using FWHM
  half_widths <- vapply(peaks, function(p) estimate_halfwidth_fwhm(chrom, p), integer(1))
  
  ord <- order(peaks)
  data.frame(
    peak_id = seq_along(peaks),
    center = peaks[ord],
    height = chrom[peaks[ord]],
    prominence = prom[ord],
    half_width = half_widths[ord]
  )
}

#' Quantify peak using local sum with adaptive or fixed boundaries
#' 
#' @param chromatogram Chromatogram vector
#' @param center Reference peak center
#' @param quant_window Fixed window size (used if adaptive_window = FALSE or half_width is NA)
#' @param position_tolerance Max distance to search for actual peak
#' @param half_width FWHM-based half-width from reference peak (for adaptive boundaries)
#' @param adaptive_window If TRUE, use half_width * width_multiplier instead of fixed quant_window
#' @param width_multiplier Multiplier for half_width (default 1.5 ≈ 90% of Gaussian area)
#' @return List with local_sum, max_height, actual_center, peak_detected, used_window
quantify_peak_local_sum <- function(chromatogram, center, quant_window = 5, 
                                     position_tolerance = 3, half_width = NA,
                                     adaptive_window = TRUE, width_multiplier = 1.5) {
  chrom <- as.numeric(chromatogram)
  n <- length(chrom)
  
  if (is.na(center) || center < 1 || center > n) {
    return(list(local_sum = NA_real_, max_height = NA_real_, 
                actual_center = NA_real_, peak_detected = FALSE,
                used_window = NA_integer_))
  }
  
  # Find actual peak within tolerance
  lo <- clamp(round(center) - position_tolerance, 1, n)
  hi <- clamp(round(center) + position_tolerance, 1, n)
  idx <- lo:hi
  idx <- idx[is_finite_num(chrom[idx])]
  
  if (length(idx) == 0) {
    return(list(local_sum = NA_real_, max_height = NA_real_, 
                actual_center = NA_real_, peak_detected = FALSE,
                used_window = NA_integer_))
  }
  
  actual_center <- idx[which.max(chrom[idx])]
  
  # Determine window size: adaptive (based on FWHM) or fixed
  if (adaptive_window && !is.na(half_width) && half_width > 0) {
    # Adaptive: use FWHM-based half_width × multiplier
    # width_multiplier = 1.5 captures ~90% of Gaussian area (similar to z=1.645)
    window_size <- as.integer(ceiling(half_width * width_multiplier))
  } else {
    # Fixed window
    window_size <- as.integer(quant_window)
  }
  
  # Define window around actual center
  start_idx <- max(1, actual_center - window_size)
  end_idx <- min(n, actual_center + window_size)
  
  window_data <- chrom[start_idx:end_idx]
  valid_data <- window_data[is_finite_num(window_data)]
  
  if (length(valid_data) == 0) {
    return(list(local_sum = NA_real_, max_height = NA_real_, 
                actual_center = actual_center, peak_detected = FALSE,
                used_window = window_size))
  }
  
  list(
    local_sum = sum(valid_data),
    max_height = max(valid_data),
    actual_center = actual_center,
    peak_detected = TRUE,
    used_window = window_size
  )
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3: LIMMA STATISTICS (EXACT copy from working P3_Height)
# ═══════════════════════════════════════════════════════════════════════════════

#' IBMT: Intensity-based Moderated T-statistic
calculate_ibmt <- function(fit, coef = 1) {
  sigma <- fit$sigma
  df_residual <- fit$df.residual
  
  avg_intensity <- tryCatch({
    if (!is.null(fit$Amean)) {
      if (is.matrix(fit$Amean)) {
        rowMeans(fit$Amean, na.rm = TRUE)
      } else {
        as.numeric(fit$Amean)
      }
    } else {
      rep(mean(sigma, na.rm = TRUE), length(sigma))
    }
  }, error = function(e) {
    rep(NA, length(sigma))
  })
  
  if (all(is.na(avg_intensity))) {
    return(eBayes(fit))
  }
  
  valid <- !is.na(sigma) & !is.na(avg_intensity) & sigma > 0
  if (sum(valid) < 10) {
    return(eBayes(fit))
  }
  
  fit_loess <- tryCatch({
    loess(log(sigma[valid]^2) ~ avg_intensity[valid], span = 0.5)
  }, error = function(e) NULL)
  
  if (is.null(fit_loess)) {
    return(eBayes(fit))
  }
  
  s2_prior <- exp(predict(fit_loess, newdata = avg_intensity))
  s2_prior[is.na(s2_prior)] <- median(sigma^2, na.rm = TRUE)
  
  d0 <- 4
  s2_post <- (d0 * s2_prior + df_residual * sigma^2) / (d0 + df_residual)
  
  se <- fit$stdev.unscaled[, coef] * sqrt(s2_post)
  t_stat <- fit$coefficients[, coef] / se
  
  df_total <- d0 + df_residual
  p_value <- 2 * pt(abs(t_stat), df = df_total, lower.tail = FALSE)
  
  fit$t <- matrix(t_stat, ncol = 1)
  fit$p.value <- matrix(p_value, ncol = 1)
  fit$s2.post <- s2_post
  fit$df.total <- df_total
  
  colnames(fit$t) <- colnames(fit$p.value) <- colnames(fit$coefficients)[coef]
  
  return(fit)
}

#' Run Limma differential analysis with IBMT
run_limma_analysis <- function(value_matrix, groups, ref_condition = NULL, 
                                treat_condition = NULL, use_ibmt = TRUE) {
  
  if (!is.null(ref_condition) && !is.null(treat_condition)) {
    groups_factor <- factor(groups, levels = c(ref_condition, treat_condition))
  } else {
    groups_factor <- factor(groups)
  }
  
  design <- model.matrix(~ 0 + groups_factor)
  colnames(design) <- levels(groups_factor)
  
  fit <- lmFit(value_matrix, design)
  
  group_levels <- levels(groups_factor)
  if (length(group_levels) != 2) {
    stop("Expected exactly 2 groups, got: ", paste(group_levels, collapse = ", "))
  }
  
  contrast_formula <- paste0(group_levels[2], " - ", group_levels[1])
  message("  Contrast: ", contrast_formula)
  contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
  
  fit2 <- contrasts.fit(fit, contrast_matrix)
  
  fit2 <- tryCatch({
    if (use_ibmt) {
      calculate_ibmt(fit2)
    } else {
      eBayes(fit2)
    }
  }, error = function(e) {
    message("  Warning: IBMT failed, falling back to eBayes: ", e$message)
    eBayes(fit2)
  })
  
  results <- tryCatch({
    topTable(fit2, number = Inf, sort.by = "none")
  }, error = function(e) {
    message("  Warning: topTable failed, creating manual results: ", e$message)
    data.frame(
      logFC = fit2$coefficients[, 1],
      AveExpr = rowMeans(value_matrix, na.rm = TRUE),
      t = if (!is.null(fit2$t)) fit2$t[, 1] else NA,
      P.Value = if (!is.null(fit2$p.value)) fit2$p.value[, 1] else NA,
      adj.P.Val = if (!is.null(fit2$p.value)) p.adjust(fit2$p.value[, 1], method = "BH") else NA,
      row.names = rownames(value_matrix)
    )
  })
  
  if ("ID" %in% colnames(results)) {
    results$feature_id <- results$ID
  } else {
    results$feature_id <- rownames(results)
  }
  
  return(results)
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4: MAIN PIPELINE FUNCTION
# ═══════════════════════════════════════════════════════════════════════════════

#' Run protein-level SEC-MS pipeline
#' 
#' @param traces_list List of peptide traces with elements named "{condition}_{replicate}"
#' @param protein_map Optional: data.frame with 'peptide_id' and 'protein_id' columns
#' @param design_matrix Optional. If NULL, created from traces_list names
#' @param ref_condition Reference condition name (e.g., "ctrl")
#' @param treat_condition Treatment condition name (e.g., "ATP")
#' @param protein_col Column name for protein ID in traces (auto-detect if NULL)
#' @param min_peptides Minimum peptides required per protein
#' @param aggregation_method Method for aggregating peptides to proteins:
#'        "maxlfq" - MaxLFQ across samples (most accurate, slower)
#'        "median" - Median of all peptides
#'        "sum" - Sum of all peptides
#'        "top3", "top5", "top10" - Sum of top N peptides by total intensity
#'        "top10median" - Median of top 10 peptides (min 4 required)
#' @param pvalue_adjustment Method for p-value adjustment:
#'        "BH" - Benjamini-Hochberg (assumes independence, can be conservative for CF-MS)
#'        "qvalue" - Storey's q-value (estimates pi0, often less conservative)
#'        "none" - No adjustment (report raw p-values only)
#' @param n_permutations Number of permutations for empirical FDR (0 = skip, recommended: 100-1000)
#' @param peak_params Parameters for peak detection
#' @param quant_window Fixed window size for quantification (used if adaptive_window = FALSE)
#' @param adaptive_window If TRUE, use FWHM-based adaptive boundaries instead of fixed
#' @param width_multiplier Multiplier for FWHM half-width (1.5 ≈ 90% of peak area)
#' @param position_tolerance Maximum position difference for peak matching
#' @param detect_new_peaks Whether to detect peaks in treatment not in control
#' @param impute_value Value to use for missing peaks (enables limma for appearing/disappearing peaks)
#' @param flag_large_shifts Whether to flag large position shifts
#' @param large_shift_multiplier Multiplier for position_tolerance to flag large shifts
#' @param output_dir Output directory
#' @param n_trace_plots Number of significant proteins to plot
#' @param custom_ids Optional: vector of protein_ids to plot
#'
run_protein_level_pipeline <- function(
    traces_list,
    protein_map = NULL,
    design_matrix = NULL,
    ref_condition,
    treat_condition,
    protein_col = NULL,
    min_peptides = 3,
    aggregation_method = c("maxlfq", "median", "top3", "top5", "top10", "top10median", "sum"),
    pvalue_adjustment = c("BH", "qvalue", "none"),
    n_permutations = 0,
    peak_params = list(
      min_height_fraction = 0.10,
      min_prominence = 0.15,
      min_distance = 3
    ),
    quant_window = 5,
    adaptive_window = TRUE,
    width_multiplier = 1.5,
    position_tolerance = 3,
    detect_new_peaks = TRUE,
    impute_value = 10,
    flag_large_shifts = TRUE,
    large_shift_multiplier = 3,
    output_dir = "protein_level_results",
    n_trace_plots = 20,
    custom_ids = NULL
) {
  
  aggregation_method <- match.arg(aggregation_method)
  pvalue_adjustment <- match.arg(pvalue_adjustment)
  
  message("\n", paste(rep("═", 70), collapse = ""))
  message("PROTEIN-LEVEL SEC-MS PIPELINE")
  message("Treatment: ", treat_condition, " vs Reference: ", ref_condition)
  message(paste(rep("═", 70), collapse = ""))
  # Parse topN if needed
  topN <- NULL
  topN_median <- FALSE
  if (aggregation_method == "top10median") {
    topN <- 10L
    topN_median <- TRUE
  } else if (grepl("^top", aggregation_method)) {
    topN <- as.integer(sub("top", "", aggregation_method))
  }
  
  message("\n  This simplified pipeline:")
  if (aggregation_method == "top10median") {
    message("  - Aggregation: TOP 10 peptides (min 4) → MEDIAN per protein per fraction")
  } else if (!is.null(topN)) {
    message("  - Aggregation: TOP", topN, " peptides (by intensity) → SUM per protein per fraction")
  } else {
    message("  - Aggregation: ", toupper(aggregation_method), " per protein per fraction")
  }
  message("  - Peak boundaries: ", ifelse(adaptive_window, 
          paste0("ADAPTIVE (FWHM × ", width_multiplier, ")"), 
          paste0("FIXED (±", quant_window, " fractions)")))
  message("  - P-value adjustment: ", toupper(pvalue_adjustment), 
          ifelse(pvalue_adjustment == "qvalue", " (Storey's q-value, less conservative)", 
                 ifelse(pvalue_adjustment == "BH", " (Benjamini-Hochberg)", " (raw p-values only)")))
  if (n_permutations > 0) {
    message("  - Permutation FDR: ", n_permutations, " permutations (empirical, accounts for correlations)")
  }
  message("  - NO proteoform clustering")
  message("  - Peak detection on CONTROL protein chromatograms")
  message("  - Quantify all samples at reference peak positions")
  message("  - Shows peptide traces in plots")
  message(paste(rep("═", 70), collapse = ""), "\n")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)
  
  # ═══════════════════════════════════════════════════════════════════════════
  # STEP 1: Setup and validation
  # ═══════════════════════════════════════════════════════════════════════════
  message("STEP 1: Setting up...")
  
  message("  Available samples: ", paste(names(traces_list), collapse = ", "))
  
  if (is.null(design_matrix)) {
    design_matrix <- data.frame(
      sample_id = names(traces_list),
      stringsAsFactors = FALSE
    )
    parsed <- lapply(design_matrix$sample_id, parse_sample_name)
    design_matrix$condition <- sapply(parsed, `[[`, "condition")
    design_matrix$replicate <- sapply(parsed, `[[`, "replicate")
    
    message("  Detected conditions: ", paste(unique(design_matrix$condition), collapse = ", "))
  }
  
  ref_samples <- design_matrix$sample_id[design_matrix$condition == ref_condition]
  treat_samples <- design_matrix$sample_id[design_matrix$condition == treat_condition]
  
  message("  Reference samples (", ref_condition, "): ", 
          ifelse(length(ref_samples) > 0, paste(ref_samples, collapse = ", "), "NONE FOUND!"))
  message("  Treatment samples (", treat_condition, "): ", 
          ifelse(length(treat_samples) > 0, paste(treat_samples, collapse = ", "), "NONE FOUND!"))
  
  if (length(ref_samples) == 0 || length(treat_samples) == 0) {
    stop("Could not find samples for both conditions. Check sample names.")
  }
  
  # Get protein mapping
  if (is.null(protein_map)) {
    protein_map <- get_protein_mapping(traces_list, protein_col = protein_col)
  }
  
  if (nrow(protein_map) == 0) {
    stop("No protein mappings found.")
  }
  
  message("  Found ", nrow(protein_map), " peptide-protein mappings")
  
  # Get proteins with enough peptides
  protein_peptide_counts <- protein_map %>%
    group_by(protein_id) %>%
    summarise(n_peptides = n(), .groups = "drop") %>%
    filter(n_peptides >= min_peptides)
  
  proteins <- protein_peptide_counts$protein_id
  message("  Proteins with >= ", min_peptides, " peptides: ", length(proteins), "\n")
  
  # Create peptide list per protein
  protein_peptide_list <- split(protein_map$peptide_id, protein_map$protein_id)
  
  # ═══════════════════════════════════════════════════════════════════════════
  # STEP 2: Build protein chromatograms
  # ═══════════════════════════════════════════════════════════════════════════
  if (topN_median) {
    agg_label <- "TOP10 MEDIAN"
  } else if (!is.null(topN)) {
    agg_label <- paste0("TOP", topN, " SUM")
  } else {
    agg_label <- toupper(aggregation_method)
  }
  message("STEP 2: Building protein chromatograms (", agg_label, " per fraction)...")
  
  trace_mats <- lapply(traces_list, function(x) extract_trace_matrix(x)$matrix)
  n_fractions <- ncol(trace_mats[[1]])
  fraction_names <- colnames(trace_mats[[1]])
  sample_names <- names(traces_list)
  
  all_peptide_mats <- trace_mats
  
  # Initialize protein chromatograms for all samples
  protein_chromatograms <- list()
  for (sample_id in sample_names) {
    prot_mat <- matrix(NA_real_, nrow = length(proteins), ncol = n_fractions)
    rownames(prot_mat) <- proteins
    colnames(prot_mat) <- fraction_names
    protein_chromatograms[[sample_id]] <- prot_mat
  }
  
  if (aggregation_method == "maxlfq") {
    # MaxLFQ: Apply per protein per fraction ACROSS all samples
    pb <- txtProgressBar(min = 0, max = length(proteins), style = 3)
    
    for (prot_idx in seq_along(proteins)) {
      setTxtProgressBar(pb, prot_idx)
      prot <- proteins[prot_idx]
      prot_peptides <- protein_peptide_list[[prot]]
      
      for (frac in 1:n_fractions) {
        # Build peptide x sample matrix for this fraction
        pep_sample_mat <- matrix(NA_real_, 
                                  nrow = length(prot_peptides), 
                                  ncol = length(sample_names))
        rownames(pep_sample_mat) <- prot_peptides
        colnames(pep_sample_mat) <- sample_names
        
        for (s_idx in seq_along(sample_names)) {
          sample_id <- sample_names[s_idx]
          trace_mat <- trace_mats[[sample_id]]
          
          peps_in_sample <- prot_peptides[prot_peptides %in% rownames(trace_mat)]
          if (length(peps_in_sample) > 0) {
            pep_sample_mat[peps_in_sample, s_idx] <- trace_mat[peps_in_sample, frac]
          }
        }
        
        # Check if we have enough data
        valid_peptides <- rowSums(!is.na(pep_sample_mat) & pep_sample_mat > 0) > 0
        pep_sample_mat <- pep_sample_mat[valid_peptides, , drop = FALSE]
        
        if (nrow(pep_sample_mat) >= min_peptides) {
          # Apply MaxLFQ across samples
          protein_abundances <- maxlfq_across_samples(pep_sample_mat)
          
          for (s_idx in seq_along(sample_names)) {
            protein_chromatograms[[sample_names[s_idx]]][prot, frac] <- protein_abundances[s_idx]
          }
        }
      }
    }
    close(pb)
    
  } else {
    # All other methods: Apply per sample
    # Methods: median, sum, top3, top5, top10, top10median
    pb <- txtProgressBar(min = 0, max = length(sample_names), style = 3)
    
    for (s_idx in seq_along(sample_names)) {
      setTxtProgressBar(pb, s_idx)
      sample_id <- sample_names[s_idx]
      trace_mat <- trace_mats[[sample_id]]
      
      for (prot in proteins) {
        prot_peptides <- protein_peptide_list[[prot]]
        prot_peptides <- prot_peptides[prot_peptides %in% rownames(trace_mat)]
        
        # For top10median: require minimum 4 peptides
        min_peps_required <- if (topN_median) 4 else min_peptides
        
        if (length(prot_peptides) >= min_peps_required) {
          pep_mat <- trace_mat[prot_peptides, , drop = FALSE]
          
          # For TopN methods: rank peptides by total intensity across all fractions
          if (!is.null(topN)) {
            pep_totals <- rowSums(pep_mat, na.rm = TRUE)
            n_to_select <- min(topN, length(pep_totals))
            top_peps <- names(sort(pep_totals, decreasing = TRUE))[1:n_to_select]
            pep_mat <- pep_mat[top_peps, , drop = FALSE]
          }
          
          for (frac in 1:n_fractions) {
            frac_intensities <- pep_mat[, frac]
            valid <- is_finite_num(frac_intensities) & frac_intensities > 0
            
            # ── FIX: enforce min_peptides PER FRACTION (not just >= 1) ──
            # Previously: if (sum(valid) >= 1) → aggregated even with a single
            # peptide, making median/top10median/top5/sum all identical.
            # Now: require at least min_peps_required valid peptides in each
            # fraction, consistent with how maxlfq handles it (line 909).
            # For top10median this is 4; for all others it is min_peptides.
            if (sum(valid) >= min_peps_required) {
              if (aggregation_method == "median") {
                protein_chromatograms[[sample_id]][prot, frac] <- median(frac_intensities[valid])
              } else if (aggregation_method == "sum") {
                protein_chromatograms[[sample_id]][prot, frac] <- sum(frac_intensities[valid])
              } else if (topN_median) {
                # TopN median: median of top N peptides
                protein_chromatograms[[sample_id]][prot, frac] <- median(frac_intensities[valid])
              } else if (!is.null(topN)) {
                # TopN sum: sum of top N peptides
                protein_chromatograms[[sample_id]][prot, frac] <- sum(frac_intensities[valid])
              }
            }
          }
        }
      }
    }
    close(pb)
  }
  
  message("\n  Built chromatograms for ", length(proteins), " proteins")
  
  # Diagnostic: Show some statistics to verify aggregation method worked
  sample1 <- sample_names[1]
  prot1 <- proteins[1]
  prot1_chrom <- protein_chromatograms[[sample1]][prot1, ]
  prot1_max <- max(prot1_chrom, na.rm = TRUE)
  prot1_sum <- sum(prot1_chrom, na.rm = TRUE)
  
  # Count peptides for this protein
  prot1_peptides <- protein_peptide_list[[prot1]]
  prot1_peptides_in_sample <- prot1_peptides[prot1_peptides %in% rownames(trace_mats[[sample1]])]
  n_peps_available <- length(prot1_peptides_in_sample)
  if (!is.null(topN)) {
    n_peps_used <- min(topN, n_peps_available)
  } else {
    n_peps_used <- n_peps_available
  }
  
  message("  Verification (", prot1, " in ", sample1, "):")
  message("    Peptides available: ", n_peps_available)
  message("    Peptides used: ", n_peps_used, " (", agg_label, ")")
  message("    Max intensity: ", round(prot1_max, 1))
  message("    Sum across fractions: ", round(prot1_sum, 1))
  
  # Show overall peptide usage summary
  peps_per_prot <- sapply(proteins, function(p) {
    peps <- protein_peptide_list[[p]]
    peps_in <- peps[peps %in% rownames(trace_mats[[sample1]])]
    length(peps_in)
  })
  
  message("  Overall peptide counts (in ", sample1, "):")
  message("    Min: ", min(peps_per_prot), ", Median: ", median(peps_per_prot), 
          ", Max: ", max(peps_per_prot))
  if (!is.null(topN)) {
    pct_affected <- 100 * mean(peps_per_prot > topN)
    message("    Proteins with >", topN, " peptides (affected by TopN): ", 
            round(pct_affected, 1), "%")
  }
  if (topN_median) {
    pct_excluded <- 100 * mean(peps_per_prot < 4)
    message("    Proteins with <4 peptides (excluded by top10median): ", 
            round(pct_excluded, 1), "%")
  }
  message("")
  
  # ═══════════════════════════════════════════════════════════════════════════
  # STEP 3: Detect peaks from CONTROL samples
  # ═══════════════════════════════════════════════════════════════════════════
  message("STEP 3: Detecting reference peaks from CONTROL samples...")
  
  control_avg <- matrix(NA_real_, nrow = length(proteins), ncol = n_fractions)
  rownames(control_avg) <- proteins
  colnames(control_avg) <- fraction_names
  
  for (prot in proteins) {
    prot_chroms <- lapply(ref_samples, function(s) {
      if (prot %in% rownames(protein_chromatograms[[s]])) {
        protein_chromatograms[[s]][prot, ]
      } else {
        rep(NA, n_fractions)
      }
    })
    prot_mat <- do.call(rbind, prot_chroms)
    control_avg[prot, ] <- colMeans(prot_mat, na.rm = TRUE)
  }
  control_avg[is.nan(control_avg)] <- NA
  
  reference_peaks <- list()
  skipped_na <- 0
  skipped_error <- 0
  
  for (prot in proteins) {
    chrom <- control_avg[prot, ]
    
    if (all(is.na(chrom)) || !any(is.finite(chrom) & chrom > 0)) {
      skipped_na <- skipped_na + 1
      next
    }
    
    tryCatch({
      peaks <- detect_peaks_local_maxima(chrom, peak_params)
      
      if (nrow(peaks) > 0) {
        peaks$protein_id <- prot
        reference_peaks[[prot]] <- peaks
      }
    }, error = function(e) {
      skipped_error <<- skipped_error + 1
    })
  }
  
  if (skipped_na > 0 || skipped_error > 0) {
    message("  Skipped: ", skipped_na, " proteins with no signal, ", 
            skipped_error, " proteins with errors")
  }
  
  reference_peaks_df <- bind_rows(reference_peaks)
  message("  Detected peaks for ", length(reference_peaks), " proteins")
  message("  Peak distribution:")
  print(table(sapply(reference_peaks, nrow)))
  
  # ═══════════════════════════════════════════════════════════════════════════
  # STEP 4: Quantify ALL samples at reference peak positions
  # ═══════════════════════════════════════════════════════════════════════════
  message("\nSTEP 4: Quantifying all samples at reference peak positions...")
  if (adaptive_window) {
    message("  Using ADAPTIVE boundaries (FWHM × ", width_multiplier, ")")
  } else {
    message("  Using FIXED boundaries (±", quant_window, " fractions)")
  }
  
  all_peak_quant <- list()
  
  for (sample_id in sample_names) {
    parsed <- parse_sample_name(sample_id)
    prot_mat <- protein_chromatograms[[sample_id]]
    
    for (prot in names(reference_peaks)) {
      ref_peaks <- reference_peaks[[prot]]
      
      if (!prot %in% rownames(prot_mat)) {
        # Protein missing in this sample - keep NA (will be imputed for limma)
        for (pk_idx in seq_len(nrow(ref_peaks))) {
          all_peak_quant[[length(all_peak_quant) + 1]] <- data.frame(
            protein_id = prot,
            peak_id = ref_peaks$peak_id[pk_idx],
            sample_id = sample_id,
            condition = parsed$condition,
            replicate = parsed$replicate,
            ref_center = ref_peaks$center[pk_idx],
            ref_half_width = ref_peaks$half_width[pk_idx],
            actual_center = NA_real_,
            local_sum = NA_real_,
            max_height = NA_real_,
            used_window = NA_integer_,
            peak_detected = FALSE,
            stringsAsFactors = FALSE
          )
        }
        next
      }
      
      chrom <- prot_mat[prot, ]
      
      for (pk_idx in seq_len(nrow(ref_peaks))) {
        ref_center <- ref_peaks$center[pk_idx]
        ref_hw <- ref_peaks$half_width[pk_idx]
        
        quant <- quantify_peak_local_sum(
          chromatogram = chrom, 
          center = ref_center, 
          quant_window = quant_window, 
          position_tolerance = position_tolerance,
          half_width = ref_hw,
          adaptive_window = adaptive_window,
          width_multiplier = width_multiplier
        )
        
        all_peak_quant[[length(all_peak_quant) + 1]] <- data.frame(
          protein_id = prot,
          peak_id = ref_peaks$peak_id[pk_idx],
          sample_id = sample_id,
          condition = parsed$condition,
          replicate = parsed$replicate,
          ref_center = ref_center,
          ref_half_width = ref_hw,
          actual_center = quant$actual_center,
          local_sum = quant$local_sum,
          max_height = quant$max_height,
          used_window = quant$used_window,
          peak_detected = quant$peak_detected,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  peak_quant <- bind_rows(all_peak_quant)
  
  # Create feature_id (KEY for limma!)
  peak_quant$feature_id <- paste(peak_quant$protein_id, "peak", peak_quant$peak_id, sep = "_")
  
  message("  Quantified ", nrow(peak_quant), " protein-peak-sample combinations")
  message("  Unique features (protein_peak): ", n_distinct(peak_quant$feature_id))
  if (adaptive_window) {
    message("  Adaptive window sizes: min=", min(peak_quant$used_window, na.rm = TRUE),
            ", median=", median(peak_quant$used_window, na.rm = TRUE),
            ", max=", max(peak_quant$used_window, na.rm = TRUE))
  }
  message("")
  
  # ═══════════════════════════════════════════════════════════════════════════
  # STEP 5: Calculate differential metrics
  # ═══════════════════════════════════════════════════════════════════════════
  message("STEP 5: Calculating differential metrics...")
  
  peak_summary <- peak_quant %>%
    filter(condition %in% c(ref_condition, treat_condition)) %>%
    group_by(feature_id, protein_id, peak_id, condition) %>%
    summarise(
      mean_local_sum = mean(local_sum, na.rm = TRUE),
      mean_max_height = mean(max_height, na.rm = TRUE),
      mean_center = mean(actual_center, na.rm = TRUE),
      n_detected = sum(peak_detected, na.rm = TRUE),
      n_replicates = n(),
      .groups = "drop"
    )
  
  ref_summary <- peak_summary %>%
    filter(condition == ref_condition) %>%
    select(feature_id, protein_id, peak_id, 
           ref_sum = mean_local_sum, ref_height = mean_max_height,
           ref_center = mean_center, ref_detected = n_detected, ref_n = n_replicates)
  
  treat_summary <- peak_summary %>%
    filter(condition == treat_condition) %>%
    select(feature_id, protein_id, peak_id,
           treat_sum = mean_local_sum, treat_height = mean_max_height,
           treat_center = mean_center, treat_detected = n_detected, treat_n = n_replicates)
  
  diff_results <- ref_summary %>%
    full_join(treat_summary, by = c("feature_id", "protein_id", "peak_id")) %>%
    mutate(
      log2FC_sum = log2((treat_sum + impute_value) / (ref_sum + impute_value)),
      log2FC_height = log2((treat_height + impute_value) / (ref_height + impute_value)),
      delta_position = treat_center - ref_center,
      peak_status = case_when(
        is.na(ref_detected) | ref_detected == 0 ~ "treatment_only",
        is.na(treat_detected) | treat_detected == 0 ~ "control_only",
        flag_large_shifts & abs(delta_position) > position_tolerance * large_shift_multiplier ~ "large_shift",
        TRUE ~ "matched"
      )
    )
  
  # ═══════════════════════════════════════════════════════════════════════════
  # STEP 6: Limma analysis (EXACT method from working P3_Height)
  # ═══════════════════════════════════════════════════════════════════════════
  message("\nSTEP 6: Running Limma differential analysis...")
  
  # Create matrix using base R (avoids dplyr issues with column names)
  limma_input <- peak_quant[peak_quant$condition %in% c(ref_condition, treat_condition), 
                             c("feature_id", "sample_id", "local_sum")]
  
  # Handle duplicates
  limma_input <- aggregate(local_sum ~ feature_id + sample_id, data = limma_input, 
                            FUN = mean, na.rm = TRUE)
  
  # Pivot to wide format using base R
  quant_wide <- reshape(limma_input, 
                        idvar = "feature_id", 
                        timevar = "sample_id", 
                        direction = "wide")
  
  colnames(quant_wide) <- gsub("^local_sum\\.", "", colnames(quant_wide))
  
  rownames(quant_wide) <- quant_wide$feature_id
  quant_wide$feature_id <- NULL
  
  sample_order <- c(ref_samples, treat_samples)
  sample_order <- sample_order[sample_order %in% colnames(quant_wide)]
  quant_mat <- as.matrix(quant_wide[, sample_order, drop = FALSE])
  
  groups <- c(rep(ref_condition, sum(sample_order %in% ref_samples)),
              rep(treat_condition, sum(sample_order %in% treat_samples)))
  
  message("  Limma matrix: ", nrow(quant_mat), " features x ", ncol(quant_mat), " samples")
  message("  Groups: ", paste(groups, collapse = ", "))
  
  # Impute NA values (enables testing control-only and treatment-only peaks!)
  n_na_before <- sum(is.na(quant_mat))
  if (n_na_before > 0) {
    quant_mat[is.na(quant_mat)] <- impute_value
    message("  Imputed ", n_na_before, " NA values with impute_value = ", impute_value)
  }
  
  # Log transform
  quant_mat_log <- log2(quant_mat + impute_value)
  
  # Remove zero variance rows
  valid_rows <- apply(quant_mat_log, 1, function(x) var(x, na.rm = TRUE) > 0)
  quant_mat_log <- quant_mat_log[valid_rows, , drop = FALSE]
  
  n_zero_var <- sum(!valid_rows)
  if (n_zero_var > 0) {
    message("  Removed ", n_zero_var, " features with zero variance")
  }
  message("  After filtering: ", nrow(quant_mat_log), " features ready for Limma")
  
  # Run Limma
  if (nrow(quant_mat_log) > 0) {
    limma_results <- run_limma_analysis(quant_mat_log, groups, 
                                         ref_condition = ref_condition, 
                                         treat_condition = treat_condition,
                                         use_ibmt = TRUE)
    
    # Rename columns
    col_rename <- c("logFC" = "limma_logFC", 
                    "t" = "limma_t", 
                    "P.Value" = "limma_pvalue", 
                    "adj.P.Val" = "limma_adj_pvalue_BH")
    
    for (old_name in names(col_rename)) {
      new_name <- col_rename[old_name]
      idx <- which(colnames(limma_results) == old_name)
      if (length(idx) > 0) {
        colnames(limma_results)[idx[1]] <- new_name
      }
    }
    
    # Apply selected p-value adjustment method
    raw_pvals <- limma_results$limma_pvalue
    
    if (pvalue_adjustment == "qvalue") {
      # Storey's q-value (less conservative, estimates pi0)
      tryCatch({
        if (!requireNamespace("qvalue", quietly = TRUE)) {
          message("  NOTE: 'qvalue' package not installed. Install with:")
          message("        BiocManager::install('qvalue')")
          message("  Falling back to BH adjustment.")
          limma_results$limma_adj_pvalue <- p.adjust(raw_pvals, method = "BH")
          pvalue_adjustment <- "BH"  # Update for reporting
        } else {
          valid_pvals <- raw_pvals[!is.na(raw_pvals) & raw_pvals > 0 & raw_pvals <= 1]
          message("  Q-value input: ", length(valid_pvals), " valid p-values")
          message("  P-value range: ", round(min(valid_pvals), 6), " - ", round(max(valid_pvals), 3))
          
          qval_result <- qvalue::qvalue(valid_pvals)
          
          # Map q-values back to full results
          limma_results$limma_adj_pvalue <- NA_real_
          valid_idx <- !is.na(raw_pvals) & raw_pvals > 0 & raw_pvals <= 1
          limma_results$limma_adj_pvalue[valid_idx] <- qval_result$qvalues
          
          # Report pi0 estimate and q-value distribution
          message("  Storey's q-value: estimated pi0 = ", round(qval_result$pi0, 3),
                  " (", round((1 - qval_result$pi0) * 100, 1), "% estimated true positives)")
          
          qval_range <- range(qval_result$qvalues, na.rm = TRUE)
          message("  Q-value range: ", round(qval_range[1], 4), " - ", round(qval_range[2], 3))
          
          qval_quantiles <- quantile(qval_result$qvalues, c(0.01, 0.05, 0.1, 0.5), na.rm = TRUE)
          message("  Q-value quantiles: 1%=", round(qval_quantiles[1], 4),
                  ", 5%=", round(qval_quantiles[2], 4),
                  ", 10%=", round(qval_quantiles[3], 4),
                  ", median=", round(qval_quantiles[4], 3))
        }
      }, error = function(e) {
        message("  Warning: q-value estimation failed: ", e$message)
        message("  Falling back to BH adjustment.")
        limma_results$limma_adj_pvalue <<- p.adjust(raw_pvals, method = "BH")
        pvalue_adjustment <<- "BH"
      })
    } else if (pvalue_adjustment == "BH") {
      limma_results$limma_adj_pvalue <- limma_results$limma_adj_pvalue_BH
    } else {
      # No adjustment - use raw p-values
      limma_results$limma_adj_pvalue <- raw_pvals
      message("  NOTE: No p-value adjustment applied. Using raw p-values.")
    }
    
    # Permutation-based FDR (if requested)
    if (n_permutations > 0) {
      message("\n  Running permutation FDR (", n_permutations, " permutations)...")
      
      # Get observed test statistics
      observed_t <- limma_results$limma_t
      abs_obs_t <- abs(observed_t)
      n_features <- length(observed_t)
      
      # Show observed t-stat distribution
      obs_t_quantiles <- quantile(abs_obs_t, c(0.9, 0.95, 0.99, 1), na.rm = TRUE)
      message("  Observed |t| quantiles: 90%=", round(obs_t_quantiles[1], 2),
              ", 95%=", round(obs_t_quantiles[2], 2),
              ", 99%=", round(obs_t_quantiles[3], 2),
              ", max=", round(obs_t_quantiles[4], 2))
      
      # Run permutations and collect ALL null t-statistics
      null_t_list <- vector("list", n_permutations)
      
      pb_perm <- txtProgressBar(min = 0, max = n_permutations, style = 3)
      for (perm_i in 1:n_permutations) {
        setTxtProgressBar(pb_perm, perm_i)
        
        # Shuffle group labels
        perm_groups <- sample(groups)
        
        # Run limma on permuted data (suppress messages)
        perm_result <- tryCatch({
          suppressMessages(
            run_limma_analysis(quant_mat_log, perm_groups,
                              ref_condition = ref_condition,
                              treat_condition = treat_condition,
                              use_ibmt = TRUE)
          )
        }, error = function(e) NULL)
        
        if (!is.null(perm_result)) {
          null_t_list[[perm_i]] <- abs(perm_result$t)
        }
      }
      close(pb_perm)
      
      # Remove failed permutations
      null_t_list <- null_t_list[!sapply(null_t_list, is.null)]
      n_valid_perms <- length(null_t_list)
      message("\n  Valid permutations: ", n_valid_perms, "/", n_permutations)
      
      if (n_valid_perms > 0) {
        # Show null t-stat distribution (pooled across all permutations)
        all_null_t <- unlist(null_t_list)
        null_t_quantiles <- quantile(all_null_t, c(0.9, 0.95, 0.99, 1), na.rm = TRUE)
        message("  Null |t| quantiles: 90%=", round(null_t_quantiles[1], 2),
                ", 95%=", round(null_t_quantiles[2], 2),
                ", 99%=", round(null_t_quantiles[3], 2),
                ", max=", round(null_t_quantiles[4], 2))
        
        # Show how many null features exceed various thresholds (average per perm)
        for (thresh in c(2, 3, 4, 5)) {
          avg_null_above <- mean(sapply(null_t_list, function(x) sum(x >= thresh, na.rm = TRUE)))
          obs_above <- sum(abs_obs_t >= thresh, na.rm = TRUE)
          message("  At |t| >= ", thresh, ": observed=", obs_above, 
                  ", null_avg=", round(avg_null_above, 1),
                  ", ratio=", round(avg_null_above / max(obs_above, 1), 3))
        }
        
        # Calculate FDR for each feature using its |t| as threshold
        limma_results$perm_fdr <- sapply(seq_along(abs_obs_t), function(i) {
          if (is.na(abs_obs_t[i])) return(NA_real_)
          
          threshold <- abs_obs_t[i]
          
          # R(t): Count observed features at or above this threshold
          R <- sum(abs_obs_t >= threshold, na.rm = TRUE)
          
          # V(t): Average count of features at or above threshold in null
          V <- mean(sapply(null_t_list, function(null_t) {
            sum(null_t >= threshold, na.rm = TRUE)
          }))
          
          # FDR = V / R
          if (R == 0) return(1)
          fdr <- V / R
          min(fdr, 1)
        })
        
        # Monotonize FDR (ensure non-increasing with |t|, i.e., non-decreasing with p-value)
        ord <- order(abs_obs_t, decreasing = TRUE)  # Sort by |t| descending
        fdr_sorted <- limma_results$perm_fdr[ord]
        fdr_monotone <- cummax(fdr_sorted)  # Ensure non-decreasing
        limma_results$perm_fdr[ord] <- fdr_monotone
        
        message("  Permutation FDR complete.")
        n_sig_perm <- sum(limma_results$perm_fdr < 0.05, na.rm = TRUE)
        n_sig_perm_10 <- sum(limma_results$perm_fdr < 0.10, na.rm = TRUE)
        n_sig_perm_20 <- sum(limma_results$perm_fdr < 0.20, na.rm = TRUE)
        message("  Significant at perm_fdr < 0.05: ", n_sig_perm)
        message("  Significant at perm_fdr < 0.10: ", n_sig_perm_10)
        message("  Significant at perm_fdr < 0.20: ", n_sig_perm_20)
        
        # Show FDR distribution
        fdr_quantiles <- quantile(limma_results$perm_fdr, c(0, 0.01, 0.05, 0.1, 0.5), na.rm = TRUE)
        message("  FDR quantiles: min=", round(fdr_quantiles[1], 3),
                ", 1%=", round(fdr_quantiles[2], 3),
                ", 5%=", round(fdr_quantiles[3], 3),
                ", 10%=", round(fdr_quantiles[4], 3),
                ", median=", round(fdr_quantiles[5], 3))
        
        # Warning for small sample sizes
        n_per_group <- min(sum(groups == ref_condition), sum(groups == treat_condition))
        if (n_per_group <= 3) {
          message("\n  ⚠ WARNING: With only ", n_per_group, " samples per group, permutation FDR")
          message("    may be overly conservative. Consider using q-value instead.")
          message("    Unique permutations possible: ", choose(length(groups), n_per_group))
        }
      } else {
        message("  Warning: No valid permutations completed. Skipping perm_fdr.")
        limma_results$perm_fdr <- NA_real_
      }
    }
    
    # Merge with diff_results
    merge_cols <- c("feature_id", "limma_logFC", "AveExpr", "limma_t", 
                    "limma_pvalue", "limma_adj_pvalue_BH", "limma_adj_pvalue")
    if (n_permutations > 0) merge_cols <- c(merge_cols, "perm_fdr")
    merge_cols <- merge_cols[merge_cols %in% colnames(limma_results)]
    limma_to_merge <- limma_results[, merge_cols, drop = FALSE]
    
    diff_results <- merge(diff_results, limma_to_merge, by = "feature_id", all.x = TRUE)
  }
  
  # Add peptide count
  diff_results <- merge(diff_results, protein_peptide_counts, by = "protein_id", all.x = TRUE)
  
  # Save results
  write.csv(diff_results, file.path(output_dir, "differential_results.csv"), row.names = FALSE)
  write.csv(peak_quant, file.path(output_dir, "peak_quantification.csv"), row.names = FALSE)
  
  # Save BH-significant UniProt IDs as plain text list (deduplicated)
  bh_sig_ids <- unique(diff_results$protein_id[
    !is.na(diff_results$limma_adj_pvalue_BH) &
    diff_results$limma_adj_pvalue_BH < 0.05
  ])
  writeLines(bh_sig_ids, file.path(output_dir, "BH_significant_uniprot_ids.txt"))
  
  # Report significance with selected method
  n_sig <- sum(diff_results$limma_adj_pvalue < 0.05, na.rm = TRUE)
  n_sig_raw <- sum(diff_results$limma_pvalue < 0.05, na.rm = TRUE)
  n_sig_bh <- sum(diff_results$limma_adj_pvalue_BH < 0.05, na.rm = TRUE)
  
  message("\n  Significant features:")
  message("    Raw p < 0.05: ", n_sig_raw)
  message("    BH adj.P < 0.05: ", n_sig_bh)
  if (pvalue_adjustment == "qvalue") {
    message("    q-value < 0.05: ", n_sig, " (SELECTED)")
  } else if (pvalue_adjustment == "none") {
    message("    (Using raw p-values, no adjustment)")
  } else {
    message("    Selected (", toupper(pvalue_adjustment), "): ", n_sig)
  }
  if (n_permutations > 0 && "perm_fdr" %in% colnames(diff_results)) {
    n_sig_perm <- sum(diff_results$perm_fdr < 0.05, na.rm = TRUE)
    message("    Permutation FDR < 0.05: ", n_sig_perm)
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # STEP 7: Visualization and QC plots
  # ═══════════════════════════════════════════════════════════════════════════
  message("\nSTEP 7: Creating plots and QC diagnostics...")
  
  plot_dir <- file.path(output_dir, "plots")
  qc_dir <- file.path(output_dir, "qc")
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --- Volcano plots for all p-value adjustment methods ---
  message("  Creating volcano plots...")
  
  # Helper to report volcano stats
  report_volcano <- function(p, name, pval_col, data) {
    if (!is.null(p)) {
      n_sig <- attr(p, "n_significant")
      # Get y-axis range
      pvals <- data[[pval_col]]
      valid_pvals <- pvals[!is.na(pvals) & pvals > 0 & pvals <= 1]
      if (length(valid_pvals) > 0) {
        y_range <- range(-log10(valid_pvals))
        message("    ", name, ": sig=", n_sig, ", y-range=[", 
                round(y_range[1], 2), ", ", round(y_range[2], 2), "]")
      } else {
        message("    ", name, ": sig=", n_sig)
      }
    }
  }
  
  # 1. Raw p-value volcano
  tryCatch({
    p_raw <- create_volcano_plot(
      data = diff_results,
      pval_col = "limma_pvalue",
      pval_label = "raw P-value",
      title_suffix = "Raw P-values",
      treat_condition = treat_condition,
      ref_condition = ref_condition,
      agg_label = agg_label
    )
    if (!is.null(p_raw)) {
      ggsave(file.path(plot_dir, "volcano_raw.pdf"), p_raw, width = 14, height = 12)
      ggsave(file.path(plot_dir, "volcano_raw.png"), p_raw, width = 14, height = 12, dpi = 150)
      message("  Saved: volcano_raw.pdf/png")
      report_volcano(p_raw, "Raw", "limma_pvalue", diff_results)
    }
  }, error = function(e) message("  Warning: Could not create raw p-value volcano: ", e$message))
  
  # 2. BH-adjusted volcano (always available)
  tryCatch({
    p_bh <- create_volcano_plot(
      data = diff_results,
      pval_col = "limma_adj_pvalue_BH",
      pval_label = "BH adj. P-value",
      title_suffix = "Benjamini-Hochberg",
      treat_condition = treat_condition,
      ref_condition = ref_condition,
      agg_label = agg_label
    )
    if (!is.null(p_bh)) {
      ggsave(file.path(plot_dir, "volcano_BH.pdf"), p_bh, width = 14, height = 12)
      ggsave(file.path(plot_dir, "volcano_BH.png"), p_bh, width = 14, height = 12, dpi = 150)
      message("  Saved: volcano_BH.pdf/png")
      report_volcano(p_bh, "BH", "limma_adj_pvalue_BH", diff_results)
    }
  }, error = function(e) message("  Warning: Could not create BH volcano: ", e$message))
  
  # 3. Q-value volcano (if q-value was used)
  if (pvalue_adjustment == "qvalue" && "limma_adj_pvalue" %in% colnames(diff_results)) {
    tryCatch({
      p_qval <- create_volcano_plot(
        data = diff_results,
        pval_col = "limma_adj_pvalue",
        pval_label = "q-value",
        title_suffix = "Storey's q-value",
        treat_condition = treat_condition,
        ref_condition = ref_condition,
        agg_label = agg_label
      )
      if (!is.null(p_qval)) {
        ggsave(file.path(plot_dir, "volcano_qvalue.pdf"), p_qval, width = 14, height = 12)
        ggsave(file.path(plot_dir, "volcano_qvalue.png"), p_qval, width = 14, height = 12, dpi = 150)
        message("  Saved: volcano_qvalue.pdf/png")
        report_volcano(p_qval, "Q-value", "limma_adj_pvalue", diff_results)
      }
    }, error = function(e) message("  Warning: Could not create q-value volcano: ", e$message))
  }
  
  # 4. Permutation FDR volcano (if permutations were run)
  if (n_permutations > 0 && "perm_fdr" %in% colnames(diff_results)) {
    tryCatch({
      p_perm <- create_volcano_plot(
        data = diff_results,
        pval_col = "perm_fdr",
        pval_label = "Permutation FDR",
        title_suffix = paste0("Permutation FDR (n=", n_permutations, ")"),
        treat_condition = treat_condition,
        ref_condition = ref_condition,
        agg_label = agg_label
      )
      if (!is.null(p_perm)) {
        ggsave(file.path(plot_dir, "volcano_permFDR.pdf"), p_perm, width = 14, height = 12)
        ggsave(file.path(plot_dir, "volcano_permFDR.png"), p_perm, width = 14, height = 12, dpi = 150)
        message("  Saved: volcano_permFDR.pdf/png")
        report_volcano(p_perm, "PermFDR", "perm_fdr", diff_results)
      }
    }, error = function(e) message("  Warning: Could not create permutation FDR volcano: ", e$message))
  }
  
  # 5. Selected method volcano (main volcano.pdf for backwards compatibility)
  tryCatch({
    selected_col <- if (pvalue_adjustment == "qvalue") "limma_adj_pvalue" 
                    else if (pvalue_adjustment == "none") "limma_pvalue"
                    else "limma_adj_pvalue_BH"
    selected_label <- if (pvalue_adjustment == "qvalue") "q-value"
                      else if (pvalue_adjustment == "none") "raw P-value"
                      else "BH adj. P-value"
    
    p_selected <- create_volcano_plot(
      data = diff_results,
      pval_col = selected_col,
      pval_label = selected_label,
      title_suffix = paste0("Selected: ", toupper(pvalue_adjustment)),
      treat_condition = treat_condition,
      ref_condition = ref_condition,
      agg_label = agg_label
    )
    if (!is.null(p_selected)) {
      ggsave(file.path(plot_dir, "volcano.pdf"), p_selected, width = 14, height = 12)
      ggsave(file.path(plot_dir, "volcano.png"), p_selected, width = 14, height = 12, dpi = 150)
      message("  Saved: volcano.pdf/png (SELECTED method: ", toupper(pvalue_adjustment), ")")
    }
  }, error = function(e) message("  Warning: Could not create selected method volcano: ", e$message))
  
  # --- Peak status distribution ---
  tryCatch({
    status_counts <- as.data.frame(table(diff_results$peak_status))
    colnames(status_counts) <- c("Status", "Count")
    
    p_status <- ggplot(status_counts, aes(x = Status, y = Count, fill = Status)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = Count), vjust = -0.5, size = 5, fontface = "bold") +
      scale_fill_manual(values = c("matched" = "steelblue", 
                                   "control_only" = "orange",
                                   "treatment_only" = "forestgreen",
                                   "large_shift" = "purple")) +
      theme_minimal() +
      theme(legend.position = "none",
            axis.text = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 14, face = "bold"),
            plot.title = element_text(size = 16, face = "bold")) +
      labs(title = "Peak Status Distribution",
           subtitle = paste0(treat_condition, " vs ", ref_condition),
           x = "Peak Status", y = "Count")
    
    ggsave(file.path(qc_dir, "peak_status_distribution.pdf"), p_status, width = 10, height = 8)
    ggsave(file.path(qc_dir, "peak_status_distribution.png"), p_status, width = 10, height = 8, dpi = 150)
    message("  Saved: peak_status_distribution.pdf/png")
  }, error = function(e) message("  Warning: Could not create status plot: ", e$message))
  
  # --- Position shift histogram ---
  tryCatch({
    shift_data <- diff_results[!is.na(diff_results$delta_position), ]
    
    if (nrow(shift_data) > 0) {
      p_shift <- ggplot(shift_data, aes(x = delta_position)) +
        geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
        geom_vline(xintercept = c(-position_tolerance, position_tolerance), 
                   linetype = "dotted", color = "orange", linewidth = 0.8) +
        theme_minimal() +
        theme(axis.text = element_text(size = 14, face = "bold"),
              axis.title = element_text(size = 14, face = "bold"),
              plot.title = element_text(size = 16, face = "bold")) +
        labs(x = "Position Shift (fractions)", 
             y = "Count",
             title = "Peak Position Shift Distribution",
             subtitle = paste0("Orange lines: ±", position_tolerance, " fraction tolerance"))
      
      ggsave(file.path(qc_dir, "position_shift_histogram.pdf"), p_shift, width = 10, height = 8)
      ggsave(file.path(qc_dir, "position_shift_histogram.png"), p_shift, width = 10, height = 8, dpi = 150)
      message("  Saved: position_shift_histogram.pdf/png")
    }
  }, error = function(e) message("  Warning: Could not create shift histogram: ", e$message))
  
  # --- FC comparison: sum vs height ---
  tryCatch({
    fc_data <- diff_results[!is.na(diff_results$log2FC_sum) & !is.na(diff_results$log2FC_height), ]
    
    if (nrow(fc_data) > 0) {
      p_fc_compare <- ggplot(fc_data, aes(x = log2FC_sum, y = log2FC_height)) +
        geom_point(alpha = 0.3, size = 1.5) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
        geom_smooth(method = "lm", se = TRUE, color = "blue") +
        theme_minimal() +
        theme(axis.text = element_text(size = 14, face = "bold"),
              axis.title = element_text(size = 14, face = "bold"),
              plot.title = element_text(size = 16, face = "bold")) +
        labs(x = "log2FC (Local Sum)", 
             y = "log2FC (Peak Height)",
             title = "Fold Change Comparison: Sum vs Height",
             subtitle = "Red: y=x line, Blue: linear fit")
      
      ggsave(file.path(qc_dir, "fc_sum_vs_height.pdf"), p_fc_compare, width = 10, height = 10)
      ggsave(file.path(qc_dir, "fc_sum_vs_height.png"), p_fc_compare, width = 10, height = 10, dpi = 150)
      message("  Saved: fc_sum_vs_height.pdf/png")
    }
  }, error = function(e) message("  Warning: Could not create FC comparison plot: ", e$message))
  
  # --- Peaks per protein histogram ---
  tryCatch({
    peaks_per_prot <- as.data.frame(table(diff_results$protein_id))
    colnames(peaks_per_prot) <- c("protein_id", "n_peaks")
    
    p_peaks_per_prot <- ggplot(peaks_per_prot, aes(x = n_peaks)) +
      geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
      theme_minimal() +
      theme(axis.text = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 14, face = "bold"),
            plot.title = element_text(size = 16, face = "bold")) +
      labs(x = "Number of Peaks", 
           y = "Number of Proteins",
           title = "Peaks per Protein Distribution")
    
    ggsave(file.path(qc_dir, "peaks_per_protein.pdf"), p_peaks_per_prot, width = 10, height = 8)
    ggsave(file.path(qc_dir, "peaks_per_protein.png"), p_peaks_per_prot, width = 10, height = 8, dpi = 150)
    message("  Saved: peaks_per_protein.pdf/png")
  }, error = function(e) message("  Warning: Could not create peaks per protein plot: ", e$message))
  
  # --- P-value histograms ---
  tryCatch({
    pval_data <- diff_results[!is.na(diff_results$limma_pvalue) & 
                               diff_results$limma_pvalue > 0 & 
                               diff_results$limma_pvalue <= 1, ]
    
    if (nrow(pval_data) > 0) {
      # Raw p-value
      p_pval_raw <- ggplot(pval_data, aes(x = limma_pvalue)) +
        geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.7) +
        geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
        theme_minimal() +
        theme(axis.text = element_text(size = 14, face = "bold"),
              axis.title = element_text(size = 14, face = "bold"),
              plot.title = element_text(size = 16, face = "bold")) +
        labs(x = "P-value (raw)", y = "Count",
             title = "Raw P-value Distribution",
             subtitle = "GOOD: Enrichment at low p-values = real signal | FLAT: No signal")
      
      ggsave(file.path(qc_dir, "pvalue_raw_histogram.pdf"), p_pval_raw, width = 10, height = 8)
      ggsave(file.path(qc_dir, "pvalue_raw_histogram.png"), p_pval_raw, width = 10, height = 8, dpi = 150)
      
      # Adjusted p-value
      pval_adj_data <- diff_results[!is.na(diff_results$limma_adj_pvalue) & 
                                     diff_results$limma_adj_pvalue > 0 & 
                                     diff_results$limma_adj_pvalue <= 1, ]
      
      if (nrow(pval_adj_data) > 0) {
        p_pval_adj <- ggplot(pval_adj_data, aes(x = limma_adj_pvalue)) +
          geom_histogram(bins = 50, fill = "coral", color = "white", alpha = 0.7) +
          geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
          theme_minimal() +
          theme(axis.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 14, face = "bold"),
                plot.title = element_text(size = 16, face = "bold")) +
          labs(x = "Adjusted P-value (BH)", y = "Count",
               title = "Adjusted P-value Distribution",
               subtitle = "GOOD: Enrichment at low p-values = real signal surviving FDR correction")
        
        ggsave(file.path(qc_dir, "pvalue_adj_histogram.pdf"), p_pval_adj, width = 10, height = 8)
        ggsave(file.path(qc_dir, "pvalue_adj_histogram.png"), p_pval_adj, width = 10, height = 8, dpi = 150)
      }
      message("  Saved: pvalue_*_histogram.pdf/png")
    }
  }, error = function(e) message("  Warning: Could not create p-value histogram: ", e$message))
  
  # --- Replicate Correlations ---
  message("  Running replicate correlation analysis...")
  tryCatch({
    correlation_results <- list()
    
    for (cond in c(ref_condition, treat_condition)) {
      cond_data <- peak_quant[peak_quant$condition == cond & peak_quant$peak_detected, ]
      
      if (nrow(cond_data) < 10) next
      
      for (metric in c("local_sum", "max_height")) {
        if (!metric %in% colnames(cond_data)) next
        
        wide_data <- cond_data %>%
          select(feature_id, replicate, !!sym(metric)) %>%
          pivot_wider(names_from = replicate, values_from = !!sym(metric), values_fn = mean)
        
        rep_cols <- setdiff(colnames(wide_data), "feature_id")
        if (length(rep_cols) < 2) next
        
        cor_mat <- cor(wide_data[, rep_cols], use = "pairwise.complete.obs")
        mean_cor <- mean(cor_mat[lower.tri(cor_mat)], na.rm = TRUE)
        
        correlation_results[[paste0(metric, "_", cond)]] <- list(
          metric = metric, condition = cond, mean_correlation = mean_cor
        )
      }
    }
    
    if (length(correlation_results) > 0) {
      corr_summary_df <- data.frame(
        comparison = names(correlation_results),
        mean_correlation = sapply(correlation_results, function(x) x$mean_correlation),
        metric = sapply(correlation_results, function(x) x$metric),
        condition = sapply(correlation_results, function(x) x$condition),
        stringsAsFactors = FALSE
      )
      
      p_corr_summary <- ggplot(corr_summary_df, 
                                aes(x = metric, y = mean_correlation, fill = condition)) +
        geom_col(position = "dodge", alpha = 0.8) +
        geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
        geom_text(aes(label = round(mean_correlation, 2)), 
                  position = position_dodge(width = 0.9), vjust = -0.5, size = 4, fontface = "bold") +
        scale_fill_brewer(palette = "Set1") +
        theme_minimal() +
        theme(axis.text = element_text(size = 14, face = "bold"),
              axis.title = element_text(size = 14, face = "bold"),
              plot.title = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 12, face = "bold")) +
        ylim(0, 1.1) +
        labs(x = "Metric", y = "Mean Replicate Correlation",
             title = "Replicate Correlation Summary",
             subtitle = "Red line: R = 0.9 threshold | Higher is better")
      
      ggsave(file.path(qc_dir, "replicate_correlation_summary.pdf"), p_corr_summary, width = 10, height = 8)
      ggsave(file.path(qc_dir, "replicate_correlation_summary.png"), p_corr_summary, width = 10, height = 8, dpi = 150)
      write.csv(corr_summary_df, file.path(qc_dir, "replicate_correlations.csv"), row.names = FALSE)
      message("  Saved: replicate_correlation_summary.pdf/png")
    }
  }, error = function(e) message("  Warning: Could not create replicate correlation plots: ", e$message))
  
  # --- CV Analysis ---
  message("  Running CV analysis...")
  tryCatch({
    cv_data <- peak_quant %>%
      filter(peak_detected) %>%
      group_by(feature_id, protein_id, peak_id, condition) %>%
      summarise(
        mean_val = mean(local_sum, na.rm = TRUE),
        sd_val = sd(local_sum, na.rm = TRUE),
        cv = ifelse(mean_val > 0, sd_val / mean_val * 100, NA),
        n_reps = n(),
        .groups = "drop"
      ) %>%
      filter(!is.na(cv), is.finite(cv), cv >= 0)
    
    if (nrow(cv_data) > 0) {
      p_cv <- ggplot(cv_data, aes(x = cv, fill = condition)) +
        geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
        geom_vline(xintercept = 20, linetype = "dashed", color = "red") +
        scale_fill_brewer(palette = "Set1") +
        theme_minimal() +
        theme(axis.text = element_text(size = 14, face = "bold"),
              axis.title = element_text(size = 14, face = "bold"),
              plot.title = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 12, face = "bold")) +
        labs(x = "CV (%)", y = "Count",
             title = "CV Distribution: Local Sum",
             subtitle = "Red line: 20% CV threshold | Lower CV = more reproducible") +
        facet_wrap(~ condition, ncol = 1)
      
      ggsave(file.path(qc_dir, "cv_distribution.pdf"), p_cv, width = 10, height = 10)
      ggsave(file.path(qc_dir, "cv_distribution.png"), p_cv, width = 10, height = 10, dpi = 150)
      
      # CV summary
      cv_summary <- cv_data %>%
        group_by(condition) %>%
        summarise(
          median_cv = median(cv, na.rm = TRUE),
          mean_cv = mean(cv, na.rm = TRUE),
          pct_below_20 = 100 * mean(cv < 20, na.rm = TRUE),
          pct_below_30 = 100 * mean(cv < 30, na.rm = TRUE),
          n_features = n(),
          .groups = "drop"
        )
      write.csv(cv_summary, file.path(qc_dir, "cv_summary.csv"), row.names = FALSE)
      message("  Saved: cv_distribution.pdf/png, cv_summary.csv")
    }
  }, error = function(e) message("  Warning: Could not create CV analysis: ", e$message))
  
  # --- MA Plot ---
  message("  Creating MA plot...")
  tryCatch({
    ma_data <- diff_results %>%
      filter(!is.na(log2FC_sum), !is.na(ref_sum), !is.na(treat_sum)) %>%
      mutate(
        avg_expression = log2((ref_sum + treat_sum) / 2 + 1),
        significant = limma_adj_pvalue < 0.05 & abs(log2FC_sum) > 1
      )
    
    if (nrow(ma_data) > 0) {
      p_ma <- ggplot(ma_data, aes(x = avg_expression, y = log2FC_sum, color = significant)) +
        geom_point(alpha = 0.4, size = 1.5) +
        geom_hline(yintercept = 0, linetype = "solid", color = "gray40") +
        geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray60") +
        geom_smooth(aes(group = 1), method = "loess", color = "blue", se = FALSE, linewidth = 0.8) +
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray50"),
                           labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
        theme_minimal() +
        theme(axis.text = element_text(size = 14, face = "bold"),
              axis.title = element_text(size = 14, face = "bold"),
              plot.title = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 12, face = "bold")) +
        labs(x = "Average Expression (log2)", 
             y = "log2 Fold Change",
             title = "MA Plot",
             subtitle = "Blue line: LOESS trend | Should be centered at 0 if properly normalized",
             color = "")
      
      ggsave(file.path(qc_dir, "ma_plot.pdf"), p_ma, width = 12, height = 10)
      ggsave(file.path(qc_dir, "ma_plot.png"), p_ma, width = 12, height = 10, dpi = 150)
      message("  Saved: ma_plot.pdf/png")
    }
  }, error = function(e) message("  Warning: Could not create MA plot: ", e$message))
  
  # --- Adaptive window size distribution (if used) ---
  if (adaptive_window) {
    tryCatch({
      window_data <- peak_quant %>%
        filter(!is.na(used_window)) %>%
        distinct(protein_id, peak_id, .keep_all = TRUE)
      
      if (nrow(window_data) > 0) {
        p_window <- ggplot(window_data, aes(x = used_window)) +
          geom_histogram(binwidth = 1, fill = "steelblue", color = "white", alpha = 0.7) +
          geom_vline(xintercept = quant_window, linetype = "dashed", color = "red") +
          theme_minimal() +
          theme(axis.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 14, face = "bold"),
                plot.title = element_text(size = 16, face = "bold")) +
          labs(x = "Adaptive Window Size (fractions)", y = "Count",
               title = "Adaptive Window Size Distribution",
               subtitle = paste0("Red line: fixed quant_window (", quant_window, ") for comparison"))
        
        ggsave(file.path(qc_dir, "adaptive_window_distribution.pdf"), p_window, width = 10, height = 8)
        ggsave(file.path(qc_dir, "adaptive_window_distribution.png"), p_window, width = 10, height = 8, dpi = 150)
        message("  Saved: adaptive_window_distribution.pdf/png")
      }
    }, error = function(e) message("  Warning: Could not create window distribution plot: ", e$message))
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # STEP 8: Trace plots with individual peptides
  # ═══════════════════════════════════════════════════════════════════════════
  message("\nSTEP 8: Creating trace plots with peptide traces...")
  
  trace_plot_dir <- file.path(output_dir, "trace_plots")
  dir.create(trace_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  sig_proteins <- diff_results %>%
    filter(!is.na(limma_adj_pvalue)) %>%
    arrange(limma_adj_pvalue) %>%
    head(n_trace_plots) %>%
    pull(protein_id) %>%
    unique()
  
  if (!is.null(custom_ids)) {
    sig_proteins <- unique(c(sig_proteins, intersect(custom_ids, proteins)))
  }
  
  if (length(sig_proteins) > 0) {
    message("  Plotting ", length(sig_proteins), " proteins...")
    
    for (prot in sig_proteins) {
      tryCatch({
        prot_peptides <- protein_peptide_list[[prot]]
        prot_results <- diff_results %>% filter(protein_id == prot)
        
        for (pk_idx in seq_len(nrow(prot_results))) {
          peak_info <- prot_results[pk_idx, ]
          peak_center <- peak_info$ref_center
          
          peptide_data <- list()
          protein_data <- list()
          
          for (sample_id in sample_names) {
            parsed <- parse_sample_name(sample_id)
            if (!parsed$condition %in% c(ref_condition, treat_condition)) next
            
            trace_mat <- all_peptide_mats[[sample_id]]
            peps_in_sample <- prot_peptides[prot_peptides %in% rownames(trace_mat)]
            
            for (pep in peps_in_sample) {
              pep_trace <- trace_mat[pep, ]
              if (!all(is.na(pep_trace))) {
                peptide_data[[length(peptide_data) + 1]] <- data.frame(
                  fraction = as.numeric(fraction_names),
                  intensity = as.numeric(pep_trace),
                  peptide_id = pep,
                  sample_id = sample_id,
                  condition = parsed$condition,
                  replicate = parsed$replicate,
                  stringsAsFactors = FALSE
                )
              }
            }
            
            prot_mat <- protein_chromatograms[[sample_id]]
            if (prot %in% rownames(prot_mat)) {
              prot_trace <- prot_mat[prot, ]
              protein_data[[length(protein_data) + 1]] <- data.frame(
                fraction = as.numeric(fraction_names),
                intensity = as.numeric(prot_trace),
                sample_id = sample_id,
                condition = parsed$condition,
                replicate = parsed$replicate,
                stringsAsFactors = FALSE
              )
            }
          }
          
          if (length(protein_data) == 0) next
          
          peptide_df <- bind_rows(peptide_data)
          protein_df <- bind_rows(protein_data)
          
          protein_df <- protein_df %>%
            mutate(
              condition = factor(condition, levels = c(ref_condition, treat_condition)),
              group = paste(condition, replicate, sep = "_")
            )
          
          if (nrow(peptide_df) > 0) {
            peptide_df <- peptide_df %>%
              mutate(
                condition = factor(condition, levels = c(ref_condition, treat_condition)),
                peptide_group = paste(peptide_id, sample_id, sep = "_")
              )
          }
          
          title_main <- paste0(prot, " (peak ", peak_info$peak_id, ")")
          title_metrics <- sprintf("log2FC=%.2f, adj.P=%.2e, n=%d peptides",
                                   peak_info$limma_logFC %||% NA, 
                                   peak_info$limma_adj_pvalue %||% NA,
                                   peak_info$n_peptides %||% NA)
          
          n_ref <- length(unique(protein_df$replicate[protein_df$condition == ref_condition]))
          n_treat <- length(unique(protein_df$replicate[protein_df$condition == treat_condition]))
          
          ref_colors <- colorRampPalette(c("#A6CEE3", "#1F78B4", "#08519C"))(max(n_ref, 1))
          treat_colors <- colorRampPalette(c("#FCBBA1", "#FB6A4A", "#CB181D"))(max(n_treat, 1))
          
          all_groups <- unique(protein_df$group)
          color_map <- character()
          for (g in all_groups) {
            parts <- strsplit(g, "_")[[1]]
            cond <- paste(parts[1:(length(parts)-1)], collapse = "_")
            rep_num <- as.integer(parts[length(parts)])
            if (cond == ref_condition && rep_num <= length(ref_colors)) {
              color_map[g] <- ref_colors[rep_num]
            } else if (cond == treat_condition && rep_num <= length(treat_colors)) {
              color_map[g] <- treat_colors[rep_num]
            } else {
              color_map[g] <- "gray50"
            }
          }
          
          p <- ggplot()
          
          if (nrow(peptide_df) > 0) {
            # Add sample_group column to peptide_df for color matching
            peptide_df <- peptide_df %>%
              mutate(sample_group = paste(condition, replicate, sep = "_"))
            
            p <- p + geom_line(
              data = peptide_df,
              aes(x = fraction, y = intensity, group = peptide_group, color = sample_group),
              linewidth = 0.3,
              alpha = 0.25
            )
          }
          
          p <- p + 
            geom_line(
              data = protein_df,
              aes(x = fraction, y = intensity, color = group),
              linewidth = 1.2,
              alpha = 0.9
            ) +
            geom_point(
              data = protein_df,
              aes(x = fraction, y = intensity, color = group),
              size = 2,
              alpha = 0.9
            )
          
          if (!is.na(peak_center)) {
            p <- p + geom_vline(xintercept = peak_center, linetype = "dashed", 
                                color = "gray30", linewidth = 0.8)
          }
          
          p <- p +
            scale_color_manual(
              values = color_map,
              name = "Sample"
            ) +
            labs(
              x = "Fraction",
              y = "Intensity",
              title = title_main,
              subtitle = title_metrics
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = 18, face = "bold"),
              plot.subtitle = element_text(size = 16, face = "bold"),
              legend.position = "right",
              legend.text = element_text(size = 18, face = "bold"),
              legend.title = element_text(size = 18, face = "bold"),
              axis.text = element_text(size = 18, face = "bold"),
              axis.title = element_text(size = 18, face = "bold"),
              panel.grid.minor = element_blank()
            )
          
          safe_prot <- gsub("[^A-Za-z0-9_]", "_", prot)
          ggsave(file.path(trace_plot_dir, paste0("trace_", safe_prot, "_pk", peak_info$peak_id, ".pdf")), 
                 p, width = 14, height = 10)
          ggsave(file.path(trace_plot_dir, paste0("trace_", safe_prot, "_pk", peak_info$peak_id, ".png")), 
                 p, width = 14, height = 10, dpi = 150)
        }
        
      }, error = function(e) {
        message("    Warning: Could not plot ", prot, ": ", e$message)
      })
    }
    
    message("  Saved trace plots to: ", trace_plot_dir)
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Summary
  # ═══════════════════════════════════════════════════════════════════════════
  message("\n", paste(rep("═", 70), collapse = ""))
  message("PROTEIN-LEVEL PIPELINE COMPLETE")
  message(paste(rep("═", 70), collapse = ""))
  message("  Total proteins analyzed: ", length(proteins))
  message("  Proteins with peaks: ", length(reference_peaks))
  message("  Total features (protein_peak): ", n_distinct(diff_results$feature_id))
  message("  Significant (adj.P < 0.05): ", n_sig)
  message("  Output directory: ", output_dir)
  message(paste(rep("═", 70), collapse = ""))
  
  results <- list(
    protein_chromatograms = protein_chromatograms,
    reference_peaks = reference_peaks_df,
    peak_quantification = peak_quant,
    differential_results = diff_results,
    protein_map = protein_map,
    peptide_traces = all_peptide_mats,
    proteins = proteins,
    parameters = list(
      ref_condition = ref_condition,
      treat_condition = treat_condition,
      min_peptides = min_peptides,
      aggregation_method = aggregation_method,
      pvalue_adjustment = pvalue_adjustment,
      n_permutations = n_permutations,
      peak_params = peak_params,
      quant_window = quant_window,
      adaptive_window = adaptive_window,
      width_multiplier = width_multiplier,
      position_tolerance = position_tolerance,
      impute_value = impute_value
    )
  )
  
  return(results)
}


# ═══════════════════════════════════════════════════════════════════════════════
# USAGE EXAMPLE
# ═══════════════════════════════════════════════════════════════════════════════
#
# # Default: MaxLFQ + BH adjustment
# results <- run_protein_level_pipeline(
#   traces_list = pepTracesList_filtered,
#   ref_condition = "ctrl",
#   treat_condition = "ATP",
#   aggregation_method = "maxlfq",
#   pvalue_adjustment = "BH",
#   output_dir = "protein_results_maxlfq"
# )
#
# # With Storey's q-value (less conservative, good for CF-MS!)
# # Requires: BiocManager::install("qvalue")
# results <- run_protein_level_pipeline(
#   traces_list = pepTracesList_filtered,
#   ref_condition = "ctrl",
#   treat_condition = "ATP",
#   aggregation_method = "top10median",
#   pvalue_adjustment = "qvalue",
#   output_dir = "protein_results_qvalue"
# )
#
# # With permutation FDR (gold standard for correlated data)
# results <- run_protein_level_pipeline(
#   traces_list = pepTracesList_filtered,
#   ref_condition = "ctrl",
#   treat_condition = "ATP",
#   aggregation_method = "top10median",
#   pvalue_adjustment = "BH",
#   n_permutations = 100,  # 100-1000 recommended
#   output_dir = "protein_results_permFDR"
# )
#
# # No adjustment (exploratory, use with caution!)
# results <- run_protein_level_pipeline(
#   traces_list = pepTracesList_filtered,
#   ref_condition = "ctrl",
#   treat_condition = "ATP",
#   aggregation_method = "top10median",
#   pvalue_adjustment = "none",
#   output_dir = "protein_results_raw"
# )
#
# # Top10 median (robust, like Hi3 but with median)
# results <- run_protein_level_pipeline(
#   traces_list = pepTracesList_filtered,
#   ref_condition = "ctrl",
#   treat_condition = "ATP",
#   aggregation_method = "top10median",
#   min_peptides = 3,
#   output_dir = "protein_results_top10median"
# )
#
# # Top3 peptides (like iBAQ)
# results <- run_protein_level_pipeline(
#   traces_list = pepTracesList_filtered,
#   ref_condition = "ctrl",
#   treat_condition = "ATP",
#   aggregation_method = "top3",
#   min_peptides = 3,
#   output_dir = "protein_results_top3"
# )
#
# # Median (robust, fast)
# results <- run_protein_level_pipeline(
#   traces_list = pepTracesList_filtered,
#   ref_condition = "ctrl",
#   treat_condition = "ATP",
#   aggregation_method = "median",
#   min_peptides = 3,
#   output_dir = "protein_results_median"
# )
#
