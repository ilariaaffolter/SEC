# ═══════════════════════════════════════════════════════════════════════════════
# STANDALONE FUNCTION: Plot Peptide Traces for Proteoforms
# ═══════════════════════════════════════════════════════════════════════════════
#
# This function should be run AFTER run_p3_height_pipeline() completes.
# It plots individual peptide traces for proteoforms, with multiple ranking options.
#
# RANKING OPTIONS:
#   A) "correlation" - Best sibling peptide correlation (most coherent traces)
#   B) "pvalue"      - Most significant p-value (strongest differential signal)
#   C) "foldchange"  - Largest absolute fold change (biggest effect)
#   D) "custom"      - Plot specific proteoforms you specify
#
# Usage:
#   source("plot_proteoform_peptide_traces.R")
#   
#   # Option A: Best correlation
#   plot_proteoform_peptides(traces_list, results, ranking = "correlation", n_top = 10)
#   
#   # Option B: Most significant
#   plot_proteoform_peptides(traces_list, results, ranking = "pvalue", n_top = 10)
#   
#   # Option C: Largest fold change
#   plot_proteoform_peptides(traces_list, results, ranking = "foldchange", n_top = 10)
#   
#   # Option D: Custom selection
#   plot_proteoform_peptides(traces_list, results, ranking = "custom", 
#                            custom_ids = c("P0A6F5_A", "P06611_B"))
#
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
})

# ═══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTION: Calculate sibling peptide correlation for ranking option A
# ═══════════════════════════════════════════════════════════════════════════════

calculate_proteoform_clustering_quality <- function(traces_list, proteoform_mapping, design_matrix = NULL) {
  
  message("  Calculating sibling peptide correlations...")
  
  # Get sample names
  sample_names <- names(traces_list)
  
  # Create design matrix if not provided
  if (is.null(design_matrix)) {
    design_matrix <- data.frame(
      sample_id = sample_names,
      condition = sapply(strsplit(sample_names, "_"), function(x) paste(x[-length(x)], collapse = "_")),
      replicate = sapply(strsplit(sample_names, "_"), function(x) x[length(x)]),
      stringsAsFactors = FALSE
    )
  }
  
  # Get unique proteoforms
  proteoforms <- unique(proteoform_mapping$proteoform_id)
  
  quality_results <- list()
  total_zero_var_peptides <- 0
  
  pb <- txtProgressBar(min = 0, max = length(proteoforms), style = 3)
  
  for (i in seq_along(proteoforms)) {
    setTxtProgressBar(pb, i)
    pf <- proteoforms[i]
    
    # Get peptides for this proteoform
    pf_peptides <- proteoform_mapping$peptide_id[proteoform_mapping$proteoform_id == pf]
    
    if (length(pf_peptides) < 2) {
      next
    }
    
    # Collect traces across all samples
    all_traces <- list()
    
    for (sample_id in sample_names) {
      traces_obj <- traces_list[[sample_id]]
      trace_dt <- traces_obj$traces
      
      fraction_cols <- grep("^[0-9]+$", colnames(trace_dt), value = TRUE)
      fraction_cols <- fraction_cols[order(as.numeric(fraction_cols))]
      
      available_peptides <- intersect(pf_peptides, trace_dt$id)
      
      if (length(available_peptides) > 0) {
        if (is.data.table(trace_dt)) {
          pep_traces <- as.matrix(trace_dt[id %in% available_peptides, ..fraction_cols])
          rownames(pep_traces) <- trace_dt[id %in% available_peptides, id]
        } else {
          pep_traces <- as.matrix(trace_dt[trace_dt$id %in% available_peptides, fraction_cols])
          rownames(pep_traces) <- trace_dt$id[trace_dt$id %in% available_peptides]
        }
        all_traces[[sample_id]] <- pep_traces
      }
    }
    
    if (length(all_traces) == 0) next
    
    sample_correlations <- c()
    
    for (sample_id in names(all_traces)) {
      pep_mat <- all_traces[[sample_id]]
      
      if (nrow(pep_mat) >= 2) {
        row_vars <- apply(pep_mat, 1, var, na.rm = TRUE)
        valid_rows <- which(row_vars > 0 & !is.na(row_vars))
        total_zero_var_peptides <- total_zero_var_peptides + (nrow(pep_mat) - length(valid_rows))
        
        if (length(valid_rows) >= 2) {
          pep_mat_filtered <- pep_mat[valid_rows, , drop = FALSE]
          cor_mat <- cor(t(pep_mat_filtered), use = "pairwise.complete.obs")
          upper_tri <- cor_mat[upper.tri(cor_mat)]
          upper_tri <- upper_tri[!is.na(upper_tri) & is.finite(upper_tri)]
          
          if (length(upper_tri) > 0) {
            sample_correlations <- c(sample_correlations, upper_tri)
          }
        }
      }
    }
    
    if (length(sample_correlations) > 0) {
      quality_results[[pf]] <- data.frame(
        proteoform_id = pf,
        n_peptides = length(pf_peptides),
        mean_sib_corr = mean(sample_correlations, na.rm = TRUE),
        median_sib_corr = median(sample_correlations, na.rm = TRUE),
        min_sib_corr = min(sample_correlations, na.rm = TRUE),
        max_sib_corr = max(sample_correlations, na.rm = TRUE),
        sd_sib_corr = sd(sample_correlations, na.rm = TRUE),
        n_correlations = length(sample_correlations),
        stringsAsFactors = FALSE
      )
    }
  }
  
  close(pb)
  
  if (total_zero_var_peptides > 0) {
    message("\n  Filtered out ", total_zero_var_peptides, " peptide-sample combinations with zero variance")
  }
  
  quality_df <- bind_rows(quality_results)
  quality_df <- quality_df %>% arrange(desc(mean_sib_corr))
  
  message("  Calculated quality for ", nrow(quality_df), " proteoforms")
  
  return(quality_df)
}


# ═══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTION: Plot peptide traces for a single proteoform
# ═══════════════════════════════════════════════════════════════════════════════

plot_single_proteoform_peptides <- function(proteoform_id, traces_list, proteoform_mapping,
                                             design_matrix = NULL, highlight_fractions = NULL,
                                             show_mean = TRUE, alpha_peptides = 0.4,
                                             title_suffix = "") {
  
  sample_names <- names(traces_list)
  
  if (is.null(design_matrix)) {
    design_matrix <- data.frame(
      sample_id = sample_names,
      condition = sapply(strsplit(sample_names, "_"), function(x) paste(x[-length(x)], collapse = "_")),
      replicate = sapply(strsplit(sample_names, "_"), function(x) x[length(x)]),
      stringsAsFactors = FALSE
    )
  }
  
  pf_peptides <- proteoform_mapping$peptide_id[proteoform_mapping$proteoform_id == proteoform_id]
  
  if (length(pf_peptides) == 0) {
    warning("No peptides found for proteoform: ", proteoform_id)
    return(NULL)
  }
  
  all_data <- list()
  
  for (sample_id in sample_names) {
    traces_obj <- traces_list[[sample_id]]
    trace_dt <- traces_obj$traces
    
    fraction_cols <- grep("^[0-9]+$", colnames(trace_dt), value = TRUE)
    fraction_cols <- fraction_cols[order(as.numeric(fraction_cols))]
    
    available_peptides <- intersect(pf_peptides, trace_dt$id)
    
    if (length(available_peptides) > 0) {
      for (pep_id in available_peptides) {
        if (is.data.table(trace_dt)) {
          pep_trace <- as.numeric(trace_dt[id == pep_id, ..fraction_cols])
        } else {
          pep_trace <- as.numeric(trace_dt[trace_dt$id == pep_id, fraction_cols])
        }
        
        all_data[[length(all_data) + 1]] <- data.frame(
          sample_id = sample_id,
          peptide_id = pep_id,
          fraction = as.numeric(fraction_cols),
          intensity = pep_trace,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(all_data) == 0) {
    warning("No trace data found for proteoform: ", proteoform_id)
    return(NULL)
  }
  
  plot_data <- bind_rows(all_data)
  plot_data <- plot_data %>% left_join(design_matrix, by = "sample_id")
  plot_data$trace_id <- paste(plot_data$sample_id, plot_data$peptide_id, sep = "___")
  
  mean_traces <- plot_data %>%
    group_by(condition, fraction) %>%
    summarise(
      mean_intensity = mean(intensity, na.rm = TRUE),
      sd_intensity = sd(intensity, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  n_peptides <- length(unique(plot_data$peptide_id))
  n_samples <- length(unique(plot_data$sample_id))
  conditions <- unique(plot_data$condition)
  
  condition_colors <- setNames(
    scales::hue_pal()(length(conditions)),
    conditions
  )
  
  p <- ggplot() +
    geom_line(
      data = plot_data,
      aes(x = fraction, y = intensity, group = trace_id, color = condition),
      alpha = alpha_peptides,
      linewidth = 0.3
    )
  
  if (show_mean) {
    p <- p +
      geom_line(
        data = mean_traces,
        aes(x = fraction, y = mean_intensity, color = condition),
        linewidth = 1.5,
        alpha = 0.9
      ) +
      geom_point(
        data = mean_traces,
        aes(x = fraction, y = mean_intensity, color = condition),
        size = 2,
        alpha = 0.8
      )
  }
  
  if (!is.null(highlight_fractions)) {
    p <- p +
      geom_vline(
        xintercept = highlight_fractions,
        linetype = "dashed",
        color = "gray40",
        linewidth = 0.5
      )
  }
  
  p <- p +
    scale_color_manual(values = condition_colors) +
    labs(
      x = "Fraction",
      y = "Intensity",
      title = paste0("Peptide Traces: ", proteoform_id),
      subtitle = paste0(n_peptides, " peptides × ", n_samples, " samples | ",
                        "Thin lines = individual peptides | Thick lines = condition mean",
                        title_suffix),
      color = "Condition"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  return(p)
}


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN FUNCTION: Plot peptide traces with multiple ranking options
# ═══════════════════════════════════════════════════════════════════════════════

#' Plot peptide traces for proteoforms with multiple ranking options
#' 
#' @param traces_list The pepTracesList_filtered used in the pipeline
#' @param pipeline_results Results from run_p3_height_pipeline()
#' @param ranking Ranking method: "correlation", "pvalue", "foldchange", or "custom"
#' @param n_top Number of top proteoforms to plot (ignored if ranking = "custom")
#' @param custom_ids Vector of proteoform_ids to plot (required if ranking = "custom")
#' @param output_dir Directory to save plots
#' @param file_format "pdf", "png", or "both"
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return Data frame with ranking information for plotted proteoforms
#' @export
plot_proteoform_peptides <- function(traces_list,
                                      pipeline_results,
                                      ranking = c("correlation", "pvalue", "foldchange", "custom"),
                                      n_top = 10,
                                      custom_ids = NULL,
                                      output_dir = "peptide_trace_plots",
                                      file_format = "both",
                                      width = 12,
                                      height = 6) {
  
  ranking <- match.arg(ranking)
  
  message("\n", paste(rep("=", 70), collapse = ""))
  message("PLOTTING PEPTIDE TRACES FOR PROTEOFORMS")
  message("Ranking method: ", toupper(ranking))
  message(paste(rep("=", 70), collapse = ""))
  
  # Validate inputs
  if (ranking == "custom" && is.null(custom_ids)) {
    stop("custom_ids must be provided when ranking = 'custom'")
  }
  
  # Create output directory
  subdir <- switch(ranking,
                   "correlation" = "by_correlation",
                   "pvalue" = "by_pvalue",
                   "foldchange" = "by_foldchange",
                   "custom" = "custom_selection")
  
  full_output_dir <- file.path(output_dir, subdir)
  if (!dir.exists(full_output_dir)) {
    dir.create(full_output_dir, recursive = TRUE)
    message("Created output directory: ", full_output_dir)
  }
  
  # Get proteoform mapping from pipeline results
  if (!"proteoform_mapping" %in% names(pipeline_results)) {
    stop("proteoform_mapping not found in pipeline_results")
  }
  proteoform_mapping <- pipeline_results$proteoform_mapping
  
  # Get differential results for ranking options B and C
  if (!"differential_results" %in% names(pipeline_results)) {
    stop("differential_results not found in pipeline_results")
  }
  diff_results <- pipeline_results$differential_results
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SELECT PROTEOFORMS BASED ON RANKING METHOD
  # ═══════════════════════════════════════════════════════════════════════════
  
  if (ranking == "correlation") {
    # Option A: Best sibling peptide correlation
    message("\nOption A: Ranking by sibling peptide correlation (best clustering quality)")
    
    quality_df <- calculate_proteoform_clustering_quality(
      traces_list = traces_list,
      proteoform_mapping = proteoform_mapping
    )
    
    # Save quality metrics
    write.csv(quality_df, file.path(full_output_dir, "ranking_by_correlation.csv"), row.names = FALSE)
    
    proteoforms_to_plot <- head(quality_df$proteoform_id, n_top)
    ranking_info <- quality_df %>%
      filter(proteoform_id %in% proteoforms_to_plot) %>%
      mutate(rank = row_number())
    
    message("\n  Top ", n_top, " proteoforms by correlation:")
    print(head(ranking_info[, c("rank", "proteoform_id", "n_peptides", "mean_sib_corr")], n_top))
    
  } else if (ranking == "pvalue") {
    # Option B: Most significant p-value
    message("\nOption B: Ranking by most significant adjusted p-value")
    
    # Get best (minimum) p-value per proteoform
    pvalue_ranking <- diff_results %>%
      filter(!is.na(limma_adj_pvalue), limma_adj_pvalue > 0) %>%
      group_by(proteoform_id) %>%
      summarise(
        min_adj_pvalue = min(limma_adj_pvalue, na.rm = TRUE),
        best_log2FC = log2FC_sum[which.min(limma_adj_pvalue)],
        best_peak_id = peak_id[which.min(limma_adj_pvalue)],
        n_peaks = n(),
        .groups = "drop"
      ) %>%
      arrange(min_adj_pvalue) %>%
      mutate(rank = row_number())
    
    # Get peptide counts
    peptide_counts <- proteoform_mapping %>%
      group_by(proteoform_id) %>%
      summarise(n_peptides = n_distinct(peptide_id), .groups = "drop")
    
    pvalue_ranking <- pvalue_ranking %>%
      left_join(peptide_counts, by = "proteoform_id")
    
    # Save ranking
    write.csv(pvalue_ranking, file.path(full_output_dir, "ranking_by_pvalue.csv"), row.names = FALSE)
    
    proteoforms_to_plot <- head(pvalue_ranking$proteoform_id, n_top)
    ranking_info <- pvalue_ranking %>% filter(proteoform_id %in% proteoforms_to_plot)
    
    message("\n  Top ", n_top, " proteoforms by p-value:")
    print(head(ranking_info[, c("rank", "proteoform_id", "n_peptides", "min_adj_pvalue", "best_log2FC")], n_top))
    
  } else if (ranking == "foldchange") {
    # Option C: Largest absolute fold change
    message("\nOption C: Ranking by largest absolute fold change")
    
    # Get maximum absolute fold change per proteoform
    fc_ranking <- diff_results %>%
      filter(!is.na(log2FC_sum), is.finite(log2FC_sum)) %>%
      group_by(proteoform_id) %>%
      summarise(
        max_abs_log2FC = max(abs(log2FC_sum), na.rm = TRUE),
        best_log2FC = log2FC_sum[which.max(abs(log2FC_sum))],
        best_pvalue = limma_adj_pvalue[which.max(abs(log2FC_sum))],
        best_peak_id = peak_id[which.max(abs(log2FC_sum))],
        n_peaks = n(),
        .groups = "drop"
      ) %>%
      arrange(desc(max_abs_log2FC)) %>%
      mutate(rank = row_number())
    
    # Get peptide counts
    peptide_counts <- proteoform_mapping %>%
      group_by(proteoform_id) %>%
      summarise(n_peptides = n_distinct(peptide_id), .groups = "drop")
    
    fc_ranking <- fc_ranking %>%
      left_join(peptide_counts, by = "proteoform_id")
    
    # Save ranking
    write.csv(fc_ranking, file.path(full_output_dir, "ranking_by_foldchange.csv"), row.names = FALSE)
    
    proteoforms_to_plot <- head(fc_ranking$proteoform_id, n_top)
    ranking_info <- fc_ranking %>% filter(proteoform_id %in% proteoforms_to_plot)
    
    message("\n  Top ", n_top, " proteoforms by fold change:")
    print(head(ranking_info[, c("rank", "proteoform_id", "n_peptides", "max_abs_log2FC", "best_log2FC", "best_pvalue")], n_top))
    
  } else if (ranking == "custom") {
    # Option D: User-specified proteoforms
    message("\nOption D: Plotting user-specified proteoforms")
    
    # Validate custom_ids exist
    available_pfs <- unique(proteoform_mapping$proteoform_id)
    valid_ids <- custom_ids[custom_ids %in% available_pfs]
    invalid_ids <- custom_ids[!custom_ids %in% available_pfs]
    
    if (length(invalid_ids) > 0) {
      message("  Warning: These proteoform_ids were not found: ", paste(invalid_ids, collapse = ", "))
    }
    
    if (length(valid_ids) == 0) {
      stop("None of the specified proteoform_ids were found in the data")
    }
    
    proteoforms_to_plot <- valid_ids
    
    # Get info for these proteoforms
    peptide_counts <- proteoform_mapping %>%
      filter(proteoform_id %in% valid_ids) %>%
      group_by(proteoform_id) %>%
      summarise(n_peptides = n_distinct(peptide_id), .groups = "drop")
    
    diff_info <- diff_results %>%
      filter(proteoform_id %in% valid_ids) %>%
      group_by(proteoform_id) %>%
      summarise(
        min_adj_pvalue = min(limma_adj_pvalue, na.rm = TRUE),
        max_abs_log2FC = max(abs(log2FC_sum), na.rm = TRUE),
        n_peaks = n(),
        .groups = "drop"
      )
    
    ranking_info <- peptide_counts %>%
      left_join(diff_info, by = "proteoform_id") %>%
      mutate(rank = row_number())
    
    # Save info
    write.csv(ranking_info, file.path(full_output_dir, "custom_selection_info.csv"), row.names = FALSE)
    
    message("\n  Proteoforms to plot: ", length(valid_ids))
    print(ranking_info[, c("rank", "proteoform_id", "n_peptides", "min_adj_pvalue", "max_abs_log2FC")])
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # GET PEAK POSITIONS FOR HIGHLIGHTING
  # ═══════════════════════════════════════════════════════════════════════════
  
  reference_peaks <- NULL
  if ("reference_peaks" %in% names(pipeline_results)) {
    reference_peaks <- pipeline_results$reference_peaks
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # CREATE PLOTS
  # ═══════════════════════════════════════════════════════════════════════════
  
  message("\nGenerating peptide trace plots...")
  
  for (i in seq_along(proteoforms_to_plot)) {
    pf <- proteoforms_to_plot[i]
    message("  [", i, "/", length(proteoforms_to_plot), "] Plotting: ", pf)
    
    # Get peak positions for highlighting
    highlight_fracs <- NULL
    if (!is.null(reference_peaks) && is.data.frame(reference_peaks)) {
      pf_peaks <- reference_peaks[reference_peaks$proteoform_id == pf, ]
      if (nrow(pf_peaks) > 0 && "center" %in% colnames(pf_peaks)) {
        highlight_fracs <- pf_peaks$center
      }
    }
    
    # Get ranking info for subtitle
    pf_info <- ranking_info[ranking_info$proteoform_id == pf, ]
    
    # Build subtitle suffix based on ranking method
    title_suffix <- ""
    if (nrow(pf_info) > 0) {
      if (ranking == "correlation") {
        title_suffix <- sprintf("\nRank #%d | Mean sib. corr: %.3f", 
                                pf_info$rank[1], pf_info$mean_sib_corr[1])
      } else if (ranking == "pvalue") {
        title_suffix <- sprintf("\nRank #%d | adj.P: %.2e | log2FC: %.2f", 
                                pf_info$rank[1], pf_info$min_adj_pvalue[1], pf_info$best_log2FC[1])
      } else if (ranking == "foldchange") {
        title_suffix <- sprintf("\nRank #%d | log2FC: %.2f | adj.P: %.2e", 
                                pf_info$rank[1], pf_info$best_log2FC[1], pf_info$best_pvalue[1])
      } else if (ranking == "custom") {
        if (!is.na(pf_info$min_adj_pvalue[1]) && !is.na(pf_info$max_abs_log2FC[1])) {
          title_suffix <- sprintf("\nadj.P: %.2e | max|log2FC|: %.2f", 
                                  pf_info$min_adj_pvalue[1], pf_info$max_abs_log2FC[1])
        }
      }
    }
    
    # Create plot
    p <- plot_single_proteoform_peptides(
      proteoform_id = pf,
      traces_list = traces_list,
      proteoform_mapping = proteoform_mapping,
      highlight_fractions = highlight_fracs,
      title_suffix = title_suffix
    )
    
    if (!is.null(p)) {
      safe_pf <- gsub("[^A-Za-z0-9_]", "_", pf)
      
      # Add rank prefix to filename for sorting
      if (ranking != "custom") {
        rank_num <- pf_info$rank[1]
        filename_prefix <- sprintf("rank%02d_%s", rank_num, safe_pf)
      } else {
        filename_prefix <- safe_pf
      }
      
      if (file_format %in% c("pdf", "both")) {
        ggsave(file.path(full_output_dir, paste0(filename_prefix, ".pdf")),
               p, width = width, height = height)
      }
      if (file_format %in% c("png", "both")) {
        ggsave(file.path(full_output_dir, paste0(filename_prefix, ".png")),
               p, width = width, height = height, dpi = 150)
      }
    }
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SUMMARY
  # ═══════════════════════════════════════════════════════════════════════════
  
  message("\n", paste(rep("=", 70), collapse = ""))
  message("COMPLETE")
  message(paste(rep("=", 70), collapse = ""))
  message("  Ranking method: ", ranking)
  message("  Plots created: ", length(proteoforms_to_plot))
  message("  Output directory: ", full_output_dir)
  
  return(ranking_info)
}


# ═══════════════════════════════════════════════════════════════════════════════
# CONVENIENCE WRAPPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

#' Plot top proteoforms by best clustering quality (sibling peptide correlation)
#' @export
plot_by_correlation <- function(traces_list, pipeline_results, n_top = 10, 
                                 output_dir = "peptide_trace_plots", ...) {
  plot_proteoform_peptides(traces_list, pipeline_results, 
                           ranking = "correlation", n_top = n_top, 
                           output_dir = output_dir, ...)
}

#' Plot top proteoforms by most significant p-value
#' @export
plot_by_pvalue <- function(traces_list, pipeline_results, n_top = 10, 
                            output_dir = "peptide_trace_plots", ...) {
  plot_proteoform_peptides(traces_list, pipeline_results, 
                           ranking = "pvalue", n_top = n_top, 
                           output_dir = output_dir, ...)
}

#' Plot top proteoforms by largest fold change
#' @export
plot_by_foldchange <- function(traces_list, pipeline_results, n_top = 10, 
                                output_dir = "peptide_trace_plots", ...) {
  plot_proteoform_peptides(traces_list, pipeline_results, 
                           ranking = "foldchange", n_top = n_top, 
                           output_dir = output_dir, ...)
}

#' Plot specific proteoforms by ID
#' @export
plot_custom <- function(traces_list, pipeline_results, custom_ids, 
                         output_dir = "peptide_trace_plots", ...) {
  plot_proteoform_peptides(traces_list, pipeline_results, 
                           ranking = "custom", custom_ids = custom_ids, 
                           output_dir = output_dir, ...)
}


# ═══════════════════════════════════════════════════════════════════════════════
# USAGE EXAMPLES
# ═══════════════════════════════════════════════════════════════════════════════
#
# # After running the main pipeline:
# results_height <- run_p3_height_pipeline(...)
#
# # Source this file:
# source("plot_proteoform_peptide_traces.R")
#
# # Option A: Best clustering quality (sibling peptide correlation)
# corr_ranking <- plot_proteoform_peptides(
#   traces_list = pepTracesList_filtered,
#   pipeline_results = results_height,
#   ranking = "correlation",
#   n_top = 10,
#   output_dir = "peptide_trace_plots"
# )
# # Or use convenience wrapper:
# corr_ranking <- plot_by_correlation(pepTracesList_filtered, results_height, n_top = 10)
#
# # Option B: Most significant p-value
# pval_ranking <- plot_proteoform_peptides(
#   traces_list = pepTracesList_filtered,
#   pipeline_results = results_height,
#   ranking = "pvalue",
#   n_top = 10,
#   output_dir = "peptide_trace_plots"
# )
# # Or use convenience wrapper:
# pval_ranking <- plot_by_pvalue(pepTracesList_filtered, results_height, n_top = 10)
#
# # Option C: Largest fold change
# fc_ranking <- plot_proteoform_peptides(
#   traces_list = pepTracesList_filtered,
#   pipeline_results = results_height,
#   ranking = "foldchange",
#   n_top = 10,
#   output_dir = "peptide_trace_plots"
# )
# # Or use convenience wrapper:
# fc_ranking <- plot_by_foldchange(pepTracesList_filtered, results_height, n_top = 10)
#
# # Option D: Custom selection
# custom_info <- plot_proteoform_peptides(
#   traces_list = pepTracesList_filtered,
#   pipeline_results = results_height,
#   ranking = "custom",
#   custom_ids = c("P0A6F5_A", "P06611_B", "P37028_A"),
#   output_dir = "peptide_trace_plots"
# )
# # Or use convenience wrapper:
# custom_info <- plot_custom(pepTracesList_filtered, results_height, 
#                            custom_ids = c("P0A6F5_A", "P06611_B"))
#
# ═══════════════════════════════════════════════════════════════════════════════
