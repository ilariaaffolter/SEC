# ═══════════════════════════════════════════════════════════════════════════════
# PROTEIN TRACE PLOTTING SCRIPT
# ═══════════════════════════════════════════════════════════════════════════════
# 
# This script plots protein chromatograms with individual peptide traces.
# Can be used standalone after running the protein-level pipeline.
#
# Usage:
#   source("plot_protein_traces.R")
#   
#   # Plot specific proteins
#   plot_protein_traces(results, protein_ids = c("P12345", "Q67890"))
#   
#   # Plot significant proteins
#   plot_significant_proteins(results, max_plots = 30)
#   
#   # Plot high fold-change proteins
#   plot_high_fc_proteins(results, fc_threshold = 2, max_plots = 30)
#
# ═══════════════════════════════════════════════════════════════════════════════

library(ggplot2)
library(dplyr)
library(tidyr)

# ═══════════════════════════════════════════════════════════════════════════════
# MAIN PLOTTING FUNCTION
# ═══════════════════════════════════════════════════════════════════════════════

#' Plot protein trace with peptide traces
#' 
#' @param protein_id Protein ID to plot
#' @param protein_chromatograms List of protein chromatogram matrices (from results)
#' @param peptide_traces List of peptide trace matrices (from results)
#' @param protein_map Data frame with peptide_id and protein_id mapping
#' @param diff_results Differential results data frame (optional, for annotations)
#' @param ref_condition Reference condition name
#' @param treat_condition Treatment condition name
#' @param show_peptides Whether to show individual peptide traces (default TRUE)
#' @param peptide_alpha Alpha for peptide traces (default 0.3)
#' @param protein_linewidth Line width for protein trace (default 1.5)
#' @param output_file Optional output file path (PDF or PNG)
#' @return ggplot object
plot_single_protein_trace <- function(
    protein_id,
    protein_chromatograms,
    peptide_traces,
    protein_map,
    diff_results = NULL,
    ref_condition = "ctrl",
    treat_condition = "ATP",
    show_peptides = TRUE,
    peptide_alpha = 0.3,
    protein_linewidth = 1.5,
    output_file = NULL
) {
  
 # Get peptides for this protein
  protein_peptides <- protein_map$peptide_id[protein_map$protein_id == protein_id]
  
  # Parse sample names to get condition info
  sample_names <- names(protein_chromatograms)
  sample_info <- data.frame(
    sample_id = sample_names,
    condition = sapply(strsplit(sample_names, "_"), `[`, 1),
    replicate = sapply(strsplit(sample_names, "_"), function(x) paste(x[-1], collapse = "_")),
    stringsAsFactors = FALSE
  )
  
  n_fractions <- ncol(protein_chromatograms[[1]])
  fraction_names <- colnames(protein_chromatograms[[1]])
  
  # Build protein trace data
  protein_data <- data.frame()
  for (sample_id in sample_names) {
    if (protein_id %in% rownames(protein_chromatograms[[sample_id]])) {
      chrom <- protein_chromatograms[[sample_id]][protein_id, ]
      df <- data.frame(
        fraction = 1:n_fractions,
        intensity = as.numeric(chrom),
        sample_id = sample_id,
        condition = sample_info$condition[sample_info$sample_id == sample_id],
        replicate = sample_info$replicate[sample_info$sample_id == sample_id],
        type = "protein",
        stringsAsFactors = FALSE
      )
      protein_data <- rbind(protein_data, df)
    }
  }
  
  # Build peptide trace data
  peptide_data <- data.frame()
  if (show_peptides && length(protein_peptides) > 0) {
    for (sample_id in sample_names) {
      pep_mat <- peptide_traces[[sample_id]]
      peps_in_sample <- protein_peptides[protein_peptides %in% rownames(pep_mat)]
      
      for (pep in peps_in_sample) {
        chrom <- pep_mat[pep, ]
        df <- data.frame(
          fraction = 1:n_fractions,
          intensity = as.numeric(chrom),
          sample_id = sample_id,
          condition = sample_info$condition[sample_info$sample_id == sample_id],
          replicate = sample_info$replicate[sample_info$sample_id == sample_id],
          peptide_id = pep,
          type = "peptide",
          stringsAsFactors = FALSE
        )
        peptide_data <- rbind(peptide_data, df)
      }
    }
  }
  
  # Get stats from diff_results if available
  subtitle_text <- ""
  peak_positions <- NULL
  if (!is.null(diff_results)) {
    prot_results <- diff_results[diff_results$protein_id == protein_id, ]
    if (nrow(prot_results) > 0) {
      # Get first peak's stats for subtitle
      first_peak <- prot_results[1, ]
      fc_text <- ifelse(!is.na(first_peak$limma_logFC), 
                        sprintf("log2FC=%.2f", first_peak$limma_logFC), "")
      pval_text <- ifelse(!is.na(first_peak$limma_adj_pvalue),
                          sprintf("adj.P=%.2e", first_peak$limma_adj_pvalue), "")
      n_pep_text <- ifelse(!is.na(first_peak$n_peptides),
                           sprintf("n=%d peptides", first_peak$n_peptides), "")
      subtitle_text <- paste(c(fc_text, pval_text, n_pep_text)[c(fc_text, pval_text, n_pep_text) != ""], 
                             collapse = ", ")
      
      # Get peak positions for vertical lines
      peak_positions <- prot_results$ref_center[!is.na(prot_results$ref_center)]
    }
  }
  
  # Define color palettes: shades of blue for ctrl, shades of coral/red for treatment
  # Get unique replicates per condition
  ctrl_samples <- sample_info$sample_id[sample_info$condition == ref_condition]
  treat_samples <- sample_info$sample_id[sample_info$condition == treat_condition]
  
  # Create color palette with different shades
  n_ctrl <- length(ctrl_samples)
  n_treat <- length(treat_samples)
  
  # Blue shades for control (light to dark)
  ctrl_colors <- colorRampPalette(c("#A6CEE3", "#1F78B4", "#08519C"))(max(3, n_ctrl))[1:n_ctrl]
  names(ctrl_colors) <- ctrl_samples
  
  # Red/coral shades for treatment (light to dark)
  treat_colors <- colorRampPalette(c("#FCBBA1", "#FB6A4A", "#CB181D"))(max(3, n_treat))[1:n_treat]
  names(treat_colors) <- treat_samples
  
  # Combine colors
  all_colors <- c(ctrl_colors, treat_colors)
  
  # Create plot - single panel with all samples
  p <- ggplot()
  
  # Add peptide traces first (background) - use same colors but more transparent
  if (show_peptides && nrow(peptide_data) > 0) {
    p <- p + geom_line(
      data = peptide_data,
      aes(x = fraction, y = intensity, group = interaction(peptide_id, sample_id), 
          color = sample_id),
      alpha = peptide_alpha, linewidth = 0.5
    )
  }
  
  # Add protein traces (foreground) - lines + points
  if (nrow(protein_data) > 0) {
    p <- p + 
      geom_line(
        data = protein_data,
        aes(x = fraction, y = intensity, group = sample_id, color = sample_id),
        linewidth = protein_linewidth, alpha = 0.9
      ) +
      geom_point(
        data = protein_data,
        aes(x = fraction, y = intensity, color = sample_id),
        size = 2, alpha = 0.9
      )
  }
  
  # Add peak position lines
  if (!is.null(peak_positions) && length(peak_positions) > 0) {
    for (pk in peak_positions) {
      p <- p + geom_vline(xintercept = pk, linetype = "dashed", color = "gray40", linewidth = 0.5)
    }
  }
  
  # Create legend labels that group by condition
  legend_labels <- sample_names
  names(legend_labels) <- sample_names
  
  # Styling - NO faceting, single panel
  p <- p +
    scale_color_manual(
      values = all_colors,
      labels = legend_labels,
      name = "Sample"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      axis.text = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 18, face = "bold"),
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14),
      legend.text = element_text(size = 18, face = "bold"),
      legend.title = element_text(size = 18, face = "bold")
    ) +
    labs(
      x = "Fraction",
      y = "Intensity",
      title = paste0(protein_id, " (peak ", 
                     ifelse(!is.null(diff_results) && nrow(diff_results[diff_results$protein_id == protein_id, ]) > 0,
                            diff_results[diff_results$protein_id == protein_id, ]$peak_id[1], "1"), ")"),
      subtitle = subtitle_text
    )
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 12, height = 8, dpi = 150)
    message("  Saved: ", output_file)
  }
  
  return(p)
}


# ═══════════════════════════════════════════════════════════════════════════════
# BATCH PLOTTING FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

#' Plot traces for a list of protein IDs
#' 
#' @param results Results object from run_protein_level_pipeline()
#' @param protein_ids Vector of protein IDs to plot
#' @param output_dir Output directory for plots (default: "custom_trace_plots")
#' @param show_peptides Whether to show peptide traces
#' @param format Output format: "pdf", "png", or "both"
#' @return Invisibly returns list of ggplot objects
plot_protein_traces <- function(
    results,
    protein_ids,
    output_dir = "custom_trace_plots",
    show_peptides = TRUE,
    format = "both"
) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Extract components from results
  protein_chromatograms <- results$protein_chromatograms
  peptide_traces <- results$peptide_traces
  protein_map <- results$protein_map
  diff_results <- results$differential_results
  ref_condition <- results$parameters$ref_condition
  treat_condition <- results$parameters$treat_condition
  
  # Filter to valid proteins
  valid_proteins <- protein_ids[protein_ids %in% results$proteins]
  if (length(valid_proteins) == 0) {
    message("No valid protein IDs found in results.")
    return(invisible(list()))
  }
  
  message("Plotting ", length(valid_proteins), " proteins...")
  plots <- list()
  
  for (prot in valid_proteins) {
    tryCatch({
      p <- plot_single_protein_trace(
        protein_id = prot,
        protein_chromatograms = protein_chromatograms,
        peptide_traces = peptide_traces,
        protein_map = protein_map,
        diff_results = diff_results,
        ref_condition = ref_condition,
        treat_condition = treat_condition,
        show_peptides = show_peptides
      )
      
      # Save plots
      base_name <- gsub("[^A-Za-z0-9_]", "_", prot)
      if (format %in% c("pdf", "both")) {
        ggsave(file.path(output_dir, paste0("trace_", base_name, ".pdf")), 
               p, width = 14, height = 8)
      }
      if (format %in% c("png", "both")) {
        ggsave(file.path(output_dir, paste0("trace_", base_name, ".png")), 
               p, width = 14, height = 8, dpi = 150)
      }
      
      plots[[prot]] <- p
      message("  Plotted: ", prot)
      
    }, error = function(e) {
      message("  Warning: Could not plot ", prot, ": ", e$message)
    })
  }
  
  message("Done! Plots saved to: ", output_dir)
  return(invisible(plots))
}


#' Plot traces for significant proteins
#' 
#' @param results Results object from run_protein_level_pipeline()
#' @param pval_threshold Adjusted p-value threshold (default 0.05)
#' @param fc_threshold Absolute log2FC threshold (default 1)
#' @param max_plots Maximum number of plots to create
#' @param output_dir Output directory
#' @param show_peptides Whether to show peptide traces
#' @param format Output format
#' @return Invisibly returns list of ggplot objects
plot_significant_proteins <- function(
    results,
    pval_threshold = 0.05,
    fc_threshold = 1,
    max_plots = 30,
    output_dir = "significant_protein_traces",
    show_peptides = TRUE,
    format = "both"
) {
  
  diff_results <- results$differential_results
  
  # Find significant proteins
  sig_proteins <- diff_results %>%
    filter(!is.na(limma_adj_pvalue) & limma_adj_pvalue < pval_threshold &
           !is.na(limma_logFC) & abs(limma_logFC) > fc_threshold) %>%
    arrange(limma_adj_pvalue) %>%
    distinct(protein_id) %>%
    head(max_plots) %>%
    pull(protein_id)
  
  if (length(sig_proteins) == 0) {
    message("No significant proteins found with adj.P < ", pval_threshold, " and |log2FC| > ", fc_threshold)
    return(invisible(list()))
  }
  
  message("Found ", length(sig_proteins), " significant proteins (adj.P < ", pval_threshold, 
          ", |log2FC| > ", fc_threshold, ")")
  
  plot_protein_traces(results, sig_proteins, output_dir, show_peptides, format)
}


#' Plot traces for high fold-change proteins
#' 
#' @param results Results object from run_protein_level_pipeline()
#' @param fc_threshold Absolute log2FC threshold (default 2)
#' @param max_plots Maximum number of plots to create
#' @param output_dir Output directory
#' @param show_peptides Whether to show peptide traces
#' @param format Output format
#' @return Invisibly returns list of ggplot objects
plot_high_fc_proteins <- function(
    results,
    fc_threshold = 2,
    max_plots = 30,
    output_dir = "high_fc_protein_traces",
    show_peptides = TRUE,
    format = "both"
) {
  
  diff_results <- results$differential_results
  
  # Find high FC proteins (regardless of p-value)
  high_fc_proteins <- diff_results %>%
    filter(!is.na(limma_logFC) & abs(limma_logFC) > fc_threshold) %>%
    arrange(desc(abs(limma_logFC))) %>%
    distinct(protein_id) %>%
    head(max_plots) %>%
    pull(protein_id)
  
  if (length(high_fc_proteins) == 0) {
    message("No proteins found with |log2FC| > ", fc_threshold)
    return(invisible(list()))
  }
  
  message("Found ", length(high_fc_proteins), " proteins with |log2FC| > ", fc_threshold)
  
  plot_protein_traces(results, high_fc_proteins, output_dir, show_peptides, format)
}


#' Plot traces for top proteins by various criteria
#' 
#' @param results Results object from run_protein_level_pipeline()
#' @param criterion Sorting criterion: "pvalue", "fc", "fc_abs"
#' @param n Number of top proteins to plot
#' @param output_dir Output directory
#' @param show_peptides Whether to show peptide traces
#' @param format Output format
#' @return Invisibly returns list of ggplot objects
plot_top_proteins <- function(
    results,
    criterion = c("pvalue", "fc", "fc_abs"),
    n = 20,
    output_dir = NULL,
    show_peptides = TRUE,
    format = "both"
) {
  
  criterion <- match.arg(criterion)
  diff_results <- results$differential_results
  
  # Sort by criterion
  top_proteins <- switch(criterion,
    "pvalue" = {
      diff_results %>%
        filter(!is.na(limma_adj_pvalue)) %>%
        arrange(limma_adj_pvalue)
    },
    "fc" = {
      diff_results %>%
        filter(!is.na(limma_logFC)) %>%
        arrange(desc(limma_logFC))
    },
    "fc_abs" = {
      diff_results %>%
        filter(!is.na(limma_logFC)) %>%
        arrange(desc(abs(limma_logFC)))
    }
  ) %>%
    distinct(protein_id) %>%
    head(n) %>%
    pull(protein_id)
  
  if (is.null(output_dir)) {
    output_dir <- paste0("top_", n, "_by_", criterion, "_traces")
  }
  
  message("Plotting top ", length(top_proteins), " proteins by ", criterion)
  
  plot_protein_traces(results, top_proteins, output_dir, show_peptides, format)
}


#' Combined plot: significant + high FC proteins
#' 
#' @param results Results object from run_protein_level_pipeline()
#' @param pval_threshold Adjusted p-value threshold
#' @param fc_threshold_sig FC threshold for significance
#' @param fc_threshold_high FC threshold for high FC (regardless of p-value)
#' @param max_plots Maximum total plots
#' @param output_dir Output directory
#' @param show_peptides Whether to show peptide traces
#' @param format Output format
#' @return Invisibly returns list of ggplot objects
plot_interesting_proteins <- function(
    results,
    pval_threshold = 0.05,
    fc_threshold_sig = 1,
    fc_threshold_high = 2,
    max_plots = 50,
    output_dir = "interesting_protein_traces",
    show_peptides = TRUE,
    format = "both"
) {
  
  diff_results <- results$differential_results
  
  # Find significant proteins
  sig_proteins <- diff_results %>%
    filter(!is.na(limma_adj_pvalue) & limma_adj_pvalue < pval_threshold &
           !is.na(limma_logFC) & abs(limma_logFC) > fc_threshold_sig) %>%
    arrange(limma_adj_pvalue) %>%
    distinct(protein_id) %>%
    pull(protein_id)
  
  # Find high FC proteins (not already in significant)
  high_fc_proteins <- diff_results %>%
    filter(!is.na(limma_logFC) & abs(limma_logFC) > fc_threshold_high) %>%
    filter(!protein_id %in% sig_proteins) %>%
    arrange(desc(abs(limma_logFC))) %>%
    distinct(protein_id) %>%
    pull(protein_id)
  
  # Combine: significant first, then high FC
  all_proteins <- c(sig_proteins, high_fc_proteins)
  all_proteins <- head(all_proteins, max_plots)
  
  message("Found ", length(sig_proteins), " significant proteins and ", 
          length(high_fc_proteins), " additional high FC proteins")
  message("Plotting ", length(all_proteins), " proteins total")
  
  plot_protein_traces(results, all_proteins, output_dir, show_peptides, format)
}


# ═══════════════════════════════════════════════════════════════════════════════
# USAGE EXAMPLES
# ═══════════════════════════════════════════════════════════════════════════════
#
# # After running the pipeline:
# results <- run_protein_level_pipeline(...)
#
# # 1. Plot specific proteins of interest
# plot_protein_traces(
#   results, 
#   protein_ids = c("P12345", "Q67890", "O12345"),
#   output_dir = "my_proteins",
#   show_peptides = TRUE
# )
#
# # 2. Plot all significant proteins (adj.P < 0.05, |FC| > 1)
# plot_significant_proteins(
#   results,
#   pval_threshold = 0.05,
#   fc_threshold = 1,
#   max_plots = 30,
#   output_dir = "significant_traces"
# )
#
# # 3. Plot high fold-change proteins (|FC| > 2, regardless of p-value)
# plot_high_fc_proteins(
#   results,
#   fc_threshold = 2,
#   max_plots = 30,
#   output_dir = "high_fc_traces"
# )
#
# # 4. Plot top 20 by p-value
# plot_top_proteins(results, criterion = "pvalue", n = 20)
#
# # 5. Plot top 20 by absolute fold change
# plot_top_proteins(results, criterion = "fc_abs", n = 20)
#
# # 6. Plot all "interesting" proteins (significant + high FC)
# plot_interesting_proteins(
#   results,
#   pval_threshold = 0.05,
#   fc_threshold_sig = 1,
#   fc_threshold_high = 2,
#   max_plots = 50,
#   output_dir = "interesting_traces"
# )
#
# # 7. Plot a single protein and get the ggplot object
# p <- plot_single_protein_trace(
#   protein_id = "P12345",
#   protein_chromatograms = results$protein_chromatograms,
#   peptide_traces = results$peptide_traces,
#   protein_map = results$protein_map,
#   diff_results = results$differential_results,
#   ref_condition = "ctrl",
#   treat_condition = "ATP"
# )
# print(p)  # Display in RStudio
#
# ═══════════════════════════════════════════════════════════════════════════════

message("Trace plotting functions loaded!")
message("Available functions:")
message("  - plot_protein_traces(results, protein_ids, ...)")
message("  - plot_significant_proteins(results, ...)")
message("  - plot_high_fc_proteins(results, ...)")
message("  - plot_top_proteins(results, criterion, n, ...)")
message("  - plot_interesting_proteins(results, ...)")
message("  - plot_single_protein_trace(...)")
