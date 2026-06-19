# ═══════════════════════════════════════════════════════════════════════════════
# DIAGNOSTIC: Check per-fraction valid peptide counts
# Run this AFTER your pipeline to understand why different aggregation methods
# produce identical traces.
# ═══════════════════════════════════════════════════════════════════════════════
#
# Usage:
#   source("diagnose_aggregation.R")
#   diagnose_aggregation_sparsity(
#     traces_list = pepTracesList_filtered_label,
#     protein_map = results_median$protein_map,
#     proteins   = results_median$proteins,
#     custom_ids = c("Q04781", "Q05468", "Q12532", "P25694", "P53044",
#                    "P33755", "Q04311", "P37840", "P42212")
#   )
# ═══════════════════════════════════════════════════════════════════════════════

diagnose_aggregation_sparsity <- function(traces_list,
                                           protein_map,
                                           proteins,
                                           custom_ids = NULL,
                                           min_peptides = 3) {
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(data.table)
  })
  
  # --- helper (same as pipeline) ---
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
  
  is_finite_num <- function(x) !is.na(x) & is.finite(x)
  
  protein_peptide_list <- split(protein_map$peptide_id, protein_map$protein_id)
  
  # Use first sample for diagnostics
  sample_id  <- names(traces_list)[1]
  trace_mat  <- extract_trace_matrix(traces_list[[sample_id]])$matrix
  n_fractions <- ncol(trace_mat)
  
  # --- 1) Global distribution of valid peptides per fraction ---
  cat("\n══════════════════════════════════════════════════════════════\n")
  cat("AGGREGATION SPARSITY DIAGNOSTIC\n")
  cat("Sample: ", sample_id, "\n")
  cat("══════════════════════════════════════════════════════════════\n\n")
  
  all_valid_counts <- c()
  
  for (prot in proteins) {
    prot_peptides <- protein_peptide_list[[prot]]
    prot_peptides <- prot_peptides[prot_peptides %in% rownames(trace_mat)]
    if (length(prot_peptides) < min_peptides) next
    
    pep_mat <- trace_mat[prot_peptides, , drop = FALSE]
    
    for (frac in 1:n_fractions) {
      vals <- pep_mat[, frac]
      n_valid <- sum(is_finite_num(vals) & vals > 0)
      all_valid_counts <- c(all_valid_counts, n_valid)
    }
  }
  
  cat("── Distribution of VALID peptides per protein×fraction ──\n")
  cat("  (This is what the aggregation function actually sees)\n\n")
  tbl <- table(all_valid_counts)
  pct <- round(100 * prop.table(tbl), 1)
  for (i in seq_along(tbl)) {
    cat(sprintf("  %2s valid peptides: %6d entries  (%5.1f%%)\n",
                names(tbl)[i], tbl[i], pct[i]))
  }
  cat(sprintf("\n  >>> Entries with 0-1 valid peptide: %.1f%%\n",
              100 * mean(all_valid_counts <= 1)))
  cat(sprintf("  >>> Entries with 0-2 valid peptides: %.1f%%\n",
              100 * mean(all_valid_counts <= 2)))
  
  if (mean(all_valid_counts <= 1) > 0.5) {
    cat("\n  *** WARNING: >50% of protein×fraction entries have 0-1 valid peptides!\n")
    cat("  *** This explains why median, top10median, top5, etc. give identical traces.\n")
    cat("  *** When only 1 peptide has signal, median == sum == that single value.\n")
  }
  
  # --- 2) Per-protein detail for custom_ids ---
  check_prots <- if (!is.null(custom_ids)) {
    intersect(custom_ids, proteins)
  } else {
    head(proteins, 5)
  }
  
  cat("\n\n── Per-protein breakdown ──\n")
  for (prot in check_prots) {
    prot_peptides <- protein_peptide_list[[prot]]
    prot_peptides <- prot_peptides[prot_peptides %in% rownames(trace_mat)]
    n_total <- length(prot_peptides)
    
    if (n_total < min_peptides) {
      cat(sprintf("\n  %s: only %d peptides (< min_peptides=%d), skipped\n",
                  prot, n_total, min_peptides))
      next
    }
    
    pep_mat <- trace_mat[prot_peptides, , drop = FALSE]
    
    valid_per_frac <- sapply(1:n_fractions, function(f) {
      vals <- pep_mat[, f]
      sum(is_finite_num(vals) & vals > 0)
    })
    
    cat(sprintf("\n  %s: %d peptides mapped, valid per fraction:\n", prot, n_total))
    cat("    Frac: ", sprintf("%3d", 1:n_fractions), "\n")
    cat("    Valid:", sprintf("%3d", valid_per_frac), "\n")
    
    pct_1 <- round(100 * mean(valid_per_frac <= 1), 1)
    cat(sprintf("    → Fractions with ≤1 valid peptide: %d/%d (%.1f%%)\n",
                sum(valid_per_frac <= 1), n_fractions, pct_1))
    
    # Show what the different methods WOULD produce for the peak fraction
    peak_frac <- which.max(valid_per_frac)
    vals_peak <- pep_mat[, peak_frac]
    valid_mask <- is_finite_num(vals_peak) & vals_peak > 0
    valid_vals <- vals_peak[valid_mask]
    
    if (length(valid_vals) > 0) {
      # top10 subset
      pep_totals <- rowSums(pep_mat, na.rm = TRUE)
      n_top10 <- min(10, length(pep_totals))
      top10_peps <- names(sort(pep_totals, decreasing = TRUE))[1:n_top10]
      top10_vals <- pep_mat[top10_peps, peak_frac]
      top10_valid <- top10_vals[is_finite_num(top10_vals) & top10_vals > 0]
      
      # top5 subset
      n_top5 <- min(5, length(pep_totals))
      top5_peps <- names(sort(pep_totals, decreasing = TRUE))[1:n_top5]
      top5_vals <- pep_mat[top5_peps, peak_frac]
      top5_valid <- top5_vals[is_finite_num(top5_vals) & top5_vals > 0]
      
      cat(sprintf("    → At PEAK fraction %d (%d valid peptides):\n", peak_frac, length(valid_vals)))
      cat(sprintf("        median (all %d):      %.1f\n", length(valid_vals), median(valid_vals)))
      cat(sprintf("        top10median (%d):      %.1f\n", length(top10_valid), 
                  if (length(top10_valid) > 0) median(top10_valid) else NA))
      cat(sprintf("        top5 sum (%d):         %.1f\n", length(top5_valid),
                  if (length(top5_valid) > 0) sum(top5_valid) else NA))
      cat(sprintf("        sum (all %d):          %.1f\n", length(valid_vals), sum(valid_vals)))
      
      if (abs(median(valid_vals) - median(top10_valid)) < 0.01 * median(valid_vals)) {
        cat("        ⚠ median ≈ top10median (same peptide pool or near-identical subset)\n")
      }
    }
  }
  
  cat("\n══════════════════════════════════════════════════════════════\n")
  cat("DONE. If most fractions have ≤1-2 valid peptides, all methods\n")
  cat("collapse to the same value. See suggestions below.\n")
  cat("══════════════════════════════════════════════════════════════\n\n")
  
  invisible(all_valid_counts)
}
