# ═══════════════════════════════════════════════════════════════════════════════
# ELUTION-AWARE PROTEOFORM CLUSTERING v2
# ═══════════════════════════════════════════════════════════════════════════════
# 
# This function clusters peptides based on INTRA-PROTEIN elution consistency.
# 
# Logic:
#   - Peptides from the same protein should co-elute (same complex)
#   - Calculate log2(peptide_i / reference) per fraction for each peptide
#   - If all peptides track together → same proteoform (protein_A)
#   - If some peptides diverge (ratio differs by ≥ threshold for ≥ N consecutive 
#     fractions) → different proteoform (protein_B, protein_C, ...)
#
# This detects peptides that elute at different positions within the same protein,
# which indicates they belong to different complexes/proteoforms.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

#' Detect proteoforms by intra-protein peptide ratio divergence
#' 
#' @param traces_list List of trace objects (pepTracesList_filtered)
#' @param protein_map Data frame with peptide_id and protein_id columns
#' @param divergence_threshold Minimum log2 ratio difference to be considered divergent (default 1.0)
#' @param min_consecutive_fractions Minimum consecutive fractions with divergence (default 3)
#' @param min_peptides Minimum peptides required per proteoform (default 3)
#' @param min_fractions_with_signal Minimum fractions with signal for a peptide (default 5)
#' @param reference_method How to calculate reference: "median" (default) or "most_abundant"
#' @param verbose Print progress messages (default TRUE)
#' 
#' @return Data frame with peptide_id, protein_id, proteoform_id, cluster_id, is_main columns
#' 
detect_proteoforms_elution_aware <- function(traces_list,
                                              protein_map,
                                              divergence_threshold = 1.0,
                                              min_consecutive_fractions = 3,
                                              min_peptides = 3,
                                              min_fractions_with_signal = 5,
                                              reference_method = c("most_abundant", "median"),
                                              verbose = TRUE) {
  
  reference_method <- match.arg(reference_method)
  
  if (verbose) {
    message("\n", paste(rep("═", 70), collapse = ""))
    message("ELUTION-AWARE PROTEOFORM CLUSTERING (Intra-protein)")
    message(paste(rep("═", 70), collapse = ""))
    message("  Method: Compare peptide ratios WITHIN each protein")
    message("  Reference: ", reference_method, " peptide intensity per fraction")
    message("  Divergence threshold: ", divergence_threshold, " log2 units")
    message("  Min consecutive fractions: ", min_consecutive_fractions)
    message("  Min peptides per proteoform: ", min_peptides)
    message(paste(rep("─", 70), collapse = ""))
  }
  
  # ─────────────────────────────────────────────────────────────────────────────
  # STEP 1: Average traces across all samples (or use one condition)
  # ─────────────────────────────────────────────────────────────────────────────
  
  if (verbose) message("\n  Step 1: Building averaged peptide trace matrix...")
  
  # Helper to extract trace matrix
  extract_trace_matrix_local <- function(trace_obj) {
    trace_dt <- trace_obj$traces
    fraction_cols <- grep("^[0-9]+$", colnames(trace_dt), value = TRUE)
    fraction_cols <- fraction_cols[order(as.numeric(fraction_cols))]
    
    if (is.data.table(trace_dt)) {
      trace_mat <- as.matrix(trace_dt[, ..fraction_cols])
      rownames(trace_mat) <- trace_dt$id
    } else {
      trace_mat <- as.matrix(trace_dt[, fraction_cols, drop = FALSE])
      rownames(trace_mat) <- trace_dt$id
    }
    colnames(trace_mat) <- fraction_cols
    return(trace_mat)
  }
  
  # Get trace matrices from all samples
  trace_mats <- lapply(traces_list, extract_trace_matrix_local)
  
  # Get union of all peptides and fraction info
  all_peptides <- unique(unlist(lapply(trace_mats, rownames)))
  n_fractions <- ncol(trace_mats[[1]])
  fraction_names <- colnames(trace_mats[[1]])
  
  if (verbose) message("    Found ", length(all_peptides), " peptides across ", n_fractions, " fractions")
  
  # Average across all samples
  avg_mat <- matrix(NA_real_, nrow = length(all_peptides), ncol = n_fractions)
  rownames(avg_mat) <- all_peptides
  colnames(avg_mat) <- fraction_names
  
  for (pep in all_peptides) {
    pep_values <- list()
    for (mat in trace_mats) {
      if (pep %in% rownames(mat)) {
        pep_values[[length(pep_values) + 1]] <- mat[pep, ]
      }
    }
    if (length(pep_values) > 0) {
      pep_mat <- do.call(rbind, pep_values)
      avg_mat[pep, ] <- colMeans(pep_mat, na.rm = TRUE)
    }
  }
  avg_mat[is.nan(avg_mat)] <- NA_real_
  
  # Filter peptides with insufficient signal
  peptides_with_signal <- rowSums(!is.na(avg_mat) & avg_mat > 0) >= min_fractions_with_signal
  avg_mat <- avg_mat[peptides_with_signal, , drop = FALSE]
  
  if (verbose) message("    Peptides with >= ", min_fractions_with_signal, " valid fractions: ", nrow(avg_mat))
  
  # ─────────────────────────────────────────────────────────────────────────────
  # STEP 2: Map peptides to proteins and initialize results
  # ─────────────────────────────────────────────────────────────────────────────
  
  if (verbose) message("\n  Step 2: Mapping peptides to proteins...")
  
  peptide_protein <- protein_map %>%
    filter(peptide_id %in% rownames(avg_mat)) %>%
    select(peptide_id, protein_id) %>%
    distinct()
  
  # Initialize results
  results <- data.frame(
    peptide_id = rownames(avg_mat),
    stringsAsFactors = FALSE
  ) %>%
    left_join(peptide_protein, by = "peptide_id") %>%
    filter(!is.na(protein_id)) %>%
    mutate(
      cluster_id = 1,
      proteoform_id = paste0(protein_id, "_A"),
      is_main = TRUE
    )
  
  proteins <- unique(results$protein_id)
  
  if (verbose) message("    Proteins to process: ", length(proteins))
  
  # ─────────────────────────────────────────────────────────────────────────────
  # STEP 3: For each protein, detect divergent peptides
  # ─────────────────────────────────────────────────────────────────────────────
  
  if (verbose) {
    message("\n  Step 3: Detecting intra-protein divergence...")
    pb <- txtProgressBar(min = 0, max = length(proteins), style = 3)
  }
  
  proteins_with_multiple <- 0
  
  for (i in seq_along(proteins)) {
    if (verbose) setTxtProgressBar(pb, i)
    
    prot <- proteins[i]
    prot_peptides <- results$peptide_id[results$protein_id == prot]
    prot_peptides <- prot_peptides[prot_peptides %in% rownames(avg_mat)]
    
    # Need at least 2 peptides to compare
    if (length(prot_peptides) < 2) next
    
    # Get peptide traces for this protein
    prot_traces <- avg_mat[prot_peptides, , drop = FALSE]
    
    # ─────────────────────────────────────────────────────────────────────────
    # Calculate reference trace (median or most abundant peptide)
    # ─────────────────────────────────────────────────────────────────────────
    
    if (reference_method == "median") {
      # Median intensity across peptides per fraction
      ref_trace <- apply(prot_traces, 2, function(x) median(x, na.rm = TRUE))
    } else {
      # Most abundant peptide (by total intensity)
      total_int <- rowSums(prot_traces, na.rm = TRUE)
      ref_peptide <- names(which.max(total_int))
      ref_trace <- prot_traces[ref_peptide, ]
    }
    
    # Handle zeros/NAs in reference
    ref_trace[ref_trace == 0 | is.na(ref_trace)] <- NA
    
    # ─────────────────────────────────────────────────────────────────────────
    # Calculate log2 ratio of each peptide vs reference
    # ─────────────────────────────────────────────────────────────────────────
    
    ratio_matrix <- matrix(NA_real_, nrow = nrow(prot_traces), ncol = ncol(prot_traces))
    rownames(ratio_matrix) <- rownames(prot_traces)
    colnames(ratio_matrix) <- colnames(prot_traces)
    
    for (pep in rownames(prot_traces)) {
      pep_trace <- prot_traces[pep, ]
      # Add small pseudocount to avoid log(0)
      ratio_matrix[pep, ] <- log2((pep_trace + 1) / (ref_trace + 1))
    }
    
    # ─────────────────────────────────────────────────────────────────────────
    # Detect divergence: |ratio| >= threshold for >= N consecutive fractions
    # ─────────────────────────────────────────────────────────────────────────
    
    # Function to find longest run of TRUE values
    find_longest_run <- function(x) {
      if (all(is.na(x)) || !any(x, na.rm = TRUE)) return(0)
      x[is.na(x)] <- FALSE
      rle_result <- rle(x)
      true_runs <- rle_result$lengths[rle_result$values]
      if (length(true_runs) == 0) return(0)
      return(max(true_runs))
    }
    
    # For each peptide, check if it diverges from reference
    divergent_peptides <- character(0)
    
    for (pep in rownames(ratio_matrix)) {
      ratios <- ratio_matrix[pep, ]
      
      # Check for consecutive fractions where |ratio| >= threshold
      is_divergent <- abs(ratios) >= divergence_threshold
      longest_run <- find_longest_run(is_divergent)
      
      if (longest_run >= min_consecutive_fractions) {
        divergent_peptides <- c(divergent_peptides, pep)
      }
    }
    
    # If no divergent peptides, all belong to same proteoform
    if (length(divergent_peptides) == 0) next
    
    # ─────────────────────────────────────────────────────────────────────────
    # Cluster divergent peptides (they might form multiple sub-groups)
    # ─────────────────────────────────────────────────────────────────────────
    
    # Main group: non-divergent peptides
    main_peptides <- setdiff(prot_peptides, divergent_peptides)
    
    # If we have divergent peptides, cluster them by their ratio profiles
    if (length(divergent_peptides) >= 2) {
      # Cluster divergent peptides by correlation of their ratio profiles
      div_ratios <- ratio_matrix[divergent_peptides, , drop = FALSE]
      
      # Handle NAs for correlation
      div_ratios_clean <- div_ratios
      div_ratios_clean[is.na(div_ratios_clean)] <- 0
      
      if (nrow(div_ratios_clean) >= 2) {
        cor_mat <- cor(t(div_ratios_clean), use = "pairwise.complete.obs")
        cor_mat[is.na(cor_mat)] <- 0
        dist_mat <- as.dist(1 - cor_mat)
        
        # NOTE (fix applied during handover): assign the tryCatch() RESULT so the
        # fallback also takes effect on error. Previously div_clusters was set
        # only inside the error handler's own scope, so an hclust/cutree failure
        # left div_clusters unset (or stale from a previous protein). Success
        # path is unchanged.
        div_clusters <- tryCatch({
          hc <- hclust(dist_mat, method = "complete")
          # Use correlation threshold of 0.5 (distance of 0.5) for sub-clustering
          cutree(hc, h = 0.5)
        }, error = function(e) {
          # On error, put all divergent in same cluster
          setNames(rep(1, length(divergent_peptides)), divergent_peptides)
        })
      } else {
        div_clusters <- setNames(1, divergent_peptides)
      }
    } else {
      # Single divergent peptide
      div_clusters <- setNames(1, divergent_peptides)
    }
    
    # ─────────────────────────────────────────────────────────────────────────
    # Assign proteoform labels
    # ─────────────────────────────────────────────────────────────────────────
    
    # Main cluster (non-divergent) gets "A"
    # Divergent clusters get "B", "C", etc.
    
    prot_idx <- which(results$protein_id == prot)
    
    for (idx in prot_idx) {
      pep <- results$peptide_id[idx]
      
      if (pep %in% main_peptides) {
        results$cluster_id[idx] <- 1
        results$proteoform_id[idx] <- paste0(prot, "_A")
        results$is_main[idx] <- TRUE
      } else if (pep %in% names(div_clusters)) {
        cluster_num <- div_clusters[pep]
        results$cluster_id[idx] <- cluster_num + 1  # +1 because main is 1
        results$proteoform_id[idx] <- paste0(prot, "_", LETTERS[cluster_num + 1])
        results$is_main[idx] <- FALSE
      }
    }
    
    # Count if we actually created multiple proteoforms
    prot_proteoforms <- unique(results$proteoform_id[prot_idx])
    if (length(prot_proteoforms) > 1) {
      proteins_with_multiple <- proteins_with_multiple + 1
    }
  }
  
  if (verbose) close(pb)
  
  # ─────────────────────────────────────────────────────────────────────────────
  # STEP 4: Filter proteoforms with too few peptides
  # ─────────────────────────────────────────────────────────────────────────────
  
  if (verbose) message("\n  Step 4: Filtering proteoforms with < ", min_peptides, " peptides...")
  
  proteoform_counts <- results %>%
    group_by(proteoform_id) %>%
    summarise(n_peptides = n(), .groups = "drop")
  
  valid_proteoforms <- proteoform_counts %>%
    filter(n_peptides >= min_peptides) %>%
    pull(proteoform_id)
  
  # For invalid proteoforms, reassign peptides back to main proteoform
  invalid_pf <- setdiff(results$proteoform_id, valid_proteoforms)
  
  if (length(invalid_pf) > 0) {
    for (pf in invalid_pf) {
      # Get the protein ID
      prot <- gsub("_[A-Z]$", "", pf)
      main_pf <- paste0(prot, "_A")
      
      # Reassign to main proteoform
      idx <- which(results$proteoform_id == pf)
      results$proteoform_id[idx] <- main_pf
      results$cluster_id[idx] <- 1
      results$is_main[idx] <- TRUE
    }
  }
  
  # Now filter again - remove proteins where even main proteoform has < min_peptides
  proteoform_counts <- results %>%
    group_by(proteoform_id) %>%
    summarise(n_peptides = n(), .groups = "drop")
  
  valid_proteoforms <- proteoform_counts %>%
    filter(n_peptides >= min_peptides) %>%
    pull(proteoform_id)
  
  n_before <- nrow(results)
  results <- results %>% filter(proteoform_id %in% valid_proteoforms)
  n_after <- nrow(results)
  
  if (verbose && n_before > n_after) {
    message("    Filtered: ", n_before, " → ", n_after, " peptides")
  }
  
  # ─────────────────────────────────────────────────────────────────────────────
  # STEP 5: Add log2_ratio column for downstream compatibility
  # ─────────────────────────────────────────────────────────────────────────────
  
  # Calculate median ratio vs protein reference for each peptide
  results$log2_ratio <- NA_real_
  
  # Recalculate using the final assignments
  for (prot in unique(results$protein_id)) {
    prot_peps <- results$peptide_id[results$protein_id == prot]
    prot_peps <- prot_peps[prot_peps %in% rownames(avg_mat)]
    
    if (length(prot_peps) > 0) {
      prot_traces <- avg_mat[prot_peps, , drop = FALSE]
      ref_trace <- apply(prot_traces, 2, median, na.rm = TRUE)
      ref_trace[ref_trace == 0 | is.na(ref_trace)] <- NA
      
      for (pep in prot_peps) {
        pep_trace <- prot_traces[pep, ]
        ratios <- log2((pep_trace + 1) / (ref_trace + 1))
        results$log2_ratio[results$peptide_id == pep] <- median(ratios, na.rm = TRUE)
      }
    }
  }
  
  # ─────────────────────────────────────────────────────────────────────────────
  # Summary
  # ─────────────────────────────────────────────────────────────────────────────
  
  n_proteoforms <- n_distinct(results$proteoform_id)
  n_proteins <- n_distinct(results$protein_id)
  n_with_multiple <- results %>%
    group_by(protein_id) %>%
    summarise(n_pf = n_distinct(proteoform_id)) %>%
    filter(n_pf > 1) %>%
    nrow()
  
  if (verbose) {
    message("\n", paste(rep("═", 70), collapse = ""))
    message("ELUTION-AWARE CLUSTERING COMPLETE")
    message(paste(rep("═", 70), collapse = ""))
    message("  Total peptides: ", nrow(results))
    message("  Total proteins: ", n_proteins)
    message("  Total proteoforms: ", n_proteoforms)
    message("  Proteins with multiple proteoforms: ", n_with_multiple)
    message(paste(rep("═", 70), collapse = ""))
  }
  
  return(results)
}


# ═══════════════════════════════════════════════════════════════════════════════
# DIAGNOSTIC: Visualize intra-protein divergence
# ═══════════════════════════════════════════════════════════════════════════════

#' Plot peptide traces and ratios for a protein to visualize divergence
#' 
#' @param traces_list List of trace objects
#' @param protein_map Peptide-protein mapping
#' @param protein_id Protein to visualize
#' @param proteoform_mapping Optional: result from detect_proteoforms_elution_aware
#' @param reference_method "median" or "most_abundant"
#' 
plot_protein_divergence <- function(traces_list, protein_map, protein_id,
                                    proteoform_mapping = NULL,
                                    reference_method = "median") {
  
  suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
  })
  
  # Extract and average traces
  extract_trace_matrix_local <- function(trace_obj) {
    trace_dt <- trace_obj$traces
    fraction_cols <- grep("^[0-9]+$", colnames(trace_dt), value = TRUE)
    fraction_cols <- fraction_cols[order(as.numeric(fraction_cols))]
    if (is.data.table(trace_dt)) {
      trace_mat <- as.matrix(trace_dt[, ..fraction_cols])
      rownames(trace_mat) <- trace_dt$id
    } else {
      trace_mat <- as.matrix(trace_dt[, fraction_cols, drop = FALSE])
      rownames(trace_mat) <- trace_dt$id
    }
    return(trace_mat)
  }
  
  trace_mats <- lapply(traces_list, extract_trace_matrix_local)
  all_peptides <- unique(unlist(lapply(trace_mats, rownames)))
  n_fractions <- ncol(trace_mats[[1]])
  
  avg_mat <- matrix(NA_real_, nrow = length(all_peptides), ncol = n_fractions)
  rownames(avg_mat) <- all_peptides
  colnames(avg_mat) <- 1:n_fractions
  
  for (pep in all_peptides) {
    pep_values <- list()
    for (mat in trace_mats) {
      if (pep %in% rownames(mat)) {
        pep_values[[length(pep_values) + 1]] <- mat[pep, ]
      }
    }
    if (length(pep_values) > 0) {
      avg_mat[pep, ] <- colMeans(do.call(rbind, pep_values), na.rm = TRUE)
    }
  }
  
  # Get peptides for this protein
  prot_peptides <- protein_map$peptide_id[protein_map$protein_id == protein_id]
  prot_peptides <- prot_peptides[prot_peptides %in% rownames(avg_mat)]
  
  if (length(prot_peptides) == 0) {
    warning("No peptides found for protein: ", protein_id)
    return(NULL)
  }
  
  prot_traces <- avg_mat[prot_peptides, , drop = FALSE]
  
  # Calculate reference
  if (reference_method == "median") {
    ref_trace <- apply(prot_traces, 2, median, na.rm = TRUE)
  } else {
    total_int <- rowSums(prot_traces, na.rm = TRUE)
    ref_trace <- prot_traces[names(which.max(total_int)), ]
  }
  ref_trace[ref_trace == 0 | is.na(ref_trace)] <- NA
  
  # Calculate ratios
  ratio_matrix <- log2((prot_traces + 1) / (matrix(ref_trace, nrow = nrow(prot_traces), 
                                                    ncol = ncol(prot_traces), byrow = TRUE) + 1))
  
  # Prepare plot data
  trace_long <- as.data.frame(prot_traces) %>%
    mutate(peptide_id = rownames(prot_traces)) %>%
    pivot_longer(-peptide_id, names_to = "fraction", values_to = "intensity") %>%
    mutate(fraction = as.numeric(fraction))
  
  ratio_long <- as.data.frame(ratio_matrix) %>%
    mutate(peptide_id = rownames(ratio_matrix)) %>%
    pivot_longer(-peptide_id, names_to = "fraction", values_to = "log2_ratio") %>%
    mutate(fraction = as.numeric(fraction))
  
  # Add proteoform info if available
  if (!is.null(proteoform_mapping)) {
    pf_info <- proteoform_mapping %>%
      filter(protein_id == !!protein_id) %>%
      select(peptide_id, proteoform_id)
    
    trace_long <- trace_long %>% left_join(pf_info, by = "peptide_id")
    ratio_long <- ratio_long %>% left_join(pf_info, by = "peptide_id")
  } else {
    trace_long$proteoform_id <- protein_id
    ratio_long$proteoform_id <- protein_id
  }
  
  # Plot 1: Raw traces
  p1 <- ggplot(trace_long, aes(x = fraction, y = intensity, 
                                group = peptide_id, color = proteoform_id)) +
    geom_line(alpha = 0.7) +
    labs(title = paste("Peptide traces:", protein_id),
         x = "Fraction", y = "Intensity") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 32),
      legend.title = element_text(size = 32),
      axis.text = element_text(size = 44),
      axis.title = element_text(size = 44)
    )
  
  # Plot 2: Ratios vs reference
  p2 <- ggplot(ratio_long, aes(x = fraction, y = log2_ratio,
                                group = peptide_id, color = proteoform_id)) +
    geom_line(alpha = 0.7) +
    geom_hline(yintercept = c(-1, 0, 1), linetype = c("dashed", "solid", "dashed"),
               color = c("red", "gray50", "red")) +
    labs(title = paste("log2(peptide/reference):", protein_id),
         subtitle = "Dashed lines = ±1 divergence threshold",
         x = "Fraction", y = "log2 ratio") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 32),
      legend.title = element_text(size = 32),
      axis.text = element_text(size = 44),
      axis.title = element_text(size = 44)
    )
  
  return(p1 / p2)
}
