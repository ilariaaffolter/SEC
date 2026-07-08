# scripts/overlap_complex_between_metabolites.R
# =============================================================================
# CROSS-METABOLITE OVERLAP of the per-metabolite COMPLEX-level hit sets.
#
# The COMPLEX twin of scripts/overlap_between_metabolites.R. It reads each metabolite's
# complex_hits.rds (written by the "COMPLEX level analysis" section when the report is rendered with
# run_complex_analysis = TRUE - e.g. via scripts/run_complex_all_metabolites.R) and overlaps the
# differential complexes across metabolites.
#
# WHEN TO RUN:
#   After scripts/run_complex_all_metabolites.R has produced
#   output/PCM_ctrl_vs_<metabolite>/rdata/complex_hits.rds for the metabolites you care about.
#   run_complex_all_metabolites.R also calls this script automatically at the end.
#   To run on its own (project open):
#     source(here::here("scripts", "overlap_complex_between_metabolites.R"))
#
# WHAT IT COMPARES (three hit types, from complex_hits.rds):
#   up       - complexes MORE abundant in the metabolite vs ctrl (treatment-relative medianLog2FC > 0)
#   down     - complexes LESS abundant in the metabolite vs ctrl (medianLog2FC < 0)
#   changed  - any significant complex-level change (pBHadj < 0.05 & |medianLog2FC| > 1)
#
# OUTPUTS (written to output/overlap_complex_between_metabolites/):
#   <type>_presence_matrix.csv        complex x metabolite 0/1 table (+ n_metabolites)
#   <type>_overlap_size_summary.csv   how many complexes are hits in exactly k metabolites
#   <type>_shared_ge2_metabolites.csv complexes that are hits in >=2 metabolites
#   <type>_pairwise_jaccard.pdf       metabolite x metabolite overlap (Jaccard + shared count)
#   <type>_upset.pdf                  UpSet plot (only if the UpSetR package is installed)
#   <type>_shared_heatmap.pdf         fallback presence heatmap of the shared complexes
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(data.table); library(ggplot2)
})

out <- here("output", "overlap_complex_between_metabolites")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

## ---- 1. discover which metabolites have been run --------------------------------------------
comp_dirs   <- list.dirs(here("output"), recursive = FALSE)
comp_dirs   <- comp_dirs[grepl("^PCM_ctrl_vs_", basename(comp_dirs))]
metabolites <- sub("^PCM_ctrl_vs_", "", basename(comp_dirs))

## ---- 2. gather hit sets: hit_type -> (metabolite -> character vector of complex ids) ---------
hit_types <- c("up", "down", "changed")
sets      <- setNames(lapply(hit_types, function(x) list()), hit_types)

for (m in metabolites) {
  f <- here("output", paste0("PCM_ctrl_vs_", m), "rdata", "complex_hits.rds")
  if (file.exists(f)) {
    h <- readRDS(f)
    for (ht in hit_types) sets[[ht]][[m]] <- unique(as.character(h[[ht]]))
  } else {
    message("  [", m, "] complex_hits.rds missing -> run scripts/run_complex_all_metabolites.R for this metabolite.")
  }
}

have <- sort(unique(unlist(lapply(sets, names))))
if (length(have) < 2) {
  stop("Found complex_hits.rds for ", length(have), " metabolite(s). Need >= 2 to compute an overlap - ",
       "run scripts/run_complex_all_metabolites.R first.")
}
message("Cross-metabolite COMPLEX overlap over: ", paste(have, collapse = ", "))

## ---- 3. per hit type: presence matrix, overlap summary, shared tables, plots ----------------
overlap_one <- function(type) {
  s <- sets[[type]]
  s <- s[vapply(s, length, 1L) > 0]                 # drop metabolites with no hits of this type
  if (length(s) < 2) {
    message("[", type, "] fewer than 2 metabolites have hits -> nothing to overlap."); return(invisible(NULL))
  }
  mets     <- names(s)
  universe <- sort(unique(unlist(s)))
  M <- vapply(s, function(v) universe %in% v, logical(length(universe)))   # complexes x metabolites
  rownames(M) <- universe

  # presence/absence table (+ how many metabolites share it)
  tab <- data.table(complex_id = universe, n_metabolites = as.integer(rowSums(M)))
  tab <- cbind(tab, as.data.table(matrix(as.integer(M), nrow(M), ncol(M), dimnames = dimnames(M))))
  setorder(tab, -n_metabolites, complex_id)
  fwrite(tab, file.path(out, paste0(type, "_presence_matrix.csv")))

  # how many complexes are hits in exactly k metabolites
  summ <- tab[, .(n_complexes = .N), by = n_metabolites][order(-n_metabolites)]
  fwrite(summ, file.path(out, paste0(type, "_overlap_size_summary.csv")))

  # complexes shared by >=2 metabolites
  shared <- tab[n_metabolites >= 2]
  fwrite(shared, file.path(out, paste0(type, "_shared_ge2_metabolites.csv")))

  # pairwise Jaccard (fill) + shared count (label)
  K <- length(mets)
  grid <- CJ(i = seq_len(K), j = seq_len(K))
  grid[, shared := mapply(function(i, j) length(intersect(s[[i]], s[[j]])), i, j)]
  grid[, uni    := mapply(function(i, j) length(union(s[[i]], s[[j]])),     i, j)]
  grid[, jaccard := ifelse(uni == 0, 0, shared / uni)]
  grid[, m1 := factor(mets[i], levels = mets)][, m2 := factor(mets[j], levels = mets)]
  g <- ggplot(grid, aes(m1, m2, fill = jaccard)) +
    geom_tile(color = "grey92") + geom_text(aes(label = shared), size = 3) +
    scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
    labs(title = paste0("Pairwise overlap of ", type, " complexes"),
         subtitle = "fill = Jaccard index, number = shared complex count", x = NULL, y = NULL) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(out, paste0(type, "_pairwise_jaccard.pdf")), g, width = 6, height = 5)

  # UpSet if the package is available; otherwise a presence heatmap of the shared complexes
  made_upset <- FALSE
  if (requireNamespace("UpSetR", quietly = TRUE)) {
    df <- as.data.frame(matrix(as.integer(M), nrow(M), ncol(M), dimnames = dimnames(M)))
    tryCatch({
      grDevices::pdf(file.path(out, paste0(type, "_upset.pdf")), width = 8, height = 5)
      print(UpSetR::upset(df, nsets = K, nintersects = NA, order.by = "freq"))
      grDevices::dev.off(); made_upset <- TRUE
    }, error = function(e) { try(grDevices::dev.off(), silent = TRUE)
      message("  UpSet plot failed for ", type, " (", conditionMessage(e), ") - using heatmap instead.") })
  }
  if (!made_upset && nrow(shared)) {
    Msh <- M[shared$complex_id, , drop = FALSE]
    hm  <- as.data.table(as.table(matrix(as.integer(Msh), nrow(Msh), ncol(Msh), dimnames = dimnames(Msh))))
    setnames(hm, c("complex_id", "metabolite", "hit"))
    g2 <- ggplot(hm, aes(metabolite, complex_id, fill = factor(hit))) +
      geom_tile(color = "grey90") +
      scale_fill_manual(values = c(`0` = "white", `1` = "firebrick"), guide = "none") +
      labs(title = paste0(type, ": complexes shared by >=2 metabolites"), x = NULL, y = NULL) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(out, paste0(type, "_shared_heatmap.pdf")), g2,
           width = 6, height = max(3, min(30, 0.12 * nrow(shared) + 2)), limitsize = FALSE)
  }

  message("[", type, "] ", length(universe), " complexes across ", K,
          " metabolites; ", nrow(shared), " shared by >=2.")
  invisible(tab)
}

invisible(lapply(hit_types, overlap_one))
message("Done. Cross-metabolite COMPLEX overlap written to: ", out)
