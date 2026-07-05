# scripts/overlap_between_metabolites.R
# =============================================================================
# CROSS-METABOLITE OVERLAP of the per-metabolite hit sets.
#
# WHY THIS IS A SEPARATE SCRIPT (not a chunk in the report):
#   DiffAnalysis_Ecoli_PCM.Rmd is rendered ONCE PER METABOLITE (each render is an
#   isolated R session that only ever holds one metabolite's data). An overlap
#   *between* the 6 metabolites therefore cannot live inside that report - it has
#   to run afterwards and read every metabolite's saved outputs. That is exactly
#   what this script does.
#
# WHEN TO RUN:
#   After run_all_metabolites.R has knitted the metabolites you care about, so
#   their outputs exist under output/PCM_ctrl_vs_<metabolite>/rdata/.
#   run_all_metabolites.R also calls this script automatically at the end.
#   To run on its own (project open):  source(here::here("scripts","overlap_between_metabolites.R"))
#
# WHAT IT COMPARES (three hit types, per your request):
#   more_assembled  - proteins MORE assembled in the metabolite vs ctrl
#   less_assembled  - proteins LESS assembled in the metabolite vs ctrl
#                     (both read from rdata/assembly_hits.rds)
#   protein_diff    - protein-level differential-abundance hits
#                     (rdata/protein_DiffExprProtein_list.RData; same cutoff as the
#                      report's diffProteins: pBHadj < 0.05 & |medianLog2FC| > 1)
#
# OUTPUTS (written to output/overlap_between_metabolites/):
#   <type>_presence_matrix.csv       protein x metabolite 0/1 table (+ gene, n_metabolites)
#   <type>_overlap_size_summary.csv  how many proteins are hits in exactly k metabolites
#   <type>_shared_ge2_metabolites.csv proteins that are hits in >=2 metabolites
#   <type>_pairwise_jaccard.pdf      metabolite x metabolite overlap (Jaccard + shared count)
#   <type>_upset.pdf                 UpSet plot (only if the UpSetR package is installed)
#   <type>_shared_heatmap.pdf        fallback presence heatmap of the shared proteins
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(data.table); library(ggplot2)
})

out <- here("output", "overlap_between_metabolites")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

## ---- 1. discover which metabolites have been run --------------------------------------------
comp_dirs   <- list.dirs(here("output"), recursive = FALSE)
comp_dirs   <- comp_dirs[grepl("^PCM_ctrl_vs_", basename(comp_dirs))]
metabolites <- sub("^PCM_ctrl_vs_", "", basename(comp_dirs))
if (length(metabolites) < 2) {
  stop("Found ", length(metabolites), " metabolite output folder(s) under output/PCM_ctrl_vs_*. ",
       "Need at least 2 to compute an overlap - run scripts/run_all_metabolites.R first.")
}
message("Cross-metabolite overlap over: ", paste(sort(metabolites), collapse = ", "))

## ---- 2. gather hit sets: hit_type -> (metabolite -> character vector of protein ids) --------
hit_types  <- c("more_assembled", "less_assembled", "protein_diff")
sets       <- setNames(lapply(hit_types, function(x) list()), hit_types)
gene_annot <- list()   # feature_id -> gene_names, to annotate the shared tables

for (m in metabolites) {
  rdir <- here("output", paste0("PCM_ctrl_vs_", m), "rdata")

  # assembly hits (more / less)
  f_asm <- file.path(rdir, "assembly_hits.rds")
  if (file.exists(f_asm)) {
    a <- readRDS(f_asm)
    sets$more_assembled[[m]] <- unique(as.character(a$more_assembled))
    sets$less_assembled[[m]] <- unique(as.character(a$less_assembled))
  } else {
    message("  [", m, "] assembly_hits.rds missing -> re-knit this metabolite to include assembly hits.")
  }

  # protein-level differential hits
  f_prot <- file.path(rdir, "protein_DiffExprProtein_list.RData")
  if (file.exists(f_prot)) {
    e <- new.env(); load(f_prot, envir = e)
    pdp <- as.data.table(e$protein_DiffExprProtein)
    if (all(c("feature_id", "pBHadj", "medianLog2FC") %in% names(pdp))) {
      sets$protein_diff[[m]] <- unique(as.character(
        pdp[pBHadj < 0.05 & abs(medianLog2FC) > 1]$feature_id))
      if ("gene_names" %in% names(pdp)) gene_annot[[m]] <- unique(pdp[, .(feature_id, gene_names)])
    } else {
      message("  [", m, "] protein_DiffExprProtein missing expected columns -> skipped.")
    }
  } else {
    message("  [", m, "] protein_DiffExprProtein_list.RData missing -> skipped.")
  }
}

# protein_id -> gene_names lookup (best-effort; NA if unavailable)
gene_lu <- if (length(gene_annot)) unique(rbindlist(gene_annot, fill = TRUE)) else
           data.table(feature_id = character(), gene_names = character())
gene_of <- function(ids) if (!nrow(gene_lu)) rep(NA_character_, length(ids)) else
                         gene_lu$gene_names[match(ids, gene_lu$feature_id)]

## ---- 3. per hit type: presence matrix, overlap summary, shared tables, plots ----------------
overlap_one <- function(type) {
  s <- sets[[type]]
  s <- s[vapply(s, length, 1L) > 0]                 # drop metabolites with no hits of this type
  if (length(s) < 2) {
    message("[", type, "] fewer than 2 metabolites have hits -> nothing to overlap."); return(invisible(NULL))
  }
  mets     <- names(s)
  universe <- sort(unique(unlist(s)))
  M <- vapply(s, function(v) universe %in% v, logical(length(universe)))   # proteins x metabolites
  rownames(M) <- universe

  # presence/absence table (+ gene + how many metabolites share it)
  tab <- data.table(protein_id = universe,
                    gene = gene_of(universe),
                    n_metabolites = as.integer(rowSums(M)))
  tab <- cbind(tab, as.data.table(matrix(as.integer(M), nrow(M), ncol(M), dimnames = dimnames(M))))
  setorder(tab, -n_metabolites, protein_id)
  fwrite(tab, file.path(out, paste0(type, "_presence_matrix.csv")))

  # how many proteins are hits in exactly k metabolites
  summ <- tab[, .(n_proteins = .N), by = n_metabolites][order(-n_metabolites)]
  fwrite(summ, file.path(out, paste0(type, "_overlap_size_summary.csv")))

  # proteins shared by >=2 metabolites
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
    labs(title = paste0("Pairwise overlap of ", type, " hits"),
         subtitle = "fill = Jaccard index, number = shared protein count", x = NULL, y = NULL) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(out, paste0(type, "_pairwise_jaccard.pdf")), g, width = 6, height = 5)

  # UpSet if the package is available; otherwise a presence heatmap of the shared proteins
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
    Msh <- M[shared$protein_id, , drop = FALSE]
    hm  <- as.data.table(as.table(matrix(as.integer(Msh), nrow(Msh), ncol(Msh), dimnames = dimnames(Msh))))
    setnames(hm, c("protein_id", "metabolite", "hit"))
    g2 <- ggplot(hm, aes(metabolite, protein_id, fill = factor(hit))) +
      geom_tile(color = "grey90") +
      scale_fill_manual(values = c(`0` = "white", `1` = "firebrick"), guide = "none") +
      labs(title = paste0(type, ": proteins shared by >=2 metabolites"), x = NULL, y = NULL) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(out, paste0(type, "_shared_heatmap.pdf")), g2,
           width = 6, height = max(3, min(30, 0.12 * nrow(shared) + 2)), limitsize = FALSE)
  }

  message("[", type, "] ", length(universe), " proteins across ", K,
          " metabolites; ", nrow(shared), " shared by >=2.")
  invisible(tab)
}

invisible(lapply(hit_types, overlap_one))
message("Done. Cross-metabolite overlap written to: ", out)
