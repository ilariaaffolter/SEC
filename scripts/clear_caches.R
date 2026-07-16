# scripts/clear_caches.R
# =============================================================================
# One-call helper to delete per-metabolite CACHE files so a re-render recomputes them.
#
# Each metabolite's report writes cache files under output/PCM_ctrl_vs_<m>/ that later runs reload to
# skip slow steps. After the UniProt-mass fix you want the mass-dependent ones rebuilt on a correct,
# fully-covered fetch - but you do NOT want to throw away the expensive mass-INDEPENDENT caches
# (normalization, protein features, grid search, the differential test), which would make re-runs slow
# again for no reason. This helper clears exactly what you ask for and keeps the rest.
#
# WHAT DEPENDS ON THE UNIPROT MASS (clear these to fix the 0%-monomer issue):
#   uniprot   -> rdata/uniprot.RData                 (the fetch itself; clearing forces a fresh,
#                                                      retry-enabled fetch)
#   assembly  -> rdata/diffAssemblyState.RData,
#                rdata/assembly_hits.rds             (the assembly-state test + its hit sets)
# MASS-INDEPENDENT (kept by default; only removed with caches="all" or by naming them):
#   normalized, proteinFeatures, grid, diff, complex
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "clear_caches.R"))
#   clear_caches()                        # ALL metabolites: drop uniprot + assembly caches (the fix)
#   clear_caches("PEP")                   # just one metabolite
#   clear_caches(c("PEP","ADP","NAD"))    # a few
#   clear_caches(caches = "assembly")     # keep the fetch, drop only the assembly caches
#   clear_caches(caches = "all")          # full reset (also normalization/features/grid/diff/complex)
#   clear_caches(dry_run = TRUE)          # show what WOULD be deleted, delete nothing
# =============================================================================

clear_caches <- function(metabolites = NULL,
                         caches      = c("uniprot", "assembly"),
                         dry_run     = FALSE) {
  suppressPackageStartupMessages(library(here))

  # cache-type -> file path(s), relative to output/PCM_ctrl_vs_<m>/ (note: grid is NOT under rdata/)
  cache_files <- list(
    uniprot         = "rdata/uniprot.RData",
    normalized      = "rdata/pepTracesNormalized.RData",
    proteinFeatures = "rdata/proteinFeatures.RData",
    grid            = "ProtGridSearch.RData",
    diff            = "rdata/protein_DiffExpr_all.RData",
    assembly        = c("rdata/diffAssemblyState.RData", "rdata/assembly_hits.rds"),
    complex         = c("rdata/complex_hits.rds", "rdata/complex_featureVals.RData")
  )

  if (identical(caches, "all")) caches <- names(cache_files)
  unknown <- setdiff(caches, names(cache_files))
  if (length(unknown))
    stop("Unknown cache type(s): ", paste(unknown, collapse = ", "),
         ".\n  Valid types: ", paste(names(cache_files), collapse = ", "), ", or \"all\".")

  # discover the metabolites that have an output folder, if not given explicitly
  if (is.null(metabolites)) {
    dirs        <- basename(list.dirs(here("output"), recursive = FALSE))
    metabolites <- sub("^PCM_ctrl_vs_", "", dirs[grepl("^PCM_ctrl_vs_", dirs)])
  }
  if (length(metabolites) == 0) {
    message("No PCM_ctrl_vs_* output folders found - nothing to clear.")
    return(invisible(character(0)))
  }

  rel     <- unlist(cache_files[caches], use.names = FALSE)
  targets <- unlist(lapply(metabolites, function(m)
    file.path(here("output", paste0("PCM_ctrl_vs_", m)), rel)))
  # the UniProt fetch also has a project-level SHARED cache (one file, reused by every metabolite)
  if ("uniprot" %in% caches) targets <- c(targets, here("output", "uniprot_annotation_shared.RData"))
  present <- targets[file.exists(targets)]

  message("Metabolites: ", paste(metabolites, collapse = ", "),
          "  |  cache types: ", paste(caches, collapse = ", "))
  if (length(present) == 0) {
    message("Nothing to delete (none of those cache files exist yet).")
    return(invisible(character(0)))
  }
  message(if (dry_run) "[dry run] would delete these files:" else "Deleting:")
  for (f in present) message("  ", f)
  if (!dry_run) unlink(present)
  message(if (dry_run) "[dry run] " else "", length(present), " file(s) ",
          if (dry_run) "would be removed." else "removed - re-render to recompute them.")
  invisible(present)
}
