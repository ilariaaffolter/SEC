# scripts/run_complex_all_metabolites.R
# =============================================================================
# Second-pass COMPLEX-level analysis, built ON TOP of run_all_metabolites.R.
#
# WHAT THIS DOES:
#   Re-renders the report once per metabolite with run_complex_analysis = TRUE, which turns ON the
#   "COMPLEX level analysis" section (gated OFF by default). Because the heavy PROTEIN-level steps
#   (findProteinFeatures, cyclic-loess normalization, the differential test) are already CACHED from
#   your run_all_metabolites.R run, they load in seconds and this pass only pays for the
#   complex-specific work (findComplexFeatures + the complex differential test + aggregation). Each
#   render writes output/PCM_ctrl_vs_<m>/rdata/complex_hits.rds (up / down / changed complexes).
#   It then runs the cross-metabolite COMPLEX overlap.
#
# PREREQUISITE:
#   Run scripts/run_all_metabolites.R FIRST, so the protein-level caches + outputs already exist.
#   The report's cache flags (load_proteinFeatures.RData / load_normalized.RData /
#   load_protein_DiffExprPep_and_Protein) must be T in SETUP (they are, by default on this branch)
#   for the protein-level work to be reused; otherwise it recomputes (correct, just slower).
#   The COMPLEX section also needs the Complex Portal file referenced by `complex_portal_file` in
#   SETUP to be present in data/raw/.
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "run_complex_all_metabolites.R"))
# =============================================================================

suppressPackageStartupMessages({
  library(rmarkdown)
  library(here)
})

rmd <- here("analysis", "DiffAnalysis_Ecoli_PCM.Rmd")

# Keep this list in sync with scripts/run_all_metabolites.R.
metabolites <- c("ATP", "PEP", "ADP", "aKG", "NAD", "PGP")

for (m in metabolites) {
  cmp <- paste0("PCM_ctrl_vs_", m)
  out <- here("output", cmp)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  message("=== Complex-level render: ", cmp, " ===")
  rmarkdown::render(
    input       = rmd,
    params      = list(metabolite = m, run_complex_analysis = TRUE),
    output_file = paste0("report_complex_", cmp, ".html"),  # separate name; keeps the protein-level report
    output_dir  = out,
    envir       = new.env(parent = globalenv())
  )
}
message("Done. Complex reports + complex_hits.rds under output/PCM_ctrl_vs_*/")

# Cross-metabolite overlap of the complex-level hit sets (up / down / changed) across all metabolites
# just rendered. Wrapped so a hiccup here never hides the successful per-metabolite renders.
message("=== Cross-metabolite COMPLEX overlap analysis ===")
tryCatch(
  source(here("scripts", "overlap_complex_between_metabolites.R")),
  error = function(e) message("Complex overlap step skipped: ", conditionMessage(e))
)
