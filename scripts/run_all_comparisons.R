# scripts/run_all_comparisons.R
# Render the yeast_QTL report once per pairwise strain comparison, without re-knitting by hand.
# Usage (RStudio console, project open):
#   source(here::here("scripts", "run_all_comparisons.R"))
#
# Each comparison needs (a) an entry in the `comparisons` registry in the .Rmd and
# (b) its .tsv present in data/raw/ whose filename contains that entry's `pattern`
# (the report finds it automatically). Report + results for each go to output/<comparison_id>/.

suppressPackageStartupMessages({
  library(rmarkdown)
  library(here)
})

rmd <- here("analysis", "DiffAnalysis_yeast_QTL.Rmd")

# Comparisons to run (must match names in the `comparisons` registry in the .Rmd).
# 31 and 83 are the mutant strains; BY and RM the parents (the reference in each pair).
comparisons_to_run <- c("QTL_31_BY", "QTL_31_RM", "QTL_83_BY", "QTL_83_RM", "QTL_RM_BY")

for (cmp in comparisons_to_run) {
  out <- here("output", cmp)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  message("=== Rendering ", cmp, " ===")
  rmarkdown::render(
    input       = rmd,
    params      = list(comparison_id = cmp),
    output_file = paste0("report_", cmp, ".html"),
    output_dir  = out,
    envir       = new.env(parent = globalenv())
  )
}
message("Done. Reports + results under output/QTL_*/")

# Cross-comparison overlap of the hit sets (more/less assembled + protein-level diff) across all
# comparisons just rendered. Wrapped so a hiccup here never hides the successful per-comparison runs.
message("=== Cross-comparison overlap analysis ===")
tryCatch(
  source(here("scripts", "overlap_between_comparisons.R")),
  error = function(e) message("Overlap step skipped: ", conditionMessage(e))
)
