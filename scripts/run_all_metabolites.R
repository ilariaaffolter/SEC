# scripts/run_all_metabolites.R
# Render the Ecoli_PCM report once per metabolite, without re-knitting by hand.
# Usage (RStudio console, project open):
#   source(here::here("scripts", "run_all_metabolites.R"))
#
# Each metabolite needs (a) an entry in `metabolite_files` in the .Rmd and
# (b) its .tsv in data/raw/PCM_ctrl_vs_<metabolite>/. Report + results for each
# go to output/PCM_ctrl_vs_<metabolite>/.

suppressPackageStartupMessages({
  library(rmarkdown)
  library(here)
})

rmd <- here("analysis", "DiffAnalysis_Ecoli_PCM.Rmd")  # rename per branch if you wish

# Metabolites to run (must match names in `metabolite_files` in the .Rmd).
# These are the 6 currently uncommented there. If your "8" includes pyr and Phe,
# uncomment them in metabolite_files (lines ~95-96) and add "pyr", "Phe" here.
metabolites <- c("ATP", "PEP", "ADP", "aKG", "NAD", "PGP")

for (m in metabolites) {
  cmp <- paste0("PCM_ctrl_vs_", m)
  out <- here("output", cmp)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  message("=== Rendering ", cmp, " ===")
  rmarkdown::render(
    input       = rmd,
    params      = list(metabolite = m),
    output_file = paste0("report_", cmp, ".html"),
    output_dir  = out,
    envir       = new.env(parent = globalenv())
  )
}
message("Done. Reports + results under output/PCM_ctrl_vs_*/")
