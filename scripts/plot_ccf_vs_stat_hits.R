# scripts/plot_ccf_vs_stat_hits.R
# =============================================================================
# Visualise the COMPLEMENT hits: the proteins that cross-correlation (CCF) flags as elution-SHIFTED,
# split into the ones the statistical differential test ALSO flags ("mixed") and the ones it MISSES
# ("shift only"). For each such protein it draws the same trace panel as plot_protein_traces()
# (metabolite on top, ctrl below; peptides grey, protein trace red), one protein per PDF page - so you
# can eyeball what CCF catches that CCprofiler's abundance test does not.
#
# INPUT: output/PCM_ctrl_vs_<m>/tables/crosscorrelation_vs_stat.txt  (written by section 7f of the
#        report). Categories are derived from its `shifted` and `stat_sig` columns, so it works with
#        tables you already have. Traces come from that metabolite's *_for_plotting.RData.
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "plot_ccf_vs_stat_hits.R"))
#   plot_ccf_vs_stat_hits()                                  # all metabolites, shift-only + mixed
#   plot_ccf_vs_stat_hits("PEP")                             # one metabolite
#   plot_ccf_vs_stat_hits("PEP", categories = "shift_only")  # only what CCprofiler misses
#   plot_ccf_vs_stat_hits("PEP", max_per_cat = 100)          # raise the per-category page cap
#
# OUTPUT (per metabolite): output/PCM_ctrl_vs_<m>/ccf_vs_stat_traces/
#   <m>_shift_only_traces.pdf   one page per protein CCF flags but the stat test does not
#   <m>_mixed_traces.pdf        one page per protein both methods flag
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })

# reuse the exact trace-panel builder so these plots match plot_protein_traces()
source(here::here("scripts", "plot_protein_traces.R"))

plot_ccf_vs_stat_hits <- function(metabolites = NULL,
                                  categories   = c("shift_only", "mixed"),
                                  max_per_cat  = 40,
                                  x_axis       = c("fraction", "mw"),
                                  aggregate    = c("condition", "replicate"),
                                  out_subdir   = "ccf_vs_stat_traces") {
  x_axis    <- match.arg(x_axis)
  aggregate <- match.arg(aggregate)
  categories <- match.arg(categories, choices = c("shift_only", "mixed"), several.ok = TRUE)

  if (is.null(metabolites)) {
    dirs        <- basename(list.dirs(here("output"), recursive = FALSE))
    metabolites <- sub("^PCM_ctrl_vs_", "", dirs[grepl("^PCM_ctrl_vs_", dirs)])
  }
  if (!length(metabolites)) stop("No PCM_ctrl_vs_* output folders found.")

  for (m in metabolites) {
    f <- here("output", paste0("PCM_ctrl_vs_", m), "tables", "crosscorrelation_vs_stat.txt")
    if (!file.exists(f)) {
      message("[", m, "] no crosscorrelation_vs_stat.txt - render this metabolite (section 7f) first; skipping."); next
    }
    tab <- fread(f)
    if (!all(c("protein_id", "shifted", "stat_sig") %in% names(tab))) {
      message("[", m, "] crosscorrelation_vs_stat.txt lacks protein_id/shifted/stat_sig; skipping."); next
    }
    tab[, `:=`(shifted_b = as.logical(shifted), statsig_b = as.logical(stat_sig))]
    sets <- list(
      shift_only = unique(as.character(tab[shifted_b %in% TRUE & !(statsig_b %in% TRUE)]$protein_id)),
      mixed      = unique(as.character(tab[shifted_b %in% TRUE &   statsig_b %in% TRUE ]$protein_id)))

    for (cat in categories) {
      ids <- sets[[cat]]; ids <- ids[!is.na(ids) & nzchar(ids)]
      if (!length(ids)) { message("[", m, "] no '", cat, "' proteins."); next }
      if (length(ids) > max_per_cat) {
        message("[", m, "] '", cat, "': ", length(ids), " proteins - plotting the first ", max_per_cat,
                " (raise max_per_cat for more).")
        ids <- head(ids, max_per_cat)
      }
      outdir <- here("output", paste0("PCM_ctrl_vs_", m), out_subdir)
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      pdf_file <- file.path(outdir, paste0(m, "_", cat, "_traces.pdf"))
      grDevices::pdf(pdf_file, width = 6, height = 6)
      np <- 0L
      for (pid in ids) {
        ok <- tryCatch({
          g <- plot_protein_traces(pid, metabolites = m, x_axis = x_axis, aggregate = aggregate,
                                   save_pdf = FALSE, print_plot = FALSE)
          if (!is.null(g)) print(g)
          !is.null(g)
        }, error = function(e) { message("   ", pid, ": ", conditionMessage(e)); FALSE })
        if (isTRUE(ok)) np <- np + 1L
      }
      grDevices::dev.off()
      message("[", m, "] '", cat, "': ", np, " protein page(s) -> ", pdf_file)
    }
  }
  invisible(NULL)
}
