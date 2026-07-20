# scripts/plot_protein_traces.R
# =============================================================================
# MANUAL QC: plot the peptide and protein SEC elution traces of a chosen UniProt id,
# ctrl vs metabolite, so you can SEE whether the chromatographic profile really changes between
# conditions (the ground truth behind an assembly / local-vs-global / diff hit).
#
# It reads each metabolite's
#   output/PCM_ctrl_vs_<m>/RData_for_further_plotting_and_analysis/*_for_plotting.RData
# (written at the end of the report), which now contains the peptide traces (pepTracesList_filtered),
# the protein traces (protein_traces_list) and the design_matrix. Replicates are collapsed to one
# profile PER CONDITION with CCprofiler's own integrateTraceIntensities(..., integrate_within =
# "Condition"), exactly as the cross-correlation section does - so nothing heavy is recomputed.
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "plot_protein_traces.R"))
#   plot_protein_traces("P0A6C1")                          # every metabolite that has been run
#   plot_protein_traces("P0A6C1", metabolites = "ATP")     # one metabolite
#   plot_protein_traces("P0A6C1", metabolites = c("ATP","PEP"))
#   plot_protein_traces("P0A6C1", x_axis = "mw")           # apparent MW (kDa) instead of fraction number
#   plot_protein_traces("P0A6C1", aggregate = "replicate") # overlay each replicate instead of averaging
#
# Each call prints the figure and writes output/trace_plots/<protein>__<metabolites>.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(data.table); library(ggplot2); library(CCprofiler)
})

# traces_obj$traces -> long data.table (id, fraction, intensity), optionally subset to `ids`.
# Mirrors the report's .get_trace_matrix (fraction columns are the integer-named columns).
.trace_to_long <- function(traces_obj, ids = NULL) {
  dt        <- as.data.table(traces_obj$traces)
  frac_cols <- grep("^[0-9]+$", colnames(dt), value = TRUE)
  frac_cols <- frac_cols[order(as.numeric(frac_cols))]
  if (!length(frac_cols)) return(NULL)
  if (!is.null(ids)) dt <- dt[id %in% ids]
  if (!nrow(dt)) return(NULL)
  long <- data.table::melt(dt, id.vars = "id", measure.vars = frac_cols,
                           variable.name = "fraction", value.name = "intensity")
  long[, fraction := as.numeric(as.character(fraction))]
  long[]
}

# fraction_number -> apparent MW (kDa) map from a traces object's fraction_annotation, if available.
.fraction_mw_map <- function(traces_obj) {
  fa <- tryCatch(as.data.table(traces_obj$fraction_annotation), error = function(e) NULL)
  if (is.null(fa) || !"molecular_weight" %in% names(fa)) return(NULL)
  key <- intersect(c("fraction_number", "id", "fraction"), names(fa))[1]
  if (is.na(key)) return(NULL)
  setNames(as.numeric(fa$molecular_weight), as.character(fa[[key]]))
}

plot_protein_traces <- function(protein_id,
                                metabolites = NULL,
                                x_axis      = c("fraction", "mw"),
                                aggregate   = c("condition", "replicate"),
                                out_subdir  = "trace_plots",
                                save_pdf    = TRUE,
                                print_plot  = TRUE) {
  x_axis    <- match.arg(x_axis)
  aggregate <- match.arg(aggregate)
  stopifnot(length(protein_id) == 1L)
  .pid <- as.character(protein_id)

  # discover metabolites from the output folders if not given
  if (is.null(metabolites)) {
    dirs        <- basename(list.dirs(here("output"), recursive = FALSE))
    metabolites <- sub("^PCM_ctrl_vs_", "", dirs[grepl("^PCM_ctrl_vs_", dirs)])
  }
  if (length(metabolites) == 0)
    stop("No PCM_ctrl_vs_* output folders found - render at least one metabolite first.")

  collected <- list()
  for (m in metabolites) {
    fdir <- here("output", paste0("PCM_ctrl_vs_", m), "RData_for_further_plotting_and_analysis")
    f    <- list.files(fdir, pattern = "_for_plotting\\.RData$", full.names = TRUE)
    if (!length(f)) {
      message("[", m, "] no *_for_plotting.RData yet - render this metabolite first; skipping."); next
    }
    e <- new.env(); load(f[1], envir = e)
    need <- c("pepTracesList_filtered", "protein_traces_list", "design_matrix")
    miss <- setdiff(need, ls(e))
    if (length(miss)) {
      message("[", m, "] file lacks ", paste(miss, collapse = ", "),
              " - re-render this metabolite (the save now includes design_matrix); skipping."); next
    }

    # collapse replicates to one profile per condition (CCprofiler's own aggregator), unless replicate view
    agg <- function(tl) {
      if (aggregate == "replicate") return(tl)
      tryCatch(integrateTraceIntensities(tl, design_matrix = e$design_matrix,
                                         integrate_within = "Condition", aggr_fun = "sum"),
               error = function(err) {
                 message("  [", m, "] integrateTraceIntensities failed (", conditionMessage(err),
                         ") - showing per-sample traces for this one."); tl })
    }
    pep_tl  <- agg(e$pepTracesList_filtered)
    prot_tl <- agg(e$protein_traces_list)

    # peptide ids of this protein, read from the peptide trace annotation (has protein_id + id)
    ann     <- as.data.table(pep_tl[[1]]$trace_annotation)
    pep_ids <- if ("protein_id" %in% names(ann)) ann[protein_id == .pid]$id else character(0)

    gather <- function(tl, level, ids) {
      rbindlist(lapply(names(tl), function(g) {
        lg <- .trace_to_long(tl[[g]], ids = ids)
        if (is.null(lg)) return(NULL)
        lg[, `:=`(group = g, level = level, metabolite = m)]
        lg
      }), use.names = TRUE)
    }
    pep_long  <- if (length(pep_ids)) gather(pep_tl,  "peptides", pep_ids) else NULL
    prot_long <- gather(prot_tl, "protein", .pid)   # protein-trace id == protein_id

    md <- rbindlist(list(pep_long, prot_long), use.names = TRUE)
    if (is.null(md) || !nrow(md)) {
      message("[", m, "] protein ", .pid, " not found in the traces; skipping."); next
    }
    # Map each trace's group (condition; or sample when aggregate="replicate") to treatment vs control,
    # so the plot puts the metabolite on top and its matching ctrl below. Control = the condition named
    # like ctrl/control/ref, else the first factor level of design_matrix$Condition (CCprofiler convention).
    .dm    <- as.data.table(e$design_matrix)
    .conds <- as.character(unique(.dm$Condition))
    .ctrl  <- .conds[grepl("ctrl|control|ref", .conds, ignore.case = TRUE)][1]
    if (is.na(.ctrl)) .ctrl <- if (is.factor(.dm$Condition)) as.character(levels(.dm$Condition))[1] else .conds[1]
    .g2c   <- c(setNames(as.character(.dm$Condition), as.character(.dm$Condition)),      # condition -> itself
                setNames(as.character(.dm$Condition), as.character(.dm$Sample_name)))    # sample    -> its condition
    md[, condition_type := ifelse(.g2c[as.character(group)] == .ctrl, "control", "treatment")]

    # attach apparent MW (kDa) if requested and available
    if (x_axis == "mw") {
      mwmap <- .fraction_mw_map(prot_tl[[1]])
      if (!is.null(mwmap)) md[, mw := mwmap[as.character(fraction)]]
    }
    collected[[m]] <- md
  }

  dat <- rbindlist(collected, use.names = TRUE, fill = TRUE)
  if (!nrow(dat)) stop("Nothing to plot for ", .pid, " in: ", paste(metabolites, collapse = ", "),
                       " (protein absent, or those metabolites not rendered with the updated report).")

  # choose x axis (fall back to fraction if MW couldn't be mapped)
  if (x_axis == "mw" && "mw" %in% names(dat) && any(is.finite(dat$mw))) {
    dat[, xval := mw]; xlab <- "apparent MW (kDa)"
  } else {
    if (x_axis == "mw") message("apparent MW not available on these traces - using fraction number instead.")
    dat[, xval := fraction]; xlab <- "fraction"
  }

  dat[, condition_type := factor(condition_type, levels = c("treatment", "control"))]   # metabolite on top, ctrl below
  npep <- length(unique(dat[level == "peptides"]$id))

  # rows = condition (metabolite / ctrl), columns = metabolite; peptides grey, protein trace bold red on top
  p <- ggplot(mapping = aes(x = xval, y = intensity, group = interaction(id, group))) +
    geom_line(data = dat[level == "peptides"], colour = "grey55",    linewidth = 0.35, alpha = 0.7) +
    geom_line(data = dat[level == "protein"],  colour = "firebrick", linewidth = 1.1) +
    facet_grid(condition_type ~ metabolite, scales = "free_y") +
    labs(title    = paste0("SEC traces: ", .pid),
         subtitle = paste0("top = metabolite (treatment), bottom = ctrl  |  ",
                           if (aggregate == "condition") "replicates aggregated per condition" else "one line per replicate",
                           "  |  ", npep, " peptide(s); grey = peptides, red = protein trace"),
         x = xlab, y = "intensity") +
    theme_bw()

  if (save_pdf) {
    outdir <- here("output", out_subdir)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    fn <- file.path(outdir, paste0(.pid, "__", paste(metabolites, collapse = "-"), ".pdf"))
    ggsave(fn, p, width = 3 + 2.4 * length(unique(dat$metabolite)), height = 6, limitsize = FALSE)
    message("Wrote ", fn)
  }
  if (print_plot) print(p)
  invisible(p)   # return the ggplot so callers (e.g. plot_ccf_vs_stat_hits.R) can compose multi-page PDFs
}
