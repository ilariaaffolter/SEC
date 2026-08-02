# scripts/plot_protein_traces.R  (yeast_QTL branch)
# =============================================================================
# PLOT INDIVIDUAL PROTEIN SEC ELUTION TRACES ACROSS COMPARISONS, so you can SEE how one protein's
# profile (position / assembly state) differs between the strain comparisons - the ground truth behind a
# differential / assembly / CCF hit, and a quick way to eyeball a protein of interest everywhere at once.
#
# This is the individual-protein companion to scripts/plot_complex_protein_traces.R (which plots whole
# complexes). It reads each comparison's
#   output/<comparison_id>/RData_for_further_plotting_and_analysis/<comparison_id>_for_plotting.RData
# which holds protein_traces_list (and pepTracesList_filtered for the optional peptide view). The yeast
# save does NOT include design_matrix, so Condition/Replicate are reconstructed from the tracesList sample
# names (Sample_name = "<condition>_<replicate>", split on the LAST underscore) and PRINTED to check; pass
# `design_matrices =` to override.
#
# TWO VIEWS
#   show_peptides = FALSE (default): one panel PER COMPARISON (rows), the protein trace coloured by
#       Condition (the two strains overlaid) - the direct cross-comparison view. Replicates are averaged
#       per condition unless aggregate = "replicate".
#   show_peptides = TRUE: facet_grid(Condition ~ comparison); peptides in grey, the protein trace bold
#       red on top - the QC view (does the protein trace summarise its peptides; do peptides agree).
#
# The MW (kDa) top axis is drawn from fraction_annotation$molecular_weight (the shared SEC calibration),
# exactly as plotFeatures / plot_complex_protein_traces do.
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "plot_protein_traces.R"))
#   plot_protein_traces("P00925")                                  # one protein, every comparison
#   plot_protein_traces(c("P00925","P00560"))                      # several proteins (one page each)
#   plot_protein_traces("P00925", comparisons = c("QTL_83_RM","QTL_83_RM_EtOH"))
#   plot_protein_traces("P00925", show_peptides = TRUE)            # QC view: peptides + protein
#   plot_protein_traces("P00925", aggregate = "replicate")        # overlay replicates, don't average
#
# OUTPUT (output/trace_plots/)
#   protein_traces__<ids>.pdf   one page per protein (facets = comparisons; colour = condition)
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })

.ppt_plotting_file <- function(cmp)
  here("output", cmp, "RData_for_further_plotting_and_analysis", paste0(cmp, "_for_plotting.RData"))

# traces_obj$traces -> long (id, fraction, intensity), optionally subset to `ids`. Fraction columns are
# the integer-named columns (mirrors plot_protein_traces.R on the Ecoli branch).
.ppt_trace_to_long <- function(traces_obj, ids = NULL) {
  dt <- as.data.table(traces_obj$traces)
  frac_cols <- grep("^[0-9]+$", colnames(dt), value = TRUE)
  frac_cols <- frac_cols[order(as.numeric(frac_cols))]
  if (!length(frac_cols) || !"id" %in% names(dt)) return(NULL)
  if (!is.null(ids)) dt <- dt[id %in% ids]
  if (!nrow(dt)) return(NULL)
  long <- melt(dt, id.vars = "id", measure.vars = frac_cols, variable.name = "fraction", value.name = "intensity")
  long[, fraction := as.numeric(as.character(fraction))]
  long[, intensity := as.numeric(intensity)]
  long[]
}

# Condition / Replicate from the tracesList sample names (design_matrix is not saved on this branch).
.ppt_design_from_names <- function(ptl) {
  sn <- names(ptl); if (is.null(sn) || !all(nzchar(sn))) sn <- paste0("sample_", seq_along(ptl))
  data.table(Sample_name = sn, Condition = sub("_[^_]+$", "", sn), Replicate = sub("^.*_([^_]+)$", "\\1", sn))
}

# gather traces for `ids` across all samples of a tracesList, tagged with comparison + Condition/Replicate.
# protein_of (peptide trace_annotation) maps peptide id -> its protein_id for the `protein` column.
.ppt_gather <- function(tl, ids, level, design, cmp, protein_of = NULL) {
  rows <- rbindlist(lapply(names(tl), function(sn) {
    lg <- .ppt_trace_to_long(tl[[sn]], ids = ids); if (is.null(lg)) return(NULL)
    lg[, Sample_name := sn]; lg
  }), use.names = TRUE, fill = TRUE)
  if (!nrow(rows)) return(NULL)
  rows[, `:=`(level = level, comparison = cmp)]
  rows <- merge(rows, design[, .(Sample_name, Condition, Replicate)], by = "Sample_name", all.x = TRUE)
  if (!is.null(protein_of) && all(c("id", "protein_id") %in% names(protein_of)))
    rows[, protein := protein_of$protein_id[match(id, protein_of$id)]]
  else rows[, protein := id]
  rows[]
}

# add the MW (kDa) top axis from a fraction_annotation (id, molecular_weight); else a plain fraction axis.
.ppt_mw_axis <- function(p, frac_ann) {
  fa <- tryCatch(as.data.table(frac_ann), error = function(e) NULL)
  idc <- if (!is.null(fa)) intersect(c("id", "fraction_number", "fraction"), names(fa))[1] else NA
  if (is.null(fa) || is.na(idc) || !"molecular_weight" %in% names(fa))
    return(p + scale_x_continuous())
  fa <- unique(fa[, .(id = as.numeric(get(idc)), molecular_weight = as.numeric(molecular_weight))])[order(id)]
  fa <- fa[is.finite(id)]
  br <- seq(min(fa$id, na.rm = TRUE), max(fa$id, na.rm = TRUE), by = 10)
  mw <- round(fa$molecular_weight[match(br, fa$id)])
  p + scale_x_continuous(breaks = br, sec.axis = dup_axis(trans = ~ ., breaks = br, labels = mw, name = "MW (kDa)"))
}

plot_protein_traces <- function(proteins,
                                comparisons     = NULL,
                                aggregate       = c("condition", "replicate"),
                                show_peptides   = FALSE,
                                design_matrices = NULL,      # optional named list to override reconstruction
                                out_subdir      = "trace_plots",
                                width = NULL, height = NULL,
                                verbose = TRUE) {
  aggregate <- match.arg(aggregate)
  proteins  <- unique(toupper(trimws(as.character(proteins))))
  if (!length(proteins)) stop("Give at least one protein id (UniProt accession).", call. = FALSE)

  if (is.null(comparisons)) {
    dd <- list.dirs(here("output"), recursive = FALSE)
    comparisons <- basename(dd)[vapply(basename(dd), function(c) file.exists(.ppt_plotting_file(c)), logical(1))]
    if (!length(comparisons)) stop("No comparison has a *_for_plotting.RData under output/<id>/RData_for_further_plotting_and_analysis/.", call. = FALSE)
  }
  message("Proteins: ", paste(proteins, collapse = ", "))
  message("Comparisons: ", paste(comparisons, collapse = ", "))

  per_cmp <- list(); frac_ann_ref <- NULL
  for (cmp in comparisons) {
    f <- .ppt_plotting_file(cmp)
    if (!file.exists(f)) { if (verbose) message("  [", cmp, "] no _for_plotting.RData -> skipped."); next }
    e <- new.env(); load(f, envir = e)
    ptl <- if ("protein_traces_list" %in% ls(e)) e$protein_traces_list else NULL
    if (is.null(ptl) || !length(ptl)) { if (verbose) message("  [", cmp, "] no protein_traces_list -> skipped."); next }
    design <- if (!is.null(design_matrices) && !is.null(design_matrices[[cmp]])) as.data.table(design_matrices[[cmp]])
              else if ("design_matrix" %in% ls(e)) as.data.table(e$design_matrix)
              else .ppt_design_from_names(ptl)
    if (!all(c("Sample_name", "Condition", "Replicate") %in% names(design))) design <- .ppt_design_from_names(ptl)
    if (verbose) { message("  [", cmp, "] design (reconstructed unless supplied):"); print(design) }

    prot_long <- .ppt_gather(ptl, proteins, "protein", design, cmp)

    pep_long <- NULL
    if (show_peptides && "pepTracesList_filtered" %in% ls(e)) {
      ptt <- e$pepTracesList_filtered
      ann <- tryCatch(as.data.table(ptt[[1]]$trace_annotation), error = function(x) NULL)
      pep_ids <- if (!is.null(ann) && "protein_id" %in% names(ann)) ann[protein_id %in% proteins]$id else character(0)
      if (length(pep_ids)) pep_long <- .ppt_gather(ptt, pep_ids, "peptides", design, cmp, protein_of = ann)
    }

    if (is.null(frac_ann_ref)) frac_ann_ref <- tryCatch(as.data.table(ptl[[1]]$fraction_annotation), error = function(x) NULL)
    md <- rbindlist(list(prot_long, pep_long), use.names = TRUE, fill = TRUE)
    if (!is.null(md) && nrow(md)) per_cmp[[cmp]] <- md
    else if (verbose) message("  [", cmp, "] none of the requested proteins were quantified here.")
  }

  dat <- rbindlist(per_cmp, use.names = TRUE, fill = TRUE)
  if (!nrow(dat)) stop("Nothing to plot: none of ", paste(proteins, collapse = ", "),
                       " were found in ", paste(comparisons, collapse = ", "), ".", call. = FALSE)

  # average replicates per condition unless the replicate overlay was requested
  if (aggregate == "condition") {
    dat <- dat[, .(intensity = mean(intensity, na.rm = TRUE)),
               by = .(protein, id, level, comparison, Condition, fraction)]
    dat[, grp := interaction(id, Condition, comparison, drop = TRUE)]
  } else {
    dat[, grp := interaction(id, Sample_name, comparison, drop = TRUE)]
  }
  dat[, comparison := factor(comparison, levels = comparisons[comparisons %in% unique(dat$comparison)])]

  outdir <- here("output", out_subdir); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ncmp <- uniqueN(dat$comparison); ncond <- uniqueN(dat$Condition)
  tag  <- if (length(proteins) <= 4) paste(proteins, collapse = "-") else paste0(length(proteins), "proteins")
  fn   <- file.path(outdir, paste0("protein_traces__", tag, ".pdf"))
  W <- if (!is.null(width))  width  else if (show_peptides) 3 + 2.6 * ncmp else 8
  H <- if (!is.null(height)) height else if (show_peptides) 2 + 2.0 * ncond else 2 + 1.7 * ncmp

  one_plot <- function(pp) {
    d <- dat[protein == pp]
    if (!nrow(d)) return(NULL)
    sub <- paste0(if (aggregate == "condition") "replicates averaged per condition" else "one line per replicate",
                  "  |  diamonds/top axis = SEC MW calibration")
    if (show_peptides) {
      np <- uniqueN(d[level == "peptides"]$id)
      p <- ggplot(d, aes(fraction, intensity, group = grp)) +
        geom_line(data = d[level == "peptides"], colour = "grey60", linewidth = 0.3, alpha = 0.6, na.rm = TRUE) +
        geom_line(data = d[level == "protein"],  colour = "firebrick", linewidth = 1.0, na.rm = TRUE) +
        facet_grid(Condition ~ comparison, scales = "free_y") +
        labs(title = paste0("SEC traces: ", pp), subtitle = paste0(np, " peptide(s); grey = peptides, red = protein  |  ", sub),
             x = "fraction", y = "intensity")
    } else {
      p <- ggplot(d[level == "protein"], aes(fraction, intensity, colour = Condition, group = grp)) +
        geom_line(linewidth = 0.9, na.rm = TRUE) +
        facet_grid(comparison ~ ., scales = "free_y") +
        labs(title = paste0("SEC trace: ", pp), subtitle = paste0("colour = condition (strain)  |  ", sub),
             x = "fraction", y = "intensity", colour = "condition")
    }
    .ppt_mw_axis(p, frac_ann_ref) + theme_bw()
  }

  plots <- Filter(Negate(is.null), lapply(proteins, one_plot))
  if (!length(plots)) stop("None of the requested proteins produced a panel.", call. = FALSE)
  grDevices::pdf(fn, width = W, height = H); for (p in plots) print(p); grDevices::dev.off()
  message("Wrote ", fn, "  (", length(plots), " protein page(s))")
  if (verbose) {
    seen <- dat[, .(quantified_in = paste(sort(unique(as.character(comparison))), collapse = ", ")), by = protein]
    message("Where each protein was found:"); print(seen)
  }
  invisible(plots[[1]])
}
