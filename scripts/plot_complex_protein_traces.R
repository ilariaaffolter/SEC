# scripts/plot_complex_protein_traces.R
# =============================================================================
# PLOT THE SUBUNIT PROTEIN TRACES OF A CHOSEN COMPLEX, in the SAME visual style as the
# "plot all_complexes" chunk of analysis/DiffAnalysis_yeast_QTL.Rmd - but WITHOUT needing a
# detected complex feature, so it works in EVERY comparison, including ones where the complex
# was not scored (e.g. Complex III in QTL_83_RM_EtOH).
#
# WHY THIS SCRIPT EXISTS
#   The chunk plots complexes with
#       plotFeatures(feature_table = scoredDataAll, traces = protein_traces_list, feature_id = <complex>, ...)
#   plotFeatures draws a *detected feature*: it needs a row for that complex in scoredDataAll. If
#   findComplexFeatures did not detect a co-eluting feature for a complex in a given comparison, that
#   complex is simply absent from scoredDataAll and plotFeatures draws nothing - which is exactly why
#   Complex III plots in QTL_83_RM but not in QTL_83_RM_EtOH. "No feature" is a real result (the subunits
#   did not co-elute tightly enough to be called a feature), but you still want to LOOK at the subunit traces to
#   judge co-elution by eye and compare conditions. This script does that: it plots the raw subunit
#   protein traces themselves, feature or no feature, so the two comparisons are directly comparable.
#
# IT REPRODUCES plotFeatures' look faithfully (verified against the differential-fork source):
#   * facet_grid(Condition ~ Replicate)                      - conditions in rows, replicates in columns
#   * MW (kDa) top axis via dup_axis, labelled from fraction_annotation$molecular_weight every 10 fractions
#   * one coloured line per subunit (protein_id), from protein_traces_list$traces
#   * a filled DIAMOND at the top of every facet at each subunit's MONOMER MW position
#       (trace_annotation$protein_mw -> fraction via calibration$MWtoFraction, or interpolated from the
#        per-fraction MW ladder if no calibration object is supplied) - the monomer_MW = T marker
#   * OPTIONAL apex (solid black) + peak boundary (grey) overlay IF you pass a scored feature table
#     (e.g. scoredDataAll) for that comparison - the only part that requires a detected feature.
#
# INPUTS IT READS (per comparison, no re-render needed)
#   output/<comparison_id>/RData_for_further_plotting_and_analysis/<comparison_id>_for_plotting.RData
#       -> protein_traces_list (a CCprofiler tracesList: one traces object per sample, carrying
#          $traces, $trace_annotation[protein_id, protein_mw], $fraction_annotation[id, molecular_weight]).
#   design_matrix is NOT saved there, so it is reconstructed from the tracesList sample names
#   (Sample_name = "<condition>_<replicate>", split on the LAST underscore) and PRINTED for you to check;
#   pass `design_matrices =` to override.
#   Complex -> subunit protein_ids: parsed from data/raw/<complex_portal_file> (the same Complex Portal
#   export the report uses), or pass an explicit vector of protein_ids as `complex`.
#
# USAGE (RStudio console, project open)
#   source(here::here("scripts", "plot_complex_protein_traces.R"))
#   # by name (case-insensitive substring of the Complex Portal name / alias / description):
#   plot_complex_protein_traces("cytochrome bc1", comparisons = c("QTL_83_RM","QTL_83_RM_EtOH"))
#   plot_complex_protein_traces("complex III",    comparisons = c("QTL_83_RM","QTL_83_RM_EtOH"))
#   # by Complex Portal id (with or without the CPX- prefix):
#   plot_complex_protein_traces("CPX-1620")
#   # by an explicit subunit list:
#   plot_complex_protein_traces(c("P07143","P08067","P00128","P00163","Q02761","P07257","P08525","P22289","P37299","P00127"))
#   # overlay the detected feature (apex+boundaries) where it exists, using your in-session scoredDataAll:
#   plot_complex_protein_traces("cytochrome bc1", feature_tables = list(QTL_83_RM = scoredDataAll))
#   # compare elution SHAPE rather than absolute intensity:
#   plot_complex_protein_traces("cytochrome bc1", normalize = "max")
#
# OUTPUT (output/complex_traces/)
#   <complex_tag>__<comparison_id>.pdf   one faithful plotFeatures-style panel per comparison
#   <complex_tag>__COMBINED.pdf          all comparisons' conditions stacked (rows) x replicates (cols),
#                                        the direct side-by-side comparison
#   <complex_tag>__subunits.csv          the subunits used, their monomer MW, and per comparison whether
#                                        each was quantified + whether a complex feature was detected
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })

# ---- resolve a comparison's saved protein_traces_list -----------------------------------------------
.pct_plotting_file <- function(cmp) {
  here("output", cmp, "RData_for_further_plotting_and_analysis", paste0(cmp, "_for_plotting.RData"))
}
.pct_load_traces <- function(cmp, verbose = TRUE) {
  f <- .pct_plotting_file(cmp)
  if (!file.exists(f)) { if (verbose) message("  [", cmp, "] missing ", f, " -> skipped."); return(NULL) }
  e <- new.env(); load(f, envir = e)
  if (!"protein_traces_list" %in% ls(e)) { if (verbose) message("  [", cmp, "] no protein_traces_list in the file -> skipped."); return(NULL) }
  ptl <- e$protein_traces_list
  if (!length(ptl)) { if (verbose) message("  [", cmp, "] protein_traces_list is empty -> skipped."); return(NULL) }
  ptl
}

# ---- design_matrix reconstructed from the tracesList sample names -----------------------------------
# Sample_name = paste0(condition_id, "_", replicate_id) (Rmd chunk 4a), so split on the LAST underscore.
.pct_design_from_names <- function(ptl) {
  sn <- names(ptl)
  if (is.null(sn) || !all(nzchar(sn))) sn <- paste0("sample_", seq_along(ptl))
  data.table(Sample_name = sn,
             Condition   = sub("_[^_]+$", "", sn),
             Replicate   = sub("^.*_([^_]+)$", "\\1", sn))
}

# ---- MW <-> fraction --------------------------------------------------------------------------------
# prefer a supplied calibration object (exactly as plotFeatures does: calibration$MWtoFraction);
# otherwise interpolate on the per-fraction MW ladder in fraction_annotation (the same calibration curve).
.pct_mw_to_fraction <- function(mw, frac_ann, calibration = NULL) {
  if (!is.null(calibration) && is.list(calibration) && is.function(calibration$MWtoFraction)) {
    v <- tryCatch(as.numeric(calibration$MWtoFraction(mw)), error = function(e) NULL)
    if (!is.null(v) && length(v) == length(mw)) return(v)
  }
  fa <- unique(frac_ann[is.finite(molecular_weight) & molecular_weight > 0, .(id, molecular_weight)])
  if (nrow(fa) < 2) return(rep(NA_real_, length(mw)))
  stats::approx(x = log10(fa$molecular_weight), y = fa$id, xout = log10(mw), rule = 2, ties = mean)$y
}

# ---- one comparison -> long table of the subunit traces --------------------------------------------
# mirrors plotFeatures' tracesList branch: toLongFormat(traces$traces) -> (id, fraction, intensity),
# tagged with Sample_name, then Condition/Replicate from the design matrix.
.pct_long <- function(ptl, subunit_ids, design, normalize = "none") {
  parts <- lapply(seq_along(ptl), function(k) {
    tr <- ptl[[k]]
    m  <- as.data.table(tr$traces)
    idc <- if ("id" %in% names(m)) "id" else names(m)[vapply(m, function(v) !is.numeric(v), logical(1))][1]
    if (is.na(idc)) idc <- names(m)[1]
    long <- melt(m, id.vars = idc, variable.name = "fraction", value.name = "intensity")
    setnames(long, idc, "id")
    long[, fraction := suppressWarnings(as.numeric(as.character(fraction)))]
    long <- long[id %in% subunit_ids & is.finite(fraction)]
    if (!nrow(long)) return(NULL)
    long[, Sample_name := names(ptl)[k]]
    long[]
  })
  L <- rbindlist(parts, fill = TRUE)
  if (!nrow(L)) return(L)
  L[, intensity := as.numeric(intensity)]
  L <- merge(L, design, by = "Sample_name", all.x = TRUE)
  if (identical(normalize, "max")) {
    L[, mx := max(intensity, na.rm = TRUE), by = .(id, Sample_name)]
    L[is.finite(mx) & mx > 0, intensity := intensity / mx][, mx := NULL]
  }
  L[]
}

# ---- the plotFeatures-style figure ------------------------------------------------------------------
.pct_plot <- function(L, frac_ann, subunitMW, title, feature = NULL, legend = FALSE,
                      free_y = FALSE, condition_order = NULL, y_lab = "intensity") {
  if (!is.null(condition_order)) L[, Condition := factor(Condition, levels = condition_order)]
  fa <- as.data.table(frac_ann)
  has_mw <- "molecular_weight" %in% names(fa) && any(is.finite(fa$molecular_weight))
  fa <- if (has_mw) unique(fa[, .(id, molecular_weight)])[order(id)] else unique(fa[, .(id)])[order(id)]
  breaks_frac <- seq(min(fa$id, na.rm = TRUE), max(fa$id, na.rm = TRUE), by = 10)
  scales_arg  <- if (isTRUE(free_y)) "free_y" else "fixed"

  p <- ggplot(L, aes(fraction, intensity, colour = id, group = interaction(id, Sample_name)))
  # optional detected-feature overlay (the only part needing a scored feature)
  if (!is.null(feature) && nrow(feature)) {
    if (all(c("left_pp", "right_pp") %in% names(feature)))
      p <- p + geom_rect(data = feature, inherit.aes = FALSE,
                         aes(xmin = left_pp, xmax = right_pp, ymin = -Inf, ymax = Inf),
                         fill = "grey70", alpha = 0.25)
    if ("apex" %in% names(feature))
      p <- p + geom_vline(data = feature, inherit.aes = FALSE, aes(xintercept = apex),
                          colour = "black", linetype = "solid")
  }
  p <- p +
    geom_line(na.rm = TRUE) +
    { if (nrow(subunitMW)) geom_point(data = subunitMW, inherit.aes = FALSE,
                                      aes(x = fraction, y = Inf, colour = id),
                                      shape = 18, size = 4, alpha = 0.5, na.rm = TRUE) } +
    facet_grid(Condition ~ Replicate, scales = scales_arg) +
    labs(title = title, x = "fraction", y = y_lab, colour = "protein_id") +
    theme_bw()
  if (has_mw) {
    breaks_MW <- round(fa$molecular_weight[match(breaks_frac, fa$id)])
    p <- p + scale_x_continuous(breaks = breaks_frac,
                                sec.axis = dup_axis(trans = ~ ., breaks = breaks_frac,
                                                    labels = breaks_MW, name = "MW (kDa)"))
  } else {
    p <- p + scale_x_continuous(breaks = breaks_frac)
  }
  if (!legend) p <- p + theme(legend.position = "none")
  p
}

# ---- complex -> subunit protein_ids (from the Complex Portal export) --------------------------------
.PCT_ACC_RE <- "^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
.pct_subunits <- function(complex, complex_portal_file = "Complex_portal_559292.tsv", verbose = TRUE) {
  # (1) explicit protein_id vector?
  if (length(complex) > 1L || (length(complex) == 1L && grepl(.PCT_ACC_RE, toupper(trimws(complex))))) {
    ids <- unique(toupper(trimws(as.character(complex))))
    return(list(ids = ids, label = if (length(ids) <= 3) paste(ids, collapse = "_") else paste0(length(ids), "proteins"),
                matched = data.table(complex_id = NA_character_, complex_name = "(explicit protein_id list)")))
  }
  q <- trimws(as.character(complex))
  f <- complex_portal_file
  if (!file.exists(f)) { f2 <- here("data", "raw", complex_portal_file); if (file.exists(f2)) f <- f2 else
    stop("Complex Portal file not found: '", complex_portal_file, "' (also tried data/raw/). ",
         "Pass complex_portal_file= or give an explicit vector of protein_ids.", call. = FALSE) }
  cp <- as.data.table(utils::read.csv(f, sep = "\t", header = TRUE, check.names = FALSE))
  # janitor::clean_names-equivalent: collapse non-alphanumerics to "_" THEN lowercase.
  # (must be [^A-Za-z0-9] - a lowercase-only class would eat every uppercase letter, incl. the leading one.)
  nm <- tolower(gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", names(cp))))
  setnames(cp, names(cp), nm)
  ac  <- grep("complex_ac", nm, value = TRUE)[1]
  nmc <- grep("recommended_name", nm, value = TRUE)[1]
  alc <- grep("aliases", nm, value = TRUE)[1]
  dsc <- grep("^description", nm, value = TRUE)[1]
  # participant/accession column: prefer expanded_participant_list (as the report does), else the
  # "identifiers and stoichiometry of molecules" column. NB avoid bare "identifier" - that also matches
  # taxonomy_identifier, which is NOT the participant column.
  plc <- grep("expanded_participant_list|participant", nm, value = TRUE)[1]
  if (is.na(plc)) plc <- grep("stoichiometry|molecules_in_complex", nm, ignore.case = TRUE, value = TRUE)[1]
  if (is.na(ac) || is.na(plc)) stop("Could not find the complex-id / participant columns in ", basename(f),
                                    ". Columns: ", paste(nm, collapse = ", "), call. = FALSE)
  # match query against id (with/without CPX-), name, alias, description
  hay <- tolower(paste(cp[[ac]], if (!is.na(nmc)) cp[[nmc]] else "", if (!is.na(alc)) cp[[alc]] else "",
                       if (!is.na(dsc)) cp[[dsc]] else ""))
  ql  <- tolower(q); ql_nocpx <- sub("^cpx-?", "", ql)
  hit <- grepl(ql, hay, fixed = TRUE)                              # fixed = safe against regex metachars
  if (grepl("^cpx-?\\d+$", ql) || grepl("^\\d+$", ql))            # only do the id-boundary match for CPX ids
    hit <- hit | grepl(paste0("(^|[^0-9])", ql_nocpx, "([^0-9]|$)"), tolower(cp[[ac]]))
  if (!any(hit)) stop("No complex matched '", q, "' in ", basename(f),
                      ". Try a shorter/different substring of the Complex Portal name.", call. = FALSE)
  sub <- cp[hit]
  if (verbose) { message("  matched ", nrow(sub), " complex(es):")
    for (i in seq_len(nrow(sub))) message("    ", sub[[ac]][i], "  ", if (!is.na(nmc)) sub[[nmc]][i] else "") }
  # parse participant accessions: "P07143(1)|P08067(1)|..." -> P07143 P08067 ... (strip stoichiometry, -PRO_)
  toks <- unlist(strsplit(paste(sub[[plc]], collapse = "|"), "\\|"))
  ids  <- sub("\\(.*$", "", toks)
  ids  <- sub("-PRO_.*$", "", ids)
  ids  <- unique(toupper(trimws(ids[nzchar(ids)])))
  ids  <- ids[grepl(.PCT_ACC_RE, ids)]          # keep protein accessions; drop CHEBI/RNA/other participants
  if (!length(ids)) stop("Matched a complex but none of its participants look like protein accessions ",
                         "(only ligands/RNA?). Pass an explicit vector of protein_ids.", call. = FALSE)
  nm_lab <- if (!is.na(nmc) && nrow(sub) == 1) gsub("[^A-Za-z0-9]+", "_", sub[[nmc]][1]) else gsub("[^A-Za-z0-9]+", "_", q)
  nm_lab <- gsub("^_|_$", "", nm_lab)
  list(ids = ids, label = nm_lab,
       matched = data.table(complex_id = sub[[ac]], complex_name = if (!is.na(nmc)) sub[[nmc]] else NA_character_))
}

# ---- is a complex feature detected in a comparison? (diagnostic only) -------------------------------
.pct_feature_detected <- function(cmp, complex_ids) {
  if (!length(complex_ids) || all(is.na(complex_ids))) return(NA)
  cid <- unique(sub("^CPX-?", "", toupper(as.character(complex_ids))))
  for (fn in c("complex_featureValsFilled.RData", "complexFeatures.RData")) {
    f <- here("output", cmp, "rdata", fn)
    if (!file.exists(f)) next
    e <- new.env(); load(f, envir = e)
    obj <- if (exists("complex_featureValsFilled", envir = e)) e$complex_featureValsFilled else
           if (exists("complexFeatures", envir = e)) e$complexFeatures else NULL
    ci <- tryCatch(unique(sub("^CPX-?", "", toupper(as.character(as.data.table(obj)$complex_id)))), error = function(x) NULL)
    if (!is.null(ci)) return(any(cid %in% ci))
  }
  NA
}

# =====================================================================================================
plot_complex_protein_traces <- function(
    complex,
    comparisons        = NULL,
    design_matrices    = NULL,
    calibration        = NULL,
    feature_tables     = NULL,
    complex_portal_file = "Complex_portal_559292.tsv",
    normalize          = c("none", "max"),
    legend             = FALSE,
    free_y             = NA,
    condition_order    = NULL,
    out_subdir         = "complex_traces",
    width = 12, height = 7, verbose = TRUE) {

  normalize <- match.arg(normalize)
  y_lab <- if (normalize == "max") "intensity (scaled to per-sample max)" else "intensity"
  outdir <- here("output", out_subdir); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # which comparisons: default = every one that has a saved plotting file
  if (is.null(comparisons)) {
    dd <- list.dirs(here("output"), recursive = FALSE)
    comparisons <- basename(dd)[vapply(basename(dd), function(c) file.exists(.pct_plotting_file(c)), logical(1))]
    if (!length(comparisons)) stop("No comparison has a *_for_plotting.RData under output/<id>/RData_for_further_plotting_and_analysis/.", call. = FALSE)
  }
  message("Comparisons: ", paste(comparisons, collapse = ", "))

  # subunits
  su <- .pct_subunits(complex, complex_portal_file, verbose)
  message("Complex label: ", su$label, "  (", length(su$ids), " annotated subunit protein_ids)")

  # per-comparison long tables + a shared subunit-MW / fraction-annotation reference
  longs <- list(); subunitMW_ref <- NULL; frac_ann_ref <- NULL; diag_rows <- list()
  for (cmp in comparisons) {
    ptl <- .pct_load_traces(cmp, verbose); if (is.null(ptl)) next
    design <- if (!is.null(design_matrices) && !is.null(design_matrices[[cmp]])) as.data.table(design_matrices[[cmp]])
              else .pct_design_from_names(ptl)
    if (verbose) { message("  [", cmp, "] design (reconstructed unless supplied):"); print(design) }

    ta  <- as.data.table(ptl[[1]]$trace_annotation)
    fa  <- as.data.table(ptl[[1]]$fraction_annotation)
    if (!"molecular_weight" %in% names(fa)) message("  [", cmp, "] WARNING: fraction_annotation has no molecular_weight -> MW top axis will be blank.")
    frac_ann_ref <- if (is.null(frac_ann_ref)) fa else frac_ann_ref

    # subunit monomer-MW diamonds (protein_mw -> fraction)
    idcol <- if ("protein_id" %in% names(ta)) "protein_id" else "id"
    mwcol <- intersect(c("protein_mw", "Mass", "mass", "molecular_weight", "monomer_mw"), names(ta))[1]
    if (!is.na(mwcol)) {
      smw <- unique(ta[get(idcol) %in% su$ids, .(id = get(idcol), mw = as.numeric(get(mwcol)))])
      smw <- smw[is.finite(mw) & mw > 0]
      if (nrow(smw)) { smw[, fraction := .pct_mw_to_fraction(mw, fa, calibration)]
                       if (is.null(subunitMW_ref)) subunitMW_ref <- smw }
    }

    L <- .pct_long(ptl, su$ids, design, normalize)
    quantified <- if (nrow(L)) uniqueN(L$id) else 0L
    detected <- .pct_feature_detected(cmp, su$matched$complex_id)
    diag_rows[[cmp]] <- data.table(comparison = cmp, subunits_quantified = quantified,
                                   subunits_annotated = length(su$ids),
                                   complex_feature_detected = detected)
    message(sprintf("  [%s] %d/%d subunits quantified | complex feature detected: %s",
                    cmp, quantified, length(su$ids), ifelse(is.na(detected), "unknown", detected)))
    if (!nrow(L)) { message("  [", cmp, "] none of the subunits were quantified -> no panel."); next }
    L[, comparison := cmp]
    longs[[cmp]] <- L

    # per-comparison panel (faithful plotFeatures style; fixed y unless overridden)
    feat <- NULL
    if (!is.null(feature_tables) && !is.null(feature_tables[[cmp]])) {
      ft <- as.data.table(feature_tables[[cmp]])
      cid <- unique(sub("^CPX-?", "", toupper(as.character(su$matched$complex_id))))
      if ("complex_id" %in% names(ft)) feat <- ft[sub("^CPX-?", "", toupper(as.character(complex_id))) %in% cid]
    }
    smw_cmp <- if (!is.null(subunitMW_ref)) subunitMW_ref[id %in% unique(L$id)] else data.table()
    p <- .pct_plot(copy(L), fa, if (is.null(smw_cmp)) data.table() else smw_cmp,
                   title = paste0(su$label, "  -  ", cmp),
                   feature = feat, legend = legend,
                   free_y = isTRUE(free_y), y_lab = y_lab)
    ggsave(file.path(outdir, paste0(su$label, "__", cmp, ".pdf")), p, width = width, height = height)
  }

  if (!length(longs)) stop("No comparison yielded any quantified subunit of '", su$label, "'.", call. = FALSE)

  # combined cross-comparison figure: all conditions (rows) x replicates (cols)
  ALL <- rbindlist(longs, fill = TRUE)
  if (is.null(condition_order)) condition_order <- unique(ALL$Condition)   # order = first appearance
  smw_all <- if (!is.null(subunitMW_ref)) subunitMW_ref[id %in% unique(ALL$id)] else data.table()
  combined_free <- if (is.na(free_y)) TRUE else isTRUE(free_y)             # compare shape by default
  pc <- .pct_plot(copy(ALL), frac_ann_ref, if (is.null(smw_all)) data.table() else smw_all,
                  title = paste0(su$label, "  -  ", paste(names(longs), collapse = " vs "),
                                 if (combined_free) "  (y free per row)" else ""),
                  feature = NULL, legend = legend, free_y = combined_free,
                  condition_order = condition_order, y_lab = y_lab)
  ggsave(file.path(outdir, paste0(su$label, "__COMBINED.pdf")), pc,
         width = width, height = max(height, 2.2 * length(condition_order)))

  # subunit / diagnostic table
  DIAG <- rbindlist(diag_rows, fill = TRUE)
  smw_out <- if (!is.null(subunitMW_ref)) subunitMW_ref[, .(protein_id = id, monomer_mw_kDa = round(mw, 1),
                                                            monomer_fraction = round(fraction, 1))] else
             data.table(protein_id = su$ids)
  fwrite(smw_out, file.path(outdir, paste0(su$label, "__subunits.csv")))
  fwrite(DIAG,    file.path(outdir, paste0(su$label, "__detection.csv")))

  message("\nDone. Written to: ", outdir)
  message("  ", su$label, "__<comparison>.pdf  (one faithful panel each)")
  message("  ", su$label, "__COMBINED.pdf      (direct side-by-side)")
  message("  ", su$label, "__subunits.csv / __detection.csv")
  if (verbose) { message("Detection summary:"); print(DIAG) }
  invisible(list(long = ALL, subunits = su$ids, matched = su$matched, detection = DIAG, plot = pc))
}
