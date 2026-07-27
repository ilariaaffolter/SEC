# scripts/globularity_check.R
# =============================================================================
# BASELINE SANITY CHECK, run ON TOP of a rendered report (reads the cached *_for_plotting.RData per
# metabolite - nothing heavy is recomputed).
#
# Question: in the CONTROL condition alone, how many proteins elute where a compact, globular species of
# their molecular weight should elute? Literature expectation is that the large majority of a soluble
# proteome behaves globularly (~95%); proteins that do NOT are either genuinely assembled, genuinely
# extended/disordered, aggregated, degraded - or a chromatography artifact. This script quantifies that
# baseline population BEFORE any metabolite interpretation, so a treatment "shift" is read against a
# known background rather than assumed.
#
# HOW (per protein, CONTROL samples only, replicates averaged):
#   apparent_mw  = the calibrated apparent MW (kDa) at the protein's APEX fraction (its dominant peak).
#                  Comes from the SEC MW calibration already stored in the traces (fraction_annotation$
#                  molecular_weight, built by calibrateMW() from the kDa standards).
#   expected_mw  = the protein's monomer MW (kDa) from UniProt (trace_annotation$protein_mw).
#   ratio        = apparent_mw / expected_mw
#   f/f0         = ratio^(1/3). For a globular protein the Stokes radius scales as M^(1/3) and the MW
#                  calibration is built on globular standards, so the cube root of the MW ratio is the
#                  APPARENT FRICTIONAL RATIO - the standard shape axis: ~1 = as compact as the globular
#                  calibrants, ~1.3 = moderately extended, >1.5 = extended / disordered / large complex.
#                  Two versions are reported:
#                    ffo_vs_monomer = ratio^(1/3)        - assumes the protein is a monomer
#                    ffo_vs_state   = (ratio/n)^(1/3)    - relative to its ASSIGNED oligomer state n
#
# CLASSIFICATION: a protein counts as "globular as expected" when its ratio matches SOME clean oligomer
# state n (1..max_oligomer) within a tolerance factor (default 1.5x, i.e. SEC is coarse) - a compact
# dimer is still globular, just assembled. Everything else is "anomalous":
#   sub_monomer        ratio << 1        (elutes smaller than its own monomer: degradation? interaction
#                                         with the column? mis-annotated MW?)
#   above_range        ratio >> max_oligomer  (far larger than any allowed oligomer: big complex,
#                                         extended/IDR, or aggregate)
#   void               apex in the void fraction(s) - excluded volume, MW cannot be inferred there
#   between_states     matches no n within tolerance (only possible for a tight tolerance)
#
# IMPORTANT CAVEAT (do not over-read): SEC apparent MW conflates SHAPE and ASSEMBLY. Eluting larger than
# the monomer can be a genuine complex OR an extended monomer - this script cannot separate them. The
# "anomalous" set is therefore a list of CANDIDATES for non-globular behaviour, not proof. Separating the
# two needs an independent shape reference (theoretical f/f0 from a structure, e.g. HYDROPRO - planned as
# a follow-up script) or orthogonal evidence.
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "globularity_check.R"))
#   globularity_check()                       # all metabolites' ctrl samples + a pooled summary
#   globularity_check("ATP")                  # one
#   globularity_check(tolerance = 1.3)        # stricter "matches a clean oligomer" window
#   globularity_check(max_oligomer = 6)       # allow up to hexamer as "globular as expected"
#   globularity_check(min_intensity = 300)    # drop low-signal proteins
#
# OUTPUT (per metabolite):
#   tables/globularity_check.txt                    per protein: expected/apparent MW, ratio, f/f0, class
#   figures/globularity_apparent_vs_expected.pdf    log-log scatter + monomer/oligomer reference lines
#   figures/globularity_ffo_distribution.pdf        proteome-wide apparent f/f0 distribution
# OUTPUT (pooled, output/globularity/):
#   globularity_summary.csv                         one row per metabolite + the headline percentages
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })

# traces_obj$traces -> numeric matrix (proteins x fractions), integer-named fraction columns
.get_mat <- function(traces_obj) {
  dt <- as.data.table(traces_obj$traces)
  fc <- grep("^[0-9]+$", colnames(dt), value = TRUE); fc <- fc[order(as.numeric(fc))]
  m  <- as.matrix(dt[, ..fc]); rownames(m) <- as.character(dt$id); m
}

# fraction_number -> apparent MW (kDa) from a traces object's fraction_annotation (written by
# annotateMolecularWeight(); calibrateMW() is built from the kDa standards, so these are kDa).
.fraction_mw_map <- function(traces_obj) {
  fa <- tryCatch(as.data.table(traces_obj$fraction_annotation), error = function(e) NULL)
  if (is.null(fa) || !"molecular_weight" %in% names(fa)) return(NULL)
  key <- intersect(c("fraction_number", "id", "fraction"), names(fa))[1]
  if (is.na(key)) return(NULL)
  setNames(as.numeric(fa$molecular_weight), as.character(fa[[key]]))
}

# protein_id -> monomer MW (kDa) from a protein traces object's trace_annotation
.protein_mw_map <- function(traces_obj) {
  ta <- tryCatch(as.data.table(traces_obj$trace_annotation), error = function(e) NULL)
  if (is.null(ta) || !"protein_mw" %in% names(ta)) return(NULL)
  key <- if ("protein_id" %in% names(ta)) "protein_id" else "id"
  setNames(suppressWarnings(as.numeric(ta$protein_mw)), as.character(ta[[key]]))
}

globularity_check <- function(metabolites   = NULL,
                              tolerance      = 1.5,   # a ratio counts as oligomer n if within this factor of n
                              max_oligomer   = 4L,    # clean oligomer states allowed as "globular as expected"
                              min_intensity  = 0,     # drop proteins below this summed ctrl intensity
                              void_fractions = 1L,    # apex here = excluded volume, MW not interpretable
                              expected_globular_pct = 95,   # literature expectation, for the headline only
                              out_subdir     = "globularity") {
  if (is.null(metabolites)) {
    dirs        <- basename(list.dirs(here("output"), recursive = FALSE))
    metabolites <- sub("^PCM_ctrl_vs_", "", dirs[grepl("^PCM_ctrl_vs_", dirs)])
  }
  if (!length(metabolites)) stop("No PCM_ctrl_vs_* output folders found.")

  summary_rows <- list()
  for (m in metabolites) {
    fdir <- here("output", paste0("PCM_ctrl_vs_", m), "RData_for_further_plotting_and_analysis")
    f    <- list.files(fdir, pattern = "_for_plotting\\.RData$", full.names = TRUE)
    if (!length(f)) { message("[", m, "] no *_for_plotting.RData - render this metabolite first; skipping."); next }
    e <- new.env(); load(f[1], envir = e)
    if (!all(c("protein_traces_list", "design_matrix") %in% ls(e))) {
      message("[", m, "] file lacks protein_traces_list/design_matrix; skipping."); next
    }
    tl <- e$protein_traces_list; samples <- names(tl)

    # CONTROL samples only (ctrl/control/ref-named condition, else the first factor level)
    dm    <- as.data.table(e$design_matrix)
    cond  <- as.character(dm$Condition[match(samples, as.character(dm$Sample_name))])
    conds <- unique(cond)
    ctrl  <- conds[grepl("ctrl|control|ref", conds, ignore.case = TRUE)][1]
    if (is.na(ctrl)) ctrl <- if (is.factor(dm$Condition)) as.character(levels(dm$Condition))[1] else conds[1]
    ctrl_idx <- which(cond == ctrl)
    if (!length(ctrl_idx)) { message("[", m, "] no control samples identified; skipping."); next }

    mats   <- lapply(ctrl_idx, function(i) .get_mat(tl[[i]]))
    common <- Reduce(intersect, lapply(mats, rownames))
    if (length(common) < 10) { message("[", m, "] <10 shared proteins in ctrl; skipping."); next }
    mats <- lapply(mats, function(M) { M <- M[common, , drop = FALSE]; M[is.na(M)] <- 0; M })
    ctrl_mean <- Reduce(`+`, mats) / length(mats)          # mean ctrl profile per protein

    mwmap <- .fraction_mw_map(tl[[ctrl_idx[1]]])
    if (is.null(mwmap)) { message("[", m, "] traces carry no MW calibration (fraction_annotation$molecular_weight); skipping."); next }
    pmw   <- .protein_mw_map(tl[[ctrl_idx[1]]])
    if (is.null(pmw)) { message("[", m, "] traces carry no protein_mw (UniProt monomer mass); skipping."); next }

    fracs     <- as.numeric(colnames(ctrl_mean))
    total_int <- rowSums(ctrl_mean)
    apex_col  <- max.col(ctrl_mean, ties.method = "first")
    apex_frac <- fracs[apex_col]
    apparent  <- unname(mwmap[as.character(apex_frac)])
    expected  <- unname(pmw[common])

    d <- data.table(protein_id = common, apex_fraction = apex_frac, ctrl_intensity = total_int,
                    expected_mw_kDa = expected, apparent_mw_kDa = apparent)
    n_all <- nrow(d)
    d <- d[is.finite(expected_mw_kDa) & expected_mw_kDa > 0 &
           is.finite(apparent_mw_kDa) & apparent_mw_kDa > 0 &
           total_int > min_intensity & is.finite(total_int) & total_int > 0]
    if (!nrow(d)) { message("[", m, "] no protein has both a monomer mass and a usable ctrl peak; skipping."); next }
    d[, ratio := apparent_mw_kDa / expected_mw_kDa]

    # Unit sanity check: both should be kDa. A median ratio near 1000x / 0.001x means a Da-vs-kDa mix-up,
    # which would make every classification meaningless - report it rather than silently classifying.
    .med <- stats::median(d$ratio, na.rm = TRUE)
    if (.med > 100 || .med < 0.01)
      warning(sprintf("[%s] median apparent/expected MW ratio is %.3g - apparent and monomer MW may be in DIFFERENT UNITS (expect both kDa). Interpret with care.", m, .med))

    # assign the closest clean oligomer state n, then test whether it is within the tolerance factor
    ns    <- as.numeric(seq_len(max_oligomer))   # numeric: vapply below is type-strict
    n_best <- vapply(d$ratio, function(r) ns[which.min(abs(log(r / ns)))], numeric(1))
    d[, oligomer_state := n_best]
    d[, dev_from_state := abs(log(ratio / oligomer_state))]          # log-distance to that state
    d[, globular_as_expected := dev_from_state <= log(tolerance)]
    # void apex: excluded volume, apparent MW is not interpretable there
    d[apex_fraction %in% void_fractions, globular_as_expected := FALSE]

    d[, class := data.table::fcase(
      apex_fraction %in% void_fractions,                    "void",
      globular_as_expected & oligomer_state == 1,           "monomer",
      globular_as_expected & oligomer_state >  1,           paste0("oligomer_", oligomer_state, "x"),
      ratio < 1 / tolerance,                                "sub_monomer",
      ratio > max_oligomer * tolerance,                     "above_range",
      default =                                             "between_states")]

    # apparent frictional ratio: cube root of the MW ratio (see header)
    d[, ffo_vs_monomer := ratio^(1/3)]
    d[, ffo_vs_state   := (ratio / oligomer_state)^(1/3)]
    setorder(d, -ratio)

    tab_dir <- here("output", paste0("PCM_ctrl_vs_", m), "tables")
    fig_dir <- here("output", paste0("PCM_ctrl_vs_", m), "figures")
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(d, file.path(tab_dir, "globularity_check.txt"), sep = "\t")

    n_test <- nrow(d); n_glob <- sum(d$globular_as_expected); pct_anom <- 100 * (1 - n_glob / n_test)
    message(sprintf("[%s] ctrl globularity: %d/%d proteins tested (%d had no monomer mass or no usable peak).",
                    m, n_test, n_all, n_all - n_test))
    message(sprintf("[%s]   globular as expected (monomer or clean oligomer <=%dx, within %.2gx): %d (%.1f%%)",
                    m, max_oligomer, tolerance, n_glob, 100 * n_glob / n_test))
    message(sprintf("[%s]   ANOMALOUS: %d (%.1f%%)  vs ~%.0f%% globular expected from literature (i.e. ~%.0f%% anomalous)",
                    m, n_test - n_glob, pct_anom, expected_globular_pct, 100 - expected_globular_pct))
    print(d[, .N, by = class][order(-N)])
    message(sprintf("[%s]   median apparent f/f0 (vs monomer) = %.2f | vs assigned state = %.2f",
                    m, stats::median(d$ffo_vs_monomer, na.rm = TRUE), stats::median(d$ffo_vs_state, na.rm = TRUE)))

    summary_rows[[m]] <- data.table(
      metabolite = m, n_tested = n_test, n_globular = n_glob,
      pct_globular = round(100 * n_glob / n_test, 1), pct_anomalous = round(pct_anom, 1),
      n_monomer      = sum(d$class == "monomer"),
      n_oligomer     = sum(grepl("^oligomer_", d$class)),
      n_sub_monomer  = sum(d$class == "sub_monomer"),
      n_above_range  = sum(d$class == "above_range"),
      n_void         = sum(d$class == "void"),
      median_ffo_vs_monomer = round(stats::median(d$ffo_vs_monomer, na.rm = TRUE), 3),
      median_ffo_vs_state   = round(stats::median(d$ffo_vs_state,   na.rm = TRUE), 3))

    # ---- plot 1: apparent vs expected MW (log-log) with monomer / oligomer reference lines ----
    ref <- data.table(n = ns, label = paste0(ns, "x"))
    g1 <- ggplot(d, aes(expected_mw_kDa, apparent_mw_kDa, colour = class)) +
      # on log10-log10 axes, y = n*x is a line of slope 1 with intercept log10(n)
      geom_abline(data = ref, aes(slope = 1, intercept = log10(n)), linetype = 2, colour = "grey60", inherit.aes = FALSE) +
      geom_point(alpha = 0.5, size = 1) +
      scale_x_log10() + scale_y_log10() +
      annotation_logticks(sides = "bl", colour = "grey70") +
      labs(title = paste0("Control elution vs expected monomer MW - PCM_ctrl_vs_", m),
           subtitle = paste0("dashed = 1x (monomer) to ", max_oligomer, "x (clean oligomers); points off those lines elute anomalously"),
           x = "expected monomer MW (kDa, UniProt)", y = "apparent MW at apex (kDa, SEC calibration)", colour = NULL) +
      theme_bw() + theme(legend.position = "right")
    ggsave(file.path(fig_dir, "globularity_apparent_vs_expected.pdf"), g1, width = 7.5, height = 5.5)

    # ---- plot 2: proteome-wide apparent frictional ratio distribution ----
    vlines <- data.table(x = ns^(1/3), label = paste0(ns, "x"))
    g2 <- ggplot(d, aes(ffo_vs_monomer)) +
      geom_histogram(bins = 80, fill = "steelblue", colour = "white") +
      geom_vline(data = vlines, aes(xintercept = x), linetype = 2, colour = "grey40", inherit.aes = FALSE) +
      geom_text(data = vlines, aes(x = x, y = Inf, label = label), vjust = 1.4, size = 3, colour = "grey30", inherit.aes = FALSE) +
      scale_x_log10() +
      labs(title = paste0("Apparent frictional ratio f/f0 (ctrl) - PCM_ctrl_vs_", m),
           subtitle = "f/f0 = (apparent MW / monomer MW)^(1/3); ~1 = globular monomer. Dashed = clean oligomer states.\nNOTE: assembly and elongation both raise f/f0 - this axis cannot separate them.",
           x = "apparent f/f0 (assuming monomer)", y = "proteins") +
      theme_bw()
    ggsave(file.path(fig_dir, "globularity_ffo_distribution.pdf"), g2, width = 7, height = 5)
  }

  if (length(summary_rows)) {
    S <- rbindlist(summary_rows, use.names = TRUE)
    out <- here("output", out_subdir); dir.create(out, recursive = TRUE, showWarnings = FALSE)
    fwrite(S, file.path(out, "globularity_summary.csv"))
    cat("\n==== control-condition globularity summary ====\n"); print(S)
    cat(sprintf("\nPooled: median %.1f%% of tested proteins elute as expected for a globular species; median %.1f%% anomalous (literature expectation ~%.0f%% globular).\n",
                stats::median(S$pct_globular), stats::median(S$pct_anomalous), expected_globular_pct))
    cat("Reminder: SEC apparent MW conflates assembly and shape - 'anomalous' proteins are CANDIDATES for\n",
        "non-globular behaviour, not proof. A theoretical f/f0 reference (HYDROPRO) would separate the two.\n", sep = "")
    invisible(S)
  } else { message("Nothing computed."); invisible(NULL) }
}
