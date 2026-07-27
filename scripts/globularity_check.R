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
# The classification is run SEPARATELY FOR EACH CONDITION (control and treatment) on the same protein set
# and the same calibration, so the two are directly comparable; the results table carries a `condition`
# column and the summary one row per metabolite x condition. The plots and the pooled pies use the
# CONTROL, which is the baseline the literature expectation refers to.
#
# HOW (per protein, per condition, replicates averaged):
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
# CLASSIFICATION: a protein counts as "globular as expected" when it elutes where SOME clean oligomer
# state n (1..max_oligomer) should - a compact dimer is still globular, just assembled.
# The window is set by `method`:
#   "fractions" (DEFAULT) - the apex may sit up to `tolerance_fractions` (default 1) away from the
#       fraction where that state would elute. This is the natural resolution unit: the apex IS a discrete
#       fraction, and on this calibration one fraction step spans a LARGE MW factor (printed per run;
#       ~1.6-1.7x here). A fixed MW-fold window narrower than one fraction step would flag every protein
#       that is a single fraction off as "anomalous" purely because of the fractionation grid - which is
#       what inflated an earlier version of this check to ~60% anomalous while the median f/f0 was ~1.0.
#   "mwfold" - the older behaviour: the apparent/expected MW ratio must be within `tolerance` (1.5x) of a
#       clean state. Kept for comparison; warns if that window is narrower than one fraction step.
# Everything else is "anomalous":
#   sub_monomer        ratio << 1        (elutes smaller than its own monomer: degradation? interaction
#                                         with the column? mis-annotated MW?)
#   above_range        ratio >> max_oligomer  (far larger than any allowed oligomer: big complex,
#                                         extended/IDR, or aggregate)
#   void               apex in the void fraction(s) - excluded volume, MW cannot be inferred there
#   between_states     matches no n within tolerance (only possible for a tight tolerance)
#   beyond_calibration apparent MW above the LARGEST calibration standard (default 670 kDa). The
#                      calibration is a log-linear fit to the standards, so above the top standard it is
#                      EXTRAPOLATED - in the earliest fractions it runs to physically impossible values
#                      (1e4-1e7 kDa). Those numbers are not measurements. They are reported separately,
#                      and the headline anomalous % is ALSO given for the calibrated range only, which is
#                      the defensible figure.
#
# CALIBRATION STANDARDS are drawn on the scatter (gold diamonds, labelled) by reading each standard's
# observed elution fraction back through the calibration curve. Because the curve was FITTED to those
# points, their offset from the 1x line is the log-linear fit RESIDUAL (a check of how well the fit
# describes the standards, especially at the top end) - not an independent validation.
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
#   globularity_check(calibration_max_kDa = 300)          # treat above 300 kDa as extrapolated
#   globularity_check(calibration_file = "PCM17_calibration_table.xlsx")   # if auto-find fails
#   globularity_check(standards = data.frame(name = c("Thyroglobulin","IgG"), mw_kDa = c(670,150)))
#
# OUTPUT (per metabolite):
#   tables/globularity_check.txt                    per protein: expected/apparent MW, ratio, f/f0, class
#   tables/globularity_standards_check.txt          each standard: expected vs recovered MW (fit residual)
#   figures/globularity_apparent_vs_expected.pdf    log-log scatter + monomer/oligomer reference lines
#   figures/globularity_ffo_distribution.pdf        proteome-wide apparent f/f0 distribution
#   figures/globularity_pies.pdf                    pie 1: globular-as-expected vs anomalous (headline)
#                                                   pie 2: full elution-class breakdown
# OUTPUT (pooled, output/globularity/):
#   globularity_summary.csv                         one row per metabolite + the headline percentages
#   globularity_pies_pooled.pdf                     the same two pies, pooled over all control sets
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })

# ---- pie-chart helpers ----------------------------------------------------------------------------
# canonical class order (globular states first, then the anomalous subtypes) and a colour per class:
# blues = compact/globular (monomer -> higher oligomers get lighter), warm/other = anomalous.
.class_levels <- function(max_oligomer)
  c("monomer", if (max_oligomer >= 2) paste0("oligomer_", 2:max_oligomer, "x"),
    "sub_monomer", "above_range", "between_states", "void", "beyond_calibration",
    "below_calibration_expected", "below_calibration_anomalous")
.class_colours <- function(max_oligomer) {
  olig <- if (max_oligomer >= 2) grDevices::colorRampPalette(c("#6B93C0", "#C6DBEF"))(max_oligomer - 1) else character(0)
  cols <- c("#2C5F8A", olig, "#F28E2B", "#E15759", "#B07AA1", "#9C755F", "#BAB0AC", "#DDD5CE", "#8C6D62")
  setNames(cols, .class_levels(max_oligomer))
}

# ---- SEC calibration standards --------------------------------------------------------------------
# Default = the standards of the kit used here (see SETUP calibration_location). Only used to draw the
# standards on the plot and to set the calibration's upper limit; nothing is re-fitted.
.default_standards <- function() data.table(
  name   = c("Thyroglobulin", "IgA", "IgG", "Ovalbumin", "Myoglobin", "Uridine"),
  mw_kDa = c(670, 300, 150, 44, 17, 0.244))

# Read the calibration table (std_weights_kDa + std_elu_fractions, the columns calibrateMW expects) so
# each standard's OBSERVED elution fraction is known. Auto-finds *calibration*.xlsx in data/raw unless a
# file is given. Returns NULL if unavailable - the standards are then drawn without observed fractions.
.read_calibration <- function(calibration_file = NULL) {
  f <- calibration_file
  if (is.null(f)) {
    cand <- list.files(here("data", "raw"), pattern = "calibration.*\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
    if (!length(cand)) return(NULL)
    f <- cand[1]
  } else if (!file.exists(f)) {
    f2 <- here("data", "raw", f); if (!file.exists(f2)) return(NULL); f <- f2
  }
  if (!requireNamespace("readxl", quietly = TRUE)) return(NULL)
  ct <- tryCatch(as.data.table(readxl::read_excel(f)), error = function(e) NULL)
  if (is.null(ct) || !all(c("std_weights_kDa", "std_elu_fractions") %in% names(ct))) return(NULL)
  ct[, .(mw_kDa = as.numeric(std_weights_kDa), fraction = as.numeric(std_elu_fractions))][is.finite(mw_kDa) & is.finite(fraction)]
}

# apparent MW (kDa) at an arbitrary (possibly non-integer) fraction, from the fraction->MW map already
# stored in the traces. The calibration is log-linear, so interpolating log10(MW) vs fraction is exact.
.mw_at_fraction <- function(mwmap, x) {
  fr <- as.numeric(names(mwmap)); mw <- as.numeric(mwmap)
  ok <- is.finite(fr) & is.finite(mw) & mw > 0
  if (sum(ok) < 2) return(rep(NA_real_, length(x)))
  10^stats::approx(fr[ok], log10(mw[ok]), xout = x, rule = 2)$y
}
# inverse: the fraction at which a given MW (kDa) would elute (log-linear calibration, so linear here)
.fraction_at_mw <- function(mwmap, mw_kDa) {
  fr <- as.numeric(names(mwmap)); mw <- as.numeric(mwmap)
  ok <- is.finite(fr) & is.finite(mw) & mw > 0
  if (sum(ok) < 2) return(rep(NA_real_, length(mw_kDa)))
  o <- order(log10(mw[ok]))
  stats::approx(log10(mw[ok])[o], fr[ok][o], xout = log10(mw_kDa), rule = 2)$y
}
# MW factor spanned by ONE fraction step - the resolution limit of this classification
.mw_per_fraction <- function(mwmap) {
  fr <- as.numeric(names(mwmap)); mw <- as.numeric(mwmap)
  ok <- is.finite(fr) & is.finite(mw) & mw > 0
  if (sum(ok) < 2) return(NA_real_)
  10^abs(stats::coef(stats::lm(log10(mw[ok]) ~ fr[ok]))[2])
}
# write a plot but never let a locked/undeletable file abort the run (OneDrive / an open PDF viewer)
.safe_save <- function(expr, path) {
  tryCatch({ force(expr); TRUE },
           error = function(e) { message("   !! could not write ", basename(path), ": ", conditionMessage(e),
                                         " (is it open in a PDF viewer? continuing)"); FALSE })
}
# one pie from a table of (label, n); slices below label_min_pct are left unlabelled to avoid clutter
.pie <- function(dt, fill_col, cols, title, subtitle, label_min_pct = 3) {
  dt <- copy(dt); data.table::setnames(dt, fill_col, "grp")
  dt[, pct := 100 * n / sum(n)]
  ggplot(dt, aes(x = "", y = n, fill = grp)) +
    geom_col(width = 1, colour = "white", linewidth = 0.3) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = cols, drop = FALSE, name = NULL) +
    geom_text(aes(label = ifelse(pct >= label_min_pct, sprintf("%.1f%%\n(%d)", pct, n), "")),
              position = position_stack(vjust = 0.5), size = 3, colour = "white", fontface = "bold") +
    labs(title = title, subtitle = subtitle) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 8.5), legend.position = "right")
}

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
                              method         = c("fractions", "mwfold"),
                              tolerance_fractions = 1,# "fractions" mode: apex may sit this many fractions
                                                      # off the position expected for a clean oligomer
                              tolerance      = 1.5,   # "mwfold" mode: a ratio counts as oligomer n if within this factor of n
                              max_oligomer   = 4L,    # clean oligomer states allowed as "globular as expected"
                              min_intensity  = 0,     # drop proteins below this summed ctrl intensity
                              void_fractions = 1L,    # apex here = excluded volume, MW not interpretable
                              calibration_max_kDa = NULL,   # NULL = the largest standard; above this the
                                                            # calibration is EXTRAPOLATED, not measured
                              calibration_min_kDa = NULL,   # NULL = the smallest standard >= 1 kDa; BELOW
                                                            # this the calibration is equally extrapolated
                              standards      = NULL,  # data.frame(name, mw_kDa); NULL = the kit defaults
                              calibration_file = NULL,# NULL = auto-find data/raw/*calibration*.xlsx
                              expected_globular_pct = 95,   # literature expectation, for the headline only
                              out_subdir     = "globularity") {
  method <- match.arg(method)
  std <- if (is.null(standards)) .default_standards() else as.data.table(standards)
  if (!all(c("name", "mw_kDa") %in% names(std))) stop("`standards` needs columns name, mw_kDa.")
  cal_obs <- .read_calibration(calibration_file)      # observed elution fraction per standard (or NULL)
  if (is.null(calibration_max_kDa)) calibration_max_kDa <- max(std$mw_kDa, na.rm = TRUE)
  # The LOW end is extrapolated just as much as the high end. Standards below ~1 kDa (e.g. uridine) sit in
  # the total-volume peak, outside the column's resolving range, so they do not constrain the curve -
  # default to the smallest standard that is a real protein.
  if (is.null(calibration_min_kDa)) {
    .sm <- std$mw_kDa[is.finite(std$mw_kDa) & std$mw_kDa >= 1]
    calibration_min_kDa <- if (length(.sm)) min(.sm) else min(std$mw_kDa, na.rm = TRUE)
    if (any(std$mw_kDa < 1, na.rm = TRUE))
      message(sprintf("Standards below 1 kDa (%s) are outside the resolving range and are NOT used as the lower calibration bound; using %.3g kDa.",
                      paste(std$name[std$mw_kDa < 1], collapse = ", "), calibration_min_kDa))
  }
  message(sprintf("Calibrated MW interval: %.3g - %.3g kDa. Outside it the log-linear fit is EXTRAPOLATED and the apparent MW is not a measurement.",
                  calibration_min_kDa, calibration_max_kDa))
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

    # CONDITIONS: the classification is run SEPARATELY for the control and for the treatment, on the same
    # protein set and the same calibration, so the two are directly comparable (the results table carries
    # a `condition` column). The control is always first and drives the plots below.
    dm    <- as.data.table(e$design_matrix)
    cond  <- as.character(dm$Condition[match(samples, as.character(dm$Sample_name))])
    conds <- unique(stats::na.omit(cond))
    ctrl  <- conds[grepl("ctrl|control|ref", conds, ignore.case = TRUE)][1]
    if (is.na(ctrl)) ctrl <- if (is.factor(dm$Condition)) as.character(levels(dm$Condition))[1] else conds[1]
    cond_names <- c(ctrl, setdiff(conds, ctrl))
    if (!length(cond_names)) { message("[", m, "] no conditions identified; skipping."); next }

    # protein set shared by ALL samples, so ctrl and treatment are classified on the same proteins
    all_mats <- lapply(seq_along(samples), function(i) .get_mat(tl[[i]]))
    common   <- Reduce(intersect, lapply(all_mats, rownames))
    if (length(common) < 10) { message("[", m, "] <10 shared proteins; skipping."); next }
    all_mats <- lapply(all_mats, function(M) { M <- M[common, , drop = FALSE]; M[is.na(M)] <- 0; M })

    mwmap <- .fraction_mw_map(tl[[1]])
    if (is.null(mwmap)) { message("[", m, "] traces carry no MW calibration (fraction_annotation$molecular_weight); skipping."); next }
    pmw   <- .protein_mw_map(tl[[1]])
    if (is.null(pmw)) { message("[", m, "] traces carry no protein_mw (UniProt monomer mass); skipping."); next }
    fracs <- as.numeric(colnames(all_mats[[1]]))
    ns    <- as.numeric(seq_len(max_oligomer))   # numeric: vapply below is type-strict

    # RESOLUTION: one fraction step spans this MW factor. If the mwfold tolerance is SMALLER than this,
    # a protein whose apex is a single fraction off is auto-flagged anomalous - which is a property of the
    # fractionation grid, not of the protein. Hence the default "fractions" method.
    mw_per_frac <- .mw_per_fraction(mwmap)
    if (is.finite(mw_per_frac))
      message(sprintf("[%s]   resolution: one fraction step = %.2fx in apparent MW (= %.3fx in f/f0).",
                      m, mw_per_frac, mw_per_frac^(1/3)))
    # the MW-fold equivalent of the accepted window, used for the anomalous SUBTYPE thresholds so they
    # follow whichever method was chosen
    tol_fold <- if (method == "fractions" && is.finite(mw_per_frac)) mw_per_frac^tolerance_fractions else tolerance
    if (method == "fractions" && is.finite(mw_per_frac))
      message(sprintf("[%s]   window: apex within %.2g fraction(s) of the expected position (~%.2fx in MW).",
                      m, tolerance_fractions, tol_fold))
    if (method == "mwfold" && is.finite(mw_per_frac) && tolerance < mw_per_frac)
      warning(sprintf("[%s] tolerance (%.2gx) is NARROWER than one fraction step (%.2fx): proteins one fraction off are flagged anomalous by construction. Use method='fractions' or raise tolerance.",
                      m, tolerance, mw_per_frac))

    # classify one condition group -> its own results table
    .classify_one <- function(cn) {
      idx <- which(cond == cn)
      if (!length(idx)) return(NULL)
      gm        <- Reduce(`+`, all_mats[idx]) / length(idx)     # mean profile per protein in this condition
      total_int <- rowSums(gm)
      apex_frac <- fracs[max.col(gm, ties.method = "first")]
      d <- data.table(condition = cn, protein_id = common, apex_fraction = apex_frac,
                      group_intensity = total_int,
                      expected_mw_kDa = unname(pmw[common]),
                      apparent_mw_kDa = unname(mwmap[as.character(apex_frac)]))
      n_all <- nrow(d)
      d <- d[is.finite(expected_mw_kDa) & expected_mw_kDa > 0 &
             is.finite(apparent_mw_kDa) & apparent_mw_kDa > 0 &
             group_intensity > min_intensity & is.finite(group_intensity) & group_intensity > 0]
      if (!nrow(d)) return(NULL)
      d[, ratio := apparent_mw_kDa / expected_mw_kDa]

      # Unit sanity check: both should be kDa. A median ratio near 1000x / 0.001x means a Da-vs-kDa
      # mix-up, which would make every classification meaningless - report rather than silently classify.
      .med <- stats::median(d$ratio, na.rm = TRUE)
      if (.med > 100 || .med < 0.01)
        warning(sprintf("[%s/%s] median apparent/expected MW ratio is %.3g - apparent and monomer MW may be in DIFFERENT UNITS (expect both kDa).", m, cn, .med))

      if (method == "fractions") {
        # distance, IN FRACTIONS, from the apex to where each clean oligomer state n would elute
        devF <- vapply(ns, function(n) abs(d$apex_fraction - .fraction_at_mw(mwmap, n * d$expected_mw_kDa)),
                       numeric(nrow(d)))
        if (is.null(dim(devF))) devF <- matrix(devF, nrow = nrow(d))
        best_i <- max.col(-devF, ties.method = "first")
        d[, oligomer_state := ns[best_i]]
        d[, dev_fractions  := devF[cbind(seq_len(nrow(d)), best_i)]]
        d[, globular_as_expected := dev_fractions <= tolerance_fractions]
      } else {
        n_best <- vapply(d$ratio, function(r) ns[which.min(abs(log(r / ns)))], numeric(1))
        d[, oligomer_state := n_best]
        d[, dev_fractions  := NA_real_]
        d[, dev_from_state := abs(log(ratio / oligomer_state))]      # log-distance to that state
        d[, globular_as_expected := dev_from_state <= log(tolerance)]
      }
      # void apex: excluded volume, apparent MW is not interpretable there
      d[apex_fraction %in% void_fractions, globular_as_expected := FALSE]
      # BEYOND the largest calibration standard the log-linear fit is EXTRAPOLATED, so the "apparent MW"
      # there is not a measurement (it runs to physically impossible values in the earliest fractions).
      d[, beyond_calibration := apparent_mw_kDa > calibration_max_kDa]
      d[, below_calibration  := apparent_mw_kDa < calibration_min_kDa]
      d[beyond_calibration == TRUE | below_calibration == TRUE, globular_as_expected := FALSE]

      d[, class := data.table::fcase(
        apex_fraction %in% void_fractions,                    "void",
        beyond_calibration == TRUE,                           "beyond_calibration",
        # a protein whose MONOMER is itself below the smallest standard BELONGS below the calibrated
        # range: it is uncalibrated, but its position is not evidence of anomaly. Only a larger protein
        # eluting there is a genuine late-elution anomaly.
        below_calibration & expected_mw_kDa <  calibration_min_kDa, "below_calibration_expected",
        below_calibration & expected_mw_kDa >= calibration_min_kDa, "below_calibration_anomalous",
        globular_as_expected & oligomer_state == 1,           "monomer",
        globular_as_expected & oligomer_state >  1,           paste0("oligomer_", oligomer_state, "x"),
        ratio < 1 / tol_fold,                                 "sub_monomer",
        ratio > max_oligomer * tol_fold,                      "above_range",
        default =                                             "between_states")]

      # proteins whose apex falls INSIDE the calibrated range - the defensible denominator
      d[, in_calibrated_range := !(apex_fraction %in% void_fractions) & !beyond_calibration & !below_calibration]
      # apparent frictional ratio: cube root of the MW ratio (see header)
      d[, ffo_vs_monomer := ratio^(1/3)]
      d[, ffo_vs_state   := (ratio / oligomer_state)^(1/3)]
      setorder(d, -ratio)
      attr(d, "n_all") <- n_all
      d
    }

    d_list <- lapply(cond_names, .classify_one); names(d_list) <- cond_names
    d_list <- d_list[!vapply(d_list, is.null, logical(1))]
    if (!length(d_list)) { message("[", m, "] no condition produced usable data; skipping."); next }
    d_all <- rbindlist(d_list, use.names = TRUE, fill = TRUE)
    d     <- d_list[[1]]          # control: drives the plots below

    tab_dir <- here("output", paste0("PCM_ctrl_vs_", m), "tables")
    fig_dir <- here("output", paste0("PCM_ctrl_vs_", m), "figures")
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(d_all, file.path(tab_dir, "globularity_check.txt"), sep = "\t")

    for (cn in names(d_list)) {
      dc <- d_list[[cn]]; n_all_c <- attr(dc, "n_all")
      role <- if (identical(cn, ctrl)) "control" else "treatment"
      n_test <- nrow(dc); n_glob <- sum(dc$globular_as_expected); pct_anom <- 100 * (1 - n_glob / n_test)
      n_rng  <- sum(dc$in_calibrated_range); n_glob_rng <- sum(dc$globular_as_expected & dc$in_calibrated_range)
      pct_anom_rng <- if (n_rng) 100 * (1 - n_glob_rng / n_rng) else NA_real_
      message(sprintf("[%s | %s (%s)] %d/%d proteins tested (%d had no monomer mass or no usable peak).",
                      m, cn, role, n_test, n_all_c, n_all_c - n_test))
      message(sprintf("   globular as expected (monomer or clean oligomer <=%dx, window: %s): %d (%.1f%%)",
                      max_oligomer,
                      if (method == "fractions") sprintf("+/-%.2g fraction(s) ~ %.2fx MW", tolerance_fractions, tol_fold)
                      else sprintf("%.2gx MW", tolerance),
                      n_glob, 100 * n_glob / n_test))
      message(sprintf("   ANOMALOUS (all tested): %d (%.1f%%)  vs ~%.0f%% globular expected from literature",
                      n_test - n_glob, pct_anom, expected_globular_pct))
      message(sprintf("   >> DEFENSIBLE headline - within the calibrated range (apex <= %.0f kDa, non-void): %d protein(s), ANOMALOUS %.1f%%",
                      calibration_max_kDa, n_rng, pct_anom_rng))
      message(sprintf("   (%d protein(s) elute outside the calibrated interval or in the void: apparent MW there is EXTRAPOLATED, not measured)",
                      n_test - n_rng))
      .bce <- sum(dc$class == "below_calibration_expected"); .bca <- sum(dc$class == "below_calibration_anomalous")
      if (.bce + .bca > 0)
        message(sprintf("   of the %d below the smallest standard: %d are proteins whose MONOMER is itself under %.3g kDa (they belong there - uncalibrated, not anomalous) and %d are larger proteins eluting late (a genuine anomaly).",
                        .bce + .bca, .bce, calibration_min_kDa, .bca))
      print(dc[, .N, by = class][order(-N)])
      message(sprintf("   median apparent f/f0 (vs monomer) = %.2f | vs assigned state = %.2f",
                      stats::median(dc$ffo_vs_monomer, na.rm = TRUE), stats::median(dc$ffo_vs_state, na.rm = TRUE)))

      summary_rows[[paste(m, cn)]] <- data.table(
        metabolite = m, condition = cn, role = role,
        n_tested = n_test, n_globular = n_glob,
        pct_globular = round(100 * n_glob / n_test, 1), pct_anomalous = round(pct_anom, 1),
        # in-range = apex inside the calibrated MW range (the defensible denominator)
        n_in_calibrated_range = n_rng, pct_anomalous_in_range = round(pct_anom_rng, 1),
        n_beyond_calibration  = sum(dc$class == "beyond_calibration"),
        n_below_calibration_expected  = sum(dc$class == "below_calibration_expected"),
        n_below_calibration_anomalous = sum(dc$class == "below_calibration_anomalous"),
        n_globular_in_range   = n_glob_rng,
        n_monomer      = sum(dc$class == "monomer"),
        n_oligomer     = sum(grepl("^oligomer_", dc$class)),
        n_sub_monomer  = sum(dc$class == "sub_monomer"),
        n_above_range  = sum(dc$class == "above_range"),
        n_between_states = sum(dc$class == "between_states"),
        n_void         = sum(dc$class == "void"),
        median_ffo_vs_monomer = round(stats::median(dc$ffo_vs_monomer, na.rm = TRUE), 3),
        median_ffo_vs_state   = round(stats::median(dc$ffo_vs_state,   na.rm = TRUE), 3))
    }

    # ---- plot 1: apparent vs expected MW (log-log) with monomer / oligomer reference lines ----
    # Standards overlay: a standard's apparent MW is its OBSERVED elution fraction read back through the
    # calibration curve. Because the curve was fitted to these points, the offset from the 1x line is the
    # log-linear fit RESIDUAL - it shows how well the calibration describes the standards (especially at
    # the top end), it is not an independent validation.
    std_pts <- copy(std)
    if (!is.null(cal_obs) && nrow(cal_obs)) {
      std_pts[, fraction := cal_obs$fraction[match(round(mw_kDa, 3), round(cal_obs$mw_kDa, 3))]]
      std_pts[, apparent := .mw_at_fraction(mwmap, fraction)]
    } else {
      std_pts[, `:=`(fraction = NA_real_, apparent = NA_real_)]
    }
    std_ok <- std_pts[is.finite(apparent) & is.finite(mw_kDa) & mw_kDa > 0]
    if (!nrow(std_ok))
      message("[", m, "]   (calibration table not found/readable - standards drawn on the x axis only; ",
              "pass calibration_file= to locate it. Needs columns std_weights_kDa + std_elu_fractions.)")

    ref <- data.table(n = ns, label = paste0(ns, "x"))
    g1 <- ggplot(d, aes(expected_mw_kDa, apparent_mw_kDa, colour = class)) +
      # everything above the largest standard is extrapolated calibration, not a measurement
      annotate("rect", xmin = 0, xmax = Inf, ymin = calibration_max_kDa, ymax = Inf,
               fill = "grey70", alpha = 0.22) +
      geom_hline(yintercept = calibration_max_kDa, linetype = 3, colour = "grey30") +
      # on log10-log10 axes, y = n*x is a line of slope 1 with intercept log10(n)
      geom_abline(data = ref, aes(slope = 1, intercept = log10(n)), linetype = 2, colour = "grey60", inherit.aes = FALSE) +
      geom_point(alpha = 0.5, size = 1) +
      scale_x_log10() + scale_y_log10() +
      annotation_logticks(sides = "bl", colour = "grey70") +
      labs(title = paste0("Control elution vs expected monomer MW - PCM_ctrl_vs_", m),
           subtitle = paste0("dashed = 1x (monomer) to ", max_oligomer, "x (clean oligomers); points off those lines elute anomalously.\n",
                             "Shaded = above the largest standard (", calibration_max_kDa,
                             " kDa): the MW calibration is EXTRAPOLATED there, so those values are not measurements.\n",
                             "Diamonds = calibration standards (offset from 1x = the log-linear fit residual)."),
           x = "expected monomer MW (kDa, UniProt)", y = "apparent MW at apex (kDa, SEC calibration)", colour = NULL)
    if (nrow(std_ok)) {
      g1 <- g1 +
        geom_point(data = std_ok, aes(mw_kDa, apparent), shape = 23, size = 3.2,
                   fill = "gold", colour = "black", stroke = 0.7, inherit.aes = FALSE) +
        geom_text(data = std_ok, aes(mw_kDa, apparent, label = name), inherit.aes = FALSE,
                  vjust = -1.1, size = 2.9, colour = "black")
    }
    g1 <- g1 + theme_bw() + theme(legend.position = "right")
    .f1 <- file.path(fig_dir, "globularity_apparent_vs_expected.pdf")
    .safe_save(ggsave(.f1, g1, width = 8.5, height = 6), .f1)

    # standards check table: expected vs recovered MW (fit residual per standard)
    if (nrow(std_ok)) {
      chk <- copy(std_ok)[, .(standard = name, expected_kDa = mw_kDa, elution_fraction = fraction,
                              recovered_kDa = round(apparent, 3))]
      chk[, ratio_recovered_expected := round(recovered_kDa / expected_kDa, 3)]
      fwrite(chk, file.path(tab_dir, "globularity_standards_check.txt"), sep = "\t")
      message("[", m, "]   calibration standards (recovered vs expected MW):"); print(chk)
    }

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
    .f2 <- file.path(fig_dir, "globularity_ffo_distribution.pdf")
    .safe_save(ggsave(.f2, g2, width = 7, height = 5), .f2)

    # ---- plot 3: pie charts - headline split, and the full class breakdown ----
    .cols <- .class_colours(max_oligomer)
    hl <- data.table(group = factor(c("globular as expected", "anomalous"),
                                    levels = c("globular as expected", "anomalous")),
                     n = c(n_glob, n_test - n_glob))[n > 0]
    p_head <- .pie(hl, "group", c("globular as expected" = "#2C5F8A", "anomalous" = "#E15759"),
                   paste0("Control globularity - PCM_ctrl_vs_", m),
                   sprintf("%d proteins tested | literature expectation ~%.0f%% globular / ~%.0f%% not",
                           n_test, expected_globular_pct, 100 - expected_globular_pct))
    cls <- d[, .(n = .N), by = class]
    cls[, class := factor(class, levels = .class_levels(max_oligomer))]
    setorder(cls, class)
    p_class <- .pie(cls, "class", .cols,
                    paste0("Elution class breakdown - PCM_ctrl_vs_", m),
                    paste0("blue = compact/globular states (monomer to ", max_oligomer,
                           "x); warm = anomalous.\nSlices under 3% are left unlabelled."))
    # the DEFENSIBLE pie: restricted to proteins whose apex falls inside the calibrated MW interval,
    # where the apparent MW is an interpolation rather than an extrapolation of the standards curve
    dr    <- d[in_calibrated_range == TRUE]
    n_rr  <- nrow(dr); n_gr <- sum(dr$globular_as_expected)
    p_rng <- if (n_rr) .pie(
      data.table(group = factor(c("globular as expected", "anomalous"),
                                levels = c("globular as expected", "anomalous")),
                 n = c(n_gr, n_rr - n_gr))[n > 0],
      "group", c("globular as expected" = "#2C5F8A", "anomalous" = "#E15759"),
      paste0("WITHIN the calibrated range - PCM_ctrl_vs_", m),
      sprintf("%d of %d proteins elute between %.3g and %.3g kDa, where the standards curve is interpolated.\nOutside that interval the apparent MW is extrapolated, so this is the defensible figure.",
              n_rr, nrow(d), calibration_min_kDa, calibration_max_kDa)) else NULL

    .f3 <- file.path(fig_dir, "globularity_pies.pdf")
    .ok3 <- tryCatch({ grDevices::pdf(.f3, width = 6.5, height = 5.5)
                       print(p_head); if (!is.null(p_rng)) print(p_rng); print(p_class); grDevices::dev.off(); TRUE },
                     error = function(e) { message("   !! could not write ", basename(.f3), ": ",
                                                   conditionMessage(e), " (open in a PDF viewer? continuing)")
                                           try(grDevices::dev.off(), silent = TRUE); FALSE })
    if (.ok3) message("[", m, "]   pies -> ", .f3)
  }

  if (length(summary_rows)) {
    S <- rbindlist(summary_rows, use.names = TRUE)
    out <- here("output", out_subdir); dir.create(out, recursive = TRUE, showWarnings = FALSE)
    fwrite(S, file.path(out, "globularity_summary.csv"))   # all metabolite x condition rows

    # pooled pies across all metabolites' CONTROL sets only (the treatment rows are in the summary table
    # and in globularity_check.txt; counts are summed over metabolites, so a protein tested in several
    # control sets is counted once per set - a pooled view, not a de-duplicated one)
    S_all <- S
    if ("role" %in% names(S)) S <- S[role == "control"]
    if (!nrow(S)) S <- S_all
    .cols <- .class_colours(max_oligomer)
    hlP <- data.table(group = factor(c("globular as expected", "anomalous"),
                                     levels = c("globular as expected", "anomalous")),
                      n = c(sum(S$n_globular), sum(S$n_tested) - sum(S$n_globular)))[n > 0]
    .plev <- c("monomer", "oligomer (2-Nx)", "sub_monomer", "above_range", "between_states", "void",
               "beyond_calibration", "below_calibration_expected", "below_calibration_anomalous")
    .g <- function(nm) if (nm %in% names(S)) sum(S[[nm]]) else 0
    pooled_class <- data.table(
      class = factor(.plev, levels = .plev),
      n = c(.g("n_monomer"), .g("n_oligomer"), .g("n_sub_monomer"),
            .g("n_above_range"), .g("n_between_states"), .g("n_void"),
            .g("n_beyond_calibration"),
            .g("n_below_calibration_expected"), .g("n_below_calibration_anomalous")))[n > 0]
    pooled_cols <- c("monomer" = "#2C5F8A", "oligomer (2-Nx)" = "#8CB3D9", "sub_monomer" = "#F28E2B",
                     "above_range" = "#E15759", "between_states" = "#B07AA1", "void" = "#9C755F",
                     "beyond_calibration" = "#BAB0AC",
                     "below_calibration_expected" = "#DDD5CE", "below_calibration_anomalous" = "#8C6D62")
    # pooled, restricted to the calibrated interval - the number to quote
    hlR <- if ("n_globular_in_range" %in% names(S)) data.table(
      group = factor(c("globular as expected", "anomalous"), levels = c("globular as expected", "anomalous")),
      n = c(sum(S$n_globular_in_range), sum(S$n_in_calibrated_range) - sum(S$n_globular_in_range)))[n > 0] else NULL
    .fp <- file.path(out, "globularity_pies_pooled.pdf")
    tryCatch({
      grDevices::pdf(.fp, width = 6.5, height = 5.5)
      print(.pie(hlP, "group", c("globular as expected" = "#2C5F8A", "anomalous" = "#E15759"),
                 "Control globularity - all metabolites pooled (ALL tested)",
                 sprintf("%d protein-observations across %d control set(s) | literature expectation ~%.0f%% globular.\nIncludes proteins eluting OUTSIDE the calibrated MW interval, where the apparent MW is extrapolated - see the next pie.",
                         sum(S$n_tested), nrow(S), expected_globular_pct)))
      if (!is.null(hlR) && nrow(hlR))
        print(.pie(hlR, "group", c("globular as expected" = "#2C5F8A", "anomalous" = "#E15759"),
                   "Pooled, WITHIN the calibrated range",
                   sprintf("%d of %d protein-observations elute between %.3g and %.3g kDa, where the standards curve is interpolated.\nThis is the defensible figure to quote.",
                           sum(S$n_in_calibrated_range), sum(S$n_tested), calibration_min_kDa, calibration_max_kDa)))
      print(.pie(pooled_class, "class", pooled_cols,
                 "Elution class breakdown - all metabolites pooled",
                 "blue = compact/globular states; warm = anomalous. Slices under 3% are left unlabelled."))
      grDevices::dev.off()
    }, error = function(e) { message("!! could not write ", basename(.fp), ": ", conditionMessage(e),
                                     " (open in a PDF viewer? continuing)")
                             try(grDevices::dev.off(), silent = TRUE) })

    cat("\n==== globularity summary (one row per metabolite x condition) ====\n"); print(S_all)
    cat("Pooled pies -> ", file.path(out, "globularity_pies_pooled.pdf"), "\n", sep = "")
    cat(sprintf("\nPooled: median %.1f%% of tested proteins elute as expected for a globular species; median %.1f%% anomalous (literature expectation ~%.0f%% globular).\n",
                stats::median(S$pct_globular), stats::median(S$pct_anomalous), expected_globular_pct))
    cat("Reminder: SEC apparent MW conflates assembly and shape - 'anomalous' proteins are CANDIDATES for\n",
        "non-globular behaviour, not proof. A theoretical f/f0 reference (HYDROPRO) would separate the two.\n", sep = "")
    invisible(S)
  } else { message("Nothing computed."); invisible(NULL) }
}
