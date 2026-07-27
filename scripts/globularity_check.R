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
    "sub_monomer", "above_range", "between_states", "void", "beyond_calibration")
.class_colours <- function(max_oligomer) {
  olig <- if (max_oligomer >= 2) grDevices::colorRampPalette(c("#6B93C0", "#C6DBEF"))(max_oligomer - 1) else character(0)
  cols <- c("#2C5F8A", olig, "#F28E2B", "#E15759", "#B07AA1", "#9C755F", "#BAB0AC")
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
                              tolerance      = 1.5,   # a ratio counts as oligomer n if within this factor of n
                              max_oligomer   = 4L,    # clean oligomer states allowed as "globular as expected"
                              min_intensity  = 0,     # drop proteins below this summed ctrl intensity
                              void_fractions = 1L,    # apex here = excluded volume, MW not interpretable
                              calibration_max_kDa = NULL,   # NULL = the largest standard; above this the
                                                            # calibration is EXTRAPOLATED, not measured
                              standards      = NULL,  # data.frame(name, mw_kDa); NULL = the kit defaults
                              calibration_file = NULL,# NULL = auto-find data/raw/*calibration*.xlsx
                              expected_globular_pct = 95,   # literature expectation, for the headline only
                              out_subdir     = "globularity") {
  std <- if (is.null(standards)) .default_standards() else as.data.table(standards)
  if (!all(c("name", "mw_kDa") %in% names(std))) stop("`standards` needs columns name, mw_kDa.")
  cal_obs <- .read_calibration(calibration_file)      # observed elution fraction per standard (or NULL)
  if (is.null(calibration_max_kDa)) calibration_max_kDa <- max(std$mw_kDa, na.rm = TRUE)
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
    # BEYOND the largest calibration standard the log-linear fit is EXTRAPOLATED, so the "apparent MW"
    # there is not a measurement (it runs to physically impossible values in the earliest fractions).
    # Flag those separately instead of counting them as a quantitative "much larger than expected".
    d[, beyond_calibration := apparent_mw_kDa > calibration_max_kDa]
    d[beyond_calibration == TRUE, globular_as_expected := FALSE]

    d[, class := data.table::fcase(
      apex_fraction %in% void_fractions,                    "void",
      beyond_calibration == TRUE,                           "beyond_calibration",
      globular_as_expected & oligomer_state == 1,           "monomer",
      globular_as_expected & oligomer_state >  1,           paste0("oligomer_", oligomer_state, "x"),
      ratio < 1 / tolerance,                                "sub_monomer",
      ratio > max_oligomer * tolerance,                     "above_range",
      default =                                             "between_states")]

    # proteins whose apex falls INSIDE the calibrated range - the defensible denominator
    d[, in_calibrated_range := !(apex_fraction %in% void_fractions) & !beyond_calibration]

    # apparent frictional ratio: cube root of the MW ratio (see header)
    d[, ffo_vs_monomer := ratio^(1/3)]
    d[, ffo_vs_state   := (ratio / oligomer_state)^(1/3)]
    setorder(d, -ratio)

    tab_dir <- here("output", paste0("PCM_ctrl_vs_", m), "tables")
    fig_dir <- here("output", paste0("PCM_ctrl_vs_", m), "figures")
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(d, file.path(tab_dir, "globularity_check.txt"), sep = "\t")

    n_test <- nrow(d); n_glob <- sum(d$globular_as_expected); pct_anom <- 100 * (1 - n_glob / n_test)
    n_rng  <- sum(d$in_calibrated_range); n_glob_rng <- sum(d$globular_as_expected & d$in_calibrated_range)
    pct_anom_rng <- if (n_rng) 100 * (1 - n_glob_rng / n_rng) else NA_real_
    message(sprintf("[%s] ctrl globularity: %d/%d proteins tested (%d had no monomer mass or no usable peak).",
                    m, n_test, n_all, n_all - n_test))
    message(sprintf("[%s]   globular as expected (monomer or clean oligomer <=%dx, within %.2gx): %d (%.1f%%)",
                    m, max_oligomer, tolerance, n_glob, 100 * n_glob / n_test))
    message(sprintf("[%s]   ANOMALOUS (all tested): %d (%.1f%%)  vs ~%.0f%% globular expected from literature",
                    m, n_test - n_glob, pct_anom, expected_globular_pct))
    message(sprintf("[%s]   >> DEFENSIBLE headline - within the calibrated range (apex <= %.0f kDa, non-void): %d protein(s), ANOMALOUS %.1f%%",
                    m, calibration_max_kDa, n_rng, pct_anom_rng))
    message(sprintf("[%s]   (%d protein(s) elute beyond the largest standard / in the void: apparent MW there is EXTRAPOLATED, not measured)",
                    m, n_test - n_rng))
    print(d[, .N, by = class][order(-N)])
    message(sprintf("[%s]   median apparent f/f0 (vs monomer) = %.2f | vs assigned state = %.2f",
                    m, stats::median(d$ffo_vs_monomer, na.rm = TRUE), stats::median(d$ffo_vs_state, na.rm = TRUE)))

    summary_rows[[m]] <- data.table(
      metabolite = m, n_tested = n_test, n_globular = n_glob,
      pct_globular = round(100 * n_glob / n_test, 1), pct_anomalous = round(pct_anom, 1),
      # in-range = apex inside the calibrated MW range (the defensible denominator)
      n_in_calibrated_range = n_rng, pct_anomalous_in_range = round(pct_anom_rng, 1),
      n_beyond_calibration  = sum(d$class == "beyond_calibration"),
      n_monomer      = sum(d$class == "monomer"),
      n_oligomer     = sum(grepl("^oligomer_", d$class)),
      n_sub_monomer  = sum(d$class == "sub_monomer"),
      n_above_range  = sum(d$class == "above_range"),
      n_between_states = sum(d$class == "between_states"),
      n_void         = sum(d$class == "void"),
      median_ffo_vs_monomer = round(stats::median(d$ffo_vs_monomer, na.rm = TRUE), 3),
      median_ffo_vs_state   = round(stats::median(d$ffo_vs_state,   na.rm = TRUE), 3))

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
    ggsave(file.path(fig_dir, "globularity_apparent_vs_expected.pdf"), g1, width = 8.5, height = 6)

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
    ggsave(file.path(fig_dir, "globularity_ffo_distribution.pdf"), g2, width = 7, height = 5)

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
    grDevices::pdf(file.path(fig_dir, "globularity_pies.pdf"), width = 6.5, height = 5.5)
    print(p_head); print(p_class); grDevices::dev.off()
    message("[", m, "]   pies -> ", file.path(fig_dir, "globularity_pies.pdf"))
  }

  if (length(summary_rows)) {
    S <- rbindlist(summary_rows, use.names = TRUE)
    out <- here("output", out_subdir); dir.create(out, recursive = TRUE, showWarnings = FALSE)
    fwrite(S, file.path(out, "globularity_summary.csv"))

    # pooled pies across all metabolites' control sets (counts summed over metabolites; a protein tested
    # in several ctrl sets is counted once per set - this is a pooled view, not a de-duplicated one)
    .cols <- .class_colours(max_oligomer)
    hlP <- data.table(group = factor(c("globular as expected", "anomalous"),
                                     levels = c("globular as expected", "anomalous")),
                      n = c(sum(S$n_globular), sum(S$n_tested) - sum(S$n_globular)))[n > 0]
    .plev <- c("monomer", "oligomer (2-Nx)", "sub_monomer", "above_range", "between_states", "void", "beyond_calibration")
    pooled_class <- data.table(
      class = factor(.plev, levels = .plev),
      n = c(sum(S$n_monomer), sum(S$n_oligomer), sum(S$n_sub_monomer),
            sum(S$n_above_range), sum(S$n_between_states), sum(S$n_void),
            sum(S$n_beyond_calibration)))[n > 0]
    pooled_cols <- c("monomer" = "#2C5F8A", "oligomer (2-Nx)" = "#8CB3D9", "sub_monomer" = "#F28E2B",
                     "above_range" = "#E15759", "between_states" = "#B07AA1", "void" = "#9C755F",
                     "beyond_calibration" = "#BAB0AC")
    grDevices::pdf(file.path(out, "globularity_pies_pooled.pdf"), width = 6.5, height = 5.5)
    print(.pie(hlP, "group", c("globular as expected" = "#2C5F8A", "anomalous" = "#E15759"),
               "Control globularity - all metabolites pooled",
               sprintf("%d protein-observations across %d control set(s) | literature expectation ~%.0f%% globular",
                       sum(S$n_tested), nrow(S), expected_globular_pct)))
    print(.pie(pooled_class, "class", pooled_cols,
               "Elution class breakdown - all metabolites pooled",
               "blue = compact/globular states; warm = anomalous. Slices under 3% are left unlabelled."))
    grDevices::dev.off()

    cat("\n==== control-condition globularity summary ====\n"); print(S)
    cat("Pooled pies -> ", file.path(out, "globularity_pies_pooled.pdf"), "\n", sep = "")
    cat(sprintf("\nPooled: median %.1f%% of tested proteins elute as expected for a globular species; median %.1f%% anomalous (literature expectation ~%.0f%% globular).\n",
                stats::median(S$pct_globular), stats::median(S$pct_anomalous), expected_globular_pct))
    cat("Reminder: SEC apparent MW conflates assembly and shape - 'anomalous' proteins are CANDIDATES for\n",
        "non-globular behaviour, not proof. A theoretical f/f0 reference (HYDROPRO) would separate the two.\n", sep = "")
    invisible(S)
  } else { message("Nothing computed."); invisible(NULL) }
}
