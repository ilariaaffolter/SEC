# scripts/gradseq_ffo.R
# =============================================================================
# ORTHOGONAL hydrodynamics: join this project's SEC data to a published GLYCEROL-GRADIENT complexome
# (Grad-seq) to estimate the frictional ratio f/f0 from a second, independent physical principle - and,
# where a protein is measured in both, to solve for native mass and shape SEPARATELY.
#
# WHY THE TWO TECHNIQUES TOGETHER ARE MORE THAN THE SUM:
#   SEC            separates by Stokes radius:      R_s  ~ (f/f0) * M^(1/3)
#   sedimentation  separates by sedimentation coeff: s   ~ M^(2/3) / (f/f0)
# Each alone confounds mass with shape, but in OPPOSITE directions: taking a protein to be a monomer when
# it is an n-mer INFLATES the SEC-derived f/f0 by n^(1/3) and DEFLATES the sedimentation-derived one by
# n^(2/3). Combining them removes both unknowns (Siegel & Monty 1966):
#   M_native = 6 pi eta N_A R_s s / (1 - vbar rho)        <- shape-independent native mass
#   f/f0     = R_s / R_min(M_native),  R_min = (3 M vbar / (4 pi N_A))^(1/3)
#
# SCALES - READ THIS BEFORE COMPARING NUMBERS:
#   The SEC "apparent f/f0" written by globularity_check.R is (M_app/M_monomer)^(1/3), where M_app comes
#   from a curve fitted to GLOBULAR STANDARDS. A protein behaving exactly like those standards therefore
#   scores 1.0 BY CONSTRUCTION - it is a ratio RELATIVE to the calibrants, not an absolute frictional
#   ratio. Real compact proteins have an ABSOLUTE f/f0 of about 1.2. The sedimentation route below yields
#   an ABSOLUTE f/f0. `globular_ffo` (default 1.2) converts the SEC scale to absolute; it is stated in
#   every plot, and every comparison is made on the absolute scale.
#
# DATA REQUIRED (not bundled - download it yourself):
#   Hor et al. 2020, "Grad-seq shines light on unrecognized RNA and protein complexes in the model
#   bacterium Escherichia coli", Nucleic Acids Research 48:9301-9319 (doi:10.1093/nar/gkaa676).
#   Linear 10-40% glycerol gradient, 20 fractions + pellet, protein quantified by MS (protein coverage
#   is about 49% of the proteome - the frequently quoted ~85% refers to TRANSCRIPTS).
#   Save the supplementary protein table (xlsx/tsv/csv) somewhere and pass its path.
#
# HOW s IS OBTAINED: the gradient is calibrated from RIBOSOMAL anchors - the small subunit (30S protein
# profiles, s = 30) and the large subunit (50S, s = 50), optionally the 70S monosome - located from the
# peak fractions of rps*/rpl* proteins. A linear s ~ fraction model is fitted (rate-zonal migration in a
# linear gradient is approximately proportional to s). With only two anchors the fit is exact and cannot
# be validated; with three the residual is reported. THE CALIBRATION IS THE WEAKEST LINK - inspect the
# diagnostic plot before trusting any absolute number.
#
# WHAT IS DEFENSIBLE AND WHAT IS NOT:
#   population level  comparing the two f/f0 DISTRIBUTIONS is a fair sanity check;
#   per protein       the Siegel-Monty combination is EXPLORATORY ONLY here, because the two datasets
#                     come from different labs, buffers, growth conditions and lysis protocols, so the
#                     complexes present are not guaranteed to be the same. Treat per-protein output as
#                     hypothesis-generating, never as a measurement.
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "gradseq_ffo.R"))
#   gs <- gradseq_load("data/raw/Hor2020_gradseq_proteins.xlsx")   # inspect what was parsed
#   cal <- gradseq_calibrate(gs)                                   # LOOK at the calibration plot
#   ff  <- gradseq_ffo(gs, cal)
#   gradseq_vs_sec(ff, metabolite = "ATP")
#   gradseq_all("data/raw/Hor2020_gradseq_proteins.xlsx", metabolite = "ATP")   # all of the above
#
# OUTPUT (output/gradseq/):
#   gradseq_profiles.csv          parsed, normalised sedimentation profiles + peak fraction
#   gradseq_calibration.pdf       ribosomal anchors, the s ~ fraction fit and its residuals
#   gradseq_ribosome_check.pdf    rps* vs rpl* profiles - do the subunits separate as published?
#   gradseq_ffo.csv               per protein: peak fraction, s, absolute f/f0 (monomer assumption)
#   gradseq_ffo_distribution.pdf  the f/f0 distribution with 1.0 / 1.2 / 1.5 / 2.0 reference lines
#   sec_vs_gradseq_distribution.pdf   THE HEADLINE: both distributions on one absolute axis
#   sec_vs_gradseq_perprotein.pdf     opposite-bias diagnostic + Siegel-Monty native mass
#   sec_vs_gradseq.csv            per-protein join, native mass, n_implied, true f/f0
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })
# the banner-row detector lives in secseq_compare.R; source it so either script works standalone
if (!exists(".find_header_skip")) try(source(here::here("scripts", "secseq_compare.R")), silent = TRUE)

# physical constants (cgs, water at 20 C - the anchors are s20,w values so the reference solvent is water)
.NA_AVO <- 6.02214076e23
.ETA_W  <- 0.01002      # poise = g/(cm s)
.RHO_W  <- 0.998        # g/cm3
.VBAR   <- 0.73         # cm3/g, typical protein partial specific volume

.gs_dir <- function(...) here("output", "gradseq", ...)

# minimal-sphere radius (cm) for a mass in g/mol
.r_min <- function(M_gmol, vbar = .VBAR) (3 * M_gmol * vbar / (4 * pi * .NA_AVO))^(1/3)

# absolute f/f0 from a sedimentation coefficient (Svedberg) and a mass (g/mol)
.ffo_from_s <- function(s_svedberg, M_gmol, vbar = .VBAR, eta = .ETA_W, rho = .RHO_W) {
  s <- s_svedberg * 1e-13
  f  <- M_gmol * (1 - vbar * rho) / (.NA_AVO * s)     # friction coefficient
  f0 <- 6 * pi * eta * .r_min(M_gmol, vbar)
  f / f0
}
# Siegel & Monty: native mass (g/mol) from Stokes radius (cm) and s (Svedberg) - shape-independent
.M_native <- function(R_s_cm, s_svedberg, vbar = .VBAR, eta = .ETA_W, rho = .RHO_W)
  6 * pi * eta * .NA_AVO * R_s_cm * (s_svedberg * 1e-13) / (1 - vbar * rho)

# ---- 1. load the published Grad-seq protein table --------------------------------------------------
# Deliberately format-tolerant: supplementary tables differ in column naming between versions.
gradseq_load <- function(file, id_col = NULL, fraction_cols = NULL, sheet = 1, skip = NULL) {
  if (!file.exists(file)) { f2 <- here(file); if (file.exists(f2)) file <- f2 else stop("File not found: ", file) }
  ext <- tolower(tools::file_ext(file))
  # supplementary workbooks often carry banner rows above the real header - reuse the detector from
  # secseq_compare.R when it is loaded, else fall back to no skip
  if (is.null(skip)) {
    skip <- if (exists(".find_header_skip")) tryCatch(.find_header_skip(file, sheet), error = function(e) 0L) else 0L
    if (skip > 0) message("Detected ", skip, " banner row(s) above the header - skipping them. Override with skip = <n>.")
  }
  X <- if (ext %in% c("xlsx", "xls", "xlsm")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("install.packages('readxl') to read ", ext)
    suppressMessages(as.data.table(readxl::read_excel(file, sheet = sheet, skip = skip)))
  } else as.data.table(data.table::fread(file, skip = skip))
  if (!nrow(X)) stop("No rows read from ", file)

  # protein id column: the one richest in UniProt-like accessions, unless named
  if (is.null(id_col)) {
    acc_rx <- "^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
    score  <- vapply(names(X), function(cc) mean(grepl(acc_rx, trimws(as.character(X[[cc]]))), na.rm = TRUE), numeric(1))
    id_col <- names(X)[which.max(score)]
    if (max(score, na.rm = TRUE) < 0.3)
      stop("Could not find a UniProt-accession column automatically. Pass id_col = '<name>'. Columns: ",
           paste(names(X), collapse = ", "))
    message("Using '", id_col, "' as the protein id (", round(100 * max(score)), "% accession-like).")
  }
  # fraction columns: numeric columns whose names carry a fraction number
  if (is.null(fraction_cols)) {
    cand <- names(X)[grepl("^([Ff]raction[ _.]*)?[0-9]{1,2}$|^F[0-9]{1,2}$|fraction", names(X), ignore.case = TRUE)]
    cand <- cand[vapply(cand, function(cc) is.numeric(X[[cc]]) || !anyNA(suppressWarnings(as.numeric(as.character(X[[cc]])))), logical(1))]
    fraction_cols <- cand
    if (length(fraction_cols) < 5)
      stop("Found only ", length(fraction_cols), " fraction column(s). Pass fraction_cols = c(...). Columns: ",
           paste(names(X), collapse = ", "))
    message("Using ", length(fraction_cols), " fraction column(s): ", paste(utils::head(fraction_cols, 25), collapse = ", "))
  }
  num <- function(v) suppressWarnings(as.numeric(as.character(v)))
  M <- as.matrix(as.data.frame(lapply(X[, ..fraction_cols], num)))
  M[!is.finite(M)] <- 0
  # fraction number parsed from the column name; the pellet (if present) is dropped for peak calling
  fno <- suppressWarnings(as.numeric(gsub("[^0-9]", "", fraction_cols)))
  ok  <- is.finite(fno)
  if (any(!ok)) message("Ignoring ", sum(!ok), " column(s) without a parsable fraction number (pellet?).")
  M <- M[, ok, drop = FALSE]; fno <- fno[ok]
  o <- order(fno); M <- M[, o, drop = FALSE]; fno <- fno[o]

  ids <- trimws(as.character(X[[id_col]]))
  # a row may list several accessions - keep the first
  ids <- sub("[;,].*$", "", ids)
  keep <- nzchar(ids) & rowSums(M) > 0
  M <- M[keep, , drop = FALSE]; ids <- ids[keep]
  rownames(M) <- make.unique(ids)

  Mn <- M / rowSums(M)                              # normalise each protein to its own profile
  peak <- fno[max.col(Mn, ties.method = "first")]
  com  <- as.vector(Mn %*% fno)                     # centre of mass, a smoother position estimate
  D <- data.table(protein_id = ids, peak_fraction = peak, com_fraction = com, total = rowSums(M))
  dir.create(.gs_dir(), recursive = TRUE, showWarnings = FALSE)
  fwrite(cbind(D, as.data.table(Mn)), .gs_dir("gradseq_profiles.csv"))
  message("Parsed ", nrow(D), " protein profile(s) over fractions ", min(fno), "-", max(fno), ".")
  invisible(list(profiles = Mn, fractions = fno, meta = D, id_col = id_col))
}

# ---- 2. calibrate fraction -> sedimentation coefficient --------------------------------------------
# Anchors: the ribosomal subunits. Their peak fractions are found from the median peak of the rps*/rpl*
# protein groups, which is far more robust than any single protein.
gradseq_calibrate <- function(gs, anchors = NULL, gene_map = NULL, save_plots = TRUE) {
  D <- gs$meta
  # protein -> gene symbol, from the shared UniProt cache unless supplied
  if (is.null(gene_map)) {
    sf <- here("output", "uniprot_annotation_shared.RData")
    if (file.exists(sf)) { e <- new.env(); load(sf, envir = e)
      u <- as.data.table(e$.uniprot_all)
      if (all(c("input_id", "gene_names") %in% names(u)))
        gene_map <- setNames(tolower(sub(" .*$", "", as.character(u$gene_names))), as.character(u$input_id)) }
  }
  if (is.null(gene_map)) stop("No gene symbols available - render a comparison first (UniProt cache), or pass gene_map.")
  D[, gene := unname(gene_map[protein_id])]

  small <- D[grepl("^rps[a-z]$", gene)]      # 30S proteins
  large <- D[grepl("^rpl[a-z]$", gene)]      # 50S proteins
  message("Ribosomal anchors found: ", nrow(small), " small-subunit (rps*), ", nrow(large), " large-subunit (rpl*) protein(s).")
  if (nrow(small) < 3 || nrow(large) < 3)
    stop("Too few ribosomal proteins matched to calibrate. Check the gene mapping, or pass explicit anchors.")

  if (is.null(anchors)) anchors <- data.table(
    name = c("30S", "50S"),
    s    = c(30, 50),
    fraction = c(stats::median(small$peak_fraction), stats::median(large$peak_fraction)))
  else anchors <- as.data.table(anchors)
  if (nrow(anchors) < 2) stop("Need at least two anchors.")
  print(anchors)
  if (anchors[name == "50S"]$fraction <= anchors[name == "30S"]$fraction)
    warning("The 50S anchor does not sediment further than the 30S anchor - the gradient orientation or the ",
            "fraction numbering may be reversed. Inspect gradseq_ribosome_check.pdf before continuing.")

  fit <- stats::lm(s ~ fraction, data = anchors)
  r2  <- if (nrow(anchors) > 2) summary(fit)$r.squared else NA_real_
  message(sprintf("Calibration: s = %.3f + %.3f * fraction%s",
                  stats::coef(fit)[1], stats::coef(fit)[2],
                  if (is.finite(r2)) sprintf("  (R2 = %.3f over %d anchors)", r2, nrow(anchors))
                  else "  (2 anchors: exact fit, NOT validatable - add a 70S anchor to test linearity)"))
  s_at_1 <- unname(stats::predict(fit, data.frame(fraction = min(gs$fractions))))
  if (s_at_1 < -5) warning(sprintf("The fit implies s = %.1f at the top fraction, which is unphysical: the linear model is poor. Treat absolute s values with caution.", s_at_1))

  if (save_plots) {
    dir.create(.gs_dir(), recursive = TRUE, showWarnings = FALSE)
    gp <- ggplot(anchors, aes(fraction, s)) +
      geom_abline(slope = stats::coef(fit)[2], intercept = stats::coef(fit)[1], colour = "steelblue") +
      geom_point(size = 3) + geom_text(aes(label = name), vjust = -1) +
      labs(title = "Gradient calibration: sedimentation coefficient vs fraction",
           subtitle = paste0("Anchors are the ribosomal subunit peaks (median over rps*/rpl* proteins).\n",
                             "Rate-zonal migration in a linear gradient is ~proportional to s; with only two anchors this fit\n",
                             "is exact and cannot be validated. THIS IS THE WEAKEST STEP - all absolute s and f/f0 depend on it."),
           x = "fraction", y = "s (Svedberg)") + theme_bw()
    tryCatch(ggsave(.gs_dir("gradseq_calibration.pdf"), gp, width = 7, height = 5), error = function(e) NULL)

    L <- rbindlist(list(cbind(small, subunit = "30S (rps*)"), cbind(large, subunit = "50S (rpl*)")))
    gr <- ggplot(L, aes(peak_fraction, fill = subunit)) +
      geom_histogram(binwidth = 1, position = "identity", alpha = 0.6) +
      geom_vline(data = anchors, aes(xintercept = fraction), linetype = 2) +
      labs(title = "Ribosomal subunit peak fractions (sanity check)",
           subtitle = "The two subunits should peak in clearly distinct fractions, as in the published A260 profile.",
           x = "peak fraction", y = "proteins", fill = NULL) + theme_bw()
    tryCatch(ggsave(.gs_dir("gradseq_ribosome_check.pdf"), gr, width = 7, height = 5), error = function(e) NULL)
  }
  invisible(list(fit = fit, anchors = anchors, gene_map = gene_map))
}

# ---- 3. per-protein absolute f/f0 from sedimentation -----------------------------------------------
gradseq_ffo <- function(gs, cal, position = c("com", "peak"), mass_map = NULL, vbar = .VBAR, save_plots = TRUE) {
  position <- match.arg(position)
  D <- copy(gs$meta)
  if (is.null(mass_map)) {
    sf <- here("output", "uniprot_annotation_shared.RData")
    if (!file.exists(sf)) stop("No UniProt cache for monomer masses - render a comparison first, or pass mass_map.")
    e <- new.env(); load(sf, envir = e); u <- as.data.table(e$.uniprot_all)
    if (!all(c("input_id", "mass") %in% names(u))) stop("UniProt cache has no input_id/mass.")
    mass_map <- setNames(as.numeric(u$mass), as.character(u$input_id))
  }
  D[, mw_Da := unname(mass_map[protein_id])]
  D[, frac_used := if (position == "com") com_fraction else as.numeric(peak_fraction)]
  D[, s_svedberg := unname(stats::predict(cal$fit, data.frame(fraction = frac_used)))]
  D <- D[is.finite(mw_Da) & mw_Da > 0 & is.finite(s_svedberg) & s_svedberg > 0]
  if (!nrow(D)) stop("No protein has both a monomer mass and a positive s - check the calibration.")
  D[, ffo_gradseq := .ffo_from_s(s_svedberg, mw_Da, vbar)]
  setorder(D, -ffo_gradseq)
  fwrite(D, .gs_dir("gradseq_ffo.csv"))
  message(sprintf("Sedimentation-derived ABSOLUTE f/f0 for %d protein(s): median %.2f (IQR %.2f-%.2f).",
                  nrow(D), stats::median(D$ffo_gradseq, na.rm = TRUE),
                  stats::quantile(D$ffo_gradseq, .25, na.rm = TRUE), stats::quantile(D$ffo_gradseq, .75, na.rm = TRUE)))
  message("   Reference: ~1.2 = compact globular (hydrated), ~1.5 = moderately elongated, >2 = extended/disordered.")
  message("   NOTE the monomer assumption DEFLATES this estimate for oligomers (by n^(2/3)), so a high value is a strong claim.")
  if (save_plots) {
    refs <- data.table(x = c(1.0, 1.2, 1.5, 2.0), lab = c("sphere", "globular", "elongated", "extended"))
    g <- ggplot(D[is.finite(ffo_gradseq) & ffo_gradseq > 0], aes(ffo_gradseq)) +
      geom_histogram(bins = 80, fill = "darkorange", colour = "white") +
      geom_vline(data = refs, aes(xintercept = x), linetype = 2, colour = "grey35") +
      geom_text(data = refs, aes(x = x, y = Inf, label = lab), vjust = 1.3, size = 3, colour = "grey25") +
      scale_x_log10() +
      labs(title = "Absolute f/f0 from the glycerol gradient (monomer assumption)",
           subtitle = "f/f0 = M(1-vbar*rho) / (N_A * s * 6*pi*eta*R_min). Oligomers are UNDER-estimated here by n^(2/3),\nso this is a conservative view of how extended the proteome is.",
           x = "absolute f/f0", y = "proteins") + theme_bw()
    tryCatch(ggsave(.gs_dir("gradseq_ffo_distribution.pdf"), g, width = 7, height = 5), error = function(e) NULL)
  }
  invisible(D)
}

# ---- 4. join to this project's SEC data ------------------------------------------------------------
gradseq_vs_sec <- function(ff, metabolite, condition = NULL, globular_ffo = 1.2, vbar = .VBAR, save_plots = TRUE) {
  gf <- here("output", paste0("PCM_ctrl_vs_", metabolite), "tables", "globularity_check.txt")
  if (!file.exists(gf)) stop("No globularity_check.txt for ", metabolite, " - run globularity_check() first.")
  G <- fread(gf)
  if ("condition" %in% names(G)) {
    cn <- unique(as.character(G$condition))
    cc <- if (!is.null(condition)) condition else { x <- cn[grepl("ctrl|control|ref", cn, ignore.case = TRUE)][1]; if (is.na(x)) cn[1] else x }
    G <- G[condition == cc]; message("Using the '", cc, "' rows of the SEC table.")
  }
  if (!all(c("protein_id", "apparent_mw_kDa", "ffo_vs_monomer") %in% names(G)))
    stop("globularity_check.txt lacks the expected columns.")

  # SEC apparent MW -> Stokes radius. The calibration reports the mass of an equivalent GLOBULAR protein,
  # so R_s = R_min(M_app) * globular_ffo, where globular_ffo is the absolute f/f0 of the calibrants
  # (~1.2). This is exactly the relative->absolute conversion described in the header.
  G[, R_s_cm := .r_min(apparent_mw_kDa * 1000, vbar) * globular_ffo]
  G[, ffo_sec_absolute := ffo_vs_monomer * globular_ffo]

  J <- merge(G[, .(protein_id, apparent_mw_kDa, expected_mw_kDa, ffo_vs_monomer, ffo_sec_absolute, R_s_cm)],
             ff[, .(protein_id, s_svedberg, mw_Da, ffo_gradseq, frac_used)], by = "protein_id")
  if (!nrow(J)) stop("No protein is present in both datasets - check the accession formats.")
  message("Proteins measured in BOTH datasets: ", nrow(J), " (SEC ", nrow(G), ", Grad-seq ", nrow(ff), ").")

  # Siegel & Monty: mass without a shape assumption, then the true f/f0 and the implied oligomeric state
  J[, M_native_Da := .M_native(R_s_cm, s_svedberg, vbar)]
  J[, n_implied   := M_native_Da / mw_Da]
  J[, ffo_true    := R_s_cm / .r_min(M_native_Da, vbar)]
  fwrite(J, .gs_dir("sec_vs_gradseq.csv"))

  msec <- stats::median(J$ffo_sec_absolute, na.rm = TRUE)
  mgrd <- stats::median(J$ffo_gradseq,      na.rm = TRUE)
  mtru <- stats::median(J$ffo_true,         na.rm = TRUE)
  message(sprintf("Median absolute f/f0 - SEC: %.2f | Grad-seq: %.2f | Siegel-Monty combined: %.2f",
                  msec, mgrd, mtru))
  message(sprintf("Median implied oligomeric state (M_native / M_monomer): %.2f", stats::median(J$n_implied, na.rm = TRUE)))
  message("Per-protein values are EXPLORATORY: different lab, buffer, growth condition and lysis, so the ",
          "complexes present need not match.")

  if (save_plots) {
    L <- rbindlist(list(
      data.table(ffo = J$ffo_sec_absolute, src = paste0("SEC (this study, x", globular_ffo, ")")),
      data.table(ffo = J$ffo_gradseq,      src = "Grad-seq (sedimentation)"),
      data.table(ffo = J$ffo_true,         src = "combined (Siegel-Monty)")))[is.finite(ffo) & ffo > 0]
    refs <- data.table(x = c(1.2, 1.5, 2.0), lab = c("globular", "elongated", "extended"))
    g1 <- ggplot(L, aes(ffo, fill = src)) +
      geom_density(alpha = 0.4, colour = NA) +
      geom_vline(data = refs, aes(xintercept = x), linetype = 2, colour = "grey35") +
      geom_text(data = refs, aes(x = x, y = Inf, label = lab), vjust = 1.4, size = 3, colour = "grey25", inherit.aes = FALSE) +
      scale_x_log10() +
      labs(title = paste0("How compact is the proteome? Two independent techniques (", metabolite, " control)"),
           subtitle = paste0("All on the ABSOLUTE scale. The SEC curve is the relative f/f0 multiplied by ", globular_ffo,
                             " (the absolute f/f0 of the globular calibrants).\nOligomers inflate the SEC estimate (n^1/3) and deflate the sedimentation one (n^2/3); the combined curve removes both."),
           x = "absolute f/f0", y = "density", fill = NULL) +
      theme_bw() + theme(legend.position = "top")
    tryCatch(ggsave(.gs_dir("sec_vs_gradseq_distribution.pdf"), g1, width = 8, height = 5.5), error = function(e) NULL)

    P <- J[is.finite(ffo_sec_absolute) & is.finite(ffo_gradseq) & ffo_sec_absolute > 0 & ffo_gradseq > 0]
    rho <- if (nrow(P) > 10) suppressWarnings(stats::cor(log(P$ffo_sec_absolute), log(P$ffo_gradseq), method = "spearman")) else NA_real_
    g2 <- ggplot(P, aes(ffo_gradseq, ffo_sec_absolute, colour = pmin(pmax(n_implied, 0.5), 8))) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
      geom_point(alpha = 0.5, size = 1) +
      scale_x_log10() + scale_y_log10() + scale_colour_viridis_c(option = "C", trans = "log10", name = "n implied") +
      labs(title = "Opposite-bias diagnostic: SEC vs sedimentation f/f0",
           subtitle = sprintf("Agreement on the dashed line means the monomer assumption holds. Systematic offsets indicate assembly:\noligomers sit ABOVE the line (SEC inflated by n^1/3, sedimentation deflated by n^2/3). Spearman rho(log) = %.2f.", rho),
           x = "f/f0 from sedimentation", y = "f/f0 from SEC (absolute)") + theme_bw()
    g3 <- ggplot(J[is.finite(M_native_Da) & M_native_Da > 0], aes(mw_Da / 1000, M_native_Da / 1000)) +
      geom_abline(slope = 1, intercept = log10(c(1, 2, 4)), linetype = 2, colour = "grey60") +
      geom_point(alpha = 0.45, size = 0.9, colour = "steelblue") +
      scale_x_log10() + scale_y_log10() +
      labs(title = "Native mass without a shape assumption (Siegel & Monty)",
           subtitle = "M_native = 6*pi*eta*N_A*R_s*s/(1-vbar*rho). Dashed = 1x, 2x, 4x the monomer mass.\nEXPLORATORY: the two measurements come from different labs and conditions.",
           x = "monomer mass (kDa, UniProt)", y = "native mass (kDa, combined)") + theme_bw()
    tryCatch({ grDevices::pdf(.gs_dir("sec_vs_gradseq_perprotein.pdf"), width = 7.5, height = 5.5)
               print(g2); print(g3); grDevices::dev.off() },
             error = function(e) try(grDevices::dev.off(), silent = TRUE))
  }
  invisible(J)
}

gradseq_all <- function(file, metabolite, id_col = NULL, fraction_cols = NULL, globular_ffo = 1.2) {
  gs  <- gradseq_load(file, id_col = id_col, fraction_cols = fraction_cols)
  cal <- gradseq_calibrate(gs)
  ff  <- gradseq_ffo(gs, cal)
  gradseq_vs_sec(ff, metabolite = metabolite, globular_ffo = globular_ffo)
}
