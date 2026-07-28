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
# HOW s IS OBTAINED: the gradient is calibrated from RIBOSOMAL anchors - the small subunit (30S, s = 30),
# the large subunit (50S, s = 50) and the 70S monosome (s = 70), located from the profiles of rps*/rpl*
# proteins. Because s is NOT additive (30 + 50 = 80, but the monosome sediments at 70), the 70S anchor is
# genuinely independent, and it is the first anchor that lets the model be tested at all: each anchor is
# predicted from the others (leave-one-out) and the errors are reported.
# The curve is forced through s = 0 at the LOAD ZONE, because a particle that does not sediment does not
# move. Leaving that constraint out - the old behaviour, still available as model = "free_linear" - gives
# a large negative intercept, over-estimates s for ordinary proteins, and was the direct cause of the
# impossible f/f0 < 1 seen earlier.
# THE CALIBRATION REMAINS THE WEAKEST LINK: every anchor is at s = 30-70 while ordinary proteins are at
# s = 2-10, so all protein values are extrapolations. Inspect the diagnostic plot first.
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
#   gradseq_selftest()                                             # verify the hydrodynamics first
#   gs <- gradseq_load("data/raw/Hor2020_gradseq_proteins.xlsx")   # inspect what was parsed
#   cal <- gradseq_calibrate(gs)                                   # LOOK at the calibration plot
#   # 70S is detected automatically; override or disable it, and switch the model, with:
#   #   gradseq_calibrate(gs, anchor_70S = 18.5)      known fraction (e.g. from the published A260 trace)
#   #   gradseq_calibrate(gs, anchor_70S = "none")    two anchors only
#   #   gradseq_calibrate(gs, model = "power")        or "free_linear" (the old, unconstrained fit)
#   gradseq_compare_models(gs, cal)                                # which model does the proteome reject?
#   # ribosomal anchors alone do NOT calibrate s for ordinary proteins (see gradseq_ffo's physical check);
#   # when that fails, use the calibration-free comparison instead:
#   gradseq_vs_sec_deviation(gs, metabolite = "ATP")
#   ff  <- gradseq_ffo(gs, cal)
#   gradseq_vs_sec(ff, metabolite = "ATP")
#   gradseq_all("data/raw/Hor2020_gradseq_proteins.xlsx", metabolite = "ATP")   # all of the above
#
# OUTPUT (output/gradseq/):
#   gradseq_profiles.csv          parsed, normalised sedimentation profiles + peak fraction
#   gradseq_calibration.pdf       all three fraction -> s models, the anchors, the load zone and where
#                                 compact proteins of known mass should elute under the model in use
#   gradseq_ribosome_check.pdf    rps* vs rpl* peak fractions, plus the co-migration trace used to place 70S
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

# inverse of .ffo_from_s: the s (Svedberg) a particle of mass M and frictional ratio f/f0 would show
.s_from_ffo <- function(M_gmol, ffo, vbar = .VBAR, eta = .ETA_W, rho = .RHO_W) {
  f0 <- 6 * pi * eta * .r_min(M_gmol, vbar)
  (M_gmol * (1 - vbar * rho) / (.NA_AVO * ffo * f0)) / 1e-13
}

# ---- 0. SELF-TEST of the hydrodynamics -------------------------------------------------------------
# Every absolute number in this file rests on four equations. They are checked here against proteins
# whose s20,w AND Stokes radius were BOTH measured independently, which makes the test non-circular:
#   f/f0 computed from s     must equal    f/f0 computed from R_s
# and the Siegel-Monty combination of the two must return the known mass. If these agree, the algebra
# and the units are right; if they do not, nothing downstream can be trusted.
# The values are the classical gel-filtration / analytical-ultracentrifugation calibration set (the
# proteins on every SEC standards vial). CHECK THEM against your own source before quoting them - they
# are transcribed here, not measured here, and different tables differ by a few per cent.
.HYDRO_STANDARDS <- data.table::data.table(
  name   = c("ribonuclease A", "chymotrypsinogen A", "ovalbumin", "BSA",
             "aldolase", "catalase", "ferritin", "thyroglobulin"),
  mw_kDa = c(13.7,  25.0,  43.5,  66.0,  158.0, 232.0, 440.0, 669.0),
  s20w   = c(1.78,  2.54,  3.55,  4.31,  7.35,  11.30, 17.60, 19.40),
  Rs_nm  = c(1.64,  2.09,  3.05,  3.55,  4.81,  5.22,  6.10,  8.50))

gradseq_selftest <- function(vbar = .VBAR, verbose = TRUE) {
  S <- data.table::copy(.HYDRO_STANDARDS)
  M <- S$mw_kDa * 1000
  S[, r_min_nm    := .r_min(M, vbar) * 1e7]
  S[, ffo_from_Rs := Rs_nm / r_min_nm]                       # geometry route
  S[, ffo_from_s  := .ffo_from_s(s20w, M, vbar)]             # hydrodynamic route
  S[, disagree_pct := 100 * abs(ffo_from_s - ffo_from_Rs) / ffo_from_Rs]
  S[, M_siegel_monty_kDa := .M_native(Rs_nm * 1e-7, s20w, vbar) / 1000]
  S[, mass_error_pct := 100 * (M_siegel_monty_kDa - mw_kDa) / mw_kDa]
  S[, s_roundtrip := .s_from_ffo(M, ffo_from_s, vbar)]       # must return s20w exactly
  if (verbose) {
    message("HYDRODYNAMIC SELF-TEST - the two independent routes to f/f0, and the mass they jointly imply:")
    print(S[, .(name, mw_kDa, s20w, Rs_nm,
                ffo_from_Rs = round(ffo_from_Rs, 3), ffo_from_s = round(ffo_from_s, 3),
                disagree_pct = round(disagree_pct, 1),
                M_SM_kDa = round(M_siegel_monty_kDa, 1), mass_err_pct = round(mass_error_pct, 1))])
    message(sprintf("Round-trip s -> f/f0 -> s: max error %.2g%% (this is pure algebra; anything above 1e-8%% is a coding bug).",
                    100 * max(abs(S$s_roundtrip - S$s20w) / S$s20w)))
    message(sprintf("Agreement of the two f/f0 routes: median %.1f%%, worst %.1f%% (%s).",
                    stats::median(S$disagree_pct), max(S$disagree_pct), S$name[which.max(S$disagree_pct)]))
    message(sprintf("Siegel-Monty mass recovery: median error %+.1f%%, worst %+.1f%%.",
                    stats::median(S$mass_error_pct), S$mass_error_pct[which.max(abs(S$mass_error_pct))]))
    message(sprintf("Empirical f/f0 of these globular standards: median %.2f (range %.2f-%.2f).\n   -> `globular_ffo` is the number used to put the SEC scale on the absolute one; the default 1.2 is at the LOW end of this set.",
                    stats::median(S$ffo_from_s), min(S$ffo_from_s), max(S$ffo_from_s)))
    message("   The two smallest standards disagree most: at that size the hydration shell is a large fraction of the particle,\n   and tabulated R_s for small proteins tends to run low. Above ~40 kDa the routes agree to a few per cent.")
  }
  invisible(S)
}

# ---- 1. load the published Grad-seq protein table --------------------------------------------------
# Deliberately format-tolerant: supplementary tables differ in column naming between versions.
# List the sheets of a workbook - a macro workbook (.xlsm) usually holds several, and the protein table
# is rarely the first. Run this before gradseq_load() if the parse looks wrong.
gradseq_sheets <- function(file) {
  if (!file.exists(file)) { f2 <- here(file); if (file.exists(f2)) file <- f2 else stop("File not found: ", file) }
  if (!requireNamespace("readxl", quietly = TRUE)) stop("install.packages('readxl')")
  s <- readxl::excel_sheets(file)
  message("Sheets in ", basename(file), ":"); print(data.frame(index = seq_along(s), sheet = s))
  invisible(s)
}

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
# Anchors: the ribosomal particles. Their positions come from the median peak of the rps*/rpl* protein
# groups, which is far more robust than any single protein.
#
# THE 70S MONOSOME IS A GENUINELY INDEPENDENT THIRD ANCHOR, not a redundant one: sedimentation
# coefficients are NOT additive (30 + 50 = 80, but the assembled particle sediments at 70), because
# joining the subunits buries surface and the friction per unit mass drops. So 70S carries information
# the two subunits do not, and it is the first anchor that lets the fraction -> s model be TESTED.
# It is located as the fraction where rps* AND rpl* proteins co-migrate most strongly, below the 50S
# peak - free 30S has no rpl*, free 50S has no rps*, only the monosome has both.
#
# THE MODEL MATTERS MORE THAN THE ANCHORS. Rate-zonal migration starts at the load zone: a particle with
# s = 0 does not move, so the curve MUST pass through s = 0 at the top of the gradient. An unconstrained
# straight line through two anchors at s = 30 and 50 ignores that, and comes out with a large negative
# intercept - which is precisely what produced impossible f/f0 < 1 for ordinary proteins. Fitting
# through the load zone is both physically required and far better behaved on extrapolation.
#   model = "zero_anchored"  s = b * (fraction - load_fraction)          standard rate-zonal treatment
#   model = "power"          s = c * (fraction - load_fraction)^p        allows the gradient to compress
#   model = "free_linear"    s = a + b * fraction                        the old, unconstrained fit
#   model = "auto"           <- DEFAULT: chosen by leave-one-out when there are three anchors,
#                               otherwise "zero_anchored" (the only one carrying a physical constraint)
# With three anchors each model is validated by leave-one-out and the errors are reported; because each
# such fit uses only two points, confirm the choice with gradseq_compare_models(), which asks the harder
# question - how much of the proteome does this model push below the physical floor f/f0 = 1? Even so, all
# anchors sit at s = 30-70 while ordinary proteins are at s = 2-10, so every protein value is an
# EXTRAPOLATION below the calibrated range. That limitation is not removed by adding 70S; it is only
# measured. gradseq_vs_sec_deviation() avoids it entirely.
gradseq_calibrate <- function(gs, anchors = NULL, gene_map = NULL,
                              anchor_70S = "auto", s_70S = 70, extra_anchors = NULL,
                              load_fraction = NULL,
                              model = c("auto", "zero_anchored", "power", "free_linear"),
                              save_plots = TRUE) {
  model <- match.arg(model)
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

  # row indices, because gs$profiles and gs$meta share row order (profile rownames are made unique)
  i_small <- which(grepl("^rps[a-z]$", D$gene))
  i_large <- which(grepl("^rpl[a-z]$", D$gene))
  small <- D[i_small]; large <- D[i_large]
  message("Ribosomal anchors found: ", nrow(small), " small-subunit (rps*), ", nrow(large), " large-subunit (rpl*) protein(s).")
  if (nrow(small) < 3 || nrow(large) < 3)
    stop("Too few ribosomal proteins matched to calibrate. Check the gene mapping, or pass explicit anchors.")

  fr  <- gs$fractions
  f30 <- stats::median(small$peak_fraction)
  f50 <- stats::median(large$peak_fraction)
  # mean normalised profile of each subunit group, each scaled to its own maximum
  ps <- colMeans(gs$profiles[i_small, , drop = FALSE]); ps <- ps / max(ps)
  pl <- colMeans(gs$profiles[i_large, , drop = FALSE]); pl <- pl / max(pl)
  co <- sqrt(pmax(ps, 0) * pmax(pl, 0))    # co-migration: high ONLY where both subunits are present
  # independent cross-check of the two subunit anchors: free 30S is where rps* exceeds rpl*, and vice versa
  message(sprintf("Subunit anchors: 30S at fraction %g, 50S at %g (median peak). Cross-check from the differential profiles: %g and %g.",
                  f30, f50, fr[which.max(ps - pl)], fr[which.max(pl - ps)]))

  # ---- locate the 70S monosome -----------------------------------------------------------------
  f70 <- NA_real_; how70 <- ""
  if (is.numeric(anchor_70S)) { f70 <- as.numeric(anchor_70S)[1]; how70 <- "supplied by the user" }
  else if (identical(anchor_70S, "auto")) {
    cand <- which(fr > f50 + 0.5)                       # the monosome sediments FURTHER than the 50S
    isloc <- if (length(cand)) vapply(cand, function(j) {
      lo <- if (j > 1L) co[j - 1L] else -Inf
      hi <- if (j < length(co)) co[j + 1L] else -Inf
      is.finite(co[j]) && co[j] >= lo && co[j] >= hi }, logical(1)) else logical(0)
    pick <- cand[isloc]
    pick <- if (length(pick)) pick[which.max(co[pick])] else integer(0)
    if (length(pick) == 1L && co[pick] >= 0.2 * max(co, na.rm = TRUE)) {
      w <- which(abs(fr - fr[pick]) <= 1.5)             # sub-fraction refinement: local centroid
      f70 <- sum(fr[w] * co[w]) / sum(co[w])
      how70 <- sprintf("detected automatically (co-migration peak at fraction %g, centroid %.2f, score %.2f of max)",
                       fr[pick], f70, co[pick] / max(co, na.rm = TRUE))
    } else {
      message("70S NOT detected: no clear rps*/rpl* co-migration peak below the 50S anchor. ",
              "Either the monosome was dissociated, or the subunits are not resolved. ",
              "Pass anchor_70S = <fraction> if you can read it off the published A260 profile.")
    }
  }
  if (is.finite(f70)) message("70S monosome anchor: fraction ", round(f70, 2), " - ", how70, ".")

  # ---- assemble and sanity-check the anchor set -------------------------------------------------
  if (is.null(anchors)) {
    anchors <- data.table(name = c("30S", "50S"), s = c(30, 50), fraction = c(f30, f50))
    if (is.finite(f70)) anchors <- rbind(anchors, data.table(name = "70S", s = s_70S, fraction = f70))
  } else anchors <- as.data.table(anchors)
  if (!is.null(extra_anchors)) anchors <- rbind(anchors, as.data.table(extra_anchors), fill = TRUE)
  anchors <- anchors[is.finite(s) & is.finite(fraction) & s > 0]
  if (nrow(anchors) < 2) stop("Need at least two usable anchors.")
  setorder(anchors, s)

  if (all(c("30S", "50S") %in% anchors$name) && anchors[name == "50S"]$fraction <= anchors[name == "30S"]$fraction)
    warning("The 50S anchor does not sediment further than the 30S anchor - the gradient orientation or the ",
            "fraction numbering may be reversed. Inspect gradseq_ribosome_check.pdf before continuing.", call. = FALSE)
  if ("70S" %in% anchors$name) {
    if (anchors[name == "70S"]$fraction <= anchors[name == "50S"]$fraction)
      warning("The 70S anchor does not sediment further than the 50S - it cannot be the monosome. Drop it with anchor_70S = 'none'.", call. = FALSE)
    # if the 'free 50S' median peak sits on top of the monosome, that anchor is mislabelled
    if (abs(anchors[name == "50S"]$fraction - anchors[name == "70S"]$fraction) < 1)
      warning("The 50S median peak coincides with the detected 70S position: most rpl* protein is probably in the ",
              "MONOSOME, so the '50S' anchor is really 70S. Set the 50S fraction explicitly via `anchors`.", call. = FALSE)
  }

  # ---- the physical zero: a particle with s = 0 stays in the load zone ---------------------------
  if (is.null(load_fraction)) load_fraction <- min(fr) - 0.5
  anchors[, d := fraction - load_fraction]
  if (any(anchors$d <= 0)) stop("An anchor sits at or above the load zone (fraction <= ", load_fraction,
                                ") - check load_fraction.")
  message(sprintf("Load zone taken as fraction %.2f (top boundary of the first collected fraction); migration distance d = fraction - %.2f.",
                  load_fraction, load_fraction))
  message("   This choice matters: it is the physical s = 0 point, and proteins peaking near the top are most sensitive to it. Override with load_fraction = <n>.")

  fit_zero <- stats::lm(s ~ 0 + d, data = anchors)
  fit_free <- stats::lm(s ~ fraction, data = anchors)
  fit_pow  <- stats::lm(log(s) ~ log(d), data = anchors)
  bz <- unname(stats::coef(fit_zero)[1]); cp <- stats::coef(fit_pow); cf <- stats::coef(fit_free)
  message(sprintf("Candidate calibrations over %d anchor(s):", nrow(anchors)))
  message(sprintf("   zero_anchored : s = %.3f * d                       (forced through s = 0 at the load zone)", bz))
  message(sprintf("   power         : s = %.3f * d^%.3f                 (exponent %s 1 => the gradient %s with depth)",
                  exp(cp[1]), cp[2], if (cp[2] > 1) ">" else "<",
                  if (cp[2] > 1) "compresses" else "expands"))
  message(sprintf("   free_linear   : s = %.3f + %.3f * fraction        (unconstrained; s = %.1f at the load zone)",
                  cf[1], cf[2], cf[1] + cf[2] * load_fraction))
  if (cf[1] + cf[2] * load_fraction < -2)
    message("   -> free_linear puts a strongly NEGATIVE s at the load zone, which is unphysical and over-estimates s for ",
            "ordinary proteins (hence f/f0 < 1). This is why it is no longer the default.")

  # ---- leave-one-out validation: the whole point of having a third anchor ------------------------
  loo <- NULL
  if (nrow(anchors) >= 3) {
    loo <- rbindlist(lapply(seq_len(nrow(anchors)), function(i) {
      A <- anchors[-i]; B <- anchors[i]
      pz <- unname(stats::coef(stats::lm(s ~ 0 + d, data = A))[1]) * B$d
      k  <- stats::coef(stats::lm(log(s) ~ log(d), data = A)); pp <- exp(k[1]) * B$d^k[2]
      pf <- unname(stats::predict(stats::lm(s ~ fraction, data = A), B))
      data.table(anchor = B$name, s_true = B$s,
                 zero_anchored = pz, power = unname(pp), free_linear = pf)
    }))
    err <- loo[, lapply(.SD, function(p) 100 * mean(abs(p - s_true) / s_true)),
               .SDcols = c("zero_anchored", "power", "free_linear")]
    message("Leave-one-out validation (predict each anchor from the others) - mean |error|:")
    print(loo[, .(anchor, s_true, zero_anchored = round(zero_anchored, 1),
                  power = round(power, 1), free_linear = round(free_linear, 1))])
    best <- names(err)[which.min(unlist(err))]
    message(sprintf("   zero_anchored %.1f%% | power %.1f%% | free_linear %.1f%%  -> best generalisation: %s",
                    err$zero_anchored, err$power, err$free_linear, best))
    if (identical(model, "auto")) {
      model <- best
      message("   model = 'auto': using '", model, "'. Note that with three anchors each fit uses two points, so ",
              "this test is suggestive rather than decisive - confirm with gradseq_compare_models().")
    } else if (best != model) {
      message("   You asked for model = '", model, "'. Re-run with model = '", best,
              "' to use the best-generalising one, and compare the resulting f/f0 distributions.")
    }
  } else {
    message("Only ", nrow(anchors), " anchors: every model fits them exactly and NONE can be validated. ",
            "Supply a 70S fraction (anchor_70S = <n>) to make the fit testable.")
    if (identical(model, "auto")) {
      model <- "zero_anchored"
      message("   model = 'auto' with fewer than three anchors: falling back to 'zero_anchored', the only one carrying a physical constraint.")
    }
  }

  predict_s <- switch(model,
    zero_anchored = function(f) pmax(bz * (f - load_fraction), 0),
    power = function(f) { d <- f - load_fraction; out <- rep(0, length(d))
                          ok <- is.finite(d) & d > 0; out[ok] <- exp(cp[1]) * d[ok]^cp[2]; out },
    free_linear = function(f) unname(stats::predict(fit_free, data.frame(fraction = f))))

  # ---- plausibility: where SHOULD ordinary globular proteins land under this calibration? --------
  # An external check that uses no ribosome: a compact protein of known mass has a predictable s
  # (verified against published standards by gradseq_selftest()). Invert the calibration to ask which
  # fraction it should peak in, then look at where such proteins actually are.
  gridf <- seq(min(fr), max(fr), by = 0.02); grids <- predict_s(gridf)
  chk <- data.table(mw_kDa = c(20, 50, 100, 250, 500))
  chk[, s_if_globular := .s_from_ffo(mw_kDa * 1000, 1.25)]
  if (all(diff(grids) > 0)) {                       # invertible only if s increases with fraction
    chk[, fraction_predicted := suppressWarnings(stats::approx(grids, gridf, xout = s_if_globular, rule = 2)$y)]
    message("Plausibility check - where a COMPACT (f/f0 = 1.25) protein should peak under this calibration:")
    print(chk[, .(mw_kDa, s_if_globular = round(s_if_globular, 2), fraction_predicted = round(fraction_predicted, 1))])
    message("   Compare with where proteins of that mass actually peak in gradseq_profiles.csv. If the predicted fractions are ",
            "EARLIER (nearer the top) than the observed ones, the calibration over-estimates s and f/f0 will come out too small.")
  } else {
    chk[, fraction_predicted := NA_real_]
    message("Plausibility check skipped: s is not monotonically increasing across the gradient under this model.")
  }

  if (save_plots) {
    dir.create(.gs_dir(), recursive = TRUE, showWarnings = FALSE)
    curves <- rbindlist(list(
      data.table(fraction = gridf, s = pmax(bz * (gridf - load_fraction), 0), model = "zero_anchored"),
      data.table(fraction = gridf, s = ifelse(gridf > load_fraction, exp(cp[1]) * pmax(gridf - load_fraction, 0)^cp[2], 0), model = "power"),
      data.table(fraction = gridf, s = cf[1] + cf[2] * gridf, model = "free_linear")))
    mdl <- model                                  # the column is also called `model`; keep them apart
    curves[, in_use := model == mdl]
    ylo <- max(-20, min(-5, min(curves$s, na.rm = TRUE)))
    gp <- ggplot(curves, aes(fraction, s, colour = model, linewidth = in_use)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2, ymax = 10, fill = "grey60", alpha = 0.18) +
      annotate("text", x = min(fr), y = 10, hjust = 0, vjust = -0.4, size = 2.8, colour = "grey30",
               label = "where ordinary proteins actually live (s = 2-10) - NO anchor is here") +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_line() +
      geom_point(data = data.table(fraction = load_fraction, s = 0), aes(fraction, s),
                 colour = "black", shape = 4, size = 3, inherit.aes = FALSE) +
      geom_point(data = anchors, aes(fraction, s), colour = "black", size = 3, inherit.aes = FALSE) +
      geom_text(data = anchors, aes(fraction, s, label = name), vjust = -1, colour = "black", inherit.aes = FALSE) +
      geom_point(data = chk[is.finite(fraction_predicted)], aes(fraction_predicted, s_if_globular),
                 shape = 21, fill = "white", size = 2.2, inherit.aes = FALSE) +
      geom_text(data = chk[is.finite(fraction_predicted)],
                aes(fraction_predicted, s_if_globular, label = paste0(mw_kDa, " kDa")),
                hjust = -0.15, size = 2.5, colour = "grey25", inherit.aes = FALSE) +
      scale_linewidth_manual(values = c(`FALSE` = 0.4, `TRUE` = 1.1), guide = "none") +
      labs(title = "Gradient calibration: sedimentation coefficient vs fraction",
           subtitle = paste0("Anchors (filled) are the ribosomal particles; the cross at fraction ", round(load_fraction, 2),
                             " is the physical s = 0 point - a particle that does not sediment stays in the load zone.\n",
                             "Open circles mark where a COMPACT protein of that mass should peak. The bold curve is the model in use ('", model, "').\n",
                             "EVERY protein value is an extrapolation into the shaded band, where there is no anchor at all - this remains the weakest step."),
           x = "fraction", y = "s (Svedberg)", colour = NULL) +
      coord_cartesian(ylim = c(ylo, max(anchors$s) * 1.15)) + theme_bw() + theme(legend.position = "top")
    tryCatch(ggsave(.gs_dir("gradseq_calibration.pdf"), gp, width = 8, height = 5.5), error = function(e) NULL)

    L <- rbindlist(list(cbind(small, subunit = "30S (rps*)"), cbind(large, subunit = "50S (rpl*)")))
    gr <- ggplot(L, aes(peak_fraction, fill = subunit)) +
      geom_histogram(binwidth = 1, position = "identity", alpha = 0.6) +
      geom_vline(data = anchors, aes(xintercept = fraction), linetype = 2) +
      geom_text(data = anchors, aes(x = fraction, y = Inf, label = name), vjust = 1.4, size = 3,
                colour = "grey20", inherit.aes = FALSE) +
      labs(title = "Ribosomal subunit peak fractions (sanity check)",
           subtitle = "The two subunits should peak in clearly distinct fractions, as in the published A260 profile.\nThe 70S anchor should sit below both - it is where rps* and rpl* co-migrate.",
           x = "peak fraction", y = "proteins", fill = NULL) + theme_bw()
    gc2 <- ggplot(data.table(fraction = rep(fr, 3),
                             value = c(ps, pl, co),
                             trace = rep(c("30S proteins (rps*)", "50S proteins (rpl*)", "co-migration = sqrt(rps* x rpl*)"), each = length(fr))),
                  aes(fraction, value, colour = trace)) +
      geom_line(linewidth = 0.8) +
      geom_vline(data = anchors, aes(xintercept = fraction), linetype = 2, colour = "grey40") +
      geom_text(data = anchors, aes(x = fraction, y = Inf, label = name), vjust = 1.4, size = 3, colour = "grey20", inherit.aes = FALSE) +
      labs(title = "Locating the 70S monosome",
           subtitle = "Free 30S carries no rpl*, free 50S carries no rps*; only the monosome carries both, so the co-migration\ntrace peaks at 70S. Each subunit trace is the mean normalised profile, scaled to its own maximum.",
           x = "fraction", y = "relative signal", colour = NULL) + theme_bw() + theme(legend.position = "top")
    tryCatch({ grDevices::pdf(.gs_dir("gradseq_ribosome_check.pdf"), width = 7.5, height = 5)
               print(gr); print(gc2); grDevices::dev.off() },
             error = function(e) try(grDevices::dev.off(), silent = TRUE))
  }
  invisible(list(fit = fit_free, predict_s = predict_s, model = model, anchors = anchors,
                 load_fraction = load_fraction, loo = loo, gene_map = gene_map,
                 fits = list(zero_anchored = fit_zero, power = fit_pow, free_linear = fit_free)))
}

# ---- 2b. which calibration model does the PROTEOME reject? -----------------------------------------
# Leave-one-out on three anchors is a weak test - each fit uses two points. This is the decisive one, and
# it needs no extra data: f/f0 >= 1 is a hard physical floor, so the fraction of proteins a model pushes
# below 1 is the fraction of answers it is definitely getting wrong. A model that makes a quarter of the
# proteome sub-spherical has been falsified, whatever its anchor residuals look like.
# It cuts both ways: values above ~4 are implausible for folded proteins, so a model producing many of
# those is under-estimating s just as badly.
gradseq_compare_models <- function(gs, cal, position = c("com", "peak"), mass_map = NULL, vbar = .VBAR) {
  position <- match.arg(position)
  D <- copy(gs$meta)
  if (is.null(mass_map)) {
    sf <- here("output", "uniprot_annotation_shared.RData")
    if (!file.exists(sf)) stop("No UniProt cache for monomer masses - render a comparison first, or pass mass_map.")
    e <- new.env(); load(sf, envir = e); u <- as.data.table(e$.uniprot_all)
    mass_map <- setNames(as.numeric(u$mass), as.character(u$input_id))
  }
  D[, mw_Da := unname(mass_map[protein_id])]
  D[, frac_used := if (position == "com") com_fraction else as.numeric(peak_fraction)]
  D <- D[is.finite(mw_Da) & mw_Da > 0 & is.finite(frac_used)]
  if (is.null(cal$fits)) stop("This `cal` predates model comparison - re-run gradseq_calibrate().")
  lf <- cal$load_fraction
  bz <- unname(stats::coef(cal$fits$zero_anchored)[1])
  cp <- stats::coef(cal$fits$power); cf <- stats::coef(cal$fits$free_linear)
  preds <- list(
    zero_anchored = pmax(bz * (D$frac_used - lf), 0),
    power         = ifelse(D$frac_used > lf, exp(cp[1]) * pmax(D$frac_used - lf, 0)^cp[2], 0),
    free_linear   = cf[1] + cf[2] * D$frac_used)
  R <- rbindlist(lapply(names(preds), function(nm) {
    s <- preds[[nm]]; ok <- is.finite(s) & s > 0
    f <- .ffo_from_s(s[ok], D$mw_Da[ok], vbar)
    data.table(model = nm, n_usable = sum(ok), pct_no_s = round(100 * mean(!ok), 1),
               median_ffo = round(stats::median(f, na.rm = TRUE), 2),
               pct_impossible = round(100 * mean(f < 1, na.rm = TRUE), 1),
               pct_over_4     = round(100 * mean(f > 4, na.rm = TRUE), 1))
  }))
  R[, verdict := data.table::fcase(
    pct_impossible > 25, "FALSIFIED - a quarter of the proteome below the physical floor",
    pct_impossible > 10, "poor - a large impossible population",
    pct_over_4     > 25, "poor - implausibly extended population (s under-estimated near the load zone)",
    default             = "usable")]
  message("Which calibration model survives contact with the proteome? (f/f0 < 1 is impossible; > 4 is implausible for folded proteins)")
  print(R)
  message("   In use: '", cal$model, "'. Choose the model with the smallest impossible AND implausible population,")
  message("   then re-run gradseq_calibrate(gs, model = '<name>'). If none is usable, the anchors cannot calibrate")
  message("   this mass range at all - say so, and use gradseq_vs_sec_deviation(), which needs no calibration.")
  invisible(R)
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
  # use the calibration's own predictor (zero-anchored / power / free linear); older cal objects only
  # carried an lm, so fall back to that
  D[, s_svedberg := if (is.function(cal$predict_s)) cal$predict_s(frac_used)
                    else unname(stats::predict(cal$fit, data.frame(fraction = frac_used)))]
  n_before <- nrow(D[is.finite(mw_Da) & mw_Da > 0])
  D <- D[is.finite(mw_Da) & mw_Da > 0 & is.finite(s_svedberg) & s_svedberg > 0]
  if (!nrow(D)) stop("No protein has both a monomer mass and a positive s - check the calibration.")
  if (n_before - nrow(D) > 0)
    message(sprintf("   %d of %d protein(s) (%.0f%%) were dropped for a NON-POSITIVE s (they sit at or above the load zone, or the fit runs negative there).",
                    n_before - nrow(D), n_before, 100 * (n_before - nrow(D)) / n_before))
  # Proteins peaking within a fraction of the load zone have barely sedimented: their position is set by
  # the width of the load band and by diffusion, not by s. Flag them rather than pretend to measure them.
  lf <- if (is.numeric(cal$load_fraction)) cal$load_fraction else NA_real_
  D[, near_load_zone := is.finite(lf) & (frac_used - lf) < 1]
  if (any(D$near_load_zone))
    message(sprintf("   %d protein(s) (%.0f%%) peak within one fraction of the load zone; their s is set by the load band and diffusion, not by sedimentation. Column `near_load_zone` flags them - exclude before quoting anything.",
                    sum(D$near_load_zone), 100 * mean(D$near_load_zone)))
  D[, ffo_gradseq := .ffo_from_s(s_svedberg, mw_Da, vbar)]

  # PHYSICAL SANITY CHECK. f/f0 cannot be below 1: the sphere is the minimum-friction shape. A sizeable
  # sub-1 population therefore proves the s calibration is wrong, not that the proteins are compact.
  .imposs <- mean(D$ffo_gradseq < 1, na.rm = TRUE)
  if (.imposs > 0.05)
    warning(sprintf(paste0("%.0f%% of proteins get f/f0 < 1, which is PHYSICALLY IMPOSSIBLE (a sphere is the minimum). ",
                           "The calibration over-estimates s for ordinary proteins: the ribosomal anchors sit at s = 30-70 ",
                           "while ordinary proteins are at s = 2-10, so the fit is extrapolated far below its anchors. ",
                           "Try the other models (gradseq_calibrate(gs, model = 'power') or 'zero_anchored'), compare their ",
                           "leave-one-out errors, add a low-s anchor of known s via extra_anchors, or use the calibration-free ",
                           "route gradseq_vs_sec_deviation()."), 100 * .imposs), call. = FALSE)
  # the mirror-image failure: forcing the curve through the load zone can under-estimate s near the top
  .huge <- mean(D$ffo_gradseq > 4, na.rm = TRUE)
  if (.huge > 0.10)
    message(sprintf("   %.0f%% of proteins get f/f0 > 4, which is extreme even for disordered chains. Near the load zone the fitted s tends to zero, which inflates f/f0 without limit - check `near_load_zone` before reading anything into the tail.",
                    100 * .huge))
  setorder(D, -ffo_gradseq)
  fwrite(D, .gs_dir("gradseq_ffo.csv"))
  message(sprintf("Sedimentation-derived ABSOLUTE f/f0 for %d protein(s): median %.2f (IQR %.2f-%.2f).",
                  nrow(D), stats::median(D$ffo_gradseq, na.rm = TRUE),
                  stats::quantile(D$ffo_gradseq, .25, na.rm = TRUE), stats::quantile(D$ffo_gradseq, .75, na.rm = TRUE)))
  message("   Reference: ~1.2 = compact globular (hydrated), ~1.5 = moderately elongated, >2 = extended/disordered.")
  message("   NOTE the monomer assumption DEFLATES this estimate for oligomers (by n^(2/3)), so a high value is a strong claim.")
  message(sprintf("   Physical check: %.0f%% of values are below 1 (impossible). Above ~5%% the absolute scale should not be used - prefer gradseq_vs_sec_deviation().",
                  100 * .imposs))
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
gradseq_vs_sec <- function(ff, metabolite, condition = NULL, globular_ffo = 1.2, vbar = .VBAR,
                           restrict_to_calibrated = TRUE, save_plots = TRUE) {
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
  # The Stokes radius below is derived from this project's apparent MW, so outside the calibrated
  # interval it would be built on an extrapolation of the standards curve - and every downstream
  # quantity (f/f0, native mass, n_implied) would inherit that.
  if (restrict_to_calibrated) {
    if ("in_calibrated_range" %in% names(G)) {
      n0 <- nrow(G); G <- G[in_calibrated_range %in% TRUE]
      message("Restricted to the calibrated MW interval: ", nrow(G), " of ", n0,
              " proteins (", round(100 * nrow(G) / n0), "%). Set restrict_to_calibrated = FALSE to use all.")
    } else message("No in_calibrated_range column - re-run globularity_check() to get it; using all proteins.")
  }
  if (!nrow(G)) stop("No proteins left after restricting to the calibrated range.")

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
    # f/f0 < 1 is physically impossible (the sphere is the minimum-friction shape), so a curve sitting
    # below 1 is proof that its calibration is wrong - stamp that on the figure rather than let it be read
    # as "very compact".
    .bad <- c(SEC = stats::median(J$ffo_sec_absolute, na.rm = TRUE),
              `Grad-seq` = mgrd, combined = mtru)
    .bad <- names(.bad)[is.finite(.bad) & .bad < 1]
    refs <- data.table(x = c(1.0, 1.2, 1.5, 2.0), lab = c("SPHERE (hard floor)", "globular", "elongated", "extended"))
    g1 <- ggplot(L, aes(ffo, fill = src)) +
      annotate("rect", xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf, fill = "grey55", alpha = 0.22) +
      geom_density(alpha = 0.4, colour = NA) +
      geom_vline(data = refs, aes(xintercept = x), linetype = 2, colour = "grey35", inherit.aes = FALSE) +
      geom_text(data = refs, aes(x = x, y = Inf, label = lab), vjust = 1.4, hjust = -0.05, angle = 90,
                size = 2.7, colour = "grey25", inherit.aes = FALSE) +
      scale_x_log10() +
      labs(title = paste0("How compact is the proteome? Two independent techniques (", metabolite, " control)"),
           subtitle = paste0("All on the ABSOLUTE scale. The SEC curve is the relative f/f0 multiplied by ", globular_ffo,
                             " (the absolute f/f0 of the globular calibrants).\n",
                             "Oligomers inflate the SEC estimate (n^1/3) and deflate the sedimentation one (n^2/3); the combined curve removes both.\n",
                             "SHADED REGION IS PHYSICALLY IMPOSSIBLE: f/f0 >= 1 always, since a sphere has the least friction for a given mass.",
                             if (length(.bad)) paste0("\n*** ", paste(.bad, collapse = " and "),
                                                      " sits below 1 -> that calibration is WRONG, not the proteins. Do not interpret it; use gradseq_vs_sec_deviation(). ***") else ""),
           x = "absolute f/f0", y = "density", fill = NULL) +
      theme_bw() + theme(legend.position = "top")
    if (length(.bad))
      message("!! ", paste(.bad, collapse = " and "), " median f/f0 is below 1, which is physically impossible - ",
              "the sedimentation calibration is unusable. Use gradseq_vs_sec_deviation() instead.")
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

# ---- 5. CALIBRATION-FREE alternative -----------------------------------------------------------------
# The absolute route above needs a trustworthy fraction -> s calibration, and ribosomal anchors alone do
# not provide one for ordinary proteins (they sit at s = 30-50, the proteome at s = 2-10, so the fit is
# extrapolated far below its anchors and returns impossible f/f0 < 1). This route avoids s entirely.
#
# In BOTH datasets, reduce each protein to how far it migrates relative to what its monomer mass predicts:
#     deviation = log10(mass expected at this position) - log10(monomer mass)
# with the expectation taken from a robust fit WITHIN that dataset. Positive = migrates as though heavier.
# Sedimentation and SEC weight mass and shape differently (s ~ M^(2/3)/(f/f0) versus R_s ~ (f/f0)M^(1/3)),
# so the two deviations are NOT the same quantity and their magnitudes should not be equated - but a
# protein that is anomalous for its mass should be anomalous in both, and that is what is tested here.
gradseq_vs_sec_deviation <- function(gs, metabolite, condition = NULL, position = c("com", "peak"),
                                     restrict_to_calibrated = TRUE, dev_cut = log10(2),
                                     mass_map = NULL, save_plots = TRUE) {
  position <- match.arg(position)
  D <- copy(gs$meta)
  if (is.null(mass_map)) {
    sf <- here("output", "uniprot_annotation_shared.RData")
    if (!file.exists(sf)) stop("No UniProt cache for monomer masses.")
    e <- new.env(); load(sf, envir = e); u <- as.data.table(e$.uniprot_all)
    mass_map <- setNames(as.numeric(u$mass), as.character(u$input_id))
  }
  D[, mw_kDa := unname(mass_map[protein_id]) / 1000]
  D[, pos := if (position == "com") com_fraction else as.numeric(peak_fraction)]
  D <- D[is.finite(mw_kDa) & mw_kDa > 0 & is.finite(pos)]
  fit <- if (requireNamespace("MASS", quietly = TRUE)) MASS::rlm(log10(mw_kDa) ~ pos, data = D) else stats::lm(log10(mw_kDa) ~ pos, data = D)
  D[, deviation_log10_gradseq := stats::predict(fit, D) - log10(mw_kDa)]

  gf <- here("output", paste0("PCM_ctrl_vs_", metabolite), "tables", "globularity_check.txt")
  if (!file.exists(gf)) stop("No globularity_check.txt for ", metabolite, " - run globularity_check() first.")
  G <- fread(gf)
  if ("condition" %in% names(G)) {
    cn <- unique(as.character(G$condition))
    cc <- if (!is.null(condition)) condition else { x <- cn[grepl("ctrl|control|ref", cn, ignore.case = TRUE)][1]; if (is.na(x)) cn[1] else x }
    G <- G[condition == cc]
  }
  if (restrict_to_calibrated && "in_calibrated_range" %in% names(G)) {
    n0 <- nrow(G); G <- G[in_calibrated_range %in% TRUE]
    message("Restricted to the calibrated MW interval: ", nrow(G), " of ", n0, " proteins.")
  }
  G[, deviation_log10_ours := log10(ratio)]
  J <- merge(G[, .(protein_id, deviation_log10_ours, class)],
             D[, .(protein_id, mw_kDa, pos, deviation_log10_gradseq)], by = "protein_id")
  J <- J[is.finite(deviation_log10_ours) & is.finite(deviation_log10_gradseq)]
  if (!nrow(J)) stop("No shared proteins.")
  ct <- suppressWarnings(stats::cor.test(J$deviation_log10_ours, J$deviation_log10_gradseq, method = "spearman"))
  message("Proteins in both datasets: ", nrow(J))
  message(sprintf("Deviation agreement (SEC vs sedimentation), Spearman rho = %+.3f (p = %.3g).",
                  unname(ct$estimate), ct$p.value))
  message("   The two techniques weight mass and shape differently, so only the AGREEMENT is meaningful, not the magnitudes.")
  J[, reproducible := abs(deviation_log10_ours) > dev_cut & abs(deviation_log10_gradseq) > dev_cut &
                      sign(deviation_log10_ours) == sign(deviation_log10_gradseq)]
  message(sprintf("Anomalous (>%.1f-fold) and in the SAME direction in both: %d protein(s) (%.1f%%).",
                  10^dev_cut, sum(J$reproducible), 100 * mean(J$reproducible)))
  dir.create(.gs_dir(), recursive = TRUE, showWarnings = FALSE)
  fwrite(J, .gs_dir("sec_vs_gradseq_deviation.csv"))
  if (save_plots) {
    g <- ggplot(J, aes(deviation_log10_gradseq, deviation_log10_ours)) +
      geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
      geom_point(aes(colour = reproducible), alpha = 0.45, size = 0.9) +
      scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#E15759"), name = "anomalous in both") +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "black", linewidth = 0.6) +
      labs(title = paste0("SEC vs sedimentation: do the same proteins migrate anomalously?  (", metabolite, ")"),
           subtitle = sprintf("Calibration-free: each dataset is referenced to its own bulk trend of mass against migration position.\nSpearman rho = %+.3f (p = %.3g, n = %d). Magnitudes are NOT comparable between techniques - only the agreement is.",
                              unname(ct$estimate), ct$p.value, nrow(J)),
           x = "deviation, Grad-seq (sedimentation)", y = "deviation, this study (SEC)") + theme_bw()
    tryCatch(ggsave(.gs_dir("sec_vs_gradseq_deviation.pdf"), g, width = 7.5, height = 5.5), error = function(e) NULL)
  }
  invisible(J)
}

gradseq_all <- function(file, metabolite, id_col = NULL, fraction_cols = NULL, globular_ffo = 1.2,
                        anchor_70S = "auto", model = "auto", load_fraction = NULL) {
  gradseq_selftest()
  gs  <- gradseq_load(file, id_col = id_col, fraction_cols = fraction_cols)
  cal <- gradseq_calibrate(gs, anchor_70S = anchor_70S, model = model, load_fraction = load_fraction)
  try(gradseq_compare_models(gs, cal), silent = TRUE)
  ff  <- gradseq_ffo(gs, cal)
  gradseq_vs_sec(ff, metabolite = metabolite, globular_ffo = globular_ffo)
}
