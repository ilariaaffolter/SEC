# scripts/secseq_compare.R
# =============================================================================
# CROSS-LAB VALIDATION of the elution behaviour, against a published E. coli SEC-seq complexome
# (Chihara et al., same separation principle as this project: size-exclusion chromatography, 20 fractions).
#
# WHY THIS AND NOT THE GRADIENT DATA: SEC-seq measures the SAME physical quantity as this project - the
# Stokes radius - so it CANNOT give an independent frictional ratio (both share the identical mass/shape
# degeneracy). What it can do, and what a sedimentation dataset cannot, is test REPRODUCIBILITY: if a
# protein elutes away from the position expected for its monomer mass here AND does so in an independent
# lab, with a different column, buffer and growth condition, then the anomaly is a property of the
# protein, not of this chromatography. That argument needs no calibration transfer and no assumption
# about the absolute f/f0 of the standards, which makes it considerably harder to attack than any
# cross-technique comparison. Use scripts/gradseq_ffo.R with a glycerol-gradient dataset for the
# orthogonal, absolute f/f0.
#
# THE MEASUREMENT USED HERE IS CALIBRATION-FREE. Rather than converting the published fractions into
# molecular weights (their standards are not in the table, and anchoring only on the ribosome would mean
# extrapolating from ~1 MDa down to ~10 kDa), each dataset is reduced to a DEVIATION:
#      deviation = log10( mass expected at this elution position ) - log10( monomer mass )
# where "mass expected at this elution position" comes from a robust fit of monomer mass against elution
# position WITHIN that dataset. Positive = the protein elutes as though it were heavier than its monomer
# (assembly and/or an extended shape); negative = it elutes late for its mass. Because each dataset is
# referenced to its own bulk behaviour, the two are directly comparable without sharing a calibration.
#
# ORIENTATION IS DETERMINED FROM THE DATA, NOT ASSUMED. Fraction numbering runs in opposite directions in
# the two datasets. The script settles it two independent ways - the ribosomal proteins (rps*/rpl*, the
# largest species present, must sit at the high-mass end) and the sign of the global mass-vs-position
# correlation - and refuses to continue if they disagree.
#
# INPUT: the supplementary protein table, e.g. Chihara_2023_RNA_SEC_Seq_Supplemental_Table_S2.xlsx, with
#        columns "Protein IDs", "Gene names", "Molecular weight [kDa]" and SEC-Fraction-1..20.
#        Values are normalised to 1 at each protein's own maximum, which is all this analysis needs.
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "secseq_compare.R"))
#   sq <- secseq_load("data/raw/Chihara_2023_RNA_SEC_Seq_Supplemental_Table_S2.xlsx")
#   secseq_orientation(sq)                       # check this before anything else
#   secseq_selfcheck(sq)                         # how anomalous is elution in THEIR data alone?
#   secseq_vs_sec(sq, metabolite = "ATP")        # the cross-lab comparison
#   secseq_all("data/raw/Chihara_2023_...xlsx", metabolite = "ATP")
#
# OUTPUT (output/secseq/):
#   secseq_profiles.csv             parsed profiles, peak and centre-of-mass fraction, monomer mass
#   secseq_orientation.pdf          the two orientation tests
#   secseq_selfcheck.pdf            their monomer mass vs elution position, with the fitted bulk trend
#   secseq_deviation.csv            per protein: elution deviation in their data
#   secseq_vs_sec.csv               per-protein join with this project's SEC
#   secseq_vs_sec.pdf               position agreement, deviation agreement, and the reproducible set
#   secseq_reproducible_anomalies.csv   proteins that elute anomalously in BOTH datasets
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })

.sq_dir <- function(...) here("output", "secseq", ...)

# Supplementary workbooks often carry one or more banner rows above the real header (e.g. merged group
# labels such as "Features" / "Normalized by maximum value"). Scan the first rows for the one that holds
# the actual column names and return how many rows to skip.
.find_header_skip <- function(file, sheet = 1, hints = c("protein.?id", "gene name", "fraction", "accession"),
                              max_scan = 12) {
  ext <- tolower(tools::file_ext(file))
  probe <- if (ext %in% c("xlsx", "xls", "xlsm")) {
    if (!requireNamespace("readxl", quietly = TRUE)) return(0L)
    suppressMessages(as.data.frame(readxl::read_excel(file, sheet = sheet, col_names = FALSE,
                                                      n_max = max_scan, .name_repair = "minimal")))
  } else {
    ln <- readLines(file, n = max_scan, warn = FALSE)
    as.data.frame(do.call(rbind, lapply(strsplit(ln, "\t|,"), function(x) { length(x) <- max(lengths(strsplit(ln, "\t|,"))); x })),
                  stringsAsFactors = FALSE)
  }
  if (!nrow(probe)) return(0L)
  for (i in seq_len(nrow(probe))) {
    txt <- paste(as.character(unlist(probe[i, ])), collapse = " | ")
    if (sum(vapply(hints, function(h) grepl(h, txt, ignore.case = TRUE), logical(1))) >= 2) return(i - 1L)
  }
  0L
}

# ---- 1. load ---------------------------------------------------------------------------------------
# the most recent parse, so the later steps can be called without passing it back in every time
.sq_cache <- new.env(parent = emptyenv())

secseq_load <- function(file, sheet = 1, id_col = NULL, mw_col = NULL, gene_col = NULL,
                        fraction_prefix = "SEC-Fraction", skip = NULL) {
  if (!file.exists(file)) { f2 <- here(file); if (file.exists(f2)) file <- f2 else stop("File not found: ", file) }
  if (is.null(skip)) {
    skip <- .find_header_skip(file, sheet)
    if (skip > 0) message("Detected ", skip, " banner row(s) above the header - skipping them. Override with skip = <n>.")
  }
  X <- if (tolower(tools::file_ext(file)) %in% c("xlsx", "xls", "xlsm")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("install.packages('readxl') to read this file")
    suppressMessages(as.data.table(readxl::read_excel(file, sheet = sheet, skip = skip)))
  } else as.data.table(data.table::fread(file, skip = skip))
  nm <- names(X)
  pick <- function(given, patterns, what) {
    if (!is.null(given)) return(given)
    for (p in patterns) { hit <- grep(p, nm, ignore.case = TRUE, value = TRUE); if (length(hit)) return(hit[1]) }
    stop("Could not find the ", what, " column. Pass it explicitly. Columns: ", paste(nm, collapse = ", "))
  }
  id_col   <- pick(id_col,   c("^Protein IDs$", "protein.?id", "accession", "uniprot"), "protein id")
  mw_col   <- pick(mw_col,   c("^Molecular weight", "molecular.?weight", "\\bmw\\b", "mass"), "molecular weight")
  gene_col <- tryCatch(pick(gene_col, c("^Gene names$", "gene"), "gene name"), error = function(e) NA_character_)

  fr <- grep(paste0("^", fraction_prefix), nm, value = TRUE)
  if (length(fr) < 5) {
    fr <- grep("fraction", nm, ignore.case = TRUE, value = TRUE)
    if (length(fr) < 5) stop("Fewer than 5 fraction columns found. Columns: ", paste(nm, collapse = ", "))
  }
  fno <- suppressWarnings(as.numeric(gsub("[^0-9]", "", fr)))
  ok  <- is.finite(fno); fr <- fr[ok]; fno <- fno[ok]
  o   <- order(fno); fr <- fr[o]; fno <- fno[o]
  message("Using ", length(fr), " fraction column(s): ", fr[1], " ... ", fr[length(fr)])

  num <- function(v) suppressWarnings(as.numeric(as.character(v)))
  M <- as.matrix(as.data.frame(lapply(X[, ..fr], num))); M[!is.finite(M)] <- 0
  ids <- sub("[;,].*$", "", trimws(as.character(X[[id_col]])))          # first accession of a group
  mw  <- num(X[[mw_col]])
  gene <- if (!is.na(gene_col)) tolower(sub("[;,].*$", "", trimws(as.character(X[[gene_col]])))) else NA_character_

  keep <- nzchar(ids) & is.finite(mw) & mw > 0 & rowSums(M) > 0
  M <- M[keep, , drop = FALSE]; ids <- ids[keep]; mw <- mw[keep]
  gene <- if (length(gene) > 1) gene[keep] else rep(NA_character_, length(ids))

  Mn <- M / rowSums(M)
  D <- data.table(protein_id = ids, gene = gene, mw_kDa = mw,
                  peak_fraction = fno[max.col(Mn, ties.method = "first")],
                  com_fraction  = as.vector(Mn %*% fno))
  dir.create(.sq_dir(), recursive = TRUE, showWarnings = FALSE)
  fwrite(cbind(D, as.data.table(Mn)), .sq_dir("secseq_profiles.csv"))
  message("Parsed ", nrow(D), " protein profile(s) over fractions ", min(fno), "-", max(fno), ".")
  out <- list(meta = D, profiles = Mn, fractions = fno, file = file)
  .sq_cache$sq <- out          # so secseq_vs_sec() can be called without re-passing it
  invisible(out)
}

# ---- 2. which way round is the fraction axis? ------------------------------------------------------
secseq_orientation <- function(sq, position = c("peak", "com"), save_plots = TRUE) {
  position <- match.arg(position)
  D <- copy(sq$meta)
  # ONE position measure for both tests. Using the peak for one and the centre of mass for the other made
  # them answer subtly different questions, and a spurious disagreement here aborts the whole comparison.
  D[, pos := if (position == "com") com_fraction else as.numeric(peak_fraction)]
  rib <- D[grepl("^rp[sl][a-z]$", gene)]
  # test 1: the ribosome (~1-2.5 MDa) must lie at the high-mass end of the axis
  t1 <- if (nrow(rib) >= 5) {
    med_rib <- stats::median(rib$pos); med_all <- stats::median(D$pos)
    list(ok = TRUE, ribosome_fraction = med_rib, overall_fraction = med_all,
         high_mass_at = if (med_rib < med_all) "low fractions" else "high fractions", n = nrow(rib))
  } else list(ok = FALSE, n = nrow(rib))
  # test 2: sign of the global monomer-mass vs position correlation
  rho <- suppressWarnings(stats::cor(log10(D$mw_kDa), D$pos, method = "spearman", use = "complete.obs"))
  t2_high_mass_at <- if (rho < 0) "low fractions" else "high fractions"

  if (isTRUE(t1$ok)) {
    message(sprintf("Orientation test 1 (ribosome, n = %d): rps*/rpl* peak at fraction %.1f vs %.1f overall -> high mass at %s.",
                    t1$n, t1$ribosome_fraction, t1$overall_fraction, t1$high_mass_at))
  } else message("Orientation test 1 skipped: only ", t1$n, " ribosomal protein(s) matched (need >= 5). Check the gene column.")
  message(sprintf("Orientation test 2 (mass vs position): Spearman rho = %+.3f -> high mass at %s.", rho, t2_high_mass_at))
  if (isTRUE(t1$ok) && t1$high_mass_at != t2_high_mass_at)
    stop("The two orientation tests DISAGREE. Inspect secseq_orientation.pdf before going further - ",
         "something is wrong with the fraction parsing or the gene mapping.")
  high_mass_at <- if (isTRUE(t1$ok)) t1$high_mass_at else t2_high_mass_at
  message("=> Adopted orientation: HIGH MASS at ", high_mass_at,
          if (high_mass_at == "high fractions") "  (opposite to this project's SEC, where fraction 1 is the void)" else "  (same as this project's SEC)")

  if (save_plots) {
    dir.create(.sq_dir(), recursive = TRUE, showWarnings = FALSE)
    g1 <- ggplot(D, aes(pos, mw_kDa)) +
      geom_point(alpha = 0.25, size = 0.7, colour = "grey45") +
      { if (nrow(rib)) geom_point(data = rib, aes(pos, mw_kDa), colour = "firebrick", size = 1.4) else NULL } +
      scale_y_log10() +
      labs(title = "Orientation of the published SEC fraction axis",
           subtitle = sprintf("Red = ribosomal proteins (the ~1-2.5 MDa particle, so they mark the high-mass end).\nSpearman rho(log mass, position) = %+.3f  =>  high mass at %s.", rho, high_mass_at),
           x = paste0(position, " fraction (published data)"), y = "monomer mass (kDa)") + theme_bw()
    tryCatch(ggsave(.sq_dir("secseq_orientation.pdf"), g1, width = 7, height = 5), error = function(e) NULL)
  }
  invisible(list(high_mass_at = high_mass_at, rho = rho, ribosome = rib))
}

# ---- 3. deviation from the bulk trend, within the published data alone -----------------------------
# A robust fit of log10(monomer mass) against elution position describes how a TYPICAL protein of this
# proteome elutes on that column. The residual is then "how much heavier the protein would have to be to
# elute here", i.e. exactly the quantity globularity_check reports as log(apparent/expected).
secseq_selfcheck <- function(sq, orientation = NULL, position = c("peak", "com"), save_plots = TRUE) {
  # DEFAULT "peak": this project's side of the comparison uses apex_fraction, which is a peak. Comparing a
  # peak against a centre of mass compares two different quantities - a centre of mass is pulled along the
  # axis by the tail of the profile - so the two datasets must be summarised the same way.
  position <- match.arg(position)
  if (is.null(orientation)) orientation <- secseq_orientation(sq, position = position, save_plots = FALSE)
  D <- copy(sq$meta)
  D[, pos := if (position == "com") com_fraction else as.numeric(peak_fraction)]
  D <- D[is.finite(pos) & is.finite(mw_kDa) & mw_kDa > 0]
  fit <- if (requireNamespace("MASS", quietly = TRUE)) MASS::rlm(log10(mw_kDa) ~ pos, data = D)
         else stats::lm(log10(mw_kDa) ~ pos, data = D)
  D[, expected_log10_mw := stats::predict(fit, D)]
  # positive = elutes as though HEAVIER than its monomer (assembled and/or extended)
  D[, deviation_log10 := expected_log10_mw - log10(mw_kDa)]
  setorder(D, -deviation_log10)
  fwrite(D[, .(protein_id, gene, mw_kDa, pos, expected_log10_mw, deviation_log10)], .sq_dir("secseq_deviation.csv"))
  message(sprintf("Published data: median |deviation| = %.2f log10 units (%.1f-fold); %.1f%% of proteins deviate by more than 2-fold.",
                  stats::median(abs(D$deviation_log10), na.rm = TRUE),
                  10^stats::median(abs(D$deviation_log10), na.rm = TRUE),
                  100 * mean(abs(D$deviation_log10) > log10(2), na.rm = TRUE)))
  if (save_plots) {
    g <- ggplot(D, aes(pos, mw_kDa)) +
      geom_point(alpha = 0.25, size = 0.7, colour = "grey45") +
      geom_line(aes(y = 10^expected_log10_mw), colour = "steelblue", linewidth = 0.9) +
      scale_y_log10() +
      labs(title = "Published SEC-seq: monomer mass vs elution position",
           subtitle = "Blue = robust fit describing how a typical protein of this proteome elutes.\nPoints far ABOVE the line elute late for their mass; far BELOW, they elute as though much heavier.",
           x = paste0(position, " fraction"), y = "monomer mass (kDa)") + theme_bw()
    tryCatch(ggsave(.sq_dir("secseq_selfcheck.pdf"), g, width = 7, height = 5), error = function(e) NULL)
  }
  invisible(D)
}

# ---- 4. the cross-lab comparison -------------------------------------------------------------------
secseq_vs_sec <- function(sq = NULL, metabolite, condition = NULL, dev_cut = log10(2),
                          position = c("peak", "com"),
                          restrict_to_calibrated = FALSE, save_plots = TRUE) {
  position <- match.arg(position)
  # restrict_to_calibrated defaults to FALSE because the deviation below is computed from a fit of mass
  # against elution POSITION and never touches the standards curve - so proteins outside the calibrated
  # MW interval are perfectly usable here, and excluding them would only discard data.
  if (is.null(sq)) {
    if (is.null(.sq_cache$sq))
      stop("No SEC-seq data loaded. Run  sq <- secseq_load(\"data/raw/<the SEC-seq table>.xlsx\")  first, ",
           "then either secseq_vs_sec(sq, metabolite = \"...\") or just secseq_vs_sec(metabolite = \"...\").")
    sq <- .sq_cache$sq
    message("Using the SEC-seq table loaded earlier: ",
            if (is.null(sq$file)) "<cached>" else basename(as.character(sq$file)), ".")
  }
  S <- secseq_selfcheck(sq, position = position, save_plots = FALSE)
  gf <- here("output", paste0("PCM_ctrl_vs_", metabolite), "tables", "globularity_check.txt")
  if (!file.exists(gf)) stop("No globularity_check.txt for ", metabolite, " - run globularity_check() first.")
  G <- fread(gf)
  if ("condition" %in% names(G)) {
    cn <- unique(as.character(G$condition))
    cc <- if (!is.null(condition)) condition else { x <- cn[grepl("ctrl|control|ref", cn, ignore.case = TRUE)][1]; if (is.na(x)) cn[1] else x }
    G <- G[condition == cc]; message("Using the '", cc, "' rows of this project's SEC table.")
  }
  if (!all(c("protein_id", "expected_mw_kDa", "apex_fraction") %in% names(G)))
    stop("globularity_check.txt lacks protein_id/expected_mw_kDa/apex_fraction.")
  # Restrict to proteins whose apex falls INSIDE the calibrated MW interval. Outside it this project's
  # apparent MW - and therefore its deviation - is an extrapolation of the standards curve, so including
  # those proteins would compare a measurement against an extrapolation.
  if (restrict_to_calibrated) {
    if ("in_calibrated_range" %in% names(G)) {
      n0 <- nrow(G); G <- G[in_calibrated_range %in% TRUE]
      message("Restricted to the calibrated MW interval: ", nrow(G), " of ", n0,
              " proteins (", round(100 * nrow(G) / n0), "%). Set restrict_to_calibrated = FALSE to use all.")
    } else message("No in_calibrated_range column - re-run globularity_check() to get it; using all proteins.")
  }
  if (!nrow(G)) stop("No proteins left after restricting to the calibrated range.")

  # THE DEVIATION MUST BE COMPUTED THE SAME WAY ON BOTH SIDES. It is the residual from a robust fit of
  # monomer mass against elution position WITHIN each dataset - so each is referenced to its own bulk
  # behaviour and NO calibration is used on either side. (Using log10(apparent/expected) here instead
  # would import this project's standards calibration onto one axis only, making the two axes different
  # quantities: the calibration-derived ratio runs to 10^6 on extrapolated proteins while a residual is
  # bounded, and correlating the two is meaningless.)
  # This project's fraction axis runs the other way (fraction 1 = void), which the fitted slope absorbs;
  # the sign convention - positive = migrates as though HEAVIER than its monomer - is preserved.
  # BOTH FITS ARE DONE ON THE SHARED SET, AGAINST ONE SINGLE MASS COLUMN. Previously each deviation was
  # fitted on its own population against its own mass annotation (expected_mw_kDa here, mw_kDa there), so
  # the "shared" -log10(mass) term was only APPROXIMATELY shared - and partialling on one of the two could
  # not fully clean both axes. Fitting both here, on the same proteins and the same masses, makes the term
  # removed literally identical on both sides, which is what makes the partial correlation exact.
  G <- G[is.finite(expected_mw_kDa) & expected_mw_kDa > 0 & is.finite(apex_fraction)]
  if ("ratio" %in% names(G)) G[, calibration_ratio_log10 := log10(ratio)]   # kept for reference only
  keepcols <- intersect(c("protein_id", "apex_fraction", "calibration_ratio_log10",
                          "expected_mw_kDa", "apparent_mw_kDa", "class"), names(G))
  J <- merge(G[, ..keepcols],
             S[, .(protein_id, gene, mw_kDa, pos_published = pos)],
             by = "protein_id")
  if (!nrow(J)) stop("No shared proteins - check the accession formats in both tables.")
  message("Proteins measured in BOTH SEC datasets: ", nrow(J),
          " (this study ", nrow(G), ", published ", nrow(S), ").")
  # sanity: do the two labs agree on the monomer mass they used?
  mwdiff <- abs(J$expected_mw_kDa - J$mw_kDa) / pmax(J$mw_kDa, 1)
  message(sprintf("Monomer mass agreement between the two annotation sources: %.1f%% within 5%% - using expected_mw_kDa for BOTH axes.",
                  100 * mean(mwdiff < 0.05, na.rm = TRUE)))
  J[, lgm := log10(expected_mw_kDa)]
  .rob <- function(f, dat) if (requireNamespace("MASS", quietly = TRUE)) MASS::rlm(f, data = dat) else stats::lm(f, data = dat)
  fitG <- .rob(lgm ~ apex_fraction, J)
  fitS <- .rob(lgm ~ pos_published, J)
  J[, deviation_log10_ours      := stats::predict(fitG, J) - lgm]
  J[, deviation_log10_published := stats::predict(fitS, J) - lgm]

  P <- J[is.finite(deviation_log10_ours) & is.finite(deviation_log10_published)]
  if (!nrow(P)) stop("No shared proteins with finite deviations.")
  rho_pos <- suppressWarnings(stats::cor(P$apex_fraction, P$pos_published, method = "spearman"))
  ct  <- suppressWarnings(stats::cor.test(P$deviation_log10_ours, P$deviation_log10_published, method = "spearman"))
  # How much positional signal does each dataset actually carry? This is what sets the scale of each
  # deviation axis - a deviation is (that dataset's mass-position slope) x (residual position) - so a
  # dataset whose position barely tracks mass yields a compressed axis and a mass-adjusted deviation that
  # is close to noise. Measure it and say so, rather than plotting it as though it were a measurement.
  .r2f <- function(f) suppressWarnings(stats::cor(stats::predict(f, P), P$lgm, use = "complete.obs")^2)
  r2_ours <- .r2f(fitG); r2_pub <- .r2f(fitS)
  message(sprintf("Positional signal - share of the variance in monomer mass explained by elution position: this study %.1f%%, published %.1f%% (slopes %.4f and %.4f log10 units per fraction).",
                  100 * r2_ours, 100 * r2_pub, stats::coef(fitG)[2], stats::coef(fitS)[2]))
  .weak <- c("this study", "the published data")[c(r2_ours, r2_pub) < 0.05]
  if (length(.weak))
    warning("In ", paste(.weak, collapse = " and "), " elution position explains under 5% of the variance in monomer ",
            "mass, so that axis carries almost no positional information and its mass-adjusted deviation is close to ",
            "noise. Treat the partial rho as an upper bound near zero, not as a measurement of agreement.",
            call. = FALSE, immediate. = TRUE)

  # PARTIAL correlation controlling for monomer mass - THE statistic that matters here.
  # Both deviations are of the form fitted(position) - log10(monomer mass), so they SHARE the
  # -log10(mass) term. Whenever position predicts mass poorly - which is the very thing under
  # investigation - each deviation collapses towards -(log10 mass - mean) and the raw correlation
  # approaches 1 for that reason alone, telling us nothing about whether the two EXPERIMENTS agree.
  # Removing the shared mass leaves the question actually being asked: does a protein sit in the same
  # relative position in both runs, beyond what its mass already dictates?
  .pcor <- function(x, y, z) {
    ok <- is.finite(x) & is.finite(y) & is.finite(z); n <- sum(ok)
    if (n < 10) return(list(rho = NA_real_, p = NA_real_, n = n))
    rx <- rank(x[ok]); ry <- rank(y[ok]); rz <- rank(z[ok])
    rxy <- stats::cor(rx, ry); rxz <- stats::cor(rx, rz); ryz <- stats::cor(ry, rz)
    den <- sqrt((1 - rxz^2) * (1 - ryz^2))
    if (!is.finite(den) || den <= 0) return(list(rho = NA_real_, p = NA_real_, n = n))
    r <- (rxy - rxz * ryz) / den
    tt <- r * sqrt((n - 3) / max(1 - r^2, .Machine$double.eps))
    list(rho = r, p = 2 * stats::pt(-abs(tt), df = n - 3), n = n)
  }
  pc <- .pcor(P$deviation_log10_ours, P$deviation_log10_published, log10(P$expected_mw_kDa))
  # the same adjustment, in linear form, so it can be plotted: what is left of each deviation
  # once the shared monomer-mass term has been regressed out of it.
  .resid_on <- function(y, z) {           # full-length residuals, NA-safe
    out <- rep(NA_real_, length(y)); ok <- is.finite(y) & is.finite(z)
    if (sum(ok) >= 3L) out[ok] <- stats::residuals(stats::lm(y[ok] ~ z[ok]))
    out
  }
  # lgm was set on J above and inherited here; the control is that SAME column, so the term removed from
  # the two axes is identical rather than merely similar.
  P[, dev_ours_adj := .resid_on(deviation_log10_ours,      lgm)]
  P[, dev_pub_adj  := .resid_on(deviation_log10_published, lgm)]
  message(sprintf("Elution POSITION agreement (Spearman, note the axes run opposite ways): rho = %+.3f", rho_pos))
  message(sprintf("DEVIATION agreement - do both labs flag the SAME proteins as eluting off their monomer mass?\n   raw Spearman rho = %+.3f (p = %.3g, n = %d)",
                  unname(ct$estimate), ct$p.value, nrow(P)))
  message(sprintf("   PARTIAL Spearman rho, monomer mass held constant = %+.3f (p = %.3g)  <- this is the honest number",
                  pc$rho, pc$p))
  # Both deviations contain the term -log10(monomer mass) by construction. If elution position is a
  # poor predictor of mass, each deviation is dominated by that shared term and the RAW correlation
  # is high whatever the experiments did. Only the partial correlation speaks to the experiments.
  if (is.finite(pc$rho) && is.finite(unname(ct$estimate))) {
    if (abs(unname(ct$estimate)) > 0.6 && abs(pc$rho) < 0.3) {
      warning("The raw deviation correlation (", sprintf("%+.3f", unname(ct$estimate)),
              ") is largely an artefact: both deviations share the term -log10(monomer mass), and with mass ",
              "held constant the agreement drops to ", sprintf("%+.3f", pc$rho),
              ". Report the PARTIAL value, not the raw one.", call. = FALSE, immediate. = TRUE)
      message("   => Do NOT quote the raw rho. The two datasets agree mainly because they are plotted against the same monomer masses.")
    } else if (abs(pc$rho) >= 0.3) {
      message("   => The agreement survives removing the shared mass term, so it reflects reproducible chromatographic behaviour.")
    }
  }
  # how much of each deviation is just the mass term? (R^2 of deviation on log10 mass)
  r2o <- suppressWarnings(stats::cor(P$deviation_log10_ours,      P$lgm, use = "complete.obs")^2)
  r2p <- suppressWarnings(stats::cor(P$deviation_log10_published, P$lgm, use = "complete.obs")^2)
  message(sprintf("   (monomer mass alone explains %.1f%% of the deviation here and %.1f%% in the published data)",
                  100 * r2o, 100 * r2p))
  if (min(r2_ours, r2_pub) < 0.05)
    message("   => NOT INTERPRETABLE: one dataset carries almost no positional signal, so its mass-adjusted deviation ",
            "is close to noise. Read the partial rho as an upper bound near zero.")
  # Read that the other way round and it is the interpretable number: because the deviation IS the
  # residual of log10(mass) ~ position, the share of the deviation NOT explained by mass is exactly the
  # share of log10(mass) that elution position does explain. Report it directly, and say plainly that
  # restricting the mass range (restrict_to_calibrated) deflates it - a classic range-restriction effect,
  # not a measure of how badly the column performs.
  message(sprintf("   Equivalently: elution position explains %.0f%% of the variance in monomer mass here and %.0f%% in the published data.",
                  100 * (1 - r2o), 100 * (1 - r2p)))
  if (isTRUE(restrict_to_calibrated))
    message("   NOTE you restricted the mass range, which mechanically lowers both figures. Re-run with restrict_to_calibrated = FALSE to see them on the full range.")

  # the reproducible set: anomalous, in the same direction, in both datasets
  P[, anom_ours := abs(deviation_log10_ours) > dev_cut]
  P[, anom_pub  := abs(deviation_log10_published) > dev_cut]
  P[, same_direction := sign(deviation_log10_ours) == sign(deviation_log10_published)]
  P[, reproducible := anom_ours & anom_pub & same_direction]
  # ...and the same set after removing the shared mass term, so the shortlist cannot be an artefact
  # of a protein simply being small or large in both datasets.
  # THE CUTOFF CANNOT BE dev_cut HERE. The mass-adjusted deviation is no longer a mass ratio: it is
  # (slope of the mass-position fit) x (residual elution position), so its scale is compressed by
  # however weak that fit is, and a fold-change threshold would simply never be met - it silently
  # returned an empty set. Each axis therefore gets its own robust spread-based cutoff, and the
  # fold-equivalent is reported so the number stays interpretable.
  so <- stats::mad(P$dev_ours_adj, na.rm = TRUE); sp <- stats::mad(P$dev_pub_adj, na.rm = TRUE)
  if (!is.finite(so) || so <= 0) so <- stats::sd(P$dev_ours_adj, na.rm = TRUE)
  if (!is.finite(sp) || sp <= 0) sp <- stats::sd(P$dev_pub_adj, na.rm = TRUE)
  cut_o <- 2 * so; cut_p <- 2 * sp
  P[, reproducible_massadj := is.finite(dev_ours_adj) & is.finite(dev_pub_adj) &
        abs(dev_ours_adj) > cut_o & abs(dev_pub_adj) > cut_p &
        sign(dev_ours_adj) == sign(dev_pub_adj)]
  message(sprintf("Anomalous (>%.1f-fold) here: %.1f%% | in the published data: %.1f%% | REPRODUCIBLE in both, same direction: %d protein(s) (%.1f%% of the shared set).",
                  10^dev_cut, 100 * mean(P$anom_ours), 100 * mean(P$anom_pub),
                  sum(P$reproducible), 100 * mean(P$reproducible)))
  message(sprintf("   after removing the shared monomer-mass term: cutoff = 2 x robust SD, i.e. %.3f here and %.3f in the published data (log10 units, NOT fold-change);\n   %d protein(s) (%.1f%%) are beyond it in both datasets and in the same direction - use THIS set.",
                  cut_o, cut_p, sum(P$reproducible_massadj), 100 * mean(P$reproducible_massadj)))
  message(sprintf("   For scale: the mass-adjusted deviations span %.3f to %.3f (this study), so the fold-change cutoff used above (%.2f-fold) does not apply on this axis.",
                  min(P$dev_ours_adj, na.rm = TRUE), max(P$dev_ours_adj, na.rm = TRUE), 10^dev_cut))
  fwrite(P, .sq_dir("secseq_vs_sec.csv"))
  fwrite(P[reproducible == TRUE][order(-abs(deviation_log10_ours))], .sq_dir("secseq_reproducible_anomalies.csv"))
  fwrite(P[reproducible_massadj == TRUE][order(-abs(dev_ours_adj))], .sq_dir("secseq_reproducible_anomalies_massadj.csv"))

  if (save_plots) {
    g1 <- ggplot(P, aes(deviation_log10_published, deviation_log10_ours)) +
      geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45") +
      geom_point(aes(colour = reproducible), alpha = 0.45, size = 0.9) +
      scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#E15759"), name = "anomalous in both") +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "black", linewidth = 0.6) +
      labs(title = paste0("Do two independent SEC experiments agree on which proteins elute anomalously?  (", metabolite, " control)"),
           subtitle = sprintf("Deviation = log10(mass expected at the elution position) - log10(monomer mass); positive = elutes as though heavier.\nBoth axes are residuals from a within-dataset fit - neither uses a standards calibration.\nRaw Spearman rho = %+.3f (p = %.3g, n = %d)  BUT both axes contain the same -log10(monomer mass) term,\nso read the PARTIAL rho instead: %+.3f (p = %.3g). See the next panel.",
                              unname(ct$estimate), ct$p.value, nrow(P), pc$rho, pc$p),
           x = "deviation, published SEC-seq", y = "deviation, this study") + theme_bw()
    # the same comparison with the shared monomer-mass term regressed out of both axes
    g1b <- ggplot(P[is.finite(dev_ours_adj) & is.finite(dev_pub_adj)], aes(dev_pub_adj, dev_ours_adj)) +
      geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45") +
      geom_point(aes(colour = reproducible_massadj), alpha = 0.45, size = 0.9) +
      scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#E15759"), name = "anomalous in both\n(mass-adjusted)") +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "black", linewidth = 0.6) +
      labs(title = "The same comparison with monomer mass removed from both axes",
           subtitle = sprintf("Each axis is the deviation after regressing out log10(monomer mass): what the chromatography adds beyond what\nthe protein's mass already dictates. THE UNITS ARE NO LONGER FOLD-CHANGES - each axis is (slope of that\ndataset's mass-position fit) x (residual elution position), so the two scales are not comparable and only the\nagreement is. Partial Spearman rho = %+.3f (p = %.3g, n = %d); red = beyond 2 robust SD in both, same direction.\nThe parallel diagonal streaks are the fraction grid: proteins sharing an elution fraction in both datasets\nfall on a line as their mass varies.",
                              pc$rho, pc$p, pc$n),
           x = "mass-adjusted deviation, published SEC-seq", y = "mass-adjusted deviation, this study") + theme_bw()
    g2 <- ggplot(melt(P[, .(protein_id, `this study` = deviation_log10_ours, `published` = deviation_log10_published)],
                      id.vars = "protein_id", variable.name = "dataset", value.name = "deviation"),
                 aes(deviation, fill = dataset)) +
      geom_density(alpha = 0.45, colour = NA) +
      geom_vline(xintercept = c(-dev_cut, 0, dev_cut), linetype = c(3, 1, 3), colour = "grey40") +
      labs(title = "Distribution of elution deviation in both datasets",
           subtitle = paste0("Dotted lines mark a ", round(10^dev_cut, 1), "-fold deviation. If the two distributions have a similar spread,\nthe anomalous population is a property of the proteome rather than of one column."),
           x = "log10(expected mass at position / monomer mass)", y = "density", fill = NULL) +
      theme_bw() + theme(legend.position = "top")
    tryCatch({ grDevices::pdf(.sq_dir("secseq_vs_sec.pdf"), width = 8, height = 6)
               print(g1); print(g1b); print(g2); grDevices::dev.off() },
             error = function(e) try(grDevices::dev.off(), silent = TRUE))
  }
  attr(P, "rho_raw")     <- unname(ct$estimate)
  attr(P, "rho_partial") <- pc$rho
  invisible(P)
}

secseq_all <- function(file, metabolite, ...) {
  sq <- secseq_load(file, ...)
  secseq_orientation(sq)
  secseq_selfcheck(sq)
  secseq_vs_sec(sq, metabolite = metabolite)
}
