# scripts/reconcile_calibration_claims.R
# =============================================================================
# Checks every number in the "the calibration only covers half the proteins" argument, against the data,
# before any of them goes into a draft. Run it, read the console, quote the numbers it prints.
#
# It answers four questions:
#   1. Do 53% / 80.9% / 44% actually agree with each other? They are over-determined: if out-of-range
#      proteins are forced non-globular, then
#           pooled %globular = (share in range) x (%globular in range)
#      and the three reported figures do not satisfy that identity. One of them is being computed
#      differently (pooled vs median-across-metabolites, or void counted separately). This finds out which.
#   2. What is the HONEST globularity rate outside the calibrated range? globularity_check.R forces
#      `globular_as_expected := FALSE` for every out-of-range protein AFTER computing dev_fractions for it.
#      So the reported 0% is by fiat. But dev_fractions is measured at n x the protein's OWN monomer mass,
#      which for an E. coli proteome is an INTERPOLATION even for a protein eluting in fraction 1 - so the
#      underlying measurement is sound and can simply be read back out.
#   3. Is the doubling of the surface-hydrophobicity rho explained by range restriction alone?
#      Restricting to a window of the dependent variable attenuates a correlation by
#           rho_restricted = rho*u / sqrt(1 - rho^2*(1 - u^2)),   u = sd(dev|in-range) / sd(dev|all)
#      An exact doubling on releasing the restriction needs u ~= 0.49, and that threshold is almost flat
#      in rho (0.498 at rho=0.10, 0.468 at rho=0.40). So u is the whole test.
#   4. Does f/f0 fall below the physical floor of 1 more often outside the calibrated range than inside?
#      This is the ONE test that is not circular: f/f0 >= 1 is imposed by physics, not by the calibration,
#      so a floor violation outside and not inside is the calibration failing where it is extrapolated.
#
# USAGE:
#   source(here::here("scripts", "reconcile_calibration_claims.R"))
#   reconcile_calibration_claims()                 # all metabolites found under output/
#   reconcile_calibration_claims("ATP")            # one
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table) })

reconcile_calibration_claims <- function(metabolites = NULL, condition = NULL,
                                         surface_csv = NULL, globular_ffo = 1.25,
                                         tolerance_fractions = 1, max_oligomer = 4L,
                                         verbose = TRUE) {
  if (is.null(metabolites)) {
    d <- basename(list.dirs(here("output"), recursive = FALSE))
    metabolites <- sub("^PCM_ctrl_vs_", "", d[grepl("^PCM_ctrl_vs_", d)])
  }
  if (!length(metabolites)) stop("No PCM_ctrl_vs_* directories under output/.")

  per <- rbindlist(lapply(metabolites, function(m) {
    f <- here("output", paste0("PCM_ctrl_vs_", m), "tables", "globularity_check.txt")
    if (!file.exists(f)) { message("[", m, "] no globularity_check.txt - skipped"); return(NULL) }
    G <- fread(f)
    if ("condition" %in% names(G)) {
      cn <- unique(as.character(G$condition))
      cc <- if (!is.null(condition)) condition else { x <- cn[grepl("ctrl|control|ref", cn, ignore.case = TRUE)][1]; if (is.na(x)) cn[1] else x }
      G <- G[condition == cc]
    }
    if (!nrow(G)) return(NULL)
    need <- c("in_calibrated_range", "globular_as_expected")
    if (!all(need %in% names(G))) { message("[", m, "] missing ", paste(setdiff(need, names(G)), collapse = "/"), " - re-run globularity_check()"); return(NULL) }

    inr <- G$in_calibrated_range %in% TRUE
    # the honest out-of-range rate: dev_fractions is measured at n x the protein's own mass, so it is
    # meaningful even out of range - the override is what zeroes it, not the measurement
    honest_out <- if ("dev_fractions" %in% names(G)) mean(G$dev_fractions[!inr] <= 1, na.rm = TRUE) else NA_real_
    honest_in  <- if ("dev_fractions" %in% names(G)) mean(G$dev_fractions[inr]  <= 1, na.rm = TRUE) else NA_real_
    # THE NON-CIRCULAR TEST - with the RIGHT threshold. `ffo_vs_monomer` is (apparent/expected)^(1/3),
    # which is RELATIVE to the globular calibrants: a protein behaving exactly like them scores 1.0 BY
    # CONSTRUCTION. It is NOT an absolute frictional ratio, so `ffo_vs_monomer < 1` is not a physics
    # violation at all - it just means "more compact than a typical globular standard", which is possible.
    # The absolute ratio is ffo_vs_monomer * globular_ffo (~1.25 measured for the classical standards by
    # gradseq_selftest), so the PHYSICALLY IMPOSSIBLE region is
    #        ffo_vs_monomer < 1 / globular_ffo  =  0.80
    # Testing against 1.0 instead of 0.80 counts a large, ordinary population as "impossible" and destroys
    # the argument. That was a bug in the first version of this script.
    ffo <- if ("ffo_vs_monomer" %in% names(G)) G$ffo_vs_monomer else rep(NA_real_, nrow(G))
    ffo_floor <- 1 / globular_ffo
    data.table(
      metabolite        = m,
      n                 = nrow(G),
      pct_in_range      = 100 * mean(inr),
      pct_glob_in       = 100 * mean(G$globular_as_expected[inr] %in% TRUE),
      pct_glob_out      = 100 * mean(G$globular_as_expected[!inr] %in% TRUE),
      pct_glob_pooled   = 100 * mean(G$globular_as_expected %in% TRUE),
      honest_glob_in    = 100 * honest_in,
      honest_glob_out   = 100 * honest_out,
      pct_impossible_in  = 100 * mean(ffo[inr]  < ffo_floor, na.rm = TRUE),
      pct_impossible_out = 100 * mean(ffo[!inr] < ffo_floor, na.rm = TRUE),
      # kept for comparison: the WRONG threshold the first version used
      pct_below1_in      = 100 * mean(ffo[inr]  < 1, na.rm = TRUE),
      pct_below1_out     = 100 * mean(ffo[!inr] < 1, na.rm = TRUE))
  }))
  if (!nrow(per)) stop("Nothing could be read.")

  cat("\n=========== PER METABOLITE ===========\n")
  print(per[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 1) else z)])

  # ---- 1. the identity check -----------------------------------------------------------------------
  per[, implied_pooled := pct_in_range / 100 * pct_glob_in + (1 - pct_in_range / 100) * pct_glob_out]
  per[, identity_gap := pct_glob_pooled - implied_pooled]
  cat("\n=========== 1. DO THE THREE NUMBERS AGREE? ===========\n")
  cat("Identity: pooled %globular = share_in_range x %glob_in + share_out x %glob_out\n")
  print(per[, .(metabolite, pct_in_range = round(pct_in_range, 1), pct_glob_in = round(pct_glob_in, 1),
                pct_glob_out = round(pct_glob_out, 1), reported_pooled = round(pct_glob_pooled, 1),
                implied = round(implied_pooled, 1), gap = round(identity_gap, 2))])
  if (max(abs(per$identity_gap), na.rm = TRUE) < 0.5)
    cat("=> Consistent. The per-metabolite figures satisfy the identity, so any mismatch with a quoted\n",
        "   headline number is a POOLED-vs-MEDIAN issue, not an error. Below:\n", sep = "")
  else
    cat("=> INCONSISTENT by up to", round(max(abs(per$identity_gap)), 2),
        "points. Something is counted differently (void? below vs beyond?). Resolve before quoting.\n")
  cat(sprintf("\n   pooled over ALL metabolites : in-range %.1f%% | glob-in %.1f%% | glob-pooled %.1f%%\n",
              100 * sum(per$n * per$pct_in_range / 100) / sum(per$n),
              100 * sum(per$n * per$pct_in_range / 100 * per$pct_glob_in / 100) / sum(per$n * per$pct_in_range / 100),
              100 * sum(per$n * per$pct_glob_pooled / 100) / sum(per$n)))
  cat(sprintf("   MEDIAN across metabolites   : in-range %.1f%% | glob-in %.1f%% | glob-pooled %.1f%%\n",
              median(per$pct_in_range), median(per$pct_glob_in), median(per$pct_glob_pooled)))
  cat("   -> quote ONE of these two rows and say which. They are different statistics.\n")

  # ---- 2. the honest out-of-range rate --------------------------------------------------------------
  cat("\n=========== 2. HONEST GLOBULARITY OUTSIDE THE CALIBRATED RANGE ===========\n")
  if (all(is.na(per$honest_glob_out))) {
    cat("dev_fractions not in the table - re-run globularity_check() to get it.\n")
  } else {
    cat("`globular_as_expected` is FORCED FALSE out of range by globularity_check.R, so the reported\n")
    cat("out-of-range rate is 0% by construction. dev_fractions is measured at n x the protein's OWN\n")
    cat("monomer mass - an interpolation - so it survives the override and can be read back:\n\n")
    print(per[, .(metabolite,
                  reported_glob_out = round(pct_glob_out, 1),
                  honest_glob_out   = round(honest_glob_out, 1),
                  honest_glob_in    = round(honest_glob_in, 1))])
    cat(sprintf("\n=> Out of range, %.1f%% of proteins still elute within one fraction of a clean oligomer state\n",
                median(per$honest_glob_out, na.rm = TRUE)))
    cat(sprintf("   (in range: %.1f%%). Quote THIS contrast, not 80.9%% vs 0%%.\n",
                median(per$honest_glob_in, na.rm = TRUE)))
  }

  # ---- 2b. THE GEOMETRIC NULL -----------------------------------------------------------------------
  # The objection that matters most: is "80.9% globular inside the calibrated range" impressive, or is it
  # what you would get from ANY set of elution positions falling in that window?
  # A protein counts as globular if its apex is within `tolerance_fractions` of ANY of its 1x..max_oligomer
  # expected positions. Those positions are only log10(max_oligomer)/|slope| fractions apart, and each
  # carries +-tolerance. So the accepting band has a FIXED width, and the calibrated window has a fixed
  # width, and their ratio is the share that would pass by geometry alone with no biology involved.
  cat("\n=========== 2b. HOW MUCH OF 80.9%% IS JUST GEOMETRY? ===========\n")
  cal <- tryCatch({
    f <- here("output", paste0("PCM_ctrl_vs_", per$metabolite[1]), "tables", "globularity_standards_check.txt")
    if (file.exists(f)) fread(f) else NULL }, error = function(e) NULL)
  if (is.null(cal) || !all(c("elution_fraction", "expected_kDa") %in% names(cal))) {
    cat("globularity_standards_check.txt not found - cannot compute the geometric null.\n")
  } else {
    ok <- is.finite(cal$elution_fraction) & is.finite(cal$expected_kDa) & cal$expected_kDa > 0
    sl <- stats::coef(stats::lm(log10(cal$expected_kDa[ok]) ~ cal$elution_fraction[ok]))[2]
    win <- diff(range(cal$elution_fraction[ok]))
    band <- log10(max_oligomer) / abs(sl) + 2 * tolerance_fractions
    null_pct <- 100 * min(band / win, 1)
    obs <- 100 * sum(per$n * per$pct_in_range / 100 * per$pct_glob_in / 100) / sum(per$n * per$pct_in_range / 100)
    cat(sprintf("calibration slope        = %.3f decades/fraction (%.2fx MW per fraction)\n", sl, 10^abs(sl)))
    cat(sprintf("calibrated window        = %.2f fractions wide\n", win))
    cat(sprintf("1x..%dx positions span    = %.2f fractions; +-%g tolerance each side\n",
                max_oligomer, log10(max_oligomer)/abs(sl), tolerance_fractions))
    cat(sprintf("=> accepting band        = %.2f fractions = %.1f%% of the window\n", band, null_pct))
    cat(sprintf("\n   GEOMETRIC NULL  %.1f%%   vs   OBSERVED  %.1f%%   ->  excess %+.1f points\n",
                null_pct, obs, obs - null_pct))
    if (obs - null_pct < 15)
      cat("   Most of the in-range globularity rate is the accepting band being nearly as wide as the\n",
          "   calibrated window. Do NOT quote 80.9%% as though it were a measurement of the proteome.\n",
          "   Quote the EXCESS over this null, or narrow tolerance_fractions / max_oligomer and re-run.\n", sep="")
    else
      cat("   A substantial excess over geometry - the in-range rate carries real information.\n")
  }

  # ---- 3. range restriction ------------------------------------------------------------------------
  cat("\n=========== 3. IS THE rho DOUBLING JUST RANGE RESTRICTION? ===========\n")
  scsv <- if (!is.null(surface_csv)) surface_csv else here("output", "surface", "surface_vs_elution_allproteins.csv")
  if (!file.exists(scsv)) {
    cat("Need", scsv, "- run surface_vs_elution(restrict_to_calibrated = FALSE) first.\n")
    u <- NA_real_
  } else {
    S <- fread(scsv)
    gm <- metabolites[1]
    gf <- here("output", paste0("PCM_ctrl_vs_", gm), "tables", "globularity_check.txt")
    G <- fread(gf)
    if ("condition" %in% names(G)) {
      cn <- unique(as.character(G$condition))
      cc <- if (!is.null(condition)) condition else { x <- cn[grepl("ctrl|control|ref", cn, ignore.case = TRUE)][1]; if (is.na(x)) cn[1] else x }
      G <- G[condition == cc]
    }
    S <- merge(S, G[, .(protein_id, in_calibrated_range)], by = "protein_id", all.x = TRUE)
    sd_in  <- stats::sd(S[in_calibrated_range %in% TRUE]$dev_log2, na.rm = TRUE)
    sd_all <- stats::sd(S$dev_log2, na.rm = TRUE)
    u <- sd_in / sd_all
    cat(sprintf("sd(dev_log2 | in range) = %.3f\nsd(dev_log2 | all)      = %.3f\nu = %.3f\n", sd_in, sd_all, u))
    cat("\nThreshold for range restriction to explain an EXACT doubling: u ~= 0.49\n")
    cat(if (abs(u - 0.49) < 0.05)
      "=> u is at the threshold. Range restriction explains the doubling COMPLETELY. Do NOT present the\n   doubling as evidence about the calibration - there is nothing left to attribute to biology.\n"
      else if (u > 0.54)
      "=> u is ABOVE the threshold, so restriction UNDER-explains the doubling. A residual effect remains\n   and is worth reporting - but verify it with deviation = 'residual' before claiming it.\n"
      else
      "=> u is BELOW the threshold, so restriction OVER-explains it: the association is genuinely WEAKER\n   outside the calibrated range. That is the opposite of the reading you had.\n")
  }

  # ---- 4. the non-circular test --------------------------------------------------------------------
  cat("\n=========== 4. THE ONE NON-CIRCULAR TEST: f/f0 BELOW THE PHYSICAL FLOOR ===========\n")
  if (all(is.na(per$pct_impossible_out))) {
    cat("ffo_vs_monomer not in the table - re-run globularity_check().\n")
  } else {
    cat("f/f0 >= 1 is imposed by PHYSICS, not by the calibration: a sphere has the least friction for a\n")
    cat("given mass. So a floor violation OUTSIDE and not INSIDE is the calibration failing exactly where\n")
    cat("it is extrapolated - and unlike the deviation argument, nothing about it is circular.\n\n")
    cat(sprintf("Threshold: ffo_vs_monomer < 1/globular_ffo = %.2f. (Testing against 1.0 is WRONG -\n", 1/globular_ffo))
    cat("ffo_vs_monomer is RELATIVE to the calibrants, which score 1.0 by construction.)\n\n")
    print(per[, .(metabolite,
                  impossible_in  = round(pct_impossible_in, 1),
                  impossible_out = round(pct_impossible_out, 1),
                  `below1_in(wrong)`  = round(pct_below1_in, 1),
                  `below1_out(wrong)` = round(pct_below1_out, 1))])
    cat(sprintf("\n=> median: %.1f%% impossible in range vs %.1f%% out of range.\n",
                median(per$pct_impossible_in, na.rm = TRUE), median(per$pct_impossible_out, na.rm = TRUE)))
    cat(if (median(per$pct_impossible_in, na.rm = TRUE) > 10)
      "   WARNING: a large impossible population INSIDE the calibrated range too. That is not the\n   calibration failing where extrapolated - it is failing everywhere, and this test cannot then be\n   used to argue that the extrapolated region is special.\n"
      else if (median(per$pct_impossible_out, na.rm = TRUE) > median(per$pct_impossible_in, na.rm = TRUE) + 5)
      "   THIS IS YOUR STRONGEST NUMBER. Physics falsifying the calibration precisely where it is\n   extrapolated, with no definitional circularity anywhere in it.\n"
      else "   No clear contrast, so this line of argument is not available. Rely on the coverage figure instead.\n")
  }

  dir.create(here("output", "surface"), recursive = TRUE, showWarnings = FALSE)
  fwrite(per, here("output", "surface", "calibration_claims_reconciled.csv"))
  cat("\nTable -> ", here("output", "surface", "calibration_claims_reconciled.csv"), "\n", sep = "")
  invisible(per)
}
