# scripts/shift_vs_charge.R
# =============================================================================
# ARTEFACT CONTROL for the DIFFERENTIAL result: is the ctrl -> metabolite elution shift correlated with
# a protein's charge?
#
# WHY: adding a charged metabolite changes the ionic strength of the sample. If that, rather than
# binding, drove the elution changes, the shift would track a protein's CHARGE across the proteome -
# basic proteins moving one way, acidic proteins the other. A correlation between signed shift and
# pI / net charge is the signature of that artefact; its ABSENCE defends the hits as specific.
# (The companion test in globularity_category_annotation.R asks the same of the BASELINE elution
# position in the control; this one asks it of the ctrl->treatment CHANGE, which is the sharper test.)
#
# SHIFT STATISTICS (per protein, replicates averaged per condition):
#   delta_com    centre-of-mass difference, treatment - control, IN FRACTIONS. The primary readout:
#                continuous, signed and directly interpretable (positive = elutes LATER = smaller
#                apparent MW; negative = elutes EARLIER = larger apparent MW).
#   signed_emd   Earth Mover's Distance with a sign (negative when the treatment mass sits at lower
#                fractions). Sensitive to partial / bimodal redistribution that a single shift misses.
#   best_lag     the cross-correlation peak shift, for continuity with ccf_differential_test.R.
#
# TESTS:
#   * PROTEOME-WIDE Spearman correlation of each shift statistic against pI and net charge. This is the
#     UNBIASED test - it uses every tested protein, so it cannot be inflated by selection.
#   * The same correlation restricted to the TOP-N candidates (smallest q from ccf_fdr_results.txt /
#     emd_fdr_results.txt). Because the permutation tests return no FDR-significant hits, this is an
#     EXPLORATORY ranked screen, not a hit list - and a correlation computed inside a selected set is
#     itself selection-biased, so it is reported for description only. Judge the artefact question on
#     the proteome-wide number.
#   Interpretation: |rho| is what matters, not p. With a few thousand proteins even rho ~ 0.04 reaches
#   p < 0.05 while explaining <1% of the variance.
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "shift_vs_charge.R"))
#   shift_vs_charge()                                  # all metabolites
#   shift_vs_charge("ATP", top_n = 100)                # one, top-100 candidates highlighted
#   shift_vs_charge("ATP", min_intensity = 300, buffer_pH = 7.4)
#
# OUTPUT (per metabolite):
#   tables/shift_vs_charge.txt          per protein: delta_com, signed_emd, best_lag, pI, net charge, rank
#   tables/shift_vs_charge_stats.txt    Spearman rho + p for every shift x charge pair, proteome-wide
#                                       and within the top-N
#   figures/shift_vs_charge.pdf         shift vs pI and vs net charge, top-N highlighted
# OUTPUT (pooled): output/shift_vs_charge/shift_vs_charge_summary.csv + rho_overview.pdf
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })

# pI / net-charge helpers (.charge_props, .PKA) live in the category-annotation script; that file only
# defines functions, so sourcing it here is side-effect free and avoids duplicating the pKa model.
source(here::here("scripts", "globularity_category_annotation.R"))

.svc_mat <- function(traces_obj) {
  dt <- as.data.table(traces_obj$traces)
  fc <- grep("^[0-9]+$", colnames(dt), value = TRUE); fc <- fc[order(as.numeric(fc))]
  m  <- as.matrix(dt[, ..fc]); rownames(m) <- as.character(dt$id); m
}
.svc_best_lag <- function(x, y, lag_max = 5L) {
  x[is.na(x)] <- 0; y[is.na(y)] <- 0
  if (stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  x <- x - mean(x); y <- y - mean(y); n <- length(x); lags <- (-lag_max):lag_max
  cc <- vapply(lags, function(L) if (L >= 0) sum(x[(1L + L):n] * y[1:(n - L)])
                                 else        sum(x[1:(n + L)] * y[(1L - L):n]), numeric(1))
  lags[which.max(cc)]
}
.svc_emd_signed <- function(x, y) {
  x[is.na(x)] <- 0; y[is.na(y)] <- 0; x[x < 0] <- 0; y[y < 0] <- 0
  sx <- sum(x); sy <- sum(y); if (sx <= 0 || sy <= 0) return(NA_real_)
  px <- x / sx; py <- y / sy
  d <- sum(abs(cumsum(px) - cumsum(py)))
  mux <- sum(seq_along(px) * px); muy <- sum(seq_along(py) * py)
  if (muy < mux) -d else d
}
# Spearman with the effect size kept in front: returns rho, p and n
.rho <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) return(list(rho = NA_real_, p = NA_real_, n = sum(ok)))
  ct <- tryCatch(suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman")), error = function(e) NULL)
  if (is.null(ct)) return(list(rho = NA_real_, p = NA_real_, n = sum(ok)))
  list(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

shift_vs_charge <- function(metabolites   = NULL,
                            top_n          = 100,
                            rank_by        = c("emd", "ccf"),  # which ranked screen supplies the top-N
                            min_intensity  = 0,
                            lag_max        = 5L,
                            buffer_pH      = 7.4,
                            out_subdir     = "shift_vs_charge") {
  rank_by <- match.arg(rank_by)
  if (is.null(metabolites)) {
    dirs        <- basename(list.dirs(here("output"), recursive = FALSE))
    metabolites <- sub("^PCM_ctrl_vs_", "", dirs[grepl("^PCM_ctrl_vs_", dirs)])
  }
  if (!length(metabolites)) stop("No PCM_ctrl_vs_* output folders found.")

  # pI / net charge for every cached protein (computed once)
  U <- .load_uniprot()
  if (is.null(U)) stop("No output/uniprot_annotation_shared.RData - render at least one comparison first.")
  if (!"sequence" %in% names(U)) stop("The UniProt cache has no `sequence` column - re-render so it is fetched.")
  message("Computing pI / net charge (pH ", buffer_pH, ") for ", nrow(U), " cached protein(s) ...")
  ch <- lapply(U$sequence, function(s) if (is.na(s) || !nzchar(s)) list(pI = NA_real_, net_charge = NA_real_, charge_per_res = NA_real_) else .charge_props(s, buffer_pH))
  CH <- data.table(protein_id = as.character(U$accession),
                   pI             = vapply(ch, function(x) x$pI, numeric(1)),
                   net_charge     = vapply(ch, function(x) x$net_charge, numeric(1)),
                   charge_per_res = vapply(ch, function(x) x$charge_per_res, numeric(1)))

  stat_rows <- list()
  for (m in metabolites) {
    fdir <- here("output", paste0("PCM_ctrl_vs_", m), "RData_for_further_plotting_and_analysis")
    f    <- list.files(fdir, pattern = "_for_plotting\\.RData$", full.names = TRUE)
    if (!length(f)) { message("[", m, "] no *_for_plotting.RData - render this metabolite first; skipping."); next }
    e <- new.env(); load(f[1], envir = e)
    if (!all(c("protein_traces_list", "design_matrix") %in% ls(e))) {
      message("[", m, "] file lacks protein_traces_list/design_matrix; skipping."); next }
    tl <- e$protein_traces_list; samples <- names(tl)
    mats   <- lapply(samples, function(s) .svc_mat(tl[[s]]))
    common <- Reduce(intersect, lapply(mats, rownames))
    if (length(common) < 10) { message("[", m, "] <10 shared proteins; skipping."); next }
    mats <- lapply(mats, function(M) { M <- M[common, , drop = FALSE]; M[is.na(M)] <- 0; M })

    dm    <- as.data.table(e$design_matrix)
    cond  <- as.character(dm$Condition[match(samples, as.character(dm$Sample_name))])
    conds <- unique(stats::na.omit(cond))
    ctrl  <- conds[grepl("ctrl|control|ref", conds, ignore.case = TRUE)][1]
    if (is.na(ctrl)) ctrl <- if (is.factor(dm$Condition)) as.character(levels(dm$Condition))[1] else conds[1]
    treat <- setdiff(conds, ctrl)[1]
    if (is.na(treat)) { message("[", m, "] only one condition present; skipping."); next }
    A <- Reduce(`+`, mats[which(cond == ctrl)])  / sum(cond == ctrl)     # mean ctrl profile
    B <- Reduce(`+`, mats[which(cond == treat)]) / sum(cond == treat)    # mean treatment profile
    fr <- as.numeric(colnames(A))

    com <- function(M) { s <- rowSums(M); ifelse(s > 0, as.vector(M %*% fr) / s, NA_real_) }
    D <- data.table(protein_id = common,
                    int_ctrl   = rowSums(A), int_treat = rowSums(B),
                    com_ctrl   = com(A),     com_treat = com(B))
    D[, delta_com := com_treat - com_ctrl]        # + = elutes LATER (smaller apparent MW)
    D[, signed_emd := vapply(seq_len(nrow(D)), function(i) .svc_emd_signed(A[i, ], B[i, ]), numeric(1))]
    D[, best_lag   := vapply(seq_len(nrow(D)), function(i) .svc_best_lag(A[i, ], B[i, ], lag_max), numeric(1))]
    D <- D[int_ctrl > min_intensity & int_treat > min_intensity]
    if (!nrow(D)) { message("[", m, "] no protein passes the intensity filter; skipping."); next }
    D <- merge(D, CH, by = "protein_id", all.x = TRUE)

    # ranked screen: smallest q from the permutation tests (exploratory - these have no FDR-significant hits)
    tab_dir <- here("output", paste0("PCM_ctrl_vs_", m), "tables")
    fig_dir <- here("output", paste0("PCM_ctrl_vs_", m), "figures")
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    rf <- file.path(tab_dir, if (rank_by == "emd") "emd_fdr_results.txt" else "ccf_fdr_results.txt")
    D[, rank_q := NA_integer_]
    if (file.exists(rf)) {
      R <- fread(rf)
      scol <- if ("emd" %in% names(R)) "emd" else if ("abs_lag" %in% names(R)) "abs_lag" else NA_character_
      if (all(c("protein_id", "qval") %in% names(R))) {
        R <- R[is.finite(qval)]
        ord <- if (!is.na(scol)) order(R$qval, -R[[scol]]) else order(R$qval)
        top_ids <- unique(as.character(R[ord]$protein_id))[seq_len(min(top_n, nrow(R)))]
        D[protein_id %in% top_ids, rank_q := match(protein_id, top_ids)]
      }
    } else message("[", m, "] no ", basename(rf), " - run the ", toupper(rank_by),
                   " test for the top-N overlay (the proteome-wide test still runs).")
    D[, is_top := !is.na(rank_q)]
    fwrite(D[order(rank_q, na.last = TRUE)], file.path(tab_dir, "shift_vs_charge.txt"), sep = "\t")

    # ---- correlations: proteome-wide (unbiased) and within the top-N (descriptive only) ----
    shifts  <- c("delta_com", "signed_emd", "best_lag")
    charges <- c("pI", "net_charge", "charge_per_res")
    rows <- list()
    for (s in shifts) for (cc in charges) {
      a <- .rho(D[[s]], D[[cc]])
      b <- if (any(D$is_top)) .rho(D[is_top == TRUE][[s]], D[is_top == TRUE][[cc]]) else list(rho = NA_real_, p = NA_real_, n = 0L)
      rows[[paste(s, cc)]] <- data.table(
        metabolite = m, shift = s, charge = cc,
        rho_all = round(a$rho, 4), p_all = signif(a$p, 3), n_all = a$n, r2_all_pct = round(100 * a$rho^2, 2),
        rho_top = round(b$rho, 4), p_top = signif(b$p, 3), n_top = b$n)
    }
    ST <- rbindlist(rows, use.names = TRUE)
    fwrite(ST, file.path(tab_dir, "shift_vs_charge_stats.txt"), sep = "\t")
    stat_rows[[m]] <- ST
    message("\n[", m, "] shift vs charge - proteome-wide Spearman (|rho| is what matters, not p):")
    print(ST[, .(shift, charge, rho_all, p_all, r2_all_pct, n_all, rho_top, n_top)])
    .key <- ST[shift == "delta_com" & charge == "pI"]
    if (nrow(.key) && is.finite(.key$rho_all))
      message(sprintf("[%s] headline: delta_com vs pI rho = %.3f (%.2f%% of variance). %s",
                      m, .key$rho_all, .key$r2_all_pct,
                      if (abs(.key$rho_all) < 0.15) "No meaningful charge dependence -> the shifts are not an ionic-strength artefact."
                      else "NOTE: a non-trivial charge dependence - inspect before interpreting the hits as binding."))

    # ---- plots: shift vs pI and vs net charge, top-N highlighted ----
    mk <- function(cc, lab) {
      dd <- D[is.finite(delta_com) & is.finite(get(cc))]
      if (!nrow(dd)) return(NULL)
      r <- .rho(dd$delta_com, dd[[cc]])
      ggplot(dd, aes(.data[[cc]], delta_com)) +
        geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
        geom_point(data = dd[is_top == FALSE], colour = "grey75", alpha = 0.45, size = 0.8) +
        geom_point(data = dd[is_top == TRUE], colour = "#E15759", alpha = 0.9, size = 1.6) +
        geom_smooth(method = "loess", se = TRUE, colour = "black", linewidth = 0.6, formula = y ~ x) +
        labs(title = paste0("Elution shift vs ", lab, " - PCM_ctrl_vs_", m),
             subtitle = sprintf("y = centre-of-mass shift (treatment - control), + = elutes later. Spearman rho = %.3f (p = %.3g, n = %d; %.2f%% of variance).\nRed = top %d candidates of the %s ranked screen (exploratory: no FDR-significant hits).\nA charge-dependent trend would indicate an ionic-strength artefact rather than binding.",
                                r$rho, r$p, r$n, 100 * r$rho^2, top_n, toupper(rank_by)),
             x = lab, y = "delta centre of mass (fractions)") +
        theme_bw()
    }
    p1 <- mk("pI", "predicted pI"); p2 <- mk("net_charge", paste0("net charge at pH ", buffer_pH))
    .fp <- file.path(fig_dir, "shift_vs_charge.pdf")
    tryCatch({ grDevices::pdf(.fp, width = 7.5, height = 5.5)
               if (!is.null(p1)) print(p1); if (!is.null(p2)) print(p2); grDevices::dev.off() },
             error = function(e) { message("   !! could not write ", basename(.fp), ": ", conditionMessage(e))
                                   try(grDevices::dev.off(), silent = TRUE) })
  }

  if (length(stat_rows)) {
    S <- rbindlist(stat_rows, use.names = TRUE)
    out <- here("output", out_subdir); dir.create(out, recursive = TRUE, showWarnings = FALSE)
    fwrite(S, file.path(out, "shift_vs_charge_summary.csv"))
    K <- S[shift == "delta_com" & charge == "pI"]
    cat("\n==== shift vs pI across metabolites (centre-of-mass shift) ====\n")
    print(K[, .(metabolite, rho_all, p_all, r2_all_pct, n_all)])
    cat(sprintf("\nLargest |rho| = %.3f (%s). %s\n",
                max(abs(K$rho_all), na.rm = TRUE),
                K$metabolite[which.max(abs(K$rho_all))],
                if (max(abs(K$rho_all), na.rm = TRUE) < 0.15)
                  "No metabolite shows a meaningful charge dependence: the differential shifts are not explained by an ionic-strength / charge artefact."
                else "At least one metabolite shows a non-trivial charge dependence - inspect that one before interpreting its hits."))
    g <- ggplot(S[shift == "delta_com"], aes(metabolite, rho_all, fill = charge)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.72) +
      geom_hline(yintercept = c(-0.15, 0.15), linetype = 3, colour = "grey45") +
      geom_hline(yintercept = 0, colour = "grey30") +
      coord_cartesian(ylim = c(-1, 1)) +
      labs(title = "Charge dependence of the elution shift, per metabolite",
           subtitle = "Spearman rho of the centre-of-mass shift against charge, all tested proteins.\nDotted lines mark |rho| = 0.15; inside that band charge explains <2.3% of the variance.",
           x = NULL, y = "Spearman rho", fill = NULL) +
      theme_bw() + theme(legend.position = "top")
    .fo <- file.path(out, "rho_overview.pdf")
    tryCatch(ggsave(.fo, g, width = 7.5, height = 5),
             error = function(e) message("   !! could not write ", basename(.fo), ": ", conditionMessage(e)))
    invisible(S)
  } else { message("Nothing computed."); invisible(NULL) }
}
