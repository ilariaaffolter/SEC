# scripts/ccf_differential_test.R
# =============================================================================
# A PROPER, FDR-CONTROLLED cross-correlation differential test - separate from the descriptive
# threshold screen in the report's section 7 (that one just flags |best_lag| >= threshold on the mean
# profiles, with no p-value). Here the peak-shift statistic gets a permutation p-value and a BH q-value,
# so "shifted" becomes FDR-controlled and directly comparable to the CCprofiler limma/beta-regression
# hits. At the end it overlaps the CCF-FDR hits with the statistical (protein_DiffExprProtein) hits and
# plots what the CCF-FDR test ADDS over CCprofiler's limma/beta-regression statistics.
#
# TEST (per metabolite, per protein):
#   statistic  = |best_lag|, the cross-correlation peak shift between the ctrl and treatment group-mean
#                elution profiles (same argmax-lag as the report's CCF, computed here per group).
#   null       = permute the condition labels across the individual replicate profiles (all label
#                splits, or a random sample if there are many), recompute |best_lag| for each split.
#   p-value    = fraction of the permutation null (>= observed |lag|). The null is POOLED across proteins
#                to get usable resolution from few replicates, but STRATIFIED: proteins are binned by
#                abundance (int_ctrl + int_treat) and each protein is scored against the null built from
#                its OWN bin. Pooling assumes the null shift is exchangeable across proteins; that breaks
#                when noisy low-signal proteins have heavier-tailed null shifts than clean ones, so binning
#                only requires exchangeability WITHIN a bin. min_intensity removes the noisiest proteins up
#                front; n_bins calibrates the rest (set n_bins = 1 for a single SAM-style pooled null). A
#                bin too thin to give a stable null (< min_null_per_bin values) falls back to the pooled
#                null. With very few replicates the test is still conservative - a calibrated screen,
#                cross-checked with the trace viewer.
#   q-value    = Benjamini-Hochberg over ALL tested proteins. Hit = q < fdr AND |best_lag| >= threshold.
#
# INPUT : each metabolite's *_for_plotting.RData (protein_traces_list per replicate + design_matrix +
#         protein_DiffExprProtein), written by the report.
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "ccf_differential_test.R"))
#   ccf_differential_test()                       # all metabolites
#   ccf_differential_test("PEP")                  # one
#   ccf_differential_test("PEP", fdr = 0.1, shift_threshold = 1, min_intensity = 300)
#   ccf_differential_test("PEP", min_intensity = 300, n_bins = 4)  # drop noise + stratified null
#   ccf_differential_test("PEP", n_bins = 1)                       # single pooled null (SAM-style)
#
# OUTPUT (per metabolite):
#   tables/ccf_fdr_results.txt          per-protein best_lag, abundance bin, pval, qval, ccf_hit
#   tables/ccf_fdr_vs_stat_overlap.txt  per-protein CCF-FDR vs statistical hit + category
#   figures/ccf_fdr_vs_stat_upset.pdf   overlap of the two hit sets
#   figures/ccf_fdr_vs_stat_scatter.pdf -log10(stat p) vs -log10(CCF q), coloured by category
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })

# traces_obj$traces -> numeric matrix (proteins x fractions), integer-named fraction columns
.get_mat <- function(traces_obj) {
  dt <- as.data.table(traces_obj$traces)
  fc <- grep("^[0-9]+$", colnames(dt), value = TRUE); fc <- fc[order(as.numeric(fc))]
  m  <- as.matrix(dt[, ..fc]); rownames(m) <- as.character(dt$id); m
}

# signed best cross-correlation lag between two profiles (argmax over +/- lag_max). |.| is the shift
# statistic. Uses the demeaned cross-covariance, whose argmax equals stats::ccf()'s (the sd/n scaling
# is constant across lags), so it matches the report's CCF but is fast enough for the permutations.
.best_lag <- function(x, y, lag_max = 5L) {
  x[is.na(x)] <- 0; y[is.na(y)] <- 0
  if (stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  x <- x - mean(x); y <- y - mean(y); n <- length(x); lags <- (-lag_max):lag_max
  cc <- vapply(lags, function(L) if (L >= 0) sum(x[(1L + L):n] * y[1:(n - L)])
                                 else        sum(x[1:(n + L)] * y[(1L - L):n]), numeric(1))
  lags[which.max(cc)]
}
# per-protein best lag for two group-mean matrices (proteins x fractions)
.best_lag_rows <- function(A, B, lag_max) vapply(seq_len(nrow(A)),
  function(i) .best_lag(A[i, ], B[i, ], lag_max), numeric(1))

# integer abundance-bin labels (1..k) by quantiles of x; collapses to 1 bin if k<=1 or x can't be split
# (e.g. too many ties at the quantile edges). Quantile bins are invariant to monotone rescaling of x.
.make_bins <- function(x, k) {
  if (k <= 1L) return(rep(1L, length(x)))
  br <- unique(stats::quantile(x, probs = seq(0, 1, length.out = k + 1L), na.rm = TRUE))
  if (length(br) < 3L) return(rep(1L, length(x)))                 # need >= 2 usable bins
  as.integer(cut(x, breaks = br, include.lowest = TRUE, labels = FALSE))
}

ccf_differential_test <- function(metabolites   = NULL,
                                  lag_max        = 5L,
                                  shift_threshold = 1L,   # a hit needs |best_lag| >= this AND q < fdr
                                  fdr            = 0.05,
                                  min_intensity  = 0,     # require > this summed intensity in BOTH conditions
                                  n_bins         = 4L,    # abundance strata for the permutation null (1 = pooled)
                                  min_null_per_bin = 100L,# a thinner bin falls back to the pooled null
                                  n_perm_max     = 2000L, # cap random label permutations if replicates are many
                                  seed           = 1L) {
  set.seed(seed)
  if (is.null(metabolites)) {
    dirs        <- basename(list.dirs(here("output"), recursive = FALSE))
    metabolites <- sub("^PCM_ctrl_vs_", "", dirs[grepl("^PCM_ctrl_vs_", dirs)])
  }
  if (!length(metabolites)) stop("No PCM_ctrl_vs_* output folders found.")

  for (m in metabolites) {
    fdir <- here("output", paste0("PCM_ctrl_vs_", m), "RData_for_further_plotting_and_analysis")
    f    <- list.files(fdir, pattern = "_for_plotting\\.RData$", full.names = TRUE)
    if (!length(f)) { message("[", m, "] no *_for_plotting.RData - render this metabolite first; skipping."); next }
    e <- new.env(); load(f[1], envir = e)
    if (!all(c("protein_traces_list", "design_matrix") %in% ls(e))) {
      message("[", m, "] file lacks protein_traces_list/design_matrix; skipping."); next
    }
    tl      <- e$protein_traces_list
    samples <- names(tl)
    mats    <- lapply(samples, function(s) .get_mat(tl[[s]]))
    common  <- Reduce(intersect, lapply(mats, rownames))
    if (length(common) < 10) { message("[", m, "] <10 shared proteins; skipping."); next }
    mats <- lapply(mats, function(M) { M <- M[common, , drop = FALSE]; M[is.na(M)] <- 0; M })

    # condition of each sample; control = ctrl/control/ref-named, else first factor level
    dm    <- as.data.table(e$design_matrix)
    cond  <- as.character(dm$Condition[match(samples, as.character(dm$Sample_name))])
    conds <- unique(cond)
    ctrl  <- conds[grepl("ctrl|control|ref", conds, ignore.case = TRUE)][1]
    if (is.na(ctrl)) ctrl <- if (is.factor(dm$Condition)) as.character(levels(dm$Condition))[1] else conds[1]
    is_treat <- cond != ctrl
    if (sum(is_treat) < 2 || sum(!is_treat) < 2) {
      message("[", m, "] need >= 2 replicates per condition for the permutation test; skipping."); next
    }
    grp_mean <- function(idx) Reduce(`+`, mats[idx]) / length(idx)

    # observed shift + intensity filter
    A_obs   <- grp_mean(which(!is_treat)); B_obs <- grp_mean(which(is_treat))
    obs_lag <- .best_lag_rows(A_obs, B_obs, lag_max); obs_abs <- abs(obs_lag)
    int_ctrl <- rowSums(A_obs, na.rm = TRUE); int_treat <- rowSums(B_obs, na.rm = TRUE)
    keep <- is.finite(obs_lag) & int_ctrl > min_intensity & int_treat > min_intensity
    if (!any(keep)) { message("[", m, "] no proteins pass the intensity/validity filter; skipping."); next }

    # permutation null (exclude the observed label split)
    n <- length(samples); nt <- sum(is_treat)
    all_splits <- combn(n, nt)
    obs_col    <- apply(all_splits, 2, function(col) setequal(col, which(is_treat)))
    splits     <- all_splits[, !obs_col, drop = FALSE]
    if (ncol(splits) > n_perm_max) splits <- splits[, sample(ncol(splits), n_perm_max), drop = FALSE]

    # |best_lag| under every label permutation, for the kept proteins (rows) x splits (cols)
    kept_idx <- which(keep)
    message("[", m, "] permutation CCF test: ", length(kept_idx), " proteins x ", ncol(splits), " label permutations ...")
    null_mat <- vapply(seq_len(ncol(splits)), function(k) {
      ti <- splits[, k]; ci <- setdiff(seq_len(n), ti)
      abs(.best_lag_rows(grp_mean(ci), grp_mean(ti), lag_max))[kept_idx]
    }, numeric(length(kept_idx)))
    if (is.null(dim(null_mat))) null_mat <- matrix(null_mat, nrow = length(kept_idx))   # 1-split guard

    # STRATIFIED null: bin the kept proteins by total signal (int_ctrl + int_treat) and build a SEPARATE
    # permutation null per bin, so each protein is judged against proteins of similar abundance. Pooling
    # everything assumes the null |best_lag| is exchangeable across ALL proteins; that fails when noisy
    # low-signal proteins have heavier-tailed null shifts than clean ones (making clean proteins look
    # conservative and noisy ones anti-conservative). Binning only needs exchangeability WITHIN a bin.
    # min_intensity has already dropped the noisiest proteins; a bin too thin for a stable null falls back
    # to the pooled null. |lag| is integer 0..lag_max, so we precompute a tail probability per level per bin.
    lev       <- 0:lag_max
    bin_stat  <- (int_ctrl + int_treat)[kept_idx]
    bin_kept  <- .make_bins(bin_stat, n_bins)
    nbin      <- max(bin_kept)
    null_all  <- null_mat[is.finite(null_mat)]                     # pooled fallback for thin bins
    obs_abs_kept <- obs_abs[kept_idx]
    pval_kept <- rep(NA_real_, length(kept_idx))
    bin_info  <- vector("list", nbin)
    for (b in seq_len(nbin)) {
      rows <- which(bin_kept == b)
      nb   <- as.vector(null_mat[rows, , drop = FALSE]); nb <- nb[is.finite(nb)]
      use_pool <- length(nb) < min_null_per_bin
      nvec  <- if (use_pool) null_all else nb
      tprob <- vapply(lev, function(v) (1 + sum(nvec >= v)) / (1 + length(nvec)), numeric(1))
      names(tprob) <- as.character(lev)
      pval_kept[rows] <- tprob[as.character(obs_abs_kept[rows])]
      bin_info[[b]] <- data.table(bin = b, n_proteins = length(rows),
                                  min_signal = round(min(bin_stat[rows]), 1),
                                  max_signal = round(max(bin_stat[rows]), 1),
                                  null_n = length(nvec), pooled_fallback = use_pool)
    }
    if (nbin > 1L) { message("[", m, "] stratified null across ", nbin, " abundance bin(s):"); print(rbindlist(bin_info)) }

    pval    <- rep(NA_real_,    length(obs_lag)); pval[kept_idx]    <- pval_kept
    qval    <- rep(NA_real_,    length(obs_lag)); qval[kept_idx]    <- p.adjust(pval_kept, method = "BH")
    bin_col <- rep(NA_integer_, length(obs_lag)); bin_col[kept_idx] <- bin_kept

    res <- data.table(protein_id = common, best_lag = obs_lag, abs_lag = obs_abs,
                      int_ctrl = int_ctrl, int_treat = int_treat, bin = bin_col, pval = pval, qval = qval)
    res[, ccf_hit := !is.na(qval) & qval < fdr & abs_lag >= shift_threshold]
    tab_dir <- here("output", paste0("PCM_ctrl_vs_", m), "tables")
    fig_dir <- here("output", paste0("PCM_ctrl_vs_", m), "figures")
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(res[order(qval, -abs_lag)], file.path(tab_dir, "ccf_fdr_results.txt"), sep = "\t")
    message("[", m, "] CCF-FDR hits (q < ", fdr, ", |lag| >= ", shift_threshold, "): ", sum(res$ccf_hit, na.rm = TRUE),
            " / ", sum(keep), " tested.")

    # ---- overlap: what the CCF-FDR test ADDS over the CCprofiler limma/beta-reg statistics ----
    if (!("protein_DiffExprProtein" %in% ls(e))) {
      message("[", m, "] protein_DiffExprProtein not in the file - skipping the overlap plot."); next
    }
    stat  <- as.data.table(e$protein_DiffExprProtein)
    pcol  <- intersect(c("pBHadj", "qVal", "pVal", "global_pBHadj"), names(stat))[1]
    idcol <- intersect(c("protein_id", "feature_id"), names(stat))[1]
    if (is.na(pcol) || is.na(idcol)) { message("[", m, "] no stat p-value/id column - skipping overlap."); next }
    stat_best <- stat[, .(stat_p = min(get(pcol), na.rm = TRUE)), by = c(idcol)]
    data.table::setnames(stat_best, idcol, "protein_id")
    mg <- merge(res, stat_best, by = "protein_id", all = TRUE)
    mg[, ccf_hit  := !is.na(ccf_hit) & ccf_hit]
    mg[, stat_hit := !is.na(stat_p)  & stat_p < 0.05]
    mg[, category := data.table::fcase(
      ccf_hit &  stat_hit, "both",
      ccf_hit & !stat_hit, "CCF-FDR only (adds)",
      !ccf_hit & stat_hit, "stat only",
      default =            "neither")]
    fwrite(mg, file.path(tab_dir, "ccf_fdr_vs_stat_overlap.txt"), sep = "\t")
    print(table(CCF_FDR = mg$ccf_hit, stat = mg$stat_hit, useNA = "ifany"))
    message("[", m, "] CCF-FDR adds ", sum(mg$category == "CCF-FDR only (adds)"),
            " protein(s) the stat test misses; ", sum(mg$category == "both"), " shared, ",
            sum(mg$category == "stat only"), " stat-only.")

    # UpSet of the two hit sets
    sets_df <- data.frame(protein_id = mg$protein_id,
                          CCF_FDR = as.integer(mg$ccf_hit), stat = as.integer(mg$stat_hit))
    if (requireNamespace("UpSetR", quietly = TRUE) && (sum(sets_df$CCF_FDR) + sum(sets_df$stat) > 0)) {
      grDevices::pdf(file.path(fig_dir, "ccf_fdr_vs_stat_upset.pdf"), width = 6, height = 4)
      print(UpSetR::upset(sets_df, sets = c("CCF_FDR", "stat"), order.by = "freq",
                          mainbar.y.label = "proteins", sets.x.label = "flagged proteins"))
      grDevices::dev.off()
    }
    # scatter: statistical significance (x) vs CCF-FDR significance (y)
    plt <- mg[!is.na(stat_p) & !is.na(qval)]
    if (nrow(plt)) {
      plt[, `:=`(sx = -log10(pmax(stat_p, .Machine$double.xmin)),
                 sy = -log10(pmax(qval,   .Machine$double.xmin)))]
      g <- ggplot(plt, aes(sx, sy, colour = category)) +
        geom_hline(yintercept = -log10(fdr),  linetype = 2, colour = "grey65") +
        geom_vline(xintercept = -log10(0.05), linetype = 2, colour = "grey65") +
        geom_point(alpha = 0.6, size = 1.1) +
        labs(title = paste0("CCF-FDR vs statistical test - ", "PCM_ctrl_vs_", m),
             subtitle = "above the horizontal line = CCF-FDR significant; left of vertical = NOT stat-significant (what CCF adds)",
             x = "-log10(BH p) - CCprofiler limma/beta-reg", y = "-log10(BH q) - permutation CCF", colour = NULL) +
        theme_bw() + theme(legend.position = "bottom")
      ggsave(file.path(fig_dir, "ccf_fdr_vs_stat_scatter.pdf"), g, width = 7, height = 5.5)
    }
  }
  invisible(NULL)
}
