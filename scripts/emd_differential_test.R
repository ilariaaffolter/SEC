# scripts/emd_differential_test.R
# =============================================================================
# A THIRD differential test, built on the EARTH MOVER'S DISTANCE (EMD, = 1-Wasserstein distance)
# between the ctrl and treatment elution profiles - a companion to the |best_lag| cross-correlation
# test in ccf_differential_test.R. Same FDR-controlled permutation framework (stratified null,
# min_intensity filter, BH q-value); only the STATISTIC differs:
#
#   |best_lag| (ccf_differential_test.R) = how many fractions you slide the WHOLE profile to line it
#       up best. Integer, captures rigid translation only; blind to partial / bimodal shifts (e.g. 40%
#       of a protein assembling into a complex while 60% stays monomeric -> the dominant peak does not
#       move -> best_lag ~ 0).
#   EMD  (this script)                   = the minimum "work" (mass x distance) to morph one normalised
#       profile into the other. Continuous, in fraction units; sees HOW MUCH mass moved HOW FAR, so it
#       catches those partial / multi-peak redistributions a single lag misses. In 1D it is just
#       sum |cumsum(p) - cumsum(q)| after normalising each profile to sum 1 (the area between the CDFs).
#       Note: normalising removes abundance, so EMD is orthogonal to CCprofiler's abundance test - it
#       asks only WHERE the mass sits, not HOW MUCH there is.
#
# It then makes a THREE-WAY comparison of the hit sets:
#   EMD-FDR (here)  vs  CCF-FDR / best_lag (ccf_differential_test.R)  vs  CCprofiler limma/beta-reg
#       (protein_DiffExprProtein). You get the overlap counts, an UpSet, the per-protein UniProt ids with
#       which method(s) flagged them, and scatters showing what EMD adds. plot_emd_hits() then draws the
#       actual SEC traces of any category so you can eyeball them.
#
# TEST (per metabolite, per protein):
#   statistic = EMD between the ctrl and treatment group-mean profiles (replicates averaged per group).
#   null      = permute the condition labels across the replicate profiles, recompute EMD per split.
#   p-value   = fraction of the permutation null (>= observed EMD), POOLED within an abundance bin
#               (int_ctrl + int_treat), same stratified scheme as the CCF test. EMD is continuous, so the
#               p-value is a direct empirical tail probability (no integer-level lookup).
#   q-value   = Benjamini-Hochberg over the tested proteins. Hit = q < fdr AND EMD >= emd_threshold.
#
# INPUT : each metabolite's *_for_plotting.RData (protein_traces_list + design_matrix +
#         protein_DiffExprProtein), and - for the three-way - tables/ccf_fdr_results.txt written by
#         ccf_differential_test.R. Run ccf_differential_test() first (with the SAME min_intensity / n_bins)
#         so the comparison is apples-to-apples; without that file you still get EMD vs CCprofiler.
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "emd_differential_test.R"))
#   emd_differential_test()                                   # all metabolites
#   emd_differential_test("PEP", min_intensity = 300, n_bins = 4)
#   emd_differential_test("PEP", emd_threshold = 1)           # also require >= ~1 fraction of mass movement
#   plot_emd_hits("PEP")                                      # trace plots of what EMD adds (see below)
#
# OUTPUT (per metabolite):
#   tables/emd_fdr_results.txt              per-protein emd, bin, pval, qval, emd_hit
#   tables/emd_vs_ccf_vs_stat_overlap.txt   per-protein EMD / CCF / stat hit flags + UniProt id + category
#   figures/emd_three_way_upset.pdf         overlap of the three hit sets
#   figures/emd_vs_bestlag_scatter.pdf      EMD vs |best_lag|, coloured by category (what EMD adds)
#   figures/emd_fdr_vs_stat_scatter.pdf     -log10(stat p) vs -log10(EMD q)
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })

# traces_obj$traces -> numeric matrix (proteins x fractions), integer-named fraction columns
.get_mat <- function(traces_obj) {
  dt <- as.data.table(traces_obj$traces)
  fc <- grep("^[0-9]+$", colnames(dt), value = TRUE); fc <- fc[order(as.numeric(fc))]
  m  <- as.matrix(dt[, ..fc]); rownames(m) <- as.character(dt$id); m
}

# 1D Earth Mover's / Wasserstein-1 distance between two elution profiles over equally-spaced fractions:
# normalise each to a distribution (sum 1) and integrate |CDF_x - CDF_y|, which in 1D is
# sum |cumsum(px) - cumsum(py)| in fraction units. NA -> 0; tiny negatives (which should not occur in
# these intensity traces) are clamped to 0 so each profile is a valid distribution. Returns NA if a
# profile carries no signal (dropped by the finite filter, like the CCF test's sd == 0 case).
.emd <- function(x, y) {
  x[is.na(x)] <- 0; y[is.na(y)] <- 0
  x[x < 0]    <- 0; y[y < 0]    <- 0
  sx <- sum(x); sy <- sum(y)
  if (sx <= 0 || sy <= 0) return(NA_real_)
  sum(abs(cumsum(x / sx) - cumsum(y / sy)))
}
# per-protein EMD for two group-mean matrices (proteins x fractions)
.emd_rows <- function(A, B) vapply(seq_len(nrow(A)), function(i) .emd(A[i, ], B[i, ]), numeric(1))

# integer abundance-bin labels (1..k) by quantiles of x; collapses to 1 bin if k<=1 or x can't be split
# (e.g. too many ties at the quantile edges). Quantile bins are invariant to monotone rescaling of x.
.make_bins <- function(x, k) {
  if (k <= 1L) return(rep(1L, length(x)))
  br <- unique(stats::quantile(x, probs = seq(0, 1, length.out = k + 1L), na.rm = TRUE))
  if (length(br) < 3L) return(rep(1L, length(x)))                 # need >= 2 usable bins
  as.integer(cut(x, breaks = br, include.lowest = TRUE, labels = FALSE))
}

# empirical upper-tail p-value of each obs against a (continuous) null vector, with the +1 correction.
# ge = sum(null >= obs), computed ties-exact via findInterval on the sorted null (left.open counts
# strictly-smaller elements): p = (1 + ge) / (1 + N).
.emp_p <- function(obs, null) {
  N <- length(null)
  if (N == 0L) return(rep(NA_real_, length(obs)))
  ge <- N - findInterval(obs, sort(null), left.open = TRUE)
  (1 + ge) / (1 + N)
}

emd_differential_test <- function(metabolites   = NULL,
                                  min_intensity  = 0,     # require > this summed intensity in BOTH conditions
                                  n_bins         = 4L,    # abundance strata for the permutation null (1 = pooled)
                                  min_null_per_bin = 100L,# a thinner bin falls back to the pooled null
                                  emd_threshold  = 0,     # a hit also needs EMD >= this (fraction units); 0 = q only
                                  fdr            = 0.05,
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

    # observed EMD + intensity filter
    A_obs   <- grp_mean(which(!is_treat)); B_obs <- grp_mean(which(is_treat))
    emd_obs <- .emd_rows(A_obs, B_obs)
    int_ctrl <- rowSums(A_obs, na.rm = TRUE); int_treat <- rowSums(B_obs, na.rm = TRUE)
    keep <- is.finite(emd_obs) & int_ctrl > min_intensity & int_treat > min_intensity
    if (!any(keep)) { message("[", m, "] no proteins pass the intensity/validity filter; skipping."); next }

    # permutation null (exclude the observed label split)
    n <- length(samples); nt <- sum(is_treat)
    all_splits <- combn(n, nt)
    obs_col    <- apply(all_splits, 2, function(col) setequal(col, which(is_treat)))
    splits     <- all_splits[, !obs_col, drop = FALSE]
    if (ncol(splits) > n_perm_max) splits <- splits[, sample(ncol(splits), n_perm_max), drop = FALSE]

    # EMD under every label permutation, for the kept proteins (rows) x splits (cols)
    kept_idx <- which(keep)
    message("[", m, "] permutation EMD test: ", length(kept_idx), " proteins x ", ncol(splits), " label permutations ...")
    null_mat <- vapply(seq_len(ncol(splits)), function(k) {
      ti <- splits[, k]; ci <- setdiff(seq_len(n), ti)
      .emd_rows(grp_mean(ci), grp_mean(ti))[kept_idx]
    }, numeric(length(kept_idx)))
    if (is.null(dim(null_mat))) null_mat <- matrix(null_mat, nrow = length(kept_idx))   # 1-split guard

    # STRATIFIED null by abundance bin (see ccf_differential_test.R for the full rationale): each protein
    # is scored against the null built from proteins of similar total signal, so exchangeability of the
    # null EMD only has to hold WITHIN a bin. min_intensity has already removed the noisiest proteins; a
    # bin too thin for a stable null falls back to the pooled null. EMD is continuous -> empirical tail p.
    bin_stat  <- (int_ctrl + int_treat)[kept_idx]
    bin_kept  <- .make_bins(bin_stat, n_bins)
    nbin      <- max(bin_kept)
    null_all  <- null_mat[is.finite(null_mat)]                     # pooled fallback for thin bins
    emd_obs_kept <- emd_obs[kept_idx]
    pval_kept <- rep(NA_real_, length(kept_idx))
    bin_info  <- vector("list", nbin)
    for (b in seq_len(nbin)) {
      rows <- which(bin_kept == b)
      nb   <- as.vector(null_mat[rows, , drop = FALSE]); nb <- nb[is.finite(nb)]
      use_pool <- length(nb) < min_null_per_bin
      nvec <- if (use_pool) null_all else nb
      pval_kept[rows] <- .emp_p(emd_obs_kept[rows], nvec)
      bin_info[[b]] <- data.table(bin = b, n_proteins = length(rows),
                                  min_signal = round(min(bin_stat[rows]), 1),
                                  max_signal = round(max(bin_stat[rows]), 1),
                                  null_n = length(nvec), pooled_fallback = use_pool)
    }
    if (nbin > 1L) { message("[", m, "] stratified null across ", nbin, " abundance bin(s):"); print(rbindlist(bin_info)) }

    pval    <- rep(NA_real_,    length(emd_obs)); pval[kept_idx]    <- pval_kept
    qval    <- rep(NA_real_,    length(emd_obs)); qval[kept_idx]    <- p.adjust(pval_kept, method = "BH")
    bin_col <- rep(NA_integer_, length(emd_obs)); bin_col[kept_idx] <- bin_kept

    res <- data.table(protein_id = common, emd = emd_obs,
                      int_ctrl = int_ctrl, int_treat = int_treat, bin = bin_col, pval = pval, qval = qval)
    res[, emd_hit := !is.na(qval) & qval < fdr & emd >= emd_threshold]
    tab_dir <- here("output", paste0("PCM_ctrl_vs_", m), "tables")
    fig_dir <- here("output", paste0("PCM_ctrl_vs_", m), "figures")
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(res[order(qval, -emd)], file.path(tab_dir, "emd_fdr_results.txt"), sep = "\t")
    message("[", m, "] EMD-FDR hits (q < ", fdr, if (emd_threshold > 0) paste0(", EMD >= ", emd_threshold) else "", "): ",
            sum(res$emd_hit, na.rm = TRUE), " / ", sum(keep), " tested.")

    # ---- three-way overlap: EMD vs CCF-FDR (best_lag) vs CCprofiler limma/beta-reg ----
    # CCprofiler statistical hits (protein_DiffExprProtein): best BH p per protein
    stat_best <- NULL
    if ("protein_DiffExprProtein" %in% ls(e)) {
      stat  <- as.data.table(e$protein_DiffExprProtein)
      pcol  <- intersect(c("pBHadj", "qVal", "pVal", "global_pBHadj"), names(stat))[1]
      idcol <- intersect(c("protein_id", "feature_id"), names(stat))[1]
      if (!is.na(pcol) && !is.na(idcol)) {
        stat_best <- stat[, .(stat_p = min(get(pcol), na.rm = TRUE)), by = c(idcol)]
        data.table::setnames(stat_best, idcol, "protein_id")
      }
    }
    if (is.null(stat_best)) message("[", m, "] no CCprofiler stat p-values in the file - stat column will be NA.")

    # CCF-FDR (best_lag) hits, from ccf_differential_test.R
    ccf_file <- file.path(tab_dir, "ccf_fdr_results.txt"); ccf <- NULL
    if (file.exists(ccf_file)) {
      ccf   <- fread(ccf_file)
      keepc <- intersect(c("protein_id", "best_lag", "ccf_hit"), names(ccf))
      ccf   <- ccf[, ..keepc]
      if ("ccf_hit" %in% names(ccf)) ccf[, ccf_hit := as.logical(ccf_hit)]
    } else {
      message("[", m, "] tables/ccf_fdr_results.txt not found - run ccf_differential_test('", m,
              "') first for the full three-way. Doing EMD vs CCprofiler only.")
    }

    mg <- res[, .(protein_id, emd, qval, emd_hit)]
    if (!is.null(stat_best)) mg <- merge(mg, stat_best, by = "protein_id", all = TRUE) else mg[, stat_p := NA_real_]
    if (!is.null(ccf))       mg <- merge(mg, ccf,       by = "protein_id", all = TRUE) else mg[, `:=`(best_lag = NA_real_, ccf_hit = NA)]
    mg[, emd_hit  := !is.na(emd_hit) & emd_hit]
    mg[, ccf_hit  := !is.na(ccf_hit) & ccf_hit]
    mg[, stat_hit := !is.na(stat_p)  & stat_p < 0.05]
    lab <- function(e_, c_, s_) { p <- c("EMD", "CCF", "stat")[c(e_, c_, s_)]; if (!length(p)) "none" else paste(p, collapse = "+") }
    mg[, category := mapply(lab, emd_hit, ccf_hit, stat_hit)]
    setcolorder(mg, intersect(c("protein_id","category","emd","qval","best_lag","stat_p","emd_hit","ccf_hit","stat_hit"), names(mg)))
    fwrite(mg[order(-emd_hit, -ccf_hit, -stat_hit, qval)], file.path(tab_dir, "emd_vs_ccf_vs_stat_overlap.txt"), sep = "\t")

    # overlap counts (this is the "number of hits overlapping"; the file above is the per-protein id list)
    nE <- sum(mg$emd_hit); nC <- sum(mg$ccf_hit); nS <- sum(mg$stat_hit)
    message(sprintf("[%s] hits - EMD:%d  CCF/best_lag:%d  CCprofiler:%d", m, nE, nC, nS))
    message(sprintf("[%s]   EMD-only (neither other): %d | EMD not in CCprofiler: %d | EMD not in CCF: %d | all three: %d",
                    m, sum(mg$emd_hit & !mg$ccf_hit & !mg$stat_hit),
                    sum(mg$emd_hit & !mg$stat_hit), sum(mg$emd_hit & !mg$ccf_hit),
                    sum(mg$emd_hit & mg$ccf_hit & mg$stat_hit)))
    print(table(mg$category))

    # UpSet of the three hit sets. UpSetR errors ("undefined columns selected") on an all-empty set, so
    # pass only the non-empty ones and plot only when >= 2 remain (an intersection needs >= 2 populated sets).
    sets_df  <- data.frame(protein_id = mg$protein_id,
                           EMD = as.integer(mg$emd_hit), CCF = as.integer(mg$ccf_hit), stat = as.integer(mg$stat_hit))
    nonempty <- c("EMD", "CCF", "stat")[c(nE > 0, nC > 0, nS > 0)]
    if (requireNamespace("UpSetR", quietly = TRUE) && length(nonempty) >= 2) {
      grDevices::pdf(file.path(fig_dir, "emd_three_way_upset.pdf"), width = 6.5, height = 4.5)
      tryCatch(print(UpSetR::upset(sets_df, sets = nonempty, order.by = "freq",
                                   mainbar.y.label = "proteins", sets.x.label = "flagged proteins")),
               error = function(err) message("[", m, "] UpSet skipped: ", conditionMessage(err)))
      grDevices::dev.off()
    } else {
      message("[", m, "] UpSet skipped - need >= 2 non-empty hit sets (have ", length(nonempty),
              "); see the overlap table and scatters instead.")
    }

    # scatter: EMD vs |best_lag| - the high-EMD / low-lag corner is exactly what the continuous EMD adds
    # over the discrete lag (partial / bimodal shifts). Only drawn when the CCF file is present.
    if (!is.null(ccf)) {
      ps <- mg[!is.na(emd) & !is.na(best_lag)]
      if (nrow(ps)) {
        g1 <- ggplot(ps, aes(emd, abs(best_lag), colour = category)) +
          geom_jitter(width = 0, height = 0.15, alpha = 0.6, size = 1.1) +
          labs(title = paste0("EMD vs |best_lag| - PCM_ctrl_vs_", m),
               subtitle = "high EMD + low |best_lag| = partial / bimodal shifts EMD catches that the lag misses",
               x = "EMD (fraction units)", y = "|best_lag| (fractions)", colour = NULL) +
          theme_bw() + theme(legend.position = "bottom")
        ggsave(file.path(fig_dir, "emd_vs_bestlag_scatter.pdf"), g1, width = 7, height = 5.5)
      }
    }
    # scatter: CCprofiler significance (x) vs EMD-FDR significance (y)
    ps2 <- mg[!is.na(stat_p) & !is.na(qval)]
    if (nrow(ps2)) {
      ps2[, `:=`(sx = -log10(pmax(stat_p, .Machine$double.xmin)),
                 sy = -log10(pmax(qval,   .Machine$double.xmin)))]
      g2 <- ggplot(ps2, aes(sx, sy, colour = category)) +
        geom_hline(yintercept = -log10(fdr),  linetype = 2, colour = "grey65") +
        geom_vline(xintercept = -log10(0.05), linetype = 2, colour = "grey65") +
        geom_point(alpha = 0.6, size = 1.1) +
        labs(title = paste0("EMD-FDR vs CCprofiler - PCM_ctrl_vs_", m),
             subtitle = "above horizontal = EMD-FDR significant; left of vertical = NOT stat-significant (what EMD adds)",
             x = "-log10(BH p) - CCprofiler limma/beta-reg", y = "-log10(BH q) - permutation EMD", colour = NULL) +
        theme_bw() + theme(legend.position = "bottom")
      ggsave(file.path(fig_dir, "emd_fdr_vs_stat_scatter.pdf"), g2, width = 7, height = 5.5)
    }
  }
  invisible(NULL)
}

# -----------------------------------------------------------------------------------------------
# VISUAL CHECK: draw the SEC traces (plot_protein_traces) of a chosen hit category, one protein per PDF
# page, so you can eyeball what EMD flags. Reads tables/emd_vs_ccf_vs_stat_overlap.txt (written above),
# so run emd_differential_test() for the metabolite first.
#
#   plot_emd_hits("PEP")                                  # default: emd_not_ccf + emd_only
#   plot_emd_hits("PEP", categories = "emd_not_stat")     # what EMD adds over CCprofiler
#   plot_emd_hits("PEP", categories = "all_three")        # highest-confidence hits
#   plot_emd_hits("PEP", max_per_cat = 100)               # raise the per-category page cap
# Categories:
#   emd_only      EMD flags, neither CCF nor CCprofiler do
#   emd_not_ccf   EMD flags, best_lag CCF does not   (partial / bimodal shifts the lag misses)
#   emd_not_stat  EMD flags, CCprofiler does not     (what EMD adds over the abundance test)
#   all_three     flagged by EMD, CCF and CCprofiler
# OUTPUT: output/PCM_ctrl_vs_<m>/emd_hit_traces/<m>_<category>_traces.pdf
# -----------------------------------------------------------------------------------------------
plot_emd_hits <- function(metabolites = NULL,
                          categories  = c("emd_not_ccf", "emd_only"),
                          max_per_cat = 40,
                          x_axis      = c("fraction", "mw"),
                          aggregate   = c("condition", "replicate"),
                          out_subdir  = "emd_hit_traces") {
  x_axis     <- match.arg(x_axis)
  aggregate  <- match.arg(aggregate)
  categories <- match.arg(categories, choices = c("emd_only", "emd_not_ccf", "emd_not_stat", "all_three"),
                          several.ok = TRUE)
  source(here::here("scripts", "plot_protein_traces.R"))   # reuse the exact trace-panel builder

  if (is.null(metabolites)) {
    dirs        <- basename(list.dirs(here("output"), recursive = FALSE))
    metabolites <- sub("^PCM_ctrl_vs_", "", dirs[grepl("^PCM_ctrl_vs_", dirs)])
  }
  if (!length(metabolites)) stop("No PCM_ctrl_vs_* output folders found.")

  for (m in metabolites) {
    f <- here("output", paste0("PCM_ctrl_vs_", m), "tables", "emd_vs_ccf_vs_stat_overlap.txt")
    if (!file.exists(f)) {
      message("[", m, "] no emd_vs_ccf_vs_stat_overlap.txt - run emd_differential_test('", m, "') first; skipping."); next
    }
    tab <- fread(f)
    for (col in c("emd_hit", "ccf_hit", "stat_hit")) if (col %in% names(tab)) tab[[col]] <- as.logical(tab[[col]])
    if (!"emd_hit" %in% names(tab)) { message("[", m, "] table lacks emd_hit; skipping."); next }
    E <- tab$emd_hit %in% TRUE
    C <- if ("ccf_hit"  %in% names(tab)) tab$ccf_hit  %in% TRUE else rep(FALSE, nrow(tab))
    S <- if ("stat_hit" %in% names(tab)) tab$stat_hit %in% TRUE else rep(FALSE, nrow(tab))
    sets <- list(emd_only     = tab$protein_id[E & !C & !S],
                 emd_not_ccf  = tab$protein_id[E & !C],
                 emd_not_stat = tab$protein_id[E & !S],
                 all_three    = tab$protein_id[E &  C &  S])

    for (cat in categories) {
      ids <- unique(as.character(sets[[cat]])); ids <- ids[!is.na(ids) & nzchar(ids)]
      if (!length(ids)) { message("[", m, "] no '", cat, "' proteins."); next }
      if (length(ids) > max_per_cat) {
        message("[", m, "] '", cat, "': ", length(ids), " proteins - plotting the first ", max_per_cat,
                " (raise max_per_cat).")
        ids <- head(ids, max_per_cat)
      }
      outdir <- here("output", paste0("PCM_ctrl_vs_", m), out_subdir)
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      pdf_file <- file.path(outdir, paste0(m, "_", cat, "_traces.pdf"))
      grDevices::pdf(pdf_file, width = 6, height = 6); np <- 0L
      for (pid in ids) {
        ok <- tryCatch({
          g <- plot_protein_traces(pid, metabolites = m, x_axis = x_axis, aggregate = aggregate,
                                   save_pdf = FALSE, print_plot = FALSE)
          if (!is.null(g)) print(g)
          !is.null(g)
        }, error = function(err) { message("   ", pid, ": ", conditionMessage(err)); FALSE })
        if (isTRUE(ok)) np <- np + 1L
      }
      grDevices::dev.off()
      message("[", m, "] '", cat, "': ", np, " protein page(s) -> ", pdf_file)
    }
  }
  invisible(NULL)
}
