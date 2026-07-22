# scripts/diagnose_power_ceiling.R
# =============================================================================
# Two small diagnostics that read the results tables already written by ccf_differential_test.R and
# emd_differential_test.R (nothing is recomputed):
#
#   diagnose_power_ceiling()  - per metabolite, per test, the MINIMUM ACHIEVABLE q-value and why. This
#       answers "was 0 hits inevitable, or did the data just fall short?". min_q is literally the smallest
#       fdr at which ANY protein would be called (BH q is monotone in p), so if min_q = 0.21 nothing can be
#       significant at 0.05 no matter the signal. The companion columns explain the ceiling:
#         n_distinct_p       how many distinct p-values the statistic can produce. |best_lag| is integer
#                            0..lag_max, so this is tiny (<= 6) -> massive ties -> BH is brutal. EMD is
#                            continuous, so it is ~n_tested -> BH can actually rank proteins.
#         min_p              smallest p the permutation null resolution allowed the data to reach.
#         n_needed_at_floor  ceil(min_p * n_tested / fdr): the MINIMUM number of proteins that must sit at
#                            the null floor for BH to reject ANY of them. If that many can't realistically
#                            reach the floor (e.g. 478 proteins all at |lag| = max), 0 hits is structural.
#       CCprofiler is shown as a reference row (its pBHadj is already adjusted; the floor columns are NA).
#
#   plot_top_q()              - plot the SEC traces (plot_protein_traces) of the top-N best-q proteins of a
#       test, one protein per PDF page, so you can eyeball the strongest candidates even when none clear
#       the FDR. Ties at the best q are broken by the effect size (|best_lag| / EMD), so you see the most-
#       shifted of the best-q group first. Also writes a small CSV listing those proteins and their q.
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "diagnose_power_ceiling.R"))
#   diagnose_power_ceiling()                       # all metabolites, all tests found
#   diagnose_power_ceiling("NAD")                  # one metabolite
#   diagnose_power_ceiling("NAD", fdr = 0.1)       # judge the ceiling against a different fdr
#   venn_top_candidates("NAD")                     # Venn of the top-10 candidates from CCF / EMD / CCprofiler
#   venn_top_candidates("NAD", n_top = 20)         # top-20 instead
#   plot_top_q("NAD", test = "emd")                # top 5 EMD proteins by q
#   plot_top_q("NAD", test = "ccf", n_top = 10)    # top 10 |best_lag| proteins by q
#   plot_top_q("NAD", test = "stat")               # top 5 CCprofiler proteins by adjusted p
#
# OUTPUT (plot_top_q, per metabolite): output/PCM_ctrl_vs_<m>/top_q_traces/
#   <m>_<test>_top<N>_traces.pdf   one page per protein
#   <m>_<test>_top<N>.csv          protein_id, q, effect size, rank
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table) })

.discover_metabolites <- function() {
  dirs <- basename(list.dirs(here("output"), recursive = FALSE))
  sub("^PCM_ctrl_vs_", "", dirs[grepl("^PCM_ctrl_vs_", dirs)])
}
.read_res <- function(path) if (file.exists(path)) data.table::fread(path) else NULL

# one summary row for a test's p/q vectors. is_perm = TRUE for our permutation tests (CCF/EMD), where the
# discreteness + floor columns are meaningful; FALSE for CCprofiler (pBHadj is already an adjusted p).
.summ_row <- function(m, test, pval, qval, fdr, is_perm = TRUE) {
  q <- qval[is.finite(qval)]; n <- length(q)
  if (!n) return(data.table(metabolite = m, test = test, n_tested = 0L, n_distinct_p = NA_integer_,
                            min_p = NA_real_, min_q = NA_real_, n_hits_q = 0L,
                            n_at_floor = NA_integer_, n_needed_at_floor = NA_integer_))
  p    <- pval[is.finite(pval)]
  minp <- min(p); minq <- min(q)
  data.table(metabolite = m, test = test, n_tested = n,
             n_distinct_p      = if (is_perm) length(unique(signif(p, 12))) else NA_integer_,
             min_p             = if (is_perm) signif(minp, 3) else NA_real_,   # raw permutation p (n/a for CCprofiler)
             min_q             = signif(minq, 3),
             n_hits_q          = sum(q < fdr),
             # proteins tied at the smallest p; this is what makes min_q differ between tests even when
             # min_p is identical, because min_q ~ min_p * n_tested / n_at_floor.
             n_at_floor        = if (is_perm) sum(p <= minp * (1 + 1e-9)) else NA_integer_,
             n_needed_at_floor = if (is_perm) as.integer(ceiling(minp * n / fdr)) else NA_integer_)
}

diagnose_power_ceiling <- function(metabolites = NULL, fdr = 0.05) {
  if (is.null(metabolites)) metabolites <- .discover_metabolites()
  if (!length(metabolites)) stop("No PCM_ctrl_vs_* output folders found.")

  rows <- list()
  for (m in metabolites) {
    tab_dir <- here("output", paste0("PCM_ctrl_vs_", m), "tables")

    ccf <- .read_res(file.path(tab_dir, "ccf_fdr_results.txt"))
    if (!is.null(ccf) && all(c("pval", "qval") %in% names(ccf)))
      rows[[paste0(m, "_CCF")]] <- .summ_row(m, "CCF (|best_lag|)", ccf$pval, ccf$qval, fdr, is_perm = TRUE)

    emd <- .read_res(file.path(tab_dir, "emd_fdr_results.txt"))
    if (!is.null(emd) && all(c("pval", "qval") %in% names(emd)))
      rows[[paste0(m, "_EMD")]] <- .summ_row(m, "EMD (Wasserstein)", emd$pval, emd$qval, fdr, is_perm = TRUE)

    ov <- .read_res(file.path(tab_dir, "emd_vs_ccf_vs_stat_overlap.txt"))
    if (is.null(ov)) ov <- .read_res(file.path(tab_dir, "ccf_fdr_vs_stat_overlap.txt"))
    if (!is.null(ov) && "stat_p" %in% names(ov))            # stat_p is CCprofiler's already-adjusted p
      rows[[paste0(m, "_stat")]] <- .summ_row(m, "CCprofiler (limma/betareg)", ov$stat_p, ov$stat_p, fdr, is_perm = FALSE)
  }
  if (!length(rows)) { message("No results tables found - run ccf_differential_test()/emd_differential_test() first."); return(invisible(NULL)) }

  S <- rbindlist(rows, use.names = TRUE)
  cat("\n==== power ceiling (min achievable q = smallest fdr at which any protein is a hit) ====\n")
  print(S)
  cat("\nHow to read it:\n",
      " - min_q is the best (smallest) q any protein reached. If min_q > your fdr, nothing is / can be a hit.\n",
      " - min_p is the permutation floor 1/(1+N): a DESIGN property (N = proteins-per-bin x label\n",
      "   permutations), so CCF and EMD share it whenever both have a protein maxing out the null - not a bug.\n",
      " - n_distinct_p: a discrete statistic (CCF: <= lag_max+1 values per bin) forces ties that BH cannot\n",
      "   separate; EMD is continuous (~n_tested distinct p) so BH can rank proteins - that is the power difference.\n",
      " - n_at_floor vs n_needed_at_floor: a hit exists only if n_at_floor >= n_needed_at_floor. n_at_floor is\n",
      "   how many proteins actually reached the smallest p (it sets min_q ~ min_p*n_tested/n_at_floor);\n",
      "   n_needed_at_floor is how many are REQUIRED. n_at_floor << n_needed_at_floor means 0 hits is structural.\n",
      " - CCprofiler row: pBHadj is already adjusted; the floor columns are NA (not a permutation test).\n", sep = "")
  invisible(S)
}

# draw a Venn of up to 3 id sets into pdf_file, using whichever Venn package is installed
# (ggVennDiagram -> VennDiagram -> eulerr). Returns TRUE if a figure was written. The overlap itself is
# always written as a CSV by the caller, so a missing Venn package only costs the picture.
.draw_venn <- function(sets, pdf_file, main) {
  sets <- sets[lengths(sets) > 0]
  if (length(sets) < 2) return(FALSE)
  if (requireNamespace("ggVennDiagram", quietly = TRUE)) {
    ok <- tryCatch({
      p <- ggVennDiagram::ggVennDiagram(sets, label = "count") + ggplot2::labs(title = main) +
           ggplot2::theme(legend.position = "none")
      ggplot2::ggsave(pdf_file, p, width = 6, height = 5.5); TRUE
    }, error = function(e) FALSE)
    if (ok) return(TRUE)
  }
  if (requireNamespace("VennDiagram", quietly = TRUE)) {
    ok <- tryCatch({
      cols <- c("#4E79A7", "#E15759", "#59A14F")[seq_along(sets)]
      g <- VennDiagram::venn.diagram(x = sets, filename = NULL, main = main, fill = cols, alpha = 0.4,
                                     disable.logging = TRUE, margin = 0.08, cex = 1.3, cat.cex = 1.1)
      grDevices::pdf(pdf_file, width = 6, height = 6); grid::grid.newpage(); grid::grid.draw(g); grDevices::dev.off(); TRUE
    }, error = function(e) FALSE)
    if (ok) return(TRUE)
  }
  if (requireNamespace("eulerr", quietly = TRUE)) {
    ok <- tryCatch({
      grDevices::pdf(pdf_file, width = 6, height = 6)
      print(plot(eulerr::euler(sets), quantities = TRUE, main = main)); grDevices::dev.off(); TRUE
    }, error = function(e) FALSE)
    if (ok) return(TRUE)
  }
  FALSE
}

# ---------------------------------------------------------------------------------------------
# Venn of the TOP-N candidate proteins (smallest q) from each of the three methods, per metabolite.
# "Candidates", not "hits": with 0 FDR hits this still shows whether the methods AGREE on their best
# proteins. Writes the Venn PDF (if a Venn package is present) and always the membership CSV.
#
#   venn_top_candidates("NAD")               # top 10 from CCF / EMD / CCprofiler
#   venn_top_candidates("NAD", n_top = 20)
# OUTPUT: output/PCM_ctrl_vs_<m>/figures/top<N>_candidate_venn.pdf
#         output/PCM_ctrl_vs_<m>/tables/top<N>_candidate_overlap.csv
# ---------------------------------------------------------------------------------------------
venn_top_candidates <- function(metabolites = NULL, n_top = 10) {
  if (is.null(metabolites)) metabolites <- .discover_metabolites()
  if (!length(metabolites)) stop("No PCM_ctrl_vs_* output folders found.")

  # top-n protein ids of a results file, ordered by q asc then effect size desc
  .top_ids <- function(path, qcol, scol, n) {
    if (!file.exists(path)) return(character(0))
    d <- data.table::fread(path)
    if (!all(c("protein_id", qcol) %in% names(d))) return(character(0))
    d <- d[is.finite(get(qcol))]
    if (!nrow(d)) return(character(0))
    ord <- if (!is.na(scol) && scol %in% names(d)) order(d[[qcol]], -d[[scol]]) else order(d[[qcol]])
    unique(as.character(head(d[ord]$protein_id, n)))
  }

  for (m in metabolites) {
    tab_dir <- here("output", paste0("PCM_ctrl_vs_", m), "tables")
    fig_dir <- here("output", paste0("PCM_ctrl_vs_", m), "figures")
    ov_file <- file.path(tab_dir, "emd_vs_ccf_vs_stat_overlap.txt")
    if (!file.exists(ov_file)) ov_file <- file.path(tab_dir, "ccf_fdr_vs_stat_overlap.txt")

    sets <- list(
      CCF        = .top_ids(file.path(tab_dir, "ccf_fdr_results.txt"), "qval",   "abs_lag", n_top),
      EMD        = .top_ids(file.path(tab_dir, "emd_fdr_results.txt"), "qval",   "emd",     n_top),
      CCprofiler = .top_ids(ov_file,                                   "stat_p", NA,        n_top))
    have <- names(sets)[lengths(sets) > 0]
    if (length(have) < 2) {
      message("[", m, "] fewer than 2 methods have results (", paste(have, collapse = ", "),
              ") - run the tests first; skipping."); next
    }

    # membership CSV (always) + printed overlap counts
    allids <- sort(unique(unlist(sets)))
    memb   <- data.table(protein_id = allids)
    for (nm in names(sets)) memb[[nm]] <- allids %in% sets[[nm]]
    memb[, n_methods := rowSums(as.matrix(memb[, names(sets), with = FALSE]))]
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(memb[order(-n_methods, protein_id)], file.path(tab_dir, paste0("top", n_top, "_candidate_overlap.csv")))

    message("[", m, "] top-", n_top, " candidate overlap:")
    message("   sizes: ", paste(sprintf("%s=%d", names(sets), lengths(sets)), collapse = "  "))
    if (all(c("CCF", "EMD") %in% have))        message("   CCF n EMD:        ", length(intersect(sets$CCF, sets$EMD)))
    if (all(c("CCF", "CCprofiler") %in% have)) message("   CCF n CCprofiler: ", length(intersect(sets$CCF, sets$CCprofiler)))
    if (all(c("EMD", "CCprofiler") %in% have)) message("   EMD n CCprofiler: ", length(intersect(sets$EMD, sets$CCprofiler)))
    if (length(have) == 3) message("   all three:        ", length(Reduce(intersect, sets)))

    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    pdf_file <- file.path(fig_dir, paste0("top", n_top, "_candidate_venn.pdf"))
    if (.draw_venn(sets, pdf_file, paste0("Top ", n_top, " candidates - PCM_ctrl_vs_", m)))
      message("   Venn -> ", pdf_file)
    else
      message("   (no Venn package found - install one of ggVennDiagram / VennDiagram / eulerr for the figure; ",
              "overlap is in top", n_top, "_candidate_overlap.csv)")
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------------------------
# plot the SEC traces of the top-N best-q proteins of a test.
# ---------------------------------------------------------------------------------------------
plot_top_q <- function(metabolites = NULL,
                       test        = c("emd", "ccf", "stat"),
                       n_top       = 5,
                       x_axis      = c("fraction", "mw"),
                       aggregate   = c("condition", "replicate"),
                       out_subdir  = "top_q_traces") {
  test      <- match.arg(test)
  x_axis    <- match.arg(x_axis)
  aggregate <- match.arg(aggregate)
  source(here::here("scripts", "plot_protein_traces.R"))   # reuse the exact trace-panel builder

  # (results file, q column, effect-size column) per test
  spec <- switch(test,
    ccf  = list(file = "ccf_fdr_results.txt",            qcol = "qval",   scol = "abs_lag", overlap = FALSE),
    emd  = list(file = "emd_fdr_results.txt",            qcol = "qval",   scol = "emd",     overlap = FALSE),
    stat = list(file = "emd_vs_ccf_vs_stat_overlap.txt", qcol = "stat_p", scol = NA,        overlap = TRUE))

  if (is.null(metabolites)) metabolites <- .discover_metabolites()
  if (!length(metabolites)) stop("No PCM_ctrl_vs_* output folders found.")

  for (m in metabolites) {
    tab_dir <- here("output", paste0("PCM_ctrl_vs_", m), "tables")
    f <- file.path(tab_dir, spec$file)
    if (test == "stat" && !file.exists(f)) f <- file.path(tab_dir, "ccf_fdr_vs_stat_overlap.txt")  # fallback overlap
    if (!file.exists(f)) { message("[", m, "] no ", basename(f), " - run the matching test first; skipping."); next }

    tab <- fread(f)
    if (!all(c("protein_id", spec$qcol) %in% names(tab))) {
      message("[", m, "] ", basename(f), " lacks protein_id/", spec$qcol, "; skipping."); next
    }
    tab <- tab[is.finite(get(spec$qcol))]
    if (!nrow(tab)) { message("[", m, "] no finite ", spec$qcol, " values; skipping."); next }

    # order by q ascending; break ties at the best q by the effect size (most-shifted first)
    ord <- if (!is.na(spec$scol) && spec$scol %in% names(tab)) order(tab[[spec$qcol]], -tab[[spec$scol]]) else order(tab[[spec$qcol]])
    top <- head(tab[ord], n_top)
    top_ids <- unique(as.character(top$protein_id)); top_ids <- top_ids[!is.na(top_ids) & nzchar(top_ids)]
    if (!length(top_ids)) { message("[", m, "] no usable protein_id in the top rows; skipping."); next }

    # console + CSV listing so each PDF page can be matched to its q
    keep_cols <- intersect(c("protein_id", spec$qcol, spec$scol), names(top))
    listing   <- top[, ..keep_cols][, rank := seq_len(.N)]
    message("[", m, "] top ", nrow(listing), " by ", test, " q:")
    print(listing)
    outdir <- here("output", paste0("PCM_ctrl_vs_", m), out_subdir)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    fwrite(listing, file.path(outdir, paste0(m, "_", test, "_top", n_top, ".csv")))

    pdf_file <- file.path(outdir, paste0(m, "_", test, "_top", n_top, "_traces.pdf"))
    grDevices::pdf(pdf_file, width = 6, height = 6); np <- 0L
    for (pid in top_ids) {
      ok <- tryCatch({
        g <- plot_protein_traces(pid, metabolites = m, x_axis = x_axis, aggregate = aggregate,
                                 save_pdf = FALSE, print_plot = FALSE)
        if (!is.null(g)) print(g)
        !is.null(g)
      }, error = function(err) { message("   ", pid, ": ", conditionMessage(err)); FALSE })
      if (isTRUE(ok)) np <- np + 1L
    }
    grDevices::dev.off()
    message("[", m, "] ", np, " protein page(s) -> ", pdf_file)
  }
  invisible(NULL)
}
