# scripts/proteoform_analysis.R
# =============================================================================
# PROTEOFORM / DOMAIN-RESOLVED elution analysis, run ON TOP of a rendered report (reads the cached
# *_for_plotting.RData per metabolite - nothing heavy is recomputed).
#
# Idea: a protein's peptides should co-elute if it behaves as one species. When a CONTIGUOUS stretch of
# the sequence elutes differently - a shoulder, a second apex, or a bigger ctrl->treatment shift than the
# rest of the protein - that is the signature of a proteoform (truncation / alternative form) or a domain
# moving on its own (e.g. the Lon N-terminal domain under aKG; the two apices of P09372 under NAD).
#
# The one building block is a peptide -> residue map: each peptide's stripped sequence is located in the
# protein's UniProt sequence, giving its (start, end, midpoint) residue. Then, PER PROTEIN (one test per
# protein - peptides are the data WITHIN a test, NOT separate hypotheses, so this does not inflate the
# multiple-testing burden; it sharpens each test):
#   B1 (static proteoform)  : is a peptide's APEX fraction associated with its residue position?
#   B2 (differential shift)  : is a peptide's ctrl->treatment SHIFT (best_lag / EMD) associated with its
#                              residue position - i.e. does a contiguous region move while the rest stays?
# Association is scored by the best CONTIGUOUS split (a two-group t-like statistic over peptides ordered
# by residue) plus a descriptive Spearman rho, each with a within-protein permutation p-value (shuffle
# which peptide carries which value; positions fixed). BH over the tested proteins -> q. Treat it as a
# RANKED SCREEN (sort by q / score, eyeball with plot_proteoform), same posture as the CCF/EMD tests.
#
# This file covers the mapping + Mode B. Mode A (domain-level beta-regression on load-normalised domain
# shares + the Lon-style per-domain chromatogram / fractional-share plots) builds on the SAME mapping and
# is added next, once you have confirmed the mapping rate on your data.
#
# INPUT (per metabolite):
#   output/PCM_ctrl_vs_<m>/RData_for_further_plotting_and_analysis/*_for_plotting.RData
#       -> pepTracesList_filtered, design_matrix, and (after re-rendering with the updated Rmd) pep_seq_map
#   protein sequences from output/PCM_ctrl_vs_<m>/rdata/uniprot.RData (falls back to the shared cache
#       output/uniprot_annotation_shared.RData). Both carry the fetched `sequence` column.
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "proteoform_analysis.R"))
#   map_peptides_to_residues("NAD")                 # just the mapping + its coverage/mismatch report
#   proteoform_scan("NAD")                          # Mode B screen -> tables/proteoform_scan.txt
#   proteoform_scan("NAD", statistic = "lag", min_intensity = 300)
#   plot_proteoform("P09372", "NAD")                # peptide traces coloured by residue position
#
# OUTPUT (per metabolite):
#   tables/proteoform_peptide_map.txt   peptide -> protein, residue start/end/midpoint, mapped flag
#   tables/proteoform_scan.txt          per-protein B1/B2 scores, permutation p, BH q, breakpoint residue
#   figures/proteoform/<protein>.pdf    (plot_proteoform) peptide traces coloured by residue position
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })

# ---- small statistics helpers (kept local so this script is standalone) ----------------------------
# traces_obj$traces -> numeric matrix (rows x fractions), integer-named fraction columns
.get_mat <- function(traces_obj) {
  dt <- as.data.table(traces_obj$traces)
  fc <- grep("^[0-9]+$", colnames(dt), value = TRUE); fc <- fc[order(as.numeric(fc))]
  m  <- as.matrix(dt[, ..fc]); rownames(m) <- as.character(dt$id); m
}
# signed best cross-correlation lag between two profiles (demeaned cross-covariance argmax; matches ccf)
.best_lag <- function(x, y, lag_max = 5L) {
  x[is.na(x)] <- 0; y[is.na(y)] <- 0
  if (stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  x <- x - mean(x); y <- y - mean(y); n <- length(x); lags <- (-lag_max):lag_max
  cc <- vapply(lags, function(L) if (L >= 0) sum(x[(1L + L):n] * y[1:(n - L)])
                                 else        sum(x[1:(n + L)] * y[(1L - L):n]), numeric(1))
  lags[which.max(cc)]
}
# 1D Earth Mover's / Wasserstein-1 distance between two profiles (normalise to a distribution, area
# between CDFs). Sign added: negative if treatment mass sits at LOWER fractions than ctrl, else positive,
# so a position-structured shift keeps a consistent direction. NA if a profile has no signal.
.emd_signed <- function(x, y) {
  x[is.na(x)] <- 0; y[is.na(y)] <- 0; x[x < 0] <- 0; y[y < 0] <- 0
  sx <- sum(x); sy <- sum(y)
  if (sx <= 0 || sy <= 0) return(NA_real_)
  px <- x / sx; py <- y / sy
  d  <- sum(abs(cumsum(px) - cumsum(py)))
  mux <- sum(seq_along(px) * px); muy <- sum(seq_along(py) * py)
  if (muy < mux) -d else d
}
# apex fraction (fraction number of the max) of an intensity vector with fraction-numbered names
.apex_frac <- function(v, fracs) { v[is.na(v)] <- 0; if (all(v == 0)) return(NA_real_); fracs[which.max(v)] }

# best contiguous two-group split of a position-ordered value vector: max |mean diff| / pooled SE over
# breakpoints. Returns the statistic and the breakpoint index k (split between position-ordered k and k+1).
# Vectorised over breakpoints (cumulative sums), so each call is O(n) - the permutation loop calls it a lot.
.seg_stat <- function(v_ordered) {
  n <- length(v_ordered)
  if (n < 4) return(list(stat = NA_real_, k = NA_integer_))
  ks  <- 2:(n - 2)
  cs  <- cumsum(v_ordered); css <- cumsum(v_ordered^2); tot <- cs[n]; totss <- css[n]
  na  <- ks; nb <- n - ks
  sa  <- cs[ks]; sb <- tot - sa
  ssa <- css[ks] - sa^2 / na                    # within-group sum of squares, left  (1:k)
  ssb <- (totss - css[ks]) - sb^2 / nb          # within-group sum of squares, right (k+1:n)
  sp  <- sqrt((ssa + ssb) / (n - 2))
  tt  <- abs(sa / na - sb / nb) / (sp * sqrt(1 / na + 1 / nb))
  tt[!is.finite(tt)] <- -Inf
  if (all(tt == -Inf)) return(list(stat = NA_real_, k = NA_integer_))
  i <- which.max(tt)
  list(stat = tt[i], k = ks[i])
}
# permutation p for the best-split statistic: shuffle values across positions (positions fixed).
.seg_perm <- function(v_ordered, n_perm) {
  o <- .seg_stat(v_ordered)
  if (!is.finite(o$stat)) return(list(stat = NA_real_, k = NA_integer_, p = NA_real_))
  ge <- 1L
  for (b in seq_len(n_perm)) { s <- .seg_stat(sample(v_ordered))$stat; if (is.finite(s) && s >= o$stat) ge <- ge + 1L }
  list(stat = o$stat, k = o$k, p = ge / (n_perm + 1L))
}

# ---- data loading ---------------------------------------------------------------------------------
.load_plotting <- function(m) {
  fdir <- here("output", paste0("PCM_ctrl_vs_", m), "RData_for_further_plotting_and_analysis")
  f    <- list.files(fdir, pattern = "_for_plotting\\.RData$", full.names = TRUE)
  if (!length(f)) return(NULL)
  e <- new.env(); load(f[1], envir = e); e
}
# protein_id -> UniProt sequence (per-metabolite uniprot.RData, else the shared cache)
.protein_seqs <- function(m) {
  grab <- function(obj) {
    d <- tryCatch(as.data.table(obj), error = function(e) NULL)
    if (is.null(d) || !all(c("input_id", "sequence") %in% names(d))) return(NULL)
    d <- d[!is.na(sequence) & nzchar(sequence)]
    setNames(as.character(d$sequence), as.character(d$input_id))
  }
  uf <- here("output", paste0("PCM_ctrl_vs_", m), "rdata", "uniprot.RData")
  if (file.exists(uf)) { e <- new.env(); load(uf, envir = e); if ("uniprot" %in% ls(e)) { s <- grab(e$uniprot); if (length(s)) return(s) } }
  sf <- here("output", "uniprot_annotation_shared.RData")
  if (file.exists(sf)) { e <- new.env(); load(sf, envir = e); if (".uniprot_all" %in% ls(e)) { s <- grab(e$.uniprot_all); if (length(s)) return(s) } }
  character(0)
}
# peptide_id -> stripped AA sequence (saved pep_seq_map, else the trace id itself, cleaned to letters)
.pep_seq_lookup <- function(e) {
  if ("pep_seq_map" %in% ls(e)) {
    d <- as.data.table(e$pep_seq_map)
    if (all(c("peptide_id", "stripped_sequence") %in% names(d)))
      return(list(map = setNames(toupper(gsub("[^A-Za-z]", "", d$stripped_sequence)), as.character(d$peptide_id)), from = "pep_seq_map"))
  }
  list(map = NULL, from = "trace_id")   # fall back to using the trace id as the sequence
}
.seq_of <- function(pep_id, lk) {
  s <- if (!is.null(lk$map) && pep_id %in% names(lk$map)) lk$map[[pep_id]] else pep_id
  toupper(gsub("[^A-Za-z]", "", s))
}

# ---- condition-mean peptide matrices --------------------------------------------------------------
.pep_condition_means <- function(e) {
  tl <- e$pepTracesList_filtered; samples <- names(tl)
  mats <- lapply(samples, function(s) .get_mat(tl[[s]]))
  common <- Reduce(intersect, lapply(mats, rownames))
  if (length(common) < 2) return(NULL)
  mats <- lapply(mats, function(M) { M <- M[common, , drop = FALSE]; M[is.na(M)] <- 0; M })
  dm   <- as.data.table(e$design_matrix)
  cond <- as.character(dm$Condition[match(samples, as.character(dm$Sample_name))])
  conds <- unique(cond)
  ctrl  <- conds[grepl("ctrl|control|ref", conds, ignore.case = TRUE)][1]
  if (is.na(ctrl)) ctrl <- if (is.factor(dm$Condition)) as.character(levels(dm$Condition))[1] else conds[1]
  is_treat <- cond != ctrl
  if (sum(is_treat) < 1 || sum(!is_treat) < 1) return(NULL)
  ann <- as.data.table(tl[[1]]$trace_annotation)
  list(ctrl  = Reduce(`+`, mats[which(!is_treat)]) / sum(!is_treat),
       treat = Reduce(`+`, mats[which(is_treat)]) / sum(is_treat),
       fracs = as.numeric(colnames(mats[[1]])),
       ann   = ann, ctrl_name = ctrl)
}

# ---- peptide -> residue mapping (standalone + report) ---------------------------------------------
map_peptides_to_residues <- function(metabolites = NULL, write = TRUE) {
  if (is.null(metabolites)) {
    dirs <- basename(list.dirs(here("output"), recursive = FALSE))
    metabolites <- sub("^PCM_ctrl_vs_", "", dirs[grepl("^PCM_ctrl_vs_", dirs)])
  }
  if (!length(metabolites)) stop("No PCM_ctrl_vs_* output folders found.")
  out <- list()
  for (m in metabolites) {
    e <- .load_plotting(m); if (is.null(e)) { message("[", m, "] no *_for_plotting.RData; skipping."); next }
    if (!"pepTracesList_filtered" %in% ls(e)) { message("[", m, "] no pepTracesList_filtered; skipping."); next }
    seqs <- .protein_seqs(m)
    if (!length(seqs)) { message("[", m, "] no UniProt sequences found (uniprot.RData / shared cache) - re-render with the updated Rmd; skipping."); next }
    lk  <- .pep_seq_lookup(e)
    ann <- as.data.table(e$pepTracesList_filtered[[1]]$trace_annotation)
    if (!all(c("id", "protein_id") %in% names(ann))) { message("[", m, "] peptide annotation lacks id/protein_id; skipping."); next }
    ann <- unique(ann[, .(id = as.character(id), protein_id = as.character(protein_id))])
    ann <- ann[!grepl("DECOY|^CON_|^iRT", protein_id)]

    map <- rbindlist(lapply(seq_len(nrow(ann)), function(i) {
      pid <- ann$protein_id[i]; pep <- ann$id[i]
      pseq <- .seq_of(pep, lk); prot <- seqs[pid]
      have_prot <- !is.na(prot) && nzchar(prot)
      st <- if (have_prot && nzchar(pseq)) as.integer(regexpr(pseq, prot, fixed = TRUE)) else -1L
      data.table(peptide_id = pep, protein_id = pid, pep_len = nchar(pseq),
                 have_protein_seq = have_prot,
                 start = if (st > 0) st else NA_integer_,
                 end   = if (st > 0) st + nchar(pseq) - 1L else NA_integer_,
                 midpoint = if (st > 0) (st + st + nchar(pseq) - 1L) / 2 else NA_real_,
                 mapped = st > 0)
    }), use.names = TRUE)

    n <- nrow(map); nmap <- sum(map$mapped); nprot <- sum(!map$have_protein_seq)
    message(sprintf("[%s] peptide->residue mapping via %s: %d/%d mapped (%.1f%%); %d peptide(s) had no protein sequence; %d unmatched despite a sequence.",
                    m, lk$from, nmap, n, 100 * nmap / max(n, 1), nprot, sum(!map$mapped & map$have_protein_seq)))
    if (lk$from == "trace_id")
      message("   NOTE: using the trace id as the peptide sequence (pep_seq_map not in the saved file). Re-render with the updated Rmd to persist the real stripped-sequence map if the mapping rate looks low.")
    if (write) {
      td <- here("output", paste0("PCM_ctrl_vs_", m), "tables"); dir.create(td, recursive = TRUE, showWarnings = FALSE)
      fwrite(map, file.path(td, "proteoform_peptide_map.txt"), sep = "\t")
    }
    out[[m]] <- map
  }
  invisible(out)
}

# ---- Mode B: per-protein position-structure screen ------------------------------------------------
proteoform_scan <- function(metabolites   = NULL,
                            statistic      = c("emd", "lag"),
                            min_intensity  = 0,      # drop peptides below this summed (ctrl+treat) intensity
                            min_peptides   = 4L,     # need at least this many mapped, usable peptides
                            lag_max        = 5L,
                            n_perm         = 2000L,
                            fdr            = 0.05,
                            seed           = 1L) {
  statistic <- match.arg(statistic); set.seed(seed)
  if (is.null(metabolites)) {
    dirs <- basename(list.dirs(here("output"), recursive = FALSE))
    metabolites <- sub("^PCM_ctrl_vs_", "", dirs[grepl("^PCM_ctrl_vs_", dirs)])
  }
  if (!length(metabolites)) stop("No PCM_ctrl_vs_* output folders found.")

  for (m in metabolites) {
    e <- .load_plotting(m); if (is.null(e)) { message("[", m, "] no *_for_plotting.RData; skipping."); next }
    if (!all(c("pepTracesList_filtered", "design_matrix") %in% ls(e))) { message("[", m, "] file lacks pepTracesList_filtered/design_matrix; skipping."); next }
    pm <- .pep_condition_means(e); if (is.null(pm)) { message("[", m, "] could not build condition-mean peptide matrices; skipping."); next }
    seqs <- .protein_seqs(m); if (!length(seqs)) { message("[", m, "] no UniProt sequences; re-render with the updated Rmd; skipping."); next }
    lk <- .pep_seq_lookup(e)

    ann <- unique(pm$ann[, .(id = as.character(id), protein_id = as.character(protein_id))])
    ann <- ann[!grepl("DECOY|^CON_|^iRT", protein_id) & id %in% rownames(pm$ctrl)]
    prot_ids <- unique(ann$protein_id)
    message("[", m, "] proteoform scan (", statistic, "): ", length(prot_ids), " proteins, ",
            n_perm, " permutations each ...")

    n_total_pep <- 0L; n_mapped_pep <- 0L
    rows <- vector("list", length(prot_ids))
    for (j in seq_along(prot_ids)) {
      pid  <- prot_ids[j]; prot <- seqs[pid]
      peps <- ann[protein_id == pid]$id
      n_total_pep <- n_total_pep + length(peps)
      if (is.na(prot) || !nzchar(prot)) next

      # map + intensity filter + per-peptide apex and shift
      recs <- lapply(peps, function(pep) {
        pseq <- .seq_of(pep, lk); if (!nzchar(pseq)) return(NULL)
        st <- as.integer(regexpr(pseq, prot, fixed = TRUE))
        if (st <= 0) return(NULL)
        cr <- pm$ctrl[pep, ]; tr <- pm$treat[pep, ]
        if (sum(cr) + sum(tr) <= min_intensity) return(NULL)
        sh <- if (statistic == "emd") .emd_signed(cr, tr) else .best_lag(cr, tr, lag_max)
        data.table(peptide_id = pep, midpoint = st + (nchar(pseq) - 1) / 2,
                   apex = .apex_frac(cr + tr, pm$fracs), shift = sh)
      })
      d <- rbindlist(recs[!vapply(recs, is.null, logical(1))], use.names = TRUE)
      if (!nrow(d)) next
      n_mapped_pep <- n_mapped_pep + nrow(d)
      d <- d[is.finite(midpoint)]
      if (nrow(d) < min_peptides) next
      setorder(d, midpoint)

      # B1 static (apex vs position) and B2 differential (shift vs position)
      ap <- d[is.finite(apex)]
      sh <- d[is.finite(shift)]
      b1 <- if (nrow(ap) >= min_peptides) .seg_perm(ap$apex, n_perm) else list(stat = NA_real_, k = NA_integer_, p = NA_real_)
      b2 <- if (nrow(sh) >= min_peptides) .seg_perm(sh$shift, n_perm) else list(stat = NA_real_, k = NA_integer_, p = NA_real_)
      rho_static <- if (nrow(ap) >= 3 && stats::sd(ap$apex)  > 0) suppressWarnings(stats::cor(ap$midpoint, ap$apex,  method = "spearman")) else NA_real_
      rho_diff   <- if (nrow(sh) >= 3 && stats::sd(sh$shift) > 0) suppressWarnings(stats::cor(sh$midpoint, sh$shift, method = "spearman")) else NA_real_
      brk <- if (is.finite(b2$k) && b2$k < nrow(sh)) round((sh$midpoint[b2$k] + sh$midpoint[b2$k + 1]) / 2) else NA_real_
      lo  <- if (is.finite(b2$k)) mean(sh$shift[seq_len(b2$k)]) else NA_real_
      hi  <- if (is.finite(b2$k)) mean(sh$shift[(b2$k + 1):nrow(sh)]) else NA_real_

      rows[[j]] <- data.table(protein_id = pid, n_pep = length(peps), n_used = nrow(d),
                              static_score = b1$stat, static_p = b1$p, rho_static = rho_static,
                              diff_score = b2$stat, diff_p = b2$p, rho_diff = rho_diff,
                              break_residue = brk, shift_lo = lo, shift_hi = hi)
    }
    res <- rbindlist(rows[!vapply(rows, is.null, logical(1))], use.names = TRUE)
    if (!nrow(res)) { message("[", m, "] no protein had >= ", min_peptides, " mapped, usable peptides; skipping."); next }
    res[, diff_q   := ifelse(is.na(diff_p),   NA_real_, p.adjust(diff_p,   method = "BH"))]
    res[, static_q := ifelse(is.na(static_p), NA_real_, p.adjust(static_p, method = "BH"))]
    setorder(res, diff_q, -diff_score, na.last = TRUE)

    td <- here("output", paste0("PCM_ctrl_vs_", m), "tables"); dir.create(td, recursive = TRUE, showWarnings = FALSE)
    fwrite(res, file.path(td, "proteoform_scan.txt"), sep = "\t")
    message(sprintf("[%s] mapping: %d/%d peptides mapped (%.1f%%). Tested %d protein(s).",
                    m, n_mapped_pep, n_total_pep, 100 * n_mapped_pep / max(n_total_pep, 1), nrow(res)))
    message(sprintf("[%s] differential (B2) q<%.2f: %d | static (B1) q<%.2f: %d  ->  tables/proteoform_scan.txt (ranked by diff_q)",
                    m, fdr, sum(res$diff_q < fdr, na.rm = TRUE), fdr, sum(res$static_q < fdr, na.rm = TRUE)))
    print(utils::head(res[, .(protein_id, n_used, diff_score, diff_p, diff_q, break_residue, static_score, static_p)], 10))
  }
  invisible(NULL)
}

# ---- diagnostic plot: peptide traces coloured by residue position ---------------------------------
plot_proteoform <- function(protein_id, metabolites, out_subdir = "proteoform", save_pdf = TRUE, print_plot = TRUE) {
  stopifnot(length(protein_id) == 1L, length(metabolites) == 1L)
  pid <- as.character(protein_id); m <- metabolites[1]
  e  <- .load_plotting(m); if (is.null(e)) stop("[", m, "] no *_for_plotting.RData.")
  pm <- .pep_condition_means(e); if (is.null(pm)) stop("[", m, "] could not build peptide matrices.")
  seqs <- .protein_seqs(m); prot <- seqs[pid]
  if (is.na(prot) || !nzchar(prot)) stop("No UniProt sequence for ", pid, " (re-render with the updated Rmd).")
  lk <- .pep_seq_lookup(e)
  ann <- unique(pm$ann[, .(id = as.character(id), protein_id = as.character(protein_id))])
  peps <- ann[protein_id == pid & id %in% rownames(pm$ctrl)]$id
  if (!length(peps)) stop("No peptide traces for ", pid, " in ", m, ".")

  long <- rbindlist(lapply(peps, function(pep) {
    pseq <- .seq_of(pep, lk); if (!nzchar(pseq)) return(NULL)
    st <- as.integer(regexpr(pseq, prot, fixed = TRUE))
    if (st <= 0) return(NULL)
    mid <- st + (nchar(pseq) - 1) / 2
    rbind(data.table(fraction = pm$fracs, intensity = pm$ctrl[pep, ],  condition = pm$ctrl_name, peptide = pep, residue = mid),
          data.table(fraction = pm$fracs, intensity = pm$treat[pep, ], condition = "treatment",  peptide = pep, residue = mid))
  }), use.names = TRUE)
  if (is.null(long) || !nrow(long)) stop("No peptide of ", pid, " mapped to its sequence.")
  long[, condition := factor(condition, levels = c("treatment", pm$ctrl_name))]

  p <- ggplot(long, aes(fraction, intensity, group = peptide, colour = residue)) +
    geom_line(linewidth = 0.5, alpha = 0.85) +
    scale_colour_viridis_c(option = "C") +
    facet_grid(condition ~ ., scales = "free_y") +
    labs(title = paste0("Proteoform view: ", pid, "  (", m, ")"),
         subtitle = paste0(length(unique(long$peptide)), " mapped peptides; colour = residue position (N->C)"),
         x = "fraction", y = "intensity", colour = "residue") +
    theme_bw()

  if (save_pdf) {
    od <- here("output", paste0("PCM_ctrl_vs_", m), "figures", out_subdir); dir.create(od, recursive = TRUE, showWarnings = FALSE)
    fn <- file.path(od, paste0(pid, ".pdf")); ggsave(fn, p, width = 7, height = 6); message("Wrote ", fn)
  }
  if (print_plot) print(p)
  invisible(p)
}
