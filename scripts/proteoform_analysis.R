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
# Mode A (domain-level, targeted - your Lon case) uses ANNOTATED domain boundaries (UniProt ft_domain, or
# a manual override table for cases UniProt lacks). It assigns each peptide to a domain, expresses each
# domain as its LOAD-NORMALISED share of the total protein signal per fraction, draws the per-domain
# chromatogram (ctrl vs treatment, mean +/- SD) and the fractional-share barplot over chosen fraction
# windows, and tests domain shares between conditions by BETA REGRESSION (the right model for a bounded
# proportion - the same family as CCprofiler's assembly test - falling back to a logit t-test).
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
#   # Mode A (targeted, domain-resolved). E. coli Lon = P0A9M0 (784 aa):
#   domain_analysis("P0A9M0", "aKG")                                       # domains from UniProt ft_domain
#   domain_analysis("P0A9M0", "aKG", windows = list(F14_16 = 14:16, F22 = 22))   # + fractional-share barplot & test
#   lon <- data.frame(protein_id = "P0A9M0", domain = c("NTD","AAA+","Protease"),
#                     start = c(1,246,588), end = c(245,587,784))          # manual boundaries (Lon)
#   domain_analysis("P0A9M0", "aKG", domains = lon, windows = list(F14_16 = 14:16, F22 = 22))
#
# OUTPUT (per metabolite):
#   tables/proteoform_peptide_map.txt   peptide -> protein, residue start/end/midpoint, mapped flag
#   tables/proteoform_scan.txt          per-protein B1/B2 scores, permutation p, BH q, breakpoint residue
#   figures/proteoform/<protein>.pdf    (plot_proteoform) peptide traces coloured by residue position
#   figures/domains/<protein>_<m>_chromatogram.pdf   (domain_analysis) per-domain load-normalised traces
#   figures/domains/<protein>_<m>_shares.pdf         (domain_analysis) fractional-share barplot over windows
#   tables/domains/<protein>_<m>_domain_test.txt     (domain_analysis) per-domain(/window) beta-reg/t test
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

# ===================================================================================================
# Mode A: domain-resolved, load-normalised elution (the targeted / Lon analysis)
# ===================================================================================================

# parse a UniProt ft_domain feature string -> data.table(domain, start, end). Handles the standard
# "DOMAIN 246..587; /note=\"AAA+\"; DOMAIN 588..784; /note=\"Peptidase S16\"" format (and <>-uncertain
# coordinates, whose numeric part is taken). Returns NULL if nothing parses.
.parse_ft_domain <- function(s) {
  if (is.na(s) || !nzchar(s)) return(NULL)
  chunks <- strsplit(s, "(?=DOMAIN )", perl = TRUE)[[1]]
  chunks <- chunks[grepl("^DOMAIN", chunks)]
  d <- rbindlist(lapply(chunks, function(ch) {
    coord <- regmatches(ch, regexpr("[0-9]+\\.\\.[0-9]+", ch))
    if (!length(coord)) return(NULL)
    se   <- as.integer(strsplit(coord, "\\.\\.")[[1]])
    note <- regmatches(ch, regexpr('/note="[^"]*"', ch))
    nm   <- if (length(note)) sub('/note="([^"]*)"', "\\1", note) else "domain"
    data.table(domain = nm, start = se[1], end = se[2])
  }), use.names = TRUE)
  if (!nrow(d)) return(NULL)
  d[order(start)]
}

# protein_id -> ft_domain string (per-metabolite uniprot.RData, else the shared cache)
.protein_domains <- function(m) {
  grab <- function(obj) {
    d <- tryCatch(as.data.table(obj), error = function(e) NULL)
    if (is.null(d) || !all(c("input_id", "ft_domain") %in% names(d))) return(NULL)
    d <- d[!is.na(ft_domain) & nzchar(ft_domain)]
    if (!nrow(d)) return(NULL)
    setNames(as.character(d$ft_domain), as.character(d$input_id))
  }
  uf <- here("output", paste0("PCM_ctrl_vs_", m), "rdata", "uniprot.RData")
  if (file.exists(uf)) { e <- new.env(); load(uf, envir = e); if ("uniprot" %in% ls(e)) { s <- grab(e$uniprot); if (length(s)) return(s) } }
  sf <- here("output", "uniprot_annotation_shared.RData")
  if (file.exists(sf)) { e <- new.env(); load(sf, envir = e); if (".uniprot_all" %in% ls(e)) { s <- grab(e$.uniprot_all); if (length(s)) return(s) } }
  character(0)
}

# resolve domains for one protein: manual override (data.frame/CSV path with protein_id,domain,start,end)
# wins; else UniProt ft_domain. Returns data.table(domain,start,end) ordered by start, or NULL.
.get_domains <- function(pid, m, override) {
  if (!is.null(override)) {
    ov <- if (is.character(override) && length(override) == 1L && file.exists(override)) fread(override) else as.data.table(override)
    if (all(c("protein_id", "domain", "start", "end") %in% names(ov))) {
      d <- ov[as.character(protein_id) == pid, .(domain = as.character(domain), start = as.integer(start), end = as.integer(end))]
      if (nrow(d)) return(d[order(start)])
    }
  }
  s <- .protein_domains(m)[pid]
  if (length(s) != 1 || is.na(s) || !nzchar(s)) return(NULL)
  .parse_ft_domain(s)
}

# per-sample peptide matrices + condition labels
.pep_samples <- function(e) {
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
  list(mats = mats, samples = samples, cond = cond, fracs = as.numeric(colnames(mats[[1]])),
       ann = as.data.table(tl[[1]]$trace_annotation), ctrl_name = ctrl)
}

# test a bounded share (in [0,1]) between two conditions: beta regression LRT (mirrors CCprofiler's
# assembly test on a proportion), falling back to a logit t-test. Returns p + method used.
.share_test <- function(share, cond) {
  cond <- factor(cond)
  if (nlevels(cond) != 2 || length(share) < 4) return(list(p = NA_real_, method = "n/a"))
  eps <- 1e-6; y <- pmin(pmax(share, eps), 1 - eps)
  if (requireNamespace("betareg", quietly = TRUE)) {
    r <- tryCatch({
      m1 <- betareg::betareg(y ~ cond); m0 <- betareg::betareg(y ~ 1)
      lr <- 2 * (as.numeric(stats::logLik(m1)) - as.numeric(stats::logLik(m0)))
      list(p = stats::pchisq(lr, df = 1, lower.tail = FALSE), method = "betareg-LRT")
    }, error = function(e) NULL)
    if (!is.null(r) && is.finite(r$p)) return(r)
  }
  tryCatch({ tt <- stats::t.test(stats::qlogis(y) ~ cond); list(p = tt$p.value, method = "logit-t") },
           error = function(e) list(p = NA_real_, method = "failed"))
}

domain_analysis <- function(protein_id,
                            metabolites,
                            domains     = NULL,     # NULL = UniProt ft_domain; else data.frame/CSV (protein_id,domain,start,end)
                            windows     = NULL,     # named list of fraction sets, e.g. list(F14_16 = 14:16, F22 = 22)
                            min_intensity = 0,
                            out_subdir  = "domains",
                            save_pdf    = TRUE,
                            print_plot  = TRUE) {
  stopifnot(length(protein_id) == 1L)
  pid <- as.character(protein_id)
  last <- NULL
  for (m in metabolites) {
    e <- .load_plotting(m); if (is.null(e)) { message("[", m, "] no *_for_plotting.RData; skipping."); next }
    ps <- .pep_samples(e);  if (is.null(ps)) { message("[", m, "] could not build peptide matrices; skipping."); next }
    seqs <- .protein_seqs(m); prot <- seqs[pid]
    if (is.na(prot) || !nzchar(prot)) { message("[", m, "] no UniProt sequence for ", pid, " (re-render with the updated Rmd); skipping."); next }
    doms <- .get_domains(pid, m, domains)
    if (is.null(doms) || !nrow(doms)) {
      message("[", m, "] no domains for ", pid, " (UniProt ft_domain empty / not fetched, and no override). ",
              "Re-render so ft_domain is fetched, or pass domains = <table>; skipping."); next
    }
    lk  <- .pep_seq_lookup(e)
    peps <- unique(ps$ann[as.character(protein_id) == pid & as.character(id) %in% rownames(ps$mats[[1]])]$id)
    if (!length(peps)) { message("[", m, "] no peptide traces for ", pid, "; skipping."); next }

    # assign each peptide to a domain by residue midpoint (else "unassigned"); computed once
    pdmap <- rbindlist(lapply(peps, function(pep) {
      pseq <- .seq_of(pep, lk); mid <- NA_real_
      if (nzchar(pseq)) { st <- as.integer(regexpr(pseq, prot, fixed = TRUE)); if (st > 0) mid <- st + (nchar(pseq) - 1) / 2 }
      dd <- if (is.finite(mid)) { hit <- doms[start <= mid & end >= mid]; if (nrow(hit)) hit$domain[1] else "unassigned" } else "unassigned"
      data.table(peptide_id = as.character(pep), domain = dd)
    }), use.names = TRUE)

    # per sample: domain share of total protein signal, per fraction (load-normalised)
    shares <- rbindlist(lapply(seq_along(ps$mats), function(si) {
      mat  <- ps$mats[[si]]; here_p <- intersect(pdmap$peptide_id, rownames(mat))
      if (!length(here_p)) return(NULL)
      sub  <- mat[here_p, , drop = FALSE]; total <- sum(sub)
      if (total <= min_intensity || total <= 0) return(NULL)
      dmp  <- pdmap[peptide_id %in% here_p]
      rbindlist(lapply(unique(dmp$domain), function(dd) {
        pep_d <- dmp[domain == dd]$peptide_id
        v <- colSums(sub[pep_d, , drop = FALSE])
        data.table(sample = ps$samples[si], condition = ps$cond[si], domain = dd,
                   fraction = ps$fracs, share = as.numeric(v) / total)
      }), use.names = TRUE)
    }), use.names = TRUE)
    if (is.null(shares) || !nrow(shares)) { message("[", m, "] ", pid, " has no signal above min_intensity; skipping."); next }

    dom_levels <- c(doms$domain, "unassigned")
    shares[, domain    := factor(domain, levels = dom_levels[dom_levels %in% unique(domain)])]
    shares[, condition := factor(condition, levels = c(setdiff(unique(condition), ps$ctrl_name), ps$ctrl_name))]  # treatment first

    tab_dir <- here("output", paste0("PCM_ctrl_vs_", m), "tables", out_subdir)
    fig_dir <- here("output", paste0("PCM_ctrl_vs_", m), "figures", out_subdir)
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

    # ---- plot 1: per-domain load-normalised chromatogram (mean +/- SD across replicates) ----
    agg <- shares[, .(mean = mean(share), sd = stats::sd(share)), by = .(domain, condition, fraction)]
    agg[is.na(sd), sd := 0]
    g1 <- ggplot(agg, aes(fraction, mean, colour = condition, fill = condition)) +
      geom_ribbon(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), alpha = 0.18, colour = NA) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~ domain, scales = "free_y") +
      labs(title = paste0("Per-domain chromatogram: ", pid, " (", m, ")"),
           subtitle = "load-normalised (share of total protein signal), mean +/- SD across replicates",
           x = "fraction", y = "share of total protein signal", colour = NULL, fill = NULL) +
      theme_bw() + theme(legend.position = "top")

    # ---- domain-level test + optional windows ----
    test_rows <- list()
    # whole-profile: center of mass per sample x domain, t-test between conditions
    com <- shares[, .(com = if (sum(share) > 0) sum(fraction * share) / sum(share) else NA_real_), by = .(sample, condition, domain)]
    for (dd in levels(shares$domain)) {
      cc <- com[domain == dd & is.finite(com)]
      if (length(unique(cc$condition)) == 2 && nrow(cc) >= 4) {
        tt <- tryCatch(stats::t.test(com ~ condition, data = cc), error = function(e) NULL)
        if (!is.null(tt)) test_rows[[paste0(dd, "_profile")]] <-
          data.table(domain = dd, window = "profile(center-of-mass)", method = "t-test",
                     p = tt$p.value, effect = diff(rev(tapply(cc$com, cc$condition, mean))))
      }
    }
    # windows: fractional share summed over each window, beta-regression / logit-t per domain x window
    win_long <- NULL
    if (!is.null(windows) && length(windows)) {
      win_long <- rbindlist(lapply(names(windows), function(wn) {
        fr <- as.numeric(windows[[wn]])
        shares[fraction %in% fr, .(share = sum(share)), by = .(sample, condition, domain)][, window := wn]
      }), use.names = TRUE)
      win_long[, window := factor(window, levels = names(windows))]
      for (dd in levels(shares$domain)) for (wn in names(windows)) {
        sw <- win_long[domain == dd & window == wn]
        if (nrow(sw) >= 4 && length(unique(sw$condition)) == 2) {
          tr <- .share_test(sw$share, sw$condition)
          mn <- tapply(sw$share, sw$condition, mean)
          test_rows[[paste0(dd, "_", wn)]] <- data.table(domain = dd, window = wn, method = tr$method,
                                                          p = tr$p, effect = as.numeric(diff(rev(mn))))
        }
      }
    }
    tests <- rbindlist(test_rows, use.names = TRUE)
    if (nrow(tests)) { tests[, p_BHadj := p.adjust(p, method = "BH")]; fwrite(tests[order(p)], file.path(tab_dir, paste0(pid, "_", m, "_domain_test.txt")), sep = "\t") }

    # ---- plot 2: fractional-share barplot over windows (mean +/- SD + replicate points) ----
    g2 <- NULL
    if (!is.null(win_long) && nrow(win_long)) {
      aggw <- win_long[, .(mean = mean(share), sd = stats::sd(share)), by = .(domain, window, condition)]
      aggw[is.na(sd), sd := 0]
      dodge <- position_dodge(width = 0.8)
      g2 <- ggplot(aggw, aes(window, mean, fill = condition)) +
        geom_col(position = dodge, width = 0.7, colour = "grey30", linewidth = 0.2) +
        geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), position = dodge, width = 0.2) +
        geom_point(data = win_long, aes(window, share, group = condition), position = position_dodge(width = 0.8),
                   shape = 21, colour = "grey20", fill = "white", size = 1.6) +
        facet_wrap(~ domain, scales = "free_y") +
        labs(title = paste0("Fractional share over windows: ", pid, " (", m, ")"),
             subtitle = "load-normalised share of total protein signal; points = replicates",
             x = NULL, y = "share of total protein signal", fill = NULL) +
        theme_bw() + theme(legend.position = "top")
    }

    if (save_pdf) {
      f1 <- file.path(fig_dir, paste0(pid, "_", m, "_chromatogram.pdf"))
      ggsave(f1, g1, width = 3 + 2.2 * length(unique(shares$domain)), height = 5, limitsize = FALSE); message("Wrote ", f1)
      if (!is.null(g2)) { f2 <- file.path(fig_dir, paste0(pid, "_", m, "_shares.pdf"))
        ggsave(f2, g2, width = 3 + 2.2 * length(unique(shares$domain)), height = 5, limitsize = FALSE); message("Wrote ", f2) }
    }
    if (print_plot) { print(g1); if (!is.null(g2)) print(g2) }
    if (nrow(tests)) { message("[", m, "] ", pid, " domain test (", nrow(tests), " row(s)) -> ", file.path(tab_dir, paste0(pid, "_", m, "_domain_test.txt"))); print(tests[order(p)]) }
    last <- list(shares = shares, tests = tests, chromatogram = g1, barplot = g2)
  }
  invisible(last)
}
