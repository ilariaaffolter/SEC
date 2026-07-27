# scripts/globularity_category_annotation.R
# =============================================================================
# CHARACTERISE the elution categories from scripts/globularity_check.R. For each category it answers:
#   (a) which GO terms are over-represented (hypergeometric test, same method as the report);
#   (b) what % of its proteins are subunits of a CURATED COMPLEX (Complex Portal) or have any annotated
#       binary interaction partner (UniProt cc_interaction);
#   (c) what % are annotated as FILAMENT / polymer-forming;
#   (d) what % carry intrinsically disordered regions (IDR);
#   (e) their ISOELECTRIC POINT / net charge - the column-interaction (artefact) test.
# Each percentage is also tested against the other tested proteins (Fisher), so "30% have IDRs" comes
# with "...which is / is not more than the globular set".
#
# CATEGORIES (from tables/globularity_check.txt):
#   beyond_calibration  apex above the largest calibration standard - the extrapolated region
#   sub_monomer         elutes clearly SMALLER than its own monomer
#   globular_1x_4x      the "as expected" set: monomer + clean oligomers 2x-4x
# (`categories =` accepts any class present in that table, plus the composite "globular_1x_4x".)
#
# ANNOTATION SOURCES - what is curated and what is a heuristic:
#   COMPLEX MEMBERSHIP  data/raw/<complex_portal_file>.tsv, the same Complex Portal export the report's
#                       complex analysis uses. CURATED. Accessions are parsed from the
#                       "identifiers and stoichiometry of molecules in complex" column.
#   INTERACTORS         UniProt cc_interaction (binary interactions). CURATED, but sparse for many
#                       organisms - a low % means "not annotated", not "no interactions".
#   FILAMENT            *** KEYWORD HEURISTIC, NOT CURATED ***: a regex over the GO terms and protein
#                       name (filament, polymeris(z)ation, cytoskelet*, flagell*, pilus/fimbria, ...).
#                       It will miss unannotated polymers and can over-call (e.g. "flagellar motor" is
#                       matched by 'flagell'). Treat as a screen; the matched terms are written out so
#                       every call can be checked. Tune with `filament_regex =`.
#   IDR                 default: FOLDINDEX (Prilusky et al. 2005, Bioinformatics 21:3435) computed from
#                       the UniProt sequence already in the cache - no external tool, no re-fetch:
#                           FoldIndex = 2.785 * <H> - |<R>| - 1.151
#                       with <H> mean Kyte-Doolittle hydropathy scaled to [0,1] and <R> mean net charge.
#                       Negative = predicted unfolded. Both a whole-protein value and a sliding-window
#                       "% disordered residues" are reported. It is a COARSE predictor - fine for
#                       comparing categories, not for calling any single protein.
#                       If output/hydropro/plddt_disorder.csv exists (from hydropro_ffo.R), the
#                       AlphaFold pLDDT-based disorder fraction is reported ALONGSIDE it, which is the
#                       better measure where available.
#   pI / NET CHARGE     computed from the cached sequence by Henderson-Hasselbalch titration with the
#                       EMBOSS 'iep' pKa set (pI by bisection; net charge evaluated at `buffer_pH`,
#                       default 7.4 - SET THIS TO YOUR ACTUAL RUNNING BUFFER). Predicted, not measured,
#                       but accurate enough to compare groups.
#                       WHY IT MATTERS: a sub_monomer protein elutes LATE, i.e. it is RETAINED. SEC
#                       resins carry residual negative charge, so cation-exchange-like retention acts on
#                       BASIC (net-positive) proteins. Two readouts are produced:
#                         * pI per category (violin + Wilcoxon vs the other tested proteins, plus % with
#                           pI > pI_basic_cut and % net-positive at buffer_pH, Fisher-tested);
#                         * pI vs log2(apparent/expected MW) across ALL tested proteins, with a Spearman
#                           correlation. A negative trend means basic proteins systematically elute late
#                           => charge-driven retention, i.e. a chromatography ARTEFACT rather than biology.
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "globularity_category_annotation.R"))
#   globularity_category_annotation("ATP")
#   globularity_category_annotation()                          # every metabolite that has been run
#   globularity_category_annotation("ATP", go_columns = c("go_p","go_f","go_c"))
#   globularity_category_annotation("ATP", categories = c("beyond_calibration","sub_monomer","globular_1x_4x","above_range"))
#
# OUTPUT (per metabolite, tables/globularity_categories/ and figures/globularity_categories/):
#   category_annotation_summary.txt    one row per category: n, % complex / interactors / filament / IDR
#                                      + Fisher p vs the other tested proteins
#   category_protein_annotation.txt    per protein: category + every annotation flag (auditable)
#   GOenrichment_<go_col>_<category>.txt   enriched GO terms per category
#   category_annotation_barplots.pdf   the percentages side by side across categories
#   category_pI_charge.pdf             pI per category (violin) + pI vs elution deviation (artefact test)
#   GOenrichment_<category>.pdf        GO bar plots
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })

# ---- GO over-representation: identical method to the report (hypergeometric on the ';'-split terms) --
.go_enrichment <- function(foreground, background, annotation, id_col = "accession",
                           go_col = "go_p", min_genes = 2, top_n = 15) {
  fg <- unique(stats::na.omit(unlist(foreground)))
  bg <- unique(stats::na.omit(unlist(c(background, fg))))
  if (length(fg) < min_genes || length(bg) < 5) return(NULL)
  am <- as.data.table(as.data.frame(annotation)[, c(id_col, go_col)])
  setnames(am, c(id_col, go_col), c("pid", "goann"))
  am <- unique(am[pid %in% bg & !is.na(goann) & nzchar(goann)])
  if (nrow(am) == 0) return(NULL)
  pt <- am[, .(term = trimws(unlist(strsplit(goann, ";")))), by = pid][nzchar(term)]
  N  <- uniqueN(pt$pid); n <- uniqueN(pt$pid[pt$pid %in% fg])
  if (n < min_genes) return(NULL)
  ts <- pt[, .(K = uniqueN(pid), k = uniqueN(pid[pid %in% fg])), by = term][k >= min_genes]
  if (nrow(ts) == 0) return(NULL)
  ts[, pval := stats::phyper(k - 1, K, N - K, n, lower.tail = FALSE)]
  ts[, padj := stats::p.adjust(pval, method = "BH")]
  ts[order(padj, -k)][seq_len(min(top_n, .N))]
}
.go_barplot <- function(go_tbl, title) {
  if (is.null(go_tbl) || !nrow(go_tbl)) return(NULL)
  gt <- copy(go_tbl); gt[, term_short := ifelse(nchar(term) > 60, paste0(substr(term, 1, 57), "..."), term)]
  ggplot(gt, aes(stats::reorder(term_short, -padj), -log10(padj), fill = k)) +
    geom_col() + coord_flip() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
    labs(title = title, x = NULL, y = "-log10(BH-adjusted p)", fill = "n proteins") + theme_bw()
}

# ---- FoldIndex disorder from sequence --------------------------------------------------------------
.KD <- c(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
         L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
# FoldIndex = 2.785*<H> - |<R>| - 1.151 ; <H> = Kyte-Doolittle scaled to [0,1]; <R> = net charge/residue.
.foldindex_vec <- function(aa) {
  h <- (unname(.KD[aa]) + 4.5) / 9                       # scaled hydropathy, NA for non-standard aa
  q <- ifelse(aa %in% c("R", "K"), 1, ifelse(aa %in% c("D", "E"), -1, 0))
  list(h = h, q = q)
}
# ---- isoelectric point / net charge from sequence ---------------------------------------------------
# Why this matters here: a protein eluting LATER than its own monomer (the sub_monomer class) is being
# RETAINED by the column. SEC resins carry residual negative charge, so cation-exchange-like retention
# acts on BASIC (high-pI, net-positive) proteins. If sub_monomer is enriched for high pI, the class is a
# chromatography artefact rather than biology - which is exactly what this test is for.
# pKa values: EMBOSS 'iep' defaults (the set behind most pI calculators). Standard Henderson-Hasselbalch
# titration, pI found by bisection. Cysteine/tyrosine are included as weak acids, as in EMBOSS.
.PKA <- list(Nterm = 8.6, Cterm = 3.6, K = 10.8, R = 12.5, H = 6.5, D = 3.9, E = 4.1, C = 8.5, Y = 10.1)
.net_charge_counts <- function(cnt, pH) {
  pos <- 1 / (1 + 10^(pH - .PKA$Nterm)) +
         cnt[["K"]] / (1 + 10^(pH - .PKA$K)) +
         cnt[["R"]] / (1 + 10^(pH - .PKA$R)) +
         cnt[["H"]] / (1 + 10^(pH - .PKA$H))
  neg <- 1 / (1 + 10^(.PKA$Cterm - pH)) +
         cnt[["D"]] / (1 + 10^(.PKA$D - pH)) +
         cnt[["E"]] / (1 + 10^(.PKA$E - pH)) +
         cnt[["C"]] / (1 + 10^(.PKA$C - pH)) +
         cnt[["Y"]] / (1 + 10^(.PKA$Y - pH))
  pos - neg
}
.charge_props <- function(seq, buffer_pH = 7.4) {
  aa <- strsplit(toupper(gsub("[^A-Za-z]", "", seq)), "")[[1]]
  if (!length(aa)) return(list(pI = NA_real_, net_charge = NA_real_, charge_per_res = NA_real_))
  cnt <- as.list(setNames(rep(0, 7), c("K", "R", "H", "D", "E", "C", "Y")))
  tb  <- table(aa)
  for (r in names(cnt)) if (!is.na(tb[r])) cnt[[r]] <- as.numeric(tb[r])
  lo <- 0; hi <- 14
  for (i in 1:60) {                                   # bisection: ~1e-4 pH units
    mid <- (lo + hi) / 2
    if (.net_charge_counts(cnt, mid) > 0) lo <- mid else hi <- mid
  }
  q <- .net_charge_counts(cnt, buffer_pH)
  list(pI = (lo + hi) / 2, net_charge = q, charge_per_res = q / length(aa))
}

.foldindex <- function(seq, window = 51L) {
  aa <- strsplit(toupper(gsub("[^A-Za-z]", "", seq)), "")[[1]]
  if (!length(aa)) return(list(global = NA_real_, disorder_frac = NA_real_, length = 0L))
  v <- .foldindex_vec(aa); h <- v$h; q <- v$q
  ok <- is.finite(h); if (!any(ok)) return(list(global = NA_real_, disorder_frac = NA_real_, length = length(aa)))
  glob <- 2.785 * mean(h[ok]) - abs(mean(q[ok])) - 1.151
  # sliding window -> fraction of residues sitting in a predicted-unfolded window
  n <- length(aa)
  if (n < window) return(list(global = glob, disorder_frac = as.numeric(glob < 0), length = n))
  h0 <- ifelse(is.finite(h), h, mean(h[ok]))
  ch <- cumsum(h0); cq <- cumsum(q)
  i1 <- seq_len(n - window + 1L); i2 <- i1 + window - 1L
  mh <- (ch[i2] - c(0, ch)[i1]) / window
  mq <- (cq[i2] - c(0, cq)[i1]) / window
  fi <- 2.785 * mh - abs(mq) - 1.151
  list(global = glob, disorder_frac = mean(fi < 0), length = n)
}

# ---- annotation sources ----------------------------------------------------------------------------
.load_uniprot <- function() {
  sf <- here("output", "uniprot_annotation_shared.RData")
  if (!file.exists(sf)) return(NULL)
  e <- new.env(); load(sf, envir = e)
  u <- tryCatch(as.data.table(e$.uniprot_all), error = function(x) NULL)
  if (is.null(u) || !"input_id" %in% names(u)) return(NULL)
  u[, accession := as.character(input_id)][]
}
# accessions of every protein listed as a subunit of a curated complex
.complex_portal_members <- function(file = NULL) {
  f <- file
  if (is.null(f)) {
    cand <- list.files(here("data", "raw"), pattern = "^Complex_portal.*\\.tsv$", full.names = TRUE)
    if (!length(cand)) return(NULL)
    f <- cand[1]
  } else if (!file.exists(f)) { f2 <- here("data", "raw", f); if (!file.exists(f2)) return(NULL); f <- f2 }
  cp <- tryCatch(as.data.table(read.csv(f, sep = "\t", header = TRUE, check.names = FALSE)), error = function(e) NULL)
  if (is.null(cp) || !nrow(cp)) return(NULL)
  idc <- grep("identifier", names(cp), ignore.case = TRUE, value = TRUE)
  col <- if (length(idc)) grep("molecul|stoichiom", idc, ignore.case = TRUE, value = TRUE)[1] else NA_character_
  if (is.na(col)) col <- if (length(idc)) idc[1] else NA_character_
  if (is.na(col)) { message("NOTE: Complex Portal file has no recognisable subunit-identifier column."); return(NULL) }
  acc <- unlist(regmatches(cp[[col]], gregexpr("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", cp[[col]])))
  unique(toupper(acc))
}

globularity_category_annotation <- function(
    metabolites   = NULL,
    categories    = c("beyond_calibration", "sub_monomer", "globular_1x_4x"),
    go_columns    = c("go_p", "go_f", "go_c"),
    go_min_genes  = 2, go_top_n = 15,
    filament_regex = "filament|polymeriz|polymeris|cytoskelet|flagell|pilus|pili\\b|fimbri|microtubul|actin|tubulin|z-ring|divisome",
    foldindex_window = 51L,
    buffer_pH     = 7.4,   # running-buffer pH: net charge is evaluated here (set to your actual buffer)
    pI_basic_cut  = 8.0,   # "basic" protein cutoff for the categorical pI test
    complex_portal_file = NULL,
    out_subdir    = "globularity_categories") {

  U <- .load_uniprot()
  if (is.null(U)) stop("No output/uniprot_annotation_shared.RData - render at least one comparison first.")
  cp_members <- .complex_portal_members(complex_portal_file)
  if (is.null(cp_members)) message("NOTE: no Complex Portal .tsv found in data/raw - the complex % will be NA.")
  plddt_f <- here("output", "hydropro", "plddt_disorder.csv")
  PL <- if (file.exists(plddt_f)) fread(plddt_f) else NULL

  # FoldIndex per protein (once, over the whole cached proteome)
  if (!"sequence" %in% names(U)) stop("The UniProt cache has no `sequence` column - re-render so it is fetched.")
  message("Computing FoldIndex disorder + pI/net charge for ", nrow(U), " cached protein(s) ...")
  fi <- lapply(U$sequence, function(s) if (is.na(s) || !nzchar(s)) list(global = NA_real_, disorder_frac = NA_real_, length = NA_integer_) else .foldindex(s, foldindex_window))
  U[, `:=`(foldindex_global = vapply(fi, function(x) x$global, numeric(1)),
           foldindex_disorder_frac = vapply(fi, function(x) x$disorder_frac, numeric(1)))]
  ch <- lapply(U$sequence, function(s) if (is.na(s) || !nzchar(s)) list(pI = NA_real_, net_charge = NA_real_, charge_per_res = NA_real_) else .charge_props(s, buffer_pH))
  U[, `:=`(pI              = vapply(ch, function(x) x$pI, numeric(1)),
           net_charge      = vapply(ch, function(x) x$net_charge, numeric(1)),
           charge_per_res  = vapply(ch, function(x) x$charge_per_res, numeric(1)))]
  # filament keyword hit over GO terms + protein name (recorded, so every call is auditable)
  .txtcols <- intersect(c("go_p", "go_f", "go_c", "protein_name"), names(U))
  .stxt <- Reduce(function(a, b) paste(a, b, sep = " ; "),
                  lapply(.txtcols, function(cc) { v <- as.character(U[[cc]]); v[is.na(v)] <- ""; v }))
  U[, search_text := .stxt]
  U[, filament_hit := grepl(filament_regex, search_text, ignore.case = TRUE)]
  U[, filament_terms := vapply(search_text, function(s) {
    tt <- trimws(unlist(strsplit(s, ";")))
    paste(unique(tt[grepl(filament_regex, tt, ignore.case = TRUE)]), collapse = " | ") },
    character(1), USE.NAMES = FALSE)]
  U[, has_interaction := if ("cc_interaction" %in% names(U)) !is.na(cc_interaction) & nzchar(cc_interaction) else NA]
  U[, in_complex_portal := if (is.null(cp_members)) NA else toupper(accession) %in% cp_members]

  if (is.null(metabolites)) {
    dirs        <- basename(list.dirs(here("output"), recursive = FALSE))
    metabolites <- sub("^PCM_ctrl_vs_", "", dirs[grepl("^PCM_ctrl_vs_", dirs)])
  }
  if (!length(metabolites)) stop("No PCM_ctrl_vs_* output folders found.")

  for (m in metabolites) {
    gf <- here("output", paste0("PCM_ctrl_vs_", m), "tables", "globularity_check.txt")
    if (!file.exists(gf)) { message("[", m, "] no globularity_check.txt - run globularity_check() first; skipping."); next }
    G <- fread(gf)
    if (!all(c("protein_id", "class") %in% names(G))) { message("[", m, "] globularity_check.txt lacks protein_id/class; skipping."); next }

    # category membership (globular_1x_4x = monomer + any clean oligomer)
    sets <- list()
    for (cat in categories) {
      ids <- if (cat == "globular_1x_4x") G[class == "monomer" | grepl("^oligomer_", class)]$protein_id else G[class == cat]$protein_id
      sets[[cat]] <- unique(as.character(ids))
    }
    sets <- sets[vapply(sets, length, 1L) > 0]
    if (!length(sets)) { message("[", m, "] none of the requested categories has proteins; skipping."); next }
    universe <- unique(as.character(G$protein_id))          # background = everything tested here

    tab_dir <- here("output", paste0("PCM_ctrl_vs_", m), "tables", out_subdir)
    fig_dir <- here("output", paste0("PCM_ctrl_vs_", m), "figures", out_subdir)
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

    # ---- per-protein annotation table (auditable) ----
    ann_cols <- intersect(c("accession", "protein_name", "gene_names", "in_complex_portal", "has_interaction",
                            "filament_hit", "filament_terms", "foldindex_global", "foldindex_disorder_frac",
                            "pI", "net_charge", "charge_per_res"), names(U))
    A <- merge(G[, .(protein_id, class, expected_mw_kDa, apparent_mw_kDa, ratio, ffo_vs_monomer)],
               U[, ..ann_cols], by.x = "protein_id", by.y = "accession", all.x = TRUE)
    if (!is.null(PL) && "plddt_disorder_frac" %in% names(PL))
      A <- merge(A, PL[, .(protein_id, mean_plddt, plddt_disorder_frac)], by = "protein_id", all.x = TRUE)
    # start from the raw class, then overwrite with the requested categories (fcase cannot take a
    # vector `default`, and a protein may belong to a composite category such as globular_1x_4x)
    A[, category := as.character(class)]
    for (cat in names(sets)) A[protein_id %in% sets[[cat]], category := cat]
    fwrite(A, file.path(tab_dir, "category_protein_annotation.txt"), sep = "\t")

    # ---- percentages per category, each Fisher-tested against the OTHER tested proteins ----
    .pct <- function(ids, flag) { v <- A[protein_id %in% ids][[flag]]; v <- v[!is.na(v)]
                                  if (!length(v)) return(c(n = 0, pct = NA_real_)); c(n = sum(v), pct = 100 * mean(v)) }
    .fisher <- function(ids, flag) {
      inn <- A[protein_id %in% ids][[flag]]; out <- A[!(protein_id %in% ids)][[flag]]
      inn <- inn[!is.na(inn)]; out <- out[!is.na(out)]
      if (!length(inn) || !length(out)) return(NA_real_)
      tryCatch(stats::fisher.test(matrix(c(sum(inn), length(inn) - sum(inn),
                                           sum(out), length(out) - sum(out)), nrow = 2))$p.value,
               error = function(e) NA_real_)
    }
    A[, idr_foldindex := foldindex_global < 0]                       # predicted unfolded overall
    A[, idr_highfrac  := foldindex_disorder_frac > 0.3]              # >30% of residues in unfolded windows
    if ("plddt_disorder_frac" %in% names(A)) A[, idr_plddt := plddt_disorder_frac > 0.4]
    A[, basic_pI      := pI > pI_basic_cut]                          # net-positive at neutral pH
    A[, net_positive  := net_charge > 0]                             # at the running-buffer pH

    flags <- intersect(c("in_complex_portal", "has_interaction", "filament_hit",
                         "idr_foldindex", "idr_highfrac", "idr_plddt",
                         "basic_pI", "net_positive"), names(A))
    # continuous variables compared per category with a two-sided Wilcoxon test (vs all other tested proteins)
    conts <- intersect(c("pI", "net_charge", "charge_per_res", "foldindex_disorder_frac"), names(A))
    .wilcox <- function(ids, v) {
      inn <- A[protein_id %in% ids][[v]]; out <- A[!(protein_id %in% ids)][[v]]
      inn <- inn[is.finite(inn)]; out <- out[is.finite(out)]
      if (length(inn) < 3 || length(out) < 3) return(NA_real_)
      tryCatch(stats::wilcox.test(inn, out)$p.value, error = function(e) NA_real_)
    }
    rows <- lapply(names(sets), function(cat) {
      ids <- sets[[cat]]
      r <- data.table(metabolite = m, category = cat, n_proteins = length(ids))
      for (fl in flags) {
        p <- .pct(ids, fl)
        r[[paste0("pct_", fl)]]      <- round(unname(p["pct"]), 1)
        r[[paste0("n_", fl)]]        <- unname(p["n"])
        r[[paste0("fisher_p_", fl)]] <- signif(.fisher(ids, fl), 3)
      }
      for (v in conts) {
        r[[paste0("median_", v)]]    <- round(stats::median(A[protein_id %in% ids][[v]], na.rm = TRUE), 3)
        r[[paste0("median_rest_", v)]] <- round(stats::median(A[!(protein_id %in% ids)][[v]], na.rm = TRUE), 3)
        r[[paste0("wilcox_p_", v)]] <- signif(.wilcox(ids, v), 3)
      }
      r
    })
    S <- rbindlist(rows, use.names = TRUE, fill = TRUE)
    fwrite(S, file.path(tab_dir, "category_annotation_summary.txt"), sep = "\t")
    message("\n[", m, "] category annotation (Fisher p vs all other tested proteins):")
    print(S[, c("category", "n_proteins", grep("^pct_", names(S), value = TRUE)), with = FALSE])
    message("[", m, "] pI / charge per category (Wilcoxon vs all other tested proteins):")
    print(S[, c("category", "n_proteins",
                intersect(c("median_pI", "median_rest_pI", "wilcox_p_pI",
                            "median_net_charge", "wilcox_p_net_charge"), names(S))), with = FALSE])

    # ---- bar plot of the percentages across categories ----
    L <- melt(S, id.vars = c("category", "n_proteins"),
              measure.vars = grep("^pct_", names(S), value = TRUE),
              variable.name = "metric", value.name = "pct")
    L[, metric := sub("^pct_", "", metric)]
    L[, category := factor(category, levels = names(sets))]
    gbar <- ggplot(L[is.finite(pct)], aes(category, pct, fill = category)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = sprintf("%.0f%%", pct)), vjust = -0.3, size = 3) +
      facet_wrap(~ metric, scales = "free_y") +
      labs(title = paste0("Elution-category annotation - PCM_ctrl_vs_", m),
           subtitle = paste0("in_complex_portal / has_interaction = curated | filament_hit = KEYWORD HEURISTIC | ",
                             "idr_* = FoldIndex prediction", if ("idr_plddt" %in% flags) " and AlphaFold pLDDT" else ""),
           x = NULL, y = "% of category") +
      theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
    .fb <- file.path(fig_dir, "category_annotation_barplots.pdf")
    tryCatch(ggsave(.fb, gbar, width = 9, height = 6),
             error = function(e) message("   !! could not write ", basename(.fb), ": ", conditionMessage(e)))

    # ---- pI / charge: the column-interaction (artefact) test ----
    # (1) pI distribution per category; (2) pI vs how far the protein elutes from its expected position,
    # across ALL tested proteins. A systematic trend in (2) is the signature of charge-driven retention
    # on the column rather than biology - basic proteins retained (eluting late) by a negatively charged
    # resin appear as sub_monomer.
    AP <- A[is.finite(pI)]
    if (nrow(AP) > 10) {
      AP[, category := factor(category, levels = unique(c(names(sets), setdiff(unique(category), names(sets)))))]
      g_pi <- ggplot(AP, aes(category, pI, fill = category)) +
        geom_violin(alpha = 0.45, colour = NA, scale = "width") +
        geom_boxplot(width = 0.16, outlier.size = 0.4, fill = "white") +
        geom_hline(yintercept = buffer_pH, linetype = 2, colour = "grey40") +
        annotate("text", x = 0.6, y = buffer_pH, label = paste0("buffer pH ", buffer_pH),
                 hjust = 0, vjust = -0.5, size = 3, colour = "grey35") +
        labs(title = paste0("Isoelectric point by elution category - PCM_ctrl_vs_", m),
             subtitle = paste0("Proteins above the dashed line are net POSITIVE in the running buffer.\n",
                               "sub_monomer enriched for high pI => cation-exchange-like retention on the column (artefact), not biology."),
             x = NULL, y = "predicted pI (EMBOSS pKa set)") +
        theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

      AP2 <- AP[is.finite(ratio) & ratio > 0]
      rho <- if (nrow(AP2) > 10) suppressWarnings(stats::cor(AP2$pI, log2(AP2$ratio), method = "spearman", use = "complete.obs")) else NA_real_
      rho_p <- if (nrow(AP2) > 10) tryCatch(suppressWarnings(stats::cor.test(AP2$pI, log2(AP2$ratio), method = "spearman"))$p.value, error = function(e) NA_real_) else NA_real_
      g_pi2 <- ggplot(AP2, aes(pI, log2(ratio))) +
        geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
        geom_vline(xintercept = buffer_pH, linetype = 3, colour = "grey55") +
        geom_point(aes(colour = category), alpha = 0.45, size = 0.9) +
        geom_smooth(method = "loess", se = TRUE, colour = "black", linewidth = 0.6, formula = y ~ x) +
        labs(title = paste0("pI vs elution deviation - PCM_ctrl_vs_", m),
             subtitle = sprintf("y = log2(apparent / expected MW); below 0 = elutes late (retained). Spearman rho = %.3f (p = %.3g).\nA negative trend = basic proteins retained by the resin: a charge artefact, not biology.",
                                rho, rho_p),
             x = "predicted pI", y = "log2(apparent / expected MW)", colour = NULL) +
        theme_bw() + theme(legend.position = "bottom")

      .fpi <- file.path(fig_dir, "category_pI_charge.pdf")
      tryCatch({ grDevices::pdf(.fpi, width = 8, height = 6); print(g_pi); print(g_pi2); grDevices::dev.off() },
               error = function(e) { message("   !! could not write ", basename(.fpi), ": ", conditionMessage(e))
                                     try(grDevices::dev.off(), silent = TRUE) })
      message(sprintf("[%s] pI vs log2(apparent/expected): Spearman rho = %.3f (p = %.3g) over %d protein(s).",
                      m, rho, rho_p, nrow(AP2)))
    }

    # ---- GO enrichment per category x GO namespace ----
    for (cat in names(sets)) {
      plots <- list()
      for (gc in intersect(go_columns, names(U))) {
        gt <- tryCatch(.go_enrichment(sets[[cat]], universe, U, "accession", gc, go_min_genes, go_top_n),
                       error = function(e) { message("[", m, "/", cat, "/", gc, "] GO failed: ", conditionMessage(e)); NULL })
        if (is.null(gt) || !nrow(gt)) { message("[", m, "] ", cat, " / ", gc, ": no enriched terms."); next }
        fwrite(gt, file.path(tab_dir, paste0("GOenrichment_", gc, "_", cat, ".txt")), sep = "\t")
        message("[", m, "] ", cat, " / ", gc, ": ", sum(gt$padj < 0.05), " term(s) at BH<0.05 (top: ",
                substr(gt$term[1], 1, 60), ")")
        plots[[gc]] <- .go_barplot(gt, paste0("GO ", gc, " - ", cat, " (", length(sets[[cat]]), " proteins)"))
      }
      if (length(plots)) {
        .fg <- file.path(fig_dir, paste0("GOenrichment_", cat, ".pdf"))
        tryCatch({ grDevices::pdf(.fg, width = 8, height = 5.5); for (p in plots) print(p); grDevices::dev.off() },
                 error = function(e) { message("   !! could not write ", basename(.fg), ": ", conditionMessage(e))
                                       try(grDevices::dev.off(), silent = TRUE) })
      }
    }
    message("[", m, "] -> ", tab_dir)
  }
  invisible(NULL)
}
