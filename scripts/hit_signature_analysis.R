# scripts/hit_signature_analysis.R
# =============================================================================
# SPECIFICITY TEST for the differential hits: do they look like SPECIFIC LIGAND BINDERS, or like a
# HYDROTROPE / non-specific solubility effect?
#
# THE QUESTION: a metabolite such as ATP is both a ligand and a hydrotrope. If the hits are enriched for
# nucleotide-binding folds and ligand-binding annotation, the effect is specific binding. If they are
# instead enriched for intrinsic disorder and net positive charge - the signature of proteins whose
# solubility/compaction responds to a charged small molecule - the effect is physicochemical rather than
# a binding event. This script tests both signatures side by side, hits versus the tested background.
#
# SIGNATURES COMPARED (hits vs all proteins tested in that comparison):
#   SPECIFIC        ligand_binding_hit   metabolite-relevant annotation (GO molecular function,
#                                        UniProt ft_binding / ft_act_site / cc_catalytic_activity and the
#                                        protein name) matched by a per-metabolite regex - see
#                                        `ligand_regex` below, which is a HEURISTIC you can and should
#                                        edit. Every matched term is written out so each call is auditable.
#                   nucleotide_fold_hit  generic nucleotide-binding evidence (P-loop, ATP/GTP/NAD binding).
#   HYDROTROPE      idr_foldindex        predicted intrinsically unfolded overall (FoldIndex < 0)
#                   idr_highfrac         >30% of residues in predicted-unfolded windows
#                   net_positive         net charge > 0 at the running-buffer pH
#                   basic_pI             pI above pI_basic_cut
#   plus the continuous versions (FoldIndex disorder fraction, pI, net charge, charge per residue).
#
# TESTS: Fisher's exact test for each binary flag, Wilcoxon for each continuous variable, hits vs the
# rest of the tested proteins; odds ratios reported for the binary flags. A GO molecular-function
# over-representation of the hits is run as well (same hypergeometric method as the report).
#
# READING THE RESULT: the two signatures are NOT mutually exclusive - report what is enriched rather than
# forcing a verdict. Enrichment of ligand/nucleotide annotation WITHOUT disorder/charge enrichment is the
# clean "specific binding" outcome; the reverse pattern is the hydrotrope outcome; both enriched means the
# hit set is mixed and should be split before interpretation.
#
# HITS are defined exactly as in scripts/overlap_between_metabolites.R:
#   protein_DiffExprProtein: pBHadj < 0.05 AND |medianLog2FC| > 1
# read from output/<cmp>/rdata/protein_DiffExprProtein_list.RData (fall back: the *_for_plotting.RData).
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "hit_signature_analysis.R"))
#   hit_signature_analysis()                       # every metabolite that has been run
#   hit_signature_analysis("ATP")
#   hit_signature_analysis("ATP", buffer_pH = 7.4, pI_basic_cut = 8)
#   hit_signature_analysis("ATP", ligand_regex = list(ATP = "ATP|adenyl|P-loop|kinase"))   # override
#
# OUTPUT (per metabolite, tables/ and figures/hit_signature/):
#   hit_signature_summary.txt      one row per metric: % in hits, % in background, OR, p, direction
#   hit_protein_annotation.txt     per hit protein: every flag + the matched ligand terms (auditable)
#   GOenrichment_go_f_hits.txt     molecular-function over-representation of the hits
#   hit_signature.pdf              binary flags (hits vs background) + continuous distributions
# OUTPUT (pooled): output/hit_signature/hit_signature_summary_all.csv + signature_overview.pdf
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })

# reuse the annotation helpers (.load_uniprot, .charge_props, .foldindex, .go_enrichment, .go_barplot);
# that file only defines functions, so sourcing it has no side effects
source(here::here("scripts", "globularity_category_annotation.R"))

# Per-metabolite ligand annotation patterns. HEURISTIC - edit or override via `ligand_regex`.
# The generic nucleotide pattern is applied to every metabolite as the "nucleotide fold" signature.
.default_ligand_regex <- function() list(
  ATP = "ATP|adenosine.?tri|adenyl|P-loop|kinase|ATPase|adenine nucleotide",
  ADP = "ADP|ATP|adenosine.?di|adenyl|P-loop|kinase|ATPase|adenine nucleotide",
  NAD = "NAD|nicotinamide|dinucleotide|dehydrogenase|oxidoreductase|redox",
  aKG = "2-oxoglutarate|oxoglutarate|alpha-ketoglutarate|ketoglutarate|tricarboxylic|citrate cycle|TCA",
  PEP = "phosphoenolpyruvate|pyruvate|glycoly|phosphotransferase",
  PGP = "phosphoglycer|glycoly|phosphoglycolate|bisphosphoglycerate|mutase")
.NUCLEOTIDE_REGEX <- "P-loop|nucleotide.?binding|ATP.?binding|GTP.?binding|NAD.?binding|nucleoside.?triphosphate"

hit_signature_analysis <- function(metabolites   = NULL,
                                   ligand_regex   = NULL,
                                   buffer_pH      = 7.4,
                                   pI_basic_cut   = 8.0,
                                   foldindex_window = 51L,
                                   pBHadj_cut     = 0.05,
                                   log2fc_cut     = 1,
                                   go_min_genes   = 2, go_top_n = 15,
                                   out_subdir     = "hit_signature") {

  U <- .load_uniprot()
  if (is.null(U)) stop("No output/uniprot_annotation_shared.RData - render at least one comparison first.")
  if (!"sequence" %in% names(U)) stop("The UniProt cache has no `sequence` column - re-render so it is fetched.")
  LR <- utils::modifyList(.default_ligand_regex(), if (is.null(ligand_regex)) list() else as.list(ligand_regex))

  # sequence-derived properties, once
  message("Computing FoldIndex disorder + pI/net charge for ", nrow(U), " cached protein(s) ...")
  fi <- lapply(U$sequence, function(s) if (is.na(s) || !nzchar(s)) list(global = NA_real_, disorder_frac = NA_real_) else .foldindex(s, foldindex_window))
  ch <- lapply(U$sequence, function(s) if (is.na(s) || !nzchar(s)) list(pI = NA_real_, net_charge = NA_real_, charge_per_res = NA_real_) else .charge_props(s, buffer_pH))
  U[, `:=`(foldindex_global        = vapply(fi, function(x) x$global, numeric(1)),
           foldindex_disorder_frac = vapply(fi, function(x) x$disorder_frac, numeric(1)),
           pI                      = vapply(ch, function(x) x$pI, numeric(1)),
           net_charge              = vapply(ch, function(x) x$net_charge, numeric(1)),
           charge_per_res          = vapply(ch, function(x) x$charge_per_res, numeric(1)))]
  # searchable annotation text: molecular function + binding/active sites + catalytic activity + name
  .txtcols <- intersect(c("go_f", "ft_binding", "ft_act_site", "cc_catalytic_activity", "cc_cofactor",
                          "protein_name", "go_p"), names(U))
  U[, annot_text := Reduce(function(a, b) paste(a, b, sep = " ; "),
                           lapply(.txtcols, function(cc) { v <- as.character(U[[cc]]); v[is.na(v)] <- ""; v }))]
  U[, nucleotide_fold_hit := grepl(.NUCLEOTIDE_REGEX, annot_text, ignore.case = TRUE)]

  if (is.null(metabolites)) {
    dirs        <- basename(list.dirs(here("output"), recursive = FALSE))
    metabolites <- sub("^PCM_ctrl_vs_", "", dirs[grepl("^PCM_ctrl_vs_", dirs)])
  }
  if (!length(metabolites)) stop("No PCM_ctrl_vs_* output folders found.")

  all_rows <- list()
  for (m in metabolites) {
    cmp <- paste0("PCM_ctrl_vs_", m)
    # ---- hits + background, defined exactly as in overlap_between_metabolites.R ----
    pdp <- NULL
    f1 <- here("output", cmp, "rdata", "protein_DiffExprProtein_list.RData")
    if (file.exists(f1)) { e <- new.env(); load(f1, envir = e)
      if ("protein_DiffExprProtein" %in% ls(e)) pdp <- as.data.table(e$protein_DiffExprProtein) }
    if (is.null(pdp)) {
      fdir <- here("output", cmp, "RData_for_further_plotting_and_analysis")
      f2 <- list.files(fdir, pattern = "_for_plotting\\.RData$", full.names = TRUE)
      if (length(f2)) { e <- new.env(); load(f2[1], envir = e)
        if ("protein_DiffExprProtein" %in% ls(e)) pdp <- as.data.table(e$protein_DiffExprProtein) }
    }
    if (is.null(pdp) || !all(c("feature_id", "pBHadj", "medianLog2FC") %in% names(pdp))) {
      message("[", m, "] no protein_DiffExprProtein with pBHadj/medianLog2FC - render this comparison first; skipping."); next
    }
    background <- unique(as.character(pdp$feature_id))
    hits       <- unique(as.character(pdp[pBHadj < pBHadj_cut & abs(medianLog2FC) > log2fc_cut]$feature_id))
    if (length(hits) < 5) { message("[", m, "] only ", length(hits), " hit(s) - too few to test; skipping."); next }

    # ---- annotate the tested proteins ----
    A <- U[accession %in% background]
    if (!nrow(A)) { message("[", m, "] none of the tested proteins is in the UniProt cache; skipping."); next }
    rgx <- LR[[m]]
    if (is.null(rgx)) { rgx <- .NUCLEOTIDE_REGEX
      message("[", m, "] no ligand pattern defined for this metabolite - using the generic nucleotide pattern. ",
              "Pass ligand_regex = list(", m, " = '...') to set one.") }
    A[, ligand_binding_hit := grepl(rgx, annot_text, ignore.case = TRUE)]
    A[, ligand_terms := vapply(annot_text, function(s) {
      tt <- trimws(unlist(strsplit(s, ";")))
      paste(unique(tt[grepl(rgx, tt, ignore.case = TRUE)]), collapse = " | ") }, character(1), USE.NAMES = FALSE)]
    A[, `:=`(idr_foldindex = foldindex_global < 0,
             idr_highfrac  = foldindex_disorder_frac > 0.3,
             net_positive  = net_charge > 0,
             basic_pI      = pI > pI_basic_cut,
             is_hit        = accession %in% hits)]

    tab_dir <- here("output", cmp, "tables", out_subdir)
    fig_dir <- here("output", cmp, "figures", out_subdir)
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    keep <- intersect(c("accession", "protein_name", "gene_names", "ligand_binding_hit", "ligand_terms",
                        "nucleotide_fold_hit", "idr_foldindex", "idr_highfrac", "foldindex_disorder_frac",
                        "pI", "net_charge", "charge_per_res", "net_positive", "basic_pI"), names(A))
    fwrite(A[is_hit == TRUE, ..keep], file.path(tab_dir, "hit_protein_annotation.txt"), sep = "\t")

    # ---- Fisher (binary) and Wilcoxon (continuous), hits vs the rest of the tested proteins ----
    bin_flags <- c(specific = "ligand_binding_hit", specific = "nucleotide_fold_hit",
                   hydrotrope = "idr_foldindex", hydrotrope = "idr_highfrac",
                   hydrotrope = "net_positive", hydrotrope = "basic_pI")
    rows <- list()
    for (i in seq_along(bin_flags)) {
      fl <- bin_flags[[i]]; sig <- names(bin_flags)[i]
      if (!fl %in% names(A)) next
      inn <- A[is_hit == TRUE][[fl]];  inn <- inn[!is.na(inn)]
      out <- A[is_hit == FALSE][[fl]]; out <- out[!is.na(out)]
      if (!length(inn) || !length(out)) next
      tt <- tryCatch(stats::fisher.test(matrix(c(sum(inn), length(inn) - sum(inn),
                                                 sum(out), length(out) - sum(out)), nrow = 2)),
                     error = function(e) NULL)
      rows[[fl]] <- data.table(metabolite = m, signature = sig, metric = fl, type = "binary",
                               n_hits = length(inn), pct_hits = round(100 * mean(inn), 1),
                               pct_background = round(100 * mean(out), 1),
                               odds_ratio = if (is.null(tt)) NA_real_ else round(unname(tt$estimate), 3),
                               p = if (is.null(tt)) NA_real_ else signif(tt$p.value, 3))
    }
    cont <- c(hydrotrope = "foldindex_disorder_frac", hydrotrope = "pI",
              hydrotrope = "net_charge", hydrotrope = "charge_per_res")
    for (i in seq_along(cont)) {
      v <- cont[[i]]; sig <- names(cont)[i]
      if (!v %in% names(A)) next
      inn <- A[is_hit == TRUE][[v]];  inn <- inn[is.finite(inn)]
      out <- A[is_hit == FALSE][[v]]; out <- out[is.finite(out)]
      if (length(inn) < 3 || length(out) < 3) next
      tt <- tryCatch(stats::wilcox.test(inn, out), error = function(e) NULL)
      rows[[v]] <- data.table(metabolite = m, signature = sig, metric = v, type = "continuous",
                              n_hits = length(inn), median_hits = round(stats::median(inn), 3),
                              median_background = round(stats::median(out), 3),
                              p = if (is.null(tt)) NA_real_ else signif(tt$p.value, 3))
    }
    S <- rbindlist(rows, use.names = TRUE, fill = TRUE)
    S[, p_BHadj := signif(stats::p.adjust(p, method = "BH"), 3)]
    fwrite(S, file.path(tab_dir, "hit_signature_summary.txt"), sep = "\t")
    all_rows[[m]] <- S

    message("\n[", m, "] ", length(hits), " hit(s) vs ", length(background) - length(hits),
            " background protein(s) - specificity signatures:")
    print(S[order(signature, p)])
    .sp <- S[metric %in% c("ligand_binding_hit", "nucleotide_fold_hit") & p_BHadj < 0.05 & odds_ratio > 1]
    .hy <- S[signature == "hydrotrope" & p_BHadj < 0.05]
    message("[", m, "] => specific-binding signature enriched: ", if (nrow(.sp)) paste(.sp$metric, collapse = ", ") else "none",
            " | hydrotrope signature enriched: ", if (nrow(.hy)) paste(.hy$metric, collapse = ", ") else "none")

    # ---- GO molecular-function over-representation of the hits ----
    if ("go_f" %in% names(U)) {
      gt <- tryCatch(.go_enrichment(hits, background, U, "accession", "go_f", go_min_genes, go_top_n),
                     error = function(e) NULL)
      if (!is.null(gt) && nrow(gt)) {
        fwrite(gt, file.path(tab_dir, "GOenrichment_go_f_hits.txt"), sep = "\t")
        message("[", m, "] top GO-F term: ", substr(gt$term[1], 1, 70), " (BH p = ", signif(gt$padj[1], 3), ")")
      } else message("[", m, "] no enriched GO molecular-function terms.")
    }

    # ---- figure: binary flags side by side + continuous distributions ----
    Lb <- S[type == "binary"]
    p1 <- if (nrow(Lb)) {
      Lm <- melt(Lb, id.vars = c("metric", "signature"), measure.vars = c("pct_hits", "pct_background"),
                 variable.name = "set", value.name = "pct")
      Lm[, set := factor(ifelse(set == "pct_hits", "hits", "background"), levels = c("hits", "background"))]
      ggplot(Lm, aes(metric, pct, fill = set)) +
        geom_col(position = position_dodge(width = 0.78), width = 0.72) +
        geom_text(aes(label = sprintf("%.0f%%", pct)), position = position_dodge(width = 0.78), vjust = -0.3, size = 2.8) +
        facet_grid(~ signature, scales = "free_x", space = "free_x") +
        scale_fill_manual(values = c(hits = "#E15759", background = "#9AA5B1"), name = NULL) +
        labs(title = paste0("Hit signature - PCM_ctrl_vs_", m, " (", length(hits), " hits)"),
             subtitle = paste0("specific = ligand / nucleotide-fold annotation | hydrotrope = disorder and net positive charge.\n",
                               "Ligand pattern (HEURISTIC, editable): ", substr(rgx, 1, 90)),
             x = NULL, y = "% of set") +
        theme_bw() + theme(legend.position = "top", axis.text.x = element_text(angle = 25, hjust = 1))
    } else NULL
    Ac <- melt(A[, c("is_hit", intersect(c("foldindex_disorder_frac", "pI", "net_charge"), names(A))), with = FALSE],
               id.vars = "is_hit", variable.name = "metric", value.name = "value")
    Ac <- Ac[is.finite(value)]
    p2 <- if (nrow(Ac)) {
      Ac[, set := factor(ifelse(is_hit, "hits", "background"), levels = c("hits", "background"))]
      ggplot(Ac, aes(set, value, fill = set)) +
        geom_violin(alpha = 0.5, colour = NA, scale = "width") +
        geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white") +
        facet_wrap(~ metric, scales = "free_y") +
        scale_fill_manual(values = c(hits = "#E15759", background = "#9AA5B1")) +
        labs(title = paste0("Hydrotrope-signature variables - PCM_ctrl_vs_", m),
             subtitle = "Wilcoxon p values in hit_signature_summary.txt", x = NULL, y = NULL) +
        theme_bw() + theme(legend.position = "none")
    } else NULL
    .fp <- file.path(fig_dir, "hit_signature.pdf")
    tryCatch({ grDevices::pdf(.fp, width = 8, height = 5.5)
               if (!is.null(p1)) print(p1); if (!is.null(p2)) print(p2); grDevices::dev.off() },
             error = function(e) { message("   !! could not write ", basename(.fp), ": ", conditionMessage(e))
                                   try(grDevices::dev.off(), silent = TRUE) })
  }

  if (length(all_rows)) {
    SA <- rbindlist(all_rows, use.names = TRUE, fill = TRUE)
    out <- here("output", out_subdir); dir.create(out, recursive = TRUE, showWarnings = FALSE)
    fwrite(SA, file.path(out, "hit_signature_summary_all.csv"))
    cat("\n==== hit signatures across metabolites (binary flags) ====\n")
    print(SA[type == "binary", .(metabolite, signature, metric, pct_hits, pct_background, odds_ratio, p_BHadj)])
    g <- ggplot(SA[type == "binary" & is.finite(odds_ratio) & odds_ratio > 0],
                aes(metabolite, log2(odds_ratio), fill = signature)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.72) +
      geom_hline(yintercept = 0, colour = "grey30") +
      facet_wrap(~ metric) +
      labs(title = "Specific-binding vs hydrotrope signature of the differential hits",
           subtitle = "log2 odds ratio, hits vs tested background. Above 0 = enriched in the hits.",
           x = NULL, y = "log2(odds ratio)", fill = NULL) +
      theme_bw() + theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust = 1))
    .fo <- file.path(out, "signature_overview.pdf")
    tryCatch(ggsave(.fo, g, width = 9, height = 6),
             error = function(e) message("   !! could not write ", basename(.fo), ": ", conditionMessage(e)))
    invisible(SA)
  } else { message("Nothing computed."); invisible(NULL) }
}
