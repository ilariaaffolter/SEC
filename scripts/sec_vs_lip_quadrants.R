# scripts/sec_vs_lip_quadrants.R
# =============================================================================
# THE SEC x LiP QUADRANT ANALYSIS  -  crossing this study's SEC hits (changes in
# elution profile = a QUATERNARY / assembly readout) against Piazza et al. 2018
# LiP-SMap hits (limited proteolysis = a LOCAL / conformational binding readout),
# metabolite by metabolite.
#
# WHY THE TWO READOUTS ARE ORTHOGONAL (and why "mostly non-overlapping" is the
# expected result, not a failure):
#   * SEC / elution shift  -> the protein changed SIZE or assembly state when the
#                             metabolite was added (monomer<->oligomer, complex
#                             gain/loss). It is blind to a binding event that does
#                             not move the protein on the column.
#   * LiP (Piazza 2018)    -> the protein changed its LOCAL protease accessibility
#                             when the metabolite was added (a binding pocket
#                             closing, a loop ordering). It is blind to an assembly
#                             change that leaves local structure intact.
#   A metabolite can bind and rigidify a pocket (LiP+) without changing the
#   particle's size (SEC-), and it can trigger assembly (SEC+) through an
#   interface far from any protected peptide (LiP-). So the two hit lists SHOULD
#   overlap only where binding ALSO remodels the quaternary structure. The
#   question this script answers is not "do they agree" but "is the overlap they
#   DO share more than chance would give" - and what biology sits in each corner.
#
# THE FOUR CLASSES (the 2x2 SEC+/- x LiP+/-):
#   SEC+/LiP+   binding that remodels assembly - dual-validated, the strongest
#               mechanistic candidates (a ligand pocket AND a size change).
#   SEC+/LiP-   assembly change with NO local structural signature - the corner
#               that is UNIQUE TO SEC and invisible to LiP. This is the added
#               value of the elution assay: allostery/assembly through an
#               interface Piazza's peptides never saw.
#   SEC-/LiP+   classic ligand binding: a local conformational change with no
#               size change - the LARGEST overlap-adjacent class, and exactly what
#               LiP was built to find. Expected to dominate the LiP hits.
#   SEC-/LiP-   neither assay flagged it - the background.
#
# WHAT IT COMPUTES, per metabolite AND pooled:
#   1. a 2x2 contingency table over the UNIVERSE (see below) with a Fisher exact
#      test + odds ratio + fold-enrichment + a hypergeometric p for the overlap
#      cell: is co-occurrence of SEC+ and LiP+ enriched above random?
#   2. the membership of every quadrant (auditable, gene-annotated);
#   3. a GO over-representation test per quadrant (hypergeometric + BH, the SAME
#      method as the report and globularity_category_annotation.R), background =
#      the shared universe;
#   4. figures: the 2x2 mosaic per metabolite + pooled, a quadrant-count bar, the
#      SEC volcano coloured by LiP status, and per-quadrant GO bar plots.
#
# THE UNIVERSE (this is the whole game for the enrichment p-value):
#   An overlap count on its own is meaningless without the set of proteins that
#   COULD have been a hit on BOTH axes. Two modes, chosen automatically from the
#   Piazza file you give:
#     BACKGROUND MODE (you provide Piazza's full tested table with a q/p column):
#       universe[m] = {SEC-measured proteins} INTERSECT {LiP-tested proteins for m}
#       LiP+ = q < cut ; LiP- = tested but not significant. This is the rigorous
#       mode and the enrichment p-value is trustworthy.
#     HITLIST MODE (you provide only Piazza's significant hits, no background):
#       universe[m] = {SEC-measured proteins} INTERSECT {proteins present ANYWHERE
#       in the Piazza file}. LiP+ = in Piazza's hit list for m ; LiP- = the rest.
#       *** FLAGGED LIMITATION: here "LiP-" mixes "tested-not-significant" with
#       "not-tested-for-this-metabolite". The overlap COUNT is exact, but the
#       enrichment p-value is an approximation - state it as such. ***
#
# INPUTS YOU MUST PROVIDE:
#   piazza_file = the Piazza et al. 2018 LiP-SMap supplementary table, as .csv /
#   .tsv / .txt (or .xlsx if the `readxl` package is installed). The script is
#   format-tolerant: it auto-detects the accession column and either a metabolite
#   column (LONG layout: one row per protein-metabolite) or per-metabolite columns
#   (WIDE layout), and PRINTS everything it decided so you can correct it. Piazza
#   2018 = Piazza A. et al., "A Map of Protein-Metabolite Interactions Reveals
#   Principles of Chemical Communication", Cell 172:358 (2018).
#
# METABOLITE-NAME MATCHING (Piazza's naming -> your run keys):
#   Your keys are ADP, aKG, ATP, NAD, PEP, PGP. The default synonym map is taken
#   from this repo's own .INTERACTOR_REGEX (scripts/validation_candidates.R) so it
#   is faithful to how the pipeline already defines each compound. PGP is the one
#   to double-check: in this codebase PGP = 6-phosphogluconate. Every match (and
#   every UNMATCHED key, with the full list of Piazza names) is printed; override
#   with `metabolite_map =` if the file names anything differently.
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "sec_vs_lip_quadrants.R"))
#   sec_vs_lip_quadrants(piazza_file = "~/Downloads/piazza2018_LiP.csv")
#   sec_vs_lip_quadrants(piazza_file = "...", metabolites = c("ATP","ADP","aKG"))
#   sec_vs_lip_quadrants(piazza_file = "...", piazza_sig_col = "qvalue",
#                        piazza_sig_cut = 0.01)                     # force bg mode
#   sec_vs_lip_quadrants(piazza_file = "...",
#                        metabolite_map = list(PGP = c("2,3-bisphosphoglycerate")))
#
# OUTPUT (output/sec_vs_lip_quadrants/):
#   tables/contingency_per_metabolite.csv   a,b,c,d + Fisher p, OR, fold-enrichment
#   tables/contingency_pooled.csv           pair-level and protein-union pooling
#   tables/quadrant_membership.csv          every protein-metabolite pair + class
#   tables/quadrant_counts.csv              metabolite x quadrant counts
#   tables/metabolite_match_report.csv      key -> matched Piazza name(s), n hits
#   tables/GOenrichment_<quadrant>_<go>.csv per-quadrant GO over-representation
#   quadrant_interpretation.txt             the biology of each corner, in words
#   figures/contingency_mosaic.pdf          the 2x2, per metabolite + pooled
#   figures/quadrant_counts_bar.pdf         class sizes per metabolite
#   figures/sec_volcano_by_lip.pdf          SEC volcano, points coloured by LiP
#   figures/GOenrichment_quadrants_<go>.pdf per-quadrant GO bar plots
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })

# ---- house-style helpers (copied small, per this repo's per-script convention) ---------------------

# UniProt cache: accession -> gene_names, go_p/go_f/go_c, sequence (from the shared render cache).
.svl_load_uniprot <- function() {
  sf <- here("output", "uniprot_annotation_shared.RData")
  if (!file.exists(sf)) return(NULL)
  e <- new.env(); load(sf, envir = e)
  u <- tryCatch(as.data.table(e$.uniprot_all), error = function(x) NULL)
  if (is.null(u) || !"input_id" %in% names(u)) return(NULL)
  u[, accession := as.character(input_id)][]
}

# GO over-representation - identical method to globularity_category_annotation.R / the report
# (hypergeometric on the ';'-split terms, BH-adjusted).
.svl_go_enrichment <- function(foreground, background, annotation, id_col = "accession",
                               go_col = "go_p", min_genes = 2, top_n = 15) {
  fg <- unique(stats::na.omit(unlist(foreground)))
  bg <- unique(stats::na.omit(unlist(c(background, fg))))
  if (length(fg) < min_genes || length(bg) < 5) return(NULL)
  if (!all(c(id_col, go_col) %in% names(annotation))) return(NULL)
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
  ts[, fold := (k / n) / (K / N)]
  ts[order(padj, -k)][seq_len(min(top_n, .N))]
}
.svl_go_barplot <- function(go_tbl, title) {
  if (is.null(go_tbl) || !nrow(go_tbl)) return(NULL)
  gt <- copy(go_tbl); gt[, term_short := ifelse(nchar(term) > 60, paste0(substr(term, 1, 57), "..."), term)]
  ggplot(gt, aes(stats::reorder(term_short, -padj), -log10(padj), fill = k)) +
    geom_col() + coord_flip() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
    labs(title = title, x = NULL, y = "-log10(BH-adjusted p)", fill = "n proteins") + theme_bw()
}

# find a user-supplied file across the usual places (mirrors surface_hydrophobicity.R's .resolve_file);
# on failure, LIST what is actually there so the real name/path is obvious.
.svl_resolve_file <- function(f) {
  if (is.null(f) || !nzchar(f)) stop("piazza_file is required - pass the path to the Piazza 2018 LiP table.", call. = FALSE)
  cand <- unique(c(f, path.expand(f), here(f)))
  base <- basename(f)
  homes <- unique(c(Sys.getenv("USERPROFILE"), Sys.getenv("HOME"), path.expand("~"), dirname(path.expand("~"))))
  homes <- homes[nzchar(homes) & dir.exists(homes)]
  dirs  <- unique(c(homes, file.path(homes, c("Downloads", "Desktop", "Documents", "Documents/Downloads")),
                    here("data", "raw"), getwd()))
  dirs  <- dirs[dir.exists(dirs)]
  cand  <- unique(c(cand, file.path(dirs, base)))
  hit   <- cand[file.exists(cand) & !dir.exists(cand)]
  if (length(hit)) return(hit[1])
  found <- unlist(lapply(dirs, function(d) list.files(d, pattern = "\\.(csv|tsv|txt|xlsx|xls)$", full.names = TRUE, ignore.case = TRUE)))
  stop("Could not find '", base, "'.\nTried:\n  ", paste(utils::head(cand, 12), collapse = "\n  "),
       if (length(found))
         paste0("\n\nBut these table files DO exist - pass one of these instead:\n  ",
                paste(utils::head(found, 15), collapse = "\n  "))
       else paste0("\n\nNo .csv/.tsv/.xlsx in any of:\n  ", paste(dirs, collapse = "\n  "),
                   "\nNOTE on Windows R, '~' is your DOCUMENTS folder - give the full path with forward slashes."),
       call. = FALSE)
}

# UniProt accession token(s) inside an arbitrary cell ("sp|P0A9M0|..." or a bare accession both work).
.SVL_ACC_RE <- "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
.svl_accs_in <- function(x) toupper(unlist(regmatches(x, gregexpr(.SVL_ACC_RE, toupper(as.character(x))))))
.svl_first_acc <- function(x) {
  x <- as.character(x); m <- regmatches(toupper(x), regexpr(.SVL_ACC_RE, toupper(x)))
  out <- rep(NA_character_, length(x)); out[nzchar(m)] <- toupper(m[nzchar(m)]); out
}

# default metabolite synonym map - faithful to this repo's .INTERACTOR_REGEX (validation_candidates.R).
# PGP = 6-phosphogluconate in this codebase; VERIFY against the Piazza file naming and override if needed.
.SVL_METAB_SYNONYMS <- list(
  ATP = c("atp", "adenosine triphosphate", "adenosine-5'-triphosphate", "adenosine 5'-triphosphate"),
  ADP = c("adp", "adenosine diphosphate", "adenosine-5'-diphosphate", "adenosine 5'-diphosphate"),
  NAD = c("nad", "nad+", "nadh", "nicotinamide adenine dinucleotide"),
  aKG = c("akg", "a-kg", "2-oxoglutarate", "2 oxoglutarate", "oxoglutarate", "alpha-ketoglutarate",
          "alpha ketoglutarate", "a-ketoglutarate", "2-ketoglutarate", "ketoglutarate", "2og"),
  PEP = c("pep", "phosphoenolpyruvate", "phospho-enol-pyruvate", "phosphoenol pyruvate"),
  PGP = c("pgp", "6-phosphogluconate", "6-phospho-d-gluconate", "6-phospho-gluconate", "phosphogluconate",
          "6pg", "6-pg", "gluconate-6-phosphate", "gluconate 6-phosphate", "kdpg"))

# does a free-text metabolite label match a given key's synonyms? (word-ish, case-insensitive).
# \Q...\E quotes the synonym literally (perl=TRUE) so '+' / '-' etc. are not treated as regex.
.svl_metab_match <- function(label, syns) {
  lab <- tolower(trimws(as.character(label)))
  if (!nzchar(lab)) return(FALSE)
  any(vapply(syns, function(s) {
    grepl(paste0("(^|[^a-z0-9])\\Q", tolower(s), "\\E([^a-z0-9]|$)"), lab, perl = TRUE)
  }, logical(1)))
}
# map an arbitrary label to at most one user key (first key whose synonyms match); NA if none.
.svl_label_to_key <- function(label, syn_map) {
  for (k in names(syn_map)) if (.svl_metab_match(label, syn_map[[k]])) return(k)
  NA_character_
}

# ---- read the Piazza table into a tidy long form: (accession, metabolite_key, is_hit, tested) -------
# returns list(long = data.table, mode = "background"|"hitlist", present = <all accessions in file>,
#              match_report = data.table, sig_col = <name or NA>)
.svl_read_piazza <- function(piazza_file, syn_map, metabolites,
                             acc_col = NULL, metabolite_col = NULL, sig_col = NULL,
                             sig_cut = 0.05, is_hitlist = NA, verbose = TRUE) {
  path <- .svl_resolve_file(piazza_file)
  ext  <- tolower(tools::file_ext(path))
  if (verbose) message("Reading Piazza table: ", path)
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("The Piazza file is an Excel workbook but the `readxl` package is not installed.\n",
           "Either install it (install.packages('readxl')) or re-save the relevant sheet as .csv.", call. = FALSE)
    raw <- as.data.table(readxl::read_excel(path))
  } else {
    raw <- tryCatch(fread(path), error = function(e) as.data.table(utils::read.csv(path, check.names = FALSE)))
  }
  if (!nrow(raw)) stop("The Piazza file '", basename(path), "' read as empty.", call. = FALSE)
  raw <- raw[, which(!duplicated(names(raw))), with = FALSE]
  if (verbose) message("  columns: ", paste(names(raw), collapse = " | "))

  ch <- names(raw)[vapply(raw, function(v) is.character(v) || is.factor(v), logical(1))]

  # 1. accession column: the character column with the most accession-like cells
  if (is.null(acc_col)) {
    score <- vapply(ch, function(cn) mean(!is.na(.svl_first_acc(raw[[cn]]))), numeric(1))
    if (!length(score) || max(score, na.rm = TRUE) < 0.5)
      stop("Could not auto-detect a UniProt-accession column in the Piazza file.\n",
           "Columns seen: ", paste(names(raw), collapse = ", "),
           "\nPass acc_col = \"<the column with accessions>\".", call. = FALSE)
    acc_col <- names(which.max(score))
    if (verbose) message("  accession column (auto): ", acc_col, "  (", round(100 * max(score)), "% accession-like)")
  }
  raw[, .acc := .svl_first_acc(get(acc_col))]

  # 2. significance column (optional): defines BACKGROUND mode + the hit threshold
  auto_sig <- NA_character_
  if (is.null(sig_col)) {
    sig_cand <- grep("q[._ ]?val|qvalue|adj.*p|p[._ ]?adj|padj|fdr|\\bp[._ ]?val|pvalue|\\bp\\b", names(raw),
                     ignore.case = TRUE, value = TRUE)
    sig_cand <- sig_cand[vapply(sig_cand, function(cn) is.numeric(raw[[cn]]) || suppressWarnings(!all(is.na(as.numeric(as.character(raw[[cn]]))))), logical(1))]
    if (length(sig_cand)) { auto_sig <- sig_cand[1]
      if (verbose) message("  significance column (auto): ", auto_sig, "  (candidates: ", paste(sig_cand, collapse = ", "), ")") }
  } else auto_sig <- sig_col
  use_sig <- if (!is.na(auto_sig) && auto_sig %in% names(raw)) auto_sig else NA_character_

  # 3. metabolite: LONG (a label column) or WIDE (one column per metabolite)
  key_names <- names(syn_map)
  long <- NULL; layout <- NA_character_

  if (is.null(metabolite_col)) {
    # LONG candidate: character column whose distinct values map to >=2 keys
    long_score <- vapply(ch, function(cn) {
      if (identical(cn, acc_col)) return(0L)
      vals <- unique(as.character(raw[[cn]])); vals <- vals[!is.na(vals) & nzchar(vals)]
      if (!length(vals)) return(0L)
      length(unique(stats::na.omit(vapply(vals, .svl_label_to_key, character(1), syn_map = syn_map))))
    }, integer(1))
    if (length(long_score) && max(long_score) >= 2L) { metabolite_col <- names(which.max(long_score)); layout <- "long" }
  } else layout <- "long"

  if (identical(layout, "long")) {
    if (verbose) message("  LONG layout; metabolite column: ", metabolite_col)
    d <- raw[!is.na(.acc)]
    d[, .key := vapply(as.character(get(metabolite_col)), .svl_label_to_key, character(1), syn_map = syn_map)]
    d <- d[!is.na(.key)]
    if (!is.na(use_sig)) {
      sigv <- suppressWarnings(as.numeric(as.character(d[[use_sig]])))
      d[, .is_hit := is.finite(sigv) & sigv < sig_cut]
      d[, .tested := is.finite(sigv)]
    } else { d[, .is_hit := TRUE][, .tested := TRUE] }   # hitlist: every row present is a hit
    long <- unique(d[, .(accession = .acc, metabolite = .key, is_hit = .is_hit, tested = .tested)])
  } else {
    # WIDE: columns whose NAME matches a key's synonyms; a cell is a hit per its type
    name_key <- setNames(lapply(names(raw), .svl_label_to_key, syn_map = syn_map), names(raw))
    wide_cols <- names(raw)[!vapply(name_key, is.na, logical(1)) & names(raw) != acc_col]
    if (length(unique(unlist(name_key[wide_cols]))) < 1L)
      stop("Could not find a metabolite column (LONG) or per-metabolite columns (WIDE) in the Piazza file.\n",
           "Columns seen: ", paste(names(raw), collapse = ", "),
           "\nPass metabolite_col = \"<label column>\", or rename/point the per-metabolite columns.", call. = FALSE)
    layout <- "wide"
    if (verbose) message("  WIDE layout; per-metabolite columns: ",
                         paste(sprintf("%s->%s", wide_cols, unlist(name_key[wide_cols])), collapse = ", "))
    parts <- lapply(wide_cols, function(cn) {
      k <- name_key[[cn]]; v <- raw[[cn]]
      num <- suppressWarnings(as.numeric(as.character(v)))
      if (!is.na(use_sig) || (is.numeric(v) && mean(num >= 0 & num <= 1, na.rm = TRUE) > 0.8)) {
        hit <- is.finite(num) & num < sig_cut; tested <- is.finite(num)          # q/p-value cells
      } else if (is.logical(v)) { hit <- isTRUE(v) | (!is.na(v) & v); tested <- !is.na(v)
      } else { s <- trimws(as.character(v)); hit <- !is.na(s) & nzchar(s) & !(tolower(s) %in% c("na","0","false","-","."))
               tested <- !is.na(v) }                                              # non-empty = a called hit
      data.table(accession = raw$.acc, metabolite = k, is_hit = hit, tested = tested)
    })
    long <- unique(rbindlist(parts)[!is.na(accession)])
  }

  # decide mode. BACKGROUND = the file actually tells us which proteins were tested-but-NOT-a-hit
  # (a significance column, a 0/1 or TRUE/FALSE matrix, ...). Detected data-driven: does any row
  # carry tested==TRUE while is_hit==FALSE? HITLIST = only hits are known.
  has_bg <- isTRUE(any(long$tested & !long$is_hit, na.rm = TRUE))
  mode <- if (!is.na(is_hitlist)) (if (isTRUE(is_hitlist)) "hitlist" else "background") else
          if (!is.na(use_sig) || has_bg) "background" else "hitlist"
  if (mode == "hitlist") long[, tested := is_hit]   # in hitlist mode "tested" is unknown -> only hits are known

  # keep only the requested run keys
  long <- long[metabolite %in% metabolites]

  # match report + all Piazza labels (so unmatched keys can be fixed)
  present <- unique(long$accession)
  all_present <- unique(.svl_first_acc(raw[[acc_col]])); all_present <- all_present[!is.na(all_present)]
  report <- rbindlist(lapply(metabolites, function(k) {
    data.table(key = k,
               matched = k %in% long$metabolite,
               n_lip_hits = long[metabolite == k & is_hit, uniqueN(accession)],
               n_lip_tested = if (mode == "background") long[metabolite == k & tested, uniqueN(accession)] else NA_integer_)
  }))
  attr(long, "all_present") <- all_present
  list(long = long, mode = mode, present = all_present, match_report = report,
       sig_col = use_sig, layout = layout, acc_col = acc_col,
       raw_labels = if (identical(layout, "long")) sort(unique(as.character(raw[[metabolite_col]]))) else names(raw))
}

# ---- load this study's SEC results: per metabolite, the measured set and the hit set ----------------
.svl_load_sec <- function(metabolites, pBHadj_cut, log2fc_cut, verbose = TRUE) {
  out <- list()
  for (m in metabolites) {
    f <- here("output", paste0("PCM_ctrl_vs_", m), "rdata", "protein_DiffExprProtein_list.RData")
    if (!file.exists(f)) { if (verbose) message("  [", m, "] SEC protein_DiffExprProtein_list.RData missing -> skipped."); next }
    e <- new.env(); load(f, envir = e)
    pdp <- as.data.table(e$protein_DiffExprProtein)
    if (!all(c("feature_id", "pBHadj", "medianLog2FC") %in% names(pdp))) {
      if (verbose) message("  [", m, "] protein_DiffExprProtein missing expected columns -> skipped."); next }
    pdp[, feature_id := as.character(feature_id)]
    measured <- unique(pdp$feature_id)
    hits     <- unique(pdp[pBHadj < pBHadj_cut & abs(medianLog2FC) > log2fc_cut, feature_id])
    stat     <- pdp[, .(metabolite = m, feature_id, medianLog2FC, pBHadj,
                        gene = if ("gene_names" %in% names(pdp)) gene_names else NA_character_)]
    out[[m]] <- list(measured = measured, hits = hits, stat = unique(stat))
    if (verbose) message("  [", m, "] SEC: ", length(measured), " measured, ", length(hits), " hits.")
  }
  out
}

# =====================================================================================================
sec_vs_lip_quadrants <- function(
    piazza_file      = NULL,
    metabolites      = NULL,                       # default: every PCM_ctrl_vs_* run present
    metabolite_map   = NULL,                       # override/extend the synonym map (named list)
    acc_col          = NULL,                       # Piazza accession column (auto if NULL)
    metabolite_col   = NULL,                       # Piazza metabolite label column (auto if NULL)
    sig_col          = NULL,                       # Piazza q/p column -> BACKGROUND mode (auto-detected)
    sig_cut          = 0.05,
    is_hitlist       = NA,                          # force TRUE/FALSE; NA = decide from sig column
    sec_pBHadj_cut   = 0.05,
    sec_log2fc_cut   = 1,
    enrichment_alt   = c("greater", "two.sided", "less"),
    go_columns       = c("go_p", "go_f", "go_c"),
    go_min_genes     = 2, go_top_n = 15,
    out_subdir       = "sec_vs_lip_quadrants",
    verbose          = TRUE) {

  enrichment_alt <- match.arg(enrichment_alt)
  outdir <- here("output", out_subdir)
  tdir <- file.path(outdir, "tables"); fdir <- file.path(outdir, "figures")
  dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)

  syn_map <- .SVL_METAB_SYNONYMS
  if (!is.null(metabolite_map)) for (k in names(metabolite_map)) syn_map[[k]] <- unique(c(syn_map[[k]], metabolite_map[[k]]))

  # which metabolites: default = the SEC runs present on disk
  if (is.null(metabolites)) {
    dd <- list.dirs(here("output"), recursive = FALSE)
    metabolites <- sub("^PCM_ctrl_vs_", "", basename(dd)[grepl("^PCM_ctrl_vs_", basename(dd))])
    if (!length(metabolites)) metabolites <- c("PGP", "NAD", "aKG", "ATP", "ADP", "PEP")
  }
  metabolites <- unique(metabolites)
  message("Metabolites: ", paste(sort(metabolites), collapse = ", "))

  # ---- load both assays ----------------------------------------------------------------------------
  message("Loading this study's SEC results ...")
  SEC <- .svl_load_sec(metabolites, sec_pBHadj_cut, sec_log2fc_cut, verbose)
  if (!length(SEC)) stop("No SEC results loaded - render at least one PCM_ctrl_vs_<metabolite> first.", call. = FALSE)

  P <- .svl_read_piazza(piazza_file, syn_map, metabolites, acc_col, metabolite_col, sig_col,
                        sig_cut, is_hitlist, verbose)
  LIP <- P$long
  message("Piazza mode: ", toupper(P$mode),
          if (P$mode == "hitlist")
            "  (LiP- mixes tested-not-significant with not-tested; overlap COUNTS exact, enrichment p approximate)"
          else if (!is.na(P$sig_col)) paste0("  (significance column: ", P$sig_col, ", cut < ", sig_cut, ")")
          else "  (tested/non-hit background read from the per-metabolite columns)")

  # metabolite match report
  fwrite(P$match_report, file.path(tdir, "metabolite_match_report.csv"))
  unmatched <- P$match_report[matched == FALSE, key]
  if (length(unmatched)) {
    message("\n*** METABOLITES NOT MATCHED in the Piazza file: ", paste(unmatched, collapse = ", "), " ***")
    message("    Distinct metabolite labels present in the file:")
    message("      ", paste(utils::head(P$raw_labels, 60), collapse = " | "))
    message("    Fix with metabolite_map = list(", unmatched[1], " = c(\"<name as written in the file>\")).\n")
  }
  if (verbose) { message("Metabolite match report:"); print(P$match_report) }

  # keep metabolites that exist on BOTH sides
  both <- intersect(names(SEC), P$match_report[matched == TRUE, key])
  if (!length(both)) stop("No metabolite is present in BOTH the SEC runs and the Piazza file after name-matching.\n",
                          "See tables/metabolite_match_report.csv and pass metabolite_map=.", call. = FALSE)
  message("Metabolites usable on both axes: ", paste(sort(both), collapse = ", "))

  U <- .svl_load_uniprot()
  gene_of <- function(ids) {
    g <- rep(NA_character_, length(ids))
    if (!is.null(U) && "gene_names" %in% names(U)) g <- U$gene_names[match(ids, U$accession)]
    # fall back to SEC-provided gene names
    sg <- rbindlist(lapply(SEC, function(x) x$stat[, .(feature_id, gene)]), fill = TRUE)
    sg <- unique(sg[!is.na(gene)]); g2 <- sg$gene[match(ids, sg$feature_id)]
    ifelse(is.na(g) | !nzchar(g), g2, g)
  }

  # ---- per-metabolite quadrant assignment + 2x2 ----------------------------------------------------
  lip_present_all <- P$present                       # LiP-observable proteome proxy (hitlist mode)
  qrows <- list(); crows <- list()
  quad_of <- function(sec, lip) fifelse(sec & lip, "SEC+/LiP+",
                                fifelse(sec & !lip, "SEC+/LiP-",
                                fifelse(!sec & lip, "SEC-/LiP+", "SEC-/LiP-")))

  for (m in both) {
    sec_meas <- SEC[[m]]$measured; sec_hit <- SEC[[m]]$hits
    lm <- LIP[metabolite == m]
    lip_hit    <- unique(lm[is_hit == TRUE, accession])
    lip_tested <- if (P$mode == "background") unique(lm[tested == TRUE, accession]) else lip_present_all

    universe <- if (P$mode == "background") intersect(sec_meas, lip_tested)
                else intersect(sec_meas, lip_present_all)
    if (!length(universe)) { message("  [", m, "] empty universe (no SEC-measured protein overlaps the LiP set) -> skipped."); next }

    sec <- universe %in% sec_hit
    lip <- universe %in% lip_hit
    q   <- quad_of(sec, lip)
    st  <- SEC[[m]]$stat
    qrows[[m]] <- data.table(metabolite = m, accession = universe, gene = gene_of(universe),
                             SEC_hit = sec, LiP_hit = lip, quadrant = q,
                             medianLog2FC = st$medianLog2FC[match(universe, st$feature_id)],
                             SEC_pBHadj  = st$pBHadj[match(universe, st$feature_id)])

    a <- sum(sec & lip); b <- sum(sec & !lip); c <- sum(!sec & lip); d <- sum(!sec & !lip)
    n <- a + b + c + d
    exp_a <- (a + b) * (a + c) / n
    ft <- tryCatch(stats::fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = enrichment_alt),
                   error = function(e) NULL)
    hyp <- stats::phyper(a - 1, a + c, b + d, a + b, lower.tail = FALSE)   # over-rep of the overlap cell
    crows[[m]] <- data.table(
      metabolite = m, mode = P$mode, n_universe = n,
      SECpos = a + b, LiPpos = a + c,
      a_SECp_LiPp = a, b_SECp_LiPm = b, c_SECm_LiPp = c, d_SECm_LiPm = d,
      expected_overlap = round(exp_a, 2), fold_enrichment = if (exp_a > 0) round(a / exp_a, 2) else NA_real_,
      odds_ratio = if (!is.null(ft)) round(unname(ft$estimate), 3) else NA_real_,
      fisher_p = if (!is.null(ft)) signif(ft$p.value, 3) else NA_real_,
      hyper_p  = signif(hyp, 3))
    message(sprintf("  [%s] universe=%d | SEC+ %d, LiP+ %d | overlap a=%d (exp %.1f, %.2gx) | Fisher p=%.2g",
                    m, n, a + b, a + c, a, exp_a, if (exp_a > 0) a / exp_a else NA_real_,
                    if (!is.null(ft)) ft$p.value else NA_real_))
  }
  if (!length(qrows)) stop("No metabolite produced a non-empty quadrant table.", call. = FALSE)

  QUAD <- rbindlist(qrows, fill = TRUE)
  CONT <- rbindlist(crows, fill = TRUE)
  fwrite(QUAD, file.path(tdir, "quadrant_membership.csv"))
  fwrite(CONT, file.path(tdir, "contingency_per_metabolite.csv"))

  # quadrant counts (metabolite x quadrant)
  qlev <- c("SEC+/LiP+", "SEC+/LiP-", "SEC-/LiP+", "SEC-/LiP-")
  QC <- dcast(QUAD[, .N, by = .(metabolite, quadrant)], metabolite ~ quadrant, value.var = "N", fill = 0)
  for (q in qlev) if (!q %in% names(QC)) QC[[q]] <- 0L
  setcolorder(QC, c("metabolite", qlev)); fwrite(QC, file.path(tdir, "quadrant_counts.csv"))
  if (verbose) { message("\nQuadrant counts per metabolite:"); print(QC) }

  # ---- pooling ------------------------------------------------------------------------------------
  # (i) PAIR-LEVEL: each protein-metabolite pair is one observation (the correct pooled test)
  poolA <- QUAD[, .(a = sum(SEC_hit & LiP_hit), b = sum(SEC_hit & !LiP_hit),
                    c = sum(!SEC_hit & LiP_hit), d = sum(!SEC_hit & !LiP_hit))]
  a <- poolA$a; b <- poolA$b; c <- poolA$c; d <- poolA$d; n <- a + b + c + d
  exp_a <- (a + b) * (a + c) / n
  ftp <- tryCatch(stats::fisher.test(matrix(c(a, b, c, d), 2, byrow = TRUE), alternative = enrichment_alt), error = function(e) NULL)
  pair_row <- data.table(pooling = "pairs (protein x metabolite)", n_universe = n,
                         a_SECp_LiPp = a, b_SECp_LiPm = b, c_SECm_LiPp = c, d_SECm_LiPm = d,
                         expected_overlap = round(exp_a, 2), fold_enrichment = if (exp_a > 0) round(a / exp_a, 2) else NA_real_,
                         odds_ratio = if (!is.null(ftp)) round(unname(ftp$estimate), 3) else NA_real_,
                         fisher_p = if (!is.null(ftp)) signif(ftp$p.value, 3) else NA_real_)

  # (ii) PROTEIN-LEVEL UNION: SEC+ = hit for ANY metabolite; LiP+ = Piazza hit for ANY metabolite
  #      (reproduces the study-vs-Piazza Venn number: this-study / Piazza / overlap)
  uni_universe <- unique(QUAD$accession)
  sec_any <- unique(QUAD[SEC_hit == TRUE, accession])
  lip_any <- unique(QUAD[LiP_hit == TRUE, accession])
  A <- length(intersect(sec_any, lip_any)); B <- length(setdiff(sec_any, lip_any))
  Cc <- length(setdiff(lip_any, sec_any)); D <- length(setdiff(uni_universe, union(sec_any, lip_any)))
  nU <- A + B + Cc + D; exp_A <- (A + B) * (A + Cc) / nU
  ftu <- tryCatch(stats::fisher.test(matrix(c(A, B, Cc, D), 2, byrow = TRUE), alternative = enrichment_alt), error = function(e) NULL)
  union_row <- data.table(pooling = "proteins (union over metabolites)", n_universe = nU,
                          a_SECp_LiPp = A, b_SECp_LiPm = B, c_SECm_LiPp = Cc, d_SECm_LiPm = D,
                          expected_overlap = round(exp_A, 2), fold_enrichment = if (exp_A > 0) round(A / exp_A, 2) else NA_real_,
                          odds_ratio = if (!is.null(ftu)) round(unname(ftu$estimate), 3) else NA_real_,
                          fisher_p = if (!is.null(ftu)) signif(ftu$p.value, 3) else NA_real_)
  POOL <- rbindlist(list(pair_row, union_row), fill = TRUE)
  fwrite(POOL, file.path(tdir, "contingency_pooled.csv"))
  message("\nProtein-level Venn (union over metabolites, within the shared universe):")
  message(sprintf("  SEC hits (this study) = %d | LiP hits (Piazza) = %d | overlap = %d | SEC-only = %d | LiP-only = %d",
                  length(sec_any), length(lip_any), A, B, Cc))
  message(sprintf("  overlap fold-enrichment = %.2gx (expected %.1f) | Fisher p = %.2g",
                  if (exp_A > 0) A / exp_A else NA_real_, exp_A, if (!is.null(ftu)) ftu$p.value else NA_real_))
  if (verbose) { message("Pooled contingency:"); print(POOL) }

  # ---- GO over-representation per quadrant (pooled pairs; bg = pooled universe) ---------------------
  go_written <- character(0)
  qtoken <- c("SEC+/LiP+" = "SECpos_LiPpos", "SEC+/LiP-" = "SECpos_LiPneg",
              "SEC-/LiP+" = "SECneg_LiPpos", "SEC-/LiP-" = "SECneg_LiPneg")
  if (!is.null(U)) {
    bg <- unique(QUAD$accession)
    for (go in go_columns) {
      if (!go %in% names(U)) next
      plots <- list()
      for (q in qlev) {
        fg <- unique(QUAD[quadrant == q, accession])
        tbl <- .svl_go_enrichment(fg, bg, U, id_col = "accession", go_col = go,
                                  min_genes = go_min_genes, top_n = go_top_n)
        if (!is.null(tbl) && nrow(tbl)) {
          fwrite(tbl, file.path(tdir, paste0("GOenrichment_", qtoken[[q]], "_", go, ".csv")))
          plots[[q]] <- .svl_go_barplot(tbl, paste0(q, "  -  ", go, "  (n=", length(fg), ", bg=", length(bg), ")"))
        }
      }
      if (length(plots)) {
        grDevices::pdf(file.path(fdir, paste0("GOenrichment_quadrants_", go, ".pdf")), width = 9, height = 6)
        for (p in plots) print(p); grDevices::dev.off()
        go_written <- c(go_written, go)
      }
    }
  } else message("NOTE: no output/uniprot_annotation_shared.RData -> GO enrichment skipped (render a comparison first).")

  # ---- figures -------------------------------------------------------------------------------------
  # (a) 2x2 mosaic, per metabolite + pooled-pairs, cell = count, fill = fold-enrichment
  mos <- rbindlist(lapply(both, function(m) {
    cc <- CONT[metabolite == m]; if (!nrow(cc)) return(NULL)
    data.table(panel = sprintf("%s (Fisher p=%.2g)", m, cc$fisher_p),
               SEC = c("SEC+", "SEC+", "SEC-", "SEC-"), LiP = c("LiP+", "LiP-", "LiP+", "LiP-"),
               n = c(cc$a_SECp_LiPp, cc$b_SECp_LiPm, cc$c_SECm_LiPp, cc$d_SECm_LiPm))
  }), fill = TRUE)
  mos_pool <- data.table(panel = sprintf("POOLED pairs (p=%.2g)", pair_row$fisher_p),
                         SEC = c("SEC+", "SEC+", "SEC-", "SEC-"), LiP = c("LiP+", "LiP-", "LiP+", "LiP-"),
                         n = c(a, b, c, d))
  mos <- rbindlist(list(mos, mos_pool), fill = TRUE)
  mos[, frac := n / sum(n), by = panel]
  gm <- ggplot(mos, aes(LiP, factor(SEC, levels = c("SEC-", "SEC+")), fill = frac)) +
    geom_tile(colour = "white") + geom_text(aes(label = n), size = 4) +
    facet_wrap(~ panel) + scale_fill_gradient(low = "grey95", high = "firebrick") +
    labs(title = "SEC x LiP contingency", subtitle = "cell = protein count; fill = share of that panel",
         x = NULL, y = NULL, fill = "share") + theme_bw() + theme(legend.position = "right")
  ggsave(file.path(fdir, "contingency_mosaic.pdf"), gm, width = 10, height = 8)

  # (b) quadrant-count bars per metabolite
  QL <- melt(QC, id.vars = "metabolite", variable.name = "quadrant", value.name = "n")
  QL[, quadrant := factor(quadrant, levels = qlev)]
  gb <- ggplot(QL, aes(metabolite, n, fill = quadrant)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("SEC+/LiP+" = "#b2182b", "SEC+/LiP-" = "#ef8a62",
                                 "SEC-/LiP+" = "#67a9cf", "SEC-/LiP-" = "grey80")) +
    labs(title = "Quadrant sizes per metabolite", x = NULL, y = "proteins in the shared universe") +
    theme_bw()
  ggsave(file.path(fdir, "quadrant_counts_bar.pdf"), gb, width = 9, height = 5)

  # (c) SEC volcano, coloured by LiP status - shows WHERE the LiP+ proteins sit in SEC space
  V <- QUAD[is.finite(medianLog2FC) & is.finite(SEC_pBHadj)]
  V[, negLogP := -log10(pmax(SEC_pBHadj, .Machine$double.xmin))]
  gv <- ggplot(V[order(LiP_hit)], aes(medianLog2FC, negLogP, colour = LiP_hit)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_hline(yintercept = -log10(sec_pBHadj_cut), linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = c(-sec_log2fc_cut, sec_log2fc_cut), linetype = "dashed", colour = "grey60") +
    scale_colour_manual(values = c(`FALSE` = "grey75", `TRUE` = "#b2182b"), name = "LiP hit (Piazza)") +
    facet_wrap(~ metabolite, scales = "free_y") +
    labs(title = "SEC differential (elution) volcano, coloured by LiP status",
         subtitle = "red = also a Piazza LiP hit for that metabolite; dashed = SEC hit thresholds",
         x = "SEC medianLog2FC", y = "-log10(SEC pBHadj)") + theme_bw()
  ggsave(file.path(fdir, "sec_volcano_by_lip.pdf"), gv, width = 11, height = 7)

  # ---- interpretation note (handover) --------------------------------------------------------------
  interp <- c(
    "SEC x LiP quadrant analysis - what each corner means",
    "=====================================================",
    "",
    sprintf("Mode: %s.", toupper(P$mode)),
    if (P$mode == "hitlist")
      "  HITLIST mode: only Piazza's significant hits were available, so 'LiP-' mixes 'tested-not-significant' with 'not-tested-for-this-metabolite'. The overlap COUNTS are exact; the enrichment p-values are approximate (interpret the direction and rough magnitude, not the exact value)."
    else sprintf("  BACKGROUND mode: Piazza's full tested table was used (sig column '%s', cut < %s), so LiP- = tested-but-not-significant and the enrichment p-values are trustworthy.", P$sig_col, sig_cut),
    "",
    "SEC+/LiP+  binding that ALSO remodels assembly. Dual-validated: a local structural",
    "           change (protease protection) AND a size/assembly change on the column.",
    "           The strongest mechanistic candidates.",
    "SEC+/LiP-  assembly change with NO local structural signature. UNIQUE TO SEC and",
    "           invisible to LiP - the added value of the elution assay. Assembly driven",
    "           through an interface Piazza's peptides never sampled.",
    "SEC-/LiP+  classic ligand binding: a local conformational change with no size change.",
    "           Expected to be the LARGEST hit-bearing class - it is what LiP is built for.",
    "SEC-/LiP-  neither assay flagged it (the background).",
    "",
    "Why 'mostly non-overlapping' is the expected result, not a failure:",
    "  SEC reads QUATERNARY state; LiP reads LOCAL structure. They are orthogonal, so they",
    "  agree only where a binding event also changes the particle's size. The Fisher test",
    "  asks whether the overlap they DO share (the SEC+/LiP+ corner) is larger than random",
    "  co-occurrence would give: enrichment there means binding sometimes remodels assembly,",
    "  even though most binding (SEC-/LiP+) and most assembly change (SEC+/LiP-) do not coincide.",
    "",
    "Where the pipeline already touches this descriptively (annotation-based, NOT the Piazza data):",
    "  - scripts/hit_signature_analysis.R  : flags SEC hits that are annotated metabolite-binders",
    "                                         (.INTERACTOR_REGEX) - an annotation proxy for 'LiP+'.",
    "  - scripts/validation_candidates.R   : ranks SEC/EMD/CCF candidates and marks known binders",
    "                                         as positive controls.",
    "  This script is the first to cross the SEC hits against Piazza's EXPERIMENTAL LiP data.")
  writeLines(interp, file.path(outdir, "quadrant_interpretation.txt"))

  message("\nDone. Written to: ", outdir)
  message("  tables/   contingency_per_metabolite.csv, contingency_pooled.csv, quadrant_membership.csv,")
  message("            quadrant_counts.csv, metabolite_match_report.csv, GOenrichment_*.csv")
  message("  figures/  contingency_mosaic.pdf, quadrant_counts_bar.pdf, sec_volcano_by_lip.pdf,",
          if (length(go_written)) paste0(" GOenrichment_quadrants_{", paste(go_written, collapse = ","), "}.pdf") else "")
  message("  quadrant_interpretation.txt")
  invisible(list(contingency = CONT, pooled = POOL, quadrants = QUAD, counts = QC,
                 match_report = P$match_report, mode = P$mode))
}
