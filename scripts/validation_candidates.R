# scripts/validation_candidates.R
# =============================================================================
# SHORTLIST proteins worth validating orthogonally with a PHOTOMETRIC ACTIVITY ASSAY for allosteric
# regulation by the metabolite that moved them.
#
# It pools the two independent lines of evidence already produced by this pipeline and scores every
# candidate on whether it is (a) well supported, (b) a relevant / central-metabolism enzyme, and
# (c) actually assayable in a cuvette or plate reader.
#
# EVIDENCE POOLED (per metabolite):
#   stat hit    CCprofiler protein-level differential test - pBHadj < 0.05 AND |medianLog2FC| > 1
#               (the same definition as overlap_between_metabolites.R). An ABUNDANCE change.
#   CCF / EMD   the top-N proteins of the permutation shift screens (smallest q in ccf_fdr_results.txt /
#               emd_fdr_results.txt). An ELUTION-PROFILE change. These screens returned no
#               FDR-significant hits, so they are used as a RANKED EXPLORATORY list, never as
#               "significant" - a protein appearing in both lines is the strongest case.
#
# SCORING (all weights are arguments - change them and re-run):
#   evidence     stat hit, presence/rank in the shift screens, effect size, and a bonus for BOTH lines
#   pathway      central-metabolism / relevant-enzyme annotation (GO + pathway keywords)
#   assayability the decisive practical filter: an EC-class-driven guess at a PHOTOMETRIC readout
#                (NAD(P)H at 340 nm directly or via a coupled system, pNPP at 405 nm, A240 lyase
#                assays, DTNB at 412 nm ...). Enzymes with no obvious photometric route score 0 and
#                drop out, which is the point of the exercise.
#   allosteric   a small CURATED prior of classic allosterically regulated E. coli enzymes (see
#                .ALLOSTERIC_PRIOR) - textbook knowledge, flagged as such, to be verified, not trusted.
#   heteromer    Complex Portal subunit composition: candidates that only work as heteromers are NOT
#                excluded (they are often the interesting ones) but are flagged with the partners that
#                must be co-purified and co-added, since that changes the experiment's cost.
#
# IMPORTANT - WHAT THIS IS AND IS NOT: this is an ANNOTATION-DRIVEN shortlist to start a conversation
# with the bench, not an assay protocol. The EC -> assay mapping is a heuristic; substrate availability,
# Km, coupling-enzyme interference, expression and purification feasibility are NOT modelled. Verify
# every suggested assay against the primary literature (e.g. BRENDA / EcoCyc) before ordering anything.
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "validation_candidates.R"))
#   validation_candidates()                                  # 10 per metabolite, all metabolites
#   validation_candidates("ATP", n_per_metabolite = 5)
#   validation_candidates(ccf_top_n = 10)                    # stricter shift screen
#   validation_candidates(require_photometric = FALSE)       # keep enzymes without an obvious assay
#
# OUTPUT (output/validation_candidates/):
#   validation_candidates_<metabolite>.txt   the full scored table for that metabolite
#   validation_candidates_all.csv            every candidate, every metabolite
#   validation_shortlist.md                  the readable shortlist: per metabolite, the top N with
#                                            gene, evidence, suggested assay and partner requirements
#   validation_shortlist.pdf                 score composition of the shortlisted proteins
#   nucleic_acid_binders_<metabolite>.txt    DNA/RNA binders, set aside from the shortlist
#   nucleic_acid_binders_all.csv             the same, combined
#   nucleic_acid_binders.md                  readable list of the nucleic-acid binders that shift, with
#                                            the matched terms and a note on how to follow them up
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })

# reuse .load_uniprot / .complex_portal_members (that file only defines functions)
source(here::here("scripts", "globularity_category_annotation.R"))

# ---- photometric assay routes by EC class (HEURISTIC - verify against BRENDA/EcoCyc) ---------------
# each entry: EC prefix regex -> suggested readout + feasibility weight (3 = direct/standard coupled,
# 1.5 = feasible but needs a coupling system or a less common substrate, 0 = no obvious photometric route)
.EC_ASSAY <- list(
  list(rx = "^1\\.1\\.1\\.",  assay = "NAD(P)H absorbance at 340 nm (direct dehydrogenase assay)",        w = 3),
  list(rx = "^1\\.2\\.1\\.",  assay = "NAD(P)H at 340 nm (aldehyde/oxo-acid dehydrogenase)",              w = 3),
  list(rx = "^1\\.6\\.",      assay = "NAD(P)H at 340 nm (NAD(P)H-dependent reductase)",                  w = 3),
  list(rx = "^1\\.8\\.",      assay = "DTNB at 412 nm (thiol) or NAD(P)H at 340 nm",                      w = 2),
  list(rx = "^1\\.",          assay = "oxidoreductase: NAD(P)H 340 nm, or DCPIP 600 nm if flavin-linked", w = 2),
  list(rx = "^2\\.7\\.",      assay = "kinase: pyruvate-kinase/LDH-coupled NADH oxidation at 340 nm",     w = 3),
  list(rx = "^2\\.6\\.1\\.",  assay = "aminotransferase: MDH- or LDH-coupled NADH oxidation at 340 nm",   w = 3),
  list(rx = "^2\\.3\\.1\\.",  assay = "acyltransferase: DTNB at 412 nm (free CoA-SH release)",            w = 2.5),
  list(rx = "^2\\.",          assay = "transferase: usually needs a coupled NAD(P)H or chromogenic system", w = 1.5),
  list(rx = "^3\\.1\\.3\\.",  assay = "phosphatase: pNPP hydrolysis at 405 nm",                           w = 3),
  list(rx = "^3\\.1\\.",      assay = "esterase/nuclease: chromogenic substrate (e.g. pNP-ester) at 405 nm", w = 2),
  list(rx = "^3\\.5\\.",      assay = "amidohydrolase: coupled (e.g. GDH/NADH at 340 nm) or chromogenic", w = 1.5),
  list(rx = "^3\\.6\\.",      assay = "ATPase/GTPase: PK/LDH-coupled NADH oxidation at 340 nm",           w = 3),
  list(rx = "^4\\.1\\.2\\.",  assay = "aldolase: TIM/GDH-coupled NADH oxidation at 340 nm",               w = 3),
  list(rx = "^4\\.2\\.1\\.",  assay = "hydro-lyase: direct A240 (fumarase/enolase-type) or coupled 340 nm", w = 2.5),
  list(rx = "^4\\.1\\.1\\.",  assay = "decarboxylase: coupled dehydrogenase, NAD(P)H at 340 nm",          w = 2.5),
  list(rx = "^4\\.",          assay = "lyase: direct A240 for conjugated products, else coupled 340 nm",  w = 1.5),
  list(rx = "^5\\.3\\.1\\.",  assay = "isomerase: coupled dehydrogenase, NAD(P)H at 340 nm",              w = 2.5),
  list(rx = "^5\\.",          assay = "isomerase/mutase: coupled NAD(P)H system at 340 nm",               w = 1.5),
  list(rx = "^6\\.",          assay = "ligase: ATP consumption via PK/LDH-coupled NADH at 340 nm",        w = 2.5))

# first regex match per element, or NA - regmatches(x, regexpr(...)) returns ONLY the matched elements,
# so it must be scattered back to the original positions rather than used directly.
.extract_first <- function(x, pattern) {
  x  <- as.character(x); x[is.na(x)] <- ""
  mt <- regexpr(pattern, x)
  out <- rep(NA_character_, length(x))
  if (any(mt > 0)) out[mt > 0] <- regmatches(x, mt)
  out
}

.assay_for <- function(ec, text) {
  if (!is.na(ec) && nzchar(ec)) for (a in .EC_ASSAY) if (grepl(a$rx, ec)) return(list(assay = a$assay, w = a$w))
  # no EC number: fall back on functional wording
  if (grepl("dehydrogenase|reductase|oxidoreductase", text, ignore.case = TRUE))
    return(list(assay = "NAD(P)H at 340 nm (inferred from name - EC not annotated)", w = 2))
  if (grepl("kinase|phosphotransferase", text, ignore.case = TRUE))
    return(list(assay = "PK/LDH-coupled NADH oxidation at 340 nm (inferred from name)", w = 2))
  if (grepl("phosphatase", text, ignore.case = TRUE))
    return(list(assay = "pNPP at 405 nm (inferred from name)", w = 2))
  if (grepl("synthetase|synthase|ligase|ATPase", text, ignore.case = TRUE))
    return(list(assay = "coupled ATP-consumption assay, NADH at 340 nm (inferred from name)", w = 1.5))
  list(assay = NA_character_, w = 0)
}

# central-metabolism / relevant-enzyme wording (GO + pathway text)
.CENTRAL_REGEX <- paste0("glycoly|gluconeogen|tricarboxylic|citrate cycle|TCA|pentose.?phosphate|",
                         "pyruvate metabol|oxidative phosphoryl|respiratory|electron transport|",
                         "glyoxylate|acetate metabol|fermentat|nucleotide biosynth|purine|pyrimidine|",
                         "amino.?acid biosynth|one.?carbon|folate|NAD biosynth|carbon metabol|",
                         "carbohydrate metabol|energy")

# DNA / RNA binders - set aside from the main prioritisation. A nucleic-acid-associated protein travels
# with DNA, RNA or the ribosome, so a change in its elution can reflect altered nucleic-acid association
# or ribosome engagement rather than an allosteric effect on the protein itself, and a photometric assay
# on a purified subunit is rarely meaningful for them. They are NOT discarded: they get their own list,
# since a metabolite-dependent shift of a nucleic-acid binder is interesting in its own right
# (ppGpp-like signalling, riboswitch-like behaviour, nucleoid remodelling, ribosome assembly).
# NOTE the pattern deliberately does NOT contain "nucleotide": a P-loop nucleotide-binding enzyme is not
# a nucleic-acid binder. HEURISTIC - override with `nucleic_acid_regex`; matched terms are written out.
.NUCLEIC_REGEX <- paste0("DNA.?binding|RNA.?binding|nucleic acid|\\brRNA|\\btRNA|\\bmRNA|ribosom|",
                         "translation|transcription|nucleoid|helicase|topoisomerase|polymerase|",
                         "ribonucleoprotein|sigma factor|nuclease|primase|recombinase|transposase|",
                         "chromosome|DNA repair|DNA replication|riboswitch|antitermination")

# CURATED PRIOR - classic allosterically regulated E. coli enzymes (textbook knowledge; VERIFY, do not
# treat as evidence from this dataset). Gene symbols as used by UniProt gene_names.
.ALLOSTERIC_PRIOR <- c(
  pfkA = "PFK-1: activated by ADP/GDP, inhibited by PEP",
  pfkB = "PFK-2 (minor isozyme)",
  pykF = "pyruvate kinase I: activated by fructose-1,6-bisP",
  pykA = "pyruvate kinase II: activated by AMP/ribose-5-P",
  ppc  = "PEP carboxylase: activated by acetyl-CoA and FBP, inhibited by aspartate",
  gltA = "citrate synthase: inhibited by NADH and 2-oxoglutarate",
  icd  = "isocitrate dehydrogenase: regulated by phosphorylation; NADP-linked",
  zwf  = "glucose-6-P dehydrogenase: NADPH feedback",
  fbp  = "fructose-1,6-bisphosphatase: inhibited by AMP",
  pyrB = "aspartate transcarbamoylase catalytic subunit: CTP/ATP regulated (heteromer with pyrI)",
  pyrI = "aspartate transcarbamoylase regulatory subunit (heteromer with pyrB)",
  glnA = "glutamine synthetase: cumulative feedback inhibition; adenylylation",
  gdhA = "glutamate dehydrogenase: 2-oxoglutarate/NADPH",
  purF = "PRPP amidotransferase: purine nucleotide feedback",
  carA = "carbamoyl-phosphate synthetase small subunit (heteromer with carB)",
  carB = "carbamoyl-phosphate synthetase large subunit: UMP/ornithine regulated",
  aceE = "pyruvate dehydrogenase E1: NADH/acetyl-CoA inhibition (complex with aceF, lpd)",
  aceF = "pyruvate dehydrogenase E2 (complex)",
  lpd  = "dihydrolipoyl dehydrogenase (PDH/OGDH complexes)",
  sucA = "2-oxoglutarate dehydrogenase E1 (complex with sucB, lpd)",
  mdh  = "malate dehydrogenase",
  gapA = "glyceraldehyde-3-P dehydrogenase",
  eno  = "enolase",
  adk  = "adenylate kinase: adenine nucleotide pool",
  nadK = "NAD kinase",
  ackA = "acetate kinase",
  pta  = "phosphotransacetylase")

validation_candidates <- function(metabolites   = NULL,
                                  n_per_metabolite = 10L,
                                  ccf_top_n      = 50L,
                                  rank_source    = c("emd", "ccf"),
                                  require_photometric = TRUE,
                                  exclude_nucleic_acid = TRUE,   # DNA/RNA binders -> their own list
                                  nucleic_acid_regex   = NULL,
                                  pBHadj_cut     = 0.05, log2fc_cut = 1,
                                  w_stat = 3, w_shift = 3, w_both = 2,
                                  w_central = 3, w_assay = 1, w_allosteric = 2,
                                  complex_portal_file = NULL,
                                  out_subdir     = "validation_candidates") {
  rank_source <- match.arg(rank_source)
  U <- .load_uniprot(); if (is.null(U)) stop("No output/uniprot_annotation_shared.RData - render a comparison first.")
  cp <- .complex_portal_members(complex_portal_file)

  # Complex Portal: protein -> the complexes it belongs to and their other subunits
  cp_tab <- NULL
  cpf <- complex_portal_file
  if (is.null(cpf)) { cand <- list.files(here("data", "raw"), pattern = "^Complex_portal.*\\.tsv$", full.names = TRUE)
                      if (length(cand)) cpf <- cand[1] }
  else if (!file.exists(cpf)) { f2 <- here("data", "raw", cpf); if (file.exists(f2)) cpf <- f2 else cpf <- NULL }
  if (!is.null(cpf) && file.exists(cpf)) {
    x <- tryCatch(as.data.table(read.csv(cpf, sep = "\t", header = TRUE, check.names = FALSE)), error = function(e) NULL)
    if (!is.null(x)) {
      idc <- grep("identifier", names(x), ignore.case = TRUE, value = TRUE)
      col <- if (length(idc)) grep("molecul|stoichiom", idc, ignore.case = TRUE, value = TRUE)[1] else NA_character_
      if (is.na(col) && length(idc)) col <- idc[1]
      nmc <- grep("recommended name", names(x), ignore.case = TRUE, value = TRUE)[1]
      if (!is.na(col)) {
        cp_tab <- rbindlist(lapply(seq_len(nrow(x)), function(i) {
          acc <- unique(toupper(unlist(regmatches(x[[col]][i], gregexpr("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", x[[col]][i])))))
          if (!length(acc)) return(NULL)
          data.table(protein_id = acc, complex_name = if (!is.na(nmc)) as.character(x[[nmc]][i]) else NA_character_,
                     n_subunits = length(acc), subunits = paste(acc, collapse = ","))
        }), use.names = TRUE)
      }
    }
  }

  if (is.null(metabolites)) {
    dirs        <- basename(list.dirs(here("output"), recursive = FALSE))
    metabolites <- sub("^PCM_ctrl_vs_", "", dirs[grepl("^PCM_ctrl_vs_", dirs)])
  }
  if (!length(metabolites)) stop("No PCM_ctrl_vs_* output folders found.")

  out <- here("output", out_subdir); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  all_cand <- list(); all_nucleic <- list(); md <- c("# Orthogonal validation shortlist", "",
    "Annotation-driven candidates for photometric activity assays of allosteric regulation.",
    "Evidence = CCprofiler differential abundance (pBHadj < 0.05, |log2FC| > 1) and/or the top-ranked",
    paste0("proteins of the ", toupper(rank_source), " elution-shift screen (exploratory: no FDR-significant hits)."),
    "", "**Verify every suggested assay against BRENDA / EcoCyc before use - the EC-to-assay mapping is a heuristic.**", "")

  for (m in metabolites) {
    cmp <- paste0("PCM_ctrl_vs_", m)
    # ---- evidence 1: CCprofiler differential ----
    pdp <- NULL
    f1 <- here("output", cmp, "rdata", "protein_DiffExprProtein_list.RData")
    if (file.exists(f1)) { e <- new.env(); load(f1, envir = e); if ("protein_DiffExprProtein" %in% ls(e)) pdp <- as.data.table(e$protein_DiffExprProtein) }
    if (is.null(pdp)) {
      fdir <- here("output", cmp, "RData_for_further_plotting_and_analysis")
      f2 <- list.files(fdir, pattern = "_for_plotting\\.RData$", full.names = TRUE)
      if (length(f2)) { e <- new.env(); load(f2[1], envir = e); if ("protein_DiffExprProtein" %in% ls(e)) pdp <- as.data.table(e$protein_DiffExprProtein) }
    }
    STAT <- if (!is.null(pdp) && all(c("feature_id", "pBHadj", "medianLog2FC") %in% names(pdp)))
      unique(pdp[, .(protein_id = as.character(feature_id), stat_p = pBHadj, log2fc = medianLog2FC)]) else NULL
    if (is.null(STAT)) message("[", m, "] no CCprofiler differential table found.")

    # ---- evidence 2: the ranked shift screen ----
    tabd <- here("output", cmp, "tables")
    rf <- file.path(tabd, if (rank_source == "emd") "emd_fdr_results.txt" else "ccf_fdr_results.txt")
    SHIFT <- NULL
    if (file.exists(rf)) {
      R <- fread(rf)
      scol <- if ("emd" %in% names(R)) "emd" else if ("abs_lag" %in% names(R)) "abs_lag" else NA_character_
      if (all(c("protein_id", "qval") %in% names(R))) {
        R <- R[is.finite(qval)]
        ord <- if (!is.na(scol)) order(R$qval, -R[[scol]]) else order(R$qval)
        R <- R[ord][seq_len(min(ccf_top_n, nrow(R)))]
        SHIFT <- data.table(protein_id = as.character(R$protein_id), shift_rank = seq_len(nrow(R)),
                            shift_q = R$qval, shift_size = if (!is.na(scol)) R[[scol]] else NA_real_)
      }
    } else message("[", m, "] no ", basename(rf), " - shift evidence unavailable (run the screen first).")

    ids <- unique(c(if (!is.null(STAT)) STAT[stat_p < pBHadj_cut & abs(log2fc) > log2fc_cut]$protein_id,
                    if (!is.null(SHIFT)) SHIFT$protein_id))
    if (!length(ids)) { message("[", m, "] no candidates from either evidence line; skipping."); next }

    D <- data.table(protein_id = ids)
    if (!is.null(STAT))  D <- merge(D, STAT,  by = "protein_id", all.x = TRUE)
    if (!is.null(SHIFT)) D <- merge(D, SHIFT, by = "protein_id", all.x = TRUE)
    D[, is_stat_hit  := !is.na(stat_p) & stat_p < pBHadj_cut & abs(log2fc) > log2fc_cut]
    D[, is_shift_hit := !is.na(shift_rank)]

    # ---- annotate ----
    ucols <- intersect(c("accession", "protein_name", "gene_names", "go_f", "go_p", "cc_catalytic_activity",
                         "ft_binding", "cc_cofactor"), names(U))
    D <- merge(D, U[, ..ucols], by.x = "protein_id", by.y = "accession", all.x = TRUE)
    if (!"gene_names" %in% names(D))    D[, gene_names := NA_character_]
    if (!"protein_name" %in% names(D))  D[, protein_name := NA_character_]
    D[, gene := toupper(trimws(sub(" .*$", "", as.character(gene_names))))]
    D[, gene_lc := tolower(gene)]
    # column-safe accessor: any annotation column absent from the cache becomes ""
    .col <- function(nm) if (nm %in% names(D)) { v <- as.character(D[[nm]]); v[is.na(v)] <- ""; v } else rep("", nrow(D))
    D[, annot_text := paste(.col("protein_name"), .col("go_f"), .col("go_p"),
                            .col("cc_catalytic_activity"), .col("ft_binding"), sep = " ; ")]
    D[, ec := sub("^EC ", "", .extract_first(.col("protein_name"), "EC [0-9]+\\.[0-9-]+\\.[0-9-]+\\.[0-9-]+"))]
    # EC may also appear only in the catalytic-activity text
    D[, ec_alt := .extract_first(.col("cc_catalytic_activity"), "[0-9]+\\.[0-9-]+\\.[0-9-]+\\.[0-9-]+")]
    D[is.na(ec), ec := ec_alt]
    D[, ec_alt := NULL]
    aa <- lapply(seq_len(nrow(D)), function(i) .assay_for(D$ec[i], D$annot_text[i]))
    D[, suggested_assay := vapply(aa, function(x) x$assay, character(1))]
    D[, assay_weight    := vapply(aa, function(x) x$w,     numeric(1))]
    D[, is_central := grepl(.CENTRAL_REGEX, annot_text, ignore.case = TRUE)]
    # DNA/RNA binders: flagged and set aside (see .NUCLEIC_REGEX for why), with the matched terms kept
    .nrx <- if (is.null(nucleic_acid_regex)) .NUCLEIC_REGEX else nucleic_acid_regex
    D[, nucleic_acid_binder := grepl(.nrx, annot_text, ignore.case = TRUE)]
    D[, nucleic_terms := vapply(annot_text, function(s) {
      tt <- trimws(unlist(strsplit(s, ";")))
      paste(unique(tt[grepl(.nrx, tt, ignore.case = TRUE)]), collapse = " | ") }, character(1), USE.NAMES = FALSE)]
    prior_names <- tolower(names(.ALLOSTERIC_PRIOR))
    D[, allosteric_note := ifelse(gene_lc %in% prior_names, unname(.ALLOSTERIC_PRIOR[match(gene_lc, prior_names)]), NA_character_)]
    D[, is_allosteric_prior := !is.na(allosteric_note)]
    # heteromer requirement
    if (!is.null(cp_tab)) {
      cpm <- cp_tab[protein_id %in% D$protein_id][order(-n_subunits)]
      cpm <- cpm[!duplicated(protein_id)]
      D <- merge(D, cpm, by = "protein_id", all.x = TRUE)
      D[, needs_partners := !is.na(n_subunits) & n_subunits > 1]
      D[, partners := ifelse(is.na(subunits), NA_character_,
                             vapply(seq_len(.N), function(i) paste(setdiff(strsplit(subunits[i], ",")[[1]], protein_id[i]), collapse = ", "), character(1)))]
    } else { D[, `:=`(complex_name = NA_character_, n_subunits = NA_integer_, needs_partners = NA, partners = NA_character_)] }

    # ---- score ----
    D[, score :=
        w_stat       * as.numeric(is_stat_hit) +
        w_shift      * as.numeric(is_shift_hit) +
        w_both       * as.numeric(is_stat_hit & is_shift_hit) +
        w_central    * as.numeric(is_central) +
        w_allosteric * as.numeric(is_allosteric_prior) +
        w_assay      * assay_weight +
        ifelse(is.na(log2fc), 0, pmin(abs(log2fc), 4) / 2)]
    D[, evidence := data.table::fcase(
      is_stat_hit & is_shift_hit, "abundance + elution shift",
      is_stat_hit,                "abundance change (CCprofiler)",
      is_shift_hit,               paste0("elution shift (", toupper(rank_source), " screen)"),
      default =                   "none")]
    if (require_photometric) D <- D[assay_weight > 0]
    setorder(D, -score)
    fwrite(D, file.path(out, paste0("validation_candidates_", m, ".txt")), sep = "\t")   # everything, flag included

    # split off the DNA/RNA binders - kept as their own list, not thrown away
    D_na <- D[nucleic_acid_binder == TRUE]
    if (exclude_nucleic_acid) D <- D[nucleic_acid_binder == FALSE]
    if (nrow(D_na)) {
      fwrite(D_na, file.path(out, paste0("nucleic_acid_binders_", m, ".txt")), sep = "\t")
      all_nucleic[[m]] <- copy(D_na)[, metabolite := m]
    }
    all_cand[[m]] <- copy(D)[, metabolite := m]
    if (!nrow(D)) { message("[", m, "] no candidates left after the DNA/RNA-binder split."); next }

    top <- head(D, n_per_metabolite)
    message("\n[", m, "] top ", nrow(top), " validation candidate(s)",
            if (exclude_nucleic_acid) paste0("  [", nrow(D_na), " DNA/RNA binder(s) moved to their own list]") else "", ":")
    print(top[, .(gene, protein_id, evidence, is_central, is_allosteric_prior, ec, score)])

    md <- c(md, paste0("## ", m, "  (", nrow(D), " scored candidates",
                       if (exclude_nucleic_acid && nrow(D_na)) paste0("; ", nrow(D_na), " DNA/RNA binders set aside") else "", ")"), "")
    for (i in seq_len(nrow(top))) {
      r <- top[i]
      md <- c(md,
        paste0("**", i, ". ", ifelse(is.na(r$gene) || !nzchar(r$gene), r$protein_id, r$gene), "** (", r$protein_id, ") - score ", round(r$score, 2)),
        paste0("- Protein: ", r$protein_name),
        paste0("- Evidence: ", r$evidence,
               ifelse(is.na(r$log2fc), "", sprintf("; log2FC = %.2f, BH p = %.3g", r$log2fc, r$stat_p)),
               ifelse(is.na(r$shift_rank), "", sprintf("; shift-screen rank %d (q = %.3g)", r$shift_rank, r$shift_q))),
        paste0("- Central metabolism: ", ifelse(isTRUE(r$is_central), "yes", "not by annotation"),
               ifelse(is.na(r$ec), "", paste0("; EC ", r$ec))),
        paste0("- Suggested photometric assay: ", ifelse(is.na(r$suggested_assay), "none obvious", r$suggested_assay)),
        if (isTRUE(r$needs_partners)) paste0("- HETEROMER: ", r$complex_name, " - co-purify/add: ", r$partners) else "- Assayable as a single purified protein (no curated obligate partners)",
        if (!is.na(r$allosteric_note)) paste0("- Known allosteric regulation (textbook prior, verify): ", r$allosteric_note) else NULL,
        "")
    }
  }

  if (length(all_cand)) {
    AC <- rbindlist(all_cand, use.names = TRUE, fill = TRUE)
    fwrite(AC, file.path(out, "validation_candidates_all.csv"))
    writeLines(md, file.path(out, "validation_shortlist.md"))
    message("\nShortlist written to ", file.path(out, "validation_shortlist.md"))

    # ---- the DNA/RNA-binder list: separate table + its own readable summary ----
    if (length(all_nucleic)) {
      AN <- rbindlist(all_nucleic, use.names = TRUE, fill = TRUE)
      fwrite(AN, file.path(out, "nucleic_acid_binders_all.csv"))
      nmd <- c("# Nucleic-acid binders with a metabolite-dependent elution change", "",
        "Set aside from the activity-assay shortlist: these proteins travel with DNA, RNA or the ribosome,",
        "so an elution change may reflect altered nucleic-acid association or ribosome engagement rather",
        "than an allosteric effect on the protein itself, and a photometric assay on a purified subunit is",
        "rarely informative for them.",
        "",
        "They are listed here because the observation is interesting in its own right - metabolite-dependent",
        "nucleoprotein remodelling, ppGpp-like signalling, or ribosome-association changes. Validating these",
        "needs a different orthogonal approach (EMSA, RNA/DNA pull-down, polysome profiling, ITC with the",
        "nucleic-acid ligand), not a cuvette activity assay.",
        "",
        "Classification is a KEYWORD HEURISTIC over GO / UniProt annotation - the matched terms are given",
        "per protein so every call can be checked.", "")
      for (m in names(all_nucleic)) {
        dn <- all_nucleic[[m]]; setorder(dn, -score)
        nmd <- c(nmd, paste0("## ", m, "  (", nrow(dn), " nucleic-acid binders)"), "")
        for (i in seq_len(min(nrow(dn), n_per_metabolite))) {
          r <- dn[i]
          nmd <- c(nmd,
            paste0("**", i, ". ", ifelse(is.na(r$gene) || !nzchar(r$gene), r$protein_id, r$gene), "** (", r$protein_id, ")"),
            paste0("- Protein: ", r$protein_name),
            paste0("- Evidence: ", r$evidence,
                   ifelse(is.na(r$log2fc), "", sprintf("; log2FC = %.2f, BH p = %.3g", r$log2fc, r$stat_p)),
                   ifelse(is.na(r$shift_rank), "", sprintf("; shift-screen rank %d", r$shift_rank))),
            paste0("- Matched nucleic-acid terms: ", substr(ifelse(is.na(r$nucleic_terms) | !nzchar(r$nucleic_terms), "(name/GO match)", r$nucleic_terms), 1, 200)),
            "")
        }
      }
      writeLines(nmd, file.path(out, "nucleic_acid_binders.md"))
      message("DNA/RNA binders (", nrow(AN), " rows across ", length(all_nucleic),
              " metabolite(s)) -> ", file.path(out, "nucleic_acid_binders.md"))
    }

    TOP <- AC[, head(.SD[order(-score)], n_per_metabolite), by = metabolite]
    TOP[, label := ifelse(is.na(gene) | !nzchar(gene), protein_id, gene)]
    g <- ggplot(TOP, aes(stats::reorder(label, score), score, fill = evidence)) +
      geom_col() + coord_flip() +
      facet_wrap(~ metabolite, scales = "free_y") +
      labs(title = "Orthogonal validation candidates", x = NULL, y = "priority score", fill = NULL,
           subtitle = "Score combines evidence (abundance and/or elution shift), central-metabolism annotation,\nphotometric assayability and a curated allosteric prior. Annotation-driven: verify before use.") +
      theme_bw() + theme(legend.position = "top")
    .fo <- file.path(out, "validation_shortlist.pdf")
    tryCatch(ggsave(.fo, g, width = 11, height = 8),
             error = function(e) message("   !! could not write ", basename(.fo), ": ", conditionMessage(e)))
    invisible(AC)
  } else { message("No candidates produced."); invisible(NULL) }
}

# ---------------------------------------------------------------------------------------------------
# Plot the SEC traces of the candidates for manual inspection.
#
# Reads output/validation_candidates/validation_candidates_all.csv (written by validation_candidates())
# and draws plot_protein_traces() for each candidate, ONE PROTEIN PER PDF PAGE, per metabolite. A protein
# that is a candidate for several metabolites gets one page in each of those metabolites' PDFs, which is
# the point - the trace is only meaningful against the condition it was called in.
#
# NOTE ON DUPLICATES: the candidate table has one row per protein PER METABOLITE, so a protein hit in
# several comparisons appears several times. Ids are de-duplicated WITHIN each metabolite here.
#
#   plot_validation_candidates()                          # all metabolites, stat hits only
#   plot_validation_candidates("ATP")
#   plot_validation_candidates(only_stat_hits = FALSE)    # include shift-screen-only candidates
#   plot_validation_candidates(ids = c("P15723","P0A6F5"))# an explicit list you pasted yourself
#   plot_validation_candidates(set = "nucleic_acid")      # the DNA/RNA-binder list instead
#   plot_validation_candidates(top_n = 20)                # only the 20 best-scoring per metabolite
#
# OUTPUT (per metabolite): output/PCM_ctrl_vs_<m>/validation_traces/
#   <m>_candidate_traces.pdf   one page per protein
#   <m>_candidate_ids.txt      the plotted accessions, one per line (handy for copy/paste elsewhere)
# ---------------------------------------------------------------------------------------------------
plot_validation_candidates <- function(metabolites   = NULL,
                                       ids            = NULL,
                                       set            = c("candidates", "nucleic_acid"),
                                       only_stat_hits = TRUE,
                                       top_n          = Inf,
                                       max_per_metabolite = 250,
                                       x_axis         = c("fraction", "mw"),
                                       aggregate      = c("condition", "replicate"),
                                       out_subdir     = "validation_traces") {
  x_axis <- match.arg(x_axis); aggregate <- match.arg(aggregate); set <- match.arg(set)
  source(here::here("scripts", "plot_protein_traces.R"))
  if (set == "nucleic_acid" && identical(out_subdir, "validation_traces")) out_subdir <- "nucleic_acid_traces"

  AC <- NULL
  f <- here("output", "validation_candidates",
            if (set == "nucleic_acid") "nucleic_acid_binders_all.csv" else "validation_candidates_all.csv")
  if (is.null(ids)) {
    if (!file.exists(f)) stop("No ", f, " - run validation_candidates() first, or pass ids = c(...).")
    AC <- fread(f)
    if (only_stat_hits && "is_stat_hit" %in% names(AC)) AC <- AC[as.logical(is_stat_hit) %in% TRUE]
    if (!nrow(AC)) stop("No candidates left after filtering (try only_stat_hits = FALSE).")
  }
  if (is.null(metabolites)) {
    metabolites <- if (!is.null(AC)) unique(as.character(AC$metabolite)) else {
      dirs <- basename(list.dirs(here("output"), recursive = FALSE))
      sub("^PCM_ctrl_vs_", "", dirs[grepl("^PCM_ctrl_vs_", dirs)]) }
  }

  for (m in metabolites) {
    if (!is.null(ids)) {
      pid <- unique(as.character(ids))
    } else {
      sub <- AC[metabolite == m]
      if ("score" %in% names(sub)) setorder(sub, -score)
      pid <- unique(as.character(sub$protein_id))          # de-duplicate within the metabolite
      if (is.finite(top_n)) pid <- head(pid, top_n)
    }
    pid <- pid[!is.na(pid) & nzchar(pid)]
    if (!length(pid)) { message("[", m, "] no candidates to plot."); next }
    if (length(pid) > max_per_metabolite) {
      message("[", m, "] ", length(pid), " candidates - plotting the first ", max_per_metabolite,
              " (raise max_per_metabolite, or use top_n).")
      pid <- head(pid, max_per_metabolite)
    }

    outdir <- here("output", paste0("PCM_ctrl_vs_", m), out_subdir)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    writeLines(pid, file.path(outdir, paste0(m, "_candidate_ids.txt")))

    pdf_file <- file.path(outdir, paste0(m, "_candidate_traces.pdf"))
    ok <- tryCatch({ grDevices::pdf(pdf_file, width = 6, height = 6); TRUE },
                   error = function(e) { message("   !! cannot write ", basename(pdf_file), ": ",
                                                 conditionMessage(e), " (open in a viewer?)"); FALSE })
    if (!ok) next
    np <- 0L; miss <- character(0)
    for (p in pid) {
      done <- tryCatch({
        g <- plot_protein_traces(p, metabolites = m, x_axis = x_axis, aggregate = aggregate,
                                 save_pdf = FALSE, print_plot = FALSE)
        if (!is.null(g)) print(g)
        !is.null(g)
      }, error = function(e) { miss <<- c(miss, p); FALSE })
      if (isTRUE(done)) np <- np + 1L
    }
    grDevices::dev.off()
    message("[", m, "] ", np, "/", length(pid), " protein page(s) -> ", pdf_file)
    if (length(miss)) message("   not plottable (absent from the traces): ", paste(utils::head(miss, 10), collapse = ", "),
                              if (length(miss) > 10) paste0(" ... +", length(miss) - 10, " more") else "")
  }
  invisible(NULL)
}
