# scripts/investigate_fraction1_effect.R
# =============================================================================
# INVESTIGATION ONLY - does NOT change the pipeline or any result on disk.
#
# Question: does the fraction-1 (void-volume) intensity spike that cyclic-loess normalization leaves
# behind actually bias the ASSEMBLY state? Fraction 1 is the highest apparent MW, so any inflated
# signal there is counted as "assembled". This script recomputes the monomer/assembled split
# WITH vs WITHOUT fraction 1 and reports how much it moves - so we decide from evidence whether the
# spike matters, before touching Benni's normalization.
#
# HOW: for each metabolite it loads the saved traces (output/PCM_ctrl_vs_<m>/RData_for_further_.../
# *_for_plotting.RData -> protein_traces_list, design_matrix), runs CCprofiler's own
# annotateMassDistribution() twice (original, and with the excluded fraction(s) zeroed), and compares
# the per-protein monomer fraction. With full = TRUE it additionally recomputes the assembly hit sets
# (getMassAssemblyChange, betareg - slow) and compares the more/less-assembled protein lists.
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "investigate_fraction1_effect.R"))
#   investigate_fraction1_effect()                     # all run metabolites, fast monomer-fraction check
#   investigate_fraction1_effect("PEP")                # one metabolite
#   investigate_fraction1_effect("PEP", exclude = 1:2) # exclude fractions 1 AND 2
#   investigate_fraction1_effect(full = TRUE)          # also recompute assembly HITS and compare (slow)
#
# OUTPUTS (output/fraction1_effect/):
#   <m>_monomer_shift.csv     per-protein monomer fraction with vs without the excluded fraction(s)
#   <m>_monomer_shift.pdf     scatter (with vs without) + histogram of the shift
#   fraction1_effect_summary.csv  one row per metabolite: median/max shift, # proteins moved
#   <m>_assembly_hit_diff.csv (full = TRUE only) hits gained/lost when the fraction is excluded
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(data.table); library(ggplot2); library(CCprofiler)
})

# --- helpers -----------------------------------------------------------------------------------
# integer-named fraction columns of a traces object's $traces
.frac_cols <- function(dt) {
  fc <- grep("^[0-9]+$", colnames(dt), value = TRUE)
  fc[order(as.numeric(fc))]
}

# return a fresh tracesList copy with the given fraction column(s) set to 0
.zero_fractions <- function(tl, fractions) {
  out <- tl
  for (s in names(out)) {
    dt   <- data.table::copy(as.data.table(out[[s]]$traces))
    cols <- intersect(as.character(fractions), names(dt))
    if (length(cols)) dt[, (cols) := 0]
    out[[s]]$traces <- dt
  }
  out
}

# per-protein monomer & assembled sums from an annotateMassDistribution() result (defensive on column names)
.extract_assembly <- function(assembly_tl) {
  pick <- function(v) if (!length(v)) NA_character_
                      else if (any(grepl("norm", v, ignore.case = TRUE))) grep("norm", v, ignore.case = TRUE, value = TRUE)[1]
                      else v[1]
  rbindlist(lapply(names(assembly_tl), function(s) {
    ta <- as.data.table(assembly_tl[[s]]$trace_annotation)
    mc <- pick(grep("monomer", names(ta), ignore.case = TRUE, value = TRUE))
    ac <- pick(grep("assembl", names(ta), ignore.case = TRUE, value = TRUE))
    if (is.na(mc) || is.na(ac)) return(NULL)
    idcol <- if ("protein_id" %in% names(ta)) "protein_id" else "id"
    data.table(protein   = as.character(ta[[idcol]]),
               monomer   = suppressWarnings(as.numeric(ta[[mc]])),
               assembled = suppressWarnings(as.numeric(ta[[ac]])))
  }), use.names = TRUE)
}

# per-protein monomer fraction (pooled over samples)
.monomer_fraction <- function(assembly_tl) {
  d <- .extract_assembly(assembly_tl)
  if (is.null(d) || !nrow(d)) return(NULL)
  d <- d[, .(monomer = sum(monomer, na.rm = TRUE), assembled = sum(assembled, na.rm = TRUE)), by = protein]
  d[, mono_frac := ifelse(monomer + assembled > 0, monomer / (monomer + assembled), NA_real_)]
  d[, .(protein, mono_frac)]
}

# --- main --------------------------------------------------------------------------------------
investigate_fraction1_effect <- function(metabolites = NULL,
                                         exclude     = 1,
                                         full        = FALSE,
                                         out_subdir  = "fraction1_effect") {
  out <- here("output", out_subdir); dir.create(out, recursive = TRUE, showWarnings = FALSE)

  if (is.null(metabolites)) {
    dirs        <- basename(list.dirs(here("output"), recursive = FALSE))
    metabolites <- sub("^PCM_ctrl_vs_", "", dirs[grepl("^PCM_ctrl_vs_", dirs)])
  }
  if (!length(metabolites)) stop("No PCM_ctrl_vs_* output folders found - render at least one metabolite first.")

  if (full) {
    # getMassAssemblyChange_aljazfix + the assembly cutoffs live in the fixes file
    fixes <- here("R", "ccprofiler_fixes.R")
    if (file.exists(fixes)) suppressWarnings(try(source(fixes), silent = TRUE))
    if (!exists("getMassAssemblyChange_aljazfix")) {
      message("full = TRUE requested but getMassAssemblyChange_aljazfix not available - doing the fast check only.")
      full <- FALSE
    }
  }

  summary_rows <- list()
  for (m in metabolites) {
    fdir <- here("output", paste0("PCM_ctrl_vs_", m), "RData_for_further_plotting_and_analysis")
    f    <- list.files(fdir, pattern = "_for_plotting\\.RData$", full.names = TRUE)
    if (!length(f)) { message("[", m, "] no *_for_plotting.RData - render this metabolite first; skipping."); next }

    # two INDEPENDENT loads so zeroing one can't touch the other (data.tables are by-reference)
    e1 <- new.env(); load(f[1], envir = e1)
    e2 <- new.env(); load(f[1], envir = e2)
    if (!all(c("protein_traces_list") %in% ls(e1))) {
      message("[", m, "] no protein_traces_list in the file; skipping."); next
    }
    # sanity: is the fraction we want to exclude actually present?
    present <- intersect(as.character(exclude), .frac_cols(as.data.table(e1$protein_traces_list[[1]]$traces)))
    if (!length(present)) { message("[", m, "] fraction(s) ", paste(exclude, collapse = ","), " not in the traces; skipping."); next }

    a_with    <- tryCatch(annotateMassDistribution(e1$protein_traces_list), error = function(err) NULL)
    a_without <- tryCatch(annotateMassDistribution(.zero_fractions(e2$protein_traces_list, exclude)), error = function(err) NULL)
    if (is.null(a_with) || is.null(a_without)) {
      message("[", m, "] annotateMassDistribution failed; skipping."); next
    }

    mf1 <- .monomer_fraction(a_with);    if (!is.null(mf1)) setnames(mf1, "mono_frac", "mono_frac_with")
    mf0 <- .monomer_fraction(a_without); if (!is.null(mf0)) setnames(mf0, "mono_frac", "mono_frac_without")
    if (is.null(mf1) || is.null(mf0)) { message("[", m, "] could not read monomer/assembled columns; skipping."); next }

    cmp <- merge(mf1, mf0, by = "protein")
    cmp[, shift := mono_frac_without - mono_frac_with]   # >0 => excluding the fraction ADDS monomer (i.e. it was inflating "assembled")
    fwrite(cmp[order(-abs(shift))], file.path(out, paste0(m, "_monomer_shift.csv")))

    n        <- nrow(cmp)
    med_abs  <- median(abs(cmp$shift), na.rm = TRUE)
    max_abs  <- max(abs(cmp$shift), na.rm = TRUE)
    p_gt05   <- mean(abs(cmp$shift) > 0.05, na.rm = TRUE)
    p_gt10   <- mean(abs(cmp$shift) > 0.10, na.rm = TRUE)
    summary_rows[[m]] <- data.table(metabolite = m, n_proteins = n,
                                    median_abs_shift = round(med_abs, 4), max_abs_shift = round(max_abs, 4),
                                    pct_moved_gt5 = round(100 * p_gt05, 1), pct_moved_gt10 = round(100 * p_gt10, 1))
    message(sprintf("[%s] %d proteins | median |Δmonomer| = %.1f%% | %.1f%% move >5%% | %.1f%% move >10%% (max %.1f%%)",
                    m, n, 100 * med_abs, 100 * p_gt05, 100 * p_gt10, 100 * max_abs))

    # plots: scatter (with vs without) + histogram of the shift
    g1 <- ggplot(cmp, aes(mono_frac_with, mono_frac_without)) +
      geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = 2) +
      geom_point(alpha = 0.4, size = 0.8) +
      coord_equal() +
      labs(title = paste0(m, ": monomer fraction, with vs without fraction ", paste(exclude, collapse = ",")),
           x = "monomer fraction (with)", y = "monomer fraction (without)") + theme_bw()
    g2 <- ggplot(cmp, aes(shift)) +
      geom_histogram(bins = 60, fill = "steelblue", colour = "white") +
      geom_vline(xintercept = 0, colour = "grey40") +
      labs(title = "shift in monomer fraction (without - with)", x = "Δ monomer fraction", y = "proteins") + theme_bw()
    grDevices::pdf(file.path(out, paste0(m, "_monomer_shift.pdf")), width = 6, height = 8)
    print(g1); print(g2); grDevices::dev.off()

    # OPTIONAL: recompute the actual assembly hit sets (slow) and compare
    if (full && "design_matrix" %in% ls(e1)) {
      hits <- function(assembly_tl, dm) {
        das <- tryCatch(getMassAssemblyChange_aljazfix(tracesList = assembly_tl, design_matrix = dm,
                          compare_between = "Condition", quantLevel = "protein_id", plot = FALSE, PDF = FALSE),
                        error = function(err) { message("   getMassAssemblyChange failed: ", conditionMessage(err)); NULL })
        if (is.null(das)) return(NULL)
        das <- as.data.table(das)
        if (!all(c("betaPval_BHadj", "meanDiff", "protein_id") %in% names(das))) return(das)  # return raw if shape differs
        das[betaPval_BHadj < 0.05 & abs(meanDiff) > 0.1]$protein_id
      }
      h_with    <- hits(a_with, e1$design_matrix)
      h_without <- hits(a_without, e2$design_matrix)
      if (!is.null(h_with) && !is.null(h_without) && is.character(h_with) && is.character(h_without)) {
        gained <- setdiff(h_without, h_with); lost <- setdiff(h_with, h_without)
        fwrite(rbind(data.table(change = "lost_when_excluded",   protein = lost),
                     data.table(change = "gained_when_excluded", protein = gained)),
               file.path(out, paste0(m, "_assembly_hit_diff.csv")))
        message(sprintf("   [%s] assembly hits: %d with, %d without | %d lost, %d gained when fraction excluded",
                        m, length(h_with), length(h_without), length(lost), length(gained)))
      }
    }
  }

  if (length(summary_rows)) {
    S <- rbindlist(summary_rows, use.names = TRUE)
    fwrite(S, file.path(out, "fraction1_effect_summary.csv"))
    message("\n==== fraction-", paste(exclude, collapse = ","), " effect on the assembly state ====")
    print(S)
    message("Interpretation: median_abs_shift is the typical change in a protein's monomer fraction when the ",
            "fraction is excluded. If it is small (a few %) and few proteins move, the void spike is cosmetic ",
            "and does not drive the assembly results. Written to: ", file.path(out, "fraction1_effect_summary.csv"))
    invisible(S)
  } else {
    message("Nothing computed - no metabolite had usable saved traces."); invisible(NULL)
  }
}
