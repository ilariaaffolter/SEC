# scripts/hydropro_ffo.R
# =============================================================================
# THEORETICAL frictional ratio (f/f0) per protein, to break the SEC degeneracy.
#
# WHY: SEC measures hydrodynamic size, not mass. What you recover is
#         apparent_MW / monomer_MW  =  n * (f/f0)^3
# with n = oligomeric state (ASSEMBLY) and f/f0 = shape factor (ELONGATION). One measurement, two
# unknowns - a disordered monomer with f/f0 = 2 elutes exactly where a compact octamer elutes.
# If f/f0 of the MONOMER is known independently, n follows:
#         n_implied = (apparent_MW / monomer_MW) / (f/f0_theoretical)^3
# This script computes that theoretical f/f0 from 3D structure with HYDROPRO, and compares it to the
# experimental ratio from scripts/globularity_check.R.
#
# PIPELINE (each step is a separate function, so you can stop/resume):
#   1. hydropro_fetch()    download AlphaFold models (AF-<ACC>-F1-model_v4.pdb) -> output/hydropro/structures/
#                          and compute the pLDDT-based disorder metrics from the PDB B-factor column
#                          (AlphaFold stores per-residue pLDDT there) -> also used by the annotation script
#   2. hydropro_prepare()  write one hydropro.dat per protein into its own run folder
#   3. hydropro_run()      run the HYDROPRO executable over those folders  [EXTERNAL SOFTWARE - see below]
#   4. hydropro_parse()    read each *-res.txt -> Stokes radius -> f/f0_theoretical
#   5. hydropro_compare()  merge with globularity_check.R and split ASSEMBLY vs ELONGATION
#   or  hydropro_all()     run 1-5 in order
#
# EXTERNAL SOFTWARE: HYDROPRO is not bundled and cannot be installed from here. Download it (free for
# academic use, Garcia de la Torre lab, https://leonardo.inf.um.es/macromol/programs/hydropro/hydropro.htm),
# unzip it, and pass the executable path:  hydropro_run(exe = "C:/hydropro/hydropro10-msd.exe")
# Steps 1, 2, 4, 5 work WITHOUT it; step 4 will simply find no results until you have run step 3.
#
# HONEST CAVEATS - read before interpreting:
#   * AlphaFold models intrinsically disordered regions as arbitrary extended ribbons. HYDROPRO on such a
#     model returns an INFLATED f/f0 that reflects the model's arbitrary IDR placement, not a real
#     conformational ensemble. Every output carries the pLDDT disorder metrics so these can be filtered
#     or flagged (see `plddt_disorder_frac`); treat high-disorder proteins as unreliable, not as evidence.
#   * The AlphaFold model is a MONOMER prediction; f/f0 is therefore the monomer shape factor, which is
#     exactly what is needed here, but it says nothing about the assembled form's shape.
#   * n_implied is continuous - it will not land exactly on an integer. Round with judgement, and prefer
#     co-elution with known partners (the Complex Portal analysis) as independent evidence of assembly.
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "hydropro_ffo.R"))
#   hydropro_fetch("ATP")                                   # structures + pLDDT for that metabolite's proteins
#   hydropro_fetch("ATP", ids = c("P0A6F5","P0A6Y8"))        # or an explicit id list
#   hydropro_prepare()
#   hydropro_run(exe = "C:/hydropro/hydropro10-msd.exe")     # slow: minutes per protein
#   hydropro_parse(); hydropro_compare("ATP")
#   hydropro_all("ATP", exe = "C:/hydropro/hydropro10-msd.exe", max_proteins = 50)
#
# OUTPUT (output/hydropro/):
#   structures/AF-<ACC>-F1-model_v4.pdb   cached AlphaFold models
#   plddt_disorder.csv                    per protein: mean pLDDT, disordered residue fractions
#   runs/<ACC>/hydropro.dat + results     HYDROPRO inputs / outputs
#   hydropro_ffo.csv                      per protein: Stokes radius, f/f0_theoretical
#   ffo_theory_vs_experiment.csv          + experimental ratio, n_implied, assembly-vs-shape call
#   ffo_theory_vs_experiment.pdf          scatter: theoretical vs experimental f/f0
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })

.hp_dir  <- function(...) here("output", "hydropro", ...)
.AF_URL  <- "https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v4.pdb"

# ---- 1. structures + pLDDT ------------------------------------------------------------------------
# AlphaFold stores per-residue pLDDT (0-100) in the PDB B-factor column. Low pLDDT is the standard
# proxy for disorder: <50 = very likely disordered, <70 = low confidence / likely flexible.
.plddt_from_pdb <- function(pdb_file) {
  ln <- tryCatch(readLines(pdb_file, warn = FALSE), error = function(e) character(0))
  ca <- grep("^ATOM.{9}CA ", ln, value = TRUE)          # one entry per residue
  if (!length(ca)) return(NULL)
  b <- suppressWarnings(as.numeric(substr(ca, 61, 66)))
  b <- b[is.finite(b)]
  if (!length(b)) return(NULL)
  list(n_res = length(b), mean_plddt = mean(b),
       frac_lt70 = mean(b < 70), frac_lt50 = mean(b < 50))
}

hydropro_fetch <- function(metabolites = NULL, ids = NULL, max_proteins = Inf, overwrite = FALSE) {
  sdir <- .hp_dir("structures"); dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(ids)) {
    if (is.null(metabolites)) stop("Give either ids = c(...) or metabolites = '<name>'.")
    f <- here("output", paste0("PCM_ctrl_vs_", metabolites[1]), "tables", "globularity_check.txt")
    if (!file.exists(f)) stop("No globularity_check.txt for ", metabolites[1], " - run globularity_check() first.")
    ids <- unique(as.character(fread(f)$protein_id))
  }
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) > max_proteins) { message("Limiting to the first ", max_proteins, " of ", length(ids), " proteins."); ids <- head(ids, max_proteins) }

  rows <- vector("list", length(ids)); nok <- 0L
  for (i in seq_along(ids)) {
    acc <- ids[i]; dest <- file.path(sdir, sprintf("AF-%s-F1-model_v4.pdb", acc))
    if (overwrite || !file.exists(dest) || file.size(dest) < 1000) {
      ok <- tryCatch({ utils::download.file(sprintf(.AF_URL, acc), dest, mode = "wb", quiet = TRUE); TRUE },
                     error = function(e) FALSE)
      if (!ok || !file.exists(dest) || file.size(dest) < 1000) {
        try(file.remove(dest), silent = TRUE)
        rows[[i]] <- data.table(protein_id = acc, has_structure = FALSE); next
      }
    }
    p <- .plddt_from_pdb(dest); nok <- nok + 1L
    rows[[i]] <- if (is.null(p)) data.table(protein_id = acc, has_structure = TRUE) else
      data.table(protein_id = acc, has_structure = TRUE, n_res = p$n_res,
                 mean_plddt = round(p$mean_plddt, 2),
                 plddt_disorder_frac = round(p$frac_lt70, 4),   # pLDDT < 70
                 plddt_verylow_frac  = round(p$frac_lt50, 4))   # pLDDT < 50
    if (i %% 50 == 0) message("  ... ", i, "/", length(ids), " structures")
  }
  D <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  dir.create(.hp_dir(), recursive = TRUE, showWarnings = FALSE)
  fwrite(D, .hp_dir("plddt_disorder.csv"))
  message("Structures available for ", nok, "/", length(ids), " proteins. pLDDT metrics -> ", .hp_dir("plddt_disorder.csv"))
  invisible(D)
}

# ---- 2. HYDROPRO input files ----------------------------------------------------------------------
# Standard HYDROPRO 10 shell-model input (one case per file; the trailing '*' closes it). Molecular
# weight and partial specific volume drive the reference sphere, so they must be per protein.
.write_hydropro_dat <- function(path, acc, pdb_basename, mw_Da, temp_C = 20, eta_poise = 0.01,
                                vbar = 0.702, solvent_density = 1.0, aer = 2.9) {
  writeLines(c(
    sprintf("%-24s !name of molecule", acc),
    sprintf("%-24s !name for output file", acc),
    sprintf("%-24s !structural (PDB) file", pdb_basename),
    "1,                       !type of calculation (1 = shell model, atomic level)",
    sprintf("%-24s !AER, radius of primary elements", paste0(aer, ",")),
    "-1,                      !NSIG",
    sprintf("%-24s !T (temperature, centigrade)", paste0(temp_C, ",")),
    sprintf("%-24s !ETA (viscosity of the solvent, poises)", paste0(eta_poise, ",")),
    sprintf("%-24s !RM (molecular weight)", paste0(round(mw_Da, 1), ",")),
    sprintf("%-24s !partial specific volume (cm3/g)", paste0(vbar, ",")),
    sprintf("%-24s !solvent density (g/cm3)", paste0(solvent_density, ",")),
    "-1                       !Number of values of Q",
    "-1                       !Number of intervals",
    "0,                       !Number of trials for MC calculation of covolume",
    "1                        !IDIF=1 (yes) for full diffusion tensors",
    "*                        !end of file"), path)
}

hydropro_prepare <- function(mw_table = NULL, temp_C = 20, eta_poise = 0.01, vbar = 0.702) {
  sdir <- .hp_dir("structures"); rdir <- .hp_dir("runs")
  pdbs <- list.files(sdir, pattern = "^AF-.*\\.pdb$", full.names = TRUE)
  if (!length(pdbs)) stop("No structures in ", sdir, " - run hydropro_fetch() first.")
  # monomer MW (Da): from the caller, else from the shared UniProt cache
  if (is.null(mw_table)) {
    sf <- here("output", "uniprot_annotation_shared.RData")
    if (!file.exists(sf)) stop("No MW source: pass mw_table = data.frame(protein_id, mw_Da) or fetch UniProt first.")
    e <- new.env(); load(sf, envir = e)
    u <- as.data.table(e$.uniprot_all)
    if (!all(c("input_id", "mass") %in% names(u))) stop("Shared UniProt cache has no input_id/mass.")
    mw_table <- u[!is.na(mass), .(protein_id = as.character(input_id), mw_Da = as.numeric(mass))]
  }
  mw_table <- as.data.table(mw_table); setkey(mw_table, protein_id)
  n <- 0L
  for (p in pdbs) {
    acc <- sub("^AF-(.*)-F1-model_v4\\.pdb$", "\\1", basename(p))
    mw  <- mw_table[acc]$mw_Da[1]
    if (!length(mw) || !is.finite(mw)) next                 # no mass -> cannot set the reference sphere
    d <- file.path(rdir, acc); dir.create(d, recursive = TRUE, showWarnings = FALSE)
    file.copy(p, file.path(d, basename(p)), overwrite = TRUE)
    .write_hydropro_dat(file.path(d, "hydropro.dat"), acc, basename(p), mw, temp_C, eta_poise, vbar)
    n <- n + 1L
  }
  message("Prepared ", n, " HYDROPRO run folder(s) under ", rdir)
  invisible(n)
}

# ---- 3. run HYDROPRO ------------------------------------------------------------------------------
hydropro_run <- function(exe, only = NULL, overwrite = FALSE, timeout_s = 1800) {
  if (missing(exe) || !nzchar(exe) || !file.exists(exe))
    stop("HYDROPRO executable not found. Download it (free for academic use) and pass exe = '<path to hydropro10-msd.exe>'.")
  rdir <- .hp_dir("runs"); dirs <- list.dirs(rdir, recursive = FALSE)
  if (!is.null(only)) dirs <- dirs[basename(dirs) %in% only]
  if (!length(dirs)) stop("No prepared run folders - run hydropro_prepare() first.")
  exe <- normalizePath(exe, winslash = "/", mustWork = TRUE)
  ok <- 0L
  for (d in dirs) {
    acc <- basename(d); res <- file.path(d, paste0(acc, "-res.txt"))
    if (!overwrite && file.exists(res)) { ok <- ok + 1L; next }
    wd <- setwd(d)
    r <- tryCatch(system2(exe, stdout = FALSE, stderr = FALSE, timeout = timeout_s),
                  error = function(e) { message("  [", acc, "] ", conditionMessage(e)); 1L })
    setwd(wd)
    if (file.exists(res)) ok <- ok + 1L else message("  [", acc, "] no -res.txt produced (exit ", r, ")")
  }
  message("HYDROPRO finished: ", ok, "/", length(dirs), " with results.")
  invisible(ok)
}

# ---- 4. parse results -> theoretical f/f0 ---------------------------------------------------------
# Prefer HYDROPRO's translational equivalent (Stokes) radius; else derive it from the translational
# diffusion coefficient via Stokes-Einstein  Rs = kT / (6 pi eta D).
.parse_hydropro_res <- function(res_file, temp_C = 20, eta_poise = 0.01) {
  ln <- tryCatch(readLines(res_file, warn = FALSE), error = function(e) character(0))
  if (!length(ln)) return(NULL)
  num1 <- function(pat) {
    i <- grep(pat, ln, ignore.case = TRUE)
    if (!length(i)) return(NA_real_)
    v <- regmatches(ln[i[1]], regexpr("[-+0-9.]+([Ee][-+]?[0-9]+)?", sub(pat, "", ln[i[1]], ignore.case = TRUE)))
    if (!length(v)) NA_real_ else as.numeric(v)
  }
  rs <- num1(".*translational.*radius|.*radius.*translational")   # cm
  dt <- num1(".*translational diffusion coefficient")             # cm2/s
  if (!is.finite(rs) && is.finite(dt)) {
    kB <- 1.380649e-16; Tk <- temp_C + 273.15                     # erg/K, cgs
    rs <- kB * Tk / (6 * pi * eta_poise * dt)                     # cm
  }
  if (!is.finite(rs)) return(NULL)
  list(stokes_radius_cm = rs, dt_cm2_s = dt)
}

hydropro_parse <- function(temp_C = 20, eta_poise = 0.01, vbar = 0.702, mw_table = NULL) {
  rdir <- .hp_dir("runs"); dirs <- list.dirs(rdir, recursive = FALSE)
  if (!length(dirs)) stop("No run folders - run hydropro_prepare()/hydropro_run() first.")
  if (is.null(mw_table)) {
    sf <- here("output", "uniprot_annotation_shared.RData")
    if (file.exists(sf)) { e <- new.env(); load(sf, envir = e); u <- as.data.table(e$.uniprot_all)
      if (all(c("input_id", "mass") %in% names(u)))
        mw_table <- u[!is.na(mass), .(protein_id = as.character(input_id), mw_Da = as.numeric(mass))] }
  }
  mw_table <- if (is.null(mw_table)) data.table(protein_id = character(), mw_Da = numeric()) else as.data.table(mw_table)
  setkey(mw_table, protein_id)

  rows <- lapply(dirs, function(d) {
    acc <- basename(d); res <- file.path(d, paste0(acc, "-res.txt"))
    if (!file.exists(res)) return(NULL)
    p <- .parse_hydropro_res(res, temp_C, eta_poise); if (is.null(p)) return(NULL)
    mw <- mw_table[acc]$mw_Da[1]
    # minimal (anhydrous sphere) radius for that mass: R_min = (3 M vbar / (4 pi N_A))^(1/3)
    rmin <- if (length(mw) && is.finite(mw)) (3 * mw * vbar / (4 * pi * 6.02214076e23))^(1/3) else NA_real_
    data.table(protein_id = acc, mw_Da = mw,
               stokes_radius_nm = p$stokes_radius_cm * 1e7,
               rmin_nm = rmin * 1e7,
               ffo_theoretical = if (is.finite(rmin)) p$stokes_radius_cm / rmin else NA_real_)
  })
  R <- rbindlist(rows[!vapply(rows, is.null, logical(1))], use.names = TRUE, fill = TRUE)
  if (!nrow(R)) { message("No parsable HYDROPRO results yet (run step 3)."); return(invisible(NULL)) }
  pf <- .hp_dir("plddt_disorder.csv")
  if (file.exists(pf)) R <- merge(R, fread(pf), by = "protein_id", all.x = TRUE)
  fwrite(R, .hp_dir("hydropro_ffo.csv"))
  message("Parsed ", nrow(R), " result(s) -> ", .hp_dir("hydropro_ffo.csv"))
  message(sprintf("  median theoretical f/f0 = %.2f", stats::median(R$ffo_theoretical, na.rm = TRUE)))
  invisible(R)
}

# ---- 5. compare with the SEC measurement ----------------------------------------------------------
hydropro_compare <- function(metabolites, plddt_disorder_max = 0.4) {
  hf <- .hp_dir("hydropro_ffo.csv")
  if (!file.exists(hf)) stop("No hydropro_ffo.csv - run hydropro_parse() first.")
  H <- fread(hf)
  for (m in metabolites) {
    gf <- here("output", paste0("PCM_ctrl_vs_", m), "tables", "globularity_check.txt")
    if (!file.exists(gf)) { message("[", m, "] no globularity_check.txt - run globularity_check() first; skipping."); next }
    G <- fread(gf)
    D <- merge(G, H, by = "protein_id")
    if (!nrow(D)) { message("[", m, "] no protein has both a HYDROPRO result and a SEC measurement."); next }

    # ratio = n * (f/f0)^3  ->  n_implied = ratio / (f/f0_theoretical)^3
    D[, n_implied := ratio / (ffo_theoretical^3)]
    D[, shape_explains := ffo_theoretical^3]        # ratio a MONOMER of this shape would already give
    D[, interpretation := data.table::fcase(
      !is.finite(n_implied),                        "undetermined",
      n_implied >= 1.5,                             "assembly (larger than shape alone explains)",
      n_implied <= 0.67,                            "elutes smaller than shape predicts",
      default =                                     "shape alone explains it (monomer)")]
    # AlphaFold renders IDRs as arbitrary extended ribbons -> theoretical f/f0 unreliable there
    if ("plddt_disorder_frac" %in% names(D))
      D[plddt_disorder_frac > plddt_disorder_max,
        interpretation := paste0(interpretation, " [UNRELIABLE: ", round(100 * plddt_disorder_frac), "% low-pLDDT]")]

    tab_dir <- here("output", paste0("PCM_ctrl_vs_", m), "tables")
    fig_dir <- here("output", paste0("PCM_ctrl_vs_", m), "figures")
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(D[order(-n_implied)], file.path(tab_dir, "ffo_theory_vs_experiment.csv"))
    message("[", m, "] ", nrow(D), " protein(s) with both measurements:")
    print(D[, .N, by = sub(" \\[.*", "", interpretation)][order(-N)])

    g <- ggplot(D, aes(ffo_theoretical, ffo_vs_monomer,
                       colour = if ("plddt_disorder_frac" %in% names(D)) plddt_disorder_frac else NULL)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
      geom_point(alpha = 0.7, size = 1.4) +
      { if ("plddt_disorder_frac" %in% names(D)) scale_colour_viridis_c(option = "C", name = "pLDDT<70\nfraction") else NULL } +
      labs(title = paste0("Shape vs assembly: theoretical vs experimental f/f0 - ", m),
           subtitle = paste0("dashed = shape alone explains the elution (monomer).\n",
                             "Above the line = elutes larger than its shape predicts -> ASSEMBLY.\n",
                             "High-disorder points (bright) are unreliable: AlphaFold models IDRs as extended ribbons."),
           x = "theoretical f/f0 (HYDROPRO on the AlphaFold monomer)",
           y = "experimental apparent f/f0 (SEC, assuming monomer)") +
      theme_bw()
    ggsave(file.path(fig_dir, "ffo_theory_vs_experiment.pdf"), g, width = 7.5, height = 6)
  }
  invisible(NULL)
}

hydropro_all <- function(metabolites, exe, max_proteins = Inf, ids = NULL) {
  hydropro_fetch(metabolites, ids = ids, max_proteins = max_proteins)
  hydropro_prepare()
  hydropro_run(exe = exe)
  hydropro_parse()
  hydropro_compare(metabolites)
}
