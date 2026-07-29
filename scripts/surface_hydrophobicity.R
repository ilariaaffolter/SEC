# scripts/surface_hydrophobicity.R
# =============================================================================
# SURFACE hydrophobicity from AlphaFold structures - the test GRAVY cannot do.
#
# WHY THIS EXISTS. globularity_category_annotation() correlates GRAVY against elution deviation and
# returns rho = +0.196. The SIGN is informative: hydrophobic retention by the column matrix would make
# sticky proteins elute LATE and appear SMALLER (a NEGATIVE correlation), so the observed positive sign
# excludes that artefact. But the positive sign itself CANNOT be read as "hydrophobic surfaces drive
# assembly", for one decisive reason: GRAVY is the mean hydropathy of the WHOLE SEQUENCE, and that mean
# is dominated by BURIED core residues. Assembly is driven by hydrophobicity that is EXPOSED. GRAVY and
# surface hydrophobicity are only loosely related - a protein can have a greasy core and a polar surface,
# which is in fact the normal state of a soluble protein.
#
# This script computes what the hypothesis is actually about, from structure:
#   sasa_total            solvent-accessible surface area (A^2), Shrake-Rupley
#   surface_gravy         SASA-WEIGHTED mean Kyte-Doolittle = the hydropathy of the SURFACE.
#                         Directly comparable to sequence GRAVY; the difference between them is the
#                         whole point of this script.
#   hydrophobic_sasa_frac fraction of the exposed surface contributed by hydrophobic residues
#   largest_patch_A2      area of the largest CONTIGUOUS exposed hydrophobic patch. This is the
#                         quantity interface prediction actually uses - assembly needs a patch, not a
#                         high average - so it is the sharpest test of the hypothesis.
#   patch_frac            that patch as a fraction of total SASA (size-normalised)
#
# THREE CONFOUNDS ARE HANDLED EXPLICITLY, because without them the test is worth no more than GRAVY:
#   size        bigger proteins have more surface AND assemble more often -> every correlation is
#               reported as a PARTIAL Spearman holding log10(monomer mass) constant, and the raw value
#               is shown next to it so the difference is visible.
#   disorder    AlphaFold's low-confidence tails have meaningless geometry -> residues below
#               plddt_min are dropped, and the discarded fraction is reported per protein.
#   membrane    membrane-associated proteins are hydrophobic AND often in large particles, which could
#               produce the whole correlation on its own -> the test is repeated with them excluded.
#
# WHAT A RESULT WOULD MEAN
#   surface hydrophobicity correlates but sequence GRAVY does not, or correlates more strongly
#       -> supports the assembly reading; the signal lives on the surface, where interfaces are.
#   both correlate about equally
#       -> the surface adds nothing over composition; most likely a compositional confound, not interfaces.
#   neither survives partialling on mass, or the effect dies when membrane proteins are removed
#       -> report the exclusion of the retention artefact and nothing further.
#
# REQUIRES AlphaFold structures in output/hydropro/structures/, as AF-<ACC>-F1-model_vN.pdb or .cif
# (both formats are read). Get them with alphafold_import_tar() for a whole proteome, or with
# hydropro_fetch() from scripts/hydropro_ffo.R, or download manually and use
# hydropro_import_structures(); everything here works offline once the files are present.
#
# USAGE (RStudio console, project open):
#   source(here::here("scripts", "surface_hydrophobicity.R"))
#   surface_selftest()                              # verify the SASA geometry first
#   # --- getting structures, easiest first ---
#   alphafold_import_tar("C:/Users/<you>/Downloads/UP000000625_83333_ECOLI_v6.tar")  # whole proteome, ONE download
#   surface_priority_list("ATP", n = 200)           # or: a short, defensibly chosen accession list
#   surface_hydrophobicity(metabolite = "ATP", max_proteins = 300)   # slow; cached and resumable
#   surface_vs_elution(metabolite = "ATP")          # the test, against elution deviation
#   surface_vs_elution_compare("ATP")               # calibrated vs all proteins, side by side
#
# COST: roughly 1-3 s per protein at the default 92 test points. Start with max_proteins = 300 to see
# whether anything is there; the cache is additive, so raising the cap later only computes the new ones.
#
# OUTPUT (output/surface/):
#   surface_metrics.csv          per protein: SASA, surface GRAVY, hydrophobic fraction, largest patch
#   surface_vs_elution_<set>.csv       the join with globularity_check, per protein set
#   surface_vs_elution_stats_<set>.csv the correlation table, per protein set
#   surface_vs_elution_COMPARISON.csv  calibrated vs all-proteins, side by side
#   surface_vs_elution_<set>.pdf       each metric vs elution deviation, per protein set
#                                      <set> is "calibrated" or "allproteins" - the two never overwrite
#   surface_vs_gravy.pdf         surface GRAVY vs sequence GRAVY - how different are they really?
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })

.sf_dir <- function(...) here("output", "surface", ...)
.struct_dir <- function() here("output", "hydropro", "structures")

# Kyte-Doolittle hydropathy; positive = hydrophobic
.KD_SF <- c(A =  1.8, R = -4.5, N = -3.5, D = -3.5, C =  2.5, Q = -3.5, E = -3.5, G = -0.4,
            H = -3.2, I =  4.5, L =  3.8, K = -3.9, M =  1.9, F =  2.8, P = -1.6, S = -0.8,
            T = -0.7, W = -0.9, Y = -1.3, V =  4.2)
# residues counted as hydrophobic for the PATCH analysis (the classic apolar set)
.HYDROPHOBIC <- c("A", "V", "L", "I", "M", "F", "W", "C", "Y")
.AA3 <- c(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G", HIS="H", ILE="I",
          LEU="L", LYS="K", MET="M", PHE="F", PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V")
# van der Waals radii (A)
.VDW <- c(C = 1.70, N = 1.55, O = 1.52, S = 1.80, H = 1.20, P = 1.80)

# ---- structure parsing -----------------------------------------------------------------------------
# Both AlphaFold formats are supported: the legacy fixed-column PDB, and mmCIF, which newer AlphaFold DB
# releases ship instead. In both, the per-residue pLDDT lives in the B-factor field.
.cif_atoms <- function(file) {
  ln <- tryCatch(readLines(file, warn = FALSE), error = function(e) character(0))
  if (!length(ln)) return(NULL)
  # the atom_site loop: the "_atom_site.<field>" lines give the column order for the rows beneath
  hdr <- grep("^_atom_site\\.", ln)
  if (!length(hdr)) return(NULL)
  fields <- sub("^_atom_site\\.", "", trimws(ln[hdr]))
  body <- ln[(max(hdr) + 1):length(ln)]
  body <- body[grepl("^(ATOM|HETATM)\\s", body)]
  if (!length(body)) return(NULL)
  M <- do.call(rbind, strsplit(trimws(body), "[ \t]+"))
  if (is.null(M) || ncol(M) < length(fields)) return(NULL)
  cn <- function(...) { h <- intersect(c(...), fields); if (length(h)) which(fields == h[1])[1] else NA_integer_ }
  i_el <- cn("type_symbol"); i_at <- cn("label_atom_id", "auth_atom_id")
  i_rn <- cn("label_comp_id", "auth_comp_id"); i_sq <- cn("label_seq_id", "auth_seq_id")
  i_x  <- cn("Cartn_x"); i_y <- cn("Cartn_y"); i_z <- cn("Cartn_z"); i_b <- cn("B_iso_or_equiv")
  if (anyNA(c(i_at, i_rn, i_sq, i_x, i_y, i_z))) return(NULL)
  A <- data.table(
    atom    = gsub('"', "", M[, i_at]),
    resname = M[, i_rn],
    resnum  = suppressWarnings(as.integer(M[, i_sq])),
    x       = suppressWarnings(as.numeric(M[, i_x])),
    y       = suppressWarnings(as.numeric(M[, i_y])),
    z       = suppressWarnings(as.numeric(M[, i_z])),
    plddt   = if (is.na(i_b)) NA_real_ else suppressWarnings(as.numeric(M[, i_b])),
    element = if (is.na(i_el)) "" else M[, i_el])
  A[is.finite(x) & is.finite(y) & is.finite(z)]
}

.pdb_atoms <- function(file) {
  if (grepl("\\.(cif|mmcif)$", file, ignore.case = TRUE)) {
    A <- .cif_atoms(file)
  } else {
    ln <- tryCatch(readLines(file, warn = FALSE), error = function(e) character(0))
    ln <- ln[startsWith(ln, "ATOM  ")]
    if (!length(ln)) return(NULL)
    A <- data.table(
      atom    = trimws(substr(ln, 13, 16)),
      resname = trimws(substr(ln, 18, 20)),
      resnum  = suppressWarnings(as.integer(substr(ln, 23, 26))),
      x       = suppressWarnings(as.numeric(substr(ln, 31, 38))),
      y       = suppressWarnings(as.numeric(substr(ln, 39, 46))),
      z       = suppressWarnings(as.numeric(substr(ln, 47, 54))),
      plddt   = suppressWarnings(as.numeric(substr(ln, 61, 66))),
      element = trimws(substr(ln, 77, 78)))
    A <- A[is.finite(x) & is.finite(y) & is.finite(z)]
  }
  if (is.null(A) || !nrow(A)) return(NULL)
  A[element == "" | is.na(element), element := substr(atom, 1, 1)]   # fall back to the atom name
  A[, aa := unname(.AA3[resname])]
  A <- A[!is.na(aa) & element != "H"]                                # heavy atoms of standard residues
  A[, radius := unname(.VDW[element])]
  A[is.na(radius), radius := 1.70]
  # side-chain flag: assembly patches are made of side chains, not the backbone
  A[, sidechain := !(atom %in% c("N", "CA", "C", "O", "OXT"))]
  A[]
}

# ---- Shrake-Rupley solvent-accessible surface area -------------------------------------------------
# Each atom is given a sphere of radius (vdW + probe); test points on that sphere are accessible unless
# they fall inside a neighbour's sphere. A cell list keeps the neighbour search local, so cost grows
# linearly with the number of atoms rather than quadratically.
.sphere_points <- function(n) {              # near-uniform points via the golden spiral
  i <- seq_len(n) - 0.5
  phi <- acos(1 - 2 * i / n); theta <- pi * (1 + sqrt(5)) * i
  cbind(cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi))
}

# SELF-TEST: for two overlapping spheres the accessible area is analytic - the buried spherical cap on
# sphere 1 has height h = r1 - (d^2 + r1^2 - r2^2)/(2d), so SASA = 4*pi*r1^2 - 2*pi*r1*h. Checking the
# point-counting against that verifies the geometry and the radii before any protein is touched.
# Validated at the default 92 points: per-atom error 0.2-5.7%, and under 1% at 252 points. Per-atom
# errors are random and average out over the thousands of atoms in a protein, so totals are far more
# accurate than the per-atom worst case; 92 is the classical Shrake-Rupley value.
surface_selftest <- function(n_points = c(92, 252), probe = 1.4) {
  cases <- data.table(R1 = c(1.70, 1.70, 1.70, 1.55, 1.70),
                      R2 = c(1.70, 1.52, 1.80, 1.70, 1.70),
                      d  = c(3.00, 2.40, 4.00, 3.50, 1.50))
  exact <- function(R1, R2, d) {
    r1 <- R1 + probe; r2 <- R2 + probe
    if (d >= r1 + r2) return(4 * pi * r1^2)
    h <- r1 - (d^2 + r1^2 - r2^2) / (2 * d)
    4 * pi * r1^2 - 2 * pi * r1 * h
  }
  cases[, exact_A2 := mapply(exact, R1, R2, d)]
  for (np in n_points) {
    cases[[paste0("n", np)]] <- mapply(function(R1, R2, d) {
      A <- data.table(x = c(0, d), y = 0, z = 0, radius = c(R1, R2))
      .sasa_atoms(A, probe = probe, n_points = np)[1]
    }, cases$R1, cases$R2, cases$d)
    cases[[paste0("err", np, "_pct")]] <- round(100 * abs(cases[[paste0("n", np)]] - cases$exact_A2) / cases$exact_A2, 1)
  }
  message("SASA SELF-TEST - point counting vs the analytic two-sphere solution:")
  print(cases[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 2) else z)])
  message(sprintf("Isolated carbon atom should be 4*pi*(1.70+%.1f)^2 = %.2f A^2; computed %.2f.",
                  probe, 4 * pi * (1.70 + probe)^2,
                  .sasa_atoms(data.table(x = 0, y = 0, z = 0, radius = 1.70), probe = probe, n_points = max(n_points))[1]))
  invisible(cases)
}

.sasa_atoms <- function(A, probe = 1.4, n_points = 92) {
  n <- nrow(A); if (!n) return(numeric(0))
  xyz <- as.matrix(A[, .(x, y, z)]); R <- A$radius + probe
  SP  <- .sphere_points(n_points)
  cut <- 2 * max(R)                                   # no atom beyond this can occlude another
  # cell list
  mn  <- apply(xyz, 2, min)
  cid <- floor(sweep(xyz, 2, mn) / cut)
  key <- paste(cid[, 1], cid[, 2], cid[, 3], sep = ",")
  buckets <- split(seq_len(n), key)
  offs <- as.matrix(expand.grid(-1:1, -1:1, -1:1))
  out <- numeric(n)
  for (i in seq_len(n)) {
    nb <- unlist(buckets[paste(cid[i, 1] + offs[, 1], cid[i, 2] + offs[, 2],
                               cid[i, 3] + offs[, 3], sep = ",")], use.names = FALSE)
    nb <- nb[!is.na(nb) & nb != i]
    if (length(nb)) {                                  # keep only spheres that can actually reach
      d2 <- colSums((t(xyz[nb, , drop = FALSE]) - xyz[i, ])^2)
      nb <- nb[d2 < (R[i] + R[nb])^2]
    }
    if (!length(nb)) { out[i] <- 4 * pi * R[i]^2; next }
    P <- sweep(SP * R[i], 2, xyz[i, ], "+")            # test points on atom i's sphere
    # a point is buried if it lies inside ANY neighbour sphere
    buried <- logical(nrow(P))
    for (j in nb) {
      buried <- buried | (colSums((t(P) - xyz[j, ])^2) < R[j]^2)
      if (all(buried)) break
    }
    out[i] <- 4 * pi * R[i]^2 * mean(!buried)
  }
  out
}

# ---- largest contiguous exposed hydrophobic patch ---------------------------------------------------
# Single-linkage clustering of EXPOSED apolar side-chain atoms; the patch area is the summed SASA of a
# cluster. A protein assembles through a patch, not through a high average, so this is the metric that
# best matches the hypothesis being tested.
.largest_patch <- function(A, sasa, link = 5.0, min_atom_sasa = 1.0) {
  sel <- which(A$sidechain & A$aa %in% .HYDROPHOBIC & sasa > min_atom_sasa)
  if (length(sel) < 3) return(list(area = 0, n_atoms = length(sel)))
  xyz <- as.matrix(A[sel, .(x, y, z)]); m <- length(sel)
  if (m > 4000) return(list(area = NA_real_, n_atoms = m))          # guard: implausible for one chain
  parent <- seq_len(m)
  find <- function(a) { while (parent[a] != a) { parent[a] <<- parent[parent[a]]; a <- parent[a] }; a }
  for (a in seq_len(m - 1)) {
    d2 <- colSums((t(xyz[(a + 1):m, , drop = FALSE]) - xyz[a, ])^2)
    for (b in which(d2 < link^2)) {
      ra <- find(a); rb <- find(a + b)
      if (ra != rb) parent[rb] <- ra
    }
  }
  roots <- vapply(seq_len(m), find, integer(1))
  areas <- tapply(sasa[sel], roots, sum)
  list(area = max(areas), n_atoms = m)
}

# ---- 0b. getting the structures ---------------------------------------------------------------------
# EASIEST ROUTE - ONE FILE, NO LIST. AlphaFold DB publishes whole proteomes as a single tar, so there is
# no need to pick proteins by hand at all. For E. coli K-12 MG1655 (UniProt proteome UP000000625,
# taxid 83333) that is roughly 4300 structures in one download:
#
#   https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000000625_83333_ECOLI_v6.tar
#
# Check the version suffix against the directory listing at
#   https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/
# before downloading - it advances with each AlphaFold DB release (v6 as of this writing); newer releases
# may drop the legacy PDB format and ship mmCIF only, which is why both are parsed here.
# Save the tar anywhere and point alphafold_import_tar() at it. If the tar is blocked but single files
# are not, use surface_priority_list() below to get a short, sensibly chosen accession list instead.
# Path resolution deserves its own helper because "~" is a trap on Windows: R maps it to the user's
# DOCUMENTS folder, not to the home directory, so the "~/Downloads/..." that works on macOS and Linux
# resolves to C:/Users/<you>/Documents/Downloads and fails. Rather than make the caller retype paths,
# try the obvious candidates and, on failure, say exactly what was tried and what tar files are actually
# sitting in those folders.
.resolve_file <- function(f, pattern = NULL) {
  cand <- unique(c(f, path.expand(f), here(f)))
  base <- basename(f)
  homes <- unique(c(Sys.getenv("USERPROFILE"), Sys.getenv("HOME"), path.expand("~"),
                    dirname(path.expand("~"))))
  homes <- homes[nzchar(homes) & dir.exists(homes)]
  dirs  <- unique(c(homes, file.path(homes, c("Downloads", "Desktop", "Documents",
                                              "Documents/Downloads")), here("data", "raw"), getwd()))
  dirs  <- dirs[dir.exists(dirs)]
  cand  <- unique(c(cand, file.path(dirs, base)))
  hit   <- cand[file.exists(cand) & !dir.exists(cand)]
  if (length(hit)) return(hit[1])
  # nothing matched the given name - show what IS there, so the real filename is obvious
  found <- unlist(lapply(dirs, function(d)
    list.files(d, pattern = if (is.null(pattern)) "\\.tar$" else pattern, full.names = TRUE)))
  stop("Could not find '", base, "'.\nTried:\n  ", paste(utils::head(cand, 12), collapse = "\n  "),
       if (length(found))
         paste0("\n\nBut these files DO exist - pass one of these paths instead:\n  ",
                paste(utils::head(found, 12), collapse = "\n  "))
       else paste0("\n\nNo matching file in any of:\n  ", paste(dirs, collapse = "\n  "),
                   "\nNOTE on Windows R, '~' means your DOCUMENTS folder, not your home folder - give the ",
                   "full path instead, e.g. \"C:/Users/<you>/Downloads/<file>.tar\" (forward slashes)."),
       call. = FALSE)
}

alphafold_import_tar <- function(tarfile, dest = .struct_dir(), keep_pdb_only = TRUE) {
  tarfile <- .resolve_file(tarfile)
  message("Using ", tarfile, " (", round(file.size(tarfile) / 1e9, 2), " GB).")
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  message("Unpacking ", basename(tarfile), " -> ", dest, " (this takes a few minutes)...")
  utils::untar(tarfile, exdir = dest)
  # some releases pack into a subdirectory - flatten so everything sits in `dest`
  sub <- list.files(dest, pattern = "^AF-.*\\.(pdb|cif)(\\.gz)?$", full.names = TRUE, recursive = TRUE)
  moved <- sub[dirname(sub) != normalizePath(dest, winslash = "/", mustWork = FALSE)]
  if (length(moved)) {
    file.rename(moved, file.path(dest, basename(moved)))
    message("Flattened ", length(moved), " file(s) out of subdirectories.")
  }
  # stream the decompression: guessing the uncompressed size and reading that many bytes silently
  # TRUNCATES any file that compressed better than the guess
  .gunzip <- function(gz, out) {
    ci <- gzfile(gz, "rb"); co <- file(out, "wb")
    on.exit({ try(close(ci), silent = TRUE); try(close(co), silent = TRUE) })
    repeat { b <- readBin(ci, "raw", 1e6); if (!length(b)) break; writeBin(b, co) }
  }
  gz <- list.files(dest, pattern = "^AF-.*\\.(pdb|cif)\\.gz$", full.names = TRUE)
  if (length(gz)) {
    message("Decompressing ", length(gz), " file(s)...")
    for (g in gz) { out <- sub("\\.gz$", "", g); if (!file.exists(out)) .gunzip(g, out); unlink(g) }
  }
  npdb <- length(list.files(dest, pattern = "^AF-.*\\.pdb$"))
  ncif <- length(list.files(dest, pattern = "^AF-.*\\.cif$"))
  # Only discard mmCIF when PDB is actually present. Newer AlphaFold DB releases have been dropping the
  # legacy PDB format, so deleting the mmCIF unconditionally could throw away the entire download.
  if (keep_pdb_only && npdb > 0 && ncif > 0) {
    unlink(list.files(dest, pattern = "^AF-.*\\.cif$", full.names = TRUE))
    message("Removed ", ncif, " redundant .cif file(s) (PDB copies are present).") ; ncif <- 0L
  }
  message(npdb, " PDB and ", ncif, " mmCIF structure(s) now in ", dest,
          ". Both formats are read by this script.")
  if (npdb + ncif == 0L)
    warning("Nothing usable was unpacked - check that the tar really is an AlphaFold proteome archive.",
            call. = FALSE)
  invisible(npdb + ncif)
}

# ---- 0c. a defensible short list, if you must download one protein at a time -------------------------
# A NOTE ON STUDY DESIGN, because "the 200 proteins known to have interaction partners" is the wrong
# sample for this particular question. The test correlates surface hydrophobicity against elution
# deviation. Known interactors mostly elute LARGER than their monomer - that is what having partners
# means chromatographically - so selecting on interaction status truncates the range of the very
# variable being predicted, and a correlation computed inside that truncated range is biased towards
# zero and uninterpretable either way.
#
# The default here therefore takes a sample STRATIFIED across the observed elution-deviation range, with
# interactors and non-interactors balanced inside each stratum. That buys two tests instead of one:
#   (1) the correlation, now over the full deviation range where it is meaningful;
#   (2) a cleaner group comparison - do proteins with curated partners have greasier SURFACES than
#       proteins without, at matched size? That is arguably the sharper test, and it needs both groups.
# Pass mode = "interactors_only" if you want the list you originally asked for; it is kept because it is
# still the right sample for a purely descriptive "what do interfaces look like" survey.
#
# Ribosomal proteins are excluded in every mode: they are one huge assembly, they would dominate any
# interactor set, and their elution says more about the particle than about their own surfaces.
surface_priority_list <- function(metabolite, n = 200, mode = c("stratified", "interactors_only"),
                                  condition = NULL, restrict_to_calibrated = TRUE,
                                  complex_portal_file = NULL, seed = 1) {
  mode <- match.arg(mode)
  gf <- here("output", paste0("PCM_ctrl_vs_", metabolite), "tables", "globularity_check.txt")
  if (!file.exists(gf)) stop("No globularity_check.txt for ", metabolite, " - run globularity_check() first.")
  G <- fread(gf)
  if ("condition" %in% names(G)) {
    cn <- unique(as.character(G$condition))
    cc <- if (!is.null(condition)) condition else { x <- cn[grepl("ctrl|control|ref", cn, ignore.case = TRUE)][1]; if (is.na(x)) cn[1] else x }
    G <- G[condition == cc]
  }
  if (restrict_to_calibrated && "in_calibrated_range" %in% names(G)) G <- G[in_calibrated_range %in% TRUE]
  G <- G[is.finite(apparent_mw_kDa) & is.finite(expected_mw_kDa) & expected_mw_kDa > 0 & apparent_mw_kDa > 0]
  G[, dev_log2 := log2(apparent_mw_kDa / expected_mw_kDa)]

  # annotation: gene symbol, curated complex membership, curated binary interactions
  U <- NULL; sf <- here("output", "uniprot_annotation_shared.RData")
  if (file.exists(sf)) { e <- new.env(); load(sf, envir = e); U <- as.data.table(e$.uniprot_all) }
  if (is.null(U)) stop("No UniProt cache - render a comparison first.")
  U[, acc := as.character(input_id)]
  U[, gene1 := tolower(sub(" .*$", "", as.character(gene_names)))]
  U[, has_interaction := if ("cc_interaction" %in% names(U)) !is.na(cc_interaction) & nzchar(cc_interaction) else FALSE]
  cp <- .complex_portal_members_sf(complex_portal_file)
  U[, in_complex_portal := if (is.null(cp)) FALSE else toupper(acc) %in% cp]
  keepU <- intersect(c("acc", "gene1", "protein_name", "has_interaction", "in_complex_portal"), names(U))
  J <- merge(G, U[, ..keepU], by.x = "protein_id", by.y = "acc")
  if (!nrow(J)) stop("No overlap between the SEC table and the UniProt cache.")

  # exclude ribosomal proteins
  is_rib <- grepl("^rp[slm][a-z]$", J$gene1) |
            grepl("(30S|50S|40S|60S)? ?ribosomal protein", as.character(J$protein_name), ignore.case = TRUE)
  message("Excluding ", sum(is_rib), " ribosomal protein(s).")
  J <- J[!is_rib]
  J[, any_partner := has_interaction %in% TRUE | in_complex_portal %in% TRUE]
  message(sprintf("%d protein(s) available; %d (%.0f%%) have a curated partner (Complex Portal and/or UniProt interactions).",
                  nrow(J), sum(J$any_partner), 100 * mean(J$any_partner)))

  set.seed(seed)
  if (mode == "interactors_only") {
    P <- J[any_partner == TRUE]
    setorder(P, -abs(dev_log2))                       # most informative elution behaviour first
    sel <- head(P, n)
    message("mode = 'interactors_only': ", nrow(sel), " protein(s). NOTE this sample is selected on ",
            "interaction status, which truncates the elution-deviation range - use it descriptively, ",
            "not for the correlation.")
  } else {
    # 10 strata across the deviation range, half interactors / half not within each
    J[, stratum := cut(dev_log2, breaks = unique(stats::quantile(dev_log2, seq(0, 1, 0.1), na.rm = TRUE)),
                       include.lowest = TRUE, labels = FALSE)]
    per <- ceiling(n / (2 * length(unique(J$stratum[is.finite(J$stratum)]))))
    sel <- J[is.finite(stratum), .SD[sample(.N, min(.N, per))], by = .(stratum, any_partner)]
    if (nrow(sel) > n) sel <- sel[sample(.N, n)]
    message("mode = 'stratified': ", nrow(sel), " protein(s) spread across ",
            length(unique(sel$stratum)), " deviation strata, ",
            sum(sel$any_partner), " with a curated partner and ", sum(!sel$any_partner), " without.")
  }
  setorder(sel, -any_partner, -abs(dev_log2))
  dir.create(.sf_dir(), recursive = TRUE, showWarnings = FALSE)
  out_cols <- intersect(c("protein_id", "gene1", "protein_name", "expected_mw_kDa", "apparent_mw_kDa",
                          "dev_log2", "class", "has_interaction", "in_complex_portal", "any_partner"), names(sel))
  fwrite(sel[, ..out_cols], .sf_dir("surface_priority_list.csv"))
  writeLines(sel$protein_id, .sf_dir("surface_priority_accessions.txt"))
  message("Wrote ", .sf_dir("surface_priority_list.csv"), " (annotated) and ",
          .sf_dir("surface_priority_accessions.txt"), " (one accession per line).")
  message("Download each as https://alphafold.ebi.ac.uk/files/AF-<ACCESSION>-F1-model_v4.pdb into ",
          .struct_dir(), " - or take the whole proteome in one tar, see alphafold_import_tar().")
  invisible(sel)
}

# same Complex Portal parser as globularity_category_annotation.R, kept local so this script stands alone
.complex_portal_members_sf <- function(file = NULL) {
  f <- file
  if (is.null(f)) {
    cand <- list.files(here("data", "raw"), pattern = "^Complex_portal.*\\.tsv$", full.names = TRUE)
    if (!length(cand)) { message("No Complex Portal export in data/raw - using UniProt interactions only."); return(NULL) }
    f <- cand[1]
  } else if (!file.exists(f)) { f2 <- here("data", "raw", f); if (!file.exists(f2)) return(NULL); f <- f2 }
  cp <- tryCatch(as.data.table(read.csv(f, sep = "\t", header = TRUE, check.names = FALSE)), error = function(e) NULL)
  if (is.null(cp) || !nrow(cp)) return(NULL)
  idc <- grep("identifier", names(cp), ignore.case = TRUE, value = TRUE)
  col <- if (length(idc)) grep("molecul|stoichiom", idc, ignore.case = TRUE, value = TRUE)[1] else NA_character_
  if (is.na(col)) col <- if (length(idc)) idc[1] else NA_character_
  if (is.na(col)) return(NULL)
  unique(toupper(unlist(regmatches(cp[[col]],
    gregexpr("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", cp[[col]])))))
}

# ---- 1. per-protein surface metrics -----------------------------------------------------------------
surface_hydrophobicity <- function(ids = NULL, metabolite = NULL, max_proteins = Inf,
                                   plddt_min = 70, probe = 1.4, n_points = 92,
                                   overwrite = FALSE, verbose_every = 25) {
  sdir <- .struct_dir()
  if (!dir.exists(sdir)) stop("No structure directory at ", sdir,
                              ". Fetch structures first: source('scripts/hydropro_ffo.R'); hydropro_fetch('ATP').")
  files <- list.files(sdir, pattern = "^AF-.*\\.(pdb|cif)$", full.names = TRUE)
  if (!length(files)) stop("No AF-*.pdb or AF-*.cif files in ", sdir,
                           " - see alphafold_import_tar(), hydropro_fetch() or hydropro_import_structures().")
  have <- sub("^AF-([^-]+)-.*$", "\\1", basename(files))   # AF-<ACC>-F1-model_vN.(pdb|cif)
  if (is.null(ids) && !is.null(metabolite)) {
    f <- here("output", paste0("PCM_ctrl_vs_", metabolite), "tables", "globularity_check.txt")
    if (file.exists(f)) ids <- unique(as.character(fread(f)$protein_id))
  }
  keep <- if (is.null(ids)) seq_along(files) else which(have %in% ids)
  if (!length(keep)) stop("None of the requested proteins has a structure in ", sdir, ".")
  files <- files[keep]; have <- have[keep]
  if (length(files) > max_proteins) {
    message("Limiting to the first ", max_proteins, " of ", length(files), " structures.")
    files <- head(files, max_proteins); have <- head(have, max_proteins)
  }

  cache <- .sf_dir("surface_metrics.csv")
  done <- if (!overwrite && file.exists(cache)) fread(cache) else NULL
  if (!is.null(done) && nrow(done)) {
    todo <- !(have %in% done$protein_id)
    message("Cache holds ", nrow(done), " protein(s); computing ", sum(todo), " new one(s). overwrite = TRUE to redo all.")
    files <- files[todo]; have <- have[todo]
  }
  if (!length(files)) { message("Nothing to compute - all requested proteins are cached."); return(invisible(done)) }

  message("Computing SASA for ", length(files), " structure(s) with ", n_points,
          " test points per atom. This is the slow step; results are cached.")
  rows <- vector("list", length(files)); t0 <- Sys.time()
  for (i in seq_along(files)) {
    A <- .pdb_atoms(files[i])
    if (is.null(A) || !nrow(A)) { rows[[i]] <- data.table(protein_id = have[i], ok = FALSE); next }
    n_all <- length(unique(A$resnum))
    # drop low-confidence residues: AlphaFold tails have no meaningful geometry, and including them
    # would inflate apparent exposed surface with something that is not a real surface
    if (is.finite(plddt_min)) A <- A[is.na(plddt) | plddt >= plddt_min]
    if (nrow(A) < 20) { rows[[i]] <- data.table(protein_id = have[i], ok = FALSE); next }
    n_kept <- length(unique(A$resnum))
    s  <- .sasa_atoms(A, probe = probe, n_points = n_points)
    kd <- unname(.KD_SF[A$aa]); kd[!is.finite(kd)] <- 0
    tot <- sum(s)
    pat <- .largest_patch(A, s)
    rows[[i]] <- data.table(
      protein_id            = have[i], ok = TRUE,
      n_residues_total      = n_all,
      n_residues_used       = n_kept,
      frac_residues_dropped = round(1 - n_kept / max(n_all, 1), 4),
      mean_plddt            = round(mean(A$plddt, na.rm = TRUE), 2),
      sasa_total            = round(tot, 1),
      # SASA-weighted mean hydropathy = the hydropathy of the SURFACE (compare with sequence GRAVY)
      surface_gravy         = round(sum(s * kd) / max(tot, 1e-9), 4),
      # unweighted sequence GRAVY over the SAME residues, so the two are like-for-like
      sequence_gravy_here   = round(mean(unname(.KD_SF[A[!duplicated(resnum)]$aa]), na.rm = TRUE), 4),
      hydrophobic_sasa_frac = round(sum(s[A$aa %in% .HYDROPHOBIC]) / max(tot, 1e-9), 4),
      largest_patch_A2      = round(pat$area, 1),
      patch_frac            = round(pat$area / max(tot, 1e-9), 4))
    if (i %% verbose_every == 0)
      message("  ... ", i, "/", length(files), "  (", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min)")
  }
  D <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (!is.null(done) && nrow(done)) D <- rbindlist(list(done, D), use.names = TRUE, fill = TRUE)
  dir.create(.sf_dir(), recursive = TRUE, showWarnings = FALSE)
  fwrite(D, cache)
  message("Surface metrics for ", sum(D$ok %in% TRUE), " protein(s) -> ", cache,
          "  (", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min)")
  invisible(D)
}

# ---- 2. the test ------------------------------------------------------------------------------------
# rank-based partial correlation, holding one covariate constant
.pcor_sf <- function(x, y, z) {
  ok <- is.finite(x) & is.finite(y) & is.finite(z); n <- sum(ok)
  if (n < 10) return(list(rho = NA_real_, p = NA_real_, n = n))
  rx <- rank(x[ok]); ry <- rank(y[ok]); rz <- rank(z[ok])
  rxy <- stats::cor(rx, ry); rxz <- stats::cor(rx, rz); ryz <- stats::cor(ry, rz)
  den <- sqrt((1 - rxz^2) * (1 - ryz^2))
  if (!is.finite(den) || den <= 0) return(list(rho = NA_real_, p = NA_real_, n = n))
  r <- (rxy - rxz * ryz) / den
  tt <- r * sqrt((n - 3) / max(1 - r^2, .Machine$double.eps))
  list(rho = r, p = 2 * stats::pt(-abs(tt), df = n - 3), n = n)
}

# Rank-based partial correlation on SEVERAL covariates at once: rank everything, regress the ranks of x
# and of y on the ranks of the covariates, and correlate what is left. With Z = (mass, sequence GRAVY)
# this answers the question the single-covariate version cannot - does the SURFACE carry information
# beyond what the amino-acid composition already provides? Surface and sequence hydropathy are
# correlated with each other, so comparing their two separate correlations does not settle that.
.pcor_multi <- function(x, y, Z) {
  Z <- as.matrix(Z)
  ok <- is.finite(x) & is.finite(y) & apply(is.finite(Z), 1, all); n <- sum(ok)
  k <- ncol(Z)
  if (n < 10 + k) return(list(rho = NA_real_, p = NA_real_, n = n))
  rx <- rank(x[ok]); ry <- rank(y[ok])
  RZ <- apply(Z[ok, , drop = FALSE], 2, rank)
  ex <- stats::residuals(stats::lm(rx ~ RZ)); ey <- stats::residuals(stats::lm(ry ~ RZ))
  r <- suppressWarnings(stats::cor(ex, ey))
  if (!is.finite(r)) return(list(rho = NA_real_, p = NA_real_, n = n))
  df <- n - 2 - k
  tt <- r * sqrt(df / max(1 - r^2, .Machine$double.eps))
  list(rho = r, p = 2 * stats::pt(-abs(tt), df = df), n = n)
}

surface_vs_elution <- function(metabolite, condition = NULL, restrict_to_calibrated = TRUE,
                               exclude_membrane = TRUE, save_plots = TRUE) {
  cache <- .sf_dir("surface_metrics.csv")
  if (!file.exists(cache)) stop("No surface metrics yet - run surface_hydrophobicity() first.")
  S <- fread(cache)[ok %in% TRUE]
  gf <- here("output", paste0("PCM_ctrl_vs_", metabolite), "tables", "globularity_check.txt")
  if (!file.exists(gf)) stop("No globularity_check.txt for ", metabolite, " - run globularity_check() first.")
  G <- fread(gf)
  if ("condition" %in% names(G)) {
    cn <- unique(as.character(G$condition))
    cc <- if (!is.null(condition)) condition else { x <- cn[grepl("ctrl|control|ref", cn, ignore.case = TRUE)][1]; if (is.na(x)) cn[1] else x }
    G <- G[condition == cc]; message("Using the '", cc, "' rows of the SEC table.")
  }
  if (restrict_to_calibrated && "in_calibrated_range" %in% names(G)) {
    n0 <- nrow(G); G <- G[in_calibrated_range %in% TRUE]
    message("Restricted to the calibrated MW interval: ", nrow(G), " of ", n0, " proteins.")
  }
  if (!all(c("apparent_mw_kDa", "expected_mw_kDa") %in% names(G)))
    stop("globularity_check.txt lacks apparent_mw_kDa/expected_mw_kDa.")
  G <- G[is.finite(apparent_mw_kDa) & is.finite(expected_mw_kDa) & expected_mw_kDa > 0 & apparent_mw_kDa > 0]
  # the SAME deviation the GRAVY test uses, so the two are directly comparable
  G[, dev_log2 := log2(apparent_mw_kDa / expected_mw_kDa)]
  J <- merge(S, G[, .(protein_id, dev_log2, expected_mw_kDa,
                      class = if ("class" %in% names(G)) class else NA_character_)], by = "protein_id")
  if (!nrow(J)) stop("No overlap between the structures and the SEC table.")
  J[, lgm := log10(expected_mw_kDa)]

  # membrane-associated proteins are hydrophobic AND often in large particles - the single most likely
  # way to get a spurious positive without any interface involvement
  memb <- character(0)
  sf <- here("output", "uniprot_annotation_shared.RData")
  if (exclude_membrane && file.exists(sf)) {
    e <- new.env(); load(sf, envir = e); u <- as.data.table(e$.uniprot_all)
    cc_col <- intersect(c("go_c", "subcellular_location", "ft_transmem", "keywords"), names(u))
    if (length(cc_col)) {
      hit <- Reduce(`|`, lapply(cc_col, function(k)
        grepl("membrane|transmembrane|lipoprotein", as.character(u[[k]]), ignore.case = TRUE)))
      memb <- as.character(u$input_id)[which(hit)]
      message("Membrane-associated annotation found for ", length(intersect(memb, J$protein_id)),
              " of the ", nrow(J), " proteins tested.")
    } else message("No membrane-related column in the UniProt cache - the membrane control cannot be run.")
  }
  J[, is_membrane := protein_id %in% memb]

  metrics <- intersect(c("surface_gravy", "hydrophobic_sasa_frac", "patch_frac",
                         "largest_patch_A2", "sequence_gravy_here"), names(J))
  res <- rbindlist(lapply(metrics, function(mt) {
    x <- J[[mt]]; y <- J$dev_log2
    raw <- suppressWarnings(stats::cor.test(x, y, method = "spearman"))
    pc  <- .pcor_sf(x, y, J$lgm)
    K <- J[is_membrane == FALSE]
    pcm <- if (nrow(K) >= 30) .pcor_sf(K[[mt]], K$dev_log2, K$lgm) else list(rho = NA_real_, p = NA_real_, n = nrow(K))
    # THE DECISIVE ONE: hold mass AND sequence GRAVY constant. If a surface metric still correlates
    # after the composition has been removed, the surface is carrying information the sequence does not.
    pcg <- if (mt == "sequence_gravy_here") list(rho = NA_real_, p = NA_real_, n = NA_integer_)
           else .pcor_multi(x, y, cbind(J$lgm, J$sequence_gravy_here))
    data.table(metric = mt, n = sum(is.finite(x) & is.finite(y)),
               rho_raw = round(unname(raw$estimate), 3), p_raw = raw$p.value,
               rho_partial_mass = round(pc$rho, 3), p_partial = pc$p,
               rho_partial_nomembrane = round(pcm$rho, 3), n_nomembrane = pcm$n,
               rho_beyond_composition = round(pcg$rho, 3), p_beyond_composition = pcg$p)
  }))
  dir.create(.sf_dir(), recursive = TRUE, showWarnings = FALSE)
  # Tag every output with the protein set it came from. Previously both settings wrote the SAME file
  # names, so running restrict_to_calibrated = FALSE silently overwrote the calibrated result and the two
  # could not be compared - which is the whole point of running both.
  .tag <- if (isTRUE(restrict_to_calibrated)) "calibrated" else "allproteins"
  .f_csv  <- .sf_dir(paste0("surface_vs_elution_", .tag, ".csv"))
  .f_stat <- .sf_dir(paste0("surface_vs_elution_stats_", .tag, ".csv"))
  .f_pdf  <- .sf_dir(paste0("surface_vs_elution_", .tag, ".pdf"))
  fwrite(J, .f_csv); fwrite(res, .f_stat)
  message("Protein set: ", .tag, " (n = ", nrow(J), "). Outputs carry that suffix, so the two settings do not overwrite each other.")
  if (!isTRUE(restrict_to_calibrated))
    message("CAUTION with restrict_to_calibrated = FALSE: outside the calibrated MW interval the apparent MW is an\n",
            "   EXTRAPOLATION of the standards curve, not a measurement, so dev_log2 is not a measured quantity for\n",
            "   those proteins. Treat this run as a robustness check: if the correlation holds in BOTH sets it is\n",
            "   solid; if it appears only here, it is a property of the extrapolation.")

  message("\nSurface hydrophobicity vs elution deviation (y = log2(apparent/expected); positive = elutes as though heavier):")
  print(res)
  message("\nWHAT EACH COLUMN IS:")
  message("  rho_raw                 plain Spearman, no covariates.")
  message("  rho_partial_mass        Spearman with log10(monomer mass) HELD CONSTANT: rank both variables,")
  message("                          regress each on rank(mass), correlate the residuals. Reads as: among")
  message("                          proteins of the SAME SIZE, does more surface hydrophobicity go with")
  message("                          eluting larger? Needed because size drives both.")
  message("  rho_partial_nomembrane  THE SAME mass-partialled statistic, recomputed on the non-membrane")
  message("                          subset only. It is not 'the correlation without mass control' - both")
  message("                          columns control for mass; this one just drops membrane proteins.")
  message("  rho_beyond_composition  mass AND sequence GRAVY held constant - does the SURFACE add anything")
  message("                          the amino-acid composition did not already say? This is the one that")
  message("                          separates an interface story from a compositional confound.")
  message("\nHOW TO READ THIS:")
  message("  * The RETENTION ARTEFACT predicts a NEGATIVE correlation (sticky protein elutes late, looks small).")
  message("    A positive value excludes it.")
  message("  * If an effect grows in `rho_partial_nomembrane`, membrane proteins were diluting it, not causing it.")
  message("  * `patch_frac` is size-normalised; `largest_patch_A2` is absolute, so it keeps a size component")
  message("    that the mass partial then removes - compare the two before reading anything into either.")
  .bc <- res[metric != "sequence_gravy_here" & is.finite(rho_beyond_composition)]
  if (nrow(.bc)) {
    b <- .bc[which.max(abs(rho_beyond_composition))]
    message(sprintf("\n  BEYOND COMPOSITION: %s keeps rho = %+.3f (p = %.3g) with mass AND sequence GRAVY held constant.",
                    b$metric, b$rho_beyond_composition, b$p_beyond_composition))
    message(if (abs(b$rho_beyond_composition) >= 0.15)
      "  => The surface carries information the sequence does not. That is the result the interface reading needs."
      else "  => Once composition is removed almost nothing is left: the signal is compositional, not an interface property. Report the artefact exclusion only.")
  }
  .best <- res[metric != "sequence_gravy_here"][which.max(abs(rho_partial_mass))]
  if (nrow(.best) && is.finite(.best$rho_partial_mass))
    message(sprintf("\n  Strongest surface metric: %s, partial rho = %+.3f (%.1f%% of the variance). %s",
                    .best$metric, .best$rho_partial_mass, 100 * .best$rho_partial_mass^2,
                    if (abs(.best$rho_partial_mass) < 0.2)
                      "Below |rho| = 0.2 this is too weak to carry an interface claim - report the artefact exclusion only."
                    else "Above |rho| = 0.2 and worth reporting, provided it survives the membrane control."))

  if (save_plots) {
    L <- melt(J[, c("protein_id", "dev_log2", "is_membrane", metrics), with = FALSE],
              id.vars = c("protein_id", "dev_log2", "is_membrane"),
              variable.name = "metric", value.name = "value")
    lab <- res[, .(metric, txt = sprintf("raw %+.3f | partial %+.3f", rho_raw, rho_partial_mass))]
    L <- merge(L, lab, by = "metric")
    L[, facet := paste0(metric, "\n", txt)]
    g1 <- ggplot(L[is.finite(value)], aes(value, dev_log2)) +
      geom_hline(yintercept = 0, colour = "grey60") +
      geom_point(aes(colour = is_membrane), alpha = 0.4, size = 0.8) +
      scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "#E15759"), name = "membrane-associated") +
      geom_smooth(method = "loess", formula = y ~ x, se = TRUE, colour = "black", linewidth = 0.6) +
      facet_wrap(~ facet, scales = "free_x") +
      labs(title = paste0("Surface hydrophobicity vs elution deviation  (", metabolite, " control) - ",
                          if (isTRUE(restrict_to_calibrated)) "WITHIN the calibrated MW interval"
                          else "ALL proteins, incl. extrapolated apparent MW",
                          "  [n = ", nrow(J), "]"),
           subtitle = paste0("y = log2(apparent / expected MW); BELOW 0 = elutes late, i.e. retained. A NEGATIVE trend would be the\n",
                             "column-retention artefact; a positive one is consistent with assembly but cannot prove it. Facet labels give the\n",
                             "raw and the mass-partialled Spearman rho - use the partial one, since size drives both surface area and assembly."),
           x = "metric value", y = "log2(apparent / expected MW)") +
      theme_bw() + theme(legend.position = "top")
    g2 <- ggplot(J[is.finite(surface_gravy) & is.finite(sequence_gravy_here)],
                 aes(sequence_gravy_here, surface_gravy)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
      geom_point(alpha = 0.45, size = 0.9, colour = "steelblue") +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "black", linewidth = 0.6) +
      labs(title = "Does the surface look like the sequence?",
           subtitle = sprintf("Sequence GRAVY averages over BURIED and exposed residues alike; surface GRAVY weights each residue by how\nmuch of it is actually exposed. Spearman rho between them = %.3f. Points below the dashed line are proteins whose\nsurface is more POLAR than their composition suggests - the normal state of a soluble protein, and the reason\nsequence GRAVY cannot test an interface hypothesis.",
                              suppressWarnings(stats::cor(J$sequence_gravy_here, J$surface_gravy, method = "spearman", use = "complete.obs"))),
           x = "sequence GRAVY (same residues)", y = "surface GRAVY (SASA-weighted)") + theme_bw()
    .ok <- tryCatch({ grDevices::pdf(.f_pdf, width = 9, height = 7)
                      print(g1); print(g2); grDevices::dev.off(); TRUE },
                    error = function(e) {
                      try(grDevices::dev.off(), silent = TRUE)
                      message("!! COULD NOT WRITE ", .f_pdf, ": ", conditionMessage(e))
                      message("   The usual cause on Windows is that the PDF is OPEN in a viewer, which locks the file.")
                      message("   Close it and re-run. (This used to fail silently, which is why the figure looked stale.)")
                      FALSE })
    if (.ok) message("Figures -> ", .f_pdf)
  }
  invisible(list(data = J, stats = res))
}

# ---- 3. calibrated vs all-proteins, side by side ----------------------------------------------------
# The direct answer to "what happens if I put the other ~50% back in". Runs the test on BOTH protein
# sets and lays the results next to each other.
#
# WHY IT MATTERS, and why it is a robustness check rather than a bigger analysis: outside the calibrated
# MW interval the apparent MW is an EXTRAPOLATION of the standards curve, not a measurement, so for those
# proteins dev_log2 is a derived guess. Adding them roughly doubles n, which always tightens p-values -
# that on its own means nothing. The informative comparison is the EFFECT SIZE:
#   similar rho in both sets   -> the association does not depend on the extrapolation. Strong evidence,
#                                 and you can report the calibrated number knowing it generalises.
#   much larger with all       -> the extra signal lives in the extrapolated region, i.e. it is a
#                                 property of the standards-curve fit rather than of the proteins.
#   much smaller with all      -> the extrapolated proteins are adding noise, as expected if their
#                                 dev_log2 is largely meaningless.
surface_vs_elution_compare <- function(metabolite, condition = NULL, exclude_membrane = TRUE) {
  message("=== Protein set 1: WITHIN the calibrated MW interval ===")
  a <- surface_vs_elution(metabolite, condition = condition, restrict_to_calibrated = TRUE,
                          exclude_membrane = exclude_membrane, save_plots = TRUE)
  message("\n=== Protein set 2: ALL proteins, including the extrapolated ones ===")
  b <- surface_vs_elution(metabolite, condition = condition, restrict_to_calibrated = FALSE,
                          exclude_membrane = exclude_membrane, save_plots = TRUE)
  A <- data.table::copy(a$stats)[, set := "calibrated"]
  B <- data.table::copy(b$stats)[, set := "all proteins"]
  cols <- intersect(c("metric", "n", "rho_raw", "rho_partial_mass", "rho_partial_nomembrane",
                      "rho_beyond_composition"), names(A))
  M <- merge(A[, ..cols], B[, ..cols], by = "metric", suffixes = c("_cal", "_all"))
  M[, delta_partial := round(rho_partial_mass_all - rho_partial_mass_cal, 3)]
  data.table::setorder(M, -rho_partial_mass_cal)
  fwrite(M, .sf_dir("surface_vs_elution_COMPARISON.csv"))
  message("\n=== SIDE BY SIDE (partial rho, monomer mass held constant) ===")
  print(M[, .(metric, n_cal, n_all,
              rho_cal = rho_partial_mass_cal, rho_all = rho_partial_mass_all, delta = delta_partial)])
  message("Comparison table -> ", .sf_dir("surface_vs_elution_COMPARISON.csv"))
  d <- max(abs(M$delta_partial), na.rm = TRUE)
  message(sprintf("\nLargest shift in partial rho when the extrapolated proteins are added: %.3f", d))
  message(if (d < 0.05)
    "  => The association is INDEPENDENT of the extrapolation. Report the calibrated figures; this run is the robustness check that justifies them."
    else if (any(abs(M$rho_partial_mass_all) > abs(M$rho_partial_mass_cal) + 0.05, na.rm = TRUE))
    "  => The effect is LARGER once extrapolated proteins are included. Be careful: their apparent MW is not a measurement, so the extra signal may belong to the standards-curve fit rather than to the proteins. Report the calibrated figures."
    else
    "  => The effect SHRINKS with the extrapolated proteins included, which is what you would expect if their dev_log2 is largely noise. Report the calibrated figures.")
  message("Note that n roughly doubles, so every p-value falls - that is arithmetic, not evidence. Compare rho, not p.")
  invisible(M)
}
