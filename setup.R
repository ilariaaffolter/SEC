# =============================================================================
# setup.R  -  ONE-TIME environment setup for the SEC differential-analysis project
# =============================================================================
#
# WHAT THIS DOES
#   Creates a private, reproducible R package library for this project using
#   {renv}, installs every package the analysis needs (from CRAN, Bioconductor
#   and GitHub), and writes renv.lock - the file that lets anyone recreate the
#   EXACT same environment later with renv::restore().
#
# WHEN TO RUN IT
#   * The first time you open this project on a new computer, OR
#   * if renv.lock does not exist yet.
#   After it has run once and renv.lock is committed, collaborators only need
#   renv::restore() (see README.md) - they do NOT need to run setup.R again.
#
# HOW TO RUN IT
#   Open SEC.Rproj in RStudio, then in the Console:
#       source("setup.R")
#   (or step through it line by line). It can take 20-40 min on a fresh machine.
#
# REQUIREMENTS
#   * R >= 4.2 recommended, and an internet connection.
#   * Windows: install Rtools (matching your R version).
#   * macOS:   install the Xcode command-line tools  ->  xcode-select --install
#   * Some CCprofiler dependencies need system libraries (notably Rmpfr -> GMP/MPFR).
#     See docs/SETUP_AND_GIT_WORKFLOW.md ("System libraries") if a package fails.
# =============================================================================

message("== SEC project setup starting ==")

## 0. renv ---------------------------------------------------------------------
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

## 1. Initialise renv for this project ----------------------------------------
#    bare = TRUE  -> start from an empty private library; we install exactly what
#                    we list below (nothing is guessed).
#    This also creates .Rprofile + renv/activate.R so the project library is
#    activated automatically every time you open SEC.Rproj in the future.
if (!file.exists("renv.lock")) {
  renv::init(bare = TRUE, restart = FALSE)
}

## 2. Helpers / Bioconductor + GitHub installers ------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE))      install.packages("remotes")

## 3. CRAN packages used by the analysis --------------------------------------
cran_pkgs <- c(
  "data.table", "stringr", "ggplot2", "betareg", "lmtest", "ggrepel",
  "scales", "protti", "tidyverse", "DT", "readxl", "writexl", "UpSetR",
  "here", "fs", "ggpmisc", "rlist", "progress",
  "foreach", "doParallel", "igraph"
)
install.packages(cran_pkgs)

## 4. Bioconductor packages used by the analysis (and by CCprofiler) ----------
bioc_pkgs <- c("PrInCE", "MSnbase", "limma", "preprocessCore", "qvalue")
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

## 5. CCprofiler -- official 'differential' branch, PINNED to an exact commit -
#    This is the official repository (NOT the AnnaPagotto/CCprofilerDiffAnna
#    fork). The fork only added one file of fixes (PPlabFunctions.R); those
#    fixes now live, fully documented, in R/ccprofiler_fixes.R inside THIS
#    project and are sourced on top of the package by the .Rmd. Pinning the
#    commit guarantees the analysis keeps working even if the branch moves on.
remotes::install_github(
  "CCprofiler/CCprofiler",
  ref          = "a6f4f4df4969e6fd5b649ad897941afc1b818346", # tip of 'differential', 2026-06-19
  dependencies = TRUE,
  upgrade      = "never"
)

## 6. Record the exact environment into renv.lock -----------------------------
#    type = "all" captures every installed package (even ones only used
#    indirectly), so the lockfile is a complete, restorable snapshot.
renv::snapshot(type = "all", prompt = FALSE)

message("== Setup complete. ==")
message("renv.lock has been written. COMMIT it (git add renv.lock) so anyone ",
        "can reproduce this environment with renv::restore().")
