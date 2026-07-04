# SEC differential analysis — QTL strains 31 vs 83

Reproducible RStudio project for the **size-exclusion chromatography (SEC)
differential analysis** of yeast strains 31 and 83, built on
[CCprofiler](https://github.com/CCprofiler/CCprofiler) (the `differential`
branch) plus a curated set of lab fixes.

This project uses **[renv](https://rstudio.github.io/renv/)** for a reproducible
package environment and **[here](https://here.r-lib.org/)** for robust file
paths, so it can be handed over and run on another machine with minimal fuss.

> **New to R/Git?** Read **[`docs/SETUP_AND_GIT_WORKFLOW.md`](docs/SETUP_AND_GIT_WORKFLOW.md)** —
> a step-by-step, beginner-friendly guide to installing everything, connecting
> RStudio to GitHub, and the day-to-day pull/commit/push loop.
>
> **Current status & what's next:** see **[`HANDOVER_NOTES.md`](HANDOVER_NOTES.md)**.
> The reproducible foundation (renv, `here`, CCprofiler + consolidated fixes, docs)
> is in place; the cleaned main/experiments split of the `.Rmd` is the next step
> and needs a couple of answers from you (listed there).

---

## Quick start

### A. Run it on your own machine (first time)

1. **Install** R, RStudio and Git (see the detailed guide).
2. **Get the project** — in RStudio: *File ▸ New Project ▸ Version Control ▸ Git*,
   and paste the repository URL. (Or `git clone` it, then open `SEC.Rproj`.)
3. **Open `SEC.Rproj`** (always open the project this way — it makes `here()` and
   `renv` work).
4. **Build the environment** (once per machine), in the RStudio Console:
   ```r
   source("setup.R")
   ```
   This installs every package (CRAN + Bioconductor + CCprofiler from GitHub) and
   writes `renv.lock`. Takes ~20–40 min the first time.
5. **Add your data** to `data/raw/` (see [`data/raw/README.md`](data/raw/README.md)).
6. **Knit** `analysis/DiffAnalysis_yeast_aSYN.Rmd` (the *Knit* button), or run
   it chunk by chunk.

### B. Someone hands the project to you (reproduce the environment)

```r
# in the RStudio Console, with SEC.Rproj open:
install.packages("renv")   # if you don't have it
renv::restore()            # installs the EXACT package versions from renv.lock
```
Then add the data to `data/raw/` and knit. No need to run `setup.R`.

---

## Project structure

```
SEC/
├── SEC.Rproj                  # RStudio project (open this; defines the here() root)
├── setup.R                    # ONE-TIME environment builder (renv + all packages)
├── renv.lock                  # exact package versions (created by setup.R; commit it)
├── README.md                  # this file
├── R/
│   └── ccprofiler_fixes.R     # the lab fixes (IBMT, Beni/Aljaž fixes, Gaussian fits…),
│                              #   sourced on top of CCprofiler. See header for provenance.
├── analysis/
│   ├── ORIGINAL_DiffAnalysis_Norm_QTL_31_83_2026.Rmd  # your original .Rmd, preserved unchanged
│   ├── DiffAnalysis_yeast_aSYN.Rmd   # ← MAIN, cleaned pipeline  (being built — see HANDOVER_NOTES.md)
│   └── experiments.Rmd                    # exploratory / alternative approaches (being built)
├── data/
│   └── raw/                   # your input files go here (git-ignored by default)
├── output/                    # generated figures / tables / RData (git-ignored)
│   ├── figures/  ├── tables/  └── rdata/
└── docs/
    └── SETUP_AND_GIT_WORKFLOW.md   # beginner guide: install, GitHub, renv, daily loop
```

---

## How the CCprofiler "fixes" work

Historically the lab used a fork, `AnnaPagotto/CCprofilerDiffAnna`. A line-by-line
comparison showed that fork is **identical** to the official
`CCprofiler@differential` except for **one added file of fixes**
(`PPlabFunctions.R`) and a package rename.

So this project drops the fork and instead:

1. installs the **official** `CCprofiler` `differential` branch, **pinned to an
   exact commit** (`a6f4f4d`) for reproducibility, and
2. keeps all the fixes in **[`R/ccprofiler_fixes.R`](R/ccprofiler_fixes.R)**, which
   the `.Rmd` loads on top of the package with `source(here::here("R", "ccprofiler_fixes.R"))`.

This makes every fix **visible and editable inside the project** (important for the
upcoming cross-correlation work) while staying faithful to the originals. The file
header documents where each fix came from.

---

## Reproducibility notes

- **renv** isolates this project's packages from your global R library. `renv.lock`
  is the source of truth; commit it whenever you add/upgrade a package
  (`renv::snapshot()`).
- **here** means paths work no matter your working directory — always
  `here::here("data", "raw", "file.tsv")`, never `setwd()` or absolute paths.
- **Data** is not committed by default. The recipient must place the input files in
  `data/raw/` (documented there).
- **System libraries**: a few CCprofiler dependencies need system libs (e.g.
  `Rmpfr` → GMP/MPFR). See the troubleshooting section of the setup guide.

---

## Roadmap

- **Phase 1 (this setup)** — reproducible project: renv + here, official CCprofiler
  pinned, fixes consolidated, main pipeline cleaned, experiments preserved & tidied. ✅
- **Phase 2 (next)** — add a cross-correlation–based differential analysis section as
  an alternative/complement to the statistical testing.
