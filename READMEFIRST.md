# ⚠️ READ ME FIRST — this project uses `renv`

**Do not `install.packages()` by hand and do not just start running the code.**
This project ships a *pinned* package environment with **[renv](https://rstudio.github.io/renv/)**
so that the analysis runs the same way on every machine (CCprofiler, its
Bioconductor dependencies, and dozens of CRAN packages, all at the exact versions
the analysis was validated with). Skipping this step is the #1 way to get errors
that look like bugs but are really "wrong package version".

Three steps, every time, in order:

1. **Open the project the right way.** Double-click **`SEC.Rproj`** (or in RStudio:
   *File ▸ Open Project*). Opening the folder or a single file is not enough — the
   `.Rproj` is what makes `renv` and `here()` work.
2. **Get the packages** (choose the one that applies):
   - **A `renv.lock` file exists in the project root** → in the RStudio Console run:
     ```r
     renv::restore()
     ```
     This downloads and installs the exact recorded versions. Say *yes* if it asks
     to proceed. First time can take ~20–40 min; after that it is instant.
   - **There is no `renv.lock` yet** (first-time setup on a brand-new project) → run:
     ```r
     source("setup.R")
     ```
     This installs everything (CRAN + Bioconductor + CCprofiler from GitHub) **and
     writes `renv.lock`**. Then commit that lockfile (see below).
3. **Only now** add your input files to `data/raw/` and knit the report in
   `analysis/` (open the `.Rmd`, press *Knit*), or run the helper scripts in
   `scripts/`.

---

## Handing this project over / transferring ownership

The lockfile **must be generated on a machine that already has the analysis
working** (so it captures the versions that actually work), and it cannot be
produced by anyone who doesn't have the packages installed. So, **before you
transfer the repo:**

```r
# in RStudio, with SEC.Rproj open, on the machine where the analysis runs:
source("setup.R")          # installs everything and writes renv.lock
```
then commit and push the environment files:
```
git add renv.lock .Rprofile renv/activate.R renv/settings.json
git commit -m "Lock package environment with renv"
git push
```
`renv/library/` (the actual packages) is intentionally **not** committed — it is
large and machine-specific; `.gitignore` already excludes it. `renv.lock` is the
only thing the next person needs.

Once `renv.lock` is in the repo, whoever clones it just does step 2A above
(`renv::restore()`) and gets an identical environment.

---

## If something goes wrong

- **A package won't build** (often `Rmpfr`, needs the GMP/MPFR system libraries;
  or Windows needs *Rtools*, macOS needs the *Xcode command-line tools*): see
  **`docs/SETUP_AND_GIT_WORKFLOW.md` → "System libraries"**.
- **CCprofiler**: this project uses the official `CCprofiler/CCprofiler`
  `differential` branch **pinned to an exact commit** in `setup.R`, plus the lab's
  fixes consolidated in **`R/ccprofiler_fixes.R`** (sourced on top of the package
  by the `.Rmd`). You do not need the old `AnnaPagotto/CCprofilerDiffAnna` fork.
- **New to R / Git / RStudio?** Full step-by-step guide:
  **`docs/SETUP_AND_GIT_WORKFLOW.md`**.
- **What is done and what is next:** **`HANDOVER_NOTES.md`**.

More detail on everything above is in **`README.md`**.
