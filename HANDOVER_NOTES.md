# Handover notes & status

This file records **what has been set up**, **what was found in the original
`.Rmd`**, and **what still needs your input** to finish a clean, runnable
handover. It is the honest "state of the project" document.

_Last updated: 2026-06-19._

---

## 1. What is done (Phase 1 — project foundation)

- ✅ **Reproducible project skeleton**: RStudio project (`SEC.Rproj`), `renv`
  bootstrap (`setup.R`), `here`-based folder layout (`R/`, `analysis/`, `data/`,
  `output/`, `docs/`), and `.gitignore`.
- ✅ **CCprofiler dependency decided**: install the **official**
  `CCprofiler@differential`, **pinned** to commit `a6f4f4d`, via `renv`
  (`setup.R`). The old `AnnaPagotto/CCprofilerDiffAnna` fork is **dropped** — see
  §3.
- ✅ **Lab fixes consolidated** into `R/ccprofiler_fixes.R` (21 functions), with
  provenance and flagged known issues. The `.Rmd` will load these with one
  `source(here::here("R", "ccprofiler_fixes.R"))` instead of ~1,200 lines of
  inline definitions.
- ✅ **Docs**: `README.md` (overview/quickstart) and
  `docs/SETUP_AND_GIT_WORKFLOW.md` (beginner Git/RStudio/renv guide).
- ✅ **External scripts received** (2026-06-19): the 5 still-needed "elution-aware"
  modules are in `R/elution_aware/` with a README. The other ~15 scripts the old
  `.Rmd` referenced are no longer needed and were left out.
- ✅ **Clean main-pipeline split done**:
  - `analysis/DiffAnalysis_yeast_QTL.Rmd` — the cleaned QTL `strain31` vs
    `strain83` pipeline (data import → QC → traces → align/impute/normalize →
    filter → assembly/complex differential → save). Fixes are sourced from
    `R/ccprofiler_fixes.R`; packages use `library()` only; paths use `here()`
    (inputs from `data/raw/`, outputs to `output/`); the `protein_traces_list`
    gap is flagged inline (section 6). Won't knit end-to-end until you fill that
    gap and add the data — everything else is structurally clean.
  - `analysis/experiments.Rmd` — everything exploratory/other-project/one-off,
    grouped (A trypsin-ratio norm; B elution-aware; C PrInCE/Gaussian/one-off),
    `eval=FALSE` by default (reference, not run), with notes on which external
    scripts were retained. Obvious `, ,` syntax bug fixed.
  - The unmodified original is preserved as
    `analysis/ORIGINAL_DiffAnalysis_Norm_QTL_31_83_2026.Rmd`.

### Git status note (2026-06-19)
- Pushing now works; all work is on branch **`claude/clever-ride-3ab2zt`**.
- The repo's **`main`** branch currently contains a stray file
  `SEC_project_foundation.tar` — that's the archive committed as a single binary
  blob (the manual upload added the `.tar` itself, not its extracted files). It's
  harmless but should be removed. The clean way: merge
  `claude/clever-ride-3ab2zt` into `main` (which has the real, extracted project)
  and delete the stray `.tar`. Do NOT keep working from the manual upload.

---

## 2. The fork question, answered

`AnnaPagotto/CCprofilerDiffAnna@differential` was compared line-by-line to the
official `CCprofiler@differential`. The **only** differences were:
1. one added file, `R/PPlabFunctions.R` (all the lab fixes),
2. the old broken `normalizeByCyclicLoess` commented out (its fix is in
   PPlabFunctions.R), and
3. a package rename in `DESCRIPTION`.

So the fork = "official + one file of fixes". This project reproduces that by
installing the official package and sourcing the fixes from
`R/ccprofiler_fixes.R`. **Nothing is lost** by dropping the fork, and the fixes
are now visible/editable in your project.

## 3. The fixes file (`R/ccprofiler_fixes.R`)

- Verified on 2026-06-19 that your inline `.Rmd` fix functions are **identical**
  to the fork's `PPlabFunctions.R`, except two functions where your `.Rmd`
  version is *better* for this setup (it uses `CCprofiler:::.tracesListTest`,
  required when sourcing on top of the package). Your `.Rmd` versions were kept.
- Function bodies are **verbatim** from your `.Rmd`. The only mechanical change
  was **un-commenting** Sections G (cyclic-loess normalization) and H
  (1-replicate draft), which removes the leading `# ` and changes no code.
- **Flagged, not patched** (see the file header for details): a few fixed
  functions call CCprofiler/PrInCE *internal* functions without a namespace
  prefix (e.g. `make_initial_conditions`, `fit_curve`, `.intersect2`). These are
  pre-existing and may need a `CCprofiler:::`/`PrInCE:::` prefix when you first
  run them. They were **not** silently changed, per your "stay faithful" request.

---

## 4. What I found in the original `.Rmd` (important)

The original file (4,652 lines, 92 code chunks) is a working research document
that has grown to contain **three different things mixed together**:

1. **The QTL 31/83 pipeline** you described (conditions `strain31`/`strain83`).
   This is the part that should become the clean main analysis.
2. **A large block of exploratory code from a *different* project**
   (alpha-synuclein / PCM; conditions `ctrl`/`aSYN`/`PEP`/`GFP`/`intox`). These
   chunks `source()` ~20 external `.R` scripts (`ELUTION_AWARE/*.R`,
   `CFMS_*.R`, …) that are **not present** in the file or repo, and index data by
   `$ctrl_1`, `$PEP`, etc. — which don't exist in the QTL data.
3. **A few one-off data-processing tasks for yet other projects** (Marc's
   Spectronaut report, PCM, aSYN batch correction) at the very end, using
   absolute Windows paths like `C:/Users/ailaria/...`.

### Concrete issues that would stop it knitting (objective, mechanical)
- **Duplicate chunk label `setup`** (lines 53 and 4568) — knitr refuses to knit.
- **Unquoted `install.packages()`** (lines 232–235: `install.packages(parallel)`
  etc.) — errors. (Not needed anyway; `renv`/`setup.R` handle installation.)
- **`testDifferentialExpression_1repfix_chatgpt` is called (line ~3701) but was
  only defined inside a commented-out chunk** — now defined in
  `ccprofiler_fixes.R` (Section H), so this is resolved.
- **`, ,` empty arguments** in some Gaussian-fitting calls (lines ~3951–3955).
- **`setwd()` with an absolute Windows path** inside a chunk (line ~4572).

### Variables used but never defined in the file
- `protein_traces_list` — used throughout the assembly/complex sections (from
  ~line 3038 on) and saved to disk, but **never created** in this file. The step
  that would produce it (protein quantification → a `tracesList` of protein
  traces) appears to live inside the missing external scripts. **This is the main
  gap for the QTL pipeline** (see §5).
- `Sample_id_treat` — used in two differential loops (lines ~3487, ~3689), never
  assigned.
- `proteins_to_exclude` — used (line ~4244) but only ever defined as a
  commented-out line (123).

None of these were invented/guessed — they are flagged here for you to confirm.

---

## 5. Decisions taken / what is still open

**Resolved (2026-06-19):**
1. ✅ **External scripts** — received the 5 still-needed ones (now in
   `R/elution_aware/`); the rest are not needed and were left out.
2. ✅ **Other-project chunks** — confirmed: move the `ctrl`/`aSYN`/`PEP`
   exploratory chunks and the one-off tasks (Marc/PCM/batch) into
   `experiments.Rmd` (kept, tidied, documented, not deleted).

**Still open:**
- ⚠️ **`protein_traces_list`** — the assembly/complex sections need a CCprofiler
  protein-level `tracesList`, but the step that creates it isn't in the `.Rmd`,
  and the uploaded `R/elution_aware/` scripts produce a *different* (custom
  `results`) object, not a CCprofiler `tracesList`. In the cleaned main pipeline
  this step is inserted as a clearly-marked **FLAG** with the most likely intended
  call (`proteinQuantification_sibPepCorrFix.tracesList(pepTracesList_filtered, …)`,
  commented out) for you to confirm/adjust — it is **not** silently invented.

---

## 6. Phase 2 (later, as you planned)

Add a **cross-correlation–based** differential analysis section as an
alternative/complement to the statistical testing. Not started yet — we'll do it
once the pipeline above is settled.
