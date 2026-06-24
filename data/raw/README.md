# data/raw — input files

Put the **raw input files** here. They are **not tracked by git** by default
(see `.gitignore`) because proteomics reports are large.

## Layout (multiple comparisons)

The `yeast_QTL` analysis runs several pairwise comparisons (mutant strains 31/83
vs parental strains BY/RM, plus RM vs BY). Each comparison has its **own
Spectronaut `.tsv`**, but they **share** the calibration table and the Complex
Portal table. So the layout is:

```
data/raw/
├── IA003_01_SEC_calibration_table.xlsx   <- shared MW calibration (all comparisons)
├── Complex_portal_559292.tsv             <- shared Complex Portal table (S. cerevisiae)
├── <...>_31_BY_<...>.tsv   ┐
├── <...>_31_RM_<...>.tsv   │  the five Spectronaut reports, placed FLAT here
├── <...>_83_BY_<...>.tsv   │  (the long original names are fine - no renaming)
├── <...>_83_RM_<...>.tsv   │
└── <...>_RM_BY_<...>.tsv   ┘
```

The `.Rmd` finds each comparison's file automatically by a short unique substring
of its name (the `pattern` in the registry), so you don't rename or type the long
Spectronaut filenames. Inputs are read with `here::here("data", "raw", <found file>)`.

## How it ties to the `.Rmd`

In the `# 1. SETUP` chunk there is a **comparison registry** and a selector:

```r
comparisons <- list(
  QTL_83_RM = list(pattern = "_83_RM_", ref = "RM", treat = "83"),
  ... )
comparison_id <- "QTL_83_RM"   # <- pick the comparison to run
```

To run a comparison:
1. Put its `.tsv` (flat) in `data/raw/`.
2. In the registry, set `pattern` to a short unique substring of that file's name
   (e.g. `"_83_RM_"`), and set `ref`/`treat` to the condition LABELS used inside
   the data (verify with `unique(design_matrix$Condition)` after the design-matrix
   chunk - they are NOT taken from the filename).
3. Set `comparison_id` to that entry and knit. Outputs land in
   `output/<comparison_id>/`.

Add a new comparison by adding a row to `comparisons` (and dropping its .tsv in `data/raw/`).

## Shared files

| File (default name)                        | What it is                          | Used as                |
|--------------------------------------------|-------------------------------------|------------------------|
| `IA003_01_SEC_calibration_table.xlsx`      | SEC molecular-weight calibration    | `calibration_location` |
| `Complex_portal_559292.tsv`                | Complex Portal table, *S. cerevisiae* | `complex_portal_file`  |

(On the `Ecoli_PCM` branch the Complex Portal file is `Complex_portal_83333.tsv`,
*E. coli* K-12; set via `complex_portal_file` in SETUP.)

## Including data in the handover

These files stay on your machine by default. To deliberately include a (small,
shareable) one, force-add it, e.g.:

```sh
git add -f "data/raw/IA003_01_SEC_calibration_table.xlsx"
```

For large files prefer sharing them separately or via [git-lfs](https://git-lfs.com/).
