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
├── QTL_31_BY/   └── <the 31-vs-BY .tsv>
├── QTL_31_RM/   └── <the 31-vs-RM .tsv>
├── QTL_83_BY/   └── <the 83-vs-BY .tsv>
├── QTL_83_RM/   └── <the 83-vs-RM .tsv>
└── QTL_RM_BY/   └── <the RM-vs-BY .tsv>
```

The folder name is the **`comparison_id`** (matches the registry in the `.Rmd`
SETUP). The `.Rmd` reads the input as
`here::here("data", "raw", comparison_id, data_location)`.

## How it ties to the `.Rmd`

In the `# 1. SETUP` chunk there is a **comparison registry** and a selector:

```r
comparisons <- list(
  QTL_31_BY = list(tsv = "TODO_QTL_31_BY.tsv", ref = "BY", treat = "strain31"),
  ... )
comparison_id <- "QTL_31_BY"   # <- pick the comparison to run
```

To run a comparison:
1. Put its `.tsv` in `data/raw/<comparison_id>/`.
2. In the registry, set that entry's `tsv` to the exact filename, and confirm
   `ref`/`treat` (which strain is the baseline) and that those names match the
   sample names in the file.
3. Set `comparison_id` to that entry and knit. Outputs land in
   `output/<comparison_id>/`.

Add a new comparison by adding a row to `comparisons` and a folder in `data/raw/`.

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
