# data/raw — input files

Put the **raw input files** for the analysis in this folder. They are **not
tracked by git** by default (see `.gitignore`) because proteomics reports are
often large and/or not yet shareable.

The analysis (`analysis/DiffAnalysis_Norm_QTL_31_83.Rmd`) expects these files
here, referenced via `here::here("data", "raw", ...)`:

| File (default name)                                          | What it is                                   | Set in .Rmd as        |
|-------------------------------------------------------------|----------------------------------------------|-----------------------|
| `20260414_QTL_83_31_IAA_Report_Angela Report (Normal).tsv`  | DIA quantitative report (peptide intensities) | `data_location`       |
| `IA003_01_SEC_calibration_table.xlsx`                        | SEC molecular-weight calibration table        | `calibration_location`|

If your file names differ, either rename the files to match, or edit the
`data_location` / `calibration_location` values in the `# 1. SETUP` chunk of the
`.Rmd`.

## Including data in the handover

By default these files stay on your machine only. If you want to hand the
project over **with** a (small, shareable) input file included, force-add it:

```sh
git add -f "data/raw/IA003_01_SEC_calibration_table.xlsx"
```

For large files, prefer sharing them separately (e.g. institutional storage)
or use [git-lfs](https://git-lfs.com/).
