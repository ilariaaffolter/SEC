# R/elution_aware — supporting scripts for the exploratory analyses

These are the lab's **"elution-aware" / protein-level** modules that the original
`.Rmd` loaded with `source(...)`. They support the **exploratory** sections of
the analysis, which now live in `analysis/experiments.Rmd` (not the main QTL
pipeline). They are sourced from there via `here::here("R", "elution_aware", "<file>")`.

You uploaded exactly the 5 scripts that are still needed; the other ~15 external
scripts the old `.Rmd` referenced (`CFMS_Hybrid_*`, `CFMS_PIQED_*`,
`CFMS_Comparison_*`, the various `_v2`/`best`/comparison modules, etc.) are **no
longer used** and were intentionally left out, per your note.

## The scripts

| File | Key function(s) | What it does |
|------|-----------------|--------------|
| `CFMS_Elution_Pipeline_P3_Height_PROTEIN_LEVEL_FIXED.R` | `run_protein_level_pipeline()` | Protein-level pipeline (no proteoform clustering): peptide→protein aggregation (`maxlfq`/`median`/`sum`/`top3/5/10`/`top10median`), FWHM peak detection on the reference condition, quantify all samples at reference peaks, limma + IBMT differential, flexible p-value adjustment (BH/qvalue/none/permutation), volcano + trace plots. Returns a `results` list. |
| `detect_proteoforms_elution_aware.R` | `detect_proteoforms_elution_aware()`, `plot_protein_divergence()` | Splits peptides of one protein into proteoforms by **intra-protein elution divergence** (log2 peptide/reference ratio over fractions). |
| `plot_protein_traces.R` | `plot_protein_traces()`, `plot_significant_proteins()`, `plot_high_fc_proteins()`, `plot_top_proteins()`, `plot_interesting_proteins()` | Protein chromatogram plots with individual peptide traces. |
| `plot_proteoform_peptide_traces.R` | `plot_proteoform_peptides()` (+ `plot_by_correlation/pvalue/foldchange/custom`) | Peptide-trace plots per proteoform, ranked by sibling-peptide correlation, p-value, or fold change. |
| `diagnose_aggregation.R` | `diagnose_aggregation_sparsity()` | Diagnostic: per-fraction valid-peptide counts (explains when different aggregation methods collapse to the same trace). |

## Important usage notes

- **They are parameterized.** Although the defaults/comments mention `ctrl`/`ATP`/`aSYN`
  (from the project they were written for), the condition names are **function
  arguments**. For this project, pass `ref_condition = "strain31"`,
  `treat_condition = "strain83"` (or whatever you compare).
- **Sample naming**: these scripts infer condition/replicate by splitting sample
  names on `_` (last part = replicate). So traces must be named like
  `strain31_1`, `strain31_2`, `strain83_1`, … for the auto-parsing to work.
- **Outputs**: each function writes plots/CSVs into an `output_dir` you pass —
  point these at `here::here("output", ...)` when you call them.

## Corrections applied during handover

- `detect_proteoforms_elution_aware.R`: the `tryCatch()` around the hierarchical
  clustering now assigns its **return value** to `div_clusters`, so the
  "all-in-one-cluster" fallback actually takes effect if `hclust`/`cutree` ever
  errors (previously the fallback was written only inside the handler's local
  scope). The normal (success) path is unchanged. This is the only code change;
  everything else is verbatim as you uploaded it.

> I could not execute R in this environment, so these scripts have had a **static**
> review only (no unquoted installs, no `setwd`/absolute paths, no separator typos
> were found). Please run them once locally to confirm functional correctness, and
> tell me if anything errors — I'll fix it.
