# What the extrapolated apparent MW does and does not corrupt

Source-verified (every analysis read from the actual code; each finding re-checked by a second reader who
tried to refute it). The question: 53% of protein observations elute outside the six-standard calibrated
window, where apparent MW is an extrapolation. For WHICH analyses does that matter?

## The map

| Analysis | Results affected? | Detail (file:line) |
|---|---|---|
| **CCprofiler protein-level differential** | **No** | Feature detection (SW/peak-picking, fraction space), scoreFeatures decoy q-values (estimateQvalues.R), testDifferentialExpression and aggregatePeptideTests (testDifferentialExperssion.R — 0 mentions of MW) are all calibration-free. `protein_DiffExprProtein` and the hit list `pBHadj<0.05 & abs(medianLog2FC)>1` are byte-identical under any monotone fraction→MW map. The section DOES also emit a *reported* assembled-vs-monomeric percentage (featureMethods.R:112-119) that uses apparent MW — a summary annotation, not the hit list. |
| **CCprofiler complex-level differential** | **No (bulk)** | Same: complexes tested, peak boundaries, p/q/pBHadj, medianLog2FC all in fraction space or at subunits' own UniProt masses. The one MW gate, `filterFeatures(min_monomer_distance_factor=2)` (Rmd:2236), tests `apex_mw > 2*max_monomer_mw` — equivalent to a fraction comparison at twice the subunit's OWN mass, an interpolation for every subunit 8.5–335 kDa, so it **cancels** for the bulk. Genuinely extrapolation-dependent only for the edge sets (max detected subunit >335 kDa or <8.5 kDa); a small second-order effect on the q-value null because the filter reshapes the decoy pool. |
| **Assembly-state change** (`getMassAssemblyChange`) | **Values: no. Model: yes.** | The split is at `2*protein_mw` (annotateMassDistribution.R:19-30); the fraction→MW *values* never enter the arithmetic, only their ordering, so distorting the extrapolated kDa at fixed fit coefficients changes nothing. BUT the split point is a prediction of the fitted line, and `remove_lowest_MW` (Rmd:160,448) moves the slope 0.262→0.395 dec/fraction and the boundary by 1–2 fractions, changing the more/less-assembled lists. **Do not call this "immune".** Correct wording: it uses the full dataset and books the 53% out-of-range observations without reading their kDa values, but the hit set depends on how the calibration is modelled. |
| **CCF / EMD differential** (custom scripts) | **No** | best_lag, EMD, permutation p, BH q, hit calls, overlaps, power-ceiling diagnostics — all from normalised elution profiles. No statistic touches MW. Only opt-in `x_axis="mw"` *trace plots* use apparent MW (cosmetic). |
| **Globularity classification** | **Yes — this IS the affected layer** | apparent_mw_kDa, ratio, ffo_vs_monomer, ffo_vs_state, class, and the pies. Expected: it is the interpretive MW layer by design. |
| **Surface hydrophobicity** | **No** | Reported statistic is a rank-based partial Spearman on a quantity monotone in elution fraction; verified identical under log-linear, 3× steeper, hard-capped, and no calibration. |
| **PrInCE** | **No** | Molecular weight does not appear anywhere in the package; features are Pearson/Euclidean/co-peak/co-apex in fraction space (fosterlab/PrInCE). |
| **PCprophet** | **No (default)** | RF features and GO-FDR carry no mass term; calibration used only in the non-default `-co CAL` mode and the `Estimated MW` report column (fossatiA/PCprophet). Its docs warn against extrapolating. |

## `protein_DiffExprProtein` — exact contents

One row per protein (from `aggregatePeptideTests` → `aggregateTests(level="protein")`):
`feature_id`, `Npeptides`, `medianLog2FC` (local), `pVal`, `pBHadj` (**hit column**), `global_medianLog2FC`,
`global_pVal`, `global_pBHadj`, `FCpVal`, `global_FCpVal`; then UniProt join: `protein_name`, `mass`
(**sequence mass, NOT apparent MW**), `gene_names`, `go_f/go_p/go_c`, `ft_act_site`, `ft_binding`,
`cc_cofactor`. Hit = `pBHadj < 0.05 & abs(medianLog2FC) > 1`. No column derives from the calibration.

## Paper sentence (broad audience)

> Size-exclusion chromatography was calibrated with six globular protein standards, defining a
> molecular-weight window that spanned only ~47% of detected protein elution positions; outside it,
> apparent molecular weight is an extrapolation of the calibration curve rather than a measurement. This
> affects only the interpretive molecular-weight layer — the apparent mass assigned to each protein, the
> derived frictional ratio, and the globular-versus-anomalous classification — and not the differential
> analysis: detection of protein and complex elution features, their significance testing, and the
> identification of metabolite-responsive proteins depend only on co-elution profiles in fraction space
> and on each protein's own sequence-derived mass, and are therefore unchanged by the calibration. The
> same holds for the two most widely used co-fractionation pipelines — PrInCE derives interactions purely
> from co-elution profiles and never converts fractions to mass, and PCprophet uses molecular-weight
> calibration only for optional reporting. Globular-standard calibration is thus a limitation of
> molecular-weight annotation and shape inference, not of differential complex detection.

One caveat to carry into the methods/supplement, not the main sentence: the assembly-state analysis,
while it never reads the extrapolated kDa values, does depend on the fitted calibration line, so its hit
set is sensitive to the choice of standards (see `remove_lowest_MW`).
