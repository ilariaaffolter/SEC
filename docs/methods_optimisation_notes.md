# SEC differential analysis — what was optimised, methods text, and take-homes

Branch `Ecoli_PCM`. Written as handover notes: the first section is an inventory,
the second is drop-in methods text, the third is what actually matters scientifically.

---

## 1. Points covered

**Differential testing**

- Permutation FDR with a **stratified (abundance-binned) null**, keeping the `min_intensity`
  noise filter — `scripts/ccf_differential_test.R`
- **EMD / 1-Wasserstein** statistic as an alternative to `|best_lag|`, with a three-way
  overlap against CCprofiler and CCF-FDR, shared UniProt IDs and trace plots —
  `scripts/emd_differential_test.R`
- **Power-ceiling diagnostic**: discrete statistics admit only a few distinct p-values, so
  `min_p = 1/(1+N)` sets a hard floor that BH may never clear. Reports `min_q`,
  `n_distinct_p`, `n_at_floor`, `n_needed_at_floor` — `scripts/diagnose_power_ceiling.R`
- UpSetR crash on all-empty sets; complex-feature parallelisation and caching

**Elution behaviour / globularity**

- `globularity_check()`: fraction-based tolerance window, calibration bounds at **both** ends,
  per-condition classification, three pie charts, all six standards overlaid for recovery
- Per-category annotation: GO enrichment, Complex Portal membership, filament heuristic,
  FoldIndex/IDR, pI and net charge, GRAVY — `scripts/globularity_category_annotation.R`
- The three "mundane cause" checks: low-end calibration coverage, silanol/pI retention,
  hydrophobic/GRAVY retention
- `shift_vs_charge.R` — is the differential shift confounded by charge? (No.)

**Biology-facing outputs**

- `hit_signature_analysis.R` — specific binding vs hydrotrope signature, with n, raw and BH p
- `validation_candidates.R` — 5–10 candidates per metabolite, photometric assays, central
  metabolism, heteromer-aware, DNA/RNA binders separated, known interactors excluded but 1–2
  retained as (+) positive controls, multi-metabolite enzymes bolded

**Cross-dataset validation**

- `secseq_compare.R` — Chihara SEC-seq, calibration-free deviation, orientation determined
  from the data, **partial Spearman controlling for monomer mass**
- `gradseq_ffo.R` — Hör Grad-seq, `gradseq_selftest()`, fraction→s models, physical floors,
  `gradseq_compare_models()`, calibration-free deviation route

**Bugs found and fixed (all were real)**

| # | Defect | Consequence |
|---|---|---|
| 1 | Deviation axes were different quantities (residual vs calibration ratio) | ρ = −0.022 meaningless |
| 2 | Raw correlation inflated by the shared −log10(mass) term | ρ = +0.971 was an artefact |
| 3 | `dev_cut` (a fold-change) applied to the mass-adjusted axis | reproducible set always empty |
| 4 | Calibration built on **peak** fraction, applied to **centre of mass** | f/f₀ 0.83 → 0.23 |
| 5 | Proteome fit regressed the noisy variable as predictor | slope attenuated, curve flat |
| 6 | Two deviations fitted on different populations against different mass columns | partial ρ only approximate |
| 7 | Mixed position measures inside `secseq_orientation` | spurious abort risk |
| 8 | Anchors with a non-finite fraction dropped silently | GroEL vanished from every figure |
| 9 | `getMassAssemblyChange(plot=T)` drawing one plot per protein (yeast) | pipeline stalled for hours |

**Still open**

- ⚠️ **HYDROPRO** — blocked on AlphaFold downloads through the institutional proxy;
  `hydropro_import_structures()` is the offline route
- ⚠️ **`QTL_83_Glc_vs_EtOH`** — awaiting the mixed Spectronaut export; needs `pattern` set and
  adding to `run_all_comparisons.R`
- ⚠️ **`QTL_83_RM_EtOH`** — speedups pushed, render not yet completed
- ⚠️ **Re-run `globularity_category_annotation()`** to pick up GRAVY and the pI tests on the
  corrected classification
- ⚠️ **Port the scripts to `main` and `yeast_aSYN`**
- ⚠️ **The absolute f/f₀ question is unresolved** — neither external dataset can settle it
  (see take-homes)

---

## 2. Methods — SCOPE NOTE FIRST

**This is a SUPPLEMENTARY methods block, not the main one.** It covers only the statistics
layer, the calibration/globularity assessment and the cross-dataset validation — i.e. what
was built or corrected in this round. The main Methods must still describe the core pipeline,
which is upstream of everything below:

| Stage | Implementation | Who supplies it |
|---|---|---|
| Culture, metabolite treatment, lysis, SEC run | — | **you** |
| LC-MS/MS acquisition, Spectronaut search settings | — | **you** |
| MW calibration and trace annotation | `calibrateMW`, `annotateMolecularWeight` | code |
| Protein feature finding | `findProteinFeatures`, `collapse_method = "apex_only"`, `perturb_cutoff = "5%"`; `corr_cutoff`, `window_size`, `rt_height`, `smoothing_length` chosen by grid search | code |
| Feature scoring / quantification | `scoreFeatures` (FDR 0.05; 0.15 variant), `extractFeatureVals`, `fillFeatureVals` | code |
| Differential expression | `testDifferentialExpression_beniFix` → `aggregatePeptideTests`; hit = `pBHadj < 0.05 & abs(medianLog2FC) > 1` | code |
| Assembly-state change | `getMassAssemblyChange_aljazfix` (beta regression per protein) | code |
| GO enrichment | hypergeometric (`phyper`, upper tail) with Benjamini–Hochberg | code |
| Complex-level analysis | `findComplexFeatures`, `corr_cutoff = 0.9`, `window_size = 7`, `collapse_method = "apex_network"`, target + decoy complex hypotheses | code |

The two paragraphs below slot in **after** that description.

---

Differential elution was tested by label permutation with an abundance-stratified null:
proteins were binned by median intensity and permuted within bins, so that a protein is
compared only against others of comparable abundance, and an intensity floor was applied
first to exclude noise-level features. Because `|best_lag|` is discrete and admits only a
few distinct p-values, an Earth Mover's Distance (1-Wasserstein) statistic was implemented in
parallel, and a power-ceiling diagnostic reports the smallest attainable p-value
(`1/(1+N)`), the number of proteins sitting at that floor, and the number that would have to
sit there for Benjamini–Hochberg to return any hit at a given FDR. Elution behaviour was
classified against the standards curve within a fraction-based tolerance window, with
proteins whose apex falls outside the calibrated interval labelled as such at both ends
rather than silently extrapolated. All hydrodynamic conversions were verified before use
against eight proteins whose s₂₀,w and Stokes radius were independently measured
(`gradseq_selftest()`): f/f₀ computed from s agreed with f/f₀ computed from R_s to a median
of 5.1%, and the Siegel–Monty combination recovered the known mass to within 3%.

Cross-dataset comparisons were made calibration-free. Each dataset was reduced to a
deviation — the residual of a robust fit of log₁₀(monomer mass) against elution position
*within* that dataset — so that no calibration is transferred between labs; both fits were
performed on the shared protein set against a single mass annotation, and both datasets were
summarised by the same position statistic (peak, matching `apex_fraction`), since a centre of
mass is displaced along the axis by the tail of the profile. Because both deviations contain
the term −log₁₀(monomer mass) by construction, agreement was quantified by the **partial**
Spearman correlation holding monomer mass constant; the raw correlation is inflated toward 1
whenever elution position predicts mass poorly, which is the very phenomenon under study.
The share of variance in monomer mass explained by position is reported for each dataset,
and a comparison in which either dataset falls below 5% is labelled uninterpretable, since
its mass-adjusted deviation is then close to noise. For the sedimentation data, f/f₀ ≥ 1 and
native mass ≥ monomer mass were used as hard falsification tests: any fraction→s model
placing a substantial share of the proteome below either floor is rejected outright, and
model selection is made on that criterion rather than on fit to the ribosomal anchors.

---

## 3. Take-home messages

- **The "55.8% anomalous" figure is dominated by extrapolation, not biology.** Restricted to
  the calibrated MW interval the proteome is **80.9% globular** (19.1% anomalous; 5318
  observations = 2932 monomer + 2386 oligomer). Quote the restricted number and state the
  denominator.

- **The calibration is the weak link, and it is quantifiable.** The six standards span
  F15.76 (thyroglobulin, 670 kDa) to F21.85 (myoglobin, 17 kDa) — about **6 fractions and
  1.6 decades** — while proteins elute over ~30 fractions and ~8 decades. **53% of
  observations fall outside the calibrated interval**, so roughly half of every
  globular-calibration-based claim in this field rests on extrapolation. This is the
  strongest form of the argument against globular-standard calibration: not that it is
  conceptually wrong, but that it covers under half the data.

- **f/f₀ ≥ 1 is a hard physical floor and makes an excellent falsification test.** A sphere
  has the least friction for a given mass, so any model placing a large share of the proteome
  below 1 is falsified regardless of how well it fits its anchors. This single criterion
  exposed the peak-vs-centre-of-mass bug, the extrapolation problem, and the unusable
  sedimentation calibration.

- **Beware shared-covariate artefacts.** Two quantities built from the same term will
  correlate at ~1 for that reason alone. A raw ρ of +0.971 fell to a partial ρ of +0.505
  once monomer mass was held constant — the latter is the real, quotable cross-lab result.
  Always partial out the shared term and report which number is which.

- **The Hör glycerol gradient cannot validate SEC at the protein level.** It was built to
  separate ribosomes; ordinary proteins pile into the top few fractions, migration position
  explains almost none of the variance in monomer mass, and no fraction→s calibration is
  recoverable by any route. The paper itself never converts fractions to Svedbergs — it
  annotates the axis with qualitative A260 landmarks and keeps every conclusion relative.
  **The SEC-seq comparison, not the gradient, is the cross-dataset evidence.**

- **Discrete statistics have a power ceiling that no amount of data fixes.** With `|best_lag|`
  the p-value floor is set by the permutation design, not by effect size, so BH can return
  zero hits even when real differences exist. Check `n_distinct_p` before concluding "no
  effect".

- **Metabolite-induced shifts are not a charge artefact.** Proteome-wide, no metabolite shows
  a shift-charge dependence (baseline ρ = 0.057, wrong-signed), and PEP's nucleotide
  enrichment (21% vs 21%) is a clean internal negative control against ATP's (38% vs 20%).

- **Nor a hydrophobic-retention artefact — but only the exclusion is claimable.** GRAVY vs
  log₂(apparent/expected) gives ρ = +0.196 (n = 2399). Column retention by hydrophobic
  interaction would make sticky proteins elute *late* and appear *smaller*, i.e. a **negative**
  correlation; the observed sign is the opposite, so that artefact is excluded. The positive
  sign itself should **not** be interpreted as assembly: GRAVY averages over buried core
  residues while assembly is driven by *surface* hydrophobicity, aggregation would produce the
  same sign as genuine complexes, and membrane-associated proteins are both hydrophobic and
  often in large particles. With ρ² = 3.8% and those confounds, report the exclusion only.

- **At n ≈ 2400, p-values stop being informative about effect size.** |ρ| > 0.04 already gives
  p < 0.05, so p = 4e-22 reports the sample size. Judge these tests on ρ² and on the *sign*
  predicted by the competing explanations, never on significance.
