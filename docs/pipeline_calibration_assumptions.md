# MW calibration in CCprofiler, PrInCE and PCprophet — and whether the claim holds

---

## 1. What each pipeline actually assumes

All three were read at source. Everything below is verified from code unless marked **UNVERIFIED**.

| | **CCprofiler** (github.com/CCprofiler/CCprofiler, master) | **PrInCE** (github.com/fosterlab/PrInCE, master, v1.7.1) | **PCprophet** (github.com/fossatiA/PCprophet, master) |
|---|---|---|---|
| **Does core inference require calibrated MW?** | **No.** Optional. Feature detection, correlation scoring, decoy FDR and q-values contain no MW term. `annotateComplexFeatures.R:151-178` has an explicit "neither monomer mw info nor calibration" branch that returns the identical detection columns. | **No — MW does not exist in the package at all.** No mass, no kDa, no calibration, no standards. `PrInCE()` accepts only `profiles` + `gold_standard`; there is no argument through which a calibration *could* enter. | **No.** `-cal` and `-mw_uniprot` both default to the string `'None'`; the calibration step is skipped in a bare `try/except`. The repo's own regression test (`run_test.py`) runs end-to-end with no calibration. |
| **Where MW enters** | `calibrateMW.R` fits it; `tracesMethods.R:120-124` (`annotateMolecularWeight`) is the **only** line that writes a fraction→MW value. Downstream: `findProteinFeatures.R:82-101` (`in_complex := apex_mw > 2*monomer_mw`), `annotateMassDistribution.R:12-30` (monomer/assembled split at `2*protein_mw`), `featureMethods.R:240-256` (`filterFeatures`, `min_monomer_distance_factor`), `annotateComplexFeatures.R:95-124` (`mw_diff`), plot axes. | Nowhere. The six classifier features are Pearson (raw/cleaned/p), Euclidean distance, co-peak (distance **in fractions**), co-apex (distance between Gaussian μ, **in fractions**). The x-axis is the raw fraction index throughout. | Post-processing only. `collapse.py:425-443` fits `log10(MW_Da) ~ fraction` with sklearn; used only in the non-default `-co CAL` redundancy-collapsing mode (`collapse_mincal`) and as the `"Estimated MW"` report column (`differential.py:614-641`). The RF feature vector (`generate_features_v2.py:402`) and the GO-based FDR contain no mass term. |
| **Calibration model** | `lm(log(MW_kDa) ~ fraction)` — **natural log**, despite roxygen saying `log2` and the diagnostic plot axis saying `log10`. Self-consistent (same base in fit and both inverses), only the docs are wrong. | None. | `LinearRegression()` on `log10(MW_Da)` vs fraction, then the slope sign is **overridden**: `calcfr = 10**(-abs(coef)*x + inter)/1000`. The fitted sign is discarded. |
| **Extrapolation handling** | **None whatsoever.** `calibrateMW` validates only file type, non-emptiness, equal vector lengths, and length ≥ 2 — nothing about the span. It never sees the traces object, so it cannot know the range it will be applied over. `annotateMolecularWeight` then applies the line to *every* fraction id unconditionally. No bound, no flag, no warning, no NA. A grep of all 37 files in `R/` and the vignette for `extrapolat|outside|valid range` returns one unrelated hit in `imputeMissingValues.R`. The only quality signal is an R² legend on an optional plot — which reports fit quality *inside* the standards and says nothing about outside. | Not applicable — no MW axis exists. The only bounds-checking is on Gaussian centres in fraction space (`filter_gaussians_center`, default TRUE). | **No code-level handling.** `calc_calibration` evaluates the fit at `xnew = list(range(1, 73))` unconditionally — no clipping, no NaN, no flag column, no warning. The only print is a global R². The warning exists **only as prose** in `PCprophet_instructions.md`: *"ONLY if a calibration curve covering all molecular weight range of the column is available. Extrapolation outside the standard leads to wrong molecular weight estimation."* |
| **Globularity assumption** | Never stated. Enters implicitly twice: (a) the stated validity condition is about the *device* ("standard SEC methods have a log-linear relationship"), not the analytes; (b) SEC-derived apparent MW is compared directly against UniProt **sequence** mass with no shape correction, at `2*monomer_mw` boundaries. An elongated or disordered monomer is silently pushed across the "assembled" line. | None instantiated. The only trace is prose in `R/PrInCE.R`: "separates … on the basis of their **diameter** or biochemical properties" — and nothing converts it. A non-globular protein is not mis-mapped; it elutes where it elutes and is compared to its partner at the same fraction index. | Never stated; the word "globular" does not appear. Implicit in the functional form and in `collapse_mincal`, which compares the **sum of UniProt sequence masses** of putative subunits against the calibration-derived apparent MW. |
| **Illustrative exposure** | Hard-coded SECexplorer standards `c(1398,699,300,150,44,17)` at fractions `c(19,29,37,46,54.5,61)` (`runSECexplorer.R:60-63`) — *every fraction below 19*, i.e. the entire large-complex end, is extrapolated. This block is a second re-implementation that bypasses `calibrateMW` entirely and so does not even get the length ≥ 2 check. | — | `cal.txt` is read from a hard-coded relative path (`differential.py:616`) while written to the cwd (`collapse.py:442`) — a stale file from an earlier run is silently picked up. |

**Gaps, flagged honestly:**
- **UNVERIFIED**: PrInCE on Bioconductor *devel* (the landing page returns 403 through the proxy). Conclusions are scoped to `fosterlab/PrInCE` master, v1.7.1.
- **UNVERIFIED**: PCprophet's packaged GUI binaries (`PCprophet_GUI/*.zip`) — no source, so I cannot see what defaults they ship. CLI/source path only.
- **Correction to a common belief**: `getMassAssemblyChange` **does not exist in CCprofiler master**. What you are using is `getMassAssemblyChange_aljazfix`, defined locally in `/home/user/SEC/analysis/DiffAnalysis_Ecoli_PMC.Rmd` (called at line 1190). The upstream analogue is `annotateMassDistribution`/`summarizeMassDistribution`. Say so in Methods, or a reviewer who greps the package will not find it.

**One-line summary of the comparison:** PrInCE never converts fractions to mass, so the question does not arise. PCprophet converts only for a report column and one optional collapsing mode. CCprofiler is the only one where a calibrated MW reaches anything you would call a result — and it is the only one of the three whose *documentation* does not warn you about extrapolation, while PCprophet's does.

---

## 2. Verdict on the claim

Your claim: *"a pipeline like CCprofiler … is doable and valid only for the ~50% of measured proteins that fall inside the calibration range, because the rest do not follow the assumption the pipeline makes."*

**Defensible:**
- The apparent MW **in kDa** is an interpolated measurement only between 17 and 670 kDa. Outside, it rests entirely on the untested assumption that the log-linear selectivity continues past the last standard. That is standard analytical doctrine (working range, Eurachem/ICH Q2), and CCprofiler gives you no flag, no warning and no NA.
- Extrapolated far enough, the model **falsifies itself**: 0.26 decades/fraction × ~30 fractions ≈ **8 orders of magnitude of mass**, against a real proteome span of ~2.7 decades. The earliest fractions get ≥10⁶ kDa; the latest get ~10² Da. Values that cannot exist.
- Inside the void fraction the map is **not injective** — everything above the exclusion limit co-elutes and elution position carries no size information at all. That is a genuine failure, and it is not fixable by better standards.

**Overreach, in three places:**

1. **"the rest do not follow the assumption."** Being outside the range is a property of *where you happened to put six standards on the column*, not of protein behaviour. In `globularity_check.R:336-337`, `beyond_calibration`/`below_calibration` are defined on `apparent_mw_kDa`, which is `mwmap[apex_fraction]` (line 302) — a strictly monotone function of the apex fraction alone. Two proteins in the same tube always get the same verdict whatever they are. GroEL, the 70S ribosome and pyruvate dehydrogenase are all out of range and all textbook globular; they violate nothing. Conversely a 40 kDa IDP eluting at 120 kDa apparent is squarely in range and completely misassigned — your own 19.1% in-range anomalous rate says so. **Range membership is neither necessary nor sufficient for the globularity assumption to hold.**

2. **"valid only for ~50%"** conflates the *kDa label* with the *inference*. In your own pipeline the two MW-dependent inferential steps do **not** depend on the extrapolation, for a reason worth understanding (see §3).

3. **The 53% is partly a choice you made.** `globularity_check.R:225-231` drops uridine (0.244 kDa) as the lower anchor — physically defensible, since it sits in the total-volume peak — but with uridine kept, the "calibrated interval" would cover nearly the whole axis and the same data would be "valid" for ~90% of observations. A validity criterion that moves when you change which vial in the kit you count is not a validity criterion. Say explicitly why uridine is excluded.

**The reformulation that survives scrutiny:**

> Within the standards' interval (F15.8–F21.9; 17–670 kDa) the apparent molecular weight is an interpolated, interval-scale quantity. Outside it, apparent MW degrades to a **censored ordinal bound** (`> 670 kDa-equivalent`, `< 17 kDa-equivalent`), and inside the void fraction it carries no size information at all because the fraction→size map is no longer injective. Any inference that treats the extrapolated value as a *number* — the apparent/monomer ratio, the frictional ratio f/f₀, a fold-change, a regression covariate — is undefined there. Inferences that depend only on the *ordering* of elution, or that evaluate the calibration at a protein's own monomer mass, are unaffected and use the full dataset.

That version concedes nothing you need, is expressed in ULOQ/LLOQ + censoring language every analytical reviewer already accepts, and — critically — it does not require you to defend the indefensible claim that out-of-range proteins are biologically different.

---

## 3. The circularity problem — read this before writing anything

This is the part that will get the paper into trouble, and your own code already contains the warning. `surface_hydrophobicity.R:569-578`:

> *"membership of `in_calibrated_range` is DEFINED by apparent MW leaving [min standard, max standard] — so for any protein whose expected mass sits inside that interval, 'out of range' FORCES |dev| past a fixed threshold. … 'The out-of-range proteins deviate most' is therefore largely a definition, not a finding."*

Concretely. For a protein of monomer mass *m*, in-range ⟺ apparent MW ∈ [17, 670] ⟺ `dev_log2 = log2(apparent/m)` lies in `[log2(17/m), log2(670/m)]`:

| monomer mass | in-range ⟺ dev_log2 in | i.e. apparent/expected in |
|---|---|---|
| 20 kDa | [−0.23, +5.07] | 0.85× – 33× |
| 35 kDa | [−1.04, +4.26] | 0.49× – 19× |
| 60 kDa | [−1.82, +3.48] | 0.28× – 11× |
| 100 kDa | [−2.56, +2.74] | 0.17× – 6.7× |

So "out-of-range proteins have the biggest deviation" is the *definition* of out-of-range restated. It would be true in a simulated dataset where every protein obeys one identical relationship and nothing is biologically different.

**Three consequences.**

**(a) The 80.9% / 44% contrast is partly manufactured by your own code.** `globularity_check.R:338`:

```r
d[beyond_calibration == TRUE | below_calibration == TRUE, globular_as_expected := FALSE]
```

Out-of-range proteins are **forced** non-globular *after* `dev_fractions` has already been computed for them, regardless of what it says. The out-of-range globularity rate is therefore 0% by fiat, not by measurement. It follows arithmetically that

> pooled % globular = (share in range) × (% globular in range) = 0.47 × 80.9 = **38.0%**

Your quoted 44% does not satisfy that identity (it implies 45%, not 53%, out of range). Most likely one figure is `median(S$pct_globular)` across metabolites (`globularity_check.R:585`) and the other is a pooled sum (lines 557-561), or void/below-calibration are counted differently. **Reconcile 53% / 80.9% / 44% against that identity before quoting any of them.** As stated they are over-determined and internally inconsistent, and a referee can spot that with a calculator.

**(b) The doubling of partial ρ is probably range restriction, not biology.** Selecting on a window of the *dependent* variable attenuates a correlation by

`ρ_restricted = ρ·u / √(1 − ρ²(1 − u²))`, where `u = sd(dev | in-range) / sd(dev | all)`.

An *exact* doubling on releasing the restriction requires **u ≈ 0.49** — and that threshold is essentially flat in ρ (0.498 at ρ=0.10, 0.468 at ρ=0.40). So the test is one line, on a file you already generate:

```r
J <- fread(".../surface_vs_elution_allproteins.csv")
sd(J[in_range == TRUE]$dev_log2) / sd(J$dev_log2)
```

- **u ≈ 0.5** → restriction explains the doubling completely. Nothing is left to attribute to biology, and the doubling must not be presented as evidence.
- **u > 0.5** → restriction under-explains it; you have a residual effect worth reporting.
- **u < 0.5** → restriction *over*-explains it; the association is genuinely **weaker** outside — the opposite of your reading.

Note also that range restriction only *attenuates* an existing correlation; it cannot manufacture one. So this does not say your hydrophobicity effect is fake. It says the **contrast between the two sets carries no information about calibration validity**.

**(c) Rank statistics are blind to the extrapolation anyway.** Because CCprofiler's map is a single monotone line evaluated everywhere, `dev_log2` is an affine function of the apex fraction: `dev_log2 = (a + b·f)/log10(2) − log2(m)`. Within a fixed mass, Spearman(dev_log2, −apex_fraction) = 1 by construction. A Spearman correlation therefore gives nearly the same answer whether you extrapolate log-linearly, extrapolate 3× steeper, hard-cap at 670 kDa, or drop the calibration entirely and use −apex_fraction. **Your hydrophobicity result is a calibration-free result wearing calibration's clothes.** That is good news for the result and fatal for using it as evidence that the calibration fails.

**The fix is already written.** `surface_hydrophobicity.R:579-596` implements `deviation = "residual"` — the residual of a robust within-dataset fit of log10(monomer mass) on apex fraction. No standards curve, no extrapolation, no definitional boundary. `surface_vs_elution_compare()` (line 751) runs all three sets side by side and prints exactly the right guidance at line 780: *"THIS is the number to quote when out-of-range proteins are included."* Run it and quote that column.

---

## 4. Making the argument airtight

**Do these four first — each closes a referee route, and each is hours of work.**

1. **Reconcile 53% / 80.9% / 44%** against `pooled = share_in_range × pct_globular_in_range`, and state whether each is pooled or a median across metabolites. Non-negotiable.
2. **Compute the honest out-of-range globularity rate**, bypassing the line-338 override. This matters more than you may realise: `dev_fractions` is computed as `|apex_fraction − .fraction_at_mw(mwmap, n × expected_mw_kDa)|` (`globularity_check.R:318`) — the reference position is evaluated at **n × the protein's own monomer mass**, n ≤ 4. For an E. coli proteome (median ~31 kDa) that is 31–124 kDa, i.e. an **interpolation**, even for a protein eluting in fraction 1. **Your globularity metric does not need the extrapolation; only the override does.** Tabulate `dev_fractions <= 1` among the out-of-range set and report that as the real number.
3. **Run `sd(dev_log2 | in-range) / sd(dev_log2 | all)`** and interpret against u ≈ 0.49 (§3b).
4. **Re-run the whole surface analysis with `deviation = "residual"`** and additionally with raw `−apex_fraction` as the outcome. If ρ is stable across all three, the calibration was never load-bearing and you can report the effect without any of this baggage.

**Then, to make the range argument itself positive rather than defensive:**

5. **f/f₀ ≥ 1 violation rate, in-range vs out-of-range.** A sphere has minimum friction at given mass, so f/f₀ < 1 is physically impossible. If it is ~0% in-range and substantial out-of-range, that is your single most persuasive number — physics falsifying the calibration precisely and only where it is extrapolated. You already build `ffo_vs_monomer` at `globularity_check.R:353`.
6. **Leave-one-out on the six standards.** Predict each from the other five. The two extreme standards are your only empirical test of mild extrapolation; their LOO error will be markedly worse than the interior ones. Zero new data, direct measurement of what extrapolation costs.
7. **Plot the 95% inverse-prediction band of the calibration across all 30 fractions**, with the calibrated interval shaded. The hyperbolic flare outside F15.8–F21.9 makes the entire argument in one panel.
8. **Rule out the "out-of-range proteins are just badly measured" objection** — compare peak SNR, summed intensity and peptides-per-protein in-range vs out-of-range. This is the objection most likely to be raised and the only one that could genuinely dent you.

**Scope the claim correctly in Methods.** In *your* pipeline (`/home/user/SEC/analysis/DiffAnalysis_Ecoli_PMC.Rmd`) MW enters at four places: `calibrateMW` (line 452), `annotateMolecularWeight` (826), `annotateMassDistribution` (1148) → `getMassAssemblyChange_aljazfix` (1190), and `CCprofiler::filterFeatures(..., min_monomer_distance_factor = 2)` (2229-2236).

Two things you should know about the last two, because they change what you can claim:

- **`filterFeatures` is in your inferential chain**, unlike in the CCprofiler vignette where the equivalent object is a dead end. Line 2242 scores the *filtered* object, and 2244 feeds it to `appendSecondaryComplexFeatures`. So this is not annotation in your hands. Do not describe it as such.
- **But it is not corrupted by extrapolation.** `apex_mw > 2 × max_monomer_mw` is, by strict monotonicity of the map, *exactly equivalent* to `apex_fraction ≤ fraction_at(2 × max_monomer_mw)`. The right-hand side is evaluated at twice a subunit's own mass — for E. coli, inside [17, 670] for any subunit between ~8.5 and ~335 kDa, i.e. essentially all of them. **The extrapolated value cancels out of the comparison.** The same argument applies to `annotateMassDistribution`, whose split point is `2 × protein_mw` (`annotateMassDistribution.R:19`).

So the honest scope statement is sharper and stronger than "half the data is invalid":

> Feature detection, co-elution scoring, decoy FDR, q-values, the assembly-state filter and the monomer/assembled split are unaffected by extrapolation, because they are evaluated at the protein's own monomer mass — an interpolation for essentially every E. coli protein. What extrapolation corrupts is the reported apparent MW in kDa and every quantity built from the apparent/expected ratio: f/f₀, dev_log2, and the globularity classification itself. All differential comparisons (ctrl vs metabolite, strain 31 vs 83) compare a protein to itself on the same column with the same calibration object, so any calibration error is a common factor that cancels exactly.

**Two more gaps to close.** `/home/user/SEC/output/` is empty on this checkout, so I could not verify 53% / 80.9% / 44% against data — they are taken from `/home/user/SEC/docs/methods_optimisation_notes.md:137-152`. And `surface_hydrophobicity.R` is blocked on AlphaFold structures, which sits on the critical path for the ρ argument. Note also that `surface_gravy`/`hydrophobic_sasa_frac` are the metrics to use, not sequence GRAVY (ρ² ≈ 3.8% is too weak to carry weight and is dominated by buried residues).

**Where this leaves you.** You have a real, publishable observation: a conventional globular-standard calibration constrains under half of a 30-fraction separation, the field extrapolates silently, and CCprofiler provides no warning while PCprophet's own documentation says not to. That stands on its own. What does *not* stand is the causal step — "the out-of-range proteins are the ones with the biggest deviation *because* they are outside the range" — which is true by definition and which your own script already flags as circular. Drop that step, run items 1–4, and the remainder is a clean result rather than an argument.
