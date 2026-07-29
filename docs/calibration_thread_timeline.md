# The calibration → globularity → surface → structure thread

A chronological record of what was asked, what was built, what came back, and what survives. Written for
handover: each step states the question it answers, so anyone picking this up knows why the code exists.

Cross-references to **PrInCE** and **PCprophet** come from `docs/pipeline_calibration_assumptions.md`,
where all three pipelines were read at source and each finding re-checked by a second reader.

---

## 1. Does the proteome elute where its monomer mass says it should?

**Built** `scripts/globularity_check.R` — classifies every control protein by comparing its apex fraction
against the position expected for a clean oligomer state (1×–4×), using a **fraction-based** tolerance
window rather than a fold-change one, because one fraction step is already 1.83× in mass on this gradient.

**Answer** Pooled over everything, ~44% behave as expected. That looked alarming — until the next step.

**Take-home** A tolerance expressed in fold-change is meaningless if it is narrower than one fraction
step. The code now warns when that happens.

---

## 2. Is globular-standard calibration defensible at all?

**Question asked** *"Is it true that CCprofiler and PCprophet rely on globular-standard calibration being
correct? I want evidence this is not good practice."*

**Built** Calibration coverage diagnostics inside `globularity_check.R`; all six standards overlaid on the
apparent-vs-expected figure to check recovery.

**Answer, and this is the durable one**

| | |
|---|---|
| standards span | F15.76 (thyroglobulin, 670 kDa) → F21.85 (myoglobin, 17 kDa) |
| that is | **~6 fractions, 1.6 decades of mass** |
| proteins elute over | **~30 fractions, ~8 decades** once the line is extrapolated |
| observations outside the standards' interval | **53%** |
| globular **within** the calibrated interval | **80.9%** |

**Take-home** The 44% figure was mostly extrapolation, not biology. Restricted to where the curve is an
interpolation, the proteome is 80.9% globular — which is roughly what the literature expects. **Quote the
restricted number and state the denominator.**

**vs PrInCE / PCprophet** PrInCE never converts fractions to mass — the question does not arise for it.
PCprophet does, but only for a report column and one non-default collapsing mode, and **its own
documentation warns**: *"Extrapolation outside the standard leads to wrong molecular weight estimation."*
CCprofiler applies the fitted line to every fraction unconditionally, with no bound, flag, NA or warning
anywhere in `R/` or the vignette. **CCprofiler is the only one of the three that neither needs the
calibration for its core inference nor warns you when it extrapolates.**

---

## 3. Is the anomalous elution just a column artefact?

Three mundane explanations, each tested by the **sign** the artefact would predict — which is the sharpest
thing these tests can do, because the artefact and the biology predict opposite directions.

| check | artefact predicts | observed | verdict |
|---|---|---|---|
| low-end calibration coverage | — | 53% outside | **real limitation** |
| silanol / pI retention | negative ρ | not negative | **excluded** |
| hydrophobic retention (GRAVY) | negative ρ (sticky → late → looks small) | **ρ = +0.196** | **excluded** |

**Take-home** Judge these on the **sign**, never on the p-value. At n = 2399, any \|ρ\| above 0.040 clears
p < 0.05, so p reports the sample size. ρ² = 3.8% here — real, weak, and only the exclusion is claimable.

---

## 4. Can an orthogonal technique give an *absolute* f/f₀?

**Built** `scripts/gradseq_ffo.R` against Hör et al. 2020 glycerol-gradient data: a hydrodynamics
self-test, three fraction→s models plus a proteome-fitted one, and hard physical floors.

**Answer NO — and the reason is worth keeping.** Migration position explains ~0% of the variance in
monomer mass in that dataset: the gradient was built to separate ribosomes, and ordinary proteins pile
into the top few fractions. No fraction→s calibration is recoverable by any route. The paper itself never
converts fractions to Svedbergs — it annotates the axis with qualitative A260 landmarks and keeps every
conclusion relative.

**Take-home** `f/f₀ ≥ 1` is a **hard physical floor** (a sphere has least friction for a given mass) and
makes an excellent falsification test. It caught three separate bugs on its own — a peak-vs-centre-of-mass
mismatch, an extrapolation problem, and an unusable calibration. Also verified along the way: the
hydrodynamic algebra reproduces published s₂₀,w and Stokes radii to a median of 5.1%, and Siegel–Monty
recovers known masses to within 3% (`gradseq_selftest()`).

---

## 5. Do two independent labs agree on which proteins elute anomalously?

**Built** `scripts/secseq_compare.R` against Chihara et al. 2023 SEC-seq — a calibration-free deviation
(residual of a within-dataset robust fit), orientation determined from the data rather than assumed.

**Answer YES.** Partial Spearman **ρ = +0.505** holding monomer mass constant. Two labs, different
columns, buffers and growth conditions.

**But first it was wrong twice, and both are instructive.** A raw ρ of **+0.971** was a shared-covariate
artefact: both deviations contain the term −log₁₀(monomer mass) by construction, so they correlate at ~1
for that reason alone. And the two axes had originally been *different quantities* — one a bounded
residual, the other a calibration ratio running to 10⁶.

**Take-home** **This is the cross-dataset evidence, not the gradient.** And whenever two quantities are
built from a shared term, partial out that term and report which number is which.

---

## 6. GRAVY cannot test the interface hypothesis. Can structure?

**Question asked** *"Can you implement the surface hydrophobicity / interaction interface with AlphaFold
so I can actually test the hypothesis GRAVY cannot?"*

**Why GRAVY cannot** It averages hydropathy over the **whole sequence**, and that average is dominated by
**buried core** residues. Assembly is driven by hydrophobicity that is **exposed**. A greasy core with a
polar shell is the normal state of a soluble protein.

**Built** `scripts/surface_hydrophobicity.R`

| metric | what it measures |
|---|---|
| `sasa_total` | solvent-accessible surface, Shrake–Rupley, self-tested against the analytic two-sphere solution (0.2–5.7% per-atom error at 92 points) |
| `surface_gravy` | SASA-**weighted** Kyte–Doolittle = the hydropathy of the **surface** |
| `hydrophobic_sasa_frac` | share of exposed area from apolar residues |
| `largest_patch_A2` | largest **contiguous** exposed apolar patch, in Å². Assembly needs a patch, not a high average — the sharpest form of the test. A real interface buries ~600–1000 Å² per side |

Three confounds handled rather than mentioned: **size** (partial ρ on log₁₀ mass), **disorder** (pLDDT < 70
dropped), **membrane association** (test repeated with them excluded).

**Structures** `alphafold_import_tar()` takes the whole E. coli proteome in one download (UP000000625,
taxid 83333) instead of hand-picking proteins — with mmCIF *and* PDB parsing, since v6 may ship
mmCIF only.

**Answer**

| metric | \|mass | \|mass, −membrane | ρ² |
|---|---|---|---|
| `hydrophobic_sasa_frac` | **0.269** | 0.292 | 7.2% |
| `surface_gravy` | 0.229 | 0.264 | 5.2% |
| `sequence_gravy_here` | 0.155 | 0.153 | 2.4% |

Membrane control **passed** — every effect *grew* when they were removed, so they were diluting the
signal, not causing it. Surface beats sequence by ~3× in variance.

**Still open** `rho_beyond_composition` — the metric with mass **and** sequence GRAVY held constant. Only
if that survives does the surface add anything the composition did not already say.

**Take-home** Aggregation and assembly remain **indistinguishable** here. A protein that clumps
non-physiologically also elutes early and also has exposed greasy patches.

---

## 7. Does adding the out-of-range 50% prove the calibration fails?

**Question asked** *"It basically doubles the rho — evidence that the pipeline is valid only for the 50%
inside the calibration range?"*

**Answer NO, on three independent grounds.**

**(a) Circularity.** `in_calibrated_range` is *defined* by apparent MW leaving [17, 670] kDa, and
`dev = log2(apparent/expected)`. So for any protein whose expected mass sits inside that interval, leaving
the range **forces** \|dev\| past a fixed threshold — a 35 kDa protein can only be called in-range if its
dev lies in [−1.04, +4.26]. "The out-of-range proteins deviate most" is the definition restated.

**(b) Range restriction.** On one **fixed** underlying relationship, restricting to the middle 50% of the
deviation range cut ρ from 0.274 to 0.092 — a **3× swing with nothing changing underneath**. A doubling is
comfortably inside that.

**(c) The result does not depend on the calibration at all.** Tested directly:

| calibration treatment | partial ρ |
|---|---|
| log-linear, extrapolated everywhere | +0.517 |
| extrapolated **3× steeper** | +0.521 |
| **hard-capped** at 670 kDa | +0.517 |
| **no calibration** — just −apex_fraction | +0.521 |

Spearman uses ranks, and `dev_log2` is a monotone function of apex fraction at fixed mass, so the
calibration contributes nothing to the ranks.

**Take-home, and it cuts both ways** The surface and globularity work is **not** invalidated by the
extrapolation — it is calibration-independent, which is stronger than it looked. But for exactly that
reason it **cannot** serve as evidence that the calibration fails.

**Built in response** `deviation = "residual"` in `surface_vs_elution()` — no standards curve, no
extrapolation, no definitional boundary; and `surface_vs_elution_compare()`, which runs all three sets and
says which column to believe.

---

## 8. What do the other pipelines actually assume?

| | CCprofiler | PrInCE | PCprophet |
|---|---|---|---|
| core inference needs calibrated MW | **No — optional** | **No — MW absent from the package entirely** | **No — `-cal` defaults to `'None'`** |
| where MW enters | `annotateMolecularWeight`; then `in_complex := apex_mw > 2*monomer_mw`, `filterFeatures`, mass distribution, plot axes | nowhere — all six classifier features are in **fraction space** | report column + one non-default collapsing mode |
| extrapolation handling | **none** | n/a | none in code, **but the docs warn** |

**The finding that reframes everything.** `filterFeatures(min_monomer_distance_factor = 2)` tests
`apex_mw > 2 × max_monomer_mw`. Because the map is strictly monotone, that is *exactly equivalent* to
comparing the apex fraction against the position of **twice the subunit's own mass** — an interpolation
for essentially every E. coli protein. **The extrapolated value cancels.** Same for the monomer/assembled
split at `2 × protein_mw`.

**Take-home — the defensible scope statement**

> Feature detection, co-elution scoring, decoy FDR, q-values, the assembly filter and the monomer/assembled
> split are **unaffected** by extrapolation, because they are evaluated at the protein's own monomer mass.
> What extrapolation corrupts is the reported apparent MW in kDa and every quantity built from the
> apparent/expected ratio: f/f₀, dev_log2, and the globularity classification itself. All differential
> comparisons put a protein against itself on the same column with the same calibration object, so any
> calibration error is a common factor that cancels exactly.

That is narrower than "valid for only 50% of proteins" — and it survives review.

**Note for Methods** `getMassAssemblyChange` is **not** on CCprofiler `master`. It is exported by the
`differential` and `DA_module` branches and defined in `R/annotateMassDistribution.R` there. Cite the
branch or a reviewer will not find it.

---

## 9. Open — running now

`scripts/reconcile_calibration_claims.R` settles four things before any number is quoted:

1. **Do 53% / 80.9% / 44% agree?** They are over-determined. `0.47 × 80.9 = 38.0%`, not 44.2%; for 44.2%
   you would need 45.4% outside, not 53%. Most likely a pooled-vs-median mix-up — the script reports both.
2. **The honest out-of-range globularity rate**, read from `dev_fractions`, which is measured at n × the
   protein's *own* mass and so survives the override that forces the reported rate to zero. Quote that
   contrast, not "80.9% vs 0%".
3. **`u = sd(dev|in-range)/sd(dev|all)` against 0.49.** At the threshold, range restriction explains the
   doubling completely.
4. **f/f₀ floor violations, in-range vs out-of-range.** **The only test in this whole argument with no
   circularity**, because f/f₀ ≥ 1 is imposed by physics rather than by the calibration. If it is ~0%
   inside and substantial outside, that is physics falsifying the calibration precisely where it is
   extrapolated — and it is a better argument than the one this thread started from.

---

## The five things worth carrying out of all this

- **Restricted to where the curve interpolates, the proteome is 80.9% globular.** The alarming pooled
  figure was extrapolation.
- **53% of observations fall outside the standards' interval.** A measured coverage fact that stands
  alone, needs no correlation, and is the strongest form of the critique.
- **The core inference of all three pipelines does not require calibrated MW.** PrInCE never uses it,
  PCprophet uses it for reporting, CCprofiler's MW-dependent thresholds evaluate at the protein's own mass
  and cancel the extrapolation. **Scope the claim to the MW-derived quantities, not to the pipelines.**
- **f/f₀ ≥ 1 is the falsification test to lead with.** Physics, not definition — the only part of this
  argument with no circularity in it.
- **Sign, then effect size, then p — in that order.** Every real conclusion in this thread came from a
  sign test or a physical floor. None came from a p-value.
