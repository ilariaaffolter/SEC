# The statistics in the surface-hydrophobicity test, explained from scratch

Written for the handover. No statistical background assumed. All examples use the real numbers from
`surface_vs_elution(metabolite = "ATP")`.

---

## 1. What the test is actually asking

One sentence:

> **Do proteins with a greasier outer surface elute from the SEC column as though they were bigger than
> a single copy of themselves?**

If yes, the natural explanation is that greasy surfaces make proteins stick to each other, so those
proteins are travelling as pairs or complexes and therefore come off the column early.

Every number below is a way of answering that one question while ruling out boring explanations.

---

## 2. The two things being correlated

**On the x-axis — surface hydrophobicity.** Computed from the AlphaFold 3D structure. Several versions,
described in section 6.

**On the y-axis — the elution deviation:**

```
dev_log2 = log2( apparent MW / expected MW )
```

- *expected MW* = the mass of one copy of the protein, from its sequence (UniProt).
- *apparent MW* = the mass the SEC column thinks it has, read off the standards curve at the fraction
  where it peaks.

Because it's a `log2`, the scale is in doublings:

| dev_log2 | meaning |
|---|---|
| 0 | elutes exactly where one copy should |
| +1 | elutes as though **2×** its own mass |
| +2 | elutes as though **4×** |
| −1 | elutes as though **half** its mass (late — retained, or unusually compact) |

---

## 3. Spearman's rho (ρ) — the strength of the association

### The idea

Rank all proteins by surface hydrophobicity: greasiest = 1, next = 2, and so on. Separately, rank them
by elution deviation. **ρ asks whether the two rankings resemble each other.**

| ρ | meaning |
|---|---|
| **+1** | identical rankings — the greasiest protein deviates most, second-greasiest second-most, perfectly |
| **0** | the rankings are unrelated; knowing one tells you nothing about the other |
| **−1** | perfectly reversed |

### Why ranks, not the raw values

Two reasons, both of which matter here:

1. **Outliers can't dominate.** One protein with a freakishly large patch would drag a normal
   (Pearson) correlation around. Ranked, it is simply "the biggest one" — it counts once.
2. **The relationship doesn't have to be a straight line.** Spearman only needs "more of this goes with
   more of that", not that doubling one doubles the other. We have no reason to expect a straight line.

### Turning ρ into something interpretable: ρ²

**ρ² is the fraction of the variation explained.** Your best metric:

```
rho = 0.269   ->   rho^2 = 0.072   ->   7.2%
```

Read that as: **7.2% of why proteins differ in elution deviation is captured by surface
hydrophobicity. The other 92.8% is something else** — real assembly driven by other forces, shape,
measurement noise, the calibration itself.

Rough calibration for a proteome-wide correlation:

| \|ρ\| | ρ² | verdict |
|---|---|---|
| < 0.1 | < 1% | noise, whatever the p-value says |
| 0.1 – 0.2 | 1–4% | real if reproducible, but explains almost nothing |
| **0.2 – 0.3** | **4–9%** | **a genuine, modest effect — where your result sits** |
| 0.3 – 0.5 | 9–25% | strong for biology of this kind |
| > 0.5 | > 25% | rare outside a direct mechanical relationship |

So: a real signal, and far too weak to be the main story.

---

## 4. The p-value — and why yours are misleading

### What it is

**The p-value answers exactly one question:**

> If there were truly NO relationship at all, how often would random chance alone produce an
> association at least this strong?

p = 0.001 means "once in a thousand". Small p = hard to explain by luck.

### What it is NOT

- ❌ *not* the probability the finding is wrong
- ❌ *not* a measure of how strong or important the effect is
- ❌ *not* a measure of whether the effect is biologically meaningful

### The trap in your data

The p-value depends on **how strong** the effect is *and* **how many proteins** you have. With enough
proteins, even a trivial association becomes "highly significant":

| number of proteins | any \|ρ\| above this gives p < 0.05 |
|---|---|
| 50 | 0.278 |
| 200 | 0.139 |
| **1161 (your test)** | **0.058** |
| 2399 (the GRAVY test) | 0.040 |

With 1161 proteins, **ρ = 0.06 would already be "significant"** — and ρ = 0.06 is 0.36% of the
variance, i.e. nothing.

So your `p = 1.09e-20` is *not* evidence of a strong effect. It says the effect is **very reliably
not zero**. Whether it is *big enough to matter* is answered by ρ² alone.

> **The rule for this project: p tells you it's there. ρ² tells you whether to care.**

---

## 5. Confounding, and what "mass-partialling" means

### The problem

A third variable can create an association between two things that have no direct link.

Here that third variable is **protein size**:

- bigger proteins have **more surface area** (more of everything)
- bigger proteins are **more often in complexes**

So even if surface hydrophobicity had *nothing* to do with assembly, you would still see the two move
together — purely because both follow size.

### A worked demonstration

I generated 1161 fake proteins where size drives surface area, size drives elution deviation, and
**there is no other connection whatsoever**:

```
raw Spearman rho         = +0.722    <- looks like a spectacular finding
partial rho (mass held)  = +0.003    <- the finding evaporates, correctly
```

The raw ρ of +0.72 is **entirely an artefact of size**. This is why the raw column cannot be trusted.

### What partialling does

"Holding mass constant" means: **compare only what is left after size has been accounted for.**

Mechanically, four steps:

1. Rank the proteins by surface hydrophobicity, by elution deviation, and by mass.
2. Work out how much of the surface ranking is predicted by the mass ranking, and subtract it. What
   remains — the **residual** — is "greasier or less greasy *than expected for a protein of this size*".
3. Do the same for elution deviation: "elutes larger or smaller *than expected for this size*".
4. Correlate the two residuals.

An analogy: comparing the heights of a 6-year-old and a 15-year-old tells you mostly about age. Convert
each to "tall or short **for their age**", and now you're comparing something meaningful. Partialling
converts every protein to "greasy for its size" and "elutes big for its size".

### Your columns, precisely

| column | what is held constant | computed on |
|---|---|---|
| `rho_raw` | nothing | all 1161 |
| `rho_partial_mass` | monomer mass | all 1161 |
| `rho_partial_nomembrane` | **monomer mass** (still!) | the 963 non-membrane proteins |
| `rho_beyond_composition` | monomer mass **and** sequence GRAVY | all 1161 |

**`rho_partial_nomembrane` is not "rho without the mass correction".** It is the same mass-corrected
number, recomputed after dropping membrane proteins. Both control for mass; only the protein set differs.

### Why membrane proteins get their own column

Membrane proteins are hydrophobic (they sit in a greasy lipid bilayer) **and** often travel in large
particles. They could single-handedly manufacture the correlation without any interface being involved.
Dropping them tests that.

**Your result: the correlation went *up* when they were removed** (0.269 → 0.292). They were diluting
the signal, not creating it. That confound is excluded.

### Why `rho_beyond_composition` is the decisive one

Surface hydropathy and sequence hydropathy are themselves correlated — a protein made of greasy amino
acids tends to have a somewhat greasy surface. So "surface ρ = 0.269 beats sequence ρ = 0.155" does
**not** prove the surface adds anything; both could be measuring the same underlying composition.

Holding sequence GRAVY constant as well asks: **once I know the amino-acid composition, does looking at
the actual 3D surface tell me anything more?** Only if yes is this an interface result rather than a
composition result.

---

## 6. What each metric measures

### First: SASA, and what Å² means

**SASA = Solvent-Accessible Surface Area.** Imagine rolling a ball the size of a water molecule
(radius 1.4 Å) over the protein. The surface the ball can touch is accessible; anything it cannot reach
is buried inside. Measured in **square Ångströms (Å²)**, written `A2` in the column names because
column names can't hold the Å character.

**1 Å = 0.1 nanometres** — about the width of one atom. A typical 50 kDa protein has a total SASA of
roughly 18,000–20,000 Å².

Useful yardstick: **a real protein–protein interface buries roughly 600–1000 Å² on each side.** So a
hydrophobic patch of that order is "interface-sized"; one of 100 Å² is not.

### The metrics

| metric | what it is | why it's here |
|---|---|---|
| **`sasa_total`** | total accessible surface, Å² | the denominator for everything else; also a size proxy |
| **`sequence_gravy_here`** | plain average hydropathy over all residues (Kyte–Doolittle scale; + = greasy, − = water-loving) | **the comparison baseline.** This is what GRAVY measured — it counts buried and exposed residues equally |
| **`surface_gravy`** | the same average, but each residue **weighted by how much of it is actually exposed** | the hydropathy **of the surface**. A protein with a greasy core and a polar shell scores high on sequence GRAVY and low here — which is the normal state of a soluble protein, and precisely why GRAVY couldn't test your hypothesis |
| **`hydrophobic_sasa_frac`** | fraction of the exposed surface belonging to greasy residues (A, V, L, I, M, F, W, C, Y) | "how much of the outside is greasy" — a proportion, so it doesn't grow just because the protein is big |
| **`largest_patch_A2`** | area, in Å², of the biggest **single connected** greasy region on the surface | **the most biologically direct.** Proteins don't stick together through a high average — they stick through one contiguous sticky spot. Computed by finding all exposed greasy side-chain atoms, joining any within 5 Å into clumps, and taking the total area of the biggest clump |
| **`patch_frac`** | that same patch as a fraction of `sasa_total` | the size-normalised version — "what share of the surface is one big greasy spot" |

**Why both `largest_patch_A2` and `patch_frac`?** The first is absolute and grows with the protein; the
second is relative. A big protein can have a big patch that is nonetheless a small fraction of its
surface. They answer slightly different questions, and they behave differently under mass-partialling —
which is exactly what you saw (`largest_patch_A2` barely changed, `patch_frac` nearly halved).

### Why residues below pLDDT 70 are dropped

AlphaFold reports a confidence score (pLDDT, 0–100) for every residue. Below ~70 the model doesn't know
where the residue is — typically floppy tails with no fixed position. Including them would add
"surface" that isn't a real surface. They're excluded, and the discarded fraction is recorded per
protein.

---

## 7. Reading your actual result

| metric | raw | \|mass | \|mass, −membrane | ρ² |
|---|---|---|---|---|
| `hydrophobic_sasa_frac` | 0.341 | **0.269** | 0.292 | 7.2% |
| `surface_gravy` | 0.309 | 0.229 | 0.264 | 5.2% |
| `largest_patch_A2` | 0.170 | 0.179 | 0.155 | 3.2% |
| `patch_frac` | 0.280 | 0.157 | 0.160 | 2.5% |
| `sequence_gravy_here` | 0.197 | 0.155 | 0.153 | 2.4% |

**What is established:**

1. **The column-retention artefact is excluded.** Sticking to the resin would make proteins elute
   *late* and look *smaller* — a negative ρ. Everything here is positive.
2. **It is not membrane proteins.** Removing them strengthened the effect.
3. **Size alone doesn't explain it.** The correlations survive mass-partialling (0.341 → 0.269), unlike
   the fake example in section 5 where it collapsed to 0.003.
4. **Surface metrics beat the sequence metric** — 0.269 vs 0.155, roughly 3× the variance explained.

**What is not yet established:**

- Whether the surface adds anything **beyond composition** → `rho_beyond_composition`, now implemented.
- The effect is **weak in absolute terms**: 7.2% of variance. Real and reproducible, not explanatory.
- These are **five views of one signal**, not five independent confirmations. The metrics are computed
  from the same structures and correlate with each other.
- **Aggregation is indistinguishable from assembly here.** A protein that clumps non-physiologically
  also elutes early and also tends to have exposed greasy patches. Nothing in this analysis separates
  the two.

**One sentence you can defend:**

> Proteins with a greater proportion of exposed hydrophobic surface elute at higher apparent molecular
> weight than their monomer mass predicts (Spearman ρ = 0.27 at matched monomer mass, n = 1161,
> p = 1e-20), an association that strengthens when membrane-associated proteins are excluded and that
> runs opposite in sign to the hydrophobic column-retention artefact. The effect accounts for ~7% of the
> variance in elution deviation.
