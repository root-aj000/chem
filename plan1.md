# plan1.md — Forensic Architecture Plan: CASMI 2026 Molecule ID From Mass Spectra

> A complete, algorithmic-first design. High-level picture first, then every mechanism
> taken apart to the algorithmic level: inputs, outputs, data structures, complexity,
> failure modes, and worked arithmetic.
>
> **Two tracks, deliberately separate:**
> - **TRACK A — LIBRARY SEARCH** (the *known* world): find the molecule among structures
>   and spectra that already exist, on disk or in a database.
> - **TRACK B — DATA-ONLY IDENTIFICATION** (the *unknown* world): reconstruct the molecule
>   from nothing but the spectra and physical/chemical knowledge. **No structure library,
>   no spectral library, no database lookup is permitted inside Track B.**
>
> Everything else exists to feed, check, or merge these two tracks.

---

## 0. How to read this document

| Section | What it gives you |
|---|---|
| 1 | The five problems the docs pose, and the seven design laws that fall out of them |
| 2 | The mental picture: one global architecture diagram, no detail |
| 3 | **Shared Inference Core (SIC)** — the machinery *both* tracks must have |
| 4 | **Track A** — library search, mechanism by mechanism |
| 5 | **Track B** — data-only identification, mechanism by mechanism |
| 6 | Fusion: how two ranked lists become one submission that cannot lose points |
| 7 | Problem → component map (every documented problem, and what kills it) |
| 8 | Research lineage: which field solved which analogous problem, and exactly what we take |
| 9 | Compute budget, artifacts, module layout (9 h, no internet, one notebook) |
| 10 | Validation protocol and honest experimentation |
| 11 | Build order and expected score progression |
| 12 | Risks, failure modes, kill-switches |
| 13 | Two full end-to-end worked examples (one per track) |

**A rule for reading:** anywhere you see a box labelled *WHY THIS BEATS A MODEL*, that is a
place where a neural network would be strictly worse than the arithmetic, and the reason is
stated. The plan's thesis: **an MS/MS spectrum is a noisy transmission of a molecular graph
through a partially-known physical channel (bond dissociation). Decoding it is a structured
inference problem, not a pattern-recognition problem. We build the decoder.**

---

## 1. The five hard problems (from the docs), and the laws they force

### 1.1 Problem 1 — A single spectrum is under-determined

Docs §12: four isomers of C8H8O2 give near-identical peak lists. The discriminating ions
(91, 95, 121) sit in the crowded 90–125 region where noise and co-isolation live. The data
sometimes does not contain the answer.

**Law 1 — Emit a posterior, not an answer.** The output of every stage is a *distribution*
over structures. MRR@25 is a ranked-list metric; the optimal submission for MRR@25 is the
list sorted by decreasing P(correct | all spectra). (Proof in one line: for any candidate at
rank r the expected contribution is P(correct)/r, so swapping any adjacent pair whose
probabilities are misordered strictly decreases expectation.) *Everything in this plan is
built to produce that one number: a calibrated posterior probability per candidate.*

### 1.2 Problem 2 — Molecules live in three difficulty classes, and you are not told which

Class 1 = reference spectra exist (search). Class 2 = structure known, no spectra (retrieve
then verify). Class 3 = not in any database (generate).

**Law 2 — Run both tracks, always, and merge by posterior.** Never bet on one track. A
candidate proposed by Track B must be able to beat a Track A candidate only if its posterior
is genuinely higher — never by being in a "protected slot", never by fiat.

### 1.3 Problem 3 — The training set is the wrong chemistry on the wrong instruments

`enveda-180` is 46 % of all spectra, is the only large timsTOF library, and is synthetic
drug-like screening compounds. The on-topic natural-product libraries are 10× smaller.
Random holdouts predict the leaderboard *backwards* (docs §7 / §13.4).

**Law 3 — The only valid validation is `ingest_lib == 'enveda-np-examples'` under three
masks** (visible / own-spectra-hidden / structure-removed-from-pool). Everything else is
diagnostics only. Track A additionally gets an *instrument-conditioned* similarity (a
timsTOF-only variant) and Track B is instrument-agnostic by construction, because it never
consults a spectrum library at all.

### 1.4 Problem 4 — The mass is exact but the *structure* is combinatorial

5 ppm at 350 Da = 0.00175 Da: the formula is essentially pinned. But a formula like
C15H10O7 has millions of isomers, and only a handful explain the peaks.

**Law 4 — Filter with arithmetic, then search with structure.** (a) Mass → formula by
integer arithmetic (no learning). (b) Formula + spectrum → *substructure inventory* by exact
set-cover. (c) Inventory → structures by *guided* graph search, never by exhaustive
enumeration. Exhaustive enumeration of C15H10O7 is ~10⁶–10⁷ isomers and is the single
biggest budget bomb in naive plans.

### 1.5 Problem 5 — Real spectra are noisy, multi-condition, and multi-spectrum

Median 20 peaks, median 3 spectra per molecule, 1–16 possible, collision energies merged
across [20,40,60,80], `base_peak_intensity` as a quality score, ~25 % of peaks unexplainable
by simple subtraction (rearrangements), NCE vs eV units in the training data.

**Law 5 — Every peak carries evidence *and* noise; model both.** Per-spectrum quality
weights from `base_peak_intensity`, per-peak emission probabilities, and evidence pooled as
*summed log-likelihoods* (Bayesian pooling of independent observations) with a hard
"contradiction penalty" so that one good spectrum can veto a candidate (docs §15.2's
min-rule, done properly).

### 1.6 The two laws that make it a system

**Law 6 — Never let a simulation outrank a measurement** (docs §5.3 library gate). A real
measured spectrum from a curated library beats an in-silico prediction, always, and the
system encodes this as a hard rule, not a learned weight.

**Law 7 — Remove structural ceilings** (docs §5.4). Any design that pins a class of
candidates to fixed ranks caps the score. In this plan the only fixed ranks are the ones the
posterior itself justifies.

---

## 2. The mental picture (global architecture)

Read this diagram top to bottom; it is the whole system. Everything below section 3 is one
of these boxes opened up.

```
        ┌──────────────────────────────────────────────────────────────────────┐
        │  RAW TEST DATA                                                       │
        │  ~400 molecules, ~1500 spectra, 1–16 spectra each, all timsTOF        │
        │  per spectrum: (m/z,intensity)[] + precursor_mz + adduct + CE list   │
        └───────────────────────────────┬──────────────────────────────────────┘
                                        │
        ╔═══════════════════════════════▼══════════════════════════════════════╗
        ║  SIC — SHARED INFERENCE CORE        (deterministic, no libraries)    ║
        ║                                                                      ║
        ║   S1  neutral mass (adduct inversion, cross-check across spectra)     ║
        ║   S2  molecular formula candidates   (integer mass-defect solve)      ║
        ║   S3  spectrum conditioning          (noise model, deisotope, top-N)   ║
        ║   S4  fragmentation tree             (Maximum Colorful Subtree, ILP)  ║
        ║   S5  neutral-loss lattice           (pairwise differences, ranked)   ║
        ║   S6  isotope consistency check      (element-count confirmation)     ║
        ║   S7  mixture-level joint inference  (co-occurring derivatives)       ║
        ║                                                                      ║
        ║   output per molecule: mass, formula set, loss lattice, FT, quality   ║
        ╚═══════════════════════════╤══════════════════════════════════════════╝
                                    │
             ┌──────────────────────┴───────────────────────┐
             ▼                                              ▼
┌────────────────────────────────────┐      ┌────────────────────────────────────┐
│ TRACK A — LIBRARY SEARCH           │      │ TRACK B — DATA-ONLY ID             │
│ (the known world)                  │      │ (the unknown world)                │
│                                    │      │                                    │
│ Pool: train ∪ COCONUT ∪ PubChem    │      │ FORBIDDEN: any structure or         │
│ ~1M structures, mass-indexed       │      │ spectral library, any database     │
│                                    │      │ ALLOWED: SIC output, chemistry      │
│ A1 mass-window retrieval           │      │ knowledge mined from data:          │
│ A2 LOSS-SPACE search (HSS/ΔMass)   │      │ the Substructure–Loss Lexicon,     │
│ A3 fragmentation-tree alignment    │      │ loss frequencies, NP class priors  │
│ A4 spectral entropy similarity     │      │                                    │
│ A5 motif-space substructure match  │      │ B1 lexicon lookup                  │
│ A6 kernel model channel (optional) │      │ B2 JIGSAW SOLVER: substructure     │
│                                    │      │    inventory by exact set-cover     │
│ A7 GBDT ranker + evidence algebra  │      │ B3 class prior from spectrum        │
│   (library gate, min-rule,         │      │ B4 SEED-AND-GROW graph assembler   │
│    popularity prior, calibration)  │      │    (beam search, never enumerate)  │
│                                    │      │ B5 exact forward verifier           │
│ output: ranked 25 + posterior      │      │ B6 log-likelihood pooling           │
│        (calibrated)                │      │ B7 diversity-aware top-25           │
└───────────────┬────────────────────┘      └───────────────┬────────────────────┘
                │                                          │
                └───────────────┬──────────────────────────┘
                                ▼
                ┌───────────────────────────────────────┐
                │ FUSION & SUBMISSION ENGINE            │
                │  sort by posterior · InChIKey14 dedup │
                │  hard rules (measurement > sim)       │
                │  class-aware + popularity fill to 25  │
                │  write submission.csv                 │
                └───────────────────────────────────────┘
```

**One-sentence summary of the flow:** *turn spectra into a formula and a loss lattice
(arithmetic); ask the known world in loss-space and tree-space (Track A); and independently
reconstruct the molecule from its pieces with a guided graph search and an exact forward
check (Track B); then let calibrated posteriors decide the order.*

---

## 2.1 What this plan changes relative to `plan.md` (and why it is better)

`plan.md` is a faithful transcription of the field's consensus: mass filter → 4 channels
(library / analog / MetFrag / neural fingerprint) → LightGBM → the five tricks → and for
Track B, SIRIUS fragmentation trees + MAYGEN enumeration + MetFrag-lite + MCES. That
architecture tops out at the current 0.44–0.47 cluster because it is what the whole
leaderboard already runs. This plan keeps everything in it that works and replaces the four
places where it leaves points on the table:

| # | `plan.md` | `plan1.md` | Why it is better |
|---|---|---|---|
| 1 | Analog channel = adduct-shifted **modified cosine over peaks** | **Loss-space search** (NIST Hybrid Similarity Search) + **fragmentation-tree alignment** | Peak-space similarity is dominated by shared housekeeping ions; loss-space similarity matches *what broke off*, which is what identifies a family. Tree alignment matches the *cascade shape*. Both are pure arithmetic and are the channels the field has not used. |
| 2 | Track B = enumerate all isomers of the formula, score each with a fragmenter | **Substructure inventory by exact set-cover, then seed-and-grow beam search with the spectral objective inside the loop, verified by an exact forward fragmenter** | Enumeration is exponential and budget-fatal (C15H10O7 ≈ 10⁶–10⁷ isomers); guided search is polynomial in beam width and directly optimises the spectrum objective. The NMR analogue of this design scores 67 % top-1 / 96 % top-10 where unguided methods collapse. |
| 3 | Multi-spectrum aggregation = **min** or mean | **Summed log-likelihoods** from a per-spectrum emission model + reliability weights + contradiction penalty | Min/mean are unprincipled; summed logs *is* Bayesian pooling of independent evidence and yields calibrated posteriors — which is exactly what the MRR-optimal sort (Law 1) needs. |
| 4 | Nothing beyond per-molecule inference | **Mixture-level joint inference** (ZODIAC-style consistency over the 400 co-occurring test molecules) + calibrated cross-track fusion | The test set is a real mixture: co-occurring derivatives share fragments and losses. Enforcing global consistency sharpens formula, substructure and class estimates for *every* molecule, and calibrated probabilities make Track A and Track B commensurable on one axis. |

Everything else in `plan.md` — the 8.5 ppm mass filter, the dedup, the always-25 rule, the
three-mask validation, the library gate, the popularity prior — is retained verbatim, because
it is correct, and is restated here at implementation depth.


---

## 3. The Shared Inference Core (SIC)

Both tracks consume these outputs. The SIC is **deterministic arithmetic and combinatorial
optimisation**: no gradients, no learned weights (the only learned quantities are *frequency
tables* mined by counting). It is also where most of the score-per-line-of-code comes from.

```
INPUT (per molecule)
   spectra[] = { mzs[], ints[], precursor_mz, adduct, base_peak_intensity, CE_ev[] }

S1  neutral mass M, adduct inversion, cross-spectrum agreement vote
S2  formula candidates F = { f1..fk } (mass-defect integer solve)
S3  conditioned spectrum per acquisition: cleaned peaks, quality weight q
S4  fragmentation tree T per spectrum (MWCS by ILP)
S5  loss lattice L per spectrum (pairwise differences, loss-prior ranked)
S6  isotope consistency score per formula candidate
S7  global pass: co-occurrence priors over the whole test mixture

OUTPUT (per molecule)
   M, formula set + scores, per-spectrum (peaks, q), loss lattice with prior scores,
   fragmentation tree with formula annotations, mixture-level consistency score
```

### 3.1 S1 — Neutral mass by adduct inversion

**What.** `M = precursor_mz ∓ passenger(adduct)`. Positive-mode adducts add mass (so we
subtract); negative-mode adducts remove mass (so we add back).

```
adduct              mode      passenger (Da)   rule
[M+H]+                +        1.00728          M = mz − p
[M+NH4]+              +       18.03382          M = mz − p
[M+H2O+H]+            +       19.01784          M = mz − p
[M−2H2O+H]+           +       37.02839          M = mz − p
[M+Na]+               +       22.98922          M = mz − p
[M+K]+                +       38.96316          M = mz − p
[M−H]−                −        1.00728          M = mz + p
[M−H2O−H]−            −       19.01784          M = mz + p
[M+CH2O2−H]−          −       44.99820          M = mz + p
[M+Cl]−               −       34.96885          M = mz + p
```

**Worked check (docs §5.3):** `precursor_mz = 301.0354`, adduct `[M−H]−` →
`M = 301.0354 + 1.00728 = 302.0427`.

**Three refinements that matter more than the formula itself:**

1. **Cross-spectrum vote.** All 1–16 spectra of one molecule must invert to the *same* M.
   Compute each; take the median; any spectrum disagreeing by > 10 ppm is flagged
   `MISLABELLED` and its evidence weight drops to 0.25 in all downstream pooling. This is
   free quality control that catches adduct-column errors (docs Part 2 §15.1).
2. **Adduct sanity.** If inversion gives M outside [150, 1159] Da (the test range from the
   dataset description), the adduct is wrong — try the other nine and keep the ones that land
   in range, ranked by (a) library frequency of that adduct for that mass region, (b) the
   resulting formula's fit (S2). This costs nothing and removes a whole class of hard failure.
3. **Charge.** The dataset promises z ∈ {+1,−1} almost always. We additionally test z = 2
   only when the neutral mass comes out at exactly half the expected range and the isotope
   pattern shows the 0.5-Da spacing; this is a 3-line check that fires on < 0.5 % of rows and
   is cheap insurance.

### 3.2 S2 — Molecular formula by mass-defect integer solve

**What.** Enumerate every chemically admissible formula within the precursor tolerance.

**Why not a learned formula predictor?** Docs Part 2 §5.3: at 5 ppm the formula is
*essentially always unique*, and the enumeration below is a few thousand integer operations.
A learned model (MIST-CF et al.) is needed only when accuracy degrades to tens of ppm.

**The algorithm (exact):**

```
elements  E = {C, H, N, O, S, P, F, Cl, Br, I}   (drop what the formula excludes)
masses    m(E) = exact monoisotopic masses
target    (M, σ)  with σ = 5 ppm × M / 1e6

for nN in 0..min(6, plausible):                   # natural products: few N
  for nO in 0..maxO:                             # maxO from mass bound
    for nS in 0..3:
      for halogens pattern in small set:         # ∏ over F,Cl,Br,I counts ≤ 6
        # two linear equations, two unknowns (nC, nH):
        #   12·nC + 1.007825·nH  =  M − base_term
        #   0·nC  + 0.007825·nH  =  defect − defect_term
        solve for nH then nC; keep if both are integers within tolerance
        apply hydrogen cap: nH ≤ 2·nC + 2 + nN + halogens
        apply RDBE:       rdbe = nC − (nH+X)/2 + nN/2 + 1  must be ≥ 0 and integral
        emit formula with score = −(|mass error in ppm|)
```

**The intuition in one picture:** every element has a *mass defect* (exact − nominal), and
defects add. Carbon's defect is 0, hydrogen's is +0.007825, oxygen's is −0.005085. So the
fractional part of M is a tiny integer combination of 3–5 numbers — like solving
`0.0427 = 10·(0.007825) + 7·(−0.005085)`. The two equations above are exactly that.

**Worked (docs Part 2 §5.2, target 302.0427, nN = 0):** the loop over nO = 0..10 yields a
single survivor, nO = 7 → nH = 10.000 → nC = 15.000, i.e. **C15H10O7**, verified by direct
summation 302.042655 ± 5 ppm. Output: a *set* of formula candidates with ppm-error scores,
almost always a single element.

**Output format:** `F = [ {formula: "C15H10O7", ppm: +1.2, rdbe: 11, score: 0.0}, ... ]`.
Track B consumes all candidates (not just the top one) with their scores — that is how
formula uncertainty is propagated honestly instead of being frozen at step one.

### 3.3 S3 — Spectrum conditioning (the noise model)

**What.** Turn a raw peak list into (a) a cleaned peak list and (b) a *quality weight*.

```
1. drop peaks with m/z > precursor + 2.0            (docs: a fragment cannot outweigh its precursor)
2. drop the top 1–2 peaks if they sit within
   ±0.02 of a 13C satellite of a stronger peak      (deisotoping)
3. keep top-N by intensity with N = 60,
   but never fewer than top-6 per 50-Da window       (the classic library-search window filter:
                                                     preserves the low-mass diagnostic ions that a
                                                     global top-N throws away — docs Part 1 §12.4)
4. intensity transform:  i' = sqrt(i)                (tames the heavy tail; docs curation note)
5. relative-intensity floor 0.002                     (noise cut)
```

**Quality weight.** From `base_peak_intensity` (the raw count — the only absolute-intensity
signal in the data). Absolute thresholds are fragile across sources, so use the *rank within
the test set*:

```
q = 0.2 + 0.8 · Φ⁻¹( rank(bpi) / n )          # squashed empirical CDF into [0.2, 1.0]
```

This is monotone, robust to the wildly different raw scales of the 12 training libraries, and
needs no calibration — and it is exactly what the multi-spectrum pooling in §6.6 multiplies
each spectrum's log-likelihood by.

**Collision energy.** Use `collision_energy_ev` only (docs are explicit). Bucket each
spectrum into `low` (< 35), `mid` (35–60), `high` (> 60), `merged` (multi-valued). The bucket
becomes a *feature* in Track A and a *conditioning variable* in the Track B emitter: peaks in
a low-energy bucket are far more likely to be the precursor / single losses, so the emission
model is bucket-dependent. This converts a documented data hazard into signal.

*WHY THIS BEATS A MODEL:* the cleaned peak set and q are derived from two documented,
physical facts (fragments cannot outweigh the precursor; low base-peak counts mean noise).
A learned cleaner needs labels for "noise" that do not exist, and will happily learn the
per-library hallmarks of the *training* libraries, which are the wrong instruments.

---

## 4. S4 — Fragmentation tree (the workhorse of both tracks)

**What.** A tree that explains the observed peaks as a cascade of neutral losses, computed
*from the spectrum alone*, with no knowledge of the structure (SIRIUS, Böcker & Rasche 2008;
Dührkop & Böcker 2016).

**The algorithm.**

```
Build the fragmentation graph G = (V, E):
  for each peak p (a "color"): enumerate all formulas f with
        |m(f) − p| ≤ tol(p),  f is a subformula of one of the S2 candidates
        (also allow the ion-variant: f ± H, since measured fragments are ions)
  node = (f, p);  color(f,p) = p
  edge (f,p) → (f',p') iff f' ≤ f componentwise, loss = f − f' , and m(loss) matches
        the difference p − p' within tolerance
  add a virtual root: for each formula f of the precursor peak, a zero-cost edge root → (f, precursor)

Find the maximum-weight colorful subtree:
  colorful = at most one node per color (each peak explained at most once)
  weight   = Σ_nodes log(1 + i_peak)  −  Σ_edges λ · cost(loss)
  cost(loss) = −log P_loss(loss)      # P_loss = a *counted* frequency table, §8.2
Solve exactly by ILP, or by the "speedy colorful subtree" beam heuristic
  (Dührkop et al. 2018) when the graph is large — the heuristic is ~100× faster
  and, crucially, *ranks* the true tree in a top position, which is what we need.
```

**What it gives us (both tracks use all three):**

1. **Per-peak formula annotations** → the isotope/loss checks become element-level.
2. **The loss cascade** `root → ... → leaf`, i.e. *which groups left, in which order.*
   Track B's jigsaw solver (B2) consumes exactly this.
3. **A no-database structural embedding** that CSI:FingerID/CANOPUS-style kernels and our
   own motif channel can consume.

**Complexity.** With a conditioned spectrum of ≤ 60 peaks, a few hundred formula candidates
per peak and ~10 formula candidates for the precursor, the ILP has ~10⁴ variables and solves
in well under a second on one core. For 1500 spectra this is the dominant SIC cost; budget
~25 min total on CPU (measured scale of SIRIUS is faster, but its scoring has more terms).

**Failure modes and the fix.**
- *Peaks above precursor* — already removed in S3.
- *Two formulas for one peak* — handled by the colorful constraint.
- *Rearranged losses (retro-Diels-Alder etc.)* — these are the documented ~25 % of peaks
  that simple subtraction cannot explain (docs Part 2 §7.5). We do **not** force them. They
  are marked `UNEXPLAINED` and become a *penalty term shared by all candidates*, so they
  cannot discriminate — which is the correct Bayesian treatment.

### 4.1 S5 — The neutral-loss lattice

**What.** All pairwise peak differences, labelled and ranked — the raw material for Track A's
loss-space search and Track B's jigsaw solver.

```
for each ordered pair (p_i, p_j) with i > j in intensity:
    d = mz_i − mz_j
    if d < 0.5: skip
    label d with the best-matching loss formula from a table of ~400 common losses
    record (d, formula, count, Σ intensity of the pair)
keep the top 200 by count × intensity
```

Two details that pay:
- **Chained losses** (docs Part 2 §7.3): record *paths* of 2–3 consecutive differences, not
  just pairs. A cascade `178.9986 → 151.0037 → 123.0088` losing CO twice is a far stronger
  signal than three independent CO losses, and it is exactly what the FT edge set encodes.
- **Loss plausibility prior** `P_loss(d)`, a *counted* table: for every loss formula, its
  frequency among all annotated peaks in the 2.5 M training spectra (mined offline, see
  §8.2). This is the only "learning" in the SIC and it is a histogram.

### 4.2 S6 — Isotope consistency (formula confirmation)

For each formula candidate, compute the theoretical M+1/M+2 intensities (binomial for C,
Poisson-corrected for N/O/S, exact for Cl/Br) and compare with any satellite peaks found in
the *MS2* spectrum (fragments do carry isotopes if the data is not deisotoped — determine
empirically per library on the training data; most timsTOF libraries are not deisotoped).

```
score_iso(f) = Σ_k  i(M+k) · log( p_k(theoretical for f) )
```

This is a soft score (the fragment-level isotope pattern is a weak signal) but it is *free*,
it is decisive for Cl (3:1) and Br (1:1) — the docs call this "the single highest-value 60
seconds of inspection" — and it costs 20 lines of code.

### 4.3 S7 — Mixture-level joint inference (the pass nobody runs)

**Fact:** the ~400 test molecules are one real mixture. The docs say the class distribution
was built to mimic "a real complex mixture", and ZODIAC (Ludwig et al. 2020, Nat Mach Intell)
showed the same signal in real data: compounds in a sample are biosynthetically related, so
their fragments, losses and formulas are *not independent*.

**Mechanism (offline-then-online, no learning):**

```
1. For each molecule, S4 gives an annotated loss set and formula candidates.
2. Build a graph over molecules: edge (a,b) if their loss sets share ≥ 2 rare losses
   (rare = P_loss below the 5th percentile) or their formulas differ by a common
   transformation (CH2, O, 2H, CO2, C6H10O5 ...).
3. Propagate: iterate a few rounds of
        log P(f | molecule) += α · Σ_neighbours w(a,b) · 1[neighbours' best formula agrees]
   (a mean-field / Gibbs-flavoured consistency pass; α small, ~3 rounds)
4. Output: sharpened formula posteriors + a "derivative family" label per molecule.
```

**Why this is worth real points here:** with 400 molecules, a candidate that is *globally
inconsistent* — e.g. a formula whose derivative relatives never appear anywhere in the
mixture while its neighbours all belong to a flavonoid-glycoside family — gets demoted, and
the correct family-consistent candidate gets promoted. It is pure joint inference over
counted statistics, and it is impossible to do with a per-molecule neural model.

---

## 5. TRACK A — LIBRARY SEARCH (the known world)

**The contract.** Track A answers one question: *is there something already measured or
already recorded that is this molecule?* It is the primary scoring engine because classes 1
and 2 are the majority of the test set and because, as the docs' phase logs show, adding
reference matching to a neural model was worth +0.19 MRR while the neural model itself was
worth +0.06.

```
TRACK A PIPELINE (per molecule)

  SIC output
     │
     ▼
  ┌──────────────────────────────────────────────────────────────┐
  │ A1  POOL: mass-window retrieval over ~1M structures          │
  │     (train 275k ∪ COCONUT ∪ PubChem subset), 8.5 ppm window  │
  │     → 20–400 candidates with mass scores                      │
  └───────────────────────────┬──────────────────────────────────┘
                              ▼
  ┌──────────────────────────────────────────────────────────────┐
  │ A2  LOSS-SPACE SEARCH  (Hybrid Similarity, ΔMass)            │
  │     "who else loses the same pieces, up to a constant shift?" │
  │     → analog evidence: analogy_tanimoto, ΔMass, family id     │
  ├──────────────────────────────────────────────────────────────┤
  │ A3  FRAGMENTATION-TREE ALIGNMENT                              │
  │     "whose cascade has the same shape?"                        │
  │     → tree_align (graph edit distance on trees, DP)           │
  ├──────────────────────────────────────────────────────────────┤
  │ A4  SPECTRAL ENTROPY SIMILARITY (exact structure)              │
  │     → lib_max per candidate, adduct-aligned                    │
  ├──────────────────────────────────────────────────────────────┤
  │ A5  MOTIF-SPACE SUBSTRUCTURE MATCH                            │
  │     "which recorded substructures does this spectrum declare?" │
  │     → motif_coverage vs the candidate's annotated substructures │
  └───────────────────────────┬──────────────────────────────────┘
                              ▼
  ┌──────────────────────────────────────────────────────────────┐
  │ A6  (optional) kernel model channel: spectrum → motif vector  │
  │     a linear/kernel model on FT features. ONE feature.        │
  ├──────────────────────────────────────────────────────────────┤
  │ A7  RANKER: GBDT (lambdarank, grouped by molecule_id)          │
  │     + evidence algebra: library gate, min-rule, prior          │
  │     + Platt calibration → posterior probability                │
  └───────────────────────────┬──────────────────────────────────┘
                              ▼
                 ranked 25 candidates, each with P(correct)
```

### 5.1 A1 — The candidate pool and its index

**Pool.** Every structure with a monoisotopic mass in [140, 1200] Da from:

| Source | Why | Size |
|---|---|---|
| `train.parquet` unique InChIKey14 | includes the test-instrument library (`enveda-180`, `enveda-np-examples`) | ~275 k |
| COCONUT (natural products) | the test chemistry; also the right prior for Track B | ~400 k |
| PubChem (mass-filtered dump) | class-2 recall; also the popularity-prior source | ~5–8 M |

**The index.** A single pass builds:

```
mass_index : dict nominal_mass → list of (structure_id, exact_mass, inchikey14)
formula_index : dict formula → list of structure_id
fp_index : structure_id → Morgan(2048) + a cheap "substructure set" (MAGMa-annotated motifs)
spec_index : inchikey14 → [(spectrum_id, peaks, adduct, ingest_lib, quality_weight)]
```

**Why 8.5 ppm, and the one improvement over the consensus.** The window is instrument-grade
(5 ppm) with a safety margin; it collapses a million structures to tens–hundreds. The
improvement: **index on nominal mass + defect band, not just nominal mass.** The mass defect
is a two-decimal fingerprint of the formula; a query with defect +0.0427 only ever matches
candidates whose *element set* can produce that defect. This is the same arithmetic as S2
applied to the database side, and it makes the pool smaller and, more importantly, makes the
S2 formula candidates *retrievable* rather than merely assumed.

**Class-1 emulation (validation only).** Because class labels are hidden, we approximate
class 1 as "the candidate set contains a structure whose training spectra are timsTOF and
whose best similarity exceeds 0.7" — which is also the library gate threshold. This gives us
a per-molecule *estimated* class posterior used only for diagnostics and fill strategy, never
for scoring.

### 5.2 A2 — Loss-space search: the analog channel, done in the right space

**The problem with the consensus version.** The field's analog channel is a *modified cosine
over fragment peaks* with an allowed constant mass offset. Peak-space similarity is dominated
by whatever all small molecules produce (the housekeeping ions, water/CO losses, the aromatic
series). Two unrelated molecules can share 80 % of their peaks.

**The better formulation — match in loss space (NIST Hybrid Similarity Search, Moorthy et
al., *Anal Chem* 2019; the same idea is the "ΔMass" parameter in NIST MS Search).**

```
Represent a spectrum by its neutral-loss vector:
    L(s) = { i_p·δ(mz_p) }  interpreted as masses *lost*: precursor − mz_p

Similarity:
    sim(q, s) = max over shifts δ of   cos( L(q), shift(L(s), δ) )
where the shift is a single global constant (the ΔMass), i.e.
    sim(q, s) = max_δ  Σ_p i_p · i_{p−δ}  /  (||L(q)|| · ||L(s,δ)||)
```

**Why this is the right space for this competition.**
- It is **invariant to the identity of the leaving group**: a query that is "candidate X plus
  a methyl" loses CH₃ where X loses H — a shift of exactly 14.0157 makes the loss vectors
  identical. That is precisely the analog relationship we want to detect, and it is why the
  docs call the analogous channel the single biggest idea (+0.19 MRR).
- It matches **what the structure is made of**, not what the instrument always produces.
- It is **adduct-agnostic**: comparing in loss space mostly cancels the adduct passenger
  (the docs' benzoic-acid example: [M−H]⁻ and [M+H]⁺ of the same molecule have *different*
  fragments but the same loss vocabulary).

**Implementation (the part people get wrong):**

```
1. OFFLINE: for every training spectrum of every pool structure with known formula,
   compute its loss vector L(s) (precursor-anchored, deisotoped, top-40).
   Hash it into a quantised loss index: bucket = round(loss / 0.02).
2. ONLINE per query spectrum:
   a. take the query's top 12 peaks by intensity.
   b. candidate shifts δ = { mz_lib − mz_q } for every (lib peak in the top-K of a
      retrieved set, query peak) pair — collected from a first pass over the
      adduct-and-mass-compatible spectra only.
   c. for each of the top ~40 shifts by vote count, compute the loss-space cosine.
   d. record sim, δ, and the owning structures.
3. Aggregate per candidate structure:  max over its spectra, then max over shifts.
```

**Cost.** The inverted index keeps the candidate set to a few thousand spectra per query;
each cosine is O(40). ~1500 queries × 40 shifts × a few thousand = fine on one core in
minutes.

**Two variants, because instrument match dominates** (docs §7.2):
- `AFIX`: restrict the index to `ingest_lib ∈ {enveda-180, enveda-np-examples}` (timsTOF) —
  matches the test instrument exactly.
- `BIO`: restrict to the natural-product-family libraries (`riken`, `gnps`, `massbank`,
  `mona`) — matches the test chemistry exactly.
Both are computed; both enter the ranker as separate features, and the GBDT learns how much
to trust each. (The published phase log shows blending the two beats either alone.)

### 5.3 A3 — Fragmentation-tree alignment (the channel the field does not use)

**Idea (Rasche et al., "Identifying the Unknowns by Aligning Fragmentation Trees", Anal
Chem 2012).** Two molecules in the same family fragment through *the same cascade shape*, even
when masses differ. Compare the trees, not the peaks.

```
OFFLINE: for every training spectrum, compute its fragmentation tree with the S4 algorithm
         (this is a one-time cost over 2.5 M spectra — see §10 for the budget) and store
         a compact signature:
             tree_sig = sorted multiset of (loss_formula, depth, intensity_share)
ONLINE:  compute the query's tree; compare to candidate structures' tree signatures with
         a tree-alignment DP:

    align(T_q, T_s) = max over monotone subtree mappings of
        Σ matched edges  [ 1[loss_formula_q == loss_formula_s] · w_match
                         + 1[depth_q == depth_s] · w_depth
                         + 1[leaf_intensity_share_q ≈ leaf_intensity_share_s] · w_int ]
        − gap penalties
    (O(n·m) DP on trees; n,m ≤ 60)
```

**Score produced:** `tree_align ∈ [0,1]` per candidate structure, plus the best-matching
*transformation set* — the multiset of loss-formula differences that reconciles the two
trees, which is literally a proposed structural edit (e.g. "add a methoxy at depth 2"). That
edit becomes a *feature* and a candidate-generation operator for the fill stage.

*WHY THIS BEATS A MODEL:* tree alignment is an exact DP over an object we compute ourselves.
It is instrument-tolerant (loss *formulas* rather than raw m/z), it generalises to mass
shifts for free, and it explains itself: you can print the aligned cascade for any decision.

### 5.4 A4 — Spectral entropy similarity (the exact-match channel)

For candidates whose structure has reference spectra, the match should be near-deterministic
for class 1, so this channel must be *accurate*, not fancy.

```
for each candidate c:
    best = max over c's training spectra s:
        sim = spectral_entropy_similarity(condensed(q), condensed(s))   # Böcker 2012
        # adduct alignment first: convert both to neutral-loss space / same anchor
        # weight s by its quality (library curation level, bpi) and by instrument match
    lib_max(c) = best
```

Specifics that matter:
- **Spectral entropy, not cosine** — it compares intensity *distributions* and is robust to
  the mismatched peak counts that different instruments and collision energies produce.
- **Weight by instrument.** A timsTOF reference for the same structure outranks an Orbitrap
  one; this is a fixed multiplicative weight, not a learned one.
- **Weight by library curation** (MassBank > community libraries), again fixed a priori.

**Hard rule (the library gate, docs §5.3):** if `lib_max ≥ 0.7`, the candidate is treated as
a *measurement*, all simulated evidence (A3's tree simulator, Track B's forward verifier) is
zeroed for that candidate, and no candidate with a lower `lib_max` may outrank it unless its
posterior from *measured* channels alone exceeds it. A real measurement never loses to a
simulation.

### 5.5 A5 — Motif-space substructure matching (a soft, library-wide channel)

**What.** Mine a *substructure vocabulary* and score candidates by which vocabulary entries
their structure can produce, against the motifs the spectrum declares.

This shares machinery with Track B's lexicon (§6.1) but with an important boundary: Track A
uses the lexicon *plus* full structures; Track B uses the lexicon *only*.

```
OFFLINE (all pool structures):
    for each structure s: enumerate its MAGMa-style fragment substructures
        (all connected subgraphs from ≤ 3 bond cleavages + H2O/NH3-type mobile losses)
    →  motif_set(s) : a set of (formula, mass, canonical substructure id)

ONLINE (query):
    decompose the query's loss lattice + FT annotations into lexicon motif loadings
    (dot product of the loss/peak vector with each motif's signature)
    →  declared_motifs with weights

SCORE per candidate:  motif_score(c) = Σ over declared motifs
        weight · 1[motif ∈ motif_set(c)]  · (1 − redundancy with motifs every structure has)
```

This channel is *not a fingerprint model*: it is dictionary matching over substructures, and
it is the one that carries the answer for class 2 (structure known, no spectra) where both
A2 and A4 are weak: the spectrum declares "hexose + flavonoid aglycone", and any candidate
that cannot produce a hexose and a flavonoid scores ~0.

### 5.6 A6 — The optional kernel model channel

**Where the neural/learned model is allowed in Track A, and where it is not.** The docs are
clear that a fingerprint model is worth having as *one feature*, and equally clear that it
must never be the decider. To honour "algorithmic, not pure ML", we make this channel
explicitly optional and give it the smallest possible footprint:

```
features per spectrum: fragmentation-tree loss multiset (binned), peak-count profile,
                       intensity distribution moments, bucket features (energy, adduct, instrument)
model:    kernel machine / linear model over these features → motif-presence logits
          (CSI:FingerID's construction; a small transformer is a drop-in replacement
           if GPU time is spare — it changes nothing else in the architecture)
feature:  fp_agreement(c) = ⟨motif_logits, motif_set(c)⟩
```

Rule: `fp_agreement` enters the GBDT as one column among ~40, and is zeroed by the library
gate whenever `lib_max ≥ 0.7`.

### 5.7 A7 — The ranker and the evidence algebra

**Features (per molecule × candidate), all cheap and heterogeneous:**

| group | features |
|---|---|
| mass | ppm error, defect error, formula S2 score, isotope score, rank in S2 list |
| A4 | lib_max, lib_2nd, best instrument-match lib score, n reference spectra |
| A2 | loss-cosine best, loss-cosine@best δ, δ value, ΔMass plausibility, n neighbours, Tanimoto of retrieved neighbours to candidate |
| A3 | tree_align, matched-loss fraction, depth agreement |
| A5 | motif_coverage, motif_redundancy, declared-motif count |
| A6 | fp_agreement |
| candidate properties | mass, heavy-atom count, ring count, H/C ratio, "famous-ness" (log PubChem substance count, log PubMed citations), presence in COCONUT |
| query properties | n spectra, median quality q, energy buckets, adduct |

**Model.** LightGBM with `objective=lambdarank`, grouped by `molecule_id`, ~40 features,
trained on the three-mask validation setup (§10). Reasons the docs give, restated: mixed
feature types, tiny per-group counts (tens of candidates), no normalisation assumptions,
free feature importance, ships as a text file.

**The evidence algebra — five rules that are *hard-coded*, not learned:**

1. **Library gate.** `lib_max ≥ 0.7` ⇒ zero all simulated features; the candidate competes
   only on measured evidence.
2. **Popularity prior.** `score += μ · (log(1+n_pubchem_substances) + log(1+n_pubmed))`
   with μ ≈ 0.25, applied only to candidates from a structure database. This is a legitimate
   prior over the label distribution (docs §5.2): a well-published natural product is far
   more likely to be in a practitioner-built test set than an obscure congener.
3. **Min-rule across spectra.** For every channel, a candidate's molecule-level score is the
   **minimum** over its spectra (or the 0.25 quantile), not the mean — one spectrum that
   flatly contradicts a candidate must be able to veto it (docs Part 2 §15.2). The
   log-likelihood pooling of Track B (§6.6) handles this with more nuance; Track A uses the
   simple min because its scores are not probabilities.
4. **No fixed slots.** Any candidate may occupy rank 1 if the algebra says so. (Docs §5.4:
   pinning outside-pool candidates to slots 4, 8, 12… cost +0.011 to remove.)
5. **Always 25.** Fill by posterior order; never submit fewer than 25.

**Calibration.** A Platt/isotone map on the GBDT output, fitted on the validation masks,
converts scores to `P(correct)`. This matters for fusion: only calibrated probabilities make
Track A and Track B comparable on one axis (Law 1).

**Track A output:** a list of ≤ 40 candidates per molecule with `P(correct)`, provenance
tags (which channel ranked it), and the library-gate flag.

---

## 6. TRACK B — DATA-ONLY IDENTIFICATION (the unknown world)

### 6.0 The philosophy, stated once, precisely

Track B may not look anything up. It must answer: *given only these spectra, what molecules
could exist?* The honest framing is the noisy-channel view:

```
    the molecule is a codeword written in the alphabet of atoms and bonds;
    the mass spectrometer applies a physical channel — collision-induced dissociation —
    which emits a random subset of the codeword's connected subgraphs plus rearranged
    losses; identifying the molecule is MAP decoding of that codeword:

        M̂ = argmax_mol  log P(spectra | mol)  +  log P(mol)
                         └── forward model, computable ──┘   └── chemistry prior ──┘
```

Everything in Track B is one of two halves of that equation:

| Half | Mechanism | Nature of the model |
|---|---|---|
| `log P(spectra \| mol)` | **B5** exact forward fragmenter + peak-matching assignment | combinatorial + counted emission probabilities |
| `log P(mol)` | **B3** class prior + **B2** substructure inventory plausibility | counted frequencies of substructures in natural products |
| the search | **B4** seed-and-grow beam search | combinatorial optimisation |

And the lexicon (**B1**) is the bridge: it is the only allowed source of "chemical knowledge",
and it is mined from data by *counting*, which is why it does not violate the "no library"
rule — it contains motifs of ≤ 12 heavy atoms, never whole molecules.

```
TRACK B PIPELINE (per molecule)

  SIC output (mass, formulas, loss lattice, FT, q)
        │
        ▼
  ┌─────────────────────────────────────────────────────────────┐
  │ B1  LEXICON LOOKUP: annotate every loss/peak with candidate  │
  │     substructures (multi-annotation, with confidences)        │
  ├─────────────────────────────────────────────────────────────┤
  │ B2  JIGSAW SOLVER: exact set-cover/ILP over the lexicon        │
  │     → substructure inventory Σ x_m·formula(m) = formula(M)     │
  │     → explains the loss lattice; mass must balance exactly     │
  ├─────────────────────────────────────────────────────────────┤
  │ B3  CLASS PRIOR: inventory → NP class distribution             │
  │     → hard "good list"/"bad list" of allowed motifs            │
  ├─────────────────────────────────────────────────────────────┤
  │ B4  SEED-AND-GROW ASSEMBLER: beam search over molecular graphs │
  │     score = substructure consistency + class prior +           │
  │            chemistry plausibility + partial loss coverage      │
  │     → candidate structures (NOT an enumeration)                │
  ├─────────────────────────────────────────────────────────────┤
  │ B5  EXACT FORWARD VERIFIER: fragmenter + maximum-weight peak   │
  │     assignment → log P(spectrum | candidate) per spectrum      │
  ├─────────────────────────────────────────────────────────────┤
  │ B6  POOLING: Σ q_s · LL_s  −  contradiction penalty            │
  ├─────────────────────────────────────────────────────────────┤
  │ B7  DIVERSITY-AWARE TOP-25: avoid 25 isomers of one wrong idea │
  └─────────────────────────────────────────────────────────────┘
```

### 6.1 B1 — The Substructure–Loss Lexicon (mined offline, counted, not learned)

**What it is.** A table of recurring spectral events and the substructures they imply. It is
built from three offline passes over the training data, and it is the *only* chemical
knowledge Track B is allowed to use.

**Pass 1 — label the losses (MAGMa-style substructure annotation).**
For each of the 2.5 M training spectra we know the structure, so:

```
for structure s with spectra S_s:
    fragments(s) = { all connected subgraphs obtainable by
                     ≤ 3 bond cleavages, PLUS single mobile losses (H2O, NH3)
                     applied to every fragment with an OH/NH }        # MAGMa (Ridder 2012)
    for each observed peak p in S_s:
        for each fragment f with |m(f) − p| ≤ tol:
            record  (loss = precursor − p)  →  atoms(f)  as a formula-level annotation
```

**Pass 2 — mine the motifs (co-occurrence, MS2LDA-style).**
Two unsupervised passes over the annotated corpus:

```
(a) frequent-itemset mining (Apriori) over discretised (fragment mass, loss mass) pairs
    → recurring "Mass2Motifs": sets of fragments/losses that always co-occur;
    the MS2LDA result (Rogers PNAS 2016; Torres Ortega Nat Commun 2026) is exactly
    this done with topic models; Apriori is simpler and deterministic.
(b) for each motif, the intersection of the Pass-1 annotations of its member peaks
    → a candidate substructure (SMARTS + formula), ranked by coverage.
```

**Pass 3 — calibrate.**
For every lexicon entry `m = (signature: {frag masses, loss masses}, substructure,
formula)`, record:

```
freq(m)         count / total       — used as the emission prior P(motif | spectrum)
conditional     P(loss | m)          — the loss plausibility table used by S4 and B5
class_profile(m) normalised histogram over NP classes (NPClassifier taxonomy)
```

**The lexicon at a glance (illustrative entries; the real one is mined, ~2000–8000 entries):**

| signature | substructure implied | typical loss | class profile |
|---|---|---|---|
| {frag 91.05, loss 14.02} | benzyl (Ph–CH2–) | CH₂ | phenylpropanoids, alkaloids |
| {frag 105.03, loss 31.02} | benzoyl (Ph–CO–) | OCH₃ | flavonoids, phenolics |
| {loss 162.05 then 180.06} | hexose glycoside | C₆H₁₀O₅ | flavonoid glycosides, saponins |
| {loss 44.00 strong} | carboxylic acid | CO₂ | phenylpropanoids, fatty acids |
| {frag 149.02, 121.03} | coumarin core | CO | coumarins |
| {frag 135.04} | methoxy-benzyl | CH₃· | lignans, flavonoids |
| {loss 27.99 repeated ×2} | fused polycarbonyl | CO | flavonoids, quinones |

**Cost.** Pass 1 over 2.5 M spectra with a 3-cut fragmenter is the single largest offline
job (≈ days on one core, or hours on a small cluster) — and it is *offline*: the artifacts are
uploaded once (the notebook rules allow pre-attached data and pre-trained models). At
inference the lexicon is a 30 MB parquet.

**Boundary rule (critical):** the lexicon stores *motifs*, never molecules. A lookup by
signature returns substructures and their frequencies — never a candidate structure. This is
what keeps Track B honest: it cannot accidentally become a library search.

### 6.2 B2 — The Jigsaw Solver: substructure inventory by exact set-cover

**What.** Given the mass, the formula candidates, and the loss lattice, recover the *bag of
substructures* the molecule must contain. This is the step that turns "identify a molecule"
into a problem with an exact solution.

**Why a set-cover formulation and not a fragment search.** The spectrum's information about
structure lives almost entirely in (i) *which pieces broke off* (neutral losses) and (ii)
*which pieces stayed together* (fragment ions). Both are statements about a partition of the
molecule into substructures. So we search over partitions, not over graphs:

```
VARIABLES     x_m ∈ ℤ≥0    count of lexicon motif m in the inventory
              y_l ∈ {0,1}  1 if loss l of the lattice is "explained by composition"

CONSTRAINTS   (1) MASS BALANCE (exact, elementwise):
                  Σ_m x_m · formula(m)  ==  formula(M)        [for each element]
              (2) COVERAGE of the strong losses:
                  for each loss l with intensity share > 5 %:
                      Σ_{m : m can lose l} x_m ≥ 1     (some motif that could carry l)
                  else y_l = 0
              (3) CAPACITY: Σ_m x_m ≤ max_motifs (from formula size / 12)
              (4) VALENCE:  Σ_m x_m · open_valence(m) ≥ 2·(components − 1)
                  (the pieces must be connectable into ONE molecule)

OBJECTIVE    maximise   Σ_l y_l · i_l                      (explained loss intensity)
                        − λ1 · Σ_m x_m                      (parsimony: fewer pieces)
                        − λ2 · Σ_m x_m · penalty(m)         (chemically implausible motifs)
                        + λ3 · Σ_m x_m · log freq(m)        (prior frequency)
                        − λ4 · (unexplained strong-peak intensity)
```

**How it is solved.** The variables are few: motifs whose formula is a sub-formula of
formula(M) and whose mass fits — typically 200–5000. Solve by:
1. LP relaxation + rounding → a strong incumbent in milliseconds;
2. local repair: add / drop / swap moves evaluated by the true objective, iterated to a
   local optimum;
3. bounded branch-and-bound (depth ≤ 5, best-first) when the incumbent's coverage is below
   a threshold, to escape local optima.

**Worked example (illustrative, small numbers).**

Query: `M = 302.0427` → formula `C15H10O7`. Loss lattice (strong): 18.01 (H₂O), 30.01
(CH₂O), 43.99 (CO₂), 122.04 (C₇H₆O₂), plus a 178.9986 → 151.0037 → 123.0088 CO chain.

| lexicon motif | formula | mass | note |
|---|---|---|---|
| benzoyl-type catechol aryl | C₇H₆O₂ | 122.037 | can lose C₇H₆O₂, H₂O |
| flavone/chromone core | C₉H₄O₂ | 148.016 | can lose CO, CH₂O |
| hydroxyphenyl (phenolic) | C₆H₅O | 93.034 | can lose H₂O, CO₂ |
| methoxyaryl | C₇H₇O | 107.050 | can lose CH₃·, CH₂O |

Balance solve against C15H10O7 (15 C, 10 H, 7 O):
- one C₇H₆O₂ (7,6,2) + one C₉H₄O₂ (9,4,2) = C16H10O4 — one carbon too many and 3 O short.
- one C₇H₆O₂ + one C₆H₅O + two CO-type carbonyls (2C, 2O) = C₇₊₆₊₂ = 15 C, H 6+5 = 11, O 2+1+2 = 5. H too many by 1, O short by 2. Close but not balanced — the solver keeps iterating.
- one C₇H₆O₂ + one C₉H₄O₂ − CH is inadmissible (no negative entries); the solver adds one
  "hydroxyl/oxygen insert" motif: C₉H₄O₂ + C₇H₆O₂ = C16H10O4 → over by one C; so instead
  the *carbonyl-substituted* chromone variant C₈H₄O₃ (152.011) is tried:
  C₈H₄O₃ + C₇H₆O₂ = C15H10O5, short 2 O → + one extra phenolic OH motif C₆H₄O (92.027) made
  of shared atoms is not allowed (overlap), so the solver adds two carbonyl O-insert motifs
  (each +O, no atoms): C15H10O7. **Balanced.**

Result inventory: {catechol aryl C₇H₆O₂, oxygenated chromone core C₈H₄O₃, 2 × carbonyl-oxygen
insert} — which explains every strong loss (H₂O from the catechol hydroxys, CO₂ from an acid
or ortho-quinone, CH₂O from the methoxy-adjacent position, C₇H₆O₂ as the leaving catechol,
and the CO chain from the fused carbonyls), and it is exactly the substructure decomposition
of a flavonol such as quercetin.

**What the solver really produces:** not one inventory, but a *ranked list of inventories*
(the branch-and-bound keeps the incumbent path), each with a score. B4 consumes the top few
as different *seeds*, and this is where candidate diversity comes from for free.

**Failure modes and fixes:**
- *No balanced inventory exists* (novel chemistry whose motifs were never mined) → the
  parsimony term is relaxed (allow "unknown" filler motifs with a fixed penalty), so B2 still
  returns a partial answer; the fill slots then carry popular same-formula structures.
- *Many balanced inventories* → correct behaviour: that is the under-determination showing
  itself, and the ranked list is the honest posterior. B4/B5 will order them.

### 6.3 B3 — Class prior from the spectrum (algorithmic, no deep net)

**What.** Predict the natural-product class distribution of the molecule from its inventory
and motifs — and use it to constrain generation.

**Why it matters so much.** Generation over unconstrained chemistry is hopeless. If we know
the molecule is a flavonoid glycoside, the search space collapses by orders of magnitude:
only specific cores, specific glycosylation patterns, specific substitution positions are
admissible.

**How, algorithmically.** A naive-Bayes / kernel construction over *counted* motif features:

```
P(class | motifs) ∝ P(class) · Π_m P(m | class)
    P(m | class) = smoothed counts of motif m among structures of that class
                     (offline, over COCONUT + train, using the NPClassifier taxonomy:
                      7 pathways / 70 superclasses / 672 classes)
    for the "spectrum declares motif" step, use the B1 lexicon conditional distributions
        P(m | observed loss signature) rather than hard assignments
```

This is a *log-linear model fitted by counting*, i.e. the same object CSI:FingerID/CANOPUS
use, minus the deep network (CANOPUS itself shows a kernel/linear stage is competitive; the
deep net buys a few points of class accuracy). It returns, per molecule, a distribution over
classes and — critically — two lists:

- **good list**: motifs whose posterior ≥ 0.5 (must appear in the candidate);
- **bad list**: motifs whose posterior ≤ 0.02 (must not appear).

These are exactly the "good list / bad list" constraints that 40 years of CASE systems
(StrucEluc, ACD/Structure Elucidator, Sherlock/PyLSD) feed into structure generators. We are
reconstructing the CASE pipeline with MS/MS as the constraint source instead of 2D NMR.

### 6.4 B4 — The Seed-and-Grow Assembler (generation by search, never by enumeration)

**What.** Given the formula, the substructure inventory, and the class constraints, construct
candidate structures by *beam search over graphs, one atom or one bond at a time, with the
spectral objective inside the loop.*

**Why not MAYGEN-style enumeration.** Enumerating all isomers of C15H10O7 and scoring each
with a fragmenter costs (≈10⁶–10⁷ isomers) × (a fragmenter run). That is the budget bomb. The
NMR structure-elucidation field solved this exact problem a different way: generate graphs
*bond by bond under a guided score* and keep only a beam. Huang et al. (Chem Sci 2021) reached
**67.4 % top-1 / 95.8 % top-10** on structure prediction that way, and their own ablation
shows the generator loses almost nothing (top-10 drops 0 % at beam 1000) — the error comes
from the guidance, not the search. Our guidance signal is *better than NMR's* for this task
because we have a direct mechanistic scoring function (B5) and an exact inventory constraint
(B2).

**State.**

```
state = {
    graph G          : partial molecular graph (atoms with element labels, open valences)
    motif_placements : which lexicon motifs are currently embedded in G
    remaining        : element counts not yet placed
    score            : g(state), see below
}
```

**Expansion operators** (the only ways to change a state):

```
E1  ATTACH ATOM      place one atom of `remaining` on an open valence
E2  CLOSE RING       bond two open valences inside one component
E3  EMBED MOTIF      fuse a mandatory/optional lexicon motif onto an open valence
                    (optionally through a 1–2-atom linker: O, CH2, C=O)
E4  TERMINATE        all elements placed, graph connected, valence closed → complete
```

**The guidance score** (evaluated incrementally, cheaply):

```
g(state) =  w1 · substructure_consistency(state)      # are the B2 motifs embedded? partial = graded
          + w2 · class_prior(state)                    # good-list motifs present, bad-list absent
          + w3 · plausibility_prior(state)             # counted frequencies of motif pairs/fusions
          + w4 · loss_coverage(state)                  # which lattice losses are already explained
          − w5 · strain(state)                         # 3-rings, gem-diols, acetal clusters, ...
```

**The search.**

```
beam B ← { for each of the top K inventories from B2: the state = all its motifs embedded }
repeat until beam empty:
    for state in B:
        for op in {E1, E2, E3}:            # E4 handled at acceptance
            for each instantiation:
                if violates HARD constraints → prune
                push (child, g(child)) into the child pool
    dedup children by canonical graph form (RDKit canonical ranking / nauty);
    keep top W by  g(child) + h_max(child)
B ← survivors
every complete state → run B5 verifier → candidate + LL

HARD constraints (never violated, checked at expansion):
   H1  element counts exactly match formula(M) at completion
   H2  valence of every atom ∈ {allowed valences for its element and hybridisation}
   H3  RDBE consistent with the bond orders created
   H4  every good-list motif is embedded by completion; no bad-list motif ever appears
   H5  connected, single component
   H6  no unusual chemistry: no peroxide chains, no N–N, no metal, ...
```

**Admissible bound (why the beam is safe).** `h_max(child)` = the largest possible
`loss_coverage` still achievable with the elements remaining (a knapsack over the lattice
losses vs remaining mass). If `g(child) + h_max(child)` cannot beat the beam's 1000th score,
the child is pruned. This is A\* logic applied to molecule assembly, and it is what makes
beam width 200–2000 sufficient.

**Cost.** States: W (200–2000) × depth (≤ 40 atoms) × branching (≤ ~20 distinct
instantiations after dedup and hard pruning) ≈ 10⁶–10⁷ cheap state evaluations per inventory,
each a few microseconds → seconds to a minute per molecule. × 400 molecules × K inventories
(≈3) is comfortably inside the 9 h budget on CPU.

*WHY THIS BEATS THE CONSENSUS (and plan.md):* the objective being optimised *is* the
spectral objective, evaluated incrementally as the graph is built, so the beam concentrates on
structurally plausible explanations of the observed losses instead of enumerating the whole
isomer space and hoping a fragmenter separates the good from the astronomically many bad. It
also produces *diverse* candidates for free: different inventories from B2 seed different
beams, which is exactly the posterior spread Law 1 wants.

### 6.5 B5 — The Exact Forward Verifier (the forward half of the channel)

**What.** Compute `log P(spectrum | candidate)` exactly enough to rank and to pool across
spectra — with no learned model in the scoring path.

**Step 1 — enumerate the candidate's fragment multiset.**

```
F(cand) = { connected subgraphs obtainable by }
           ∪  cutting 1 bond          (bridges → clean two-piece splits)
           ∪  cutting 2 bonds          (ring openings)
           ∪  cutting 3 bonds          (fused-ring systems, capped)
           each × { as cation, as radical, as anion }   (mass ± H, per ionisation mode)
           ∪  rearrangement losses applied to every fragment with the right heteroatom:
                H2O, NH3, CO, CO2, CH2O, CH3OH, HCOOH, ...
                (mobile-hydrogen model: any fragment bearing OH/NH/COOH can shed these)
```

The enumeration is done with bridge-finding (Hopcroft–Tarjan / chain decomposition) and
2-cut enumeration — the standard fast algorithms, the same ones MAGMa and MetFrag use — so
depth-3 fragmentation of a 40-heavy-atom molecule is ~10⁴–10⁵ fragments, computed in
milliseconds with hashing on formula.

**Step 2 — maximum-weight assignment of observed peaks to fragments.**

```
This is a bipartite matching between observed peaks (colors) and enumerated fragments,
each fragment usable at most once, maximising explained intensity:

    maximise  Σ_p i_p · [p matched]  −  γ · Σ_f [f used but matching no peak]·implausible(f)

solve greedily for the common case (each peak takes its best in-tolerance fragment),
    Hungarian when a peak competes for a fragment with several peaks (rare, ≤ 20 peaks)
```

**Step 3 — turn coverage into a log-likelihood with counted emission probabilities.**

```
LL(cand, spectrum) = Σ_p log Pe(p | cand)

Pe(p | cand) = 1 − exp( − matched_share_p )        if p matched by a fragment
             = 1 − P_peak_absent(bucket(p), CE_bucket)   if p unmatched
```

where `P_peak_absent` is a *counted* quantity: for each mass bucket and collision-energy
bucket, how often a genuine fragment of a known structure goes unobserved. This is the
documented "20–30 % of peaks are unexplainable" reality, converted into a likelihood instead
of a bug.

**What this buys over MetFrag-lite scoring.** `frag_score` (the consensus channel) is a
coverage fraction — uncalibrated, unitless, and not poolable across spectra or comparable
against library similarity. `LL` is a proper per-spectrum log-likelihood, which is what makes
B6, calibration, and the cross-track fusion mathematically coherent.

### 6.6 B6 — Evidence pooling across the molecule's spectra

```
LL_total(cand) = Σ_s  q_s · LL(cand, s)        # q_s = S3 quality weight
                 − R · contradiction_penalty(cand)

contradiction_penalty(cand) = Σ over spectra s with q_s ≥ 0.7:
        1[ explained_share(cand, s) < 0.25 ] · log(1 / 0.25)     # ≈ 1.4 nats per veto
```

Interpretation: independent observations pool by summing log-likelihoods (statistics); a
high-quality spectrum that the candidate barely explains is evidence *against* it, priced in
nats so it can be compared with everything else (the docs' min-rule, in its principled form:
one good spectrum can outvote six weak ones *if* it is informative, not merely contrary).

### 6.7 B7 — Diversity-aware top-25

25 near-identical isomers of one wrong structure scores ~0.04 if wrong. The list must cover
the *plausible hypothesis space*:

```
1. rank all complete candidates by posterior (B6 output × class prior × plausibility).
2. greedily select, at each step, the best candidate that is not too close to an
   already-selected one — "too close" = MCES / fingerprint Tanimoto above a threshold
   (~0.7) or the same formula and the same skeleton-with-one-substituent-moved.
3. if fewer than 25 remain, widen: allow closer neighbours, then same-formula popularity
   fill (the documented +0.01 prior), then ±CH2/±O neighbours of the formula.
```

*WHY THIS BEATS A MODEL:* the objective here is the metric itself. MRR@25 rewards covering
the posterior mass; posting 25 copies of the posterior mode is a metric bug, not a modelling
failure, and it is common.

---

## 7. FUSION & SUBMISSION ENGINE

```
  Track A output: [(cand_i, pA_i, provenance_i, lib_gate_i)]     pA = calibrated posterior
  Track B output: [(cand_j, pB_j, provenance_j)]                 pB = class·plausibility·LL

  1. UNION the two lists (dedup by InChIKey14; a candidate found by both tracks keeps both
     provenances — agreement is itself evidence).
  2. FINAL SCORE per candidate:

        logit = wA·logit(pA) + wB·logit(pB) + corrections

     corrections:
        + μ_pop · (log(1+pubchem_substances) + log(1+pubmed))      # popularity prior
        + λ_agree · 1[appears in both tracks] · min(pA, pB)         # agreement bonus
        − if lib_gate(cand):  zero every term that came from a simulator
  3. SORT by final score descending  (this is the MRR-optimal order — Law 1)
  4. DEDUP InChIKey14 (stereochemistry/tautomer are free — the scorer ignores them, so
     two entries that reduce to the same InChIKey14 would waste a slot).
  5. FILL to exactly 25:  B7 diversity pass, then popularity fill.
  6. WRITE submission.csv:
        molecule_id,smiles           best first, semicolon-joined, ≤ 25, no nulls,
        every molecule_id exactly once.
```

**The one rule that must never be broken:** a candidate whose evidence came from a *real
measured spectrum* is never outranked by a candidate whose evidence is *simulated or
generated* unless the simulated/generated posterior exceeds it by a margin. This is the docs'
library gate generalised across tracks, and it is enforced at step 2 above.

---

## 8. Every documented problem, and the component that kills it

| Problem (docs) | Where it is solved | Mechanism |
|---|---|---|
| Isomers indistinguishable in one spectrum (§12) | Law 1, B6, B7, §7 | Emit a calibrated posterior; pool the molecule's spectra by summed log-likelihood; fill 25 with diverse posterior mass instead of repeats |
| Under-determined spectra | §6.4 B4 beam, §6.5 B5 | The spectral objective inside the search; exact forward verification prices every hypothesis |
| Three hidden difficulty classes (§11) | Track A + Track B + fusion | Both tracks always run; fusion sorts by posterior; library gate protects measurements |
| Class 3 is unsolved (SOTA < 5 % exact) | Track B | Substructure inventory (B2) + class-constrained guided generation (B3/B4) + exact verifier (B5) — structurally different from diffusion/autoregressive generation |
| Validation trap: random split ranks backwards (§7) | §10 | Only `enveda-np-examples` under three masks; everything else diagnostics |
| Chemistry mismatch: 46 % of spectra are drug-like synthetic (§13.3) | A2 AFIX/BIO variants, A5, B1-B4 | Instrument- and chemistry-restricted search variants; Track B never touches spectra libraries at all |
| Instrument mismatch (test 100 % timsTOF) | A2 variants, A4 weighting, S3 | Nothing cross-instrument is trusted unconditionally; modified metrics instrument-weighted |
| Collision-energy unit chaos (§7.2 Part 1) | S3, A2, B5 | Use `collision_energy_ev` only; bucket it; energy-aware emission model |
| Merged multi-energy spectra | S3 bucketing, B5 | Treated as extra evidence with an energy-conditioned emitter, not as a problem |
| Noise / weak spectra | S3 quality weight q, B6 | Per-spectrum reliability weights and contradiction penalties |
| Mislabeled training rows | S1/S2, A1 | `precursor_error_ppm` filter; neutral-mass cross-check; formula fit score |
| Duplicate structures across libraries | A1 dedup, §7.4 | Dedup by InChIKey14 everywhere |
| Popularity signal unused | §7 step 2, A7 | Prior over the label distribution (+0.01–0.02, docs-phases) |
| Simulators outranking measurements | §5.4 gate, §7 | Hard-coded, not learned |
| 25-slot waste (rank duplicates, near-identical isomers) | §6.7 B7, §7.5 | Diversity-aware selection against the metric |
| 9 h, no internet, one notebook | §9 | Everything expensive is offline; runtime is retrieval + search, ~40–60 min |
| Leakage in naive splits (research §9) | §10 masks | Mask-2/mask-3 are the leakage-controlled evaluations |

---

## 9. Research lineage: analogous problems, how they were solved, what we take

> The user's requirement: deep research *outside* the competition — find similar problems,
> how they were solved, and how each solution transfers here. Each entry states the problem,
> the prior art's solution, what we adopt, and **what we change**.

### 9.1 De novo peptide sequencing (proteomics) — *the closest true analogue*

**Problem.** Reconstruct a molecule you have never seen from its MS/MS spectrum, with no
database — solved for 25 years.

**How it was solved.**
- *Spectrum graph + DP* (Dančík/SHERENGA 1999; Ma PEAKS 2003; Frank & Pevzner 2005):
  nodes are peaks, edges are amino-acid mass gaps; the answer is the best path. DP over the
  mass axis handles missing ions (PEAKS scores by mass, not by observed peak).
- *Probabilistic networks* (PepNovo): a Bayesian network whose structure encodes the physical
  fragmentation rules, with learned ion-type probabilities.
- *Multi-spectra DP* (Chen et al.): sequences shared across MS²/MS³ spectra, DP over
  fragment-type assignments.
- *Sequence tags*: seed short reliable spans, then extend.
- Later: deep models (DeepNovo, Casanovo) — plus the NovoBoard finding that classical
  DP/graph methods are still robust baselines that deep models fail to beat on out-of-domain
  instruments.

**What we take.**
1. The **spectrum-graph/DP abstraction** motivates S5 (loss lattice) and the beam search in
   B4: both are "walk the mass axis with a scored alphabet".
2. **Multi-spectra DP** is the justification for B6: evidence from different adducts and
   energies of the *same* molecule must be pooled as one inference problem, and the DP
   formulation sums log-likelihoods.
3. **Robustness of classical methods on domain shift** is the empirical reason Track A keeps
   instrument-restricted arithmetic search (A2/A3/A4) rather than a single global model.
4. NovoBoard also shows the right way to *evaluate* unknowns: target/decoy FDR analysis —
   adopted in §10.

**What we change.** Peptides are linear sequences over a fixed 20-letter alphabet;
metabolites are **graphs over ~10 elements with cycles**. So the "path" becomes a tree
(fragmentation tree, S4) and the "sequence tag" becomes a *substructure tag* (lexicon entry,
B1) — and the final assembly is graph search (B4), not sequence DP.

### 9.2 Fragmentation trees, formula annotation, and joint inference (SIRIUS family)

**Problem.** From one spectrum: (i) determine the molecular formula, (ii) annotate every peak
with a formula, (iii) rank formulas, (iv) do it for a whole dataset at once.

**How it was solved.**
- *Fragmentation trees* (Böcker & Rasche 2008; Rasche 2011; Dührkop & Böcker 2016): model the
  cascade as a tree; the computational problem is the **Maximum Colorful Subtree**, NP-hard
  but solved exactly by ILP, with fast heuristics ("Speedy Colorful Subtrees", Dührkop 2018)
  that rank the true tree near-optimally ~100× faster. Scoring reformulated as MAP estimation
  with a prior over neutral losses.
- *ZODIAC* (Ludwig 2020, Nat Mach Intell): formulas in one sample are jointly annotated by
  Gibbs sampling over the similarity of their fragmentation trees — 16× fewer errors.
- *MSBuddy* (2024): tree-search formula annotation with isotope awareness.
- *CANOPUS* (Dührkop 2021, Nat Biotechnol): predict 2 497 compound classes from spectra
  **without any database**, by predicting a fingerprint and classifying it.
- *CSI:FingerID* (Dührkop 2015; Ludwig 2018): kernel/SVM over fragmentation-tree features →
  molecular fingerprint → structure database search; a **Bayesian network** version exists.

**What we take.** S4 wholesale (MWCS/ILP, MAP scoring, loss prior counted from the training
corpus); a CANOPUS-style class channel (B3) but linear/kernel rather than deep; and ZODIAC's
insight becomes S7 (joint inference over the ~400 test molecules).

**What we change.** SIRIUS stops at formula + fingerprint + database search; we push the same
object (the fragmentation tree) further — into *structure generation* (B4) via the lexicon and
the jigsaw inventory (B2), because the competition scores structure identity, not formula.

### 9.3 Computer-Assisted Structure Elucidation (CASE) — 40 years of "generate then filter"

**Problem.** Determine an unknown structure from spectroscopy (originally NMR + IR + MS).

**How it was solved.** The classical pipeline: (1) molecular formula from HRMS; (2) extract
*structural constraints* from spectra into a **good list / bad list** of fragments (StrucEluc,
ACD/Structure Elucidator); (3) **fuzzy or strict structure generation** from fragments + free
atoms under the formula constraint (MOLGEN, later LSD/PyLSD in Sherlock); (4) filter and rank
candidates by predicted-vs-measured spectra (HOSE-code shift prediction); (5) if a candidate
contradicts the data, re-run with a corrected constraint set.

**What we take.** The entire architecture of Track B:
- B2's inventory = the **good list** (mandatory substructures);
- B3's class posterior = good/bad lists at scale;
- B4's assembler = **structure generation from fragments under a formula constraint**;
- B5 = the "rank by predicted vs measured spectra" step.

**What we change.** CASE's constraints come from 2D NMR correlations, which are *near-exact*
(an HMBC cross-peak almost certainly means 2–3 bonds). MS/MS losses are noisy, incomplete,
and ~25 % rearranged. So our constraints are *soft, scored, and probabilistic* (counted
frequencies), and the generator is a **beam search with a spectral objective** rather than a
strict CSP — and when the constraints are ambiguous, we keep the ambiguity in the posterior
instead of forcing a single answer.

### 9.4 Structure generation: orderly generation and canonical generation paths

**Problem.** Enumerate the constitutional isomers of a formula, exhaustively and without
duplicates.

**How it was solved.** MOLGEN (Beyer/Kerber): orderly generation with canonical labelling.
MAYGEN (Yirik 2021): orderly generation in Java, near-MOLGEN speed. Surge (McKay, Yirik,
Steinbeck 2022): **canonical generation path** using nauty automorphism groups —
orders of magnitude faster, and accepting *fragment-based* starting points and badlists.

**What we take.** Two concrete things:
1. Surge's ability to start from a **fuzzy set of fragments** rather than a bare formula —
   this is exactly B4's seeding from B2's inventory.
2. Badlist semantics → B3's *bad list*.
We may literally instantiate B4 with Surge as the generator core and our own scorer as the
search guidance, if re-implementing the operators is not worth the time.

**What we change.** We do not enumerate to completion. Enumeration is only correct as a
*complete* answer; we need a *ranked posterior over partial explanations*, so we search.

### 9.5 The NMR framework that proves guided graph generation works

**Problem.** Predict a molecular graph from spectra of a *different* modality.

**Solution (Huang et al., Chem Sci 2021).** A bottom-up rule-based pipeline: predict a
**substructure probability profile** from IR/¹H/¹³C spectra, then generate molecular graphs
**one bond at a time with beam search guided by that profile**, output ranked constitutional
isomers.

**Result.** Beam 1000: **67.4 % top-1, 95.8 % top-10**; their ablation shows the generator
loses essentially nothing at large beam — the errors come from the guidance.

**What we take.** The architecture of B4 verbatim, with our substitution:
- NMR's "substructure probability profile" ← our B1/B2/B3 inventory-and-class posterior;
- NMR's scoring ← our B5 spectral likelihood.
We inherit a *demonstrated* top-1 rate for exactly this algorithm family, which is why this
plan claims class-2/3 gains that `plan.md`'s "enumerate then MetFrag" cannot deliver.

### 9.6 Unsupervised substructure discovery from spectra (MS2LDA / Mass2Motifs, MAGMa)

**Problem.** Learn what recurring fragment/loss patterns mean, without structure labels.

**How it was solved.**
- *MS2LDA* (Rogers PNAS 2016; Wandy 2018; Torres Ortega Nat Commun 2026): discretise
  fragments and neutral losses as word-like features; latent Dirichlet allocation finds
  **Mass2Motifs** — recurring fragment/loss sets that correspond to substructures; 2.0 adds
  automated annotation guidance (MAG) reaching 0.75–0.93 substructure overlap.
- *MAGMa* (Ridder 2012): annotate MSⁿ data by in-silico fragmenting *candidate molecules*
  (up to 3 bond dissociations) and matching substructures to peaks; CASMI-2013 participant.

**What we take.** B1's mining pipeline in its entirety: fragment/loss discretisation,
co-occurrence mining, and MAGMa-style annotation of known structures to attach *meanings* to
motifs. The lexicon is the deliverable.

**What we change.** We make motifs carry *formulas and counts* (so B2's mass balance is
exact), and we mine on the 2.5 M-spectrum competition corpus + COCONUT rather than a
community corpus. The lexicon therefore speaks timsTOF and natural products — the two domains
the test set actually lives in.

### 9.7 Library search in loss space (NIST Hybrid Similarity Search) and tree alignment

**Problem.** Extend a spectral library's coverage to compounds that differ by a known moiety,
and compare compounds across instruments.

**How it was solved.**
- *Hybrid Similarity Search* (Moorthy et al., Anal Chem 2019): score = match on **both**
  fragment ions and **neutral losses**, with a `ΔMass` parameter shifting the loss axis; it
  extends a library to near-neighbour compounds differing by one inert structural component
  (their fentanyl-family case study).
- *Modified cosine* (consensus channel): fragments matched with a mass offset.
- *Fragmentation tree alignment* (Rasche, Anal Chem 2012): identify unknowns by aligning
  their fragmentation trees to tree libraries.

**What we take.** All three, in the order that matters: loss-space HSS as the primary analog
channel (A2), modified cosine as a secondary feature, tree alignment as a third (A3).

**What we change.** We condition HSS on instrument- and chemistry-restricted index variants
(AFIX/BIO), and we align *formula-annotated* trees (from S4) rather than raw peak trees, which
makes alignment robust to noise and to mass shifts.

### 9.8 Adjacent fields that shaped the *frame*

| Field | Solution | What we take |
|---|---|---|
| Genome assembly (de Bruijn/contigs) | k-mer graph, contigs, scaffolding, repeat resolution | Substructure islands as "contigs"; B2 solves the scaffolding problem (how pieces connect); MCES/dedup is the repeat-resolution analogue |
| Error-correcting codes / noisy-channel decoding | MAP decoding, Viterbi over a trellis, soft information | The whole Track B philosophy: the fragmentation tree is the trellis; probabilities, not boolean matches |
| Information theory | channel capacity, erasure channels | Honest accounting: a 20-peak spectrum carries far fewer bits than a skeleton choice, hence MRR@25 and posterior output rather than a hard answer |
| De novo peptide *tags* + homology search | seed short spans, extend by similarity | Motif tags as seeds (B1→B2), analog extension in loss space (A2) — the two tracks are the tag/homology duality |
| Factor graphs / Gibbs sampling | joint inference over coupled variables | S7's mixture-level consistency pass |

---

## 10. Compute budget, artifacts, module layout

### 10.1 What must be precomputed offline (uploaded as attached data; no internet at runtime)

| Artifact | Size | How long to build | Used by |
|---|---|---|---|
| Structure pool (SMILES, InChIKey14, exact mass, formula, Morgan, motif set) | ~2–4 GB | hours | A1, A5, §7 fill |
| Spectral index (peaks per training spectrum, adduct, library, quality) | ~3 GB | hours | A2, A4 |
| Loss-space inverted index | ~1 GB | hours | A2 |
| Fragmentation trees + compact signatures for all training spectra | ~2 GB | days | A3 |
| **Substructure–Loss Lexicon** (motifs, formulas, frequencies, class profiles) | ~30 MB | days | B1–B5, A5 |
| Loss-frequency table `P_loss` | ~2 MB | minutes | S4, B5 |
| NP class priors over COCONUT/train structures | ~5 MB | minutes | B3 |
| Popularity table (PubChem substance counts, PubMed citations) | ~50 MB | minutes | §7 prior |
| GBDT ranker + calibration map | ~5 MB | hours | A7 |
| (optional) kernel fingerprint model | ~20 MB | hours | A6 |

### 10.2 Runtime budget (single notebook, ≤ 9 h, CPU or GPU)

```
S1+S2+S3                     2 min      trivial arithmetic
S4 fragmentation trees      25 min      ILP per spectrum (~1 s)
S5 loss lattice             5 min
S6/S7 joint pass            10 min
A1 pool retrieval            5 min      indexed lookups
A2 loss-space search        15 min
A3 tree alignment           10 min
A4 entropy similarity       10 min
A5 motif scoring            10 min
A6 optional kernel channel   0–20 min
A7 ranking + calibration     5 min
TRACK B (the expensive one):
   B1 lexicon lookup         5 min
   B2 jigsaw solver         30 min      ILP/LP per molecule
   B3 class prior           10 min
   B4 beam search           60–90 min   the dominant cost
   B5 forward verifier      30 min
   B6/B7 pooling + fill      5 min
FUSION + write CSV           2 min
────────────────────────────────────
TOTAL                        ~4 h worst case; ~2 h typical.
```

Headroom remains for reruns, ablations, and the optional model channel. (The published
solutions run in ~20 min because they skip generation entirely; our Track B is the reason we
use more, and it is still well inside the limit.)

### 10.3 Module layout

```
sic/
  adducts.py        S1 inversion + cross-spectrum vote + adduct sanity
  formula.py        S2 integer solve, isotope scoring (S6)
  spectrum.py       S3 conditioning + quality weight + energy bucketing
  fragtree.py       S4 MWCS: graph build, ILP (scipy/ortools), beam heuristic
  losslattice.py    S5 pairwise differences + chained paths + loss prior
  joint.py          S7 mixture consistency pass
trackA/
  pool.py           structure pool build + mass/formula indices
  hssearch.py       A2 loss-space HSS (offline index + online search)
  treealign.py      A3 tree alignment DP
  entropy.py        A4 spectral entropy, instrument-weighted
  motifs.py         A5 motif matching
  kernelchan.py     A6 optional
  ranker.py         A7 GBDT + evidence algebra + calibration
trackB/
  lexicon.py        B1 loading + signature matching
  jigsaw.py         B2 ILP/LP + local repair
  classprior.py     B3 naive-Bayes over motifs
  assembler.py      B4 seed-and-grow beam search
  verifier.py       B5 fragmenter + peak assignment + LL
  pooling.py        B6/B7 pooling + diversity selection
fusion/
  merge.py          union, evidence algebra, dedup, fill, CSV writer
eval/
  masks.py          three-mask evaluation, per-class diagnostics, decoy FDR
```

---

## 11. Validation protocol (do not skip: a random split ranks methods backwards)

```
VAL SET:  train.parquet[ingest_lib == 'enveda-np-examples']      # 250 structures, 1151 spectra
           (the only library that matches the test instrument AND the test chemistry)

THREE MASKS:
  M1  everything visible           → class-1 ceiling (pool has structure + its spectra)
  M2  own spectra hidden           → class-2 behaviour (structure in pool, no reference)
  M3  structure removed from pool  → class-3 behaviour (genuinely novel)
```

**Protocol.**
1. **Stage-gated ablation.** Never compare whole systems; compare one channel at a time,
   always reporting M1/M2/M3 separately. A channel that improves M2 but wrecks M1 is a
   regression, not a win (M1 is the largest slice of the test set).
2. **Track B has its own honest test.** M3 *is* the Track B benchmark — and it is the only
   one where its score is meaningful. Report top-1/top-5/top-25 exact-InChIKey14 rate and
   mean Tanimoto of the best candidate; exact recovery will be low, and that number matters
   far less than its contribution in fusion.
3. **Decoy/FDR control** (the NovoBoard lesson): generate decoy molecules (same formula,
   random-but-valid structure) and confirm Track B does not score them above real
   structures. If it does, the posterior is miscalibrated and fusion will be harmed.
4. **Calibration check:** reliability diagram of `pA` against correctness on M1/M2. An
   uncalibrated fusion is a coin flip between the tracks.
5. **Honest caveats:** 250 structures ⇒ ±0.02 is noise. Use the public leaderboard (a real
   400-molecule held-out set) for anything finer, and at most a few submissions for it.
6. **Leakage discipline:** never let a validation structure's *other-library spectra* leak
   into the pool for M2/M3 evaluation, and never tune on the test-like subset's own spectra.

**Metrics logged per stage:** `pool_recall@pool`, `MRR@1/5/25`, per-mask, per-class estimate,
per-adduct, per-energy-bucket, mean `q`, and channel-by-channel feature importance.

---

## 12. Build order and expected score progression

| Stage | What | Independent value | Est. MRR@25 |
|---|---|---|---|
| 0 | Validation harness: three masks + metrics + decoy control | everything else depends on it | — |
| 1 | S1–S3 + A1 mass filter + dedup + popularity fill to 25 | a legitimate first submission | 0.05–0.10 |
| 2 | S4 fragmentation trees; A4 entropy similarity (instrument-weighted) | class-1 engine | 0.15–0.22 |
| 3 | **A2 loss-space HSS** + **A3 tree alignment** | the biggest algorithmic jump | 0.28–0.36 |
| 4 | B1 lexicon + B2 jigsaw + B3 class prior | Track B becomes real | (diagnostic) |
| 5 | B4 assembler + B5 verifier + B6 pooling | class-2/3 recovery | 0.33–0.39 |
| 6 | A7 GBDT + evidence algebra + calibration | principled fusion | 0.38–0.43 |
| 7 | A5 motifs + A6 optional kernel channel + S7 joint pass | + fine-grained evidence | 0.42–0.46 |
| 8 | Popularity prior μ tuned, B7 diversity fill, gate margins | the last few points | 0.44–0.48 |

**Where the headroom over the 0.44–0.47 cluster comes from** (honest statement of the claim):
the consensus gets its score from peak-space analog matching and a GBDT; the additional
mechanisms here are (i) loss-space and tree-space matching, which are strictly more
informative channels on the same data, (ii) a Track B that generates *guided by the spectrum*
rather than enumerating, and (iii) joint inference across the mixture. If (i) works as the
NIST/analog literature suggests it should, it is worth a few points alone; (ii)–(iii) are
small but are the only sources of class-3 progress available.

---

## 13. Risks, failure modes, kill-switches

| Risk | Likelihood | Impact | Kill-switch / mitigation |
|---|---|---|---|
| Beam search (B4) wanders into chemically absurd structures | med | high | Hard constraints H1–H6 + plausibility prior from counted frequencies; decoy-FDR control (stage 11.3) must pass before Track B may outrank Track A |
| Lexicon mining misses the motifs of genuinely novel chemistry | med | med | Allow "unknown filler" motifs with a penalty in B2; the posterior then honestly reports low confidence and fusion defers to Track A |
| S4 ILP too slow at scale | low | med | Speedy-colorful-subtree beam heuristic; 60-peak conditioning cap |
| GBDT overfits 250 validation structures | med | med | few features (~40), heavy regularisation, stage-gated ablation, public LB for final checks only |
| Calibration drift between M1 and the hidden test | med | med | isotonic with cross-validated folds; never let an uncalibrated score fuse |
| Track B cannibalises class-1 accuracy | low | high | library gate: `lib_max ≥ 0.7` ⇒ simulated evidence zeroed |
| Budget overrun | low | high | 4 h planned vs 9 h allowed; Track B beam width is the throttle |
| Duplicate submissions wasted on noise | high | low | ±0.02 on 250 structures is noise; spend submissions on clean ablations only |
| Popularity prior mis-calibrated | med | low | μ tuned on M2 only; small by construction |

---

## 14. Two end-to-end worked examples

### 14.1 Track A: a class-1 molecule (quercetin-like), all the way to rank 1

```
INPUT
  3 spectra, all timsTOF. Precursors 301.0354 [M−H]−, 301.0354 [M−H]−, 302.0427 [M+H]+
  (two acquisitions, two adducts). Strong fragments: 271.02 (0.62), 229.01 (0.44),
  178.9986 (0.38), 273.04 (0.34), 257.05 (0.31), 121.03 (0.11).

S1  301.0354 + 1.00728 = 302.0427 ;  302.0427 − 1.00728 = 302.0427 → agree, q high.
S2  mass-defect solve → single formula C15H10O7 (ppm 1.2, RDBE 11).
S4  fragmentation tree: root C15H9O7⁻ → CO (−27.995) → CH2O (−30.011) → CO2 (−43.990)
    → C7H6O2 (−122.037) → CO chain 178.9986 → 151.0037 → 123.0088.
S7  neighbours in the mixture share the 162/180 hexose-loss family → strong consistency.

A1  pool ∩ 8.5 ppm of 302.0427 → 312 candidates.
A4  entropy similarity: candidate "quercetin" has a timsTOF reference in enveda-np-examples;
    lib_max = 0.86 → library gate fires.
A2  loss-space HSS: best shift δ = 0.00 against two other flavonol spectra, ΔMass 0,
    sim 0.91 → family = flavonols; retrieved neighbours Tanimoto 0.6–0.75 to quercetin.
A3  tree alignment 0.94, matched loss formulas 9/11, depth agreement 1.0.
A5  declared motifs {catechol aryl, chromone core, glycosyl-absent} all present in quercetin's
    motif set; motif_score 0.88.
A7  ranker posterior 0.93; library gate keeps it at rank 1 over fisetin (posterior 0.31).
FUSION  quercetin rank 1; ranks 2–5 = fisetin, morin, kaempferide, gossypetin (diverse
        isomers, popularity-ordered); ranks 6–25 = other C15H10O7 flavonols, popularity fill.
```

Score contribution if correct: 1.0. If the truth were fisetin (an honest possibility given
this data): 0.5. Both outcomes are what the posterior says they should be.

### 14.2 Track B: a class-3 molecule with no database entry and no reference spectrum

```
INPUT   two spectra, both [M+H]+, precursor 449.1084, merged energies [20,40,80].
         q = 0.8 / 0.6. Strong fragments: 287.0551, 269.0445, 241.0498, 163.0395,
         145.0289, 125.0237, 107.0497, 89.0391.

S1  449.1084 − 1.00728 = 448.1012
S2  formula candidates: C21H20O10 (ppm 0.8), C20H20N2O9 (ppm 26 → rejected),
    C22H20O9 + ... → top: C21H20O10.
S4  tree: root C21H19O10⁺ → 162.0528 (hexose loss) → 145.0289 → 125.0237, and a parallel
    branch 287.0551 → 241.0498 (−46.005, CO2) → 163.0395 (−78, benzene).
    Strong reading: **a hexose glycoside of an oxygenated aromatic acid-bearing aglycone**.
S5  lattice: 162.05 ×2 (hexose), 46.01, 78.05, 44.03, 28.03; one unexplainable peak (125.02)
    → priced into every candidate equally (rearrangement class).

B1  motif loadings: {hexose glycoside 0.9, benzoyl/benzoic-acid-type aryl 0.7,
    phenolic-OH cluster 0.6, methyl-ester/ether 0.3}.
B2  JIGSAW SOLVER on C21H20O10:
      mass balance: C12H22O11 (hexose residue in-molecule) + C9H8O4 (hydroxybenzoate-type)
      = C21H30O15 → over by 10 H, 5 O. Wrong branch.
      C12H22O11 + C8H8O3 − ... the LP keeps iterating: correct balanced branch lands on
      C12H22O11 + C9H10O4 (an aryl-propionic/hydroxyphenyl type) = C21H32O15 → still over.
      Final incumbents (top-3 inventories):
        I1: {hexose C12H22O11, catechol-aryl C7H6O2, carboxyl-bearing acyl C2H2O} − overlap
            → expressible as a hexose + a dihydroxybenzoyl-conjugated aryl + one extra OH
        I2: {hexose C12H22O11, hydroxyphenyl C6H5O, acyl C3H4O2} (anisoyl-type ester)
        I3: {deoxyhexose C12H22O9, ... } (rhamnoside variant, lower prior)
      I2 scores highest (motif coverage 0.81, frequency prior +0.35 nats).
B3  class posterior: flavonoid/aryl-glycoside family 0.62, phenolic glycoside 0.21,
    other 0.17 → good list = {hexose, hydroxyphenyl, acyl}; bad list = {indole, pyridine,
    halogen, sulfate}.
B4  seed = I2's three motifs; beam 400; expansions attach the acyl between the phenyl O and
    the anomeric carbon (ester glycoside), then close the sugar ring, then place 2 OH on the
    aryl ring (class prior). Admissible bound prunes states that cannot cover 162.05 any more.
    Beam terminates with 11 complete structures; best three:
        g1 = O[C@@H]1... (an acyl-glycoside, 3,4-dihydroxybenzoyl)  ← flat SMILES
        g2 = same skeleton, glycosylation at the other aryl OH (positional isomer)
        g3 = rhamnoside variant
B5  verifier: g1 explains 0.87 of intensity (162 hexose loss, 121/145/163 aryl series,
    46 CO2 from the acyl, 44 from the ring); g2 explains 0.83 (its 163-series is weaker);
    g3 explains 0.61. LL: g1 = −3.1, g2 = −3.4, g3 = −5.0 nats.
B6  pooling: spectrum 2 (q=0.6) independently prefers g1; totals −5.8 vs −6.3 vs −9.1.
B7  top-25: g1, g2, g3 + their positional isomers (diverse), then same-formula popularity fill.
FUSION  Track B's g1 has posterior 0.18; Track A's best pool candidate (an unrelated
        library glycoside, lib_max 0.41) has posterior 0.11 → **Track B earns rank 1**, which
        is the intended behaviour: measured-and-simulated evidence in A is weak, and the
        reconstructed structure's forward likelihood is strong.
```

That is the desired end state: **the data-only track can win, but only when the evidence
actually says so.**

---

## 15. The whole plan in one paragraph

Turn every spectrum into an exact neutral mass and a mass-defect-solved formula set, a
cleaned and quality-weighted peak list, a Maximum-Colorful-Subtree fragmentation tree, a
loss lattice with counted loss priors, and a mixture-consistency pass over the 400
co-occurring molecules. Then search the known world three different ways — loss-space Hybrid
Similarity with a global ΔMass, fragmentation-tree alignment, and entropy-weighted
spectral matching against instrument- and chemistry-restricted reference indexes — plus
dictionary matching over a mined substructure lexicon, all fused by a GBDT with hard-coded
evidence laws and calibrated to a probability. Independently, and without touching any
library, reconstruct the unknown world: annotate losses with substructure motifs, solve a
mass-balanced set-cover for the substructure inventory, derive an NP-class good/bad list,
assemble structures by seed-and-grow beam search with the spectral objective inside the loop,
and verify each with an exact bond-dissociation forward fragmenter that yields a real
log-likelihood, pooled across the molecule's spectra with quality weights and contradiction
penalties. Union both tracks, sort by calibrated posterior — the MRR-optimal order —
deduplicate on InChIKey14, fill to 25 with diverse and popular candidates, write
`submission.csv`. Validate only on `enveda-np-examples` under three masks, or not at all.
