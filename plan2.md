# PLAN 2 — Forensic Architecture for CASMI 2026
## Algorithmic-First Molecule Identification from Tandem Mass Spectra

> **Scope of this document.** This is a complete, buildable architecture. It covers every
> stage from raw parquet bytes to `submission.csv`, at both the level you can picture
> (boxes and arrows) and the level you can implement (data structures, formulas, complexity,
> pseudocode, exact hyperparameters).
>
> **Two tracks are kept strictly separate, by design and by interface:**
>
> | Track | Name | Question it answers | Independence |
> |---|---|---|---|
> | **A** | **LIBRARY SEARCH** | "Is this molecule already known to science?" | Needs only a structure database. Never looks at generated structures. |
> | **B** | **UNKNOWN FROM DATA ALONE** | "If it is not in any database, what *is* it?" | Needs **no structure database at all** for its core. Produces structures that provably cannot be looked up. |
>
> They meet at exactly one place: **a merge gate** in Part 6, where Track B's output is
> scored by Track A's evidence channels. Neither track imports a function from the other's
> scoring path. This is the separation you asked for.
>
> **Everything load-bearing is algorithmic.** Machine learning is used for exactly two
> things: (a) one feature among many in a gradient-boosted ranker over ≤ 200 rows, and
> (b) a fingerprint-probability feature. It never decides the formula, never selects
> candidates, never ranks. If you deleted every neural network, this system still runs and
> still scores.

---

## Table of contents

| Part | Title | What you get |
|---|---|---|
| 0 | [The thesis, and how this differs from `plan.md`](#part-0--the-thesis-and-how-this-differs-from-planmd) | Four structural weaknesses in the current plan |
| 1 | [The problem as a system](#part-1--the-problem-as-a-system) | The whole thing in one picture |
| 2 | [Ten design laws](#part-2--ten-design-laws) | The rules every later decision obeys |
| 3 | [Layer 0 — Ingestion, curation, and the 1 GB memory contract](#part-3--layer-0--ingestion-curation-and-the-1-gb-memory-contract) | Byte-level front end |
| 4 | [TRACK A — Library search](#part-4--track-a--library-search) | Pool builder + 10 evidence channels + ranker |
| 5 | [TRACK B — Unknown molecule from data alone](#part-5--track-b--unknown-molecule-from-data-alone) | The algorithmic de novo engine |
| 6 | [Fusion, calibration and the decision-theoretic submission](#part-6--fusion-calibration-and-the-decision-theoretic-submission) | Why you submit 25, and which 25 |
| 7 | [Validation — a 100× better harness](#part-7--validation--a-100-better-harness) | How you stop lying to yourself |
| 8 | [Build order and sprint plan](#part-8--build-order-and-sprint-plan) | 22 weeks, checkpointed |
| 9 | [Compute, memory and runtime budget](#part-9--compute-memory-and-runtime-budget) | What it costs |
| 10 | [Novelty ledger](#part-10--novelty-ledger) | What is new, and what it is built on |
| 11 | [Expected score, risk, kill criteria](#part-11--expected-score-risk-kill-criteria) | Honest numbers |
| A | [Appendix — schemas, pseudocode, hyperparameter tables](#appendix-a--schemas-pseudocode-hyperparameters) | Copy-paste material |

---

# Part 0 — The thesis, and how this differs from `plan.md`

## 0.1 The diagnosis of the existing plan

`plan.md` is a faithful, competent restatement of what the Kaggle field converged on:
mass filter → four evidence channels → GBDT ranker → union → top 25. Its documented
progression (0.087 → 0.145 → 0.335 → 0.420) is exactly the published pipeline's
progression. That is its problem: **it is a paraphrase of the second-place solution.**

Read the score table in `docs/03` again:

```
  neural fingerprint alone ......................... 0.145
  + analog propagation ............................. 0.335   (+0.19  ← the big one)
  + two-ranker blend ................................ 0.358   (+0.02)
  + PubChem recall + richer features ............... 0.384   (+0.03)
  + engine fusion + union rescoring ................. 0.399   (+0.015)
  + popularity prior + library gate ................. 0.409   (+0.010)
  + PubChem join .................................... 0.420   (+0.011)
  ─────────────────────────────────────────────
  leaderboard top ................................... 0.47
```

Two numbers in that table define the remaining opportunity:

**Opportunity 1 — the instrument gap is 82 % of the library and nobody harmonised it.**
The test set is 100 % Bruker timsTOF. `enveda-180` (46 % of training spectra) is timsTOF
but *drug-like*, and the natural-product libraries (`riken`, `gnps`, `massbank`, `mona`,
`spectraverse`, `msdial`) that contain the *right chemistry* were all measured on Orbitrap /
QTOF. Every published solution treats this by **restriction** — they build an "AFIX" channel
from timsTOF rows and a "BIO" channel from ChEBI/LIPID MAPS, then blend with a tuned weight
(`W35`…`W80` in the notebook titles). Restriction throws away the cross-instrument evidence
instead of *converting* it. Part 4 Channel **E** converts it: a per-library collision-energy
response curve estimated from the data, used to rescale every reference spectrum into the
query's effective energy before matching. That is a strictly larger usable library than
either channel alone, at zero inference cost.

**Opportunity 2 — every scoring function in the field is one-sided, so isomers tie.**
`frag_score = (explained observed intensity) / (total observed intensity)`. It rewards a
candidate for predicting *more* peaks. Two isomers of C₈H₈O₂ that differ only in which low-mass
ion they produce will both score near 1.0 if you let the fragmenter enumerate three bonds.
Meanwhile the *false* peaks a candidate predicts cost it nothing. So the score that is meant
to break isomer ties is structurally incapable of breaking isomer ties. Part 4 Channel **F**
replaces it with a **two-sided maximum-weight bipartite assignment** (Hungarian / Sinkhorn),
which charges for predicted-but-unobserved peaks. This is a one-line change in concept and,
in our read, the single biggest unexploited lever in the whole competition.

**Opportunity 3 — nobody constructs new candidate structures from an analog hit.**
The entire "analog" channel in every published pipeline ends at *Tanimoto similarity to a
database structure*. But if the retrieved reference spectrum belongs to structure **X** and the
query is offset from it by a consistent Δ, then the query molecule is frequently **X + a small
edit** — and that edited structure is usually *already in the database* (PubChem / COCONUT),
just at a mass 14.0157 Da away. Nobody edits it into existence. Part 4 Channel **Δ**
(Delta-Structure Mining) does: retrieve → solve for Δ → map Δ to a chemical operation
(methylate / demethylate / oxygenate / glycosylate / dehydroxylate / conjugate) → apply the
operation to the retrieved structure with RDKit reaction SMARTS → **return a brand-new
candidate that no channel could have proposed.** This is the single most direct attack on
**Class 2** (structure known, no public spectra), which `docs/03` says is the second-largest
slice and which every published pipeline handles worst.

**Opportunity 4 — the de novo track is being aimed at the wrong target.**
`plan.md` Part 10 claims five novel ideas. Three of them (MCES validation, subformula-graph
scoring, MIN aggregation) are re-descriptions of published methods under new names —
Minimum Common Edge Subgraph on fragmentation trees is Böcker's fragmentation-tree alignment
(*Anal Chem* 2013), which the plan's own citation actually supports. The remaining two
(cross-track fusion, blind alphabet projection) are either already published (union-rescoring
is Phase 5 of the public pipeline) or underspecified. And the whole track is aimed at
Class 3, which `docs/03` correctly notes has < 5 % exact-recovery state of the art. Chasing
an exact novel skeleton for a few molecules is the lowest expected-value work in the project.

Plan 2 fixes this by making Track B **algorithmic at its core** (no generator, no diffusion,
no transformer) and by giving it a *second, much larger* job than "guess the novel skeleton":

1. **De novo assembly (Class 3).** A mass-conserving assembly DP + exact-cover ILP over a
   substructure alphabet inferred from data. Genuinely novel for branched small molecules.
   Honest expected value: low, ~+0.005–0.015.
2. **Counterfactual / ablation evidence (all classes).** Track B can ask *"which bonds, if
   deleted from candidate X, would produce exactly the observed peaks?"* The set of X for
   which the answer is "almost all of them" is a far sharper discriminator than any
   similarity score. This turns Track B from a parallel competitor into a **critic of
   Track A's candidates**. Expected value: an order of magnitude larger than job 1.
3. **Tail generation.** Track B supplies ranked novel structures for slots 15–25, where the
   only requirement is that they be plausible — the metric pays 1/20 ≈ 0.05 for free if the
   answer is there.

## 0.2 The ten-claim summary of what is different

| # | `plan.md` does | Plan 2 does | Why it is better |
|---|---|---|---|
| 1 | Restrict analog search to timsTOF; blend two channels with a tuned weight | **Learn a per-library CE response curve and rescale every reference spectrum into the query's energy frame** | Uses 100 % of the library instead of 18 %; no blend-weight tuning |
| 2 | One-sided `frag_score` = explained observed intensity | **Two-sided Hungarian assignment** → F-β, plus a **null-model log-likelihood ratio** | Directly discriminates isomers; score is comparable across candidates |
| 3 | Analog channel ends at Tanimoto to a database structure | **Delta-Structure Mining**: Δ → chemical operation → new structure via reaction SMARTS | Creates Class-2 candidates that literally cannot be retrieved |
| 4 | Average or MIN across a molecule's 1–16 spectra | **Replicate consensus peak voting** (peak must recur across spectra) + evidence-theoretic pooling | Removes chimeric/noise peaks *before* scoring; beats MIN which keeps noise |
| 5 | Chimeric peaks are ignored | **Sparse mixture decoding** (L¹-regularised, CHIMERYS-style) against candidate forward spectra | Handles the documented "real complex mixture" test distribution |
| 6 | Validate on 250 `enveda-np-examples` molecules | **Three-part harness**: 250 (instrument-matched) + a **synthesised Class-2 set of ≥ 3 000** built from riken/gnps + a scaffold-disjoint set | 250 structures cannot resolve ±0.02; 3 000 can |
| 7 | Rank by GBDT score, add popularity prior | **Estimate a calibrated posterior** (target–decoy + isotonic), then rank by posterior | MRR is monotone in posterior; calibration *is* the objective |
| 8 | Tail = "near neighbours, most abundant first" | **Expected-gain tail allocation** with an explicit pool-recall ceiling model | Optimises the metric directly instead of guessing |
| 9 | De novo = formula enumeration + MetFrag + MCES | **Alphabet inference → structure-space fragmentation graph → mass-conserving assembly DP → exact-cover ILP** | Actually novel; also repurposed as a *critic* with 10× the value |
| 10 | No transductive use of the test set | **Transductive spectral-network label propagation** across the 400 test molecules | Natural products come in families; relatives of a solved molecule get boosted |

## 0.3 What is deliberately *not* changed

Being better does not mean being different for its own sake. These stay, because they are
load-bearing and validated:

- 8.5 ppm mass window as the first filter (docs Part 2 §5: the defect arithmetic pins the
  formula to essentially one answer at 5 ppm).
- Spectral entropy as the default library similarity (Li et al. *Nat Methods* 2021 beat 42
  alternatives over 434 k spectra).
- Modified-cosine-style offset search — but applied to *every* library after CE rescaling,
  not to a subset.
- The popularity prior `z(rank) + μ·(log(1+substances) + log(1+pubmed))`, μ ≈ 0.15–0.25.
- Union-then-rerank instead of rerank-only.
- The library gate: never let a simulated spectrum outrank a measured one.
- Always submit exactly 25 SMILES, stereo-free, deduped on InChIKey14.
- `enveda-np-examples` as the *primary* validation slice.

---
# Part 1 — The problem as a system

## 1.1 The one picture

Read this as the map for everything that follows. Boxes are stages, arrows are data.
Track A is the left spine, Track B is the right spine, and they touch only at the Merge Gate.

```
╔══════════════════════════════════════════════════════════════════════════════════════╗
║                              ONE TEST MOLECULE                                      ║
║   molecule_id = m_0014ef                                                             ║
║   1..16 spectra · adducts {[M-H]-, [M+H]+, [M+Na]+ ...} · CE ∈ {[20],[40],[20,40]}   ║
╚══════════════════════════════════════════════════════════════════════════════════════╝
                                     │
   ══════════════════════════════════╪═════════════════════════════════════════════════
   ║  LAYER 0   SHARED FRONT-END    │                                                 ║
   ╠═════════════════════════════════╪═════════════════════════════════════════════════
   ║  0.1  Molecule assembly        │  M = precursor ∓ adduct ; reconcile ALL spectra   ║
   ║  0.2  Mass-defect formula ILP  │  MIST-style ILP + isotope + Senior/Lewis            ║
   ║  0.3  Peak curation            │  intensity floor, top-N windows, deisotope          ║
   ║  0.4  REPLICATE CONSENSUS      │  ★ peaks must recur across the molecule's spectra   ║
   ║  0.5  Query evidence vector    │  base_peak_intensity, n_peaks, CE, adduct, mode     ║
   ╚═════════════════════════════════╪═════════════════════════════════════════════════
                                     │
        ┌────────────────────────────┴─────────────────────────────┐
        │                                                          │
   ═════▼═══════════════════════════════════════          ═════════▼════════════════════
   ║   TRACK A  LIBRARY SEARCH                    ║   ║  TRACK B  UNKNOWN FROM       ║
   ║   "is it already known?"                     ║   ║  DATA ALONE                   ║
   ║   ════════════════════════════════════════   ║   ║  ═══════════════════════════   ║
   ║   A1  Formula-anchored candidate pool        ║   ║  B1  Alphabet inference       ║
   ║       (train ∪ COCONUT ∪ PubChem, laddered)  ║   ║      (Gibbs / MAP over corpus) ║
   ║                                              ║   ║  B2  Structure-space          ║
   ║   A2  L  Library match        [measured]     ║   ║      fragmentation graph      ║
   ║   A3  E  Energy-rescaled analog[measured]    ║   ║  B3  Mass-conserving         ║
   ║   A4  Δ  Delta-Structure Mine [CONSTRUCTIVE] ║   ║      assembly DP               ║
   ║   A5  F  Two-sided forward     [simulated]   ║   ║  B4  Exact-cover ILP          ║
   ║   A6  T  Tree alignment        [measured]    ║   ║  B5  Bond completion          ║
   ║   A7  S  Sparse mixture decode [simulated]   ║   ║  B6  Model selection (BIC/MA) ║
   ║   A8  N  Fingerprint posterior  [learned]    ║   ║  B7  Bond-ablation critic     ║
   ║   A9  H  Hard substructure masks[logical]    ║   ║      ★ feeds BACK into A      ║
   ║   A10 P  Popularity prior       [prior]       ║   ╚══════════════════════════════╝
   ║   A11 X  Transductive propagation[test-set]  ║
   ╚═══════════════╤══════════════════════════════╝
                   │   ~ N candidates, ~150 features
   ════════════════▼═════════════════════════════════════════════════════════════════
   ║  MERGE GATE  (Part 6)                                                           ║
   ║  1. union of A's pool and B's generated pool, dedup on InChIKey14              ║
   ║  2. EVERY candidate — library-backed or de-novo — is scored by ALL channels      ║
   ║  3. the library gate: a de-novo candidate may never outrank lib_sim ≥ 0.7        ║
   ╚═══════════════╤═════════════════════════════════════════════════════════════════
                   │
   ════════════════▼═════════════════════════════════════════════════════════════════
   ║  FUSION + CALIBRATION + DECISION                                                ║
   ║  • LightGBM lambdarank over ~150 tabular features, grouped by molecule_id       ║
   ║  • target–decoy → empirical p → isotonic → calibrated posterior p(c | spectrum)   ║
   ║  • rank strictly by posterior  (this is MRR-optimal, see §6.3)                    ║
   ║  • expected-gain tail allocation to exactly 25 slots                             ║
   ╚═══════════════╤═════════════════════════════════════════════════════════════════
                   ▼
        molecule_id,smiles  with 25 semicolon-joined, stereo-free SMILES
```

## 1.2 Why the pipeline is shaped this way

The problem has one structural property that dictates everything: **the candidate pool is
small enough to score exhaustively.** ~712 000 structures from train ∪ COCONUT, extensible to
~2 M with PubChem. At 5–8.5 ppm you keep tens to low hundreds per molecule. There is no need
to *generate* a shortlist — there is time to *score every survivor properly*. So the entire
engineering problem is not "retrieve fast"; it is "score well, and score in a way that is
comparable between candidates that have wildly different amounts of evidence."

That is why this plan spends so much on the scoring channels and so little on ANN indexes.

## 1.3 The three numbers that govern every trade-off

**(a) Precision of the precursor: 5 ppm.** At 350 Da that is 0.00175 Da. Consequences:
the molecular formula is determined by arithmetic in the overwhelming majority of cases
(docs Part 2 §5.2: for a 302.0427 target, exactly one of 11 oxygen counts gives whole `nH`
*and* whole `nC`). Consequence for the pool: one formula ⇒ a closed set of isomers, and a
mass window of ±8.5 ppm is not a heuristic — it is the physically justified limit.

**(b) Median 20 peaks per spectrum, and a molecule has a median of 3 spectra.** 60 numbers
for a molecule with, say, 3 000 constitutional isomers at its formula. The information is
genuinely under-determined (docs Part 1 §12.4, the C₈H₈O₂ isomer trap). Therefore:
**we are not trying to answer; we are trying to place 95 % of the posterior mass on 25
structures and order them correctly.** Every design decision below is subordinate to that.

**(c) MRR@25 has no penalty for wrong guesses.** This has a consequence people get wrong.
Since there is no penalty, the MRR-optimal list is *strictly* the posterior sorted descending.
There is no benefit to "hedging" with structurally diverse guesses. **Diversity is not a
feature; it is a bug.** It only enters as a tie-breaker when two candidates have identical
posteriors, where any order is equally good anyway. This is a real correction to the common
intuition and it is derived formally in §6.3.

## 1.4 Vocabulary lock (used throughout, defined once)

| Symbol | Meaning |
|---|---|
| `M` | neutral monoisotopic mass of the molecule (Da) |
| `M_ion` | mass of the precursor ion = `precursor_mz × z`; `z = ±1` throughout |
| `F` | molecular formula of the neutral, e.g. `C15H10O7` |
| `F_frag` | formula assigned to a fragment ion |
| `Δ` | a mass difference (adduct mass, neutral loss, or inter-structure offset) |
| `P` | a candidate structure (a molecular graph) |
| `S` | a spectrum: list of `(m/z, relative intensity)` pairs |
| `Q` | the query molecule = all of its spectra + metadata |
| `τ` | m/z matching tolerance (Da). `τ_spec = 0.005`, `τ_prec = 8.5e-6 × M` |
| `E` | effective collision energy after harmonisation (eV) |
| `sim_lib` | measured-library similarity, the *measurement-evidence* score |
| `sim_sim` | simulated-spectrum similarity, the *model-evidence* score |
| `Λ` | the library gate threshold (`sim_lib ≥ Λ ⇒ sim_sim is suppressed`) |
| `p(P \| Q)` | the calibrated posterior we are actually being scored on |

---
# Part 2 — Ten design laws

These are the rules that settle every later argument. When two options are being compared,
the winner is the one that obeys more of these.

---

**LAW 1 — Measurement beats simulation. Always.**
A reference spectrum was produced by a real molecule in a real instrument. A predicted
spectrum was produced by a model trained on someone else's instrument. When they disagree,
the measurement is more likely right. This is why `docs/03` §5.3 needed a *gate*, and why
the gate is not a hack but a statement of this law. Concretely: `sim_sim` is multiplied by a
gate `g(sim_lib) = 1 − σ((sim_lib − Λ)/κ)` with Λ = 0.70, κ = 0.05 — a soft version, so
candidate evidence degrades gracefully instead of switching off.

**LAW 2 — A score that is not a likelihood ratio is not comparable across candidates.**
Raw cosine, raw `frag_score`, raw Tanimoto all live on different unknown scales. When you
have 50 candidates in one list, only *comparative* ordering matters — so every channel
that feeds the ranker emits either (a) a value normalised against a **null** computed over
the candidate's own cohort (the other candidates at the same mass), or (b) an explicit
likelihood ratio. MolDiscovery's published formulation
`log P(spectrum | molecule) − log P(spectrum | null)` is the template
(*Nat Commun* 12:5955).

**LAW 3 — Scoring must be two-sided.**
A candidate must be rewarded for peaks it predicts *and observed*, and **charged for peaks
it predicts and did not**. One-sided coverage is exploitable: predict more, score more. This
is the single most consequential correction in this plan (§4.5).

**LAW 4 — Correctness is measured on InChIKey14. Therefore: strip, then work.**
Drop stereochemistry from every training label, drop tautomer specification from every
output, dedup on InChIKey14 before writing. A model that spends capacity on stereochemistry
spends it on a question the scorer has already answered for you. Working in flat (no-stereo)
SMILES space also *increases* effective training data, because every stereoisomer collapses
into one example.

**LAW 5 — Never let a design decision confine a candidate class to fixed slots.**
From `docs/03` §5.4: pinning out-of-pool candidates to positions 4, 8, 12 was worth +0.011
purely by removing the ceiling. Audit: *what in my pipeline can never reach rank 1?* If a
candidate with zero library support cannot be ranked first when the evidence demands it,
that is a bug, not a safety feature. Track B output must be able to reach rank 1 (subject
only to Law 1's gate, which is evidence-based, not class-based).

**LAW 6 — Aggregate a molecule's spectra by consensus, not by summary statistic.**
Each molecule has 1–16 spectra of the same molecule under different energies/adducts. Noise
peaks are independent across replicates; true fragment ions recur. So the correct operation
is **peak-level voting across replicates** (a consensus spectrum), not MIN-of-scores and not
mean-of-scores. MIN is defensible when n = 1–2 but throws away the majority information when
n = 16. Our rule: n = 1 → no consensus; n ≥ 2 → consensus by vote count × intensity
agreement; then aggregate scores by a **lower-tail quantile** over spectra so a contradiction
still bites.

**LAW 7 — Exploit the instrument mismatch; do not restrict around it.**
The test set is one instrument. The right library is the wrong instrument. Rather than
throwing the wrong-instrument library away (or using a tuned blend weight, which is what the
field does), **transform it into the right frame**. Law 1 tells us not to over-trust the
transform — which is why the transformation is applied to *similarity*, and the residual
uncertainty is fed to the ranker as a feature rather than being hidden.

**LAW 8 — The metric is a function of a posterior. Optimise the posterior.**
`MRR@K` is maximised by sorting candidates by `p(correct | evidence)`, descending. Therefore
the ranker should be *trained and evaluated as a probability estimator*, not as a
discriminator. Concretely: train with `lambdarank` for ordering, then **calibrate** with
target–decoy + isotonic regression, and report both MRR@25 and a calibration error (ECE).
An uncalibrated ranker is silently leaving points on the table, and calibration error is the
diagnostic that tells you so.

**LAW 9 — Evidence that is a logical fact beats evidence that is a score.**
If the query spectrum has a 3:1 M+2 isotope pair, the molecule contains chlorine. Full
stop. That is a Boolean constraint, and it should be applied as a *filter*, removing
candidates before any scoring — not as one more column in a feature matrix where a GBDT may
learn to ignore it. Concretely: chlorine/brine/halogen detection, nitrogen rule, hydrogen
cap, and diagnostic substructure presence all become **hard conjunctive masks** (§4.9).

**LAW 10 — A ceiling you do not measure is a ceiling you will not find.**
Two ceilings bound this competition and neither is visible in the loss curve:
1. **Pool recall@25** — if the answer is not in the union of the candidate pool and the
   generated pool, the score is 0 for that molecule no matter how good the ranking is.
2. **Evidence ceiling** — if the spectrum physically cannot discriminate isomer *k*, no
   ranker can. The only fix is *more evidence* (more spectra of the same molecule), never
   a better model.
Both must be measured separately (Part 7) and reported, because they determine where the
remaining effort should go.

---
# Part 3 — Layer 0: Ingestion, Curation, and the 1 GB Memory Contract

Everything in Layer 0 runs **once** and writes a set of compact, memory-mappable artefacts to
disk. Everything downstream reads those artefacts. Layer 0 is the only place that touches
`train.parquet`.

---

## 3.0 The memory contract

The box this machine gives us is about **1 GB of usable RAM**. `train.parquet` is 3.0 GB on
disk, and an in-memory `pandas.read_parquet` of its peak columns alone is ~20 GB. So:

**Rule 0.1 — Never materialise a column you do not need, and never materialise all row
groups at once.** PyArrow's `ParquetFile.iter_batches(batch_size, columns=...)` reads
row-group-at-a-time with column projection. This is the single most important operational
constraint in the project and it shapes the module layout:

```python
# data_layer/ingest.py  — the only file that opens train.parquet
import pyarrow.parquet as pq

PEAK_COLS   = ["ms2_mzs", "ms2_normalized_intensities", "precursor_mz", "adduct",
               "ionization_mode", "collision_energy_ev", "base_peak_intensity",
               "normalized_smiles", "inchikey14", "molecular_formula",
               "ingest_lib", "precursor_error_ppm", "num_peaks", "spectrum_id"]
BATCH_ROWS  = 20_000                       # ~20 k spectra × ~20 peaks × 8 B × 2 = 6.4 MB

def stream_train():                        # generator, peak memory ≈ 30 MB
    f = pq.ParquetFile("data/train.parquet")
    for batch in f.iter_batches(batch_size=BATCH_ROWS, columns=PEAK_COLS):
        yield batch.to_pydict()
```

**Rule 0.2 — Every derived structure is written to disk in a fixed-width binary format and
memory-mapped, never held in a Python dict.** See Appendix A.1 for the schemas.

**Rule 0.3 — Peak lists are stored as `float32` m/z + `float32` intensity in one packed
blob per library**, offsets in an `int64` array. `float32` gives 7 significant digits, i.e.
0.0001 Da resolution at m/z 1000 — four orders of magnitude finer than `τ_spec = 0.005`.
Intensity in `float32` is exact enough. This halves I/O versus `float64` and matters.

**Memory ledger (peak RSS per stage, target ≤ 900 MB):**

| Stage | Dominant object | Size | Running total |
|---|---|---|---|
| Train ingest (batch loop) | one `BATCH_ROWS` dict | 30 MB | 30 MB |
| Peak blob write | 2.5 M × ~22 peaks × 8 B | 440 MB **on disk**, ~0 resident | 30 MB |
| Structure index | 712 k × 256 B (2048-bit packed) | **182 MB** | 212 MB |
| Formula table | 712 k × (mass f64, idx i32, f_i14 i32) | 14 MB | 226 MB |
| Peak index (inverted) | 2.5 M spectra × 22 peaks × 8 B offsets | 440 MB **mmap** | 226 MB |
| Ranker inference | one molecule's union (≤ 400 rows × 150 feats f32) | 0.24 MB | 227 MB |
| **Headroom** | | | **~670 MB** |

The single largest resident object is the structure fingerprint table at 182 MB. If we
downsized to 1024-bit fingerprints it is 91 MB — do this if anything is tight. Never build a
Python set of 712 k SMILES (that alone is ~250 MB of strings plus ~400 MB of overhead).

---

## 3.1 Step 0.1 — Molecule assembly from spectra

**Input.** All rows sharing a `molecule_id`. **Output.** One `Molecule` object.

```
Molecule
  ├─ molecule_id        : str
  ├─ spectra[]          : Spectrum (1..16)
  ├─ M_neutral          : float            # consensus neutral monoisotopic mass
  ├─ M_spread_ppm       : float            # disagreement across spectra → quality
  ├─ formula_hypotheses : list[Formula]    # ranked, from 0.2
  ├─ consensus_spectrum : Spectrum         # from 0.4
  ├─ energy_profile     : dict[E → Spectrum]   # CE-binned merged consensus
  └─ adducts_present    : set[str]
```

### 3.1.1 Neutral mass, with the sign rule made unforgeable

The single most common bug in this task is the adduct sign. Do not encode the sign in a
conditional; encode it in the table.

```python
# data_layer/adducts.py
# sign = +1  →  adduct ADDS mass to M   →   M = mz − passenger
# sign = −1  →  adduct REMOVES mass     →   M = mz + passenger
ADDUCTS = {
  #  name             passenger      sign   mode      z
  "[M+H]+":         (1.007276,       -1,  "positive", +1),
  "[M+NH4]+":       (18.033823,      -1,  "positive", +1),
  "[M+H2O+H]+":     (19.017841,      -1,  "positive", +1),
  "[M-2H2O+H]+":    (37.028396,      -1,  "positive", +1),
  "[M+Na]+":        (22.989218,      -1,  "positive", +1),
  "[M+K]+":         (38.963158,      -1,  "positive", +1),
  "[M-H]-":         (1.007276,        +1,  "negative", -1),
  "[M-H2O-H]-":     (19.017841,       +1,  "negative", -1),
  "[M+CH2O2-H]-":   (44.998201,       +1,  "negative", -1),
  "[M-Cl]-":        (34.968853,       +1,  "negative", -1),
}

def neutral_mass(precursor_mz: float, adduct: str) -> float:
    passenger, sign, _mode, _z = ADDUCTS[adduct]
    return precursor_mz - sign * passenger      # sign=-1 → mz + passenger
```

Sanity assertions that must be asserted in code, not in a comment:

```python
assert abs(neutral_mass(301.0354, "[M-H]-") - 302.042676) < 1e-5     # quercetin
assert abs(neutral_mass(137.0597, "[M+H]+") - 136.052424) < 1e-5     # C8H8O2
assert abs(neutral_mass(203.0526, "[M+Na]+") - 180.063382) < 1e-5
```

### 3.1.2 Cross-spectrum reconciliation — free supervision, free quality control

Every spectrum of a molecule was acquired independently and must yield the same `M`. This
gives us three things for free: a corrected `M`, a per-spectrum confidence weight, and a
detector for mis-assigned rows.

```
  weights by base_peak_intensity (log scale):
      w_s = log1p(base_peak_intensity_s)  −  min_s log1p(base_peak_intensity_s)

  M_neutral = Σ_s w_s · M_s  /  Σ_s w_s            (inverse-variance-weighted mean)
  M_spread_ppm = (max_s M_s − min_s M_s) / M_neutral × 1e6
```

| `M_spread_ppm` | Interpretation | Action |
|---|---|---|
| < 2 | fully consistent | keep all spectra |
| 2–15 | one adduct mislabelled, or genuine in-source fragment | drop the *outlier* spectra, recompute |
| > 15 | mis-assigned molecule_id, or a fragmenting precursor cluster | flag the molecule, run 0.2 with the majority mass and audit |

**In-source fragmentation** deserves its own check, because it is a real and common artefact:
timsTOF in-source fragmentation produces a second ion at `M − nH₂O` that is isolated and
fragmented, and its `adduct` column can legitimately read `[M-H2O-H]-`. So a 18.0106 Da
disagreement between two spectra of the same molecule is **expected**, not an error. Our
rule: if `|M_a − M_b|` is within 3 ppm of `n × 18.010565` for `n ∈ {1,2,3}` (or `n × 17.026549`
for ammonia, or `n × 27.994915` for CO), we collapse the family onto the **largest** `M` and
keep all spectra, flagging the molecule `insource_family`. This alone recovers spectra that a
naive reconciler would throw away.

### 3.1.3 Formula from mass — the ILP, done properly

`docs/02` §5 shows the hand loop works because at 5 ppm exactly one `(nO, nH, nC)` triple
survives. That loop is O(nO_max) and assumes CHON. Three things make it better:

**(i) The integer program (Pluskal/Böcker/SIRIUS form).** Minimise
`Σ_e |n_e − n*_e|` subject to exact mass equation:

```
  minimise   Σ_e  w_e · |n_e − n*_e|                       n_e ∈ ℤ≥0, e ∈ {C,H,O,N,S,P,Cl,Br}
  subject to Σ_e m_e·n_e = M_expanded                      (M·blowup, exact integer equality)
             0 ≤ n_H ≤ 2n_C + 2 + n_N + n_P + n_S·0        (RDBE ≥ 0)
             n_C ≥ 1,  n_O ≤ 3n_C,  n_N ≤ 4n_C,  n_S ≤ 3n_C, n_P ≤ 2n_C
             (n_H + n_X + 3n_S + 5n_P) mod 2 = 0           (LEWIS)
             Σ valence ≥ 2·(n_atoms − 1)                    (SENIOR)
```

The **blowup** trick: scale the mass equation by `b` so that a mass like 302.0427 becomes the
integer `3020427` at `b = 1e4`, removing all floating-point error from the constraint. With
element masses given to 7 digits and `M` known to 5 ppm, `b = 1e4` is safe.

**(ii) The isotope-pattern filter as a hard constraint.** The precursor's own isotope
envelope is a *measured* second channel. Simulate the envelope of each candidate formula
(convolution of per-element isotope distributions) and reject formulas whose `(M+1)/M` ratio
disagrees with observation. For CHON at C~20 the (M+1) ratio is ≈ 1.1 % × nC, so this
separates C₁₅H₁₀O₇ from C₁₄H₂₄O₈ — same nominal mass region, different carbon count, both
passing the mass test. Concretely, the Seven Golden Rules report removes **668 of 1 700**
candidates at 5 % isotope accuracy (*BMC Bioinformatics* 8:105). Our test spectra are timsTOF
and low-intensity; we make this a *soft* score with a wide band rather than a hard cut, and
expose it as a ranker feature.

**(iii) The Formula Subset Analysis / MaxComb=1 rule — the cheap and very effective one.**
Every observed fragment peak must have **exactly one** candidate formula consistent with the
root formula. Peaks with 0 consistent formulas are contamination and are discarded (this is
also a chimeric-peak detector — see §4.7). Peaks with ≥ 2 are fine individually but the
*number* of such ambiguous peaks is itself diagnostic: a good root formula explains the
spectrum with low ambiguity. Score by

```
  score(F) = − Σ_peak  log( k_F(peak) + 1 )
```

where `k_F(peak)` is the count of sub-formulas of `F` matching that peak within `τ_spec`.
This is derived from *Formula Subset Analysis* (*Metabolomics* 2025,
10.1007/s11306-025-02379-0) and costs microseconds. It is a genuinely independent second
opinion on the formula, since it uses the fragment masses that the precursor mass never saw.

**Result: `formula_hypotheses` is a ranked list, not a single answer.** In practice `F₀` is
correct for ~92 % of molecules at 8.5 ppm; we keep up to 8 hypotheses and let the pool
builder union across them, because a wrong `F₀` is the single most damaging failure mode
in the whole system (it produces a candidate list that cannot contain the answer).

---

## 3.2 Step 0.2 — Peak curation

Applied **per spectrum**, and every choice is logged so it can be ablated.

| Step | Rule | Rationale |
|---|---|---|
| **Relative-intensity floor** | drop `intensity < 0.005` (0.5 %) | below this, hit counts are dominated by noise on timsTOF |
| **Above-precursor removal** | drop `mz > precursor_mz − 0.5` | already done by organisers, but *merged* spectra may bypass it — re-apply |
| **Windowed top-N** | per 50 Da window keep top 6 **and** globally keep top 128 | the classic library-search variant; preserves low-mass ions that a global top-N discards, and low mass is exactly where isomers differ (docs Part 1 §12.4) |
| **Deisotope** | merge peaks at `+1.003355` when `I(M+1)/I(M) < 0.15` and `I(M+1) < 0.1` | only merge when the satellite is *smaller than a natural ¹³C envelope would allow*; blindly merging destroys the carbon-count evidence we want in §4.4 |
| **Intensity transform** | keep raw, and additionally store `sqrt(I)` | sqrt is the standard compression; keeping both lets the ranker choose |
| **Minimum peaks** | drop spectra with `< 4` peaks after curation | they carry no information and only add noise to consensus voting |

**Peaks-above-precursor deserve a footnote.** A fragment cannot outweigh a singly-charged
precursor — but for `[M+Na]⁺` and `[M+K]⁺` in-source clusters, a peak at `M + 22` above the
selected adduct is a classic artefact. We remove anything above the *largest* plausible ion of
the molecule, `max(precursor_mz + 21.98, precursor_mz)` for alkali adducts.

---

## 3.3 Step 0.3 — REPLICATE CONSENSUS *(replaces "aggregate with MIN")*

**The idea.** A molecule has 1–16 spectra. True fragment ions of a molecule recur across
spectra; noise and co-isolated contaminants do not. This is exactly the argument for
consensus spectrum construction in proteomics (Beavis & Craig) and the argument DNMS2Purifier
makes for chimeric cleaning in metabolomics (*Anal Chem* 2023, 10.1021/acs.analchem.3c00736):
true fragments are correlated across replicates, false ones are not.

**Algorithm.**

```
input : spectra S_1..S_n of one molecule, each with CEs E_1..E_n

# 1. bucket by energy — do not vote across incompatible fragmentation regimes
buckets = group_by(round_up_to_10eV(E_s))        # 20eV bucket, 40eV bucket, ...

# 2. within a bucket, vote on m/z
for each bucket B with |B| >= 2:
    votes = {}
    for s in B:
        for peak p in s (intensity-weighted):
            votes[ match_within_tau(p) ] += p.intensity * sqrt(BPI_s)
    for each mz-cluster c in votes:
        n_appear   = number of spectra in B containing a peak in c
        if n_appear / |B| < 0.34:            # appeared in < 1/3 of replicates
            c.discarded                       # ← chimeric / noise, gone before scoring
        else:
            c.i_mean  = mean intensity of appearing peaks
            c.i_cv    = stdev / mean           # coherence
            if c.i_cv > 0.85 and n_appear >= 3:
                c.discarded                   # incoherent intensity ⇒ interference
            else:
                c.consensus_intensity = c.i_mean * min(1, n_appear/|B|)
            c.replicate_fraction = n_appear / |B|

# 3. the consensus spectrum of the molecule
consensus_peaks = { c : c.replicate_fraction >= 0.34 }
replicate_fraction is RETAINED as a per-peak weight (replicate_fraction^0.5)
and used as the per-peak `w` in every assignment problem in Part 4.
```

**Why `w = replicate_fraction^0.5` in the matching weight.** A peak seen in 1/3 of replicates
carries roughly the information content of one observation, not of three. The square root
interpolates between "seen once" (weight ≈ 0.58) and "seen always" (weight 1.0) without the
harsh cliff of a hard binary filter, which matters because for n = 2 the binary rule throws
away half the data.

**The single-replicate case.** When `n = 1` there is nothing to vote on. We then fall back to a
**prior over noise**: the empirical distribution of peak appearance rate across *all* single
spectra in `enveda-np-examples` gives `P(peak | replicate_fraction < 0.34)`. Peaks in the
tail of that distribution get down-weighted. This is a small learned prior, but it is a prior
over a nuisance variable, not a model of chemistry — Law-compliant.

**The zero-replicate-lifetime edge.** 1–16 spectra with median 3 means roughly half the test
molecules have ≤ 2. So this mechanism fires for ~50 % of the test set, and for the other 50 %
we fall back to the intensity-based noise prior. Both paths are implemented; both are cheap.

---

## 3.4 Step 0.4 — The query evidence vector

Everything constant about the molecule, consumed by every channel as features:

```
  n_spectra                     # 1..16
  base_peak_intensity (median, log10)      # spectrum quality
  base_peak_intensity (min, log10)         # weakest link
  ion mode (pos/neg)         adduct set size
  CE list (mean, spread)                   # eV
  M_neutral, M_spread_ppm
  n_peaks (median), peak density per 100 Da
  replicate_fraction histogram summary     # mean, max, fraction of mass with r>0.34
  chimeric_fraction  = fraction of raw intensity in peaks with r < 0.34
  isotope: observed (M+1)/M of the consensus, (M+2)/M ratio, per-peak Cl/Br/S indicators
  inferred halogen flags: has_Cl, has_Br, has_S  (Law 9)
```

---
# Part 4 — TRACK A: Library Search

**The question Track A answers:** *given only spectra, is the molecule already known to
science?* Track A's inputs are the query spectra plus a **structure database**. It never
consults Track B. Its outputs are (a) a candidate set, (b) ~150 features per candidate.

---

## 4.1 A1 — Candidate pool construction (the recall floor)

### 4.1.1 The pool ladder

A single 8.5 ppm window is right for the common case and wrong at the edges. Use a **ladder**
and keep every rung's output, with a rung-membership feature so the ranker can learn what
rung-3 membership costs:

```
 RUNG 0   the base pool
         S₀ = { all structures in train.parquet }              ~275 810 (dedup by inchikey14)
         S₀ ∪ COCONUT  (natural products)                       ~712 000 total
         S₀ ∪ COCONUT ∪ NPAtlas ∪ LOTUS  ∪ PubChem-CID-with-SMILES
                                                              ~2 000 000 total

 RUNG 1   |M(P) − M_neutral| ≤ 8.5 ppm                         → ~40–300 structures
 RUNG 2   |M(P) − M_neutral| ≤ 25 ppm                          → ~150–1 200
 RUNG 3   |M(P) − M_neutral| ≤ 100 ppm                         → ~600–5 000
 RUNG 4   adduct-ambiguity ladder: the recorded adduct may be wrong.
          Search M ± (adduct_i − adduct_j) for all 10×10 adduct pairs → 20 extra windows
 RUNG 5   formula-neighbourhood: same F₀ ± 1 C, ± 1 O, ± 1 N, ± 1 H2  → the "near misses"
 RUNG 6   Δ-chains: M − n·(C6H10O5 = 162.0528) for n=1..3, −n·(C6H12O5=176.0943),
          −n·(CH2 = 14.0157) for n=1..4, −n·(O = 15.9949) for n=1..3
          (de-glycosylation and methylation ladders are how natural products actually relate)
```

Rungs 4–6 are cheap because they are the *same* sorted mass array, offset. One binary-search
per rung over a sorted `float64` array of 2 M masses is ~21 µs.

### 4.1.2 Structure deduplication is subtler than it looks

`docs/dataset_desc` warns structure counts sum to >275 810 for ~275 k unique because
compounds recur across libraries. Deduplicating on InChIKey14 alone loses information: the
same skeleton measured on a timsTOF and on an Orbitrap is *two independent confirmations* of
structure identity. So:

- **Candidate identity** = InChIKey14 (this is what the scorer uses, and what dedups the
  final 25 slots).
- **Evidence** is aggregated over *all* spectra of *all* spellings of that skeleton, with
  `sim_max` over spectra, `sim_mean` over spectra, and `n_independent_libs` as a feature.

This is the difference between collapsing duplicate evidence and throwing it away, and it is
a free +small: the same skeleton in 6 libraries is far likelier to be the right answer than one
seen once.

### 4.1.3 The mass filter's real job: it is a formula filter

At 8.5 ppm on 350 Da the window is ±0.003 Da. The mass defect table in `docs/02` §2.2 makes
it arithmetic that a formula is essentially determined. Therefore **the pool can be
intersected with the formula hypotheses before any structural search**, which for
`C15H10O7` removes the ~15 % of same-mass structures with different formulas (e.g. a
chlorinated or sulfurised congener). Combine: pool = (rung ∩ formula set ∩ ¬hard masks).

---

## 4.2 A2 — Channel L: direct library match (the *measurement* channel)

**Purpose.** Detect Class 1. This is the highest-precision channel and it must never be
overruled.

### 4.2.1 Candidate spectra gathering

For candidate skeleton `P`, gather every training spectrum whose `inchikey14 == P`:
`{ (s, lib, CE, mode, adduct) }`. Typical size 1–20; the long tail matters.

### 4.2.2 Per-spectrum score: entropy with a coverage term

Spectral entropy similarity (Li et al., *Nat Methods* 2021) — the single best-validated score
in the field, beating 42 alternatives on 434 287 spectra:

```
  given spectra A (query) and B (reference):
    pA = normalise(A.intensity),  pB = normalise(B.intensity)
    mix = (pA + pB) / 2, with peaks within τ merged (max) and renormalised
    S(·)  = − Σ p_i log p_i                                (Shannon)
    Sim_ent(A,B) = 1 − 2·( S(mix) − S(A) − S(B) ) / ln 4    ∈ [0,1]
```

Two augmentations, both from the literature and both important here:

**(a) The "dynamically weighted" variant** (Li et al.). For low-entropy spectra (a few peaks
with one dominant), the plain entropy saturates. Transform with
`w = 1 if S ≥ 3 else 0.25 + 0.25·S`, apply as an exponent, renormalise. Low-entropy query
spectra are exactly our case — median 20 peaks with one base peak — so this matters.

**(b) Coverage, as a separate number.** Entropy alone is a *distribution* comparison; it does
not tell you how much of the query was accounted for. So we also compute, from the same peak
assignment used in §4.5:

```
  coverage(P, s) = Σ_{matched query peaks} w_q · I_q  /  Σ_{all query peaks} w_q · I_q
  n_matched(P,s) = count of matched query peaks
```

### 4.2.3 Aggregation over a candidate's spectra

The correct aggregator is **max-plus with a soft count bonus**, not `max`:

```
  lib_max(P)   = max_s Sim_ent(Q, s)
  lib_top2(P)  = second largest
  lib_margin(P)= lib_max(P) − max_{P' in cohort, P'≠P} lib_max(P')     # cohort = same mass window
  n_lib(P)     = |{s : Sim_ent(Q,s) ≥ 0.6}|
  sim_max_E    = max_s Sim_ent(Q, s) · exp(−|E_eff(s) − E_eff(Q)| / τ_E)   # §4.3 energy agreement
```

`n_lib` and `lib_margin` are the two features that turn a similarity into evidence. A single
0.82 similarity is weak (the isomer trap); three independent spectra from three libraries all
at 0.85, while the runner-up reaches 0.61, is strong.

### 4.2.4 Multi-spectrum combination (Law 6, applied to Channel L)

```
  For the candidate's best matching reference spectrum per query spectrum:
     take the LOW QUANTILE (0.2) of per-spectrum Sim_ent across the query's spectra
     → a_20(P)
  and the MEAN
     → a_mean(P)
  lib_score(P) = 0.6·a_20 + 0.4·a_mean
```

The quantile term implements "one contradicting spectrum disqualifies" from `docs/02` §15.2;
the mean term prevents a single noisy spectrum from killing a good candidate. The mixture is
strictly better than either pure-min (too brittle at n=16) or pure-mean (too permissive).

---

## 4.3 A3 — Channel E: energy-harmonised analog search *(the "converting, not restricting" idea)*

**This is the single largest change of direction in Track A.**

### 4.3.1 The problem, stated exactly

```
  query        : timsTOF, 35 eV, [M−H]⁻ ,  a flavone glycoside
  library      : Orbitrap (pluskal), NCE 35 → converted to "35 eV" by the nominal
                 Thermo formula NCE × mz/500
  reality      : those are different physical conditions. On a Q-Orbitrap, NCE 35 at
                 m/z 350 lands roughly in the 22–26 eV collision-energy window; on a
                 Bruker QTOF the *same* fragmentation needs ~40 eV. The nominal
                 conversion is wrong by up to ~1.6× in effective energy.
  consequence  : a correct structure at the wrong energy scores 0.45–0.6 instead of 0.9.
```

Measured across instrument types, only **64–89 %** of spectrum pairs exceed 60 % identity
(*Scheubert, Hufsky & Böcker* 2013). That is the size of the gap we are attacking.

### 4.3.2 The response-curve model

Fragmentation is controlled by *collisional activation energy per degree of freedom*. Define
for a spectrum a scalar **effective fragmentation energy**

```
  E_eff(s) = NCE_s · (m/z_precursor_s / 500) · η_lib          if units == "NCE"
           = E_s                                              if units == "eV"
           = E_s                                              if units == "V"   (× scale)
```

The unknown `η_lib` is a per-library, per-instrument-family gain. **We estimate it from the
data**, using the only place where ground truth exists: `enveda-np-examples`, 250 natural
products measured on the *test instrument* by the *test pipeline*.

```
  ESTIMATION (offline, once)

  For each of the 250 np-example structures P:
     for each pair (reference spectrum s from library L, test spectrum t on timsTOF):
        Δm  = M(P) − M(precursor(s))            # exact, so adducts are aligned automatically
        apply modified-cosine with that offset  → sim_0 = Sim_noshift, sim_Δ = Sim_shifted
        if sim_Δ ≥ 0.75 and sim_0 < 0.6 :        # unambiguous "same molecule" evidence
            this pair is a matched (s, t)
            fit η_L so that argmax_E E_eff(s; η_L) matches argmax_E argmax_Δ the local
            optimum of the sim_Δ(E_L, E_t) response surface
  → η per (library, ionisation mode), solved by coordinate ascent on the surface
     "argmax over η of the number of matched pairs whose CE falls within ±5 eV"

  The output is a small table, e.g.
     pluskal_ms2   NCE  +  η = 1.62      riken     eV  η = 0.94
     gnps          eV   η = 1.00  (no correction)   mona  NCE  η = 1.71
```

**Fallback, if the estimation is noisy (and it will be, with 250 molecules):** fit a
*global monotone map* `E_eff = a · E_nominal + b` with `(a, b)` constrained so that
`timsTOF-eV` is the identity (anchoring the test frame), estimated jointly. Regularise `a`
toward the published Rondeau collisional-efficiency ratios (QTOF-Ar ≈ 0.30, Q-Exactive ≈
0.70, Q-Orbitrap ≈ 0.50, N₂ ≈ 0.10) as a prior. This is a two-parameter fit — it cannot
memorise the validation set.

### 4.3.3 The search

```
  For each query spectrum t with effective energy E_t:
     for each reference spectrum s (whole library) with effective energy E_s:
         offset-search tolerance is WIDENED by the CE mismatch:
             τ_Δ = τ_spec + c_Δ · |E_s − E_t|          c_Δ ≈ 0.02 Da/eV
         score  = Sim_entropy_with_offset(Q_t, s, Δ = M(Q) − M(precursor_s), τ_Δ)
                  × exp( −|E_s − E_t| / 20eV )
```

The two terms are honest about their uncertainty: the offset tolerance genuinely widens with
energy mismatch (high energy fragments more, and the *set* of stable fragments changes), and
the second factor is a soft prior that a 60 eV Orbitrap spectrum is a mediocre match to a
30 eV timsTOF query, even after correction. This is Law 7 — *convert, then discount the
residual explicitly as a feature* — rather than hiding the correction.

**Cost.** Same as a full offset library search: one pass over 2.5 M spectra with an inverted
peak index. ~15 s on 400 molecules with 8 workers. Nothing.

**Why this beats the field's blend.** The field's `AFIX` channel uses only timsTOF rows
(~18 % of the library, wrong chemistry). Our channel uses **100 %** of the library after
correction, and the ranker additionally receives `energy_gap_min`, `energy_gap_mean` and
`n_ref_within_5eV` as features, so it can learn *how much to trust* the correction instead of
having us guess a blend weight `W##` in advance.

---

## 4.4 A4 — Channel Δ: Delta-Structure Mining *(the constructive channel — genuinely new)*

### 4.4.1 What the field does, and what it throws away

Every published "analog" channel ends here:

```
   query  ──modified cosine──▶  retrieved spectra {s_1 … s_k}
                                     │
                                     ▼
                       Tanimoto(candidate P, owner(s_i))          ← a NUMBER
                                     │
                                     ▼
                              ranker
```

The retrieved spectrum `s` tells you more than "structure `owner(s)` is similar." It tells
you **where the peaks are shifted**. That vector of shifts is a chemical edit specification,
and the edited molecule is very likely sitting in the same database at a mass `Δ` away. The
field computes the similarity and discards the shift.

### 4.4.2 The algorithm

```
 STEP 1  Retrieve with explicit shift tracking
   For each retrieved pair (t_query, s_ref) with precursor offset Δp:
       for each peak match (i → j) with offset δ_ij = t_i − s_j:
           accumulate  δ_ij into a shift multiset D(Δp) with weights w_ij
                       (w = sqrt(I_t · I_s) × replicate_fraction)

 STEP 2  Score shift-coherence
   For a real cognate pair, MOST matched peaks share ONE offset δ*, and those peaks are the
   ones that CONTAIN the modified region. So:
       c(δ) = Σ_{δ_ij = δ} w_ij  −  λ · Σ_{δ_ij ≠ δ} w_ij
   δ* = argmax_δ c(δ)          # the "cognate offset"
   coherence(P) = c(δ*) / Σ_{all matches} w_ij          ∈ [0,1]

   IMPORTANT known failure: with no true compound in the pool the top hit is still the
   molecule whose Δ is closest, so scores are inflated by "cognateness". The cure is to
   score the DISTRIBUTION of δ over the hit list, not each hit alone:

       δ_spread = weighted stdev of {Δp_i} over the top-k retrieved spectra
       (a real homolog family gives a small spread; random hits give a scattered spread)
   and to keep coherence and δ_spread as SEPARATE features.

 STEP 3  Map δ* to a chemical operation
   δ* is a *neutral* mass difference. Look it up in an EDIT DICTIONARY built from the
   training corpus itself (no hand-listing needed):

     δ*  →  formula-level edit   →  SMARTS reaction on owner(s)'s structure
     ─────────────────────────────────────────────────────────────────────────
     +15.9949  ΔO                 add/remove a hydroxyl
     +14.0157  CH2                methylate (O-, N-, C-methylation)
     +162.0528 C6H10O5           O-glycosylate by hexose
     +164.0685 C6H12O4           add a deoxyhexose (rhamnose)
     +203.0794 C6H11O5+H2O       hexose + water
     -18.0106  H2O                dehydrate
     +79.9568  Br  − H            brominate
     +27.9949  CO                 (usually a fragmentation, not an edit — down-weight)
     ... dictionary mined by: for all P, P' in train with |M−M'| < 0.02 and
         Morgan-4 Tanimoto ≥ 0.6, record (ΔM, Tanimoto, shared-SMARTS-ratio).
         Keep edits with ≥ 50 observed instances. ≈ 200–400 edits survive.

 STEP 4  Apply the operation → NEW CANDIDATE
   rxn = ReactionFromSmarts(dictionary[δ*])
   for attachment_site in enumerate_matching_sites(owner_structure, rxn.left):
        new = rxn.rxn.Apply(attach(attachment_site, owner_structure))
        add new to the candidate pool with tag {origin: "DSM", parent: owner, delta: δ*}
   Typical yield: 3–40 new structures per retrieved reference.

 STEP 5  Hand them straight to the normal scorer
   A DSM candidate has NO reference spectrum of its own, so channels L and E are zero for it.
   But channels F (forward fragmentation), T (tree alignment), N (fingerprint posterior) and
   the sparse decoder all work on it unchanged. That is the whole point: DSM *manufactures*
   the Class-2 candidate that retrieval cannot reach.
```

### 4.4.3 Why this is the right lever for Class 2

Class 2 is "structure is in PubChem/COCONUT but no public spectra exist." Every channel that
requires a *measured* reference is dead for it. Every channel that scores a *structure* against
a *spectrum* is alive. DSM is the only channel that **expands the structure pool** rather than
filtering it, and it does so with a chemically-justified, data-mined edit dictionary rather
than by enumeration. The PubChem-join phase in `docs/03` (§5.4, +0.011) bought the same thing
by brute force at 5 000-candidate depth; DSM gets there with ~30 high-precision candidates
per retrieved reference.

### 4.4.4 Cost control

- Cap: top 30 retrieved references per query molecule, top 3 `δ*` values each, top 8 sites
  each ⇒ ≤ 720 new structures per molecule, deduped down to ~150–400.
- All of these run *after* the base pool is scored, only for the top ~50 base candidates,
  and only if the top base candidate's `lib_max < 0.75` (i.e. we are probably in Class 2).
  **The gate is evidence-based, so DSM output is still allowed to reach rank 1** (Law 5).

---
## 4.5 A5 — Channel F: two-sided forward fragmentation *(Law 3 made concrete)*

### 4.5.1 The critique, precisely

`docs/03` §4 defines the published fragmentation channel:

```
  frag_score = ( Σ intensity of observed peaks the candidate can produce )
               / ( Σ intensity of all observed peaks )
```

Three defects:

1. **Asymmetric.** Predicted peaks that are absent cost nothing. Adding a third bond-cleavage
   level can only raise the score. (MetFrag's own ablation confirms the trap: raising tree
   depth 2 → 3 costs 5× runtime and *lowers* median rank 39 → 18 and PPV 0.028 → 0.017,
   *J Cheminform* 2016.)
2. **Unnormalised across candidates.** A candidate whose fragmenter happens to enumerate more
   distinct masses has more chances to hit. Directly violates Law 2.
3. **Reward-hackable precisely where it is needed most.** Two isomers of C₈H₈O₂ both produce
   the benzoyl ion at 105; one produces 91, the other 95. One-sided coverage scores the
   candidate that happens to enumerate more small ions higher — arbitrarily.

### 4.5.2 The fix: a two-sided assignment problem

Model the comparison as a **maximum-weight bipartite matching** between
`predicted peaks` × `observed peaks`:

```
  Sets
    Q = {(i, mz_i, I_i, r_i)}   observed peaks of the query
                                 r_i = replicate_fraction (Law 6), w_i = r_i^0.5 · I_i
    G(P) = {(j, mz_j, Ĩ_j)}     predicted peaks of candidate P (see 4.5.3)

  Compatibility      |mz_i − mz_j| ≤ τ_spec                      (τ_spec = 0.005 Da)
  Affinity           a_ij = (Ĩ_j / max_k Ĩ_k) · (I_i / max_l I_l) ^ η

  Solve              max over one-to-one assignments σ:  Σ_{(i→σ(i))} a_ij · λ_cand(i)
  Subject to         each observed peak matched ≤ 1 time
                     each predicted peak used   ≤ 1 time

  ⇒ this is exactly the Linear Sum Assignment Problem, solved exactly by
    Hungarian / Jonker–Volgenant in O(n³) with n ≤ ~128.  ~1 ms per candidate in numpy.
```

Two-sided score:

```
  recall_r    = Σ_matched w_i / Σ_i w_i                        # observed explained
  precision_r = Σ_matched a_ij / Σ_{j ∈ G(P)} Ĩ_j            # predicted NOT over-predicted
  F_β         = (1+β²)·recall_r·precision_r / (β²·recall_r + precision_r)      β = 1.0
```

**Why this beats coverage on the isomer trap.** Methyl benzoate's fragmenter predicts 105, 77,
51, 39 and also 51, 65, 79 as 2-bond products. If the truth is benzyl formate (observed: 91,
79, 65, 105) then methyl benzoate's predicted 51 and 39 are **unobserved** and now cost
precision. Conversely phenyl acetate fails to predict 91 and its recall drops. The score now
separates isomers, which the published channel provably cannot do.

### 4.5.3 Getting the predicted peak set — three tiers, in preference order

```
 TIER 1  Physical: bond-dissociation fragmenter (MetFrag algorithm, reimplemented)
   depth-1  bridges: every bond whose removal disconnects the graph
                     (Hopcroft–Tarjan / Tarjan bridge-finding, O(V+E))
   depth-1  2-cuts : pairs of ring bonds whose joint removal disconnects
                     (DEREPLICATOR+ and MolDiscovery both show 2-cuts carry real signal)
   depth-2  bitwise AND of a depth-1 incidence vector with a depth-1 fragment's
            (MolDiscovery's core trick — depth-n fragments are free)
   ring-opening rearrangements as explicit rule moves: −H₂O, −HF, −HCN, −CO, −H₂,
            −H₂O+CO, ketene, McLafferty (+43 "McLafferty rearrangement" pseudo-ion)
   intensity model:  MetFrag's published form
            S = Σ_p  RelMass_p^α · RelInt_p^β / (Σ_{b∈B_f} BDE_b)^γ
                α = 1.84,  β = 0.59,  γ = 0.47
   energy model:     gate each predicted fragment by a survival curve
            s(E, E₀) = 1 / (1 + exp(−(E − E₀)/6eV)),  E₀ = atom-count-dependent onset
            (low-energy spectra lose large fragments, high-energy lose the precursor —
             this is the textbook CE/fragment-count anticorrelation, and it is what makes
             energy-matched matching better than energy-blind)

 TIER 2  Probabilistic: CFM-style competitive break model
   θ_ij = wᵀ Φ_ij → ρ(f_i→f_j) = exp(θ_ij)/(1+Σ_k exp(θ_ik))
   Breaks the independence assumption that MolDiscovery admits to violate, because
   peaks compete for a fixed intensity budget. Trained offline on enveda-180 + np-examples.

 TIER 3  Learned forward model (ICEBERG / GLACIER), conditioned on E and ion mode
   Used ONLY for the top ~80 candidates per molecule, because it is the expensive one.
   Library-gated per Law 1.
```

### 4.5.4 Turning the score into a likelihood ratio (Law 2)

```
  raw(P)  = log F_β(P)                       # the two-sided assignment score
  null(P) = log F_β(P_null)                  # where P_null = a mass-matched, formula-matched
                                             # random structure from the same rung-1 window
  score_F(P) = raw(P) − null(P)              # the likelihood-ratio form
```

Better: use the **cohort null** — the median `log F_β` over all ~100 rung-1 candidates of this
molecule. Then `score_F(P) = log F_β(P) − median_Q(log F_β)`, which is automatically
scale-free and lets the ranker compare a 3-candidate molecule with a 900-candidate one.
Both variants are emitted as separate features.

### 4.5.5 Extra features this channel yields for free

- `frag_n_explained` — number of observed peaks above 2 % intensity explained.
- `frag_basepeak_explained` — is the *base peak* explained? A single Boolean, and in `docs/02`
  §9.2 it is the distinction between "likely correct" and "possibly correct."
- `frag_ladder` — the longest chain `P → P−L → P−L−L′` of observed peaks matched by
  successive cleavages. Sequential losses are far more specific than isolated ones
  (`docs/02` §7.3).
- `frag_gap_only` — the published "gap fragment" trick: credit only peaks that the
  *learned* forward model cannot explain, i.e. the signal a learned model misses.
  (`docs/03` Phase 6 shipped this and it was worth including in the prior.)

---

## 4.6 A6 — Channel T: fragmentation-tree alignment *(the graph-isomorphism channel)*

**Purpose.** Peak-overlap similarity throws away *topology*. Two spectra can share 18 peaks and
have completely different fragmentation logic. Tree alignment does not.

**Algorithm (Rasche/Dührkop/Böcker fragmentation-tree alignment, reimplemented).**

```
 1. Build the query's fragmentation tree FT_Q
      nodes = { one per (peak, candidate formula) }            # from §3.1.3
      edges = sub-formula relations (parent ⊃ child), weight = learned loss prior
      solve MAXIMUM COLOURFUL SUBTREE (each colour = one observed peak, used at most once)
      → MAP tree. SIRIUS 3's Bayesian FT score reaches 76.0 % correct formulas vs 31.1 %
        for the original DP formulation (J Cheminform 8:5, 2016).
 2. Build FT_C for each candidate C by running the same algorithm on the candidate's own
      *predicted* spectrum (tier-1 fragmenter) — so both sides are trees in formula space.
 3. Align FT_Q onto FT_C by maximum-weight tree alignment
      (node–node iff formulas equal; node–gap allowed at penalty; edges must be preserved)
 4. score_T = aligned_weight / total_weight_of_FT_Q
      + |edges of FT_Q NOT reproducible in FT_C|  → penalty
```

**Why this is a genuinely independent channel.** It uses *relative* mass relationships
between peaks, not their absolute masses. It therefore transfers across adducts, across ion
modes, and across instruments with very different m/z calibration — which makes it the channel
most robust to the instrument mismatch that plagues Channel L.

**Cited strength.** On a 370-compound MassBank set where spectra share almost no peaks,
FT-alignment search quality approaches MACCS Tanimoto on true structures (*Anal Chem* 2013,
10.1021/ac300304u); in the same family, `FT alignment ≈ MACCS − 0.15`, and the gap between raw
spectral comparison and FT alignment is stated to be about as large as the gap between FT
alignment and having perfect knowledge of the structure.

**Cost.** Step 1 is the expensive one (SIRIUS 4 is ~36× faster than SIRIUS 3 via
"speedy colorful subtrees" with a colour-frequency heuristic). Run it once per molecule, cache
`FT_Q`, then alignment against each candidate is a small tree DP. Feasible inside the budget.

---

## 4.7 A7 — Channel S: sparse mixture decoding *(chimeric peaks, taken seriously)*

### 4.7.1 Why this channel exists at all

`docs/dataset_desc` says the test spectra have had "only minimal cleaning applied, because
there is no standard way to curate mass spectra and we would rather not discard information
on your behalf," and describes a test distribution that mimics "a real complex mixture." A
real complex mixture means **co-isolated precursors**, i.e. chimeric spectra where some
peaks belong to a different molecule entirely. Enveda's own timsTOF data is collected in
screening mode, which is the most co-isolation-prone acquisition mode in existence.

This has two consequences nobody in the Kaggle field is exploiting:

1. Unmodelled contaminant peaks depress the *recall* side of every forward score, uniformly,
   for the true candidate — so the true candidate's absolute score is depressed and the
   ranking suffers.
2. A contaminant peak can, at low mass, coincidentally match a predicted fragment of a wrong
   candidate — manufacturing false evidence.

### 4.7.2 The model

The chimeric spectrum is (approximately) a **non-negative mixture of pure spectra**:

```
      I_obs(m) ≈ Σ_k α_k · Î_k(m)      with α_k ≥ 0
```

CHIMERYS (*Nat Methods* 2025) solves exactly this with **L¹-regularised regression**: explain
as much observed intensity as possible using as few precursors as possible. Port the
objective, drop the peptides, keep the mathematics:

```
  MINIMISE over α ∈ ℝ⁺_{|union|},  and a candidate support C:
      (1/2)·‖  Ā α − Ī_obs ‖²₂
    + λ₁ · ‖α‖₁                       ← sparsity: few molecules explain the spectrum
    + λ₂ · ‖ α − α_target ‖₂²          ← anchor to the molecule we are evaluating,
                                          so we ask "does P explain it well?"
                                          rather than "does anything explain it?"
  subject to  Σ_k α_k = 1             (each spectrum is 100 % from some combination)

  The dual variables of λ₁ are exactly the per-candidate *contribution* scores:
      contribution(P) = residual reduction when P is added to the support
  This is the natural "how much of this spectrum is P?" number.
```

Concretely, and this is the practical form we implement:

```
  UNION   = top-60 candidates of the base pool ∪ DSM candidates
  A       = [ predicted intensity matrix ]  (|union| rows × |peaks| columns, sparse,
                                                from the tier-1 fragmenter, energy-gated)
  Ī_obs  = query consensus intensities, replicate-weighted

  solve   min_{α}  ½‖Aᵀα − Ī_obs‖² + λ₁‖α‖₁ + λ₂‖α − e_P‖²

  features:
      chi_P     = α_P                              # fitted weight of candidate P
      resid_P   = ‖Ī_obs − α_P·Ā_P‖²              # unexplained residual if P is ALONE
      resid_uni = ‖Ī_obs − Aᵀα*‖²                 # unexplained residual for the best mixture
      P_gain    = resid_P − resid_uni              # how much P ADDS to the best mixture
      contam_frac= 1 − resid_uni / ‖Ī_obs‖²       # global chimeric contamination estimate
```

**`contam_frac` is itself a ranker feature.** A molecule whose spectrum is 40 % chimeric is
one where absence of a peak is weak evidence and presence of a foreign peak is strong
evidence of contamination — and the ranker should trust the forward channels less. This is
free, principled, and completely unused in the field.

**Cost.** With `|union| = 80` and `|peaks| = 128`, `A` is a 128×80 sparse matrix.
Coordinate descent on the LASSO objective, 200 iterations, ~2 ms. Run it only for the top-80
by cheap channels. Trivial.

---

## 4.8 A8 — Channel N: fingerprint posterior (the *only* ML channel)

Deliberately small, deliberately demoted. Its entire role is to be one column among 150.

```
  model    : spectrum → Morgan(2048) bits, trained on COCONUT (natural products, matched
             to test chemistry — docs/03 §4 notes this mattered)
  feature 1: n_bits·( candidate_Morgan  ∧  predicted_bit_distribution )
             (a proper probabilistic dot product: candidate bits vs the model's *predicted
              probabilities*, not a hard 0/1 fingerprint)
  feature 2: log-likelihood of the candidate's bits under the predicted Bernoulli posterior
  feature 3: Jensen–Shannon divergence between the two bit distributions, per candidate
  feature 4: bits predicted with p > 0.7 that the candidate LACKS  → hard-negative evidence
             (this one is often more informative than the positive dot product)
```

Cost of NOT making this bigger: the documented result is that a neural fingerprint model alone
scores 0.145 versus 0.335 for library matching. It is a rounding error next to the
measurement channels. The right response is to give it four good features, not to make it the
system.

---

## 4.9 A9 — Channel H: hard logical masks *(Law 9 — filters, not features)*

These are Boolean. They are applied by *removing candidates before scoring*, and only
candidates surviving all of them enter the ranker. A GBDT is free to ignore a column; it
cannot ignore a removed row.

```
  H1  NITROGEN RULE      odd nominal M ⟺ odd n_N     — exact, from the defect table
  H2  HYDROGEN CAP      n_H ≤ 2n_C + 2 + n_N          — exact
  H3  LEWIS + SENIOR    even-electron; Σvalence ≥ 2(n_atoms−1)   — exact
  H4  HALOGEN           if consensus shows a 3:1 M+2 pair  ⟹ F contains Cl;
                        if 1:1 ⟹ Br; else F contains NO halogen   (docs/02 §6.2)
  H5  SULFUR            M+2 ≈ 4.4 % ⟹ exactly one S (or none, if below detection)
  H6  RDBE CONSISTENCY  the candidate's RDBE must equal the formula's RDBE   — exact
  H7  DIAGNOSTIC IONS   conjunctive substructure requirements, e.g.
                            base peak at 179.03 → candidate MUST contain a
                            C9H7O4 (flavone/flavonol) substructure
                            peak at 91.05    → candidate MUST contain C7H7 (benzyl)
                            peak at 105.03   → candidate MUST contain C7H5O (benzoyl)
                            peak at 121.029  → candidate MUST contain C7H5O2
                            peaks 179 AND 121 AND 151 in a CO chain → flavonoid core
                        (docs/02 Tables 3 and 4 give the dictionary; it is 20 lines of code)
  H8  SUGAR EVIDENCE    a −162.0528 loss ⟹ candidate MUST contain an O-glycosidic link
  H9  ACYL EVIDENCE     a −42.0106 ketene loss ⟹ candidate MUST have an aryl–O–C(=O)–
                        or an aryl ketone; a −30.0106 loss ⟹ an aromatic methoxy or an
                        O in a ring adjacent to carbon
 H10  MASS ACCURACY     |M(P) − M_neutral| must be within rung tolerance AND consistent
                        across all spectra of the molecule (else the structure is misassigned)
```

**H7–H9 are the important ones and they are cheap.** The substructure index is a bitset over
the 712 k-structure pool: `contains[C7H7]` as a 712 k-bit bitmap is **89 kB**. Precompute the
bitmaps once for the ~40 diagnostic substructures (one pass over the pool with RDKit
`HasSubstructMatch`, parallel, ~20 min). Every subsequent mask is a bitwise AND and a
popcount. This turns "does candidate P contain benzyl" from an RDKit call into a bit test.

**Safety.** Hard masks can be wrong. Mitigations, in order:
1. Every mask is computed with a tolerance and an "inconclusive" state. `H4` returns
   `has_Cl | unknown | no_Cl` — never `no_Cl` from a spectrum with too few peaks to resolve
   M+2.
2. Every mask's firing is logged per molecule. On validation, if a mask ever excludes the
   correct answer, the mask's threshold is relaxed, never the mask's absence from the log.
3. Nothing above ~5 % prevalence gets hard status. Below that it stays a feature.

---

## 4.10 A10 — Channel P: the popularity prior

Reproduced because it is worth +0.010 and costs nothing:

```
  final_prior(P) = μ · ( log(1 + n_pubchem_substances(P)) + log(1 + n_pubmed_hits(P)) )
  final_score(P) = z_rank(ranker(P)) + final_prior(P) · 1{ P came from a structure database }
  μ = 0.15 … 0.25, tuned on the public leaderboard (NOT on the 250-structure val set)
```

**Three refinements that make it better than the published version:**

1. **Split the prior by novelty class.** Class 1 and Class 2 molecules are drawn from
   *databases*, so the "is it famous?" prior is informative. Class 3 molecules are not in any
   database, so for a candidate the system believes is Class 3, the prior should be **flat**.
   Emitting `μ_F` conditioned on the model's own class posterior
   (`P(class|Q)` from the calibrated model) is strictly better than a constant μ.
2. **Add natural-product-specific popularity.** `n_chebi`, `n_coconut`, `lifesciencedb hits`
   are closer to the test distribution than PubChem counts, which are dominated by
   synthetic screening compounds (i.e. by `enveda-180`, the wrong chemistry).
3. **Use the prior as a *calibration* target, not an additive bonus.** If the ranker is
   already calibrated, the honest place for a prior is *inside* the posterior:
   `log P(P) ∝ log p_ranker(P) + log p_popularity(P)`. An additive bonus applied after a
   calibrated model double-counts and breaks calibration. This is a real correctness
   improvement over the published formula.

---

## 4.11 A11 — Channel X: transductive propagation across the test set *(new, cheap)*

**The observation.** `docs/dataset_desc` says the test molecules resemble those found in
natural samples. Natural products come in **families** — 8 of the 250 `enveda-np-examples`
molecules are flavonoid glycosides of the same aglycone; RIKEN's plant-specialised-metabolite
library is full of such series. And the *test set is scored as a whole*, so we may look at all
400 molecules at once.

**The algorithm.**

```
 1. Build a spectral similarity graph over the 400 test molecules
      edge weight w_ij = max over spectra pairs of Sim_ent(Q_i, Q_j)
                         × exp(−|ΔM|/τM)          # homologues also share mass series
                         × exp(−|E_eff_i − E_eff_j|/20)

 2. Soft-cluster: connected components at w ≥ 0.45, then Leiden on the weighted graph.
    Typical: 400 molecules → 90–130 components; ~70 % of molecules land in a component of
    size ≥ 3. A size-3 component of C16H12O6 flavonoid glycosides is a *hypothesis*, not a
    coincidence.

 3. Confidence-weighted label propagation (Zhou-style, but on spectra not labels):
      initialise  y_i^(0) = one-hot over the candidate list from stage-1 scoring
      iterate     y_i^(k+1) = α · A_norm y^(k) + (1−α) y_i^(0),      α = 0.35
      normalise rows
    Where α·A y pulls a molecule toward the posterior of its spectral neighbours.

 4. Blend in the *aggregate* evidence of the cluster, with an explicit "family-level"
      hypothesis: the strongest shared substructure across the component's top candidates
      is evidence for all of them (Law 9's conjunction, now at cluster level).

 5. Guard: any molecule whose cluster-level promotion would override a strong library match
      (lib_max ≥ 0.8) is refused. Law 1 again.
```

**Why this is safe and why it is worth it.** It cannot invent a molecule: it only *reorders*
posteriors toward those of spectrally similar test molecules, and the seed posteriors come
from Track A. If the whole cluster is wrong, propagation makes the wrong answer more
confident — the risk is real. Hence: (a) α = 0.35 is conservative, (b) the seed never
changes, (c) the gain is applied *only* where the seed's top-1 and top-5 candidates are
mutually consistent across the cluster, (d) it is a post-hoc reorderable stage so its value
can be measured on the leaderboard with A/B submissions.

---

## 4.12 A12 — The Track A ranker

```
  groups        = molecule_id
  objective     = LightGBM lambdarank, ndcg_eval_at=[1,5,25], label_gain=[0, 1]
                 (labels derived from the calibration set, not from the raw rank — see below)
  rows          ≈ 150–400 candidates per molecule (rung1 ∩ formula ∩ ¬hard masks) ∪ DSM
  features      ≈ 150–220 columns, listed in Appendix A.3
  class balance = the hard-masked rows never enter, so no masking feature is needed —
                  instead a `rung_membership` and `origin` flag is required so the model
                  knows whether a candidate came from PubChem, DSM, or the base pool
```

**The training-label subtlety that most people get wrong.** You cannot train a ranker on
training-set labels directly, because in training data the answer almost always *has* a
reference spectrum, so the model learns "high `lib_max` ⇒ correct" and then fails on Class 2
and 3. The fix, and it is the single most important validation design decision in this plan:

```
  Build the ranker's TRAINING set from enveda-np-examples only, under three masks
  (see Part 7), so the model sees the full range of evidence-availability regimes.
  Sample Class-1 :Class-2 :Class-3 regimes at 45 : 40 : 15 — deliberately oversampling the
  hard regimes relative to the true test mix, because they are where the discriminative
  signal lives and where the model would otherwise be underfit.
```

---
# Part 5 — TRACK B: Unknown Molecule from Data Alone

## 5.0 Scope, independence, and honest expected value

**Independence requirement (yours, honoured literally).** Track B's *core* — the assembly
engine — never reads a structure database. Its inputs are:

```
   ALLOWED INPUTS                      FORBIDDEN INPUTS
   ──────────────                       ─────────────────
   the query molecule's spectra         train.parquet normalized_smiles
   the learned ALPHABET (a mass table)  COCONUT / PubChem / NPAtlas structures
   the chemical element table           any candidate list from Track A
   generic RDKit valence/aromaticity rules
```

The only permitted cross-track contact is at the **merge gate** (§6.1), where Track B's
output structures are scored by Track A's channels on equal footing. Track B never sees
Track A's reasoning, and Track A never sees Track B's reasoning.

**No neural network, no diffusion, no language model.** The entire engine is: integer
programming, dynamic programming, bipartite assignment, and a substructure dictionary. This
is a deliberate design constraint, and it is the constraint that makes it novel — every
published de-novo MS/MS method from DiffMS through MARLIN to MS-Flow is a generative model,
and they all top out in the low single digits to low twenties of exact match.

**Three modes of operation, in decreasing priority by expected value:**

| Mode | What it does | Value | Effort | Target class |
|---|---|---|---|---|
| **B-critic** | For each Track A candidate, asks *"which bonds must break for this molecule to yield exactly the peaks observed?"* and returns the fit quality | **High** — it is a sharper discriminator than any similarity score | Medium | 1, 2, **and 3** |
| **B-tail** | Emits novel structures for slots 15–25 | Low but free (0.04–0.05 × hit rate) | Low | 3 |
| **B-assembly** | Full de novo reconstruction of an unseen skeleton | Low, but it is the scientifically honest core | High | 3 |

**Why B-critic is worth more than B-assembly, quantitatively.** Track A's `docs/03` score
table shows the whole library of scoring tricks plateaus around 0.42–0.47. Every one of those
is a *similarity* between a candidate and the spectrum. The remaining errors are the isomer
ties of `docs/02` §14 — molecules that are genuinely under-determined by *any* similarity.
B-critic attacks a different question: not "how similar is P's predicted spectrum to Q?" but
**"is there a set of bonds whose cleavage reproduces Q's peaks, and is that set consistent
with P's bond types?"** Two isomers that produce similar spectra can still differ sharply in
*which bonds* produce them. This is the only signal class available that the field does not
exploit.

---
## 5.1 B1 — Alphabet inference: discovering the building blocks from the corpus

### 5.1.1 The assumption we are deleting

Every structure-search method in the literature — MetFrag, CFM-ID, MzTree, the spectrum-graph
de novo peptide assemblers, MAGMa — assumes the set of "residue masses" is *known*. For
peptides the alphabet is the 20 amino acids. For small molecules there is no fixed alphabet:
sugars, terpenoids, polyketide units, and halogenated aromatics each contribute their own.

Kreitzberg, Bern, Shu, Yang & Serang (*J Proteome Res* 2019, 18:3268) removed exactly this
assumption for glycans by **MAP alphabet inference**: find the mass set `δ` maximising
`Π_l P(D^(l) | δ) · P(δ)` over a corpus of spectra, by Gibbs sampling one `δ_k` at a time,
accepted by Metropolis–Hastings with a tempering exponent θ. Their result recovered hexose,
HexNAc, dHex and NeuAc **plus three previously unknown masses** to 4–5 digits. That is the
method we port, with three changes for small molecules.

### 5.1.2 Our corpus (this is a design decision, not an accident)

**We run alphabet inference over the 1 213 test spectra themselves** — plus, optionally, the
natural-product subsets of `train.parquet` (`riken`, `gnps`, `mona`, `massbank`,
`spectraverse`, `msdial`) used as *unlabelled* data. Labels are not used; only peak lists.
This is legitimate: it is unsupervised inference on the test distribution, and it means the
alphabet is the alphabet *of the chemistry actually present in this sample*, not the
alphabet of the 2026 literature.

```
  Preconditioning (important, and where most of the gain is)
    • deisotope (isotope satellites are a 1.003355 / 1.99705 redundancy that swamps the
      difference statistics)
    • remove precursor-adjacent peaks (mz > M_neutral − 2)
    • keep only peaks with intensity ≥ 1 %  → median ~12–15 peaks per spectrum
    • discard m/z < 40 (below that, background ions and solvent dominate)
    → ~16 000 peaks over 1 213 spectra. Gibbs sampling at 16 000 iterations is < 1 minute
      on one core.
```

### 5.1.3 Three changes for branching molecules

**(a) Hyperedges, not just edges.** In a peptide, a residue contributes an *edge* between two
peaks. In a branched molecule, a substructure contributes an **n-ary condensation**: `n ≥ 2`
fragments joining with `n−1` bond formations and a possible neutral elimination. So the
alphabet search must consider, for a candidate mass δ, whether it appears as:

```
   2-ary :  m_b − m_a = δ                              (one bond cleavage)
   3-ary :  m_a + m_b − δ = m_c   (i.e. m_a + m_b − m_c = δ)   (two bonds + one loss)
   n-ary :  any subset S ⊆ peaks with  Σ_{i∈S} m_i − m_parent = δ
```

The `n`-ary case is the exponential part. We restrict it: `|S| ≤ 4` and only consider S whose
span is < 250 Da (a substructure plus its complement). That is ~3·10⁴ candidate subsets per
spectrum — fine, and it is exactly the "combining fragments by mass" recursion used in the
assembly-index algorithms of Marshall et al. (*Nat Commun* 2021, 12:3033) and the split-branch
algorithm of the assembly-theory toolchain.

**(b) Intensity-aware membership.** Glycan alphabet projection scores edges by
`Π p_i · p_j`. We additionally require that a δ is **intensity-coherent**: a real
substructure's fragments have intensities that co-vary across spectra in which the substructure
survives. So each δ carries two weights: `freq(δ)` (how often δ appears) and `coh(δ)` (the
correlation of the involved peaks' intensities across the corpus). `coh` separates a real
building block from a coincidental mass difference, which is the dominant false-positive mode.

**(c) Canonical-mass mapping is mandatory.** The Kreitzberg paper's own critical step: test
`δ_k` **and its integer fractions `δ_k/c` and multiples**, and keep whichever generates the
most edges. Without it the sampler reports water as 36.02 instead of 18.01 when only z=1,2
gaps are present. We adopt this verbatim — it is a bug that will otherwise silently destroy
the alphabet.

### 5.1.4 The output object

```
  ALPHABET = [
     { mass: 18.010565,  tol: 0.005, freq: 0.94, coh: 0.71, arity: 2,
       kind: "elimination",   smarts: "[OX2H]-[CX4] >> [C]=[C]      # water loss
       typical_bond: "CO-C" },
     { mass: 162.052824, tol: 0.005, freq: 0.41, coh: 0.83, arity: 2,
       kind: "glycoside",    smarts: "[OX2]-[c] >> [*]O[C@@H]1O[C@H]([C@H]([C@@H](O1)CO)O)O",
       typical_bond: "O-C_ar" },
     ...
  ]
```

The `smarts` field is what makes the alphabet *usable* rather than merely interesting: once a
δ is identified as a mass, we find which chemical operations have that exact mass by
enumerating reaction SMARTS from RDKit's `rdChemReactions` + a hand-mined list of ~60
natural-product biotransformations (methylation, glycosylation, oxidation, reduction,
hydrolysis, conjugation with glycine/sulfate/glucose/methyl). This bridge from *mass* to
*chemistry* is what separates B from the published alphabets, which stop at masses.

---

## 5.2 B2 — The structure-space fragmentation graph *(nodes are substructures, not formulas)*

### 5.2.1 The upgrade

SIRIUS's fragmentation graph has **formula** nodes: for each observed peak, enumerate all
CHONPS formulas matching its mass; connect parent–child when one is a subformula of the other;
then take the maximum-weight colourful subtree. Powerful, but a formula carries no
structural information — every C₆H₅O⁺ is the same node.

We replace the node label with a **substructure**, using a fixed dictionary:

```
  SUBSTRUCTURE DICTIONARY  (built once, from chemistry, not from any test label)

  source : (a) all ring systems ≤ 8 atoms from RDKit's ring database + our own enumerator
           (b) all acyclic functional-group skeletons: ester, amide, ether, ketone,
               aldehyde, acid, amine, sulfate, phosphate, nitro, halide, nitrile, ...
           (c) every ALPHABET entry from B1, expanded into actual substructure graphs
           (d) frequency-filtered: only substructures occurring in ≥ 50 of the pool's
               structures (so the dictionary is the *chemically common* vocabulary)
  size   : ~8 000–15 000 substructures, as canonical SMARTS + InChIKey14 of the fragment

  For each observed peak i:
      N(i) = { substructure S : |m(S) − m_i| ≤ τ_spec }         # mass-consistent
             filtered by |C(S) − Ĉ(i)| ≤ 2                        # carbon-count evidence §5.2.2
             filtered by isotope evidence (§5.2.2)
             typically |N(i)| ∈ [1, 40]; median ~8
```

The carbon-count filter is the one that makes this tractable and it is free:

```
  Ĉ(i) ≈ ( I(M+1 peak of i) / I(i) ) / 0.0107      # 1.07 % ¹³C per carbon
  For a peak at m/z 179 with I(M+1)/I = 0.135  ⟹  Ĉ ≈ 12.6  ⟹ only substructures with
  11–14 carbons survive. Across a 20-peak spectrum this typically reduces the total node set
  from ~10^5 to ~10^3.
```

### 5.2.2 The graph

```
  G_B = ( V, E ) where
     V = ∪_i N(i)   ∪   { ROOT = the query molecule }
     E = { (S, S') : S' ⊂ S and both are mass-consistent with observed peaks,
                    or S' is reachable from S by one ALPHABET-consistent cleavage }
     edge weight w(S→S') = loss_prior(S ⊖ S')  ×  isotope_coherence  ×  arity_fit

  Then: find the maximum-weight SUBTREE in which each observed peak is used at most once
        (Maximum Colourful Subtree, NP-hard in general; solved by ILP, as SIRIUS does, or
         by a branch-and-bound over the ≤ 40-node-per-peak neighbourhoods).
  The winning subtree is the query's *structural fragmentation signature*.
```

**Why this beats a formula tree.** A formula tree says "peak 179 is C9H7O4, derived from the
root by losing C6H3O3." A structure-space tree says "peak 179 is the C9H7O4 fragment
`c1ccc(O)cc1C(=O)OC`, derived by cleaving the ring-B/C-ring bond of a flavonoid." The second
statement is worth 6–8 order-of-magnitude more search pruning, because the second one is
directly checkable against a *substructure dictionary* and directly usable for *assembly*.

---
## 5.3 B3 — Mass-conserving assembly DP *(the core algorithmic contribution)*

### 5.3.1 The formulation

We want a molecular graph `P` with:

```
  (A) exact composition:   atoms(P) = F₀   exactly (from §3.1.3, or over the hypothesis set)
  (B) forward consistency:
        every substructure of P with mass in the observed-peak set should be predictable,
        and the set of "large" substructures of P should be ≈ the observed peaks
  (C) backward consistency:
        every observed strong peak is produced by some cleavage of P,
        and the implied neutral loss is either 0, a mass in ALPHABET, or a documented
        rearrangement
  (D) energy consistency:
        the cleavage events that must occur are energetically accessible at the measured E
  (E) parsimony:   the number of distinct substructures that must be posited is minimal
```

(A) + (B) + (E) is a **weighted set-cover / exact-cover problem**. (A) is an *exact* cover on
atoms. This is where we depart from everything published.

### 5.3.2 The DP

State = a partial assembly = a multiset of already-built fragments whose atom multisets are
pairwise disjoint and whose union is a sub-multiset of `F₀`.

```
  DP over states, ordered by decreasing mass (so a state's predecessors are always lighter):

  state = ( atom_multiset_as_bitmask_or_counter , structure_assignment )
  We CANNOT store structure in the state (too big). Instead:

  ► We store the *mass* state exactly and the *structure* implicitly via a parent pointer,
    exactly like the spectrum-graph DP of Dančík et al. (1999) which works for peptides
    because the mass uniquely identifies the prefix. It does NOT here (isomers), so we
    keep the top-K structures per mass state by partial score.  K = 200.

  transition:  state(mass w, structure π)  ⊕  building-block B
               →  state(w + m(B) − μ(B), π ⊕ B)
      where μ(B) is the eliminated neutral (0 for a condensation, 18.010565 for a
      glycosylation-like condensation-with-loss, 27.994915 for a CO-loss condensation, …)
      SUBJECT TO:  π ⊕ B is valence-feasible (below), and  w + m(B) − μ(B) ∈ observed masses
                   ∪ { intermediate values that admit a completion }

  INIT:  states = { each single building block B : m(B) ∈ observed peaks or < 40 }
  GOAL:  states with mass = M_neutral  and  atoms = F₀ exactly

  SCORE(state) = Σ over the K best transitions of  [ cleavage_likelihood(bond formed) ]
                 − Σ over unexplained observed peaks of  gap_penalty
                 − λ · (number of distinct building blocks used)         ← parsimony (E)

  Termination:  branch-and-bound with an admissible upper bound
                UB(state) = score(state) + Σ_{atoms remaining} max possible gain
  Ranking:      keep the top 500 goal states, canonicalise each with RDKit
                (no stereo, remove salts, standardise), dedup on InChIKey14
```

### 5.3.3 The exact-cover formulation (the rigorous version)

The DP above is a heuristic over structure. The ILP below is exact over composition and is
what guarantees (A) and gives a **provable** no-solution certificate.

```
  BINARY VARIABLES
      y_b ∈ {0,1}              for each building block b in the (pruned) dictionary
      x_ij ∈ {0,1}             for each (block b, peak i) assignment of b to peak i
      z_i ∈ {0,1}              peak i is explained

  MAXIMISE   Σ_b  s(b)·y_b  +  Σ_i  w_i·z_i
             − λ·Σ_b  a(b)·y_b                     (a(b) = atom count of b)

  s.t.  (A1)  Σ_b  y_b  ≥  1                            at least one building block
        (A2)  Σ_b  atoms(b)·y_b  =  F₀                   EXACT cover of atoms
        (B1)  Σ_b  x_{b,i}·y_b  =  z_i  ∀i              peak i explained by ≥1 block
        (C1)  Σ_b  x_{b,i}  ≤  1  ∀i                    each peak used at most once
        (C2)  Σ_i  x_{b,i}·y_b  ≥  1  ∀b, b selected    a selected block must claim a peak
        (D1)  compatibility(i,b)                          from §5.2's node sets
        (E1)  Σ_b  μ(b)  =  Σ_i unaccounted neutral      mass balance
        (F1)  Σ_b  ρ(b)·y_b  ≤  1     ∀ bond class       chemical-valence compatibility of
                                                            simultaneous bonds at one atom
```

`(F1)` is what stops the ILP from producing valence-nonsense: a benzene ring, a furan and a
piperidine cannot all be bonded to the same carbon. It is a clique constraint on the
compatibility graph, i.e. a **maximum-weight clique** subproblem, and modern clique solvers
(MaxSAT / MILP with clique cuts) handle it.

**Solver.** HiGHS or CBC via `pulp`/`mip`, warm-started from the DP's best solution. For a
CHNOS formula, `|blocks| ≈ 3 000` after carbon-count and isotope filtering, `|peaks| ≈ 20`,
so ~63 000 binaries — but >95 % are fixed to 0 by (D1) before the solver starts, leaving
typically < 500 live variables. Solvable in seconds to minutes.

### 5.3.4 Cost control, and what gets cut

The combinatorial explosion is real, so cut in this order, and log each cut's cost in
`recall@assembly`:

```
  cut 1  discard blocks whose mass matches no observed peak AND whose mass > 60 Da
         (large substructures always fragment to something observed)
  cut 2  carbon-count filter (§5.2.2) — kills ~90 % of blocks
  cut 3  isotope filter — kills another ~60 %
  cut 4  require Σ atoms(blocks) ≥ F₀ in total mass (reachability), else prune
  cut 5  cap the number of blocks selected at n_max = ceil(MA_upper_bound) + 2, where
         MA is the molecular assembly index. This is a real, published complexity measure
         computable by DP over assembly states with branch-and-bound
         (Marshall et al., Nat Commun 12:3033; assembly-theory toolchain arXiv:2507.08852),
         and it bounds the depth of the assembly DP honestly.
  cut 6  energy: a bond class that cannot break at the measured E is not a valid cleavage,
         so blocks that would require it are dropped from the *cleavage* role (they may still
         appear as an intact terminal piece).
```

### 5.3.5 What the literature says has not been done here

From the research: *"there is no published generative algorithm that starts from fragment masses
and constructs a structure by constraint satisfaction."* The nearest analogues are
(i) spectrum-graph DP, which is **linear-only** and has never been extended to branching
molecules; (ii) assembly-index algorithms, which compute complexity from MSⁿ but are given
the molecule; (iii) MOLGEN/surge, which enumerate formula → structures with **no spectral
evidence at all** and are honest that they need "user experience"; (iv) FOAM's formula-
constrained genetic algorithm, which is a *stochastic search*, is seeded from formula only,
and admits it yields "structures that are valid by chemical valence rules but highly unlikely
to exist."

Our formulation is the missing joint: **spectral evidence as a hard constraint inside an exact
cover, with branching allowed and minimality as the objective.**

---

## 5.4 B4 — Bond completion and valence satisfaction

Having selected blocks, we must decide how to connect them. This is a small graph problem.

```
  1. CANDIDATE BONDS
     For each pair of selected blocks (a,b) with an open valence, propose a bond.
     Open valence types:
         aryl–O (glycosidic), O–glycosidic, O–methyl, N–methyl, C–C, C–C(aromatic),
         C(acyl)–O(ester), C(acyl)–N(amide), C–C(malonyl, +44), C–C(terpene, +136/+204),
         S–O (sulfate), O–P (phosphate), C–halide
     Each bond type carries:
         (i)  the mass it adds: 0, or a small eliminated neutral (H₂O, NH₃, H₂S, HBr, HCl)
         (ii) a valence-compatibility predicate
         (iii) an energy prior p(cleaves at E)

  2. FORMULATION: maximum-weight spanning-forest on the block graph
        maximise  Σ_edges  p(bond)  −  λ_ring · (#rings formed)  −  strain_penalty
        subject to every open valence is either filled by a bond or satisfied by H
     Rings: a spanning forest forbids cycles. Rings are added afterwards as explicit
     "ring-closure bonds" with a separate prior, and the count is penalised by λ_ring.
     Aromaticity is perceived by RDKit after the graph is built.

  3. STEREOCHEMISTRY IS NOT INFERRED AND NOT OUTPUT. InChIKey14 ignores it (Law 4). This is
     not a simplification we chose; it is a decision the scorer made for us, and it removes
     the single largest error source in structure elucidation.
```

---

## 5.5 B5 — Scoring Track B's structures: the same two-sided law

A Track B structure has no reference spectrum, so it is scored exactly like a Track A candidate
with zero measured evidence — which is precisely the Class-2/3 regime. The scorer is Channel F
(§4.5), unmodified:

```
  F_β via Hungarian assignment between the structure's predicted fragments
       (tier-1 fragmenter, energy-gated) and the query consensus peaks, with
       per-peak weight w_i = replicate_fraction_i^0.5 · I_i

  PLUS the coverage-style assembly agreement:
       assembled_blocks_matched  = fraction of the assembly blocks whose mass matches an
                                   observed peak with w_i above the 20th percentile
       unexplained_strong_peaks  = count of peaks with I > 0.05 not produced by ANY
                                   cleavage of the structure        ← the two-sided charge
       chain_consistency         = fraction of assembly-tree edges whose parent and child
                                   BOTH appear in the observed peak set
```

`unexplained_strong_peaks` is the discriminator that makes Track B's outputs compete fairly.
A structure that explains the base peak and two others but leaves nine unexplained strong
peaks is wrong, and a one-sided score would not know.

---

## 5.6 B6 — Model selection: which of 500 assemblies do we believe?

Five hundred exact-cover solutions is not a ranking problem, it is a **model selection**
problem: which assembly is the simplest hypothesis that explains the data?

```
  score(P) = log P(observed peaks | P)        = the two-sided assignment log-likelihood
           − ½ · k · log n                   ← BIC penalty on the free parameters
           where k = #distinct blocks + #ring closures (the free parameters of P)
             n = number of distinct observed peaks explained
           −  MA(P) · γ_MA                   ← assembly-index prior: prefer molecules that
                                                are not gratuitously complex. MA is the
                                                number of joining operations, so this is a
                                                direct Occam term in the language of
                                                assembly theory.
           −  strain(P) · γ_strain            ← RDKit ring-strain + exotic-substructure penalty
           −  rarity(P) · γ_rarity            ← RDKit SA score (Ertl's synthetic-accessibility
                                                proxy); this is NOT a popularity prior, it is
                                                a "does this substructure exist in nature"
                                                prior, and it is what stops the assembly from
                                                proposing furans on strained bicyclics
```

Then select the **MMD / MAP set** (minimum description length) of structures within
`ΔBIC ≤ 2` of the best, which is the standard device for reporting a *set* of acceptable
explanations rather than one. That set is exactly what the MRR@25 metric wants.

---

## 5.7 B7 — B-critic: the bond-ablation auditor

This is the highest-value mode and it is short to specify.

```
  INPUT : a candidate structure P from Track A's pool, and the query molecule Q

  QUESTION: is there a cleavage set B ⊆ bonds(P) such that cleaving P along B produces
            (a subset of) the observed peaks, with the implied neutral losses all
            belonging to ALPHABET ∪ documented rearrangements?

  ALGORITHM (per candidate, ~2–5 ms)
    1. Build the fragmentation graph of P: bridges ∪ 2-cuts (Hopcroft–Tarjan).
       For P with 25 heavy atoms this is ~40 fragments at depth 1, ~400 at depth 2.
    2. Maximum-weight bipartite assignment between fragments and observed peaks, weight
         a_ij = (Ĩ_j/max Ĩ) · w_i^obs  × p(cleaves | bond_type, E)
       — note p(cleaves | bond_type) is the MolDiscovery-style conditional table,
         estimated by counting (bond_type, fragment intensity rank) pairs in training.
    3. OUTPUT:
         ablation_fit   = fraction of observed top-8 peak intensity matched by P
         ablation_pen   = Σ over matched pairs of −log p(cleaves | bond_type, E)
         ablation_necessity = for the top-3 matched peaks, is the matching UNIQUE?
                            (a peak matched by exactly one cleavage path is highly
                             informative; a peak matched by nine is not)
         ablation_novel_bond = the set of bonds P contains that NO reference structure
                            in the pool shares in that arrangement  ← Class-2 signal
    4. FEATURES fed to the ranker: all four, plus their cohort-normalised versions.
```

**Why this is not redundant with Channel F.** Channel F asks "how many of P's fragments match
the query". B-critic asks "which of P's *bonds* would you have to cut, and is that cut
plausible at the measured energy, and are the cuts unique". The second question separates
isomers that the first cannot. Concretely: methyl benzoate and benzyl formate share nearly
the same fragment *masses*; but for methyl benzoate the 91 ion requires cutting a
C(aryl)–CH₂ bond **that does not exist in methyl benzoate**, so it is not merely low
intensity — it is structurally impossible, and B-critic assigns it an explicit hard
infeasibility rather than a small negative contribution.

---

## 5.8 What Track B honestly cannot do

Stated plainly, because a plan that pretends otherwise is not a plan.

```
  ✗  Regiochemistry and positional isomerism inside an aromatic ring. The observed peak
     set usually cannot distinguish 3-OH-quercetin from 3'-OH-quercetin. B will generate
     both and rank them by Occam, which is not evidence.

  ✗  Stereochemistry. Deliberately not attempted (Law 4). Fortunately not scored.

  ✗  Full de-novo skeletons above ~500 Da. Formula enumeration for CHNOS at 1 000 Da admits
     > 350 million formulas; with our spectral constraints the surviving block-set count is
     still typically > 10^6. We cap at 500 Da and hand everything above to Track A.

  ✗  Rearrangement-dominated spectra. ~20–30 % of real peaks do not come from a single bond
     cleavage (docs/02 §7.5). We handle the top documented rearrangements (McLafferty +44,
     retro-Diels–Alder, H₂O/CO/NH₃ losses) as explicit graph moves, but an exotic
     rearrangement is invisible to us and is charged as an unexplained peak.

  ✗  Confidence calibration on genuinely novel molecules. Our posterior over structures is
     calibrated on validation data drawn from known molecules. On a truly novel skeleton the
     posterior is over-confident by an unknown factor. Mitigation: cap Track B's rank-1
     probability at 0.35 unless ablation_fit ≥ 0.85.
```

---
# Part 6 — Fusion, Calibration and the Decision-Theoretic Submission

## 6.1 The merge gate

```
  POOL_LIBRARY   =  Track A's scored candidates (rung ∩ formula ∩ ¬hard_masks ∪ DSM)
  POOL_DENOVO    =  Track B's MMD set (≤ 25 structures) ∪ Track B's tail proposals (≤ 40)

  UNION          =  POOL_LIBRARY ∪ POOL_DENOVO,  deduped on InChIKey14
                    (dedup key computed AFTER RDKit tautomer canonicalisation, matching
                     the scorer exactly — otherwise we may "dedup" two structures the
                     scorer would count as one, or fail to dedup two it counts as one)

  RESCORE        =  every candidate in UNION is scored by EVERY channel:
                      • measured channels (L, E) are ZERO for de-novo candidates, and that
                        zero is itself informative — emit `origin ∈ {lib, dsm, denovo}`
                        as a feature so the model can learn what to do with an
                        evidence-free candidate
                      • model channels (F, T, N, S, B-critic) work identically for all

  LIBRARY GATE   =  g = 1 − σ((sim_lib − Λ)/κ),   Λ = 0.70, κ = 0.05
                    every model-channel score is multiplied by g
                    ⟹  a de-novo candidate with excellent forward evidence CAN reach rank 1
                        (Law 5), but only where measurement evidence is silent (Law 1)
```

**The gate is soft, not a switch.** A hard threshold creates a discontinuity that the
validation set cannot calibrate. The sigmoid is monotone, has two interpretable parameters,
and produces a smooth feature.

## 6.2 Fusion: one ranker, not two

Everything is scored by **one** LightGBM lambdarank model over one union. Two separate rankers
were tried by the field (`W##` blend weights in the notebook titles) and beat one general
ranker — but they were splitting on *library provenance*, which is a proxy for
*evidence availability*, and that proxy is now an explicit feature. Splitting on a proxy you
have already measured is strictly worse than feeding the measurement.

## 6.3 Why the submission is literally "sort by posterior" *(the metric derivation)*

Let `p_1 ≥ p_2 ≥ … ≥ p_N` be candidate posteriors and let the true structure occupy unknown
rank `R*`. The score for the molecule is `1/R*` if `R* ≤ 25`, else 0.

```
  E[score | p_1 ≥ … ≥ p_N]  =  Σ_{r=1..25}  P(R* = r) · 1/r
                           =  Σ_{r=1..25}  p_r · 1/r

  (because P(R* = r) = p_r under a perfectly calibrated posterior)

  d/dr [ Σ_r p_r/r ]  w.r.t. swapping adjacent r and r+1:
      = p_{r+1}/(r+1) + p_r/r  −  p_r/(r+1) − p_{r+1}/r
      = (p_r − p_{r+1}) · (1/r − 1/(r+1))  ≥ 0  for p_r ≥ p_{r+1}
  ⟹ sorting by descending p is optimal.  Strictly. Always.
```

**Three consequences, each of which contradicts common intuition:**

1. **Submit exactly 25 whenever you have ≥ 25 distinct candidates.** The marginal
   contribution of the 25th is `p_25/25 ≥ 0` — strictly positive unless `p_25 = 0`. And
   since slot 25 costs a rank only if the answer is *also* earlier, and the answer is only
   once, filling is free.
2. **Never re-order for "diversity".** A structurally diverse tail has *lower* `Σ p_r/r` than
   a probability-sorted one. The intuition that a structurally dissimilar 25th guess is
   "safer" is wrong under this metric: it is only safer if you *believe* it has high
   probability, which is what sorting by p already does.
3. **The tail is where calibration shows up.** Positions 5–25 contribute `Σ_{r=5..25} p_r/r`,
   which for a typical posterior is 0.05–0.15 — a third of the total score. A ranker that
   nails rank 1 but is uncalibrated in the tail loses that third.

## 6.4 Calibration pipeline

```
 STAGE 1  Score
         s(P) = LightGBM lambdarank output for candidate P in molecule Q

 STAGE 2  Normalise within the molecule
         s̃(P) = ( s(P) − median_Q s ) / ( MAD_Q s + 1e-9 )
         (raw GBDT scores are not comparable across molecules of different candidate counts)

 STAGE 3  Softmax with learned temperature
         p_uncal(P) = exp(s̃(P)/T) / Σ_{P'} exp(s̃(P')/T),   T fitted on validation

 STAGE 4  Target–decoy empirical p-value
         Build the decoy set: for each real candidate P, sample 20 decoys that share P's
         formula and approximate mass but are NOT P. Run them through the identical
         pipeline. This yields, for each real candidate, the null distribution of
         "score of a wrong molecule of the right kind".
         p_empirical(P) = (1 + #{decoy scores ≥ s(P)}) / (1 + n_decoy)
         (This is the same target–decoy device COSMIC uses to convert CSI:FingerID scores
          into calibrated p-values, which is precisely what let it beat library search 27:1.)

 STAGE 5  Isotonic regression
         Fit isotonic p(P) → 1{correct} on validation.
         Output the CALIBRATED posterior.
         Report expected calibration error (ECE) and a reliability diagram as a first-class
         validation metric, not a diagnostic you look at once.

 STAGE 6  Merge the popularity prior INSIDE the posterior (Law 8 / refinement 5.10.3)
         log p_final(P)  =  log p(P) + β·(log(1+substances) + log(1+pubmed))
                          + β_np·log(1+coconut_hits) + β_lib·1{has reference spectrum}
         renormalise
         (β fitted on the leaderboard; β_lib captures "Class 1 is likelier than Class 3",
          which is a prior over the label-generating process, not over chemistry)
```

## 6.5 The pool-recall ceiling and how to spend effort

```
  CEILING  =  P(true structure ∈ POOL_LIBRARY ∪ POOL_DENOVO)
            =  recall@pool of the final union

  MEASURED at three depths:   @50, @200, @800
  REPORTED per novelty class

  DECISION RULE
    if recall@50  is high but MRR@1 is low  →  the problem is RANKING. Spend on channels.
    if recall@50  is low                      →  the problem is POOL. Spend on DSM,
                                                 COCONUT/NPAtlas coverage, and Track B.
    if both are high but MRR@25 is low        →  the problem is the EVIDENCE CEILING
                                                 (Law 10.2). Spend on nothing; it is
                                                 physics. Accept it.
```

## 6.6 Submission writer

```python
def write_submission(predictions: dict[str, list[str]], path: str) -> None:
    rows = []
    for mid in EXPECTED_IDS:                       # exactly the test molecule_ids, sorted
        cands = dedup_by_inchikey14(predictions[mid])   # see below
        cands = cands[:25]
        if not cands:                               # never emit an empty cell — it is a
            cands = [global_fallback_for(mid)]      # hard rejection, so guarantee non-empty
        smiles = ";".join(cands)
        assert 1 <= len(cands) <= 25
        assert ";" not in any(cands) and "\n" not in smiles
        rows.append((mid, smiles))
    df = pd.DataFrame(rows, columns=["molecule_id", "smiles"])
    assert df.molecule_id.is_unique and df.notna().all().all()
    df.to_csv(path, index=False)
```

**Dedup must replicate the scorer exactly:**

```python
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

TAUT = rdMolStandardize.TautomerEnumerator()

def scorer_key(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return None
    mol = Chem.Mol(mol)                       # avoid mutating a cached mol
    try:
        mol = rdMolStandardize.Cleanup(mol)
        mol = TAUT.Canonicalize(mol)
    except Exception:
        pass
    return Chem.MolToInchiKey(mol).split("-")[0]     # InChIKey14
```

Three failure modes to guard explicitly, each of which silently costs slots:
1. `MolFromSmiles` returns `None` for one of the 25 → drop it, do not emit it.
2. Tautomer canonicalisation can be slow on large molecules → cache by canonical SMILES, and
   cap per-candidate time with a try/except (the scorer pins RDKit 2026.03.3; we should pin the
   same version in our environment).
3. Two different SMILES with the same InChIKey14 → keep the higher-ranked one, and keep the
   *stereo-free* spelling (smaller string, faster to write, same score).

---
# Part 7 — Validation — a 100× better harness

## 7.1 The trap, restated

`docs/03` §7 is unambiguous: a random holdout from `train.parquet` **ranked three rankers in
the opposite order** from the public leaderboard. Two compounding mismatches:

```
  TRAIN (random holdout)                  TEST
  ────────────────────                    ─────
  82 % non-timsTOF (Orbitrap/QTOF)        100 % timsTOF
  46 % enveda-180 = drug-like synthetic    natural products
  11 × more drug-like than natural         ~0 drug-like
```

A random split validates *drug fragmentation on the wrong instrument*. There is no version of
that which correlates with the board.

## 7.2 The harness has three slices, not one

### Slice A — **instrument-matched, tiny, uncontaminated** (the `enveda-np-examples` 250)

```
  val_A = train[ingest_lib == "enveda-np-examples"]        # 250 structures, 1 151 spectra
  PRIMARY USE: catch large regressions. Do NOT tune μ or temperature here — 250 structures
               cannot resolve ±0.02, and tuning on it is how you overfit 250 molecules.
  MASK 1  everything visible                         ≈ Class 1
  MASK 2  hide that structure's own spectra          ≈ Class 2
  MASK 3  remove that structure from the pool        ≈ Class 3
```

### Slice B — **chemistry-matched, large, instrument-mismatched** *(new; the workhorse)*

This is the slice that does not currently exist anywhere and it is where the discriminative
power lives. Build it from the natural-product libraries by **synthesising the Class-2
condition at scale**:

```
  val_B sources, ranked by on-topic-ness:
      riken        15 892 structures   plant specialised metabolites   ← best chemistry match
      gnps         45 750 structures   largest NP collection, heterogeneous
      mona         11 681 structures   Fiehn Lab
      spectraverse   9 631 structures   harmonised, obscure sources
      massbank      9 180 structures   curated
      msdial        9 127 structures
      ── excluded ──
      enveda-180   182 941 structures  drug-like, WRONG chemistry, dominates by count
      pluskal_ms2   46 821 structures  screening library, wrong chemistry
      drug_plus      2 539             drugs, wrong chemistry

  Construction:
    for each structure P in the selected libraries with ≥ 1 spectrum:
       • all of P's spectra are REMOVED from the search library
         ⟹ Channel L is dead for P, and P's *training-time* reference is gone  → CLASS 2
       • P is KEPT in the candidate pool (we still know its structure)
         ⟹ the pool-recall question is answered by us, not by luck
    This yields ≈ 100 000 Class-2 validation molecules instead of 250.

  Additional sub-slices inside B, which isolate specific failure modes:
    B1  LOW-ENERGY      spectra with CE_ev ≤ 25 only
    B2  HIGH-ENERGY     spectra with CE_ev ≥ 45 only
    B3  NEGATIVE-ION    negative-mode only
    B4  LOW-QUALITY     base_peak_intensity in the lowest tercile
    B5  SINGLE-SPECTRUM molecules with exactly 1 spectrum   ← mirrors the real n=1 case
    B6  SKELETON-DISJOINT  sample so that no two molecules share a Bemis–Murcko scaffold
                         (prevents scaffold leakage inflating every score)
```

**The instrument mismatch in Slice B is a feature, not a flaw.** It is precisely the
condition Channel E (§4.3) is designed for. Running the ablation *with* Slice B and seeing
Channel E's Δ is how we learn whether the η correction is real or a fiction.

### Slice C — **held-out chemistry region** *(scaffold-disjoint, for honest generalisation)*

```
  Cluster all structures by Bemis–Murcko scaffold (MurckoScaffold.GetScaffoldForMol).
  Group into clusters by scaffold + first-shell atom environment (radius-1 Morgan fingerprint
  on the scaffold) so that positional isomers land in the same cluster.
  Hold out entire clusters — never molecules.

  This is the slice that catches "the model has memorised a scaffold family". Published work
  found ~95 % of compounds still elude exact recovery even on leakage-controlled benchmarks,
  and DiffMS drops from 12.4 % (oracle scaffold) to 1.54 % (predicted scaffold) — an 8×
  collapse. Any feature that relies on scaffold memorisation shows it here and nowhere else.
```

## 7.3 What gets measured, per slice, per mask, per channel

```
  POOL METRICS
      pool_recall@50 / @200 / @800       ← the ceiling (Law 10.1)
      pool_size after each filter         ← which filter costs what
      formula_recall@1 / @5              ← §3.1.3 quality
      false_kill_rate per hard mask       ← did any mask ever remove a correct answer?

  RANK METRICS
      MRR@25, MRR@5, MRR@1, hit@1, hit@5, hit@25
      PER CHANNEL, in isolation           ← this is the ablation table that drives all
                                            subsequent engineering decisions
      RANK CORRELATION between channels    ← Kendall τ; low correlation between two channels
                                              means the pair is worth keeping (Law: combining
                                              weak independent signals beats one strong)
      ABLATION GAIN of each channel when removed from the fused ranker, on each slice

  CALIBRATION
      ECE (expected calibration error), reliability diagram
      mean |predicted p₁ − empirical accuracy@1|

  EVIDENCE CEILING (Law 10.2)
      For each validation molecule, the maximum achievable MRR if we had a PERFECT ranker over
      the pool = 1/rank of the answer under an oracle that scores all isomers equally.
      The gap between our MRR and this value is irreducible *given the spectra*.
      Measuring it tells us how much of the remaining loss is fixable at all.
```

## 7.4 Guarding against the new failure mode this harness creates

Slice B is large, so it can overfit in a way 250 cannot. Controls:

```
  1. GROUP-LEVEL cross-validation: split by molecule_id, and additionally verify that the
     scaffold distribution of the fold-0 validation molecules is within a stated overlap
     (report Jaccard of Murcko scaffolds between train folds and val fold).
  2. The ranker's own hyper-parameters are chosen on Slice A + Slice B masks 1–2 only.
     The *final* read-out is on Slice C and on the public leaderboard, both untouched until
     the last week.
  3. Every reported number carries its sample size and a bootstrap CI. On 250 molecules, the
     95 % CI on MRR@25 is roughly ±0.06 — that is the resolution limit and pretending
     otherwise is how people chase noise for a month.
  4. Explicitly forbidden: tuning `collision_energy_ev` conversion, the hard-mask thresholds,
     or the peak-intensity floor on the public leaderboard more than three times. Each such
     tune is one bit of leaderboard overfitting; three is a strategy.
```

---
# Part 8 — Build Order and Sprint Plan

Total: **22 weeks** to the 14 Dec 2026 final submission (start 14 Sep 2026), with
checkpointed deliverables. Every stage is independently submittable. **A valid `submission.csv`
exists from day 3 and is re-submitted after every stage**, so there is never a zero.

---

## Sprint 0 — Foundation (week 1)

| # | Task | Output | Gate |
|---|---|---|---|
| 0.1 | `data_layer/ingest.py` streaming reader, binary peak blobs | `peaks.lms`, `meta.parquet` | peak count parity with pyarrow on a 1 % sample |
| 0.2 | Adduct table + neutral-mass function with the three assertions in §3.1.1 | unit tests | tests pass |
| 0.3 | Structure index: dedup on InChIKey14, pack Morgan-2048 to 256 B | `fp2048.npy` (182 MB), `structs.parquet` | ≤ 950 MB RSS |
| 0.4 | Sorted mass array + all rung windows (§4.1.1) | `massidx.npy` | pool_recall@800 on Slice A computed |
| 0.5 | **Baseline submission**: rung-1 pool, popularity order, 25 slots | `submission.csv` | format validator passes |

> **Checkpoint 0.5 — SUBMIT.** Expected MRR@25 ≈ 0.05–0.10.

## Sprint 1 — Formula + curation + consensus (weeks 2–3)

| # | Task |
|---|---|
| 1.1 | Formula ILP (§3.1.3) + Seven-Golden-Rules filter + isotope scoring |
| 1.2 | Fragmentation-tree MAP (SIRIUS-style colourful subtree) for `formula_hypotheses` |
| 1.3 | Formula Subset Analysis / MaxComb=1 score |
| 1.4 | Peak curation suite (§3.2) with ablation switches |
| 1.5 | **Replicate consensus** (§3.3) — one of the highest value-per-line items in the project |
| 1.6 | Slice B construction script (the 100 k Class-2 validation set) |

> **Checkpoint 1.6.** Measure: does consensus improve MRR on Slice B5 (single-spectrum) vs
> B-multi? If consensus *hurts* on n=1 molecules (it must not — it falls back to the prior),
> there is a bug.

## Sprint 2 — Measured channels: L and E (weeks 4–6)

| # | Task |
|---|---|
| 2.1 | Spectral entropy with dynamic weighting + coverage + cohort normalisation |
| 2.2 | Peak inverted index over 2.5 M spectra, memory-mapped, 8 workers |
| 2.3 | Modified-cosine offset search with explicit Δ tracking (§4.4.2 Step 1) |
| 2.4 | **Channel E: η estimation + energy-rescaled search** (§4.3) |
| 2.5 | Ablation: `lib_max` alone vs `lib_score` alone; Slice A masks 1/2/3 |

> **Checkpoint 2.5.** Channel E's Δ on Slice B is the go/no-go for the whole §4.3 idea.
> If Δ(MRR) < +0.005 across all of Slice B, demote E to a feature and stop investing.

## Sprint 3 — Model channels: F, T, N (weeks 7–9)

| # | Task |
|---|---|
| 3.1 | MetFrag-class fragmenter with bridges/2-cuts/Hopcroft–Tarjan + rearrangements |
| 3.2 | **Two-sided Hungarian assignment + F_β + cohort LLR** (§4.5) |
| 3.3 | **Ablation F vs the published one-sided `frag_score`** — this is the headline experiment |
| 3.4 | Structure-space fragmentation tree + tree alignment (Channel T) |
| 3.5 | Fingerprint model on COCONUT → 4 features (Channel N) |
| 3.6 | Ranker v1: LightGBM over these channels; Slice A/B validation |

> **Checkpoint 3.3.** If two-sided F does not beat one-sided `frag_score` by ≥ +0.01 MRR on
> Slice B, the theory in §4.5 is wrong and the channel reverts to one-sided with a
> precision penalty. Run this experiment early; it is cheap and it de-risks a core claim.

## Sprint 4 — Constructive + critic channels: Δ and B-critic (weeks 10–12)

| # | Task |
|---|---|
| 4.1 | Edit dictionary mining from the training corpus (§4.4.2 Step 3) |
| 4.2 | DSM: apply edits via reaction SMARTS, generate candidates, integrate into the pool |
| 4.3 | B-critic: bond-ablation auditor (§5.7) with the ablation_* feature family |
| 4.4 | Hard-mask engine (Channel H) + the diagnostic substructure bitmaps |
| 4.5 | Popularity prior with the class-conditional μ |

> **Checkpoint 4.2.** DSM's value is entirely in Slice B (Class 2) and is measurable:
> `pool_recall@800` on Slice B before and after DSM. If it moves by < 2 points, DSM is dead
> weight and gets demoted to a tail-filler.

## Sprint 5 — Track B core (weeks 13–17)

| # | Task |
|---|---|
| 5.1 | Substructure dictionary (8–15 k entries) + canonical SMARTS |
| 5.2 | B1: alphabet inference by Gibbs over the test corpus; canonical-mass mapping |
| 5.3 | B2: structure-space fragmentation graph with carbon-count and isotope filters |
| 5.4 | B3a: mass-conserving assembly DP with K-best per mass state |
| 5.5 | B3b: exact-cover ILP with the clique (valence) constraint |
| 5.6 | B4: bond completion via maximum-weight spanning forest + ring closures |
| 5.7 | B5/B6: two-sided scoring + BIC/MA model selection → MMD set |
| 5.8 | B7 ablation mode is already shipped from Sprint 4; here we add Track B's *own* candidates |
| 5.9 | Channel S: sparse mixture decoding (CHIMERYS-style L¹) |

> **Checkpoint 5.5.** Report `assembly_recall@25` on Slice A mask 3 (the 250 structures with
> no pool entry). This is the honest number for Track B. Published de-novo exact-match is
> 4–20 %; anything ≥ 10 % is a good result, and we ship it regardless of its small EV,
> because it is the scientifically interesting part and it costs the leaderboard nothing
> (Law 1's gate).

## Sprint 6 — Fusion, calibration, transduction (weeks 18–19)

| # | Task |
|---|---|
| 6.1 | Merge gate with the soft library gate |
| 6.2 | Target–decoy + isotonic calibration, ECE as a tracked metric |
| 6.3 | Posterior-merged popularity prior |
| 6.4 | Channel X: transductive test-set propagation |
| 6.5 | Decision-theoretic tail allocation + writer + validator |
| 6.6 | Slice C (scaffold-disjoint) audit: which features are memorising scaffolds? |

## Sprint 7 — Packaging and the Kaggle notebook (weeks 20–22)

| # | Task |
|---|---|
| 7.1 | One `submission.csv` generator script, offline, ≤ 9 h, ≤ 20 GB RAM |
| 7.2 | All artefacts pre-uploaded as Kaggle datasets (internet is OFF) |
| 7.3 | Determinism: fixed seeds, sorted iteration, pinned RDKit `2026.03.3` |
| 7.4 | Dry-run on the sample test set, then 3 leaderboard probes max (§7.4 rule 4) |
| 7.5 | Final: ablate each channel against the last two leaderboard scores and keep what earns |

---

## 8.1 The five experiments to run in the first two weeks, before writing any ranker

These are cheap, and each one can kill a major component. Running them late is expensive.

| # | Experiment | Kills what if it fails |
|---|---|---|
| E1 | Formula ILP top-1 accuracy on Slice A + first 5 000 of Slice B | the formula stage; replace with SIRIUS's decomp + BLAST |
| E2 | Consensus peaks (§3.3) vs raw peaks, on Slice B multi-spectrum molecules | §3.3 |
| E3 | Channel E (energy-rescaled) vs plain offset search, on Slice B | §4.3 entirely |
| E4 | Two-sided F_β vs one-sided `frag_score` on Slice A mask 2 + Slice B | §4.5 |
| E5 | DSM pool recall gain on Slice B | §4.4 |

---
# Part 9 — Compute, Memory and Runtime Budget

## 9.1 Memory ledger (the hard constraint: ≤ 1 GB resident for research/iteration; the Kaggle
run may use more but the design stays within 1 GB so that local iteration is free)

```
 STAGE                                   RESIDENT     ON DISK      NOTES
 ─────────────────────────────────────   ─────────    ─────────    ──────
 Train ingest (streaming, 20 k rows)     30 MB        440 MB      never held whole
 Formula index (712 k × 20 B)           14 MB         14 MB      mass f64 + f_i14 i32
 Structure fingerprint table 2048-bit    182 MB       182 MB      packed uint8[256]
   └ fallback 1024-bit                    91 MB        91 MB      if the box is tight
 Structural hard-mask bitmaps (~40)      3.6 MB       3.6 MB      89 kB each
 Peak blob (per library)                   0          440 MB      mmap, random access
 Peak inverted index                      0          900 MB      mmap; postings as int32
 Fragment cache (top-80 candidates)      18 MB         18 MB      rebuilt per molecule
 Ranker matrix (400 rows × 200 feats)    0.3 MB             —      float32
 Substructure dictionary                 11 MB         11 MB      SMARTS + masses
 ─────────────────────────────────────   ─────────    ─────────
 PEAK RESIDENT                                       ≈ 260 MB
 HEADROOM to 1 GB                                       ≈ 740 MB
```

The design consequence: **every index is memory-mapped and read on demand.** No Python dict
larger than 10 k entries, ever. No `pandas.read_parquet` on peak columns, ever.

## 9.2 Runtime ledger (per full 400-molecule test run; 8 workers)

```
 ─────────────────────────────────────────────────────   ────────   ────────
 STAGE                                                   1 CORE     8 CORES
 ─────────────────────────────────────────────────────   ────────   ────────
 ingest train.parquet → binary blobs                     28 min     4 min
 structure dedup + InChIKey14 + fingerprints             22 min     3 min
 substructure bitmaps + diagnostic masks                  18 min     3 min
 formula ILP + isotope + FSA  (400 molecules)            40 s       8 s
 replicate consensus                                      6 s       2 s
 Channel L: entropy vs ~2.5 M spectra                    46 min     6 min
 Channel E: η fit + energy-rescaled offset search         52 min     7 min
 Channel Δ: DSM (top-30 refs × 3 δ × 8 sites)             14 min     3 min
 Channel F: fragmenter + Hungarian, top-120 candidates   38 min     6 min
 Channel T: FT build + alignment, top-120                 22 min     4 min
 Channel N: fingerprint model (GPU, offline-trained)      —         —
 Channel S: sparse mixture decode, top-80                18 s       4 s
 Channel H: hard masks                                     2 s      < 1 s
 Channel B-critic: ablation audit, top-120               26 min     4 min
 Track B: alphabet + graph + DP + ILP, top-40 molecules  71 min    12 min
 ranker + calibration + writer                            22 s      4 s
 ─────────────────────────────────────────────────────   ─────────  ────────
 TOTAL (CPU-only)                                                    ≈ 54 min
 TOTAL (incl. frozen model inference)                                  ≈ 60 min
```

Comfortably inside the 9-hour Kaggle limit. The margin is deliberate: we want room for the
full Track B on *every* molecule, not just the top 40, if the measurements say it pays.

## 9.3 What we give up if the box is smaller than 1 GB

Ordered by loss, so that if you are squeezed you drop from the bottom:

```
 DROP 1  fragment cache on disk → recompute.                    −0 runtime, 0 score
 DROP 2  fingerprint table 2048 → 1024 bits.                    −0 score (measured), −91 MB
 DROP 3  Channel T alignment → top-40 candidates only.         ~−0.003 MRR
 DROP 4  Channel E to timsTOF + np-examples only.               ~−0.02 MRR  (last resort)
 DROP 5  Track B assembly to top-20 molecules.                  ~−0.001 MRR
 NEVER DROP: formula, mass filter, Channel L, two-sided F, hard masks, dedup, the 25 slots.
```

---
# Part 10 — Novelty Ledger

Every claim of novelty in this plan, stated honestly: what is new, what is *not* new, and what
it is built on. Nothing here is claimed to be unprecedented without checking.

| # | Component | Novelty status | Built on | Prior closest work |
|---|---|---|---|---|
| N1 | **Two-sided Hungarian assignment scoring with F_β + cohort LLR** | **New for this problem.** The two-sidedness is textbook (optimal transport / LSAP) and has been evaluated in MS/MS (NIST's ξ₆ Earth-mover; WABI-2018; CHICKN) but **never as the primary candidate-scoring function in a metabolite ID pipeline**. | LSAP; MetFrag's fragmenter; MolDiscovery's LLR | `docs/03` §4 uses one-sided `frag_score` |
| N2 | **Energy-harmonised offset search** (learned per-library η) | **New.** Per-library CE response curves are published (Rondeau's collisional efficiencies; Palit & Mallard's FEI); applying a *fitted* η to make cross-instrument spectra comparable **for a ranking pipeline, with the residual gap exposed as a ranker feature** is not. | Rondeau 2020; Palit & Mallard 2009; SIRIUS linear-time modified cosine | Field uses *restrictions* (AFIX/BIO blend), not conversion |
| N3 | **Delta-Structure Mining: retrieved spectrum → mass shift → chemical edit → new structure** | **New and central.** Prior art stops at *ranking* retrieved cognates (NIST hybrid; MZmine analog search; ModiFinder localises the modification). **Nobody edits the retrieved structure into existence as a candidate.** The mass→reaction-SMARTS dictionary mined from the corpus is also new. | RDKit reactions; NIST hybrid Δm; modification-mass tables | the field's analog channel ends at Tanimoto |
| N4 | **Replicate-consensus peak voting at the molecule level** | **Adapted, not new** — consensus spectra are decades-old in proteomics, and chimeric cleaning from replicate statistics is DNMS2Purifier. **New as the primary aggregation operator for a multi-spectrum prediction task**, and it replaces the field's MIN/mean heuristic. | Beavis & Craig; DNMS2Purifier (*Anal Chem* 2023) | `docs/03` §12.7 "multi-spectrum aggregation: MIN" |
| N5 | **Sparse mixture decoding (CHIMERYS-style L¹) for chimeric natural-product spectra** | **Port, new application.** The mathematics is published for peptides; applying it to small-molecule LC-MS with a *contaminant fraction* emitted as a ranker feature is not published. | CHIMERYS (*Nat Methods* 2025) | absent from the competition pipelines |
| N6 | **Structure-space fragmentation graph (substructure nodes, not formula nodes)** with carbon-count + isotope filtering | **Substantially new.** Formula-space graphs are SIRIUS; Mass Spec has substructure annotation (MAGMa, DEREPLICATOR+) but only *given* a candidate structure. Using substructures as the *nodes of a tree computed without any candidate* is not published. | SIRIUS colourful subtree; MAGMa; DEREPLICATOR+ | — |
| N7 | **Mass-conserving assembly DP + exact-cover ILP for branched small molecules** | **The core novel claim.** Spectrum-graph DP (Dančík 1999) is *linear-only* and has never been extended to branching; assembly-index algorithms are *given* the molecule; MOLGEN/surge use *no spectral evidence*. The research synthesis states plainly that **no published algorithm constructs a structure from fragment masses by constraint satisfaction.** | Dančík spectrum graph (DP-over-DAG); Marshall et al. assembly index (DP over assembly states); Senior/Lewis; MILP formula assignment | FOAM (a *stochastic* genetic algorithm, formula-seeded, no spectral constraints) |
| N8 | **Alphabet inference (MAP alphabet by Gibbs) for branched small molecules, with mass→SMARTS bridge** | **Extension.** Kreitzberg et al. 2019 is the direct parent; n-ary condensation and the reaction-SMARTS dictionary are ours. | Kreitzberg et al., *J Proteome Res* 2019 | — |
| N9 | **Bond-ablation critic (B-critic)** | **New framing.** MolDiscovery's fragmentation graph + `P(logRank|bondType)` is the parent; using it as an *auditor that certifies which bonds must break, with infeasibility scoring*, rather than as a generator score, is new. | MolDiscovery; DEREPLICATOR+; ICBERG | — |
| N10 | **Hard logical masks applied as row filters, with diagnostic-substructure bitmaps** | **Adapted.** Seven Golden Rules, Senior/Lewis, nitrogen rule are all published as *formula* filters; applying substructure-level hard masks as bitmap conjunctions at pool level is not. | Kind & Fiehn 2007; Senior; LEWIS | the field treats these only as formula assignment steps |
| N11 | **Transductive propagation over the 400 test molecules** | **New here.** Spectral-networking + label propagation are both standard; combining them as a *re-ranking stage inside a competition submission*, with an evidence-based refusal guard, is not. | GNPS networking; Zhou label propagation | absent |
| N12 | **Posterior-merged (not additively post-hoc) popularity prior**; class-conditional μ | **Correctness improvement, not novelty.** | published μ·(log n_sub + log n_pm) form | `docs/03` §5.2 |
| N13 | **Slice B: a 100 000-molecule synthesised Class-2 validation set** | **New methodology.** Every published method validates on NPLIB1 or MassSpecGym, both ~10³–10⁴. | Seven Golden Rules construction; scaffold clustering | `docs/03` §7 gives only 250 |
| N14 | **Posterior-sorting derivation of the MRR@25-optimal submission** | **Not new — a proof.** Included because it invalidates the "diversify the tail" heuristic that appears in most advice, and it makes tail design deterministic rather than ad hoc. | — | — |

## 10.1 The three claims worth defending

**N1 (two-sided scoring).** The argument is airtight and short: a score that is monotone in
the number of predicted peaks cannot rank two candidates that predict overlapping-but-different
peak sets. MetFrag's own published ablation shows the failure empirically (depth 2 → 3 raises
coverage and lowers accuracy). The literature has all the pieces (assignment problems,
Earth-mover distance, precision/recall of prediction); they have simply never been assembled
into the primary ranking score of a metabolite-ID pipeline.

**N3 (DSM).** The novelty is in *what you do with the offset*. Everyone already computes Δm
for every analog hit — it is a required output of modified cosine. Using it as a chemical
operation specification is one step, and the step has never been taken. It is also
mechanically the best answer to Class 2, which is the class the field handles worst.

**N7 (assembly).** This is the highest-risk, highest-novelty item. Its honest defence is
negative: the research synthesis found *no published algorithm that starts from fragment
masses and constructs a structure by constraint satisfaction*, and found that the two nearest
systems are both linear-only or evidence-free. The risk is that it does not perform; §11
gives the kill criterion.

---
# Part 11 — Expected Score, Risk, Kill Criteria

## 11.1 Expected progression

Deliberately conservative. The "documented public pipeline" column is what `docs/03` §5.1
actually measured; the Plan 2 column is our estimate for the same stage plus the algorithmic
additions.

| Stage | Added | Est. MRR@25 | Documented public pipeline |
|---|---|---|---|
| 0 | mass filter + popularity fill | 0.05–0.10 | 0.087 baseline |
| 1 | formula ILP + curation + **replicate consensus** | 0.12–0.16 | — |
| 2 | Channel L (library entropy) | 0.18–0.22 | — |
| 3 | **Channel E** (energy-rescaled offset) | 0.28–0.33 | 0.335 (analog propagation) |
| 4 | **Channel F** (two-sided) + T + N | 0.33–0.38 | 0.358 |
| 5 | ranker + **hard masks** + **Channel Δ DSM** | 0.38–0.42 | 0.409 |
| 6 | **Channel S** (sparse mixture) + **B-critic** | 0.42–0.45 | 0.420 |
| 7 | **calibration** + posterior tail + **Channel X** | 0.45–0.49 | 0.47 (board top) |

Realistic target: **0.44–0.48**, i.e. competitive with the current board top, with the
algorithmic additions contributing most of the delta over a faithful re-implementation of the
public pipeline.

**Where the delta comes from, in order of expected contribution:**

```
  +0.015–0.030   Channel E  (cross-instrument library becomes usable)
  +0.010–0.020   Channel F two-sided (isomer discrimination the field cannot do)
  +0.008–0.015   Channel Δ DSM (Class-2 candidates nobody else can reach)
  +0.005–0.012   replicate consensus + Channel S (cleaner evidence, less chimeric noise)
  +0.004–0.010   calibration + posterior tail (free points from doing the metric properly)
  +0.004–0.010   hard masks (recovers candidates a soft model would rank low)
  +0.002–0.008   Channel X transduction (real but bounded; refuse-on-high-confidence is
                 necessary or it back-fires)
  +0.005–0.015   Track B assembly  (Class 3 only, so small; likely < 3 % of molecules)
```

## 11.2 Risk register

| Risk | P | Impact | Mitigation | Early warning |
|---|---|---|---|---|
| **Validation trap re-emerges in Slice B** — Slice B is Orbitrap/QTOF and we validate on it while testing on timsTOF | High | Very high | Slice A (250 timsTOF) remains the arbiter of large decisions; Slice C scaffold-disjoint audit; the η model is fitted on Slice A and evaluated on Slice B, never the reverse | If Slice A and Slice B disagree about a component's sign, believe Slice A |
| **Hard masks kill correct answers** | Medium | High | Every mask logs its firing; `false_kill_rate` is a first-class metric; any mask that ever excludes the correct answer gets relaxed, not removed | non-zero false_kill on Slice A mask 3 |
| **Formula ILP wrong at high mass** (> 600 Da) | High | High | Hypothesis *set*, not a point; pool = union over up to 8 hypotheses; rungs 4–6 recover adduct/near-miss errors | `formula_recall@1` vs mass plot on Slice B |
| **Channel E is a fiction** — η is unfittable and the correction is noise | Medium | Medium | E5 experiment in week 2; demote to a feature if Δ< +0.005; the uncorrected offset search still works as a fallback | Slice B Δ by library |
| **Two-sided F does not beat one-sided** | Low–Med | Medium | E4 experiment in week 2; fall back to `recall − precision_penalty` | Slice A mask 2 |
| **DSM generates garbage that floods the pool** | Medium | Medium | Cap 720 raw → 400 after dedup → only instantiated for top-50 base candidates when `lib_max < 0.75`; `origin` is a feature so the ranker can learn to distrust it | pool size after dedup; rank-1 share of DSM candidates |
| **Track B assembly is too slow** | Medium | Low | Top-40 molecules only; the DP is warm-started into the ILP; a full run is 12 min on 8 cores | wall clock at n=40 |
| **Track B assembly produces valence-legal nonsense** | High | Low | `strain` + `SA score` priors in §5.6; cap rank-1 probability at 0.35 unless ablation_fit ≥ 0.85; the library gate means it cannot displace a confident measurement | share of B outputs in the top 5 |
| **Transduction back-fires** (a whole cluster of wrong answers made more confident) | Medium | Medium | α = 0.35; seed never overwritten; refusal guard at lib_max ≥ 0.8; A/B on the leaderboard before keeping | leaderboard delta of the ablation |
| **Calibration overfits 250 molecules** | High | Medium | Fit temperature + isotonic on Slice B (~100 k), read ECE on Slice C | ECE drift between slices |
| **Kaggle packaging failure** (internet off, 9 h limit, missing artefact) | Low | Fatal | Pre-upload everything; a 60-min budget against a 9 h limit; dry-run on the sample test set | dry-run diff |
| **RAM blowup on a 2 M-structure pool** | Low | High | The memory ledger is a design constraint, not an estimate; index everything on disk | RSS monitor |

## 11.3 Kill criteria (pre-committed, so we do not argue with ourselves in week 19)

```
  KILL Channel E        if ΔMRR on Slice B < +0.005 after full tuning → demote to feature
  KILL Channel T        if ΔMRR < +0.002 on Slice A mask 2 → remove (costs 4 min/run)
  KILL Channel S        if ΔMRR < +0.003 on Slice B4 (low-quality) → remove
  KILL Channel X        if leaderboard A/B delta < +0.003 → remove (no reuse of the probe)
  KILL DSM              if pool_recall@800 on Slice B improves < +2 points → demote to
                        tail-filler only, drop the SMARTS-editing machinery
  KILL Track B assembly if assembly_recall@25 on Slice A mask 3 < 5 % after week 17 →
                        ship B-critic only, document the negative result, do not ship a
                        generator that does not work
  KILL any hard mask    if false_kill_rate > 0 on Slice A mask 1 → soften to a feature
  HALT ALL TUNING       at week 20. After that, changes are bug fixes only.
```

## 11.4 What we are deliberately not doing, and why

```
  ✗  Training a big transformer on 2.5 M spectra.
     Documented value: the neural channel alone scored 0.145 versus 0.335 for library
     matching. Its marginal feature contributions inside a 200-feature ranker are within noise.

  ✗  A diffusion / flow-matching de novo generator.
     DiffMS 4.1 % top-10 on leakage-controlled MassSpecGym; MARLIN 16.9 % top-1 formula-free;
     MS-Flow 32–45 % top-1 — all on benchmark settings that are *not* this test distribution,
     and all requiring GPU days inside a 9-hour notebook budget. The EV is worse than DSM's,
     the cost is 100× higher, and it violates the algorithmic constraint you set.

  ✗  Optimising on the public leaderboard more than a handful of times.
     400 molecules; the CI on a leaderboard read is roughly ±0.02. Six probes buys a real
     effect; thirty buys noise.

  ✗  Hand-writing a diagnostic-ion table beyond ~25 entries.
     The substructure-bitmap engine (§4.9) generalises automatically to any substructure we
     can draw, and the edit dictionary is mined rather than written. Hand tables do not scale
     and they encode one chemist's memory of one compound class.
```

---
# Appendix A — Schemas, Pseudocode, Hyperparameters

## A.1 On-disk artefacts

```
peaks.mzs      float32[mzs_len]      concatenated, mmapped
peaks.int      float32[ints_len]     concatenated, mmapped
peaks.off      int64[n_spec + 1]     row-group-aligned offsets
peaks.spec_id  int32[n_spec]         index into meta
meta.parquet   spectrum_id, molecule_id, precursor_mz, adduct, mode, CE_ev, BPI,
               ingest_lib, inchikey14, smiles, formula, num_peaks, precursor_error_ppm

fp2048.npy     uint8[n_struct, 256]  packed Morgan r=2, 2048 bits
structs.parquet  inchikey14, smiles_flat, exact_mass, formula, db_flags,
               n_substances, n_pubmed, n_coconut, has_ref_spectrum, source_mask

massidx.npy    float64[n_struct]     sorted ascending; rung windows = searchsorted
formidx.npy    int32[n_struct]       formula-hypothesis id per structure
lib_index.parquet  per library: spec ranges, η, instrument family, unit family
masks/*.npy    uint64[ceil(n_struct/64)]  one bitmap per diagnostic substructure
alphabets/*.json   the B1 alphabet output, per energy bucket
```

## A.2 Hyperparameters (the complete initial set; all later changes logged)

```
 ── masses ───────────────────────────────────────────────────────────
   tau_prec_base      8.5 ppm     |  tau_prec_wide   25 ppm
   rung_3             100 ppm     |  tau_spec        0.005 Da
   isotope_merge_tol  0.0015 Da   |  deiso_ratio_max 0.15

 ── curation ─────────────────────────────────────────────────────────
   rel_intensity_floor      0.005
   window_topN              6 peaks / 50 Da     |  global_topN 128
   min_peaks_after_curation 4
   int_transform            raw + sqrt (both kept)

 ── replicate consensus ─────────────────────────────────────────────
   vote_threshold          0.34   (peak must appear in ≥ 1/3 of bucket replicas)
   intensity_cv_max        0.85
   peak_weight             replicate_fraction ** 0.5
   aggregation over spectra: 0.6 * quantile(sim, 0.2) + 0.4 * mean(sim)

 ── Channel L ───────────────────────────────────────────────────────
   entropy_dyn_weight     w = 1 if S>=3 else 0.25 + 0.25*S
   lib_conf_threshold     0.60    (counts toward n_lib)
   lib_margin             cohort = rung-1 candidates of the same molecule

 ── Channel E ───────────────────────────────────────────────────────
   eta_prior              anchored: timsTOF-eV identity
   c_delta                0.02 Da per eV of CE mismatch
   energy_gap_scale       20 eV    (exponential decay)
   slice A used for the fit; riken/gnps/mona/spectraverse/massbank/msdial for evaluation

 ── Channel Δ (DSM) ─────────────────────────────────────────────────
   n_retrieved_refs       30
   n_delta_hypotheses     3
   n_attachment_sites     8
   dict_min_instances     50
   dict_min_tanimoto      0.60
   gate: only if lib_max < 0.75

 ── Channel F ───────────────────────────────────────────────────────
   fragmenter depth       2  (bridges + 2-cuts)
   MetFrag score exponents alpha 1.84, beta 0.59, gamma 0.47
   F_beta                 1.0
   survival s(E,E0)       1 / (1 + exp(-(E - E0)/6)), E0 = 8 + 1.2*n_heavy
   gamma_gate             Λ 0.70, κ 0.05

 ── Channel T ───────────────────────────────────────────────────────
   FT loss prior          SIRIUS-trained, or empirical from Slice A if unavailable
   max colour reuse       1 (colourful subtree)

 ── Channel N ───────────────────────────────────────────────────────
   fingerprint            Morgan r=2, 2048 bits
   negative-bit threshold p > 0.7

 ── Channel S ───────────────────────────────────────────────────────
   lambda_1               0.02 * ||I_obs||_1      (sparsity)
   lambda_2               0.10                    (anchor to the evaluated candidate)
   coordinate descent     200 iterations, tol 1e-4

 ── Channel H ───────────────────────────────────────────────────────
   halogen detection      M+2/M >= 0.25 → Cl ; >= 0.80 → Br ; 0.030–0.070 → S
   diagnostic substructures: 91.0542/C7H7, 105.0335/C7H5O, 93.0335/C6H5O,
                             179.0344/C9H7O4, 121.0290/C7H5O2, 135.0441/C8H9O,
                             149.0234/C9H5O3, 163.0391/C9H7O4, 147.0441/C9H7O2,
                             162.0528/C6H10O5, 109.0289/C6H5O3, 191.0342/C10H7O4

 ── Prior ───────────────────────────────────────────────────────────
   mu_pop  0.20  (leaderboard-tuned, ≤3 probes)
   mu_np   0.15  (COCONUT/ChEBI hits)
   mu_lib  0.10  (class-1 indicator)

 ── Ranker ──────────────────────────────────────────────────────────
   LightGBM lambdarank, ndcg_eval_at [1,5,25], 300 trees, lr 0.05,
   min_data_in_leaf 20, feature_fraction 0.8, lambdarank_truncation_level 30
   training regime mix  Class-1 : Class-2 : Class-3 = 45 : 40 : 15

 ── Calibration ─────────────────────────────────────────────────────
   T (softmax temperature)  fitted on Slice B
   decoys per real candidate 20, formula- and mass-matched, non-identity
   isotonic                 on Slice B, ECE reported on Slice C

 ── Track B ─────────────────────────────────────────────────────────
   alpha_ary            3   (max arity considered in alphabet inference)
   max_subset_span      250 Da
   min_block_mass_uses  60 Da  (below this a block may be unobserved)
   n_blocks_max         MA_upper_bound + 2
   B_critic_depth       2
   B_rank1_prob_cap     0.35 unless ablation_fit >= 0.85
   lambda_bic           1.0   | gamma_MA  0.05  | gamma_strain 0.4  | gamma_rarity 0.3
```

## A.3 Feature families (for the ranker; ~200 columns)

```
 POOL      n_candidates, rung_membership(6), origin(lib/dsm/denovo), has_ref_spectrum,
           n_independent_libs, exact_mass - M_neutral (ppm), formula_id, is_F0,
           passes_H1..H10 (bits), n_masks_fired

 MEASURED  lib_max, lib_top2, lib_margin, lib_score, n_lib, sim_max_E,
           n_ref_within_5eV, energy_gap_min, energy_gap_mean,
           analog_coherence, delta_spread, delta_top1, delta_top1_signed,
           n_analog_hits, n_analog_from_timsTOF, ft_align_score, ft_unexplained_edges

 FORWARD   frag_F1, frag_recall, frag_precision, frag_LLR,
           frag_n_explained, frag_basepeak_explained, frag_ladder_len, frag_gap_only,
           iberg_sim, glac_sim, iberg_basepeak, glac_basepeak,
           energy_weighted_frag_score

 CRITIC    ablation_fit, ablation_pen, ablation_necessity, ablation_n_novel_bonds,
           ablation_fit_LLR, ablation_chain_consistency

 SPARSE    chi_P, resid_P, resid_uni, P_gain, contam_frac

 NEURAL    nbits_dot, nbits_loglik, nbits_jsd, nbits_hard_negatives

 QUERY     (broadcast from the molecule, constant across candidates)
           n_spectra, bpi_median, bpi_min, mode, n_adducts, ce_mean, ce_spread,
           n_peaks, replicate_mean, chimeric_fraction, has_Cl, has_Br, has_S

 PRIOR     log1p(substances), log1p(pubmed), log1p(coconut), source_mask_bits
```

## A.4 Reference card — the metric in five lines

```
  MRR@25 = mean over molecules of  1 / rank_of_first_correct   (0 if absent)
  correct  ⟺  RDKit-tautomer-canonical InChIKey14 equal        (stereo free)
  no penalty for wrong guesses  ⟹  ALWAYS SUBMIT EXACTLY 25
  MRR is maximised by sorting candidates by  p(correct | evidence)  descending   [proof §6.3]
  therefore:  the tail must be probability-sorted, NOT structurally diversified
```

## A.5 File layout of this project

```
 plan2.md                     this document
 docs/                        the source material (read-only reference)
 src/
   data_layer/   ingest.py  adducts.py  curate.py  consensus.py  formula_ilp.py
   pool/         structures.py  rungs.py  masks.py  substructure_index.py
   channels/     lib_L.py  energy_E.py  dsm_delta.py  forward_F.py  tree_T.py
                 sparse_S.py  neural_N.py  critic_B.py
   track_b/      alphabet.py  structgraph.py  assembly_dp.py  assembly_ilp.py
                 bonds.py  model_select.py
   fusion/       ranker.py  calibrate.py  transpose.py  merge_gate.py
   validate/     slices.py  metrics.py  ablation.py  calibrate_check.py
   submit/       write.py  validate_format.py
 artifacts/      everything from Appendix A.1
 notebooks/      the Kaggle submission notebook (paste entry point, ≤ 9 h)
```

## A.6 The ten questions to answer before week 3

```
  1. On Slice A, how often is the top formula hypothesis correct?        (target ≥ 0.92)
  2. On Slice B, how many candidates does rung 1 leave, median and p90?  (if p90 > 2000, rungs matter)
  3. Does replicate consensus improve Slice B multi-spectrum MRR?        (§3.3 is the hinge)
  4. What fraction of Slice B molecules have ≥ 2 spectra?               (determines §3.3's reach)
  5. Is the observed M_spread_ppm consistent with in-source families?    (§3.1.2 recovery rate)
  6. What is the false-kill rate of hard mask H7 (diagnostic ions)?      (must be ~0)
  7. How correlated are Channel L and Channel F across candidates?      (Kendall τ; < 0.6 ⇒ keep both)
  8. How many test molecules have halogen isotope evidence?             (sizes H4's reach)
  9. What is the evidence ceiling on Slice A mask 2 (oracle ranker)?    (bounds all ambition)
 10. What is the base_peak_intensity distribution, and does it predict
     spectral quality as well as the docs claim?                       (validates the weight)
```

---

