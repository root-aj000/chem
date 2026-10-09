# Part 3 — How People Actually Solve It

> **What this document is:** the current state of the Enveda CASMI 2026 competition
> (open 14 Sep 2026, entry deadline 7 Dec 2026, final submission 14 Dec 2026), the exact
> architecture behind the top public solutions with their real numbers, the published
> research methods from the literature, and a build order you can follow.
>
> **Who this is for:** someone who has read [Part 1](01-what-is-the-problem.md) and
> [Part 2](02-solving-it-by-hand.md) and now wants to know what actually wins.

---

## Table of contents

- [1. Where the competition actually stands](#1-where-the-competition-actually-stands)
- [2. The single most important structural fact](#2-the-single-most-important-structural-fact)
- [3. The winning architecture: retrieval, not generation](#3-the-winning-architecture-retrieval-not-generation)
- [4. The four evidence channels in detail](#4-the-four-evidence-channels-in-detail)
- [5. A full public solution, phase by phase, with real scores](#5-a-full-public-solution-phase-by-phase-with-real-scores)
- [6. The five ideas that produced most of the score](#6-the-five-ideas-that-produced-most-of-the-score)
- [7. The validation trap — the most important operational fact](#7-the-validation-trap--the-most-important-operational-fact)
- [8. What the researchers are doing: de novo generation](#8-what-the-researchers-are-doing-de-novo-generation)
- [9. The leakage controversy you must know about](#9-the-leakage-controversy-you-must-know-about)
- [10. The forward direction: predict the spectrum](#10-the-forward-direction-predict-the-spectrum)
- [11. The toolbox](#11-the-toolbox)
- [12. A recommended build order](#12-a-recommended-build-order)
- [13. Common mistakes](#13-common-mistakes)
- [14. The one-paragraph summary](#14-the-one-paragraph-summary)

---

## 1. Where the competition actually stands

### 1.1 The leaderboard

Snapshot from a public standings aggregator, October 2026:

| Rank | Handle | MRR@25 | Last submission |
|---|---|---|---|
| 1 | `pikachu jzhpfn` | **0.47** | ~18 h ago |
| 2 | `MarvinTMB marvintmb` | **0.47** | ~1 day ago |
| 3 | `Nicolas Krusche nicolaskrusche` | 0.44 | ~17 h ago |
| 4 | `Bertan Sezgin bertansezgin` | 0.44 | ~19 h ago |
| 5 | `Shehab Alshehabi` | 0.44 | ~20 h ago |
| 6 | `Sho Saga shosaga` | 0.44 | ~1 day ago |
| 7 | `chopper mishrasiddhant` + team | 0.44 | ~1 day ago |
| 8 | `Arsokan ilkeryasin8` + team | 0.44 | ~17 h ago |

The board is tightly clustered at the top — **0.44 to 0.47 across eight teams.** That
clustering is informative: several teams converged on essentially the same approach, and
nobody has found a step change yet.

### 1.2 The public baselines

Not everyone is at 0.47. Public Kaggle notebooks span a wide range:

| Public notebook | Score |
|---|---|
| `⚡ [0.336+ SOTA] Enveda CASMI26 Analog Ranker` | **0.319 – 0.336** |
| `CASMI'26: 4-Channel Molecular Retrieval & Ranking` | ~0.33 |
| `CASMI26: Quad-Channel Evidence Ranker` | ~0.33 |
| `Enveda CASMI 2026 V7 Fast SOTA Bagging` | ~0.33 |
| `🧬 From Mass to Molecule: A Transparent Retrieval` | 0.17 |
| The official tutorial notebook (de novo transformer) | low |

### 1.3 Look at the notebook *titles* — they tell the whole story

```
  CASMI26: Four-Channel Ranker, Tautomer Dedup
  CASMI'26: 4-Channel Molecular Retrieval & Ranking
  CASMI26: Quad-Channel Evidence Ranker
  CASMI26 20260919 Two Ranker W80
  CASMI26 20260919 Two Ranker W65
  CASMI26 20260919 Two Ranker W50
  CASMI26 20260919 Two Ranker W35
  Analog Propagation - CASMI 2026 baseline
  CASMI26 One Engine Nine Scores
  Enveda CASMI26 Analog Ranker (community method)
```

Read those names and you can reconstruct the entire field's convergence:

- Nearly every title contains **"retrieval," "ranker," or "analog."**
- **None** contains "diffusion," "autoregressive," "SMILES decoder," or "generation."
- There are at least **six variants** of "Two Ranker W##" — people tuning a weighted blend
  of two gradient-boosted rankers.

> **The field converged on retrieval-plus-learning-to-rank within about two weeks, and the
> convergence was total.** The de novo generative approach that the official tutorial
> teaches is not competitive. Understanding *why* is more useful than any single trick.

---

## 2. The single most important structural fact

Three facts from [Part 1 §11](01-what-is-the-problem.md) combine into the conclusion that
determines the entire architecture:

1. **Class 1** (in a public spectral library) and **Class 2** (structure in PubChem or
   COCONUT but no public spectra) are *retrieval* problems. The answer **exists in a
   database**.
2. **MRR@25 with no penalty for wrong guesses** means a correct shortlist scores well even
   when imperfectly ordered.
3. The candidate pool — training structures + COCONUT — is only **~712,000 structures.**

**712,000 structures is not a large number.** A modern laptop can hold the fingerprints of
all of them in a few gigabytes. You can compute, for every single one, how well it matches
the test spectrum. That is a tractable, exhaustive problem.

**So the correct frame is:**

> This is **not** "invent a molecule from a spectrum." It is **"search 712,000 molecules
> and rank them."**

The molecules that are genuinely novel (Class 3) are the minority, and for those you will
lose points regardless of what you do — the published state of the art is under 5% exact
recovery. **Do not spend your effort there.**

---

## 3. The winning architecture: retrieval, not generation

The pipeline that the top public solutions converged on:

```
              ┌──────────────────────────────────────────┐
              │  ALL SPECTRA OF ONE MOLECULE            │
              │  (median 3, up to 16)                    │
              └────────────────────┬─────────────────────┘
                                   │
                                   ▼
              ┌──────────────────────────────────────────┐
              │  STEP 1 — DEDUCE THE MASS                │
              │  M = precursor_mz ∓ adduct               │
              │  Cross-check all spectra agree           │
              └────────────────────┬─────────────────────┘
                                   │
                                   ▼
              ┌──────────────────────────────────────────┐
              │  STEP 2 — CANDIDATE POOL (~712,000)      │
              │  training structures ∪ COCONUT           │
              │  (later + PubChem)                       │
              │  FILTER: within 8.5 ppm of target mass   │
              └────────────────────┬─────────────────────┘
                                   │  typically leaves
                                   │  tens to hundreds
                                   ▼
              ┌──────────────────────────────────────────┐
              │  STEP 3 — FOUR EVIDENCE CHANNELS         │
              │   (a) library   spectral similarity      │
              │   (b) analog    adduct-shifted search    │
              │   (c) fragment  bond-break coverage      │
              │   (d) model     fingerprint logits       │
              └────────────────────┬─────────────────────┘
                                   │  ~30–160 features
                                   ▼
              ┌──────────────────────────────────────────┐
              │  STEP 4 — GRADIENT-BOOSTED RANKER        │
              │  LightGBM / CatBoost  →  one score       │
              │  + popularity prior                      │
              │  + library gate                          │
              └────────────────────┬─────────────────────┘
                                   │
                                   ▼
              ┌──────────────────────────────────────────┐
              │  STEP 5 — DEDUP ON InChIKey14            │
              │  Write top 25 as semicolon-joined SMILES │
              └──────────────────────────────────────────┘
```

### 3.1 Why the mass filter is the load-bearing step

**8.5 ppm on a 350 Da molecule is 0.003 Da.** With that filter, 712,000 candidate
structures collapse to a few dozen. Every subsequent stage operates on a tiny list, which
is why the whole thing fits in a 9-hour notebook.

This filter is pure arithmetic. It is also, per [Part 2 Step 4](02-solving-it-by-hand.md),
almost always sufficient to pin the **molecular formula** exactly. Nothing clever happens
here — and it is the highest-value step in the pipeline.

### 3.2 Why a gradient-boosted ranker and not a neural net for the final decision

After the mass filter, the remaining problem is: *given 30–160 heterogeneous numerical
features about each of ~50 candidates, put them in the best order.*

That is a **tabular learning-to-rank problem.** It is LightGBM's home turf.

And a GBDT has decisive advantages here:

| Advantage | Why it matters |
|---|---|
| **Handles mixed feature types** | Spectral entropy scores, Tanimoto values, fragment coverage fractions, counts — no common scale needed |
| **Works on ~50 rows** | You have tiny training sets after filtering. Deep nets would not train. |
| **No normalisation assumptions** | Distances and similarities are wildly skewed; trees do not care. |
| **Feature importance for free** | You can see which evidence channel is carrying the score, which tells you what to work on next. |
| **Ships as a text file** | Fits the Kaggle notebook constraint trivially. |

---

## 4. The four evidence channels in detail

This is the part that matters most, because it is where the actual chemistry insight lives.
Each channel is a **completely independent** line of evidence, and combining four weak
signals beats one strong one.

### CHANNEL A — Library channel: "does a reference spectrum already exist?"

For each candidate structure, look up all its training spectra. Compute a **spectral
similarity score** against the query spectrum. Similarity is typically **spectral entropy
similarity**, which is much better than naive cosine at handling mismatched peak counts.

```
  lib_max  =  max over the candidate's training spectra of  similarity(train_spectrum, query)
```

- If `lib_max ≈ 1.0`, the answer is almost certainly correct. This is **Class 1** and it
  should go to rank 1.
- If `lib_max ≈ 0`, the candidate is in the pool but has no reference spectrum — Class 2.

**Why spectral entropy, not cosine?** Two spectra of the same molecule rarely have
identical peak counts. Entropy-based similarity (Böcker & Młotkowski 2012) compares the
*intensity distributions* and is far more robust to missing peaks than cosine.

### CHANNEL B — Analog channel: "is there a closely-related molecule with a reference?"

**This is the single cleverest idea in the whole pipeline and most people miss it.**

The problem: the query's structure isn't in any library, but a *similar* structure might be.

The trick: **spectra of similar molecules are often offset by a constant mass difference.**
If candidate X has a reference spectrum and the query looks like X's spectrum but shifted by
the mass of, say, a methyl group (14.0157) or an oxygen (15.9949), that shift tells you the
query is "X plus a methyl."

```
  1. Take the query spectrum.
  2. Compute modified-cosine similarity to EVERY training spectrum,
     allowing a constant mass offset.
  3. For the best matches, note which candidate structures own
     those training spectra.
  4. Compute Tanimoto similarity between the query's implied
     fingerprint and those candidates.
```

A **modified cosine** score is used because it allows for the small systematic m/z offset
between different instruments — essential when the training set spans dozens of instrument
types and the query is timsTOF.

Implementation notes from the successful versions:
- **Adduct-shifted matching** — align the precursor masses before comparing, otherwise every
  cross-adduct pair is garbage.
- Two variants proved valuable: an **"AFIX"** variant restricted to timsTOF spectra, and a
  **"BIO"** variant pulling from ChEBI and LIPID MAPS.
- Blending both beat either alone.

> **Intuition:** if you don't know what compound this is, but you *do* know that compounds
> with a 3,4-dihydroxyphenyl fragment all fragment the same way, then finding that
> fragment tells you the answer is in that family — and you rank the family's members first.

### CHANNEL C — Fragmentation channel: "can the candidate's bonds explain these peaks?"

Run a **bond-dissociation fragmenter** on the candidate and ask how much of the observed
peak intensity it can account for.

The best-known open implementation is **MetFrag**; a simplified version is
**MetFrag-lite**.

```
  frag_score  =  (sum of intensities of observed peaks that the
                  candidate's single/double bond cleavages can produce)
                 /  (total observed intensity)
```

**Why this is a genuinely independent channel:** it never looks up any reference spectrum.
It asks a purely physical question about the candidate's own bond connectivity. A molecule
can be in the pool with zero training spectra (Channel A and B both dead) and still be
strongly supported by Channel C.

### CHANNEL D — Model channel: "does a neural network's fingerprint prediction agree?"

Train a **spectrum → molecular fingerprint** transformer. Then:

```
  model_score  =  dot_product( candidate's Morgan fingerprint,
                               network's predicted fingerprint logits )
```

This is the "modern neural network" component, and note how it is used: **not to generate
structures, but to produce one more similarity number for the ranker.** That is a much
easier learning problem than generation, and it is where the neural network earns its keep.

Implementation detail that mattered: training the fingerprint model on **COCONUT** rather
than the training set. COCONUT is a natural-products database, so it matches the test
chemistry better than Enveda's drug-like data.

### 4.1 The complete public pipeline, as published

From the community pipeline released by `prvsiyan` and `haideptry` under Apache 2.0, and
documented in the `shiferaxa/casmi26` repository:

> 1. **Candidates** are every structure within **8.5 ppm** of the target mass in a pool of
>    about **712k structures**: COCONUT plus every structure in the training set.
> 2. **Four kinds of evidence** are computed per candidate:
>    - **spectral entropy similarity** to training spectra of that structure (*library*)
>    - **Tanimoto similarity** to structures whose training spectra resemble the query
>      after a mass shift (*analog*)
>    - the **fraction of peak intensity explained by one or two bond breaks** of the
>      candidate (*fragmentation*)
>    - the **dot product of the candidate fingerprint with the logits of a
>      spectrum-to-fingerprint transformer** (*model*)
> 3. A **gradient-boosted ranker** turns the evidence into a probability.
> 4. Candidates are **deduplicated by canonical key**, and the **top 25** are written.

That is four channels, a GBDT, a dedup, and a 8.5 ppm mass filter. **Everything that scores
above 0.40 on this leaderboard is a variation on that skeleton.**

---

## 5. A full public solution, phase by phase, with real scores

One participant published a complete phase log. It is the most useful single artefact for
understanding this competition, because it shows **the exact order in which the score
arrived** — and that order is itself the lesson.

### Phase 1 — Neural fingerprint model alone

```
  Morgan + MACCS fingerprint transformer, trained on COCONUT,
  decoder predicting SMILES from the predicted fingerprint.

  0.087  →  0.145
```

**The plateau is the finding.** Pure neural prediction ran out of road almost immediately,
and the reason is stated plainly: *lack of reference matching.* The model can only say "a
structure like this might fit," which is worth very little when the answer is a specific
member of a specific family.

### Phase 2 — Analog propagation

```
  Adduct-shifted modified cosine + MetFrag-lite + "Ranker A"

  0.332  →  0.335
```

**+0.19 in one step — the biggest single jump in the entire log.** Adding reference-spectrum
matching to the neural model took the score from 0.145 to 0.332.

That ratio is the whole competition in one number: **reference matching is worth about
twice what the neural network was worth.**

### Phase 3 — Two-ranker ensemble

```
  Channel BIO  = ChEBI + LIPID MAPS derived
  Channel AFIX = timsTOF derived
  Ranker A (31 features)  +  Ranker B (51 features), blended

  0.350  →  0.358
```

Two specialised rankers beat one general one. The "W##" suffixes you saw in the notebook
titles (`W80`, `W65`, `W50`, `W35`) are people tuning this blend weight.

### Phase 4 — Deep PubChem recall + richer features

```
  fe_v4 feature families + fpnet_full1 + PubChem recall N_1=5000

  0.358  →  0.384
```

### Phase 5 — Engine fusion + union rescoring

```
  The top-40 candidates from EACH engine are pooled into a union.
  Every candidate in the union is then re-scored by the forward
  spectrum predictors ICEBERG and GLACIER (MIT ms-pred).
  Final score = z(ranker) fused across engines.

  0.384  →  0.399
```

The innovation here: **most re-ranking schemes only re-score the existing list.** This one
first *unions* the lists, letting a candidate proposed by one engine be validated by
another engine's evidence.

### Phases 6–7 — Popularity prior and library gate

```
  Popularity prior mu=0.15                    0.399 → 0.404
  Library gate: ICE_LAM_LIB=0, GL_LAM_LIB=0   0.404 → 0.404  (held, but safe)
  Stronger prior POP_MU=0.25                  0.404 → 0.409
```

### Phase 8 — PubChem join

```
  PC_JOIN_N = 50, PC_JOIN_MAX_LIB = 0.7

  0.409  →  0.420
```

### 5.1 The score table

| Phase | What was added | Cumulative LB |
|---|---|---|
| 0 | baseline (neural fingerprint only) | 0.087 |
| 1 | better fingerprint transformer | 0.145 |
| 2 | **analog propagation + MetFrag-lite** | **0.335** |
| 3 | two-ranker blend, ChEBI/LIPID/timsTOF | 0.358 |
| 4 | deep PubChem recall, more feature families | 0.384 |
| 5 | engine fusion + ICEBERG/GLACIER union rescoring | 0.399 |
| 6 | popularity prior μ=0.15 + gap fragment rescoring | 0.404 |
| 7 | stronger prior μ=0.25 + library gate | 0.409 |
| 8 | PubChem join (N=50, max_lib=0.7) | **0.420** |

The leaderboard top is 0.47, so roughly 0.05 more sits in solutions not publicly
documented.

### 5.2 The popularity prior, written out

From the published solution. When candidates tie on spectroscopic evidence, break the tie
by how often the compound appears in the literature:

```
  final_score  =  z(ranker_score)
               +  μ · ( log(1 + n_pubchem_substances) + log(1 + n_pubmed_citations) )
```

with `μ = POP_MU`, tuned from `0.15` to `0.25`.

Applied only to candidates that came from a structure database.

**Why this is legitimate and not cheating:** it encodes a genuine prior over the label
distribution. Nature products that have been isolated, named, and published are far more
likely to be the answer in a test set built from a realistic mixture than an obscure
synthetic congener nobody has ever reported. The ground truth is drawn from a
chemistry-practitioner-drawn sample, and that sample is biased toward known compounds.

### 5.3 The library gate

```
  ICE_LAM_LIB  =  0.0
  GL_LAM_LIB   =  0.0
```

Meaning: **when the experimental library similarity is high (≥ 0.7), switch the forward
simulators off entirely.**

**Why:** ICEBERG and GLACIER *simulate* spectra from candidate structures. Their simulators
were trained on particular instruments and do not know your instrument is a timsTOF with a
specific fragmentation regime. So a genuinely correct library match — which came from
*real measured data* — would get **overruled by a simulation**. The gate says: never let a
simulation outrank a real measurement.

This is a nice general principle, not just a hyperparameter: **a well-sourced observation
should always beat a generated prediction.**

### 5.4 The PubChem join

```
  PC_JOIN_N       = 50
  PC_JOIN_MAX_LIB = 0.7
```

Previously, candidates found in PubChem but absent from the train∪COCONUT pool were
**pinned to fixed slots** (4th, 8th, 12th, …) in the final list. This version instead
**injects the top 50 PubChem-only candidates directly into the main engine's ranker**, where
they get all ~160 features and compete normally — with a gate protecting molecules that
already have a confident library match.

The gain: **+0.011**, because outside-pool answers could finally reach rank 1.

> **General lesson that transfers to any ranking competition:** any design decision that
> *structurally prevents* a class of candidates from ever ranking first is a hard ceiling
> on your score. Find those decisions and remove them.

### 5.5 The operational reality of shipping it

From the published run instructions:

```
  Accelerator   : GPU  T4 x2
  Internet      : Off
  Attached data : 14 datasets (including casmi26-pubchem-popularity-prior)
  Entry point   : paste kaggle_artifacts/sota_top1_submission.py into a notebook
  Runtime       : ~15–20 minutes smoke check on 12 molecules
  Output        : /kaggle/working/submission.csv
```

Two things to note:

1. **It runs in 20 minutes, not 9 hours.** The 9-hour notebook limit is generous; nobody
   needs it. Retrieval is cheap. (The *tutorial's* de novo model needs the time, because
   autoregressive generation without a KV cache is slow.)
2. **All 14 datasets must be attached as Kaggle input**, because internet is off. Anything
   you want — pre-trained weights, COCONUT, PubChem tables — has to be pre-uploaded. This
   is a real logistics constraint that bites people who plan to download anything at
   runtime.

---

## 6. The five ideas that produced most of the score

If you only remember five things:

### 1. Analog propagation via adduct-shifted modified cosine — worth ~+0.19

Search training spectra with a **mass offset allowed**, then use Tanimoto similarity on
the retrieved neighbours. Enormous. Easy to implement. Skipped by almost everyone at first.

### 2. Multiple independent evidence channels — worth ~+0.02 to +0.05

Library, analog, fragmentation, model. Four channels that fail independently. Combine them
in a GBDT.

### 3. The popularity prior — worth ~+0.01

`z(ranker) + μ·(log(1+substances) + log(1+pubmed))`, μ ≈ 0.25. Ridiculously cheap, nobody's
instinct, real points.

### 4. Never let a simulation outrank a measurement — protective

The library gate. Zero direct gain but prevents a category of regression.

### 5. Remove structural ceilings — worth ~+0.01

Let PubChem-only candidates compete for rank 1. Audit your own pipeline for anything that
confines a candidate class to fixed positions.

---

## 7. The validation trap — the most important operational fact

### 7.1 What happened

From the `shiferaxa/casmi26` repository, which tested this explicitly:

> The only split that has predicted the public leaderboard is the **250 structures the host
> measured on the test instrument** (`ingest_lib == enveda-np-examples`), evaluated under
> three masks. **Random held-out training structures predicted the opposite ordering of
> three rankers and were wrong on the board.**

Read that again. A random holdout did not merely fail to help — it **ranked the methods
backwards.** You could have optimised it perfectly and moved *away* from the answer.

### 7.2 Why it happens

Two mismatches compound:

| | Training set | Test set |
|---|---|---|
| **Instrument** | 10 of 12 libraries are Orbitrap/QTOF/etc. | **100% timsTOF** |
| **Chemistry** | `enveda-180` is 46% of all spectra, and it is **synthetic drug-like screening compounds** | **natural products** |

A random split from `train.parquet` is dominated by `enveda-180`. So you validate on
drug-like fragmentation measured on the wrong instrument, and test on natural products
measured on the right one. Those are different problems.

### 7.3 The correct validation design

```python
val = train_df[train_df.ingest_lib == 'enveda-np-examples']   # 250 structures
```

Evaluate under **three masks**, because each answers a different question:

| Mask | What it measures | Maps to |
|---|---|---|
| Everything visible | Best case — molecule is in pool, spectra available | Class 1 |
| Structure's own spectra hidden | Molecule in pool, no reference spectrum | Class 2 |
| Structure **removed from the candidate pool** | Genuinely novel structure | Class 3 |

**Read the mask-2 and mask-3 scores separately.** They tell you what fraction of the test
set you are structurally incapable of solving. That is genuinely useful information about
where your remaining points are.

### 7.4 The honest caveat

250 structures is a small validation set. Differences of ±0.02 on it are probably noise.
Use it to catch large regressions and to compare clearly different approaches — **not** to
tune μ to three decimal places. For that, prefer the public leaderboard itself, which is
a real 400-molecule held-out set.

---

## 8. What the researchers are doing: de novo generation

For Class 3 — molecules not in any database — you must generate. This is genuinely hard
research, and the honest state of the art is **low single-digit to ~20% exact match**.

### 8.1 The generational lineage

```
  1st era  fingerprint-conditioned RNNs and decoders
           ├── MSNovelist (Stravs et al. 2022)   CSI:FingerID fingerprint → SMILES
           ├── Spec2Mol    (Litsa et al. 2023)    SMILES autoencoder latent + CNN encoder
           ├── MassGenie   (Shrivastava 2021)     ~400M-param transformer, 4.7M synthetic
           ├── MS2Mol      (Butler et al. 2023)   m/z tokenisation + BPE SMILES
           ├── MS2SMILES                          explicit chemistry constraints
           ├── Mass2SMILES                        VAE + TCN decoder, 83k spectra
           └── Neuraldecipher                     invert fingerprints to structures
                   │
                   ▼
  2nd era  formula-conditioned GRAPH DIFFUSION
           ├── DiffMS  (Cooley group 2025)        the canonical reference implementation
           └── MADGEN  (Wang et al. 2025)         SCAFFOLD-conditioned
                   │
                   ▼
  3rd era  large-scale diffusion language models / posterior querying
           ├── FRIGID           diffusion LM, hundreds of millions of structures
           ├── MSFlow          flow matching on SAFE sequences + CDDD embeddings
           ├── MARLIN          formula-FREE
           └── MS-GPT          molecule-LM + LoRA + posterior band
```

### 8.2 DiffMS — the reference implementation

**Architecture:** two stages.

- **Encoder:** the **MIST formula transformer** (Goldman et al.), which treats a spectrum as
  an unordered **set** of (m/z, intensity) peaks — note it deliberately avoids treating
  m/z values as a *token sequence*, because mass spectra are permutation-invariant and a
  sequence model gets that wrong. Each peak is embedded using a **predicted formula
  assignment** from SIRIUS, and a **set transformer** implicitly models pairwise
  **neutral losses** between fragments. The embedding at the precursor position is the
  structural condition *y*.
- **Decoder:** a **discrete graph diffusion** model (DiGress) on the molecular adjacency
  matrix, **hard-constrained by the heavy-atom composition of a known chemical formula.**

**The two key insights:**

1. **Formula as a hard constraint.** Fixing the atom inventory prunes the search space
   enormously, bounds output length, and permits hard filtering and re-ranking.
2. **Pretrain the decoder on fingerprint–structure pairs.** There are only ~tens of
   thousands of real spectrum–structure pairs, but **2.8 million** fingerprint–structure
   pairs (DSSTox, HMDB, COCONUT, MOSES). Because the decoder is conditioned on a
   fingerprint-like vector either way, pretraining transfers directly. **This is the single
   most important practical trick in the field** — it decouples decoder training from the
   scarcity of labelled spectra.

**Result:** **4.1% top-10 on MassSpecGym.** And note the robustness result: dropping from
the *ground-truth* formula to a *predicted* formula (MIST-CF) costs only 4.25% → 4.1%.
The model is remarkably insensitive to formula-prediction error.

### 8.3 MARLIN — formula-free (July 2026)

The sharpest critique of the field: **almost everything assumes you already know the
molecular formula**, which is itself predicted with substantial error — and which does not
exist at all for a genuinely novel compound.

MARLIN removes the oracle entirely:

- **Self-supervised encoder** (DreaMS) predicts a fingerprint straight from raw peaks.
- **Block-diffusion language model** generates candidates conditioned only on the
  fingerprint **and the measured precursor mass**.
- A **"provably safe mass-shell constraint"** keeps every candidate on-mass *without*
  fixing the atom inventory. It cannot exclude the true molecule.
- A **symmetric noise objective** absorbs encoder error; a **candidate-diversity
  mechanism** prevents all 25 candidates collapsing to the same structure.
- Candidates accepted by exact **ppm agreement** (τ = 10).

**Results on NPLIB1, formula-unknown setting:**

| Method | Top-1 exact | Top-1 Tanimoto |
|---|---|---|
| Prior SOTA **given the ground-truth formula** (MBGen) | 12.20% | — |
| **MARLIN (DreaMS, fully formula-free)** | **16.94%** | **0.55** |
| MARLIN (MIST encoder) | 19.18% | — |
| FRIGID | — | — |

It **beats every prior system while those systems are handed the true formula.**

And a lovely bonus result: MARLIN never computes a formula, yet the correct formula emerges
in its top-1 candidate for **76.7%** of spectra (**96.9%** of spectra where it returns a
candidate at all) — comparable to dedicated formula predictor MIST-CF at 82.3%, and far
above FRIGID at 32.6%.

### 8.4 MS-GPT — posterior querying (July 2026)

Identifies a specific, general pathology:

> The decoder is trained on **oracle** fingerprints computed from molecules, but at
> inference it is queried with a noisy, **thresholded, collapsed** spectrum-induced
> fingerprint. Training-inference mismatch.

**Fix:** condition a molecule-language model on fingerprints *and* formulas, then convert
the noisy posterior into a **band of fingerprint queries near the oracle manifold** via
**active-bit density calibration**. Sample candidates across the band, then pool them and
rank by **generation-frequency consensus** (how often the model regenerates the same
structure — a cheap self-consistency signal). A **LoRA adapter** handles domain-specific
posterior bias without disturbing the pretrained molecular prior.

**Results:**

| Benchmark | Top-1 | Top-10 |
|---|---|---|
| NPLIB1 | **29.8%** | **41.1%** |
| MassSpecGym | **23.9%** | **28.7%** |

### 8.5 MSFlow — flow matching (Feb 2026)

- **Encoder:** MIST + **SIRIUS-predicted formula**, aggregating peaks with a **set
  transformer** into a fixed-length **CDDD** embedding (Winter et al.'s continuous
  autoencoder-derived molecular descriptors — more expressive than binary fingerprints).
- **Decoder:** **discrete flow matching**, a BERT-like transformer with **bidirectional**
  attention generating **SAFE** sequences (a better molecular string than SMILES, because it
  encodes substructures directly rather than non-sequentially). Trained with
  **classifier-free guidance** (drop the condition randomly so the model learns both
  conditional and unconditional generation).

**Results:**

| Benchmark | MSFlow Top-1 | DiffMS Top-1 | Ratio |
|---|---|---|---|
| NPLIB1 | **45.00%** | 8.34% | 5.4× |
| MassSpecGym | **32.00%** | 2.30% | **13.9×** |

| Benchmark | MSFlow Top-10 | MSFlow Top-1 Tanimoto | MSFlow MCES |
|---|---|---|---|
| NPLIB1 | 58.53% | 0.72 | 3.79 |
| MassSpecGym | 42.53% | 0.69 | 6.38 |

MSFlow's authors also identify the bottleneck honestly: **the spectrum encoding stage is
the primary limitation** — a large gap remains between oracle CDDD descriptors and ones
derived from experimental spectra.

### 8.6 Test-time tuning (Oct 2025)

An orthogonal axis: rather than a better architecture, **adapt a pre-trained model at test
time.**

- **BART-base** encoder–decoder. Spectrum, formula, and SMILES **all tokenised as text** —
  (m/z, intensity) pairs become strings. Zero modality-specific machinery.
- An **MLP head predicts the Morgan fingerprint** alongside the SMILES, with a
  **binary cross-entropy fingerprint-alignment loss** added to the sequence loss. This
  forces the internal representation to be chemically meaningful rather than
  string-matching-shaped.
- Pretrained on **simulated** spectra (NEIMS/RASSP-style), then **test-time tuned** on
  experimental spectra by selecting informative samples.
- **Formula-constrained decoding:** at each step, mask out tokens that would violate the
  given formula. Output is always formula-consistent.

**Result:** surpasses DiffMS by **100%** on NPLIB1 and **20%** on MassSpecGym. On MassSpecGym
— the heterogeneous, high-domain-shift benchmark — **test-time tuning is essential**, while
on NPLIB1 ordinary fine-tuning is already optimal. Gains from simulated pretraining are
explicit: zero-shot generalisation that pure experimental training lacks.

### 8.7 SpecTUS — the counterexample worth reading

SpecTUS is **not LC-MS/MS** (it is GC-EI-MS) so it does not transfer directly, but it is
the most instructive *evaluation* result in the recent literature:

- BART-like, 354M params, 12 encoder + 12 decoder blocks, 16 heads, embedding 1024,
  FF 4096. Pretrained on ~millions of synthetic spectra from custom-trained **NEIMS** and
  **RASSP**, then fine-tuned on NIST 20.
- On 28,267 held-out NIST 20 spectra: **43% single-suggestion perfect reconstruction**;
  **65% with ten suggestions.**
- **On 76% of the test set, its single suggestion has strictly better Morgan Tanimoto
  similarity than hybrid database search.** On 84% with ten suggestions.

**Why this matters:** against a curated, well-covered library like NIST, a de novo model can
*beat database search outright* on a similarity metric. Retrieval is not automatically
dominant — it is dominant **when the library has coverage**, and it is not when it doesn't.
For CASMI's Class 2 and Class 3, the library has no coverage, and this is precisely where
generation should help. It also suggests **using both and re-ranking**, which is what the
Kaggle solutions do.

---

## 9. The leakage controversy you must know about

### 9.1 The finding

A 2026 critical review in *Molecules* (31(5):769) surveyed this field and reached a
sobering conclusion:

> State-of-the-art models achieve only **4.1% top-10 accuracy** on rigorously
> leakage-controlled benchmarks like MassSpecGym. This stands in stark contrast to earlier,
> overly optimistic reports, which we attribute to **pervasive data leakage in naive data
> splits**. Over 95% of compounds in leakage-controlled benchmarks still elude exact
> structural recovery.

Specifically, the review attributes historical ~45% top-k figures (e.g. MSNovelist) to
**leakage in the original GNPS splits** — where training and test share scaffolds and
near-duplicates. GNPS spectra come from many sources with overlapping chemistry.

### 9.2 Why this bites CASMI participants specifically

CASMI has the same disease in a different form:

| Leakage vector | CASMI specifics |
|---|---|
| **Scaffold overlap** | Compounds recur across the 12 training libraries. Structure counts sum to **more than 275,810** for ~275k unique — explicit overlap. |
| **Instrument overlap** | Test is 100% timsTOF; `enveda-180` is also timsTOF with 183k structures. |
| **Near-duplicate spectra** | Exact duplicates were removed, but the *same compound measured independently by several libraries* remains. |
| **Population overlap** | `enveda-np-examples` (250 common natural products) overlaps with other libraries. |

### 9.3 The review's two most useful conclusions

**(a) Formula conditioning helps, reliably.**
Providing the correct formula consistently improves exact- and close-match accuracy for
every architecture tested. Gains are architecture-dependent but always positive.

**(b) Scaffold prediction, not structure generation, is the bottleneck for scaffold-based
methods.**

| Method | Condition | Top-10 |
|---|---|---|
| MADGEN | **oracle** scaffold | 12.4% |
| MADGEN | **predicted** scaffold | **1.54%** — collapses |
| DiffMS | oracle formula | 4.25% |
| DiffMS | predicted formula (MIST-CF) | 4.1% — barely moves |

An **8× collapse** when the scaffold must be predicted, versus a **2% drop** when the
formula is predicted. That contrast should shape your design: **formula is a safe
conditioning variable; scaffold is not.**

It also explains why the Kaggle winners do not use scaffold-conditioned generators — in a
fully predictive setting they would collapse to ~1.5%.

---

## 10. The forward direction: predict the spectrum

### 10.1 The idea

Every channel so far asks "could this candidate produce these peaks?" Fragmentation
models answer **forward**: given a structure, predict its spectrum, then compare.

### 10.2 Two generations of tools

| Tool | Type | Notes |
|---|---|---|
| **CFM-ID** | classic probabilistic | The long-standing standard |
| **MetFrag** / MetFrag-lite | bond dissociation | Fast; the `frag_score` channel |
| **ICEBERG**, **GLACIER** | modern learned | MIT `ms-pred`; the Phases 5–7 forward models |

### 10.3 Why forward models are only *re-scoring* channels, never primary

Because they are **simulators, and simulators have instrument blind spots.** A fragmenter
trained on Orbitrap libraries will not reproduce timsTOF fragmentation, especially in
collision-energy terms. The Phase 7 library gate exists precisely because letting ICEBERG
and GLACIER override real measured similarity was actively harmful.

**The correct pattern:**

```
  real measured evidence  →  primary signal
  simulated evidence      →  tie-breaker, and ONLY where real evidence is silent
```

### 10.4 Cycle consistency — the interesting frontier

FRIGID pushes this into a loop: generate candidate structures → forward-fragment them with
**ICEBERG** → find fragments inconsistent with the observed spectrum → **refine** those
parts → regenerate. This is the same self-consistency logic MS-GPT uses with
generation-frequency consensus, but using an external physical model instead of sampling
frequency.

It is expensive (one forward-fragmentation call per candidate) but it directly attacks the
failure mode that kills de novo models: **generating a structure that is chemically valid
but spectrally impossible.**

---

## 11. The toolbox

### Formula and structure annotation

| Tool | Purpose |
|---|---|
| **SIRIUS** | formula + fingerprint annotation; widely used, several papers build on its formula output |
| **MIST-CF** | modern formula predictor; ~82.3% top-1 on NPLIB1 |
| **BUDDY** | formula assignment |
| **CSI:FingerID** | fingerprint prediction (used by MSNovelist's pipeline) |

### Fragmentation / forward models

| Tool | Purpose |
|---|---|
| **CFM-ID** | classic in-silico fragmentation |
| **MetFrag** | bond dissociation; fast `frag_score` |
| **ICEBERG**, **GLACIER** | modern learned forward models (MIT ms-pred) |
| **NEIMS**, **RASSP** | neural fragmentation simulators for synthetic-spectrum pretraining |

### Pretraining and spectra infrastructure

| Tool / resource | Purpose |
|---|---|
| **DreaMS** | self-supervised spectrum embeddings + the **GeMS** corpus of curated unannotated spectra mined from MassIVE |
| **MassIVE** | billions of unannotated spectra — a natural self-supervised corpus |
| **matchms** | spectra cleaning, curation, similarity in Python |
| **FragHub** | harmonised merge of the major open MS/MS libraries |

### Benchmarks for honest evaluation

| Benchmark | Why |
|---|---|
| **NPLIB1** | natural-product focused, but historically leaky |
| **MassSpecGym** | rigorously leakage-controlled; the honest number |

### Label caveats from the official dataset description

- **GNPS propagated ("suspect") annotations** — structures transferred to unidentified
  spectra via spectral-network similarity — were **excluded** from `train.parquet` because
  the labels are inferred, not measured. If you pull them from GNPS yourself, be aware you
  are adding inferred labels.
- **`precursor_error_ppm`** ships uncleaned and doubles as a label-quality signal. Filter on
  it to drop mislabelled training rows.

---

## 12. A recommended build order

If you are building this from scratch, do it in this order. Each stage is independently
submittable and independently useful.

### Stage 0 — Validation harness (do this first, or you will waste a month)

```
  ✓ Build the enveda-np-examples split (250 structures)
  ✓ Implement all three masks
  ✓ Log: candidate-pool recall@pool, MRR@1, MRR@5, MRR@25, per mask
  ✓ Sanity checks: all three spectra of a molecule give the same neutral mass
```

**Do not skip this.** [Section 7](#7-the-validation-trap--the-most-important-operational-fact)
explains that a random split ranks methods backwards.

### Stage 1 — Mass filter (2 hours, biggest step per unit effort)

```
  ✓ Neutral mass per molecule from precursor + adduct
  ✓ Load train ∪ COCONUT fingerprints  (~712k structures)
  ✓ Filter to within 8.5 ppm
  ✓ Deduplicate on InChIKey14
  ✓ Fill remaining slots with most-abundant same-mass structures
  ✓ SUBMIT THIS — it is a valid, non-trivial baseline
```

### Stage 2 — Library channel (1 day)

```
  ✓ Spectral entropy similarity to every training spectrum of each candidate
  ✓ Aggregate per structure (max over its spectra)
  ✓ Rank by lib_max
  ✓ Weight spectra by base_peak_intensity when aggregating
```

### Stage 3 — Analog channel (2 days, biggest single jump)

```
  ✓ Precompute an index: (m/z, intensity) inverted index over all training spectra
  ✓ For each query, find matching training spectra ALLOWING A CONSTANT MASS OFFSET
  ✓ Use modified cosine (tolerates cross-instrument m/z offsets)
  ✓ Align precursor masses / handle adduct shift before comparing
  ✓ Map retrieved spectra → owning structures → Tanimoto to candidates
  ✓ Restrict one variant to timsTOF; another to ChEBI/LIPID MAPS
```

### Stage 4 — Fragmentation channel (2 days)

```
  ✓ Run MetFrag-lite on each candidate
  ✓ frag_score = explained observed intensity / total observed intensity
  ✓ Also compute a "gap" variant: only credit peaks the forward model cannot score
```

### Stage 5 — Fingerprint model (2–3 days, GPU)

```
  ✓ Spectrum → Morgan/MACCS fingerprint transformer, trained on COCONUT
  ✓ Feature = dot product of candidate fingerprint with predicted logits
  ✓ Keep it as ONE FEATURE, not as the ranker
```

### Stage 6 — Ranker (2 days)

```
  ✓ Assemble 30–160 features per candidate
  ✓ Train LightGBM with lambdarank, grouped by molecule_id
  ✓ Inspect feature importance — drop what does nothing
  ✓ Validate ONLY on the enveda-np-examples split
```

### Stage 7 — The five tricks (each is small, each is worth points)

```
  ✓ Popularity prior:  z(ranker) + μ(log(1+substances)+log(1+pubmed)), μ≈0.25
  ✓ Library gate: suppress forward-model features when lib_max ≥ 0.7
  ✓ Union, don't just rerank: pool top-40 from each engine, then rank the union
  ✓ Remove fixed-slot constraints on outside-pool candidates
  ✓ Multi-spectrum aggregation: MIN across a molecule's spectra, not mean
```

### Stage 8 — Forward rescoring (2 days)

```
  ✓ Union top-40 candidates
  ✓ Forward-fragment each with ICEBERG / GLACIER
  ✓ Re-rank the union with the forward scores as extra features
  ✓ Library-gate it so it cannot override real measurements
```

### Stage 9 — Only then, generation

```
  ✓ Take the top ~50 candidates, run de novo generation, merge into the pool
  ✓ This helps Class 3 only, which is a minority of the test set
  ✓ Do NOT let generated candidates displace high-library-confidence matches
```

---

## 13. Common mistakes

| Mistake | Consequence | Fix |
|---|---|---|
| **Validating on a random train split** | Optimises toward a score that is **negatively** correlated with the leaderboard | Use `ingest_lib == 'enveda-np-examples'` |
| **Trusting `collision_energy_orig`** | Compares eV against NCE — meaningless | Use `collision_energy_ev` |
| **Running the tutorial's de novo model as the main approach** | Loses to retrieval by ~0.15 MRR | Retrieval first |
| **Submitting fewer than 25 guesses** | Leaves free points on the table | Always 25 |
| **Writing stereochemistry** | Wasted effort; scoring ignores it | Flat SMILES |
| **Using only the strongest spectrum** | Throws away corroborating evidence | Aggregate all of a molecule's spectra |
| **Averaging evidence across spectra** | One contradicting spectrum gets outvoted | Take min / low quantile |
| **Letting a simulator outrank a measurement** | Actively harmful (Phase 7) | Library gate |
| **Pinning outside-pool candidates to fixed slots** | Hard ceiling on achievable score | Let them compete |
| **Forgetting the InChIKey14 dedup** | Two spellings of one skeleton waste slots | Dedup before writing the CSV |
| **Ignoring `base_peak_intensity`** | Learns from noise | Downweight low-quality spectra |
| **Assuming 275k unique structures** | Counts sum to >275,810 — compounds recur across libraries | Dedupe by InChIKey14 |
| **Downloading anything at runtime** | Internet is off in submissions | Pre-upload all assets as Kaggle datasets |
| **Ignoring the 9h × notebook constraint** | Submission rejected | Retrieval needs ~20 min; de novo needs hours |

---

## 14. The one-paragraph summary

Nobody on this leaderboard is generating molecules. The competition is won by **searching
~712,000 candidate structures**: filter them by an 8.5 ppm mass window (which also pins
the molecular formula exactly, by pure arithmetic), then score each survivor with **four
independent evidence channels** — spectral similarity against reference spectra, Tanimoto
similarity to structures found by an **adduct-shift-tolerant modified-cosine search** over
the training set (the single biggest idea, worth ~+0.19), the fraction of peak intensity a
bond-dissociation fragmenter can explain, and the agreement between a candidate's
fingerprint and a spectrum→fingerprint transformer's prediction — and fuse those with a
**LightGBM ranker**. On top sit five cheap ideas: a **popularity prior** from PubChem
substance and citation counts (+0.01), a **library gate** that stops simulated spectra
outranking measured ones, **union-then-rerank** instead of rerank-only, removing any fixed
slot that confines a candidate class, and aggregating a molecule's 1–16 spectra with a
**min** rather than a mean. This reaches **0.420** in one documented public pipeline and
**0.47** at the top of the board. Meanwhile, for genuinely novel structures the research
community is running formula-conditioned graph diffusion (DiffMS, 4.1% top-10 on
leakage-controlled MassSpecGym), formula-free block diffusion (MARLIN, 16.9% top-1),
flow matching on continuous molecular descriptors (MSFlow, 32–45% top-1), and
molecule-language models queried over a fingerprint posterior band (MS-GPT, 29.8% top-1) —
all of which remain far behind retrieval, because **most answers are already in a database**,
and because a single MS/MS spectrum is genuinely under-determined for novel skeletons. The
one thing every serious participant agrees on: **validate on `enveda-np-examples` or not at
all**, because a random holdout from the training set predicts the leaderboard backwards.

---

**Where to go next:**

- [Part 1 — What the Problem Is](01-what-is-the-problem.md)
- [Part 2 — Solving It By Hand](02-solving-it-by-hand.md)
- `research/casmi-denovo-tutorial-notebook.ipynb` — the official de novo baseline. Read it
  to understand the architecture; do not submit it as your main approach.
