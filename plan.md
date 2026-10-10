# Forensic-Level Architecture Plan: CASMI 2026 Molecule Identification from Mass Spectra

## Executive Summary

This plan describes a complete system for identifying chemical structures from tandem mass spectrometry (MS/MS) data. The system is organized as **two independent but complementary tracks**:

- **Track A — Library Search**: Find the molecule by matching against known spectral libraries and structure databases. This is the primary scoring engine and targets Class 1 (in public spectral libraries) and Class 2 (known structure, no public spectra) molecules.
- **Track B — De Novo Identification**: Identify unknown molecules from data alone, without any spectral library. This targets Class 3 (novel structures not in any database).

Both tracks share a common front-end (mass determination, formula inference, candidate generation) but diverge in how they score and rank candidates. The final submission merges both tracks' outputs.

The architecture is **algorithmic first, ML-assisted second**. Every scoring channel is grounded in a physical or combinatorial model of how molecules fragment, not in an end-to-end learned black box. This is deliberate: the docs show that pure ML approaches plateau at 0.087–0.145 MRR while algorithmic retrieval reaches 0.42+.

---

## Part 1: The Problem in One Picture

```
Input:  1–16 MS/MS spectra per molecule
        Each spectrum = [(m/z, intensity), ...]  (median 20 peaks)
        + precursor_mz, adduct, ionization_mode, collision_energy_ev

Output: Ranked list of up to 25 SMILES strings per molecule_id
        Scored by MRR@25 (InChIKey14 match after tautomer canonicalization)

Key constraints:
  - Stereochemistry and tautomer form are FREE (InChIKey14 ignores them)
  - No penalty for wrong guesses → always submit 25
  - Mass accuracy ~5 ppm → molecular formula is essentially determined
  - Isomers with identical formulas are the fundamental ambiguity
  - Test set: ~400 molecules, ~1500 spectra, all on Bruker timsTOF
  - Training set: ~2.5M spectra, ~275k structures, 12 heterogeneous libraries
```

---

## Part 2: High-Level Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        SHARED FRONT-END                                │
│                                                                         │
│  ┌──────────────┐   ┌──────────────┐   ┌──────────────────────────┐   │
│  │  Neutral Mass │──▶│   Molecular  │──▶│  Candidate Pool          │   │
│  │  Determination│   │   Formula    │   │  (train ∪ COCONUT ∪      │   │
│  │              │   │  Inference   │   │   PubChem, 8.5 ppm)       │   │
│  └──────────────┘   └──────────────┘   └──────────┬───────────────┘   │
│                                                    │                    │
│                    ┌───────────────────────────────┼──────────────┐     │
│                    │                               │              │     │
│                    ▼                               ▼              │     │
│  ┌─────────────────────────────┐  ┌─────────────────────────────┐ │     │
│  │   TRACK A: LIBRARY SEARCH   │  │  TRACK B: DE NOVO ID        │ │     │
│  │                             │  │                             │ │     │
│  │  Channel A1: Spectral       │  │  Channel B1: Fragmentation  │ │     │
│  │    Entropy Similarity       │  │    Tree (FT) Computation    │ │     │
│  │  Channel A2: Analog         │  │  Channel B2: Subformula     │ │     │
│  │    Propagation (mass-shift  │  │    Graph Annotation         │ │     │
│  │    modified cosine)         │  │  Channel B3: Bond-Dissoc.   │ │     │
│  │  Channel A3: Fragmentation  │  │    Fragmenter (MetFrag-     │ │     │
│  │    Coverage (MetFrag-lite)  │  │    style, algorithmic)     │ │     │
│  │  Channel A4: Fingerprint    │  │  Channel B4: Formula-       │ │     │
│  │    Model (transformer →     │  │    Constrained Structure    │ │     │
│  │    fingerprint logits)      │  │    Enumeration             │ │     │
│  │                             │  │                             │ │     │
│  │  ┌───────────────────────┐  │  │  ┌───────────────────────┐  │ │     │
│  │  │  LightGBM Ranker      │  │  │  │  Combinatorial        │  │ │     │
│  │  │  + Popularity Prior   │  │  │  │  Structure Assembler  │  │ │     │
│  │  │  + Library Gate       │  │  │  │  + MCES Validation    │  │ │     │
│  │  └───────────┬───────────┘  │  │  └───────────┬───────────┘  │ │     │
│  └──────────────┼──────────────┘  └──────────────┼──────────────┘ │     │
│                 │                                │                │     │
│                 └────────────┬───────────────────┘                │     │
│                              ▼                                    │     │
│  ┌─────────────────────────────────────────────────────────────┐ │     │
│  │              FUSION & SUBMISSION ENGINE                      │ │     │
│  │                                                             │ │     │
│  │  1. Union candidates from both tracks                       │ │     │
│  │  2. Re-score with cross-track features                      │ │     │
│  │  3. Deduplicate on InChIKey14                               │ │     │
│  │  4. Fill to 25 with popularity-ranked same-formula          │ │     │
│  │     near-neighbours                                         │ │     │
│  │  5. Write submission.csv                                    │ │     │
│  └─────────────────────────────────────────────────────────────┘ │     │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Part 3: Shared Front-End — How Each Piece Works

### 3.1 Neutral Mass Determination

**What it does**: Converts the measured precursor m/z and adduct type into the neutral molecule's exact mass.

**How it works, step by step**:

1. Look up the adduct string (e.g., `[M-H]-`) in a lookup table of 10 adducts.
2. Each adduct has a known "passenger mass" — the mass of the attached group.
3. Apply the sign rule:
   - Positive-mode adducts (`[M+H]+`, `[M+Na]+`, etc.) ADD mass to the molecule → **subtract** the passenger mass from the measured m/z.
   - Negative-mode adducts (`[M-H]-`, `[M+Cl]-`, etc.) REMOVE mass from the molecule → **add** the passenger mass to the measured m/z.
4. Sanity check: the result must be in the range 150–1159 Da (the test set's mass range).

**Why this matters**: This is the single highest-value computation in the entire pipeline. A 5 ppm mass accuracy at 350 Da means we know the mass to within 0.00175 Da. This precision is what makes the next step — formula inference — essentially deterministic.

**Multi-spectrum consistency check**: All spectra of the same molecule must yield the same neutral mass (within measurement error). If they don't, the disagreeing spectrum is mis-assigned and should be discarded. This is free quality control.

### 3.2 Molecular Formula Inference

**What it does**: Determines the exact molecular formula (e.g., C₁₅H₁₀O₇) from the neutral mass.

**How it works — the mass defect method**:

The key insight is that atoms do not weigh round numbers. Each element has a "mass defect" — the difference between its exact mass and its nominal (rounded) mass:

| Element | Exact Mass | Defect |
|---------|-----------|--------|
| H | 1.007825 | +0.007825 |
| C | 12.000000 | 0.000000 |
| N | 14.003074 | +0.003074 |
| O | 15.994915 | −0.005085 |
| S | 31.972071 | −0.027929 |
| Cl | 34.968853 | −0.031147 |
| Br | 78.918338 | −0.081662 |

The molecule's total defect is the **sum of its atoms' defects**. This is the "fingerprint" that makes formula inference possible.

**The algorithm**:

1. Split the neutral mass into nominal mass (integer part) and defect (fractional part).
2. Apply the **nitrogen rule**: odd nominal mass → odd number of nitrogens; even nominal mass → even number of nitrogens. Try 0 nitrogens first (most common for natural products).
3. Apply the **hydrogen cap**: max H = 2×C + 2 + N. Any formula violating this is impossible.
4. For each candidate oxygen count (0, 1, 2, ...), solve for hydrogen count:
   ```
   nH = (defect + 0.005085 × nO) / 0.007825
   ```
5. Check if nH is a whole number (within measurement tolerance).
6. Compute carbon count: `nC = (nominal − nH − 16×nO − 14×nN) / 12`
7. Check if nC is a whole number.
8. Exactly one combination survives → that's the formula.

**Why this is powerful**: With 5 ppm mass accuracy, the formula is essentially always determined uniquely. The mass alone — before looking at a single fragment — tells us the atom inventory. This is why the docs say "the formula is essentially always determined uniquely."

**Isotope confirmation**: Check for M+2 satellite peaks:
- 3:1 ratio → chlorine present
- 1:1 ratio → bromine present
- ~4.4% → sulfur present
- M+1 at ~1.1% per carbon → carbon count estimate

### 3.3 Candidate Pool Construction

**What it does**: Assembles the set of all known structures that could possibly be the answer.

**How it works**:

1. Load all structures from:
   - Training set (~275k unique structures)
   - COCONUT natural products database (~400k+ structures)
   - PubChem (subset within mass range, ~millions of structures)
2. For each structure, compute its exact monoisotopic mass from its molecular formula.
3. Filter to structures within 8.5 ppm of the target neutral mass.
4. Deduplicate on InChIKey14 (first block of InChIKey — the connectivity skeleton).
5. Typical result: tens to hundreds of candidates per molecule.

**Why 8.5 ppm**: The test instrument (timsTOF) has ~5 ppm accuracy. We use 8.5 ppm to have a safety margin while still being selective enough to keep the candidate pool manageable.

---

## Part 4: Track A — Library Search (Detailed)

Track A is the primary scoring engine. It asks: "Does a reference spectrum or known structure exist that matches this query?"

### 4.1 Channel A1: Spectral Entropy Similarity (Library Channel)

**The question it answers**: "Does a reference spectrum already exist for this molecule or a very close analogue?"

**How it works**:

1. For each candidate structure, look up all its training spectra (from the 12 source libraries).
2. For each training spectrum, compute **spectral entropy similarity** to the query spectrum.
3. Take the **maximum** similarity across all of the candidate's training spectra.

**Spectral entropy similarity** (Böcker & Młotkowski 2012):
- Unlike naive cosine similarity, it compares the *intensity distributions* of two spectra.
- It is robust to mismatched peak counts — two spectra of the same molecule rarely have identical peak counts due to noise, different instruments, different collision energies.
- The entropy transform down-weights the contribution of noise peaks (which have random intensities) and up-weights the contribution of true signal peaks.
- Score ranges from 0 (no similarity) to 1 (identical).

**Why this is the first channel**: If `lib_max ≈ 1.0`, the answer is almost certainly correct — this is Class 1 (the molecule has a public reference spectrum). It should go to rank 1 immediately.

**Aggregation across multiple spectra**: When a molecule has multiple query spectra, take the **minimum** similarity across spectra (not the mean). One spectrum that flatly contradicts a candidate should disqualify it — a mean would let four agreeable spectra outvote the one that knows better.

### 4.2 Channel A2: Analog Propagation (The Single Biggest Idea)

**The question it answers**: "Is there a closely-related molecule with a reference spectrum, even if this exact molecule doesn't have one?"

**Why this matters**: This is worth ~+0.19 MRR — the single largest improvement in the published solutions. Most people miss it.

**The key insight**: Spectra of similar molecules are often offset by a **constant mass difference**. If compound X has a reference spectrum, and the query looks like X's spectrum but shifted by the mass of a methyl group (14.0157 Da) or an oxygen (15.9949 Da), then the query is likely "X plus a methyl" or "X plus an oxygen."

**How it works, step by step**:

1. **Precompute an inverted index** over all training spectra: for each m/z value, store which spectra contain a peak at that m/z.
2. **For each query spectrum**:
   a. Compute **modified cosine similarity** to every training spectrum, but **allow a constant mass offset**.
   b. The modified cosine works by: for each pair of peaks (one from query, one from training), compute the mass difference. If many peak pairs share the same mass difference, that's evidence of an analogue relationship.
   c. Specifically: shift the training spectrum by Δ m/z and compute cosine similarity. Try all possible Δ values (from the peak-pair differences). The best Δ gives the modified cosine score.
3. **Map retrieved spectra to owning structures**: Each training spectrum belongs to a known structure. The retrieved spectra tell us which structures are "analogues" of the query.
4. **Compute Tanimoto similarity** between the query's implied molecular fingerprint and the fingerprints of those analogue structures.
5. **Rank candidates** by this Tanimoto similarity.

**Adduct-shifted matching**: Before comparing, align the precursor masses. A query with `[M+H]+` at m/z 300 and a training spectrum with `[M-H]-` at m/z 298 are the same molecule (neutral mass 299) but the raw m/z values differ by 2. Without adduct alignment, every cross-adduct pair would be garbage.

**Two variants that proved valuable**:
- **AFIX variant**: Restricted to timsTOF spectra only (matches the test instrument).
- **BIO variant**: Pulls from ChEBI and LIPID MAPS (biologically relevant databases).
- Blending both variants beats either alone.

### 4.3 Channel A3: Fragmentation Coverage (MetFrag-style)

**The question it answers**: "Can the candidate's bond connectivity explain the observed fragment peaks?"

**How it works**:

1. For each candidate structure, run a **bond-dissociation fragmenter**:
   a. Represent the candidate as a molecular graph (atoms = vertices, bonds = edges).
   b. Systematically break each bond (or pair of bonds for rings) and compute the mass of each resulting fragment.
   c. This is a **breadth-first enumeration**: start with the intact molecule, break all possible bonds to get depth-1 fragments, then break all possible bonds in each depth-1 fragment to get depth-2 fragments, etc.
   d. Prune: only keep fragments whose mass matches an observed peak (within tolerance).
2. **Score**: `frag_score = (sum of intensities of explained peaks) / (total observed intensity)`

**Why this is genuinely independent**: It never looks up any reference spectrum. It asks a purely physical question: "If this molecule fragmented, would it produce the peaks we see?" A molecule can be in the pool with zero training spectra (Channels A1 and A2 both dead) and still be strongly supported by Channel A3.

**Algorithmic details** (from MetFrag, Wolf et al. 2010):
- Linear bonds: cleave → two fragments.
- Ring bonds: must cleave two bonds simultaneously to open the ring.
- Rearrangement rules: a small set of rules account for known rearrangements (e.g., neutral losses of unconnected fragments like H₂O, CO₂).
- Redundancy check: eliminate duplicate fragments using graph isomorphism heuristics.
- Scoring: weighted peak count (high-mass, high-intensity peaks weighted more) minus bond dissociation energy penalty (bonds with high BDE are less likely to break).

### 4.4 Channel A4: Fingerprint Model

**The question it answers**: "Does a neural network's predicted molecular fingerprint agree with the candidate's actual fingerprint?"

**How it works**:

1. Train a **spectrum → molecular fingerprint** transformer:
   - Input: MS/MS spectrum (set of (m/z, intensity) peaks).
   - Output: logits over fingerprint bits (e.g., 2048-bit Morgan fingerprint).
   - Training data: COCONUT structures with their known fingerprints (not spectra — this decouples training from the scarcity of labelled spectra).
2. At inference:
   a. Feed the query spectrum to the transformer → predicted fingerprint logits.
   b. For each candidate, compute the dot product of its Morgan fingerprint with the predicted logits.
   c. This gives a similarity score.

**Why this is used as ONE feature, not the ranker**: The neural network is not generating structures — it's producing one more similarity number for the ranker. This is a much easier learning problem than generation, and it's where the neural network earns its keep.

**Key implementation detail**: Training on COCONUT rather than the training set. COCONUT is a natural-products database, so it matches the test chemistry better than Enveda's drug-like data.

### 4.5 The LightGBM Ranker

**What it does**: Combines all four evidence channels into a single ranking score.

**How it works**:

1. **Assemble features** per candidate (30–160 features):
   - A1: `lib_max` (max spectral entropy similarity)
   - A2: `analog_tanimoto` (Tanimoto to best analogue), `analog_modified_cosine`
   - A3: `frag_score` (fragmentation coverage), `frag_weighted_count`, `frag_bde_penalty`
   - A4: `fp_dot_product` (fingerprint model agreement)
   - Metadata: mass difference in ppm, number of spectra, base_peak_intensity
   - Structural: Tanimoto to query's implied fingerprint, molecular weight, formula features
2. **Train LightGBM** with `lambdarank` objective, grouped by `molecule_id`.
3. **Add popularity prior**: `final_score = z(ranker_score) + μ × (log(1 + n_pubchem_substances) + log(1 + n_pubmed_citations))` with μ ≈ 0.25.
4. **Library gate**: When `lib_max ≥ 0.7`, suppress forward-model features. A real measured library match should never be overruled by a simulation.

**Why LightGBM and not a neural net**: After the mass filter, we have ~50 candidates and 30–160 heterogeneous numerical features. This is a tabular learning-to-rank problem — LightGBM's home turf. It handles mixed feature types, works on tiny training sets, needs no normalization, and provides feature importance for free.

---

## Part 5: Track B — De Novo Identification (Detailed)

Track B handles the case where the molecule is NOT in any database — Class 3. This is the genuinely hard case, and the approach must be fundamentally algorithmic.

### 5.1 Channel B1: Fragmentation Tree Computation

**What it is**: A fragmentation tree (FT) is a graph-theoretical model of how a molecule broke apart, computed directly from the spectrum without any database.

**The algorithm** (Böcker & Rasche 2008, SIRIUS):

1. **Build the fragmentation graph**:
   - For each peak in the spectrum, find all molecular formulas that match its mass (within tolerance). These are the "colored" nodes — all formulas for the same peak share the same color.
   - Insert directed edges whenever one formula is a sub-formula of another (i.e., one fragment could be derived from another by losing atoms).
2. **Find the maximum-weight colorful subtree**:
   - A "colorful" subtree uses at most one node per color (each peak is explained at most once).
   - The subtree with the maximum total edge weight is the best fragmentation tree.
   - This is the **Maximum Colorful Subtree** problem — NP-hard in general, but solvable in practice via integer linear programming (ILP) or dynamic programming.
3. **Score the tree** using a Bayesian model:
   - Prior probability of each fragmentation step (learned from training data).
   - Likelihood of the observed intensities given the tree.
   - The best tree is the Maximum A Posteriori (MAP) estimate.

**How this helps identification**: The fragmentation tree tells us:
- The most likely molecular formula (the root of the tree).
- The most likely neutral losses (the edges).
- The fragmentation cascade (which fragments came from which).

This is used to:
- Rank candidate formulas (the root formula of the best FT is the most likely formula).
- Predict a molecular fingerprint (via CSI:FingerID-style ML on the FT).
- Search structure databases (compare predicted fingerprint to database fingerprints).

### 5.2 Channel B2: Subformula Graph Annotation

**What it is**: A newer, simpler approach (Li 2024) that represents the spectrum as a graph of formula-subformula relationships.

**The algorithm**:

1. For each candidate parent formula, construct the **Parent Subformula Graph (PSG)**:
   - Nodes = molecular formulas of fragments.
   - Edges = subformula relationships (fragment A's formula is a sub-formula of fragment B's formula).
2. Score the PSG using **purely topological features** (no empirical weights):
   - Vertex count (more explained peaks = better).
   - Edge density (more subformula relationships = more internally consistent).
   - Connectivity patterns.
3. Rank candidate formulas by their PSG scores.

**Why this is novel and algorithmic**: It requires no training data, no learned parameters, and no fragmentation rules. It purely asks: "Which candidate formula produces the most internally consistent network of formula-subformula relationships?"

### 5.3 Channel B3: Bond-Dissociation Fragmenter (Algorithmic MetFrag)

**What it is**: A purely combinatorial fragmenter that enumerates all possible fragments of a candidate structure by systematically breaking bonds.

**The algorithm** (improved from MetFrag, Wolf et al. 2010):

1. **Input**: Candidate molecular graph (atoms + bonds).
2. **Enumerate fragments**:
   - Depth 1: Break each single bond → pairs of fragments.
   - Depth 2: Break each single bond in each depth-1 fragment → depth-2 fragments.
   - For rings: break two bonds simultaneously.
   - Continue to a maximum depth (typically 2–3).
3. **Match to observed peaks**: For each enumerated fragment, check if its mass matches an observed peak (within 5 ppm).
4. **Score**: Weighted peak count (high-mass, high-intensity peaks weighted more) minus BDE penalty.

**Key algorithmic optimizations**:
- **Hopcroft-Tarjan bridge finding**: Efficiently find all bridges (bonds whose removal disconnects the graph) — these are the most likely cleavage points.
- **2-cut enumeration**: Find all pairs of bonds whose removal disconnects the graph — these model ring openings.
- **Bitwise AND for depth-n fragments**: Depth-n fragments can be computed by intersecting depth-(n-1) fragments with depth-1 fragments, avoiding redundant enumeration.
- **Graph isomorphism pruning**: Eliminate duplicate fragments using canonical labeling.

### 5.4 Channel B4: Formula-Constrained Structure Enumeration

**What it is**: Given a molecular formula, enumerate all chemically valid structures, then score each against the spectrum.

**The algorithm**:

1. **Generate all structures** for the formula using a structure generator (e.g., MAYGEN, SURGE, MOLGEN).
2. **Filter by chemical validity**:
   - Valence constraints (each atom has the right number of bonds).
   - Ring constraints (no impossible ring sizes).
   - Degree of unsaturation (RDBE) must match the formula.
3. **Score each surviving structure** using the bond-dissociation fragmenter (B3).
4. **Rank** by fragmentation score.

**Why this is exponential but tractable**: For a formula like C₁₅H₁₀O₇, there are millions of possible structures. But:
- The mass filter already pinned the formula.
- Chemical validity filtering eliminates >99% of candidates.
- The fragmenter is fast (milliseconds per structure).
- We only need the top ~25, not the full enumeration.

**The MCES validation step** (from Fragment Tree Reconciliation, Domschke et al. 2026):
- For each candidate structure, compute the **Maximum Common Edge Subgraph (MCES)** between the candidate and the fragmentation tree.
- The MCES tells us which bonds are preserved and which are broken.
- Candidates with larger MCES (more bonds preserved) are more consistent with the observed fragmentation.
- This is a graph-theoretic measure that requires no learning.

### 5.5 The Combinatorial Structure Assembler

**What it does**: Assembles candidate structures from the fragmentation evidence, then validates them.

**How it works**:

1. **Extract constraints** from the fragmentation tree:
   - The root formula (molecular formula).
   - The set of neutral losses (edges of the tree).
   - The set of observed fragment masses.
2. **Assemble candidate structures**:
   - Start with the root formula.
   - For each neutral loss, determine which substructure could have been lost.
   - Use a **constraint satisfaction** approach: find all structures that are consistent with all observed losses simultaneously.
3. **Validate with MCES**:
   - For each assembled candidate, compute the MCES with the fragmentation tree.
   - Rank by MCES size (larger = more consistent).
4. **Fill remaining slots** with near-neighbours of the best candidate (same formula, one substituent changed).

---

## Part 6: Fusion & Submission Engine

### 6.1 Union, Don't Just Re-Rank

**The key insight**: Most re-ranking schemes only re-score the existing candidate list. This approach first **unions** the lists from both tracks, letting a candidate proposed by one track be validated by another track's evidence.

**How it works**:

1. Take the top ~40 candidates from Track A (library search).
2. Take the top ~40 candidates from Track B (de novo).
3. Union them (deduplicate on InChIKey14).
4. Re-score every candidate in the union with **all features from both tracks**.
5. This lets a candidate that was rank 15 in Track A (because it has no library match) get boosted by Track B's fragmentation evidence.

### 6.2 Cross-Track Feature Fusion

**How it works**:

For each candidate in the union, compute features from BOTH tracks:
- Track A features: lib_max, analog_tanimoto, frag_score, fp_dot_product
- Track B features: ft_score, psg_score, bond_dissoc_score, mces_score
- Cross-features: does the candidate appear in both tracks? (agreement signal)

Train a final LightGBM ranker on these cross-track features.

### 6.3 Deduplication on InChIKey14

**Why**: The same molecule can be written many ways (different SMILES, different tautomers, different stereochemistry). The scorer uses InChIKey14 (first block of InChIKey) which encodes only connectivity — stereochemistry and tautomer form are ignored.

**How**: Before writing the submission, deduplicate all candidates on their InChIKey14. If two candidates have the same InChIKey14, keep only the higher-ranked one.

### 6.4 Fill to 25

**Why**: There is no penalty for wrong guesses. A correct guess at rank 25 scores 0.04 — better than 0.

**How**: If we have fewer than 25 unique candidates, fill the remaining slots with:
1. Near-neighbours of the best candidate (same formula, one substituent changed or moved).
2. The most database-abundant compounds with the same formula (popularity prior).
3. If still short, compounds with formulas within 1 CH₂ (14.0157 Da) of the target.

### 6.5 Submission Format

```
molecule_id,smiles
m_0014ef,CC1=CC(=O)C=CC1=O;OC(=O)c1ccccc1O;CN1C=NC2=C1C(=O)N(C)C(=O)N2C;...
```

- Up to 25 SMILES per molecule, semicolon-separated, best first.
- No stereochemistry needed (InChIKey14 ignores it).
- Every molecule_id must appear exactly once.

---

## Part 7: Validation Strategy

### 7.1 The Validation Trap

**The most important operational fact**: A random holdout from the training set will **actively mislead** you. It predicts the leaderboard ranking **backwards**.

**Why**: Two compounding mismatches:
1. **Chemistry mismatch**: `enveda-180` is 46% of all training spectra and is synthetic drug-like screening compounds. The test set is natural products.
2. **Instrument mismatch**: Test is 100% timsTOF. Only `enveda-180` and `enveda-np-examples` are timsTOF-labelled; everything else is Orbitrap, QTOF, etc.

### 7.2 The Correct Validation Design

```
Validation set = train_df[train_df.ingest_lib == 'enveda-np-examples']  # 250 structures
```

Evaluate under **three masks**:

| Mask | What it measures | Maps to |
|------|-----------------|---------|
| Everything visible | Best case — molecule is in pool, spectra available | Class 1 |
| Structure's own spectra hidden | Molecule in pool, no reference spectrum | Class 2 |
| Structure removed from candidate pool | Genuinely novel structure | Class 3 |

**Read the mask-2 and mask-3 scores separately**. They tell you what fraction of the test set you are structurally incapable of solving.

### 7.3 Metrics to Track

- `candidate_pool_recall@pool`: Is the correct structure in the candidate pool at all?
- `MRR@1`, `MRR@5`, `MRR@25`: Standard reciprocal rank metrics.
- Per-mask breakdown: Class 1 / Class 2 / Class 3 scores separately.

---

## Part 8: Build Order

Each stage is independently submittable and independently useful.

### Stage 0 — Validation Harness (do this first, or waste a month)

```
✓ Build the enveda-np-examples split (250 structures)
✓ Implement all three masks
✓ Log: candidate-pool recall@pool, MRR@1, MRR@5, MRR@25, per mask
✓ Sanity check: all spectra of a molecule give the same neutral mass
```

### Stage 1 — Mass Filter (2 hours, biggest step per unit effort)

```
✓ Neutral mass per molecule from precursor + adduct
✓ Load train ∪ COCONUT fingerprints (~712k structures)
✓ Filter to within 8.5 ppm
✓ Deduplicate on InChIKey14
✓ Fill remaining slots with most-abundant same-mass structures
✓ SUBMIT THIS — it is a valid, non-trivial baseline
```

### Stage 2 — Library Channel (1 day)

```
✓ Spectral entropy similarity to every training spectrum of each candidate
✓ Aggregate per structure (max over its spectra)
✓ Rank by lib_max
✓ Weight spectra by base_peak_intensity when aggregating
```

### Stage 3 — Analog Channel (2 days, biggest single jump)

```
✓ Precompute an index: (m/z, intensity) inverted index over all training spectra
✓ For each query, find matching training spectra ALLOWING A CONSTANT MASS OFFSET
✓ Use modified cosine (tolerates cross-instrument m/z offsets)
✓ Align precursor masses / handle adduct shift before comparing
✓ Map retrieved spectra → owning structures → Tanimoto to candidates
✓ Restrict one variant to timsTOF; another to ChEBI/LIPID MAPS
```

### Stage 4 — Fragmentation Channel (2 days)

```
✓ Run MetFrag-lite on each candidate
✓ frag_score = explained observed intensity / total observed intensity
✓ Also compute a "gap" variant: only credit peaks the forward model cannot score
```

### Stage 5 — Fingerprint Model (2–3 days, GPU)

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

### Stage 7 — The Five Tricks (each is small, each is worth points)

```
✓ Popularity prior: z(ranker) + μ(log(1+substances)+log(1+pubmed)), μ≈0.25
✓ Library gate: suppress forward-model features when lib_max ≥ 0.7
✓ Union, don't just rerank: pool top-40 from each engine, then rank the union
✓ Remove fixed-slot constraints on outside-pool candidates
✓ Multi-spectrum aggregation: MIN across a molecule's spectra, not mean
```

### Stage 8 — Forward Rescoring (2 days)

```
✓ Union top-40 candidates
✓ Forward-fragment each with ICEBERG / GLACIER
✓ Re-rank the union with the forward scores as extra features
✓ Library-gate it so it cannot override real measurements
```

### Stage 9 — De Novo Track (Track B, for Class 3)

```
✓ Fragmentation tree computation (SIRIUS-style ILP)
✓ Subformula graph annotation
✓ Bond-dissociation fragmenter on de novo candidates
✓ Formula-constrained structure enumeration for top candidates
✓ MCES validation of assembled structures
```

---

## Part 9: How Similar Problems Were Solved (Research Synthesis)

### 9.1 The Fragmentation Tree Problem

**Similar problem**: Given a set of fragment masses, find the most likely fragmentation cascade.

**How it was solved** (Böcker & Rasche 2008, SIRIUS):
- Formulated as Maximum Colorful Subtree on a fragmentation graph.
- NP-hard, but solved via ILP or dynamic programming.
- Bayesian scoring (MAP estimate) with learned parameters.
- Used for molecular formula determination and fingerprint prediction.

**How we use it**: Track B's primary tool. The FT gives us the most likely formula and fragmentation cascade without any database. We use it to:
1. Rank candidate formulas.
2. Predict a molecular fingerprint (via ML on the FT).
3. Constrain structure enumeration.

### 9.2 The Bond-Dissociation Problem

**Similar problem**: Given a candidate structure, predict its fragments.

**How it was solved** (MetFrag, Wolf et al. 2010):
- Combinatorial enumeration of all bond cleavages.
- Breadth-first with pruning (mass check, redundancy check).
- Scoring: weighted peak count minus BDE penalty.
- Rearrangement rules for known non-bond-breaking losses.

**How we use it**: Track A's Channel A3 and Track B's Channel B3. We improve on MetFrag with:
- Hopcroft-Tarjan bridge finding for efficient cleavage point identification.
- Bitwise AND for depth-n fragment computation.
- MCES-based validation against the fragmentation tree.

### 9.3 The Fingerprint Prediction Problem

**Similar problem**: Given a spectrum, predict which substructures are present.

**How it was solved** (CSI:FingerID, Dührkop et al. 2015):
- Compute fragmentation tree from spectrum.
- Use kernel methods (SVM) to predict fingerprint bits from FT + spectrum.
- Multiple kernel learning to combine different similarity measures.
- Search structure database by fingerprint Tanimoto.

**How we use it**: Track A's Channel A4. We use a transformer instead of SVM, but the principle is the same: predict a fingerprint, then search by fingerprint similarity.

### 9.4 The Database Search Problem

**Similar problem**: Given a spectrum, find the best-matching structure in a database.

**How it was solved** (MolDiscovery, Cao et al. 2021):
- Construct fragmentation graph using bridge/2-cut enumeration.
- Learn a probabilistic model P(logRank | bondType) from training data.
- Score candidates by log-likelihood ratio.
- Efficient indexing for searching millions of structures.

**How we use it**: Track A's overall architecture. We adopt the fragmentation graph + probabilistic scoring approach, but replace the learned model with a GBDT ranker that combines multiple evidence channels.

### 9.5 The De Novo Sequencing Problem

**Similar problem**: Given a spectrum of a linear molecule (peptide), reconstruct the sequence.

**How it was solved** (Kreitzberg et al., alphabet projection):
- Infer an "alphabet" of mass differences from the spectrum collection.
- Build a graph where edges connect peaks differing by an alphabet mass.
- Find large connected components → likely substructures.
- Use LSH to find recurring subgraphs across spectra.

**How we use it**: Track B's inspiration. For branched molecules (not just linear), we use the fragmentation tree instead of a simple alphabet, but the principle is the same: find the set of mass differences that best explains the observed peaks.

### 9.6 The Structure Enumeration Problem

**Similar problem**: Given a molecular formula, generate all valid structures.

**How it was solved** (MAYGEN, SURGE, MOLGEN):
- Combinatorial generation with chemical validity filtering.
- Canonical labeling to avoid duplicates.
- Constraint satisfaction to incorporate spectral evidence.

**How we use it**: Track B's Channel B4. We enumerate structures for the top candidate formulas and score each against the spectrum.

---

## Part 10: Novel Algorithmic Contributions

This plan introduces several novel algorithmic ideas that go beyond pure ML:

### 10.1 Cross-Track Feature Fusion with Union Rescoring

**The idea**: Instead of running Track A and Track B independently and taking the best, we:
1. Union the candidate lists.
2. Re-score every candidate with features from BOTH tracks.
3. This lets a candidate with no library match (Track A score = 0) get boosted by strong fragmentation evidence (Track B score = high).

**Why it's novel**: Most systems run retrieval and de novo independently. The union-then-rerank approach lets evidence from one track rescue candidates that the other track would have missed.

### 10.2 MCES-Based Structure Validation

**The idea**: Use Maximum Common Edge Subgraph (MCES) between a candidate structure and the fragmentation tree as a validation score.

**Why it's novel**: MCES is a graph-theoretic measure that requires no learning. It directly quantifies how much of the candidate's bond structure is consistent with the observed fragmentation. This is a purely algorithmic alternative to learned scoring.

### 10.3 Subformula Graph Scoring

**The idea**: Score candidate formulas by the topological properties of their Parent Subformula Graph (PSG) — no learned parameters needed.

**Why it's novel**: This is a training-free scoring method that works for any instrument, any chemistry, any mass range. It provides a baseline that doesn't depend on the training set's biases.

### 10.4 Blind Alphabet Projection for Fragmentation Patterns

**The idea**: Infer an "alphabet" of mass differences from the spectrum itself (not from a predefined list of neutral losses), then use this alphabet to build a fragmentation graph.

**Why it's novel**: Traditional fragmentation methods use a predefined list of common losses (H₂O, CO₂, etc.). The blind approach discovers the relevant losses from the data, making it applicable to novel chemistry where the common losses are unknown.

### 10.5 Multi-Spectrum Evidence Aggregation with MIN

**The idea**: When aggregating evidence across a molecule's multiple spectra, take the MINIMUM score across spectra (not the mean).

**Why it's novel**: This is a simple but important design choice. One spectrum that flatly contradicts a candidate should disqualify it. A mean would let four agreeable spectra outvote the one that knows better. This is a principled approach to evidence aggregation that doesn't require learning.

---

## Part 11: Expected Score Progression

| Stage | What was added | Expected MRR@25 |
|-------|---------------|-----------------|
| 0 | Validation harness only | — |
| 1 | Mass filter + fill to 25 | ~0.05–0.10 |
| 2 | + Library channel | ~0.15–0.20 |
| 3 | + Analog channel | ~0.30–0.35 |
| 4 | + Fragmentation channel | ~0.33–0.38 |
| 5 | + Fingerprint model | ~0.35–0.40 |
| 6 | + LightGBM ranker | ~0.38–0.42 |
| 7 | + Five tricks | ~0.40–0.44 |
| 8 | + Forward rescoring | ~0.42–0.45 |
| 9 | + De novo track (Track B) | ~0.43–0.47 |

The leaderboard top is ~0.47. This plan is designed to reach that level.

---

## Part 12: Risk Analysis

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| Validation trap (random split misleads) | High | Very high | Use enveda-np-examples split exclusively |
| Training set chemistry mismatch | High | High | Weight enveda libraries higher; use COCONUT for fingerprint training |
| Instrument mismatch | High | High | Restrict analog channel to timsTOF; use modified cosine for cross-instrument |
| Candidate pool misses the correct structure | Medium | High | Use 8.5 ppm window; add PubChem join; check recall@pool |
| Overfitting to validation set | Medium | Medium | 250 structures is small; use public leaderboard for fine-tuning |
| De novo track produces garbage | Medium | Medium | Library gate: never let de novo override high-confidence library matches |
| Submission format errors | Low | High | Validate format before submission; test with sample_submission.csv |

---

## Part 13: Summary

This plan describes a complete, algorithmic-first system for molecule identification from mass spectra. The key design principles are:

1. **Mass filter first**: The 5 ppm mass accuracy makes formula inference essentially deterministic. This is the highest-value step.
2. **Retrieval before generation**: Most answers are in databases. Search first, generate only for the remainder.
3. **Multiple independent evidence channels**: Library, analog, fragmentation, and model channels fail independently. Combining four weak signals beats one strong one.
4. **Algorithmic scoring, not black-box ML**: Every channel is grounded in a physical or combinatorial model of fragmentation. ML is used only to combine channels, not to replace them.
5. **Union, don't just rerank**: Pool candidates from both tracks, then re-score with cross-track features.
6. **Always submit 25**: No penalty for wrong guesses. Fill with popularity-ranked near-neighbours.
7. **Validate correctly**: Use enveda-np-examples or not at all. A random split predicts the leaderboard backwards.

The two tracks are:
- **Track A (Library Search)**: Spectral entropy similarity → analog propagation → fragmentation coverage → fingerprint model → LightGBM ranker. Targets Class 1 and Class 2.
- **Track B (De Novo)**: Fragmentation tree → subformula graph → bond-dissociation fragmenter → formula-constrained enumeration → MCES validation. Targets Class 3.

Together, they cover all three difficulty classes and are designed to reach the top of the leaderboard (~0.47 MRR@25).
