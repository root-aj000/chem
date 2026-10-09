# Part 1 — What the Problem Actually Is

> **Who this is for:** someone who has never heard the words "mass spectrometry" and
> thinks chemistry is just memorising the periodic table. That is fine. We start from atoms.
>
> **How to read this:** read top to bottom. Every term is defined the first time it appears.
> If a section says "you already know X" — you do, we covered it above. There is nothing to
> look up.

---

## Table of contents

1. [The story in one paragraph](#0-the-whole-thing-in-one-paragraph)
2. [Molecules: what are we even identifying?](#1-molecules-what-are-we-even-identifying)
3. [Mass: why every atom has a slightly "wrong" weight](#2-mass-why-every-atom-has-a-slightly-wrong-weight)
4. [The mass spectrometer: a machine that weighs and then breaks things](#3-the-mass-spectrometer-a-machine-that-weighs-and-then-breaks-things)
5. [m/z: the number on the x-axis](#4-mz-the-number-on-the-x-axis)
6. [Adducts: the molecule picks up a passenger](#5-adducts-the-molecule-picks-up-a-passenger)
7. [A spectrum is just a list of numbers](#6-a-spectrum-is-just-a-list-of-numbers)
8. [Collisions: why molecules shatter in a controlled way](#7-collisions-why-molecules-shatter-in-a-controlled-way)
9. [SMILES: writing a molecule as a sentence](#8-smiles-writing-a-molecule-as-a-sentence)
10. [InChIKey14: the answer key that ignores the unimportant stuff](#9-inchikey14-the-answer-key-that-ignores-the-unimportant-stuff)
11. [The metric: MRR@25, computed by hand](#10-the-metric-mrr25-computed-by-hand)
12. [The three difficulty classes](#11-the-three-difficulty-classes)
13. [Why this is genuinely hard: the isomer problem](#12-why-this-is-genuinely-hard-the-isomer-problem)
14. [What the computer is handed, precisely](#13-what-the-computer-is-handed-precisely)
15. [A summary you can recite](#14-a-summary-you-can-recite)

---

## 0. The whole thing in one paragraph

Nature makes thousands of little molecules. A lab tool called a **mass spectrometer**
weighs a molecule, breaks it into pieces, and weighs the pieces. The output is a table of
numbers: "a piece weighing 121.07 that is 30% as strong as the biggest piece, a piece
weighing 93.03 that is 12% as strong." Given only that table, your job is to say which
molecule produced it. You may guess up to 25 times, best guess first. You are scored on
how close to the front your correct guess lands.

That is the entire competition. Everything else is detail.

---

## 1. Molecules: what are we even identifying?

A molecule is a small collection of atoms stuck together by "chemical bonds."

### 1.1 Atoms, and the atoms that matter here

Every element has a mass. For this competition you only need six:

| Element | Symbol | Mass (exactly) | Notes |
|---|---|---|---|
| Carbon | C | 12.00000 | The skeleton of all organic molecules |
| Hydrogen | H | 1.00783 | The "padding" — fills empty space on carbon |
| Oxygen | O | 15.99491 | In alcohols, acids, ethers, carbonyls |
| Nitrogen | N | 14.00307 | Rare in natural products, very informative when present |
| Sulfur | S | 31.97207 | In garlic/onion/egg compounds |
| Halogens (F, Cl, Br) | — | 18.998 / 34.969 / 78.918 | Extremely easy to spot in spectra |

**Organic molecule** = a molecule built mostly from carbon, hydrogen, oxygen, nitrogen,
sulfur, halogens. Every natural product — penicillin, morphine, caffeine, quinine — is one.

### 1.2 The shape of a molecule is information

Carbon makes 4 bonds. It can be a straight chain, a branch, or a **ring**. Two molecules
can have *exactly the same atoms in exactly the same numbers* but be arranged differently —
think of a necklace: same beads, different order. Those are called **isomers**.

This matters enormously and we return to it in [Section 12](#12-why-this-is-genuinely-hard-the-isomer-problem).

### 1.3 Two names for the same molecule

A molecule can be described two ways:

- **Structure / diagram** — a picture. Clear to a human, impossible to store reliably.
- **SMILES** — a text code. This is what the competition asks you to output.

---

## 2. Mass: why every atom has a slightly "wrong" weight

This section is the single most important one in the document, because it is where all
of the mathematical power comes from.

### 2.1 The "wrong weights" are the whole trick

Here is the thing that surprises people: **atoms do not weigh round numbers.**

If you weigh one carbon atom on a perfect scale, it is *exactly* 12.000000. Boring.

Weigh one hydrogen atom: it is **1.007825**. Not 1. It is 0.007825 *heavier* than 1.

Weigh one oxygen atom: it is **15.994915**. Not 16. It is 0.005085 *lighter* than 16.

Weigh one nitrogen atom: **14.003074**. Slightly heavier than 14.

This is not a measuring error. This is real, and it is called the **mass defect** — the tiny
discrepancy between an atom's actual mass and its whole-number "nominal" mass.

### 2.2 The defect table (memorise this; it is 5 numbers)

| Element | Exact mass | **Defect** (exact − rounded) | Sign |
|---|---|---|---|
| H | 1.007825 | **+0.007825** | positive |
| C | 12.000000 | **0.000000** | zero |
| N | 14.003074 | **+0.003074** | positive |
| O | 15.994915 | **−0.005085** | negative |
| S | 31.972071 | **−0.027929** | negative |
| P | 30.973762 | −0.026238 | negative |
| F | 18.998403 | −0.001597 | slightly negative |
| Cl | 34.968853 | −0.031147 | negative |
| Br | 78.918338 | −0.081662 | negative |

**The three facts to burn into your memory:**

1. **Carbon is the reference.** It weighs exactly its whole number. Everything else has a defect.
2. **Hydrogen's defect is the biggest positive one** (+0.0078), and molecules are full of hydrogen.
3. **Oxygen's defect is negative** (−0.0051). So oxygen *cancels out* some of hydrogen's contribution.

### 2.3 Why this is useful: defects add up

Here is the miracle. The defect of a **molecule** is just the sum of the defects of its atoms.

**Worked example — caffeine, C₈H₁₀N₄O₂:**

```
8 carbons  × (+0.000000) =  0.000000
10 hydrogens ×(+0.007825) = +0.078250
 4 nitrogens ×(+0.003074) = +0.012296
 2 oxygens  ×(−0.005085) = −0.010170
                              --------
  total defect               = +0.080376
```

Now compute caffeine's actual mass: 8×12 + 10×1.007825 + 4×14.003074 + 2×15.994915 = **194.080376**.

Check: 194.080376 − 194 (the whole-number part) = **0.080376**. ✓ Exact match.

**You just predicted a molecule's exact mass to six decimal places using nothing but a
6-number table and one multiplication.** No machine required. This is exactly the
calculation a chemist does by hand.

### 2.4 Nominal mass and exact mass, side by side

For any measured mass we split it in two:

```
  194.080376
  ──┬───────
    │        └── DEFECT  = 0.080376   ("the interesting bit")
  NOMINAL MASS = 194                  ("the round bit")
```

- **Nominal mass** = the mass rounded down to a whole number. Useful for a quick rule
  (the nitrogen rule, Section 2.5).
- **Exact mass** = the full number. This is what the mass spectrometer reports, and what
  you must match.

### 2.5 The nitrogen rule (a free shortcut)

Look at the elements again. C, H, O, S, halogens — do you notice something?

| Element | Nominal mass | Odd or even? |
|---|---|---|
| Carbon | 12 | even |
| Hydrogen | 1 | **odd** |
| Oxygen | 16 | even |
| Sulfur | 32 | even |
| Chlorine | 35 | **odd** |

Only **nitrogen (14 — even!)**, hydrogen, and the halogens complicate this. Let me state
the practical rule, which holds for normal organic molecules:

> **Odd nominal mass ⟺ an odd number of nitrogens. Even nominal mass ⟺ an even number of nitrogens.**

(Strictly, the rule assumes an even number of hydrogens relative to halogens — "Rule of
13s" territory. For ordinary natural products, the simplified version above is safe and is
what practitioners use.)

**Example:** caffeine has nominal mass 194 (even) and 4 nitrogens (even). ✓
Morphine, C₁₇H₁₉NO₃, nominal 285 (odd), 1 nitrogen (odd). ✓

**Why this helps:** it halves your formula search space for free.

### 2.6 The hydrogen cap

Carbon forms 4 bonds. If you have `c` carbons and `n` nitrogens, the most hydrogens you can
possibly attach without breaking anything is:

```
  max hydrogens  =  2c + 2 + n
```

This is called the **degree of unsaturation / ring-and-double-bond** constraint. It lets
you reject most bad formula guesses instantly, with no arithmetic beyond doubling.

Formula validity check: `h ≤ 2c + 2 + n` and `h ≥ 0`. Most wrong formulas fail this.

---

## 3. The mass spectrometer: a machine that weighs and then breaks things

A mass spectrometer is a tube with a vacuum inside and magnets at the ends. Two stages
matter to you.

### 3.1 Stage A — weigh the intact molecule

1. Dissolve the sample in liquid and spray it into the machine. (Electrospray ionisation.)
2. The spray is put under huge voltage, so molecules gain or lose charge — most commonly
   they **pick up a proton** (H⁺). See Section 5.
3. The charged molecules fly down the tube. A magnet bends their path.
4. **Bending harder = heavier.** The machine measures exactly how hard the ion is pushed,
   and reports back a number: mass divided by charge, called **m/z**.
5. In this instrument (a Bruker **timsTOF**, the same one used for the test data) this is
   extremely precise — accurate to about **5 parts per million**, i.e. 0.000005 relative.

> **5 ppm on a 350 Da molecule = 0.00175 Da.** That is the resolution of the entire puzzle.
> We are not approximating; we know the mass with absurd precision.

### 3.2 Stage B — break it and weigh the pieces

1. The instrument selects **one** ion (the **precursor ion**) and isolates it — everything
   else is thrown away.
2. It fills a collision cell with inert gas (usually nitrogen or helium) and rams the
   precursor into the gas at high speed (**collision energy**, measured in electron volts,
   eV).
3. The molecule **shatters**. Bonds break. The pieces fly out.
4. The machine weighs each piece and records how many there were.

This whole process is called **tandem mass spectrometry**, abbreviated **MS/MS** or
**MS²**. "Tandem" just means "two stages in a row."

### 3.3 Why the shattering is useful — and why it is not random

Here is the key insight that separates this problem from noise:

**Molecules break in *predictable* places.** Specifically, bonds that hold on to a small
side-group tend to break. So the set of pieces you see is a fingerprint of which side
groups existed.

- A molecule with an **–OH group** (like water, or an alcohol) often loses a piece of mass
  **18.01** (a water molecule).
- A molecule with a **carboxylic acid** (–COOH, like in vinegar) often loses **43.99**
  (carbon dioxide, CO₂).
- A molecule with an **aromatic ring** (like benzene) often keeps a piece of mass **91.05**
  (a "tropylium" ion), which is a recognisable ring signature.

Once you learn the list of common "lost pieces," you can work backwards from the observed
pieces to guess what was lost — and therefore what the molecule contained.

That backwards-guessing procedure is [Part 2 of this series](02-solving-it-by-hand.md).

---

## 4. m/z: the number on the x-axis

The x-axis of a mass spectrum is **m/z**: mass divided by charge.

For the molecules in this competition, the charge is almost always **1** (written +1 or −1),
because these are small molecules. So:

> **When the charge is 1, m/z numerically equals the mass of the ion.**

So we can stop saying "m/z" and just say "mass," as the official competition description
does. Good — one less thing.

(The exception is if a molecule carries 2 or 3 charges. Then you divide. The dataset says
this is rare for small molecules, so we assume charge 1 throughout.)

---

## 5. Adducts: the molecule picks up a passenger

A neutral molecule has no charge, so the machine cannot see it. It must be **charged**
before it can fly. The charge usually comes from a passenger atom-group attaching.

Each passenger has a known mass. Write it down — you will need all of them.

### 5.1 Positive-mode adducts (molecule gains a positive passenger)

| Adduct notation | Meaning | Passenger mass |
|---|---|---|
| `[M+H]⁺` | molecule + 1 proton | **1.00728** |
| `[M+NH₄]⁺` | molecule + ammonium | **18.03382** |
| `[M+H₂O+H]⁺` | molecule + water + proton | **19.01784** |
| `[M−2H₂O+H]⁺` | molecule − 2 water + proton | **37.02839** |
| `[M+Na]⁺` | molecule + sodium | **22.98922** |
| `[M+K]⁺` | molecule + potassium | **38.96316** |

### 5.2 Negative-mode adducts (molecule gains a negative passenger)

| Adduct notation | Meaning | Passenger mass |
|---|---|---|
| `[M−H]⁻` | molecule lost a proton | **1.00728** |
| `[M−H₂O−H]⁻` | molecule − water − proton | **19.01784** |
| `[M+CH₂O₂−H]⁻` | molecule + formate anion | **44.99820** |
| `[M+Cl]⁻` | molecule + chloride | **34.96885** |

(`M` always means the **neutral** molecule — the thing with no charge. That is what we are
trying to identify.)

### 5.3 The rule for getting back to the neutral mass

Read the signs:

- **Positive adducts ADD mass to the molecule.** To undo, **subtract**.
- **Negative adducts REMOVE mass from the molecule.** To undo, **add**.

```
  POSITIVE:   M  =  m/z  −  passenger mass
  NEGATIVE:   M  =  m/z  +  passenger mass
```

**Worked example.** The dataset says a peak cluster sits at precursor m/z **301.0354**, and
the adduct column says `[M−H]⁻`. Negative mode, so we ADD:

```
  M  =  301.0354  +  1.00728  =  302.0427
```

That is the neutral molecule's exact mass: **302.0427 Da**. One subtraction. That is
Step 1 of the entire manual procedure.

---

## 6. A spectrum is just a list of numbers

Dropping the graph picture, a mass spectrum in this dataset is literally this:

```
ms2_mzs:                 [301.0354, 271.0248, 257.0456, 229.0143, 201.0194, ...]
ms2_normalized_intensities: [1.00,    0.62,     0.31,     0.44,     0.18,     ...]
```

Row 1 = the mass of each fragment. Row 2 = how strong each fragment was. Aligned: the
*k*-th m/z pairs with the *k*-th intensity.

### 6.1 Normalised intensity: why the biggest peak is 1.0

Raw spectra have wildly varying scales depending on how much sample was injected. So the
competition divides every intensity in a spectrum by the largest intensity in that same
spectrum. Result:

> **The strongest peak in every spectrum has intensity exactly 1.0.** Every other peak is
> a fraction of it: 0.62 means "62% as strong as the strongest peak."

The strongest peak is called the **base peak**.

### 6.2 Base peak intensity is a quality signal

There is a *second*, raw intensity column: `base_peak_intensity` — the un-normalised count
of the base peak. This is a **confidence score for the spectrum itself.** The dataset notes
spectra whose base peak was under 1,000 raw counts were dropped as noise. So:

- base_peak_intensity = 50,000 → strong, trustworthy spectrum
- base_peak_intensity = 900 → noisy, mostly garbage

This is one of the few free signals in the whole dataset. Use it.

### 6.3 Typical size

The competition says the median spectrum has about **20 peaks**. Some have thousands.

**Median 20 peaks is a sobering fact.** Hold onto it. Twenty numbers, and you must name the
exact molecule. Many different molecules will be consistent with those twenty numbers.

---

## 7. Collisions: why molecules shatter in a controlled way

### 7.1 Collision energy

The `collision_energy_ev` column says how hard the precursor was rammed into the gas.

- **Low energy (~20 eV):** the molecule barely breaks. Mostly the intact ion plus a few
  small pieces.
- **High energy (~80 eV):** the molecule explodes. Many small fragments, and often the
  intact skeleton is gone.

The dataset says the instrument sometimes merges several acquisitions at different energies
into one row, so a value can be a list like `[20, 40, 60, 80]`. More energy settings =
more information.

### 7.2 The training-set trap (read this, it will bite you)

The training set was collected from **12 different source libraries** with **dozens of
different instruments**. They reported collision energy differently:

| Unit | Meaning | Example |
|---|---|---|
| `eV` | absolute electron volts | `40` |
| `NCE` | "normalised collision energy" — a Thermo *percentage* scale, 0–100+ | `35` |
| `V` | raw accelerating voltage | `6V` |
| `unknown` | not recorded | `null` |

**`NCE` and `eV` are different quantities on different scales.** A training row with
`collision_energy_orig = 35` might be 35 eV (beam instrument) or 35 NCE (Orbitrap) — those
are not comparable. The competition provides the column `collision_energy_ev` as a best-
effort conversion to a common scale.

**Use `collision_energy_ev`, not `collision_energy_orig`.** The dataset says so explicitly.

---

## 8. SMILES: writing a molecule as a sentence

SMILES is a text code for a structure. You need to read at least simple ones.

### 8.1 The atoms

Any symbol that is a valid element symbol is an atom: `C` carbon, `N` nitrogen, `O`
oxygen, `S` sulfur, `Cl` chlorine, `Br` bromine.

**Lowercase means "in a ring with aromatic (delocalised) bonding":** `c`, `n`, `o`, `s`.

### 8.2 Reading the simplest possible example

```
CCO
```

Read left to right. Each letter is an atom. Consecutive letters are bonded by a single bond.

```
  C  —  C  —  O
  |     |     |
 CH3  CH2   OH
```

That is **ethanol**. The second carbon has 2 hydrogens (2 bonds used), the last oxygen has
1 hydrogen (1 bond used), the first carbon has 3 hydrogens.

### 8.3 The punctuation

| Symbol | Meaning | Example |
|---|---|---|
| `(` `)` | branch — the next part hangs off the atom just before the `(` | `CC(C)O` |
| `=` | double bond | `CC=O` (acetaldehyde) |
| `#` | triple bond | `C#N` (hydrogen cyanide) |
| `1` `2` `3`… | ring-closure labels | `C1CCCCC1` (cyclohexane) |
| `c1ccccc1` | a benzene ring | — |
| `[` `]` | bracket an unusual atom, e.g. charges | `[NH4+]` |

More examples to build intuition:

| SMILES | Structure |
|---|---|
| `CCO` | ethanol (2 carbons, 1 oxygen) |
| `CC(=O)O` | acetic acid (vinegar) |
| `c1ccccc1` | benzene (6-carbon aromatic ring) |
| `CC(=O)Oc1ccccc1C(=O)O` | aspirin |
| `CN1C=NC2=C1C(=O)N(C)C(=O)N2C` | caffeine (the example from the official docs) |

### 8.4 The starter example, unrolled

Take that last one — caffeine — and count. Caffeine is C₈H₁₀N₄O₂, exact mass 194.080376.
You should be able to find 8 `C`, 4 `N`, 2 `O` in that string. You have just read a
molecule's formula without being told it.

### 8.5 What you must output

A SMILES string. It should be **valid** — a computer must be able to parse it. It does not
need to be canonical, or pretty, or stereochemically complete (see below).

---

## 9. InChIKey14: the answer key that ignores the unimportant stuff

This is the crucial scoring detail, and getting it right removes a huge amount of work.

### 9.1 The problem

The same molecule can be written many ways. Water is `O`, `OWO`… well, `O` and `[OH2]` are
the same molecule. A benzene ring can start at any of its six carbons. So a naive
string-comparison scoring would mark correct answers wrong for silly reasons.

### 9.2 The fix: canonicalisation, then a hash

The scorer does this to **both** your answer and the real answer:

1. Parse the SMILES.
2. Apply **RDKit's tautomer canonicalisation** — pick one standard member of the family
   of structures that differ only by where a hydrogen sits on a double bond.
3. Compute the **InChIKey** — a fixed-length text fingerprint of the molecule.
4. Take only the **first block** — the first 14 characters before the first hyphen.

Example from the official docs, glucose:

```
  WQZGKKKJIJFFOK-VFUOTHSA-NSA  ←  full InChIKey
  └────┬────┘
   InChIKey14 ← THIS is what gets compared
```

### 9.3 What each block means

| Block | Encodes | Matters for scoring? |
|---|---|---|
| 1st (14 chars) | **Connectivity / skeleton** — which atoms are bonded how | **YES — this is the score** |
| 2nd | Stereochemistry, isotopes | no |
| 3rd | Version, protonation | no |

### 9.4 The enormous practical consequence

**Stereo and tautomer are free.** So this pair of predictions scores *identically* against a
glucose answer:

```
  OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O    ← stereochemistry written out
  OCC1OC(O)C(O)C(O)C1O                       ← stereochemistry omitted
```

Both reduce to `WQZGKKKJIJFFOK`.

### 9.5 What this means for how you work

| Feature | Do you need to get it right? |
|---|---|
| Which atoms are bonded to which | **YES — this is the whole game** |
| Ring sizes and ring membership | **YES** |
| Which functional groups exist | **YES** |
| Chirality / `@` / `@@` wedges | **NO** |
| Chiral cis/trans, E/Z double bonds | **NO** |
| Protonation state, salts | **NO** |
| Which tautomer you wrote | **NO** |

**You can output flat, stereo-free SMILES.** That is a real simplification: you do not
have to solve stereochemistry, a notoriously error-prone step. It also means the training
labels can be stripped of chirality, which makes learning easier.

---

## 10. The metric: MRR@25, computed by hand

### 10.1 The definition

> **MRR@25 = Mean Reciprocal Rank at 25.**
>
> For each molecule: look down your ranked list of up to 25 guesses. Find the **first**
> position where you guessed correctly. Take **1 ÷ that position**. If you never guess
> correctly, that molecule scores **0**.
>
> Average that number across all test molecules.

### 10.2 The score for each possible position

| Your correct guess is at position… | Reciprocal = 1/position | Score |
|---|---|---|
| 1 | 1/1 | **1.00** |
| 2 | 1/2 | **0.50** |
| 3 | 1/3 | 0.333 |
| 5 | 1/5 | 0.20 |
| 10 | 1/10 | 0.10 |
| 25 | 1/25 | **0.04** |
| not in your list at all | — | **0.00** |

### 10.3 Two worked calculations

**Molecule A.** Your top guess is right.
```
  position = 1  →  1/1  =  1.000
```

**Molecule B.** Your first four guesses are wrong, the fifth is right.
```
  position = 5  →  1/5  =  0.200
```

**Molecule C.** All 25 wrong.
```
  →  0.000
```

**Your final score for these three molecules:**
```
  (1.000 + 0.200 + 0.000) / 3  =  1.200 / 3  =  0.400
```

### 10.4 The four things this metric tells you

**(a) Being right at rank 1 is worth 25× being right at rank 25.**
There is a huge payoff for the very top of your list and almost nothing for the bottom.
Ranking quality, not just recall, is where the points are.

**(b) There is NO penalty for a wrong guess.**
The docs say this explicitly: "there is no penalty for a wrong guess beyond the rank it
occupies." So filling all 25 slots with junk costs you only the chance that junk displaced
a good answer. **Always submit 25 guesses.** Never submit 5 to "be safe."

**(c) Only the FIRST correct guess counts.**
Each molecule has exactly one right answer. Guessing the right structure twice, at ranks 3
and 7, scores 1/3 — no bonus.

**(d) Getting many molecules into the top 25 matters more than getting a few into the top 1.**
If you can put 25 plausible structures in front of every molecule, your floor is ~0.04
even with zero discrimination. Real top solutions score ~0.47, which means they are
genuinely identifying a large fraction of molecules outright.

### 10.5 A realistic spread

If you identified **no** molecules and just submitted 25 random database structures each
time, you'd score maybe 0.01. If you perfectly ranked everything first-try, 1.00.
Current competition leaders are around **0.47** — i.e. they get roughly half the credit
available from a perfect rank-1 guess, averaged over all molecules.

---

## 11. The three difficulty classes

The organisers split test molecules into three classes. **Which class each molecule belongs
to is hidden from you.**

| Class | Definition | How you'd have to solve it |
|---|---|---|
| **1** | The structure has publicly available reference spectra | **Search** a library for the closest spectrum. Near-deterministic if the library has it. |
| **2** | The structure is known (it's in PubChem or the natural-product database COCONUT) but **no public spectra exist** | **Retrieve** candidate structures by mass, then *predict* their spectra in software and compare. |
| **3** | The structure is not in PubChem at all | **Generate from nothing** — invent a structure consistent with the spectrum. |

In increasing difficulty: class 1 → class 2 → class 3.

### 11.1 Why this distinction dominates everything

**The correct algorithm is different for each class.** A de novo generator (which Part 2
and the tutorial notebook build) is the right tool for class 3 and nearly useless for
class 1. A library search is the reverse.

Because you don't know the class in advance, a strong solution runs **all approaches and
merges the results**. That merge, not any single model, is what scores well.

### 11.2 The design intent

The organisers say the class distribution mimics "a real complex mixture," where a small
handful of known compounds dominate and most of what you see is unknown. So expect
**class 1 to be the largest slice** and class 3 the smallest.

**Practical consequence: retrieval beats generation.** If most answers are already in a
database, then searching databases well is most of the score. This is not a hunch — it is
what the leaderboard shows. See [Part 3](03-how-people-actually-solve-it.md).

---

## 12. Why this is genuinely hard: the isomer problem

This section is the crux. If you understand one thing about this competition, understand this.

### 12.1 A concrete demonstration

Four molecules. Here they are:

```
  methyl benzoate      phenyl acetate       benzyl formate      2-hydroxyacetophenone
  C6H5-C(=O)-O-CH3     C6H5-O-C(=O)-CH3     C6H5-CH2-O-CHO    C6H5-C(=O)-OH, OH ortho
```

### 12.2 Step 1: they all have the same exact mass

```
  8 carbons × 12.000000  =  96.000000
  8 hydrogens × 1.007825 =   8.062600
  2 oxygens  × 15.994915 =  31.989830
                        ----------
  C8H8O2 exact mass      = 136.052430 Da
```

**All four: 136.0524.** The mass filter cannot separate them. Not by a single unit.

### 12.3 Step 2: they all produce a benzoyl ion at m/z 105

```
  136.0524  −  105.0340  =  31.0184  =  CH3O (methoxy radical)  → for methyl benzoate
  136.0524  −  105.0340  =  31.0184  =  CH3 (methyl radical)    → for phenyl acetate
```

Same observed peak, same arithmetic, two completely different molecules. And a library
search keyed on that peak ranks all of them essentially tied.

### 12.4 Step 3: what *does* separate them

Only the **low-mass region** and the **relative intensities**:

| Compound | Characteristic low-mass ions |
|---|---|
| methyl benzoate | 77 (C₆H₅⁺), 51, 39 |
| phenyl acetate | 94, 66 |
| benzyl formate | **91** (tropylium, C₇H₇⁺), 79, 65 |
| 2-hydroxyacetophenone | **121** (loss of CH₃•), 93, 65 |

These *are* distinguishable — but only with fine attention, and **only if the spectrum is
clean enough to resolve ions below m/z 40**. Noisy spectra lose exactly the ions that
carry the discriminating information.

### 12.5 The punchline

> **A single MS/MS spectrum is mathematically under-determined.** Many distinct molecules
> produce nearly the same peak list. The data does not always *contain* the answer.

This is why:

1. The metric is **MRR@25, not accuracy**. The organisers know you cannot always tell them
   apart, so they reward you for putting the plausible ones in the right order.
2. The test set gives **multiple spectra per molecule** (median 3). More independent
   fragmentations = more chances to catch a discriminating ion. If a molecule's three
   spectra are all consistent with methyl benzoate but only one rules out benzyl formate,
   you can still win.
3. De novo generation is genuinely unsolved at the exact-skeleton level. See Part 3.

### 12.6 The compound of the two

The manual procedure in [Part 2](02-solving-it-by-hand.md) narrows candidates down, but it
**cannot always finish the job**. What it produces is a *ranked shortlist*. Submitting that
shortlist is a legitimate, respectable solution — and, per Section 10.4(d), a surprisingly
hard one to beat on recall.

---

## 13. What the computer is handed, precisely

Here is the complete data specification. Every column, in plain language.

### 13.1 Test set — `test.parquet`, ~1,500 spectra / ~400 molecules

| Column | What it contains |
|---|---|
| `molecule_id` | Anonymous label (`m_0014ef`). **You predict per molecule_id.** |
| `spectrum_id` | Unique per spectrum. Several spectra share one `molecule_id`. |
| `ms2_mzs` | Array of fragment m/z values. |
| `ms2_normalized_intensities` | Aligned intensities, largest = 1.0. |
| `base_peak_intensity` | Raw (un-normalised) base peak count. **A quality score.** |
| `adduct` | One of the ten in Section 5. Tells you how to get the neutral mass. |
| `ionization_mode` | `positive` or `negative`. |
| `instrument_type` | **Always `timsTOF`** for test. |
| `precursor_mz` | Measured precursor m/z. Combined with `adduct` → neutral mass. |
| `collision_energy_ev` | eV, possibly a list like `[20, 40, 60, 80]`. |
| `collision_energy_orig` | Raw as recorded by the source. Lossless but inconsistent. |
| `collision_energy_orig_units` | `eV`, `NCE`, `V`, or `unknown`. |

### 13.2 Training set — `train.parquet`, ~2.5M spectra / ~275k structures

All 12 test columns **plus** these 6, which are the labels and the useful extras:

| Column | What it contains |
|---|---|
| `normalized_smiles` | **The answer.** RDKit-standardised SMILES. |
| `inchikey` | Full InChIKey of the structure. |
| `inchikey14` | **First block only — the scoring key.** |
| `molecular_formula` | e.g. `C15H10O7`. |
| `ingest_lib` | Which of the 12 source libraries. **Critically important — see 13.4.** |
| `adduct_orig` | Adduct as the source recorded it, before standardisation. |
| `precursor_error_ppm` | ppm error between measured precursor and the true structure's expected mass. **A label-quality flag.** |
| `num_peaks` | How many fragment peaks. |

### 13.3 The training set composition

| Library | Spectra | Structures | Character |
|---|---|---|---|
| `enveda-180` | 1,153,785 | 182,941 | **Same timsTOF as test** — but drug-like synthetic screening compounds, wrong chemistry |
| `pluskal_ms2` | 527,581 | 46,821 | Orbitrap, consistent protocol, multiple energies |
| `riken` | 347,171 | 15,892 | **Plant specialised metabolites — on-topic** |
| `gnps` | 220,849 | 45,750 | Largest natural-product library, most heterogeneous |
| `massbank` | 101,727 | 9,180 | Curated consortium, many instruments |
| `mona` | 92,416 | 11,681 | Fiehn Lab, community |
| `spectraverse` | 50,933 | 9,631 | Harmonised aggregation incl. obscure sources |
| `msdial` | 40,765 | 9,127 | Bundled with MS-DIAL software |
| `drug_plus` | 2,545 | 2,539 | ~2,500 drugs, ~1 spectrum each, no CE metadata |
| `enveda-np-examples` | **1,151** | **250** | **Closest to test.** Natural products, same instrument, same pipeline |
| `masaryk` | 652 | 416 | RECETOX chemical standards |

Note the tension: `enveda-180` is 46% of all training spectra and is the *only large*
library on the test instrument — but its chemistry is **synthetic drug-like screening
compounds**, explicitly described as "a different region of chemical space from the test
molecules." Meanwhile the on-topic natural-product libraries are 10× smaller.

### 13.4 The validation trap — the most important operational fact

If you hold out a random 500 structures from training to validate, your validation score
will **actively mislead you**, for two compounding reasons:

1. **Chemistry mismatch.** `enveda-180` teaches drug-like fragmentation; the test set is
   natural products.
2. **Instrument mismatch.** Test is 100% timsTOF. `enveda-180` and `enveda-np-examples` are
   the only timsTOF-labelled libraries; everything else is Orbitrap, QTOF, etc.

The one split that correlates with the real leaderboard is:

```
  ingest_lib == 'enveda-np-examples'
```

and within it, evaluate under **three masks**:

| Mask | Question it answers |
|---|---|
| Everything visible | How well do I do when the molecule *is* in my pool? |
| That structure's own spectra hidden | Can I identify a molecule with no reference spectra? (≈ class 2) |
| That structure removed from the candidate pool entirely | Can I handle a genuinely novel structure? (≈ class 3) |

This is not a guess. Teams that reported this found that a random held-out split predicted
the leaderboard ranking **backwards** — preferring exactly the model that was worst on the
real board. See Part 3.

### 13.5 A note on `precursor_error_ppm`

This is the ppm gap between the measured precursor m/z and the m/z the *labelled* structure
would produce. Because it ships uncleaned, a large value means the label is probably wrong
— the spectrum was assigned the wrong molecule. It is a free way to **drop mislabelled
training rows** and stop your model learning from errors.

---

## 14. A summary you can recite

**The task.** Given 1–16 MS/MS spectra of one molecule (each a list of ~20 fragment masses
with relative strengths, plus an accurate precursor m/z and a known adduct), rank up to 25
candidate structures, best first.

**The scoring.** Mean Reciprocal Rank at 25. Score for a molecule = 1/(rank of first correct
guess), or 0 if absent. Correctness = InChIKey14 match after RDKit tautomer canonicalisation —
so **connectivity only; stereochemistry and tautomer are free**. No penalty for wrong guesses,
so **always submit 25**.

**What you get for free.**
- Neutral mass to ~5 ppm → at 350 Da that is 0.00175 Da, which pins the **molecular formula**
  to 1–5 options.
- Fragment masses → **neutral losses**, which reveal functional groups.
- `base_peak_intensity` → per-spectrum quality.
- 1–16 spectra per molecule → repeated independent evidence.

**What makes it hard.** Isomers. Different molecules with identical formulas can generate
nearly identical spectra, because one spectrum is genuinely under-determined. Hence 25
guesses, and hence the need to aggregate evidence across a molecule's spectra.

**What actually wins.** Not the deepest neural network. A mass filter that narrows ~700,000
candidate structures down to a handful, several independent evidence channels, and a learned
ranker to order them — with a **prior favouring commonly-reported structures** to break
isomer ties. Top score is ~0.47.

**What to build first, in order.**
1. Neutral mass from precursor + adduct. (One subtraction.)
2. Candidate structures from train ∪ COCONUT ∪ PubChem within ~8.5 ppm.
3. Score each candidate against the molecule's spectra.
4. Aggregate the molecule's spectra into one ranked list.
5. Deduplicate on InChIKey14, write top 25.

**Then, and only then,** consider a neural network — and remember that a de novo generator
is aimed at the hardest third of the problem where the published state of the art still only
recovers ~4% of molecules exactly.

---

**Next:** [Part 2 — Solving It By Hand](02-solving-it-by-hand.md) — the full nine-step manual
procedure with every arithmetic step shown, worked twice from raw spectrum to candidate list.
