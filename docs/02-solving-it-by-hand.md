# Part 2 — Solving It By Hand (Pen and Paper)

> **What you need:** a pencil, paper, and a calculator. That is genuinely all.
>
> **What you do NOT need:** calculus, linear algebra, statistics, machine learning, a
> chemistry degree, or any software.
>
> **How to use this document:** Sections 1–10 are the procedure. Sections 11–14 are four
> complete worked examples where I do every single arithmetic step out loud, including the
> dead ends. Sections 15–17 are the wrap-up. If you read only one thing, read
> [Section 11](#11-worked-example-a--quercetin-c15h10o7) — it is the whole method end to end.

---

## Table of contents

- [0. The promise](#0-the-promise)
- [1. Your toolkit](#1-your-toolkit-the-only-things-you-must-memorise)
- [2. Step 1 — Get the neutral mass](#2-step-1--get-the-neutral-mass)
- [3. Step 2 — Split into nominal mass and defect](#3-step-2--split-into-nominal-mass-and-defect)
- [4. Step 3 — Apply the nitrogen rule](#4-step-3--apply-the-nitrogen-rule)
- [5. Step 4 — Find the molecular formula](#5-step-4--find-the-molecular-formula)
- [6. Step 5 — Read the isotope peaks](#6-step-5--read-the-isotope-peaks)
- [7. Step 6 — Compute neutral losses](#7-step-6--compute-neutral-losses)
- [8. Step 7 — Recognise diagnostic ions](#8-step-7--recognise-diagnostic-ions)
- [9. Step 8 — Assemble candidates and test them](#9-step-8--assemble-candidates-and-test-them)
- [10. Step 9 — Rank and submit](#10-step-9--rank-and-submit)
- [11. Worked Example A — quercetin](#11-worked-example-a--quercetin-c15h10o7)
- [12. Worked Example B — benzoic acid](#12-worked-example-b--benzoic-acid-two-adducts-two-answers)
- [13. Worked Example C — vanillin](#13-worked-example-c--vanillin-finding-the-aldehyde-and-the-methoxy)
- [14. Worked Example D — the isomer trap](#14-worked-example-d--the-isomer-trap-why-you-submit-25)
- [15. Using multiple spectra of the same molecule](#15-using-multiple-spectra-of-the-same-molecule)
- [16. What this method cannot do](#16-what-this-method-cannot-do)
- [17. One-page cheat sheet](#17-one-page-cheat-sheet)

---

## 0. The promise

A chemist with a pencil genuinely can identify a molecule from a mass spectrum. They will
usually be right. It takes practice — maybe 50–100 spectra before you are fast — but the
procedure is deterministic bookkeeping, not inspiration.

Here is the shape of it in one picture:

```
  ┌─────────────────────────────────────────────────────────────┐
  │  STEP 1   Precursor m/z  ──subtract adduct──▶  NEUTRAL MASS │
  └─────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
  ┌─────────────────────────────────────────────────────────────┐
  │  STEP 2-4  Split into 302 + 0.0427                           │
  │            Nitrogen rule, hydrogen cap, tiny search          │
  │            ──▶  MOLECULAR FORMULA   (C15H10O7)              │
  └─────────────────────────────────────────────────────────────┘
                                  │
                    ┌─────────────┴─────────────┐
                    ▼                           ▼
  ┌───────────────────────────────┐  ┌──────────────────────────────┐
  │  STEP 6   Look at the gaps     │  │  STEP 5   Look at the +2     │
  │  between peaks. Each gap is    │  │  satellites. 3:1 = Cl,      │
  │  a LOST PIECE.                 │  │  1:1 = Br, ~4% = S.         │
  │  ──▶ FUNCTIONAL GROUPS         │  │  ──▶ HALOGENS / SULFUR       │
  └───────────────────────────────┘  └──────────────────────────────┘
                    │                           │
                    └─────────────┬─────────────┘
                                  ▼
  ┌─────────────────────────────────────────────────────────────┐
  │  STEP 7-8  Recognise ring and group "fingerprints"         │
  │            Then DRAW candidates and check each one's        │
  │            predicted fragments against the observed peaks.  │
  │            ──▶  RANKED SHORTLIST OF 2–5 STRUCTURES         │
  └─────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
  ┌─────────────────────────────────────────────────────────────┐
  │  STEP 9  Submit the shortlist, plus near-neighbours        │
  │           of the best guess, up to 25 total.               │
  └─────────────────────────────────────────────────────────────┘
```

Everything below is just those five boxes filled in with arithmetic.

---

## 1. Your toolkit: the only things you must memorise

Five tables. That is the entire memorisation burden. Table 1 is the important one.

### TABLE 1 — Atom masses and defects (the master table)

| Element | Symbol | Exact mass | Defect (exact − whole) | Sign |
|---|---|---|---|---|
| Hydrogen | H | 1.007825 | **+0.007825** | + |
| Carbon | C | 12.000000 | **0.000000** | 0 |
| Nitrogen | N | 14.003074 | **+0.003074** | + |
| Oxygen | O | 15.994915 | **−0.005085** | − |
| Sulfur | S | 31.972071 | **−0.027929** | − |
| Phosphorus | P | 30.973762 | −0.026238 | − |
| Fluorine | F | 18.998403 | −0.001597 | − |
| Chlorine | Cl | 34.968853 | **−0.031147** | − |
| Bromine | Br | 78.918338 | **−0.081662** | − |
| Iodine | I | 126.904473 | −0.095527 | − |

**Three sentences to remember:**
1. Carbon is exactly 12. It is the ruler everything else is measured against.
2. Hydrogen's defect is **big and positive** (+0.0078), and molecules are full of hydrogen.
3. Oxygen's defect is **negative** (−0.0051), so it *cancels* hydrogen's contribution.

### TABLE 2 — Adducts: how to get from m/z back to the neutral mass

| Adduct | Mode | Rule | Passenger mass |
|---|---|---|---|
| `[M+H]⁺` | + | subtract | 1.00728 |
| `[M+NH₄]⁺` | + | subtract | 18.03382 |
| `[M+H₂O+H]⁺` | + | subtract | 19.01784 |
| `[M−2H₂O+H]⁺` | + | subtract | 37.02841 |
| `[M+Na]⁺` | + | subtract | 22.98922 |
| `[M+K]⁺` | + | subtract | 38.96316 |
| `[M−H]⁻` | − | **add** | 1.00728 |
| `[M−H₂O−H]⁻` | − | **add** | 19.01784 |
| `[M+CH₂O₂−H]⁻` | − | **add** | 44.99820 |
| `[M+Cl]⁻` | − | **add** | 34.96885 |

**The rule in one breath:** positive-mode adducts *add* mass, so **subtract** them.
Negative-mode adducts *remove* mass, so **add** them back.

### TABLE 3 — Neutral losses: the "lost piece" lookup

When a fragment is lighter than the precursor, the difference is a piece that broke off.
Look the difference up here.

| Loss | Exact mass | Tells you |
|---|---|---|
| H radical | 1.00728 | just a deprotonation — uninformative |
| H₂ | 2.01565 | very common background; usually ignore |
| **CH₃ radical** | **15.02347** | **a methoxy group (–OCH₃), or an N-methyl** |
| H₂O | 18.01056 | an **–OH**, an alcohol, or a **carboxylic acid** |
| HCN | 27.01090 | a nitrogen-containing ring (pyridine-like) |
| **CO** | **27.99491** | **a carbonyl; a ring opening; retro-aldol** |
| **CH₂O** | **30.01056** | **a methoxy on an aromatic ring; formaldehyde** |
| CH₃OH | 32.02621 | a **methyl ester** (–COOCH₃) |
| **CO₂** | **43.98983** | **a carboxylic acid — often the base peak** |
| CH₂CO (ketene) | 42.01056 | an aromatic ester |
| HCOOH | 46.00548 | a formate ester |
| H₂O + CO | 46.00547 | a phenolic ring opening |
| C₂H₄O₂ | 60.02113 | an acetyl group; an ester side chain |
| C₃H₄O₂ | 72.02113 | a ring fragment |
| C₆H₅O₂ | 109.02895 | loss of a dihydroxyphenyl piece |
| C₇H₆O₂ | 122.03678 | loss of a benzoyl group |

> **Important habit:** also look at losses of *44*, *43*, *18*, *30*, *60*, *72*, *42*, *120*,
> *162*. Many of these are **sugars** — glucose is 180, rhamnose 164, hexose 162. Seeing
> a 162 loss means "there is a sugar attached." Seeing 162 **and** 180 **and** 324 means
> "there are two sugars attached."

### TABLE 4 — Diagnostic fragment ions: ring and group "fingerprints"

| m/z | Formula | Means |
|---|---|---|
| 77.04 | C₆H₅⁺ | a plain phenyl ring |
| **91.05** | **C₇H₇⁺** | **tropylium → a benzyl group (Ph–CH₂–)** |
| 93.03 | C₆H₅O⁺ | an oxygen on a phenyl ring (phenol / anisole) |
| **105.03** | **C₇H₅O⁺** | **benzoyl → Ph–CO– (a benzoate or ketone)** |
| 107.05 | C₇H₇O⁺ | oxygenated benzyl |
| 121.03 | C₇H₅O₂⁺ | two oxygens on a 7-carbon ring/chain |
| 121.07 | C₈H₉O⁺ | methoxytropylium |
| **135.04 / 135.07** | — | **methoxy on an aromatic ring (methoxybenzyl)** |
| 149.02 | C₉H₅O₃⁺ | coumarin / benzoyl-derived |
| 163.04 | C₉H₇O₄⁺ | flavonoid-type core |
| 179.03 | C₉H₇O₄⁺ | flavone / flavonoid RDA ion |
| 147.04 | C₉H₇O₂⁺ | methoxyphenylpropenoyl (e.g. feruloyl) |

> **Which ones are arithmetic, which are memorised?**
> The **masses** in this table are arithmetic — you can compute any of them from Table 1.
> The **assignments** ("91 means tropylium") are things you look up and learn. That is the
> only real memorisation in the whole method, and it is a lookup list, not a derivation.

### TABLE 5 — Isotope spacings

| Offset from a peak | Meaning |
|---|---|
| +1.00336 | a ¹³C satellite |
| +0.99705 | a ¹⁵N satellite |
| +1.99840 | an ¹⁸O satellite |
| +1.99580 | a ³⁴S satellite |
| +1.99705 | a ³⁷Cl satellite |
| +1.99795 | an ⁸¹Br satellite |

---

## 2. Step 1 — Get the neutral mass

**Goal:** from the `precursor_mz` and `adduct` columns, compute the neutral molecule's mass.

**The rule:** look up the adduct in Table 2, apply the sign rule.

```
  POSITIVE-MODE adduct  →   M  =  m/z  −  passenger mass
  NEGATIVE-MODE adduct  →   M  =  m/z  +  passenger mass
```

**Worked micro-examples:**

| Given | Compute | Result |
|---|---|---|
| m/z 301.0354, `[M−H]⁻` | 301.0354 **+** 1.00728 | **302.0427** |
| m/z 153.0546, `[M+H]⁺` | 153.0546 **−** 1.00728 | **152.0473** |
| m/z 203.0526, `[M+Na]⁺` | 203.0526 **−** 22.98922 | **180.0634** |
| m/z 349.0704, `[M+CH₂O₂−H]⁻` | 349.0704 **+** 44.99820 | **394.0686** |

**Sanity check you should always do:** neutral masses of small organic molecules are
essentially always between **150 and 500 Da**. If your subtraction gives you 27 Da or
1,900 Da, you misread the adduct column. Go back.

---

## 3. Step 2 — Split into nominal mass and defect

Take the neutral mass and write it as a whole number plus a fraction.

```
  M = 302.0427

      302  │  0.0427
     ──────┴───────
      │        └── DEFECT      = +0.0427
      │
  NOMINAL MASS   = 302
```

**Why you bother:** because the whole-number part is easy to reason about with the
nitrogen rule, and the fraction part is where all the chemistry hides.

---

## 4. Step 3 — Apply the nitrogen rule

> **Odd nominal mass ⟺ odd number of nitrogens. Even nominal mass ⟺ even number of nitrogens.**

302 is **even**, so the molecule has an **even** number of nitrogens: 0, 2, 4, 6…

Natural products frequently have **zero** nitrogens (sugars, terpenes, phenolics,
polyketides). So **try 0 nitrogens first.** It is usually right.

**Second rule — the hydrogen cap:**

```
  max hydrogens  =  2 × (carbons)  +  2  +  (nitrogens)
```

If your candidate formula has more hydrogens than that, it is impossible — the molecule
would need more bonds than carbon has available. Reject it.

---

## 5. Step 4 — Find the molecular formula

This is the first place where the procedure feels like real work, but it is only
**division and subtraction**.

### 5.1 The relationship

From Table 1, the molecule's defect is the sum of its atoms' defects:

```
  defect  =  0.007825·nH  +  0.003074·nN  −  0.005085·nO  −  0.027929·nS  −  ...
             (nH = number of hydrogens, nN = nitrogens, etc.)
```

If you assume no nitrogen and no sulfur (the common case), this becomes:

```
  defect  =  0.007825·nH  −  0.005085·nO
```

**Rearrange it to solve for the hydrogen count:**

```
             defect  +  0.005085·nO
  nH  =  ─────────────────────────────────
                       0.007825
```

### 5.2 The loop you actually run

Try oxygen counts 0, 1, 2, 3, … and each time compute nH, then check whether nH is a
**whole number**. Then check the carbon count is also a whole number:

```
  nC  =  ( nominal mass  −  nH  −  16·nO  −  14·nN )  ÷  12
```

**Worked example: target = 302, defect = +0.04265, assume nN = 0**

Now run the loop. For each oxygen count, compute the hydrogen count, then the carbon count:

```
  nH = (0.04265 + 0.005085·nO) / 0.007825
  nC = (302 − nH − 16·nO) / 12
```

| nO | numerator | nH computed | nH whole? | nC computed | nC whole? |
|---|---|---|---|---|---|
| 0 | 0.04265 | 5.450 | ✗ | 24.750 | ✗ |
| 1 | 0.04774 | 6.100 | ✗ | 23.333 | ✗ |
| 2 | 0.05283 | 6.750 | ✗ | 21.917 | ✗ |
| 3 | 0.05791 | 7.400 | ✗ | 20.583 | ✗ |
| 4 | 0.06299 | 8.050 | ~ | 19.167 | ✗ |
| 5 | 0.06808 | 8.700 | ✗ | 17.750 | ✗ |
| 6 | 0.07316 | 9.350 | ✗ | 16.417 | ✗ |
| **7** | **0.07824** | **10.000** | **✓** | **15.000** | **✓ WIN** |
| 8 | 0.08333 | 10.650 | ✗ | 13.583 | ✗ |
| 9 | 0.08842 | 11.300 | ✗ | 12.250 | ✗ |
| 10 | 0.09350 | 11.950 | ~ | 10.833 | ✗ |

**Exactly one row survives: nO = 7, nH = 10, nC = 15.**

Check the hydrogen cap: 10 ≤ 2(15) + 2 = 32. ✓ Comfortably valid.

So the formula is **C₁₅H₁₀O₇**. Verify by direct sum:

```
  15 × 12.000000  = 180.000000
  10 ×  1.007825  =  10.078250
   7 × 15.994915  = 111.964405
                  ────────────
  total           = 302.042655   →  302.0427   ✓
```

### 5.3 Also try the nitrogen branches

For target 302.0427 with **2 nitrogens**, the two nitrogens contribute +0.006148 of defect, so
the remaining H/O balance must supply 0.04265 − 0.00615 = 0.03650. Then run the same loop.
It produces no whole-number answer (and the carbon count never works out). **Rejected.**

With **4 nitrogens**: remaining defect 0.03035. Same outcome, rejected.

With **1 sulfur**: sulfur's defect is −0.027929, a big negative. You would need
extraordinarily many hydrogens to climb back to +0.04265, and the hydrogen cap
(2c + 2 + n) blocks it. **Rejected.**

> **The takeaway:** because the measurement is good to 5 ppm, and because a defect is a
> 5-decimal number built from only a handful of atoms, **the formula is essentially always
> determined uniquely.** The mass alone — before you look at a single fragment — usually
> tells you the answer's atom inventory. This is why the mass filter is so powerful, and
> why Part 3's retrieval solutions put it first.

### 5.4 The "C₆H₁₂ subtraction" shortcut (saves you time)

Real analysts often do this instead of the loop. Subtract **C₆H₁₂ = 72.0688** repeatedly from
your target and see what falls out.

For 302.0427:
```
  302.0427 − 72.0688 = 229.9739
  229.9739 − 72.0688 = 157.9051
  157.9051 − 72.0688 =  85.8363     ← too small; only two CH2-subtracts valid
```
Then the residue 85.8363 gets built up from C/H/O directly. This is faster for large
molecules and less reliable for small ones. **The loop in 5.2 is safer and I recommend it.**

---

## 6. Step 5 — Read the isotope peaks

Free information. Takes 60 seconds.

### 6.1 The carbon counter

Every carbon has a ~1.1% chance of being ¹³C instead of ¹²C. So a peak at *m/z* *X* is
accompanied by a small peak at *X* + 1.00336, and the size of that satellite tells you how
many carbons are in the fragment:

```
  carbon count  ≈  100 × ( intensity of M+1 peak  ÷  intensity of M peak )
```

| If the M+1 satellite is… | The fragment has about… |
|---|---|
| 1% of M | 9 carbons |
| 5% of M | 45 carbons |
| 10% of M | 90 carbons |

### 6.2 The halogen and sulfur detectors

| You see | Ratio to the main peak | Conclusion |
|---|---|---|
| M+2 at **~3:1** (i.e. 33%) | two isotopes ³⁵Cl (75.8%) and ³⁷Cl (24.2%) | **one chlorine** |
| M+2 at **~1:1** (100%) | two isotopes ⁷⁹Br (50.7%) and ⁸¹Br (49.3%) | **one bromine** |
| M+2 at **~4.4%** | ³⁴S is 4.3% of ³²S | **one sulfur** |
| M+2 at **~18%** | two oxygens (each ~0.2% ¹⁸O) plus… | suggests 2+ oxygens, weak signal |

**Chlorine and bromine are essentially never ambiguous.** Seeing a 3:1 or 1:1 pair tells
you instantly that the molecule contains a halogen, and the ratio tells you which one. This
is the single highest-value 60 seconds of inspection in the whole method.

---

## 7. Step 6 — Compute neutral losses

**This is the heart of the method.**

### 7.1 The procedure

1. Write down the precursor m/z.
2. Take the **~20 strongest** fragment peaks. (Ignore the rest — noise lives in the tail.)
3. For each one, compute `precursor − fragment`. That number is the mass of the **piece that
   broke off**.
4. Look each difference up in Table 3.
5. Each recognisable difference is **evidence that a particular group exists in the molecule.**

### 7.2 The two directions of subtraction

You can subtract in two directions, and you should do both.

| Direction | What it answers |
|---|---|
| **precursor − fragment** | "What group broke off?" → tells you what the molecule **contains** |
| **fragment − smaller fragment** | "What was added or rearranged?" → tells you how the ring **opens** |

### 7.3 Sequential (chained) losses

Very often a peak is reached by losing two things in a row. Do the subtraction **twice**:

```
  301.0354  →  229.0143  →  201.0194

  step 1:  301.0354 − 229.0143  =  72.0211   =  C3H4O2    (a ring fragment left)
  step 2:  229.0143 − 201.0194  =  27.9949   =  CO        (a carbonyl opened)
```

Chains reveal a *sequence*, which is far more informative than isolated peaks. A molecule
that loses CO twice (→ CO twice → nothing) is telling you it has two accessible carbonyls.

### 7.4 The regular-spacing trick

Look for the **differences between all pairs of your top peaks**, not just from the
precursor. A difference that repeats many times is a repeating unit:

| Spacing | Meaning |
|---|---|
| 14.0157 | CH₂ — a saturated alkyl chain |
| 30.0106 | CH₂O — methoxy, or an O in a chain |
| **42.0106** | **C₂H₂O — an acetyl/acyl unit** |
| **162.0528** | **a hexose sugar** |
| 164.0685 | a rhamnose/deoxyhexose sugar |
| 203.0794 | hexose + water |
| 324.1056 | **two hexose sugars** → a diglycoside |

If the same spacing shows up three or four times, you are looking at a chain or a polymer.

### 7.5 When subtraction *fails* — and why that is normal

Sometimes you compute a loss and it matches nothing in Table 3. Example, from
[Section 11](#11-worked-example-a--quercetin-c15h10o7):

```
  301.0354  −  179.0344  =  122.0010
```

122.0010 is not C₇H₆O₂ (122.0368). It is not any simple loss.

**This is not a mistake — it is information.** Ion 179.0344 is the famous flavonoid ion.
It is produced by a **retro-Diels-Alder rearrangement**, where two bonds break *and* the
hydrogens redistribute around the ring. Simple subtraction assumes one bond breaks and
nothing moves; rearrangements move things.

> **Rule: about 20–30% of the peaks in a real spectrum will not be assignable by simple
> subtraction. That is normal and expected.** Do not burn time forcing them.

---

## 8. Step 7 — Recognise diagnostic ions

Check your peak list against Table 4. These are the "smoking gun" peaks.

| If you see | Immediately suspect |
|---|---|
| 91.05 | a benzyl group (a phenyl ring with a CH₂ attached) |
| 105.03 | a benzoyl group (Ph–CO–) |
| 93.03 | a phenyl ring bearing one oxygen |
| 135.04 / 135.07 | a **methoxy group on an aromatic ring** |
| 179.03 | a flavone/flavonoid skeleton |
| 163.04 | a flavonoid-type polyoxygenated core |
| 147.04 | a feruloyl / methoxycinnamoyl unit |
| strong 60-ish | an acetyl or ester side chain |
| 179 / 161 / 143 sequence | a highly oxygenated aromatic |

**Combine this with Table 3.** For example, finding 91.05 (benzyl) *and* a loss of
18.0106 (water) *and* a loss of 162.0528 (hexose) suggests: **a benzyl ring, an OH group,
and a sugar** — i.e. a benzyl glycoside. Now go look for such compounds.

---

## 9. Step 8 — Assemble candidates and test them

### 9.1 Assemble

Using formula + groups + rings, write down **every structure you can imagine** that fits.
Draw them on paper. Typical yield: **2 to 6 candidates.**

Be generous here. A candidate that later fails the test costs you five minutes; a candidate
you never considered costs you the whole molecule.

### 9.2 Test — this is the step people skip, and it's the one that matters

For each candidate, **predict its fragments** and compare to what you observed.

**The prediction method:** for each bond in the candidate that could break, compute the mass
of the piece that would fall off and see whether it is in your loss table. The candidate
that explains the most peaks — especially the *strong* peaks — is your answer.

**Scoring yourself:**

| Peaks explained | Verdict |
|---|---|
| Most, including the base peak | Likely correct |
| Most, but base peak unexplained | Possibly correct; find an alternative pathway |
| Only the small peaks | Probably wrong |
| Nothing | Wrong |

### 9.3 Fill the tail of your shortlist

Once you have 2–5 candidates, add **near-neighbours of the best one**: same formula, one
substituent changed or moved.

Examples for C₁₅H₁₀O₇ (quercetin's formula): the other flavonols of that formula include
**fisetin, morin, kaempferide, myricetin hexoside**, etc. These are the compounds a
confident-but-wrong guess should be followed by.

**Rule:** if your top candidate's formula is X, the other 20 guesses should be the 20 most
abundant compounds in databases with formula X. (Part 3 shows that this "popularity" prior
is worth real leaderboard points.)

---

## 10. Step 9 — Rank and submit

### 10.1 Build the list

1. Your best-tested candidate: rank 1.
2. Other surviving candidates, ranked by how many peaks each explains.
3. Near-neighbour isomers of rank 1.
4. Fill to 25 with the most database-abundant remaining compounds of the right mass.

### 10.2 The submission rules (from the official docs)

| Rule | Detail |
|---|---|
| Every `molecule_id` **exactly once** | no duplicates, no omissions |
| Format | `molecule_id,smiles` — semicolon-separated, best first |
| **Maximum 25** guesses | more is rejected outright |
| Header required | yes |
| No nulls | neither column |
| Filename | `submission.csv` |

Example row from the official docs:
```
m_0014ef,CC1=CC(=O)C=CC1=O;OC(=O)c1ccccc1O;CN1C=NC2=C1C(=O)N(C)C(=O)N2C
```

### 10.3 Two reminders that matter more than they look

- **Always submit 25.** There is *no penalty* for a wrong guess beyond the slot it takes.
- **No stereochemistry needed.** Because scoring is InChIKey14, `OCC1OC(O)C(O)C(O)C1O` is
  a perfectly good glucose prediction with no stereocentres written. Don't waste effort.

---

## 11. Worked Example A — quercetin, C15H10O7

**This is the whole method, start to finish, with nothing skipped.**

### 11.1 The data we're given

From a test-like row:

```
  precursor_mz                = 301.0354
  adduct                      = [M-H]-
  ionization_mode             = negative
  collision_energy_ev         = [40]

  ms2_mzs                    =  [301.0354, 283.0248, 273.0405, 271.0248, 257.0456,
                                  241.0142, 229.0143, 213.0557, 178.9986, 151.0037,
                                  123.0088, 121.0290,  95.0491,  65.9979]

  ms2_normalized_intensities = [1.00,    0.05,    0.34,    0.62,    0.31,
                                  0.04,    0.44,    0.03,    0.38,    0.09,
                                  0.02,    0.11,    0.08,    0.10]
```

**Read the intensities first.** The base peak (intensity 1.00) is the **precursor itself,
301.0354**. The strongest *fragments* are:

| m/z | intensity | note |
|---|---|---|
| 271.0248 | 0.62 | second strongest |
| 229.0143 | 0.44 | third |
| 178.9986 | 0.38 | fourth |
| 273.0405 | 0.34 | |
| 257.0456 | 0.31 | |
| 121.0290 | 0.11 | |
| 283.0248 | 0.05 | weak |

**Do your loss arithmetic on those strong peaks first.** Working on the weak tail is how
analysts waste hours. Only go to the weak peaks once the main skeleton is settled.

### 11.2 Step 1 — neutral mass

```
  adduct is [M-H]-  →  NEGATIVE mode  →  ADD the passenger

  M  =  301.0354  +  1.00728  =  302.0427
```

Sanity check: 302 is in the 150–500 range for a small natural product. ✓

### 11.3 Steps 2–4 — the formula

```
  nominal mass  =  302          defect  =  +0.0427
  302 is EVEN   →  even number of nitrogens → try 0 first
```

Run the loop from [Section 5.2](#52-the-loop-you-actually-run). **Exactly one row survives:**

```
   nO=7  →  nH = 10.000 (exact!)  →  nC = 15.000 (exact!)

   hydrogen cap:  10 ≤ 2(15) + 2 = 32  ✓
   direct check:  15(12) + 10(1.007825) + 7(15.994915) = 302.0427  ✓
```

### FORMULA = C₁₅H₁₀O₇

(For the record, the branches with 2 and 4 nitrogens and the branch with 1 sulfur were all
run and all produced no whole-number answer. Details in [Section 5.3](#53-also-try-the-nitrogen-branches).)

### 11.4 Step 5 — isotopes

Look for a peak 1.00336 above any strong peak. There is a peak at **273.0405** and one at
**273.0405 − 1.003 = 272.037**? Not in our list, so no useful satellite here. And no 3:1 or
1:1 M+2 pairs anywhere, so:

> **No chlorine. No bromine. No sulfur.** Consistent with C₁₅H₁₀O₇ — good, the isotopes
> confirm the formula rather than contradict it.

### 11.5 Step 6 — the loss table

Subtract every significant fragment from the precursor, then look up Table 3.

```
   fragment    precursor − fragment    Table 3 match      strength
  ─────────────────────────────────────────────────────────────────────
   283.0248        18.0106             H2O  (18.0106)      0.05  weak
   273.0405        27.9949             CO  (27.9949)        0.34  STRONG
   271.0248        30.0106             CH2O (30.0106)       0.62  BASE-ish
   257.0456        43.9898             CO2 (43.9898)        0.31  strong
   241.0142        60.0212             C2H4O2 (60.0211)     0.04  weak
   229.0143        72.0211             C3H4O2 (72.0211)     0.44  STRONG
   213.0557        87.9797             2 × CO2 (87.9797)    0.03  weak
   178.9986       122.0368             C7H6O2 (122.0368)    0.38  strong
    95.0491       205.9863             (no match — low-mass, see 11.6)
    65.9979       235.0375             (no match)
```

**Every single assignment above is exact to 4 decimal places.** That is not luck — it is
what makes this method work.

Now the **sequential** losses:

```
   229.0143 − 201.0194 = 27.9949 = CO    (so 201.0194 is a real peak in richer spectra)
   178.9986 − 151.0037 = 27.9949 = CO    ← our list HAS 151.0037
   151.0037 − 123.0088 = 27.9949 = CO    ← and 123.0088
```

So there is a chain: **178.9986 → 151.0037 → 123.0088, losing CO each step.** Three ions
down a path, each one carbon monoxide gone. That means the neutral contains a **sequence of
carbonyl-bearing positions that can open one after another** — a hallmark of a fused
polyoxygenated ring system.

### 11.6 Step 7 — diagnostic ions

```
   178.9986  ≈  C9H7O4+  (179.0344)   ← a FLAVONOID ion
   121.0290  =  C7H5O2+  exactly       ← a two-oxygen 7-carbon ring system
    95.0491  =  C6H7O+   exactly
```

And crucially, Table 4 flags **179.03 as the flavone/flavonoid RDA ion.** Our
`precursor − 179.0344 = 122.0010` matches nothing — which per
[Section 7.5](#75-when-subtraction-fails--and-why-that-is-normal) is exactly what a
retro-Diels-Alder rearrangement looks like. It is a **positive** confirmation that the
flavonoid hypothesis is right, not a puzzle to solve.

### 11.7 What we have concluded so far

| Evidence | Points to |
|---|---|
| Formula C₁₅H₁₀O₇ | 15 C, 10 H, 7 O |
| Loss of CO₂ | a **carboxylic acid** — *or*, see below |
| Loss of 30.0106 (CH₂O), strong | an **oxygen in a ring adjacent to carbon**, ring opening |
| Loss of CO, three times in a chain | a fused **polycyclic** system with accessible carbonyls |
| 179.03 / 121.03 | a **flavonoid / flavonol skeleton** |
| No halogens, no sulfur | C/H/O only ✓ |

**Candidate class: flavonols** — compounds built from two benzene rings joined through a
six-membered oxygen-containing ring with a carbonyl, and heavily decorated with –OH.

### 11.8 Step 8 — assemble candidates

Flavonols with formula C₁₅H₁₀O₇. Real ones:

| Compound | OH positions | Notes |
|---|---|---|
| **quercetin** | 3, 5, 7, 3′, 4′ | five hydroxyls, B-ring catechol |
| fisetin | 3, 7, 3′, 4′ | no 5-OH, no 6-ring OH |
| morin | 3, 5, 7, 2′, 4′ | 2′-OH instead of 3′ |
| kaempferide | 3, 5, 7, 4′ | four OH, no catechol on B-ring |
| herbacetin | 3, 5, 6, 7, 4′ | extra 6-OH |

**Plus the non-flavonol alternatives** you should not skip: gossypetin (a flavonol with
6 OH), and gallocatechin-type structures. Write them all down.

### 11.9 Step 8 — test each candidate

**Test quercetin.** Quercetin has **five –OH groups**. It should be able to lose water —
and indeed we see 301.0354 → 283.0248 (loss H₂O). ✓

Quercetin's B-ring is a **catechol** (two adjacent –OH on the pendant phenyl). Losing
**C₇H₆O₂ (122.0368)** — which is exactly a dihydroxyphenyl fragment — is what we observe at
178.9986. ✓✓

Quercetin's C-ring carries a **carbonyl at C4**. Losing **CO** is expected. ✓ And the
three-step CO chain (178.9986 → 151.0037 → 123.0088) is consistent with sequential opening
of the fused system. ✓

**Test fisetin** (3,7,3′,4′-tetrahydroxy). It has the catechol too, so it also explains
122.0368. It has a 5-**H** rather than a 5-**OH**, so it is slightly less able to lose
water — but water loss is a weak peak here anyway, so this is **not decisive**.

**Test kaempferide** (3,5,7,4′-tetrahydroxy). Its B-ring has only **one** OH, so it cannot
lose C₇H₆O₂ as readily — it would more likely lose C₇H₆O₂ with one fewer O. The observed
122.0368 loss has **two oxygens in the leaving group**, which favours the catechol (B-ring
with two OH). Slight evidence for quercetin over kaempferide. ✓

**Test morin.** Morin has its 2′-OH adjacent to the ring junction, which *sterically
hampers* the O on that ring from participating in loss. Slightly disfavoured.

**Test gossypetin** (3,5,6,7,8-pentahydroxyflavone — a flavonol with OH on the A-ring).
Strongly oxygenated; should show extra water losses. We only see a weak one. **Disfavoured.**

**Verdict: quercetin is the front-runner, with fisetin a close second.** The data cannot
absolutely exclude fisetin — which is exactly the situation the metric anticipates.

### 11.10 Step 9 — build the submission

```
  Rank   1  quercetin                              (best explanation)
  Rank   2  fisetin                                (same formula, same catechol loss)
  Rank   3  morin                                  (same formula, near-neighbour)
  Rank   4  kaempferide                            (same formula)
  Rank   5  gossypetin                             (same formula)
  Ranks 6–25  other C15H10O7 flavonols and
            near-neighbours, most-abundant first
```

**And now the trick that wins real points:** ranks 6–25 should be filled with the
**most commonly reported** C₁₅H₁₀O₇ compounds in PubChem — because when the spectrum
genuinely cannot discriminate, the answer is more likely to be a famous compound than an
obscure one. Part 3 quantifies this: that single idea is worth roughly **+0.01 MRR** on the
real leaderboard, which is enormous.

---

## 12. Worked Example B — benzoic acid (two adducts, two answers)

Short example, to show how **different adducts reveal different things**.

### 12.1 The data

We are given **two spectra of the same molecule**, measured two different ways.

```
  Spectrum 1:  precursor 121.0295,  adduct [M-H]-   (negative mode)
               peaks: 121.0295 (1.00), 105.0335 (0.22), 77.0397 (0.61), 50.0023 (0.03)

  Spectrum 2:  precursor 123.0441,  adduct [M+H]+   (positive mode)
               peaks: 123.0441 (1.00), 105.0335 (0.55),  77.0391 (0.09), 51.0229 (0.02)
```

### 12.2 Neutral mass — must be consistent

```
  Spectrum 1:  121.0295  +  1.00728  =  122.0368
  Spectrum 2:  123.0441  −  1.00728  =  122.0368      ✓  SAME
```

**This consistency check is free and powerful.** All of a molecule's spectra must yield the
same neutral mass. If they don't, one of them is mis-assigned — and you learn that without
any model at all.

```
  M = 122.0368  →  nominal 122, defect +0.0368
```

### 12.3 The formula

```
  defect +0.0368,  nN = 0 (122 even)

   nO   nH = (0.0368 + 0.005085·nO)/0.007825      nH whole?   nC = (122 − nH − 16·nO)/12
  ─────────────────────────────────────────────────────────────────────────────────────────
   1         0.04189/0.007825 = 5.352               ✗           10.067     ✗
   2         0.04697/0.007825 = 6.002   ≈ 6         ✓ 6          7.000     ✓✓  WIN
   3         0.05206/0.007825 = 6.652               ✗
   4         0.05714/0.007825 = 7.302               ✗
```

**nO=2, nH=6, nC=7 → C₇H₆O₂.** Verify:
```
  7(12) + 6(1.007825) + 2(15.994915) = 84 + 6.04695 + 31.98983 = 122.03678   ✓
```

### 12.4 Now the fragmentation — the *best* lesson in the document

**Spectrum 1, negative mode:**

```
   121.0295 − 77.0397 = 43.9898   =  CO2     ← Table 3: a CARBOXYLIC ACID
   121.0295 − 105.0335 = 15.9960  =  ~16     ← a deprotonated/rearranged form
```

**Spectrum 2, positive mode:**

```
   123.0441 − 105.0335 = 18.0106   =  H2O     ← Table 3: an OH / acid
   123.0441 −  77.0391 = 46.0050   =  HCOOH   ← the whole acid group, neutral
```

**This is a perfect, two-sided confirmation of one carboxylic acid:**

- In **negative** mode, the acid readily loses **CO₂** (43.9898) — the textbook behaviour.
- In **positive** mode, the acid loses **H₂O** (18.0106) and the neutral **formic acid
  moiety** (46.0055), leaving the benzoyl cation **C₇H₅O⁺ = 105.0340**.

And 105.0340 appears in **both** spectra. It is Table 4's **benzoyl** ion — definitive
evidence for a **benzene ring bearing a carbonyl**.

Combine: benzene ring + carbonyl + CO₂ loss + C₇H₆O₂ formula.

> ### Answer: **benzoic acid**, C₆H₅COOH.

**What this example teaches:** the same molecule measured two ways gives you *independent,
complementary* evidence. Negative mode told you about CO₂; positive mode told you about
water and the benzoyl ring. **Using all of a molecule's spectra together roughly doubles
your confidence.** This is why the competition gives you 1–16 spectra per molecule.

---

## 13. Worked Example C — vanillin (finding the aldehyde and the methoxy)

```
  precursor 153.0546,  adduct [M+H]+

  M = 153.0546 − 1.00728 = 152.0473        nominal 152, defect +0.0473
  152 is even → try 0 nitrogens

   nO   nH = (0.0473 + 0.005085·nO)/0.007825     nH whole?    nC = (152 − nH − 16·nO)/12
  ────────────────────────────────────────────────────────────────────────────────────────
   1         0.05239/0.007825 = 6.699              ✗            9.917    ✗
   2         0.05739/0.007825 = 7.334              ✗            9.250    ✗
   3         0.06248/0.007825 = 7.984   ≈ 8        ✓ 8           8.000    ✓✓  WIN
   4         0.06756/0.007825 = 8.634              ✗            6.750    ✗
```

### FORMULA = C₈H₈O₃   (verify: 96 + 8.0626 + 47.9847 = 152.0473 ✓)

Peaks: `153.0546 (1.00), 138.0311 (0.10), 125.0597 (0.48), 109.0284 (0.35), 107.0491 (0.14), 93.0335 (0.42), 81.0335 (0.11), 65.9979 (0.06)`

### 13.1 The losses

```
   153.0546 − 138.0311 =  15.0235  =  CH3 radical   ← a METHYL group leaving, alone
   153.0546 − 125.0597 =  27.9949  =  CO            ← a CARBONYL, ejected whole
   153.0546 − 109.0284 =  44.0262  =  C2H4O         ← acetaldehyde
   125.0597 − 107.0491 =  18.0106  =  H2O           ← (from the CO-loss ion)
   125.0597 −  93.0335 =  32.0262  =  CH3OH         ← (from the CO-loss ion)
```

### 13.2 Reading the structure out of this

**(a) The 15.0235 loss = a methyl radical leaving by itself.** A lone CH₃⁺ radical only
happens when it is bonded to something electronegative that stabilises the remaining
positive charge — classically **an –OCH₃ methoxy group**, where the oxygen can hold the
charge. Vanillin has one methoxy. ✓

**(b) The 27.9949 (CO) loss = a carbonyl ejected alone.** Aldehydes (–CHO) famously
leave as CO. So there is an **aldehyde group**. Vanillin has one. ✓

**(c) The 93.0335 ion = C₆H₅O⁺, exactly.** Table 4 says: an oxygen on a phenyl ring. And
how did we get there? By losing **CH₃OH (methanol)** from the CO-loss ion. That means the
methoxy group left **as intact methanol** — which is precisely what an aromatic methoxy
does when the adjacent ring can stabilise the charge. The resulting C₆H₅O⁺ is the classic
**methoxyphenyl / phenoxy** ion. ✓✓

**(d) The 109.0284 ion = loss of C₂H₄O.** A two-carbon-oxygen unit leaving together.
Consistent with losing HCOOH + something, or acetaldehyde from the –CH₂OH– ring position.
Vanillin's ring has a CH₂OH after the phenolic OH loses its H. Reasonable.

### 13.3 Assemble

- Formula C₈H₈O₃
- 3 oxygens
- an **aldehyde** (from the CO loss)
- a **methoxy** (from the CH₃ radical loss and the CH₃OH loss)
- an aromatic ring (from C₆H₅O⁺ at 93.03)
- the aromatic ring bears **at least two oxygen substituents** (the aldehyde contributes one
  C=O, and there's a phenolic/methoxy oxygen)

### ANSWER: **vanillin**, 4-hydroxy-3-methoxybenzaldehyde, `COc1cc(C=O)ccc1O`

**Verify against the losses:**
- methoxy → CH₃⁺ loss and CH₃OH loss ✓
- aldehyde → CO loss ✓
- phenolic OH → H₂O loss ✓
- aromatic + OMe → C₆H₅O⁺ at 93.03 ✓

**Alternatives considered and rejected:** isovanillin (3-hydroxy-4-methoxy) has the same
formula and the same functional groups and would produce nearly identical losses — **this
one genuinely cannot be resolved by losses alone.** It is a perfect rank-2 candidate, and
exactly why you submit 25.

---

## 14. Worked Example D — the isomer trap (why you submit 25)

This is the most important example in the document, and it is a *negative* result: a case
where **the manual method runs out of road.**

### 14.1 The setup

Four molecules:

```
  A  methyl benzoate        C6H5–C(=O)–O–CH3
  B  phenyl acetate         C6H5–O–C(=O)–CH3
  C  benzyl formate         C6H5–CH2–O–CHO
  D  2-hydroxyacetophenone  C6H5–C(=O)–CH2OH
```

### 14.2 First: they are the same mass

```
  8 × 12.000000  =  96.000000
  8 ×  1.007825  =   8.062600
  2 × 15.994915  =  31.989830
                 ────────────
  C8H8O2         = 136.052430 Da     ← all four, identically
```

Any mass-based filter returns all four. Any mass-based formula step returns **C₈H₈O₂** and
therefore all four.

### 14.3 Second: they all produce the same precursor

```
  [M+H]+ = 136.05243 + 1.00728 = 137.05971        for all four
```

### 14.4 Third: the losses separate them — barely

```
   A  methyl benzoate:  137.0597 − 32.0262 (CH3OH) = 105.0335  = C7H5O+  BENZOYL
   B  phenyl acetate:   137.0597 − 42.0106 (CH2CO) =  95.0491  = acylium
   C  benzyl formate:   137.0597 − 46.0055 (HCOOH) =  91.0542  = C7H7+   TROPYLIUM
   D  2-hydroxyacetophenone: 137.0597 − 16.0307 (CH4) = 121.0290 = C7H5O2+
```

Four different strong ions. **If your spectrum is clean, you can tell them apart.**

### 14.5 But now suppose the spectrum is noisy

| Reality check | Consequence |
|---|---|
| The discriminating ions are **91, 95, 121** — all in the crowded 90–125 region | Many co-isolated contaminants live there |
| The precursor column says spectra with `base_peak_intensity` under 1,000 raw counts were dropped as noise | Surviving spectra are *relatively* clean, but "relatively" is doing a lot of work |
| Test molecules have a **median of only 3 spectra** | You often have too few chances to see the one ion that discriminates |
| Real spectra have heavy-tailed peak counts | Some of your peaks are noise, and a noise peak at 95.05 looks exactly like signal |

**Now suppose 95.0491 is present at 4% intensity and 91.0542 is absent because it was
buried in the noise floor.** Now:

- Is it phenyl acetate (95) or benzyl formate (91 missing)?
- Or is the 95 a noise peak and the real answer is methyl benzoate?

**You cannot tell. The data does not contain the information.**

### 14.6 What the manual method hands you

```
  Rank 1  benzyl formate      (loss of formic acid 46.0055 is diagnostic; absence of
                               105 rules out A; the unassigned 91-sized region is weak)
  Rank 2  phenyl acetate
  Rank 3  methyl benzoate
  Rank 4  2-hydroxyacetophenone
  Rank 5–25  other C8H8O2 esters and ketones, most-abundant first
```

Which is **exactly right**, because the metric pays 1/rank and this ordering captured all
the probability mass.

### 14.7 The lesson

> **The manual method's real output is not "the answer." It is a *distribution* over
> answers.** Its job is to move 95% of the probability onto 25 candidates. Anything beyond
> that is not information, it's a coin flip.

And that is precisely why:
- the metric is **MRR@25** and not accuracy,
- **no penalty** exists for wrong guesses,
- and "always submit 25" is the single most important submission habit.

---

## 15. Using multiple spectra of the same molecule

The test set gives **1–16 spectra per molecule, median 3.** Different adducts, different
collision energies, sometimes merged.

### 15.1 How to combine them

| Situation | What to do |
|---|---|
| Same adduct, different collision energies | Run the procedure on the **highest-energy** one first (most fragmentation → most structural detail), then use the others to confirm or add groups |
| **Different adducts** | Run each separately. Different adducts give **different losses** — the benzoic acid example. **Cross-check the neutral masses agree.** |
| Any spectrum that disagrees on the neutral mass | **Discard it as mis-assigned.** Free quality control. |
| A low `base_peak_intensity` spectrum | **Downweight it.** It is mostly noise. |
| Any spectrum that contradicts your candidate | That candidate is suspect. Demote it. |

### 15.2 The aggregation rule that actually matters

Do **not** just average the lists. Do this:

> **A candidate must be consistent with EVERY spectrum of the molecule.** Score each
> candidate per spectrum, then take the **minimum** (or a low quantile) across spectra —
> not the mean. One spectrum that flatly contradicts a candidate should disqualify it,
> and a mean would let four agreeable spectra outvote the one that knows better.

This is a real design decision, and getting it backwards is a common bug.

---

## 16. What this method cannot do

Honesty, so you know where to stop.

### 16.1 Where the manual method runs out

| Situation | Why it fails |
|---|---|
| Isomers with identical fragmentation (Section 14) | The information is not in the data |
| Very large molecules (mass > 500) | Too many possible formulas and skeletons; the search explodes |
| Novel scaffolds with no known fragment rules | Table 3 and Table 4 are built from known chemistry; a new skeleton has no entry |
| Noisy spectra | The discriminating ions are below the noise floor |
| Molecules with 3+ sugars | Combinatorial explosion in how they're linked |
| Stereochemistry | Not scored, so irrelevant — lucky |

### 16.2 Where the manual method is actually *better* than a neural network

Don't assume ML wins. It often doesn't:

| Task | Manual | Neural network |
|---|---|---|
| Getting the **formula** right | **Instant, exact** | Must be learned; often wrong |
| Assigning a **known** fragment loss | Exact, interpretable | Statistical, may hallucinate |
| **Handling a novel instrument** | Needs new rules | Needs retraining |
| Explaining a decision to a human | Perfect | Opaque |
| Ruling out a candidate with one decisive fact | Instant | Needs evidence to accumulate |

> **The strongest practical systems use both.** The mass and formula come from arithmetic
> (exact). The candidate pool comes from a database (retrieval). The neural network is used
> to *rank* among near-ties, never to *decide* the formula. Part 3 shows that the
> leaderboard-leading solutions do exactly this — they start with the mass filter and never
> let a model override it.

### 16.3 The one thing you cannot do by hand

**Enumerate the candidate pool.** A human can hold a few dozen molecules in their head.
A computer can hold 700,000. Since C₁₅H₁₀O₇ has many hundreds of possible flavonols, only
a database search finds them all. **That transition — from "reason about one structure" to
"rank 700,000 structures" — is exactly the step that needs machine learning, and it is
where the leaderboard scores come from.**

---

## 17. One-page cheat sheet

```
  ┌────────────────────────────────────────────────────────────────────────────┐
  │ 1. NEUTRAL MASS     M = precursor_mz − adduct   (positive mode)            │
  │                     M = precursor_mz + adduct   (negative mode)            │
  │                     Cross-check all spectra of one molecule agree.         │
  ├────────────────────────────────────────────────────────────────────────────┤
  │ 2. SPLIT            nominal = whole part      defect = fractional part     │
  ├────────────────────────────────────────────────────────────────────────────┤
  │ 3. NITROGEN RULE    odd nominal  <-> odd  #N                                │
  │                     even nominal <-> even #N     try 0 N first              │
  │                     H cap:  nH ≤ 2·nC + 2 + nN                              │
  ├────────────────────────────────────────────────────────────────────────────┤
  │ 4. FORMULA          nH = (defect + 0.005085·nO) / 0.007825      [nN=0]     │
  │                     loop nO = 0,1,2,…  keep only whole nH AND whole nC     │
  │                     re-verify by direct summation                           │
  ├────────────────────────────────────────────────────────────────────────────┤
  │ 5. ISOTOPES         M+1 ≈ 1.1% per carbon  →  carbon count                 │
  │                     M+2 at 3:1   → Cl        M+2 at 1:1 → Br                │
  │                     M+2 at ~4.4% → S                                        │
  ├────────────────────────────────────────────────────────────────────────────┤
  │ 6. LOSSES           loss = precursor − fragment   (look up)                 │
  │   18.0106  H2O     27.9949  CO     30.0106  CH2O    43.9898  CO2          │
  │   15.0235  CH3•    32.0262  CH3OH  60.0211  C2H4O2  162.05   hexose       │
  │   then subtract AGAIN between fragments — chains are informative           │
  │   ~25% of peaks will not assign. That is normal (rearrangements).          │
  ├────────────────────────────────────────────────────────────────────────────┤
  │ 7. DIAGNOSTICS      91.05 benzyl   105.03 benzoyl   93.03 Ph–O             │
  │                     135.04 OMe-on-ring   179.03 flavonoid   147.04 feruloyl│
  ├────────────────────────────────────────────────────────────────────────────┤
  │ 8. ASSEMBLE         draw 2–6 candidates from formula + groups + rings       │
  │                     TEST each: does it predict the peaks you actually see?  │
  ├────────────────────────────────────────────────────────────────────────────┤
  │ 9. SUBMIT           rank best first, fill to 25 with same-formula near-    │
  │                     neighbours, most-abundant-in-databases first            │
  │                     NO stereochemistry needed (InChIKey14 ignores it)        │
  │                     ALWAYS 25 guesses — no penalty for wrong ones          │
  └────────────────────────────────────────────────────────────────────────────┘
```

---

**Next:** [Part 3 — How People Actually Solve It](03-how-people-actually-solve-it.md) — the
real competition leaderboard, the exact pipeline behind the top public solutions with
configurations and scores, the published researchers' methods (DiffMS, MS-GPT, MSFlow,
MARLIN and others), and a build order you can follow.
