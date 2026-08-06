# Exact within-cross additive variance for residual-heterozygous parents

**Status:** DH **and RIL-∞** derivations **complete, validated, and implemented**
(`ng_gms_additive_var_general`, wired into `ng_score_crosses` for
`parent_type="ril"` with phased haplotypes). Gives the exact analytical
within-cross additive variance for *arbitrary phased diploid parents* — fully
inbred **or** residual-heterozygous (real RILs). Reduces exactly to the inbred
`a'Ra` kernel; validated to Monte-Carlo precision against the assumption-free
pure-Haldane transmission across coupling/repulsion/coexisting-Δ,δ/both-het/random
configurations (max |err| 0.31% at 3×10⁵ reps; DH −0.24%, RIL-∞ +0.49%). Only
**finite-generation** RILs (residual *progeny* het) remain simulation-only.

## Problem

The production kernel `VPM = a'Ra`, `a_k = ½(x1_k − x2_k)β_k`,
`R_kl = 1 − 2r_kl` (DH), assumes both parents are fully homozygous
(`x ∈ {0,2}`). Real RIL parents are inbred *by descent* but carry residual
heterozygosity; at those loci the parent segregates and contributes a
Mendelian-sampling term the kernel omits, biasing the variance **low**. We want
the exact variance as a function of the *phased* parental haplotypes.

## Scope: DH and RIL are two independent axes

Heterozygosity enters at two *separate* points; conflating them is the usual
error.

1. **Progeny target** (what the cross produces). A **DH** family is fixed —
   every DH line is 100% homozygous by construction (single gamete, doubled).
   A **RIL** family (finite selfing) may still carry residual heterozygosity in
   individual lines. This axis chooses the kernel: DH `1−2r`, RIL-∞
   `(1−2r)/(1+2r)`, finite-RIL by simulation.

2. **Parent line type** (what is being crossed). **DH / fully fixed parents are
   homozygous** — a heterozygous locus in a declared-DH parent is a genotyping
   or data error and must be a **blocker (do not proceed)**. **RIL parents**
   legitimately carry residual het at a few loci and must be **allowed**, with
   the het-aware variance below.

This derivation lives at the intersection **RIL (residual-het) parents →
DH progeny**: the parents segregate (Term 2), yet every DH progeny line is
itself fully homozygous. So "residual het" here is a property of the *parents*,
never of the DH progeny. Governance: **`parent_type = c("inbred","dh","ril")`**
(default `"inbred"`). `"dh"` declares doubled haploids — 100% homozygous by
construction — and gets **zero** het-marker tolerance: *any* true het call
blocks (wrong ploidy, contamination, a RIL mislabelled DH). `"inbred"` (finished
inbred lines) blocks het but keeps a small fraction tolerance for genotyping
noise (`inbred_marker_fraction`, default 2%). `"ril"` declares RIL inputs and
**proceeds** with the het-aware kernel. (The per-dosage numeric tolerance
`inbred_tolerance`, e.g. `1.998 → 2`, applies in all cases.) The legacy boolean `assume_inbred` is
deprecated (reconciled `TRUE→inbred`, `FALSE→ril` with a one-time warning). The
two axes are orthogonal — RIL parents can produce DH progeny (this doc) or RIL
progeny (RAH simulation), and `parent_type` (audit) is independent of the
progeny `target` (kernel).

## Notation

Loci `k = 1..m`; effects `β_k`; recombination frequency `r_kl` between loci
`k, l` (net, from the map). Phased parents:

- P1 haplotypes `a = (a_k)`, `b = (b_k)`, alleles in `{0,1}`.
- P2 haplotypes `c = (c_k)`, `d = (d_k)`.

Derived quantities:

- dosages `x1_k = a_k + b_k`, `x2_k = c_k + d_k`  (∈ {0,1,2})
- allele frequencies `p1_k = ½(a_k+b_k)`, `p2_k = ½(c_k+d_k)`
- parent contrasts `Δ_k = p1_k − p2_k = ½(x1_k − x2_k)`
- **within-parent haplotype contrasts** `δ1_k = a_k − b_k`, `δ2_k = c_k − d_k`
  (∈ {−1,0,1}; **nonzero only where that parent is heterozygous**)
- DH kernel `φ_kl = 1 − 2r_kl`

## Transmission model (DH)

`P1 → gamete u` (recombination mosaic of `a,b`); `P2 → gamete v` (mosaic of
`c,d`); `F1 = (u,v)`; `F1 → gamete w` (mosaic of `u,v`); DH line `= (w,w)`,
dosage `X_k = 2w_k`. Genetic value `G = Σ_k β_k X_k = 2 Σ_k β_k w_k`, so

```
Var(G) = 4 Σ_{k,l} β_k β_l Cov(w_k, w_l).
```

`w_k` traces through three independent meioses: an F1-meiosis indicator `S_k ∈
{P1,P2}` (which F1 side), and parental indicators `T1_k, T2_k` (which parental
haplotype), each a symmetric 2-state Markov chain along the chromosome.

## Key expectations

`E[u_k] = p1_k`, `E[v_k] = p2_k`, `E[w_k] = ½(p1_k + p2_k)` ⇒
`E[X_k] = p1_k + p2_k = ½(x1_k + x2_k)` (mid-parent). ✔

For the F1-meiosis, `P(S_k=S_l=P1) = P(S_k=S_l=P2) = ½(1−r_kl)` and the
cross-side cases each have probability `½ r_kl`. Conditioning `E[w_k w_l]` on
`(S_k,S_l)`:

```
E[w_k w_l] = ½(1−r_kl) E[u_k u_l] + ½(1−r_kl) E[v_k v_l]
           + ½ r_kl (p1_k p2_l + p2_k p1_l).
```

### Within-parent covariance (the new ingredient)

A parental gamete `u` is a mosaic of `a,b` with meiosis recombination `r_kl`.
Carrying out `E[u_k u_l] − p1_k p1_l` gives, after cancellation,

```
Cov_u(k,l) = ¼ (1 − 2r_kl) (a_k − b_k)(a_l − b_l) = ¼ φ_kl δ1_k δ1_l,
Cov_v(k,l) = ¼ φ_kl δ2_k δ2_l.
```

Zero unless the parent is heterozygous at **both** loci — exactly the residual
segregation the inbred kernel drops.

## Assembling Cov(w_k, w_l)

The frequency-product part collapses (same algebra as the within-parent step) to
`¼ φ_kl Δ_k Δ_l`. The within-parent parts carry the same-side coefficient
`½(1−r_kl)`:

```
Cov(w_k,w_l) = ¼ φ_kl Δ_k Δ_l  +  ½(1−r_kl)[Cov_u(k,l) + Cov_v(k,l)]
             = ¼ φ_kl Δ_k Δ_l  +  ⅛(1−r_kl) φ_kl (δ1_k δ1_l + δ2_k δ2_l).
```

Multiplying by `4 β_k β_l` and summing:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  Var(G)_DH =  a'R a   +   ½ Σ_{k,l} β_k β_l (1−r_kl) φ_kl (δ1_k δ1_l + δ2_k δ2_l) │
└─────────────────────────────────────────────────────────────────────────────┘
```

with `a_k = ½(x1_k − x2_k)β_k` and `R_kl = φ_kl = 1 − 2r_kl`. As one quadratic
form `Var(G)_DH = β' K β`:

```
K_kl = R_kl [ ¼(x1_k−x2_k)(x1_l−x2_l) + ½(1−r_kl)(δ1_k δ1_l + δ2_k δ2_l) ].
```

- **Term 1** `a'Ra` is the classic inbred DH variance (parent-difference `Δ`).
- **Term 2** is the residual-het correction: within-parent haplotype contrasts
  `δ1, δ2`, weighted by `(1−r_kl)` (probability both loci co-inherit from the
  same F1 side) times the DH kernel `φ_kl`. Its parent-marginal
  `Var_{u,v}[E(G|F1)] = β'(L_P1+L_P2)β` is precisely what
  `ng_exact_gms_additive_var()` already computes — i.e. that orphaned function
  equals the **between-F1** variance component of this decomposition. (It is
  *not* numerically half of Term 2: the coefficients `¼φ` vs `½(1−r)φ` differ,
  so `2 × ng_exact_gms` reconstructs Term 2 only at `r = 0`.)

## Sanity checks

**Reduces to inbred `a'Ra`.** Inbred parents ⇒ `δ1 = δ2 = 0` ⇒ Term 2 vanishes,
`Var = a'Ra`. ✔

**Single locus, P1 = Aa × P2 = aa** (`x1=1, x2=0, δ1=1, δ2=0`, `r=0`, `β=1`):

```
Term 1 = [½(1−0)]² · 1 = ¼ ;  Term 2 = ½·(1−0)·(1)·(1²) = ½ ;  Var = ¾.
```

Matches the exact transmission-tree value `¾β²`, and the RAH lineage simulation
measured `0.744` for this cross (`tests/rah_residual_het_parent.R`). ✔

**Single locus, P1 = AA × P2 = aa:** `δ1=δ2=0`, Term 1 `= 1`, Var `= 1`. ✔

**Multi-locus, numerical.** 5 loci, one chromosome, effects `β`, P1 het at two
loci; DH truth by simulation:

Against the assumption-free pure-Haldane transmission MC (three validated
`ng_rah_meiosis_gamete` draws, 3×10⁵ reps):

| Configuration | Closed | MC | Rel. err |
| --- | ---: | ---: | ---: |
| P1=Aa × P2=AA (Δ and δ coexist at a locus) | 0.7500 | 0.7523 | −0.31% |
| Coupling linked het (both `1` on one homolog) | 2.5995 | 2.6008 | −0.05% |
| **Repulsion** linked het (same loci, opposite homologs) | 1.7037 | 1.7038 | −0.00% |
| Both parents heterozygous, mixed phase | 1.5141 | 1.5126 | +0.10% |
| Random het configuration (6 loci) | 3.7636 | 3.7546 | +0.24% |

Two checks matter most for correctness: (i) coupling (2.60) vs repulsion (1.70)
are *identical loci and effects* differing only in phase — the signed
`δ_kδ_l` term reproduces both, so the phase-dependence is exact (and is why
**phased** haplotypes are required, not dosages); (ii) `Δ` and `δ` coexisting
at one locus needs **no** `Δ×δ` cross-term — confirmed. AlphaSimR agrees for
*unlinked* het (−0.33%) but sits ~5% above for *linked* het because its default
gamma crossover model carries **interference**, whereas the package (and the
`a'Ra` kernel it extends) assume Haldane (no interference). The formula is
therefore exact **for the package's recombination model**; interference is a
separate, package-wide modelling assumption, not a defect of this derivation.

The closed form matches a pure-Haldane Monte-Carlo of the exact transmission
(three validated `ng_rah_meiosis_gamete` draws: P1→gamete, P2→gamete, F1→DH) to
MC precision — it is **exact for the package's Haldane (no-interference)
recombination model**. The ~5% gap against AlphaSimR is AlphaSimR's default
crossover **interference** (gamma model, fewer double-crossovers → linked loci
stay more correlated), not a formula error; the same interference offset appears
in the inbred baseline. For reference, the inbred-only `a'Ra` under-counts this
linked-het cross by −37.6%, which the correction removes.

## PMV (effect uncertainty)

`K` is exactly the progeny **genotypic covariance** `Cov(X_k, X_l)`
(`Var(G)=β'Kβ` for all `β`; diagonal check: `K_kk = ¼(x1_k−x2_k)² +
½(δ1_k²+δ2_k²)` equals the DH dosage variance, e.g. `¾` for `Aa×aa`). So PMV
carries the standard form with this `K` as the transmission covariance:
`PMV = β̂'Kβ̂ + trace(K Σ_β)`, `Σ_β` the posterior effect covariance. The
parent-homozygosity assumption lives entirely in `K`; `Σ_β` is orthogonal.

## RIL-∞ extension (derived and validated)

The RIL (infinite-selfing) case follows by the **same** law-of-total-covariance
argument; only the F1-side co-inheritance changes from `(1−r_kl)` (the single DH
meiosis) to the RIL two-point same-origin probability `(1−R_kl) = 1/(1+2r_kl)`.
Its product with the (single) parental-meiosis kernel `(1−2r_kl)` collapses:
`(1−R_kl)(1−2r_kl) = (1−2r_kl)/(1+2r_kl) = Φ*_kl`. So Term 2 reuses the **same
RIL kernel `Φ*` as Term 1**, with **no** extra `(1−r)` factor:

```
Var(G)_RIL∞ = a' Φ* a  +  1/2 β'[ Φ* ∘ (δ1δ1' + δ2δ2') ] β ,   Φ*_kl = (1−2r_kl)/(1+2r_kl)
```

**Validated** against a pure-Haldane single-seed-descent-to-fixation Monte-Carlo
(single-locus `Aa×aa → ¾`; multi-locus coupling vs repulsion split reproduced;
inbred acid-test recovers `a'Φ*a`). The **only** genuinely open case is
**finite-generation** RILs (F2:F3, F3:F4 …), which additionally carry
generation-dependent residual *progeny* heterozygosity — there the RAH lineage
simulation remains the authoritative oracle.

## Implementation (done)

`ng_gms_additive_var_general(parent1, parent2, haplo_mat, beta, recomb_decay_mat,
beta_cov = NULL, target = c("DH","RIL"))` computes `VPM` and
`PMV = β'Kβ + trace(K Σ_β)` for both targets:

- DH  : Term-2 kernel `(1−r)∘R = ½(1+R)∘R`;  RIL : Term-2 kernel `R` (= `Φ*`).
- Reduces byte-identically to `a'Ra` for inbred parents (δ = 0).
- Validated (`tests/gms_general_var.R`): single-locus closed forms, inbred
  reduction, phase-dependence, PMV, and a pure-Haldane MC (DH −0.24%, RIL +0.49%).

**Wired into `ng_score_crosses`** (and threaded through `ng_design_crosses` /
`ng_run_cross_prediction`): when `parent_type = "ril"` and `phased_haplotypes`
are supplied, het-parent crosses get the exact variance while inbred-parent
crosses keep the a'Ra path byte-identical (`tests/het_parent_correction_e2e.R`).
Requires **phased** parental haplotypes; dosage-only input recovers only the
diagonal/unlinked part of Term 2. Uses a **dense** kernel (the correction is
exact only densely; the a'Ra `window_cm` truncation does not apply to it).
