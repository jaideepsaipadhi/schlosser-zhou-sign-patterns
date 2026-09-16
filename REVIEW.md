# Independent verification and self-review

This file records how the results were re-checked independently of the original code, and what risks remain.
**No human expert has refereed the paper yet.** Everything below is machine verification plus a self-audit.

## 1. Independent re-implementation of the exact computations (PARI/GP 2.15.4)
Directory: `verification/pari/`. This implementation shares no code with the original FLINT/Python scripts and
uses a different algorithm: the original uses Descartes' rule with bisection, while this one uses integer-scaled
polynomials $n!\,c_n(\delta)$ and Sturm root counting (`polsturm`). Endpoint zeros are divided out exactly.

| check | PARI/GP result | agrees with original |
|---|---|---|
| Q8 on [8/3,4], n ≤ 700 | no violations | ✓ |
| Q8 on [β⁺, 8/3], n ≤ 700 | no violations | ✓ |
| Q8 on [−1, −0.772], n ≤ 400 (c₄ divided by δ²−7δ−6) | no violations | ✓ |
| Q12 on [2,3], n ≤ 400 | no violations | ✓ |
| Q12 on [−1, δ₁⁺], n ≤ 780 | only n=22, one root, correct sign at −1 | ✓ |
| G11 on [γ,2], n ≤ 310 | no violations; the degree-18 factor of c₂₁ has exactly one root there | ✓ |
| G7 on [2,3], n ≤ 897 | no violations | ✓ |
| G7 on [3, δ_lo], n ≤ 897 | no violations | ✓ |
| G7 on [δ_lo, δ_hi], n ≤ 897 | only n=897: one root, signs −/+; no root on [δ_hi, δ*⁺] | ✓ |
| Cusp typing, Q8, k ≤ 40 (490 cusps) | 0 mismatches | ✓ |
| Cusp typing, Q12, k ≤ 36 (396 cusps) | 0 mismatches | ✓ |

## 2. Audit of the analytic bounds against exact data
Directory: `verification/audit/`. Each analytic bound used in the proofs was compared with exact coefficients
(or exact derivatives). A mistake in a bound would show up as a ratio above 1.

| bound tested | worst ratio (must be ≤ 1) |
|---|---|
| Farey/Ford chord lemma: chord length, Re z, arc length, Re(1/z) ≥ 1 (all fractions, N ≤ 59 and N = 97, 150, 211) | 0.85, 1.00 (equality case), 0.75; min Re(1/z) = 1.000000 |
| Q8 expansion + error bound (δ = 2.7, 3.2, 3.9; 60 ≤ n ≤ 420) | 0.078 |
| G7 expansion + error bound (δ = 2.5 … 4.8735) | 0.028 |
| Q12 expansion + error bound (δ = 2.2 … 3) | 1.6·10⁻⁴ |
| G11 expansion + error bound (δ = 1.76 … 2) | 0.010 |
| G7 derivative error bound (window near δ_c) | 3.4·10⁻⁶ |
| Q8 ε-difference bound (near δ = −1) | 1.1·10⁻⁵ |

## 3. The hypothesis behind "one check at n₁ covers all n ≥ n₁"
Directory: `verification/monotone/`. Every competitor/dominant ratio has the form C·y^a·e^{−g y}, which is
decreasing for y ≥ a/g. The script checks, for every dominance and derivative checker, that g > 0 and
y(n₁) ≥ a/g, with the conservative value a = 5. All cases pass (see `n1_monotone_check.log`); the smallest
margin is y(n₁)/(a/g) ≈ 1.6, for G7 at δ = 5 in its correction-dominated residues.
The one term with g = 0 (the |a₇|·D·𝒫₇ term in the G7 window derivative) is handled separately. It only
matters for 898 ≤ n < 2500, and its value at n = 2500 bounds it there.

## 4. Corrections found during review
* **Missing harmonic factor.** The G7 arc-integral derivative bound omitted a harmonic factor
  Σ 1/k ≤ (1/7)(1+log(N/7)). It has been added, and every affected certificate was rerun; margins are
  essentially unchanged.
* **Chord length.** The chord length bound is 2√2k/(N+1), not the 2k/(N+1) quoted in the literature.
  The constants use the corrected value.
* **Transformation convention.** Schlosser–Zhou's Prop. 7 must use hh′ ≡ −1 (mod k). This is certified in
  `code/conj20/gp_branch_cert.py`.

## 5. Residual risks
1. **Unrefereed analytic arguments.** An error could survive in a step that holds numerically but is argued
   wrongly in the text. The sections to read most carefully are §7 (Q8 near δ = −1) and §§11–13
   (Conjectures 20, 23, 24).
2. **Classical input.** The eta transformation formula is used as a classical input; only its branch and
   sign convention were certified.
3. **Software.** All certified computations rely on FLINT/Arb and PARI/GP.
