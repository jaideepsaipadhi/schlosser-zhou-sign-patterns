# Sign patterns of real powers of infinite products: resolution of Conjectures 20, 21, 23, 24 of Schlosser–Zhou

**Paper:** `paper/sz_joint.pdf` (source `paper/sz_joint.tex`). Earlier split drafts are in `archive/`.

| Conjecture | product | result |
|---|---|---|
| 20 | $G_7$ | true for δ=1, 2≤δ≤3, 3<δ≤δ_c, δ=5; **false on (δ_c,5)**, δ_c=4.8735075853867342634… (root of c(897)) |
| 21 | $Q_8$ (Göllnitz–Gordon) | true for δ=±2, −1<δ≤(7−√73)/2, 8/3≤δ≤4; **false on [β,8/3)** |
| 23 | $G_{11}$ | **true** |
| 24 | $Q_{12}$ | true for δ=±1, 2≤δ≤3, −1<δ≤δ_1; **false on (δ_1,0)**, δ_1=−0.6441134808… (root of c(22)) |

Reproduce: `run_all.sh` (Conjecture 21) and `run_all_rest.sh` (Conjectures 20, 23, 24).

---
Paper and certified computations for a complete resolution of Conjecture 21 of
M. J. Schlosser and N. H. Zhou, *On the infinite Borwein product raised to a positive real power*,
Ramanujan J. **61** (2023), 515–543, Appendix 1.

Let $Q_8(q)=(q,q^7;q^8)_\infty/(q^3,q^5;q^8)_\infty$ and $Q_8(q)^\delta=\sum_n c_\delta(n)q^n$.

| part of Conjecture 21 | result |
|---|---|
| (a) δ = 2, length-16 pattern `+-++-+--+--+-++-` | **true** (certified proof; earlier proof in He–Li) |
| (b) β ≤ δ ≤ 4, pattern `+-++-+--` | **true on [8/3, 4]**, **false on [β, 8/3)**; 8/3 is sharp |
| (c) −1 < δ ≤ (7−√73)/2, pattern `+++----+` | **true** (the range −1 < δ < −0.99 was open) |
| (d) δ = −2, length-16 pattern `+++++----+++----` | **true** (certified proof; earlier proof in He–Li) |

A zero coefficient is counted as compatible with either sign. The zeros in the ranges are c₄(3), c₋₂(8), and c_b(4) at b = (7−√73)/2.

An explicit counterexample to (b) is c_{2.66448}(1 000 010) < −5.7·10³⁸⁷. It is certified twice, independently:
- analytically, in §6 of the paper;
- by direct power-series computation (`results/q8_lite_run_860000.txt`). This run finds 12 499 sign violations for 760 014 ≤ n ≤ 860 000.

## Repository layout
```
paper/     sz_joint.tex, sz_joint.pdf (joint paper); archive/ has the earlier split drafts
code/      all verification scripts (Python 3, python-flint / Arb)
results/   logs of the certified runs quoted in the paper
run_all.sh reproduces every certified computation (~20–25 min, single core)
```

## Requirements
```
pip install -r requirements.txt   # python-flint 0.9.0 (FLINT 3.6.0), mpmath, numpy
```
The runs quoted in the paper used Python 3.12.3, python-flint 0.9.0, and FLINT 3.6.0.

## What each script certifies
| paper | script | kind |
|---|---|---|
| Lemma 3.1 (sanity check) | `weil_check.py` | floating point, not used in proofs |
| Prop. 3.2–3.3 (orbit, cusp types) | `orbit.py`, `cusp_types.py` | exact, ℤ[ζ₁₂₈] |
| Lemma 3.5 (growing cusps, k ≤ 40) | `cusp_exact.py 40`, `cusp_exact_inv.py 40` | exact, mod Φ₁₆ₖ |
| Prop. 3.6 (K_NG), and K′_NG in §7 | `supK.py`, `supKinv.py` | ball arithmetic |
| Lemma 3.7 (Gauss sums) | `exact.py` | exact |
| Lemma 3.8, and the negative-δ phases in §7 | `branch_cert.py`, `inv_exact.py` | ball arithmetic |
| Remark 4.5 (error constant) | `consts.py` | ball arithmetic |
| §5.5, n ≤ 700 on [8/3, 4] | `c3.py 700` | exact polynomials in ℚ[δ] |
| §5.3, 701 ≤ n ≤ 20000 | `certify_mid.py 701 20000 16` | ball arithmetic |
| §5.4, tail | `tail_sc.py` | ball arithmetic |
| Theorem 1.3 / §6, explicit counterexample | `discert_sc.py` | ball arithmetic |
| Remark on [β⁺, 8/3], n ≤ 700 | `c3_betaplus.py 700` | exact |
| §7, part (c) | `neg1.py`, `neg_finite2.py 400`, `neg_certify.py 372 20000 0.772`, `neg_tail772.py` | exact / ball arithmetic |
| §7, parts (a), (d) | `pm2_finite.py`, `pm2_certify.py 297 20000`, `pm2_tail.py` | exact / ball arithmetic |
| independent disproof | `q8_lite.py 860000 2.66448 2600 400 8192` (~30 min, <1 GB) | ball arithmetic |
| independent disproof (older, 15–20 GB) | `q8_signcheck.py` | ball arithmetic |

Scripts `neg_certify.py 272 20000 0.99`, `neg_finite.py 300` and `neg_tail.py` handle the sub-range −1 < δ ≤ −0.99 alone. They are subsumed by the 0.772 runs.

## Trust boundary
Every certified step uses either exact integer/rational arithmetic or Arb ball arithmetic. Arb returns intervals guaranteed to contain the true value, and an inequality is accepted only when it holds on the whole interval. Floating point is used only for exploration and for dictionary keys, and every floating-point match is confirmed exactly.

The mathematical inputs from the literature are classical:
- Poisson summation;
- the fact that S and T generate SL₂(ℤ);
- Ford circles and Farey neighbours;
- integral representations of I₁;
- the Andrews–Bressoud vanishing theorem (1979), used in §7.

## Related work
- He–Li, arXiv:2509.10023, proved parts (a), (d), and (c) for δ ≥ −0.99.
- Richmond–Szekeres (1978) treated δ = 1.


## Independent verification
See `REVIEW.md`. The directory `verification/` contains:
* `pari/`: an independent PARI/GP re-implementation of the exact finite checks and of the cusp typing;
* `audit/`: checks of every analytic bound against exact data;
* `monotone/`: a programmatic check of the monotonicity hypothesis behind "one check at n₁ covers all n ≥ n₁".

## License
Released under the MIT License (see `LICENSE`).

## Conjectures 20, 23, 24 (sections 8–12 of the joint paper)
| Conjecture | result |
|---|---|
| 20 ($G_7$) | true for $\delta=1$, $2\le\delta\le3$, $3<\delta\le\delta_c$, $\delta=5$; **false on $(\delta_c,5)$**, $\delta_c=4.87350758538673426336\ldots$ (root of $c(897)$) |
| 23 ($G_{11}$) | **true** |
| 24 ($Q_{12}$) | true for $\delta=\pm1$, $2\le\delta\le3$, $-1<\delta\le\delta_1$; **false on $(\delta_1,0)$**, $\delta_1=-0.6441134808134765\ldots$ (root of $c(22)$); e.g. $c_{-1/4}(6)=545/65536$ |

Reproduce with `run_all_rest.sh`; logs are in `results/conj20`, `results/conj23`, `results/conj24`.
All certificates cited in the paper have logs in `results/`, including `results/conj24/supK12_plus.log` ($K_{12}\le550.7490$), `supK12_minus.log` ($K'_{12}\le547.2229$) and `results/conj20/pointcheck_7000.log` (no violation for n≤7000 at $\delta_{lo}$; only n=897 at $\delta_{hi}$).
The $G_7$ exact-polynomial pickle (~330 MB) is regenerated by `c897.py` and is not included.

Note: the transformation of Schlosser–Zhou [Prop. 7] must be used with $hh'\equiv-1\pmod k$; `gp_branch_cert.py` certifies this and the branch at all needed cusps.
