#!/usr/bin/env bash
# Reproduce the certificates for Conjectures 20, 23, 24 (paper2). Logs -> results/conj2x/. Approx. total 1.5-2 h, one core.
set -euo pipefail; cd "$(dirname "$0")/code"
# ---- Conjecture 24 (Q12) ----
cd conj24
python3 orbit12.py && python3 types12.py                   # orbit (48) and cusp types (exact, Z[zeta_192])
python3 q12_cusps.py 36                                    # exact typing of cusps k<=36
python3 q12_exact.py      > ../../results/conj24/exact_phases.log   # Gauss sums + certified branches
python3 supK12.py 1 > ../../results/conj24/supK12_plus.log; python3 supK12.py -1 > ../../results/conj24/supK12_minus.log                 # K12=550.75, K12'=547.23 (~5 min each)
python3 q12_finite.py 400 2/1 3/1 '+-+--+-+-++-'  > ../../results/conj24/finite_pos_400.log
python3 q12_pos.py 400 11000 8 > ../../results/conj24/pos_400_11000.log
python3 q12_pos.py 11001 20000 8 > ../../results/conj24/pos_11001_20000.log
python3 q12_pos_tail.py  > ../../results/conj24/pos_tail.log
python3 q12_pm1.py       > ../../results/conj24/pm1.log
python3 q12_neg_threshold.py 260                            # delta_1 = root of c_22
python3 q12_finite.py 780 -1/1 -6441134808/10000000000 '+++++------+' > ../../results/conj24/finite_neg_780.log
python3 q12_neg.py 774 20000 0.6441134808 > ../../results/conj24/neg_774_20000.log
python3 q12_neg_tail.py  > ../../results/conj24/neg_tail.log
cd ..
# ---- Conjecture 23 (G11) ----
cd conj23
python3 ../conj20/gp_branch_cert.py 11 11,22 > ../../results/conj23/branch_cert_p11.log
python3 g11_finite.py 310 > ../../results/conj23/finite_310.log
python3 g11_dom.py 1.7584 2 300 '+--+++-+--+' > ../../results/conj23/dom_gamma_2.log
python3 g11_point.py 600 > ../../results/conj23/points_1_3.log
cd ..
# ---- Conjecture 20 (G7) ----
cd conj20
python3 gp_branch_cert.py 7 7,14,21,28 > ../../results/conj20/branch_cert_p7.log
python3 c897.py                                            # exact c_n, n<=897 (pickle ~330 MB), root isolation of c_897
for r in "0 500" "501 600" "601 680" "681 740" "741 790" "791 830" "831 865" "866 897"; do
  python3 g7_finite.py $r 3/1 487350758538673426330055236816/100000000000000000000000000000 "+-++---"; done > ../../results/conj20/finite_upper.log
for r in "0 500" "501 600" "601 680" "681 740" "741 790" "791 830" "831 865" "866 897"; do
  python3 g7_finite.py $r 487350758538673426330055236816/100000000000000000000000000000 487350758538673426343362426758/100000000000000000000000000000 "+-++---"; done > ../../results/conj20/finite_bracket.log
for r in "0 500" "501 680" "681 790" "791 865" "866 897"; do python3 g7_finite2.py $r 2/1 3/1 "+--+++-"; done > ../../results/conj20/finite_23.log
python3 g7_tile.py 2 2.999999999 898 '+--+++-' 1e-12 > ../../results/conj20/tile_2_3.log
python3 g7_tile.py 3.000000001 4.87 898 '+-++---' 1e-12 > ../../results/conj20/tile_3_487.log
python3 g7_tile.py 3 3.000000001 898 '+-++---' 1e-12 0136 > ../../results/conj20/tile_near3_up.log
python3 g7_tile.py 2.999999999 3 898 '+--+++-' 1e-12 0136 > ../../results/conj20/tile_near3_down.log
python3 g7_tile.py 4.87 4.8735075853 898 '+-++---' 1e-14 > ../../results/conj20/tile_487_a.log
python3 g7_tile.py 4.8735075853 4.87350758538 898 '+-++---' 1e-16 > ../../results/conj20/tile_487_b.log
python3 g7_tile.py 4.87350758538 4.8735075853867342634 2500 '+-++---' 1e-22 > ../../results/conj20/window_n2500.log
python3 g7_tile.py 4.87350758538 4.8735075853867342634 898 '+-++---' 1e-20 023456 > ../../results/conj20/window_other_residues.log
python3 g7_deriv.py 4.87350758538 4.87350758538673426343362426758 898 > ../../results/conj20/window_derivative.log
python3 pointcheck.py 7000 3000 > ../../results/conj20/pointcheck_7000.log   # point values n<=7000 at delta_lo/hi
python3 g7_eps3.py  > ../../results/conj20/eps_near3.log
python3 g7_point.py 400 > ../../results/conj20/points_1_3_5.log
echo done
