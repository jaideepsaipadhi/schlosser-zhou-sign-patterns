#!/usr/bin/env bash
# Reproduces every certified computation in the paper. Logs go to results/.
# Approximate single-core runtimes are given; total ~20-25 min (excluding the optional direct runs).
set -euo pipefail
cd "$(dirname "$0")/code"
R=../results
run(){ name=$1; shift; echo ">> $name: $*"; "$@" > "$R/$name.log" 2>&1; tail -3 "$R/$name.log"; }
# --- Section 3: cusp structure (exact) -----------------------------------------
run weil_check        python3 weil_check.py                    # Lemma 3.1 sanity check
run orbit             python3 orbit.py                         # Prop 3.2 (writes orbit.pkl)
run cusp_types        python3 cusp_types.py                    # Prop 3.3
run cusp_exact        python3 cusp_exact.py 40                 # Lemma 3.5   (~1 min)
run cusp_exact_inv    python3 cusp_exact_inv.py 40             # Sec 7 analogue (~1 min)
run exact             python3 exact.py                         # Lemma 3.7 Gauss sums
run branch_cert       python3 branch_cert.py                   # Lemma 3.8
run supK              python3 supK.py                          # Prop 3.6, K_NG
run supKinv           python3 supKinv.py                       # Sec 7, K'_NG
run inv_exact         python3 inv_exact.py                     # Sec 7 exact phases/branches
# --- Theorem A: delta in [8/3,4] -----------------------------------------------
run consts            python3 consts.py
run c3_700            python3 c3.py 700                        # Sec 5.5, n<=700 exact (~1-2 min)
run certify_mid       python3 certify_mid.py 701 20000 16      # Sec 5.3 (~2 min)
run tail_sc           python3 tail_sc.py                       # Sec 5.4
# --- Theorem B: disproof on [beta, 8/3) ----------------------------------------
run discert_sc        python3 discert_sc.py                    # Sec 6 explicit counterexample
run c3_betaplus_700   python3 c3_betaplus.py 700               # Remark, n<=700 on [beta+,8/3] (~4 min)
# --- Theorem C: parts (a), (c), (d) --------------------------------------------
run neg1              python3 neg1.py                          # vanishing/derivative data at delta=-1
run neg_finite2_400   python3 neg_finite2.py 400               # part (c), n<=400 exact on [-1,-0.772] (~2 min)
run neg_certify_772   python3 neg_certify.py 372 20000 0.772   # part (c), 372<=n<=20000
run neg_tail772       python3 neg_tail772.py                   # part (c), n>20000
run pm2_finite        python3 pm2_finite.py                    # parts (a),(d), n<=400 exact
run pm2_certify       python3 pm2_certify.py 297 20000         # parts (a),(d), 297<=n<=20000
run pm2_tail          python3 pm2_tail.py                      # parts (a),(d), n>20000
echo "All certified checks completed."
# --- Optional: direct, analysis-free disproof (long; ~30 min, <1 GB RAM) --------
#   python3 q8_lite.py 860000 2.66448 2600 400 8192 > ../results/q8_lite_run_860000.txt
