# Reproducing the certificates

Tested on Linux, Python 3.11+, single core.

## 1. Install

    python3 -m venv venv
    source venv/bin/activate
    pip install -r requirements.txt

(`python-flint==0.9.0`, `mpmath>=1.3`, `numpy>=1.24`.)

## 2. Run

    chmod +x run_all.sh run_all_rest.sh
    ./run_all.sh        # Q8: Theorems A, B, C.   ~20-25 min
    ./run_all_rest.sh   # Conjectures 20, 23, 24. ~1.5-2 h

Run them from the repo root, not from `code/`. Logs are written to
`results/` and `results/conj2{0,3,4}/`.

Use a machine that will not be interrupted: `run_all_rest.sh` is a
single long serial job with no checkpointing. To survive a closed
terminal:

    nohup ./run_all_rest.sh > rest_stdout.log 2>&1 &

## 3. Disk

`code/conj20/c897.py` writes a pickle of roughly 330 MB. Everything
else is small. Budget ~1 GB free.

## 4. What success looks like

`run_all.sh` ends with `All certified checks completed.`
`run_all_rest.sh` ends with `done`.
No log should contain a Python traceback. Note that several logs
legitimately contain the strings "failures 0" and "relative error" --
grepping for "error" or "fail" gives false positives.

## 5. Values to check against the paper

These were confirmed on a clean checkout (no cached .pkl):

| quantity                        | expected              | log |
|---------------------------------|-----------------------|-----|
| sup on [8/3,4]                  | 0.375                 | tail_sc.log |
| kappa(x0)                       | 0.99412               | tail_sc.log |
| c_{2.66448}(1000010)            | < -5.7393e387         | discert_sc.log |
| exact zeros delta=4, n<=700     | [3]  (c_4(3)=0 only)  | c3_700.log |
| exact zeros delta=-2, n<=400    | [8]  (c_-2(8)=0 only) | pm2_finite.log |
| K_12 (sgn +1)                   | 550.748906            | conj24/supK12_plus.log |
| K_12 (sgn -1)                   | 547.222822            | conj24/supK12_minus.log |

Not yet reproduced from a clean checkout (the container running them
was restarted mid-job): everything in `run_all_rest.sh` from
`q12_pos.py` onward, i.e. most of Conjectures 24, 23 and 20.

The most delicate single number in the paper is delta_c, isolated to a
bracket of width 1.3e-19 by `code/conj20/g7_tile.py` and `c897.py`.
Check `results/conj20/window_n2500.log` and `window_derivative.log`
against Theorem 1.8 before relying on it.
