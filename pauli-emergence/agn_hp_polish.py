"""L-BFGS polish of a mod-3 residue code on the AGN channel with the coherent information and its gradient
in 128-bit ball arithmetic (hp_par, six charge blocks in parallel). Needed once I_c ~ 1e-10: there the
double-precision gradient is dominated by rounding in the logarithms of near-zero output eigenvalues and
agn_polish.py stalls. Saves the best state so far to out.npy after every improvement.
usage: python3 agn_hp_polish.py sigma2 init.npy out.npy Nin D iters [procs]   (sigma2 as a fraction, e.g. 33/80)"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import json
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "depolarizing-threshold"))
import numpy as np
from scipy.optimize import minimize

import hp_par
from optimize import pack, supports, unpack

if __name__ == "__main__":
    s2, path, out = Fraction(sys.argv[1]), sys.argv[2], sys.argv[3]
    Nin, D, iters = int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
    procs = int(sys.argv[7]) if len(sys.argv) > 7 else 4
    par = hp_par.ParIc(procs=procs)
    key = f"agn:{s2}"
    sup = supports([0, 1], 3, Nin)
    S = np.load(path)
    C = np.zeros((2, Nin + 1))
    C[:, : min(Nin + 1, S.shape[1])] = S[:, : Nin + 1]
    C = unpack(pack(C, sup), sup, Nin)
    C /= np.linalg.norm(C)
    t0 = time.time()
    Ic0, _ = par(key, Nin, D, C)
    print(json.dumps({"sigma2": str(s2), "start": Ic0, "sec": round(time.time() - t0)}), flush=True)
    s = 1.0 / max(abs(Ic0), 1e-40)
    best = [Ic0]
    n_eval = [0]

    def f(x):
        X = unpack(x, sup, Nin)
        Ic, g = par(key, Nin, D, X)
        n_eval[0] += 1
        if Ic > best[0]:
            best[0] = Ic
            np.save(out, X / np.linalg.norm(X))
        print(json.dumps({"eval": n_eval[0], "Ic": Ic, "best": best[0], "sec": round(time.time() - t0)}), flush=True)
        return -s * Ic, -s * pack(g, sup)

    minimize(f, pack(C, sup), jac=True, method="L-BFGS-B",
             options={"maxiter": iters, "maxcor": 30, "ftol": 1e-15, "gtol": 1e-14})
    print(json.dumps({"RESULT": str(s2), "best": best[0], "evals": n_eval[0]}), flush=True)
