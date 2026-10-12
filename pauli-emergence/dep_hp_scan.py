"""Depolarizing record regime with high-precision I_c (hp_ic): squeezed-cat residue codes
(residues 0, 1 mod 3 of D(alpha) S(r)|0>, weights w, 1 - w) at a given eta.
usage: python3 dep_hp_scan.py grid eta alphas rs w   |   python3 dep_hp_scan.py opt eta alpha r w tag"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import json
import sys
import time
from fractions import Fraction
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "depolarizing-threshold"))
import flint
import numpy as np
from scipy.optimize import minimize

import hp_ic
from agn_sqcat import squeezed_coherent
from certify import ArbChannel

NIN, D = int(os.environ.get("NIN", 220)), int(os.environ.get("DOUT", 320))
flint.ctx.prec = 128


def code(alpha, r, w):
    w = min(max(w, 1e-6), 1 - 1e-6)
    psi = squeezed_coherent(alpha, r)[: NIN + 1]
    n = np.arange(NIN + 1)
    cs = []
    for x in (0, 1):
        c = np.where(n % 3 == x, psi, 0.0)
        cs.append(c / np.linalg.norm(c))
    return np.array([np.sqrt(w) * cs[0], np.sqrt(1 - w) * cs[1]])


if __name__ == "__main__":
    mode, eta = sys.argv[1], Fraction(sys.argv[2])
    ch = ArbChannel(eta, NIN, D)
    f = lambda y: float(hp_ic.ic(ch, code(*y)).mid())
    if mode == "grid":
        w = float(sys.argv[5])
        for al in [float(x) for x in sys.argv[3].split(",")]:
            for r in [float(x) for x in sys.argv[4].split(",")]:
                t = time.time()
                print(json.dumps({"eta": str(eta), "alpha": al, "r": r, "w": w, "Ic": f((al, r, w)), "sec": round(time.time() - t)}), flush=True)
    else:
        x0 = [float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])]
        tag = sys.argv[6]
        best = [None, -1.0]

        def obj(y):
            v = f(y)
            if v > best[1]:
                best[0], best[1] = list(y), v
                np.save(f"runs/{tag}.npy", code(*y))
                print(json.dumps({"eta": str(eta), "alpha": y[0], "r": y[1], "w": y[2], "Ic": v}), flush=True)
            return -v * 1e12

        minimize(obj, x0, method="Nelder-Mead", options={"initial_simplex": [x0, [x0[0] + 0.3, x0[1], x0[2]], [x0[0], x0[1] + 0.1, x0[2]], [x0[0], x0[1], x0[2] + 0.02]],
                                                        "xatol": 1e-4, "fatol": 1e-4, "maxiter": 200})
        print(json.dumps({"RESULT": tag, "best": best}), flush=True)
