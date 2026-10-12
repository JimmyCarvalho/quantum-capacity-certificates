"""L-BFGS on residue codes with high-precision I_c and gradient (hp_ic.ic_and_grad).
For the near-classical plateau of the depolarizing record (|I_c| ~ 1e-12), where double-precision line
searches stall. usage: python3 hp_opt.py init.npy eta Nin D maxiter tag [scale] [channel]
channel: depol (default; eta = shrinking factor) or th<nu>/agn as in certify_limit (eta = p there)."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import json
import sys
import time
from fractions import Fraction

import flint
import numpy as np
from scipy.optimize import minimize

import hp_ic
from certify import ArbChannel
from optimize import pack, supports, unpack

flint.ctx.prec = 128

if __name__ == "__main__":
    path, eta, Nin, D, maxiter, tag = sys.argv[1], Fraction(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
    scale = float(sys.argv[7]) if len(sys.argv) > 7 else None
    ch = ArbChannel(eta, Nin, D)
    C0 = np.zeros((2, Nin + 1))
    S = np.load(path)
    C0[:, : min(Nin + 1, S.shape[1])] = S[:, : Nin + 1]
    sup = supports([0, 1], 3, Nin)
    C0 = unpack(pack(C0, sup), sup, Nin)
    C0 /= np.linalg.norm(C0)
    t0 = time.time()
    best = {"Ic": -np.inf}
    n = [0]

    def f(x):
        C = unpack(x, sup, Nin)
        Ic, g = hp_ic.ic_and_grad(ch, C)
        n[0] += 1
        if Ic > best["Ic"]:
            best["Ic"] = Ic
            np.save(f"runs/{tag}.npy", C / np.linalg.norm(C))
        print(json.dumps({"tag": tag, "eval": n[0], "Ic": Ic, "best": best["Ic"], "sec": round(time.time() - t0)}), flush=True)
        s = scale if scale else 1.0 / max(abs(Ic), 1e-30)
        if scale is None:
            f.scale = getattr(f, "scale", s)
            s = f.scale
        return -s * Ic, -s * pack(g, sup)

    res = minimize(f, pack(C0, sup), jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "maxcor": 30, "ftol": 1e-15, "gtol": 1e-12})
    print(json.dumps({"RESULT": tag, "Ic": best["Ic"], "message": str(res.message), "nit": res.nit}), flush=True)
