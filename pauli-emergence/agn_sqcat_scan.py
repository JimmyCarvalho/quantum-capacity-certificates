"""Squeezed-cat inputs on the AGN channel: residues (r0, r1) mod k of D(alpha) S(r)|0>, weights w, 1 - w.
Continuation in sigma^2 from s0 to s1; prints the optimal (alpha, r, w) and I_c at each point.
usage: python3 agn_sqcat_scan.py k r0,r1 s0 s1 step alpha0 r0sq w0 tag"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import json
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "depolarizing-threshold"))
import numpy as np
from scipy.optimize import minimize

from agn_sqcat import squeezed_coherent
from gausschan import EmergentChannel, charge_blocks, coherent_info
from limit_codes import emergent_TG

NIN, D = 200, 400


def code(k, res, alpha, r, w):
    w = min(max(w, 1e-6), 1 - 1e-6)
    psi = squeezed_coherent(alpha, r)[: NIN + 1]
    n = np.arange(NIN + 1)
    cs = []
    for x in res:
        c = np.where(n % k == x, psi, 0.0)
        cs.append(c / np.linalg.norm(c))
    return np.array([np.sqrt(w) * cs[0], np.sqrt(1 - w) * cs[1]])


if __name__ == "__main__":
    k, res = int(sys.argv[1]), [int(x) for x in sys.argv[2].split(",")]
    s0, s1, step = float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
    x = np.array([float(sys.argv[6]), float(sys.argv[7]), float(sys.argv[8])])
    tag = sys.argv[9]
    bl = charge_blocks(res, k, D)
    s2 = s0
    while s2 <= s1 + 1e-12:
        T, G = emergent_TG("agn", s2)
        ch = EmergentChannel(0.5, NIN, D, T=T, G=G)
        f = lambda y: -coherent_info(ch, code(k, res, *y), blocks=bl)
        best = None
        for start in [x, x + [0.3, -0.05, 0.0], x + [-0.3, 0.05, 0.0]]:
            out = minimize(f, start, method="Nelder-Mead", options={"xatol": 1e-7, "fatol": 1e-14, "maxiter": 1500})
            if best is None or out.fun < best.fun:
                best = out
        x = best.x
        psi = squeezed_coherent(x[0], x[1])
        nbar = float(np.sum(np.arange(len(psi)) * psi ** 2))
        print(json.dumps({"tag": tag, "sigma2": round(s2, 5), "Ic": -best.fun, "alpha": round(x[0], 5), "r": round(x[1], 5),
                          "w0": round(x[2], 5), "nbar": round(nbar, 2)}), flush=True)
        if -best.fun <= 1e-10:
            break
        s2 = round(s2 + step, 6)
