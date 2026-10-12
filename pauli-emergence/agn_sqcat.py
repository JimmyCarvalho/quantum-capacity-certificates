"""Three-component squeezed-cat inputs on the AGN channel: the residue-0 and residue-1 (mod 3) parts of a
displaced squeezed state D(alpha) S(r) |0>, with weights w and 1 - w. Parameters (alpha, r, w); r < 0 is
phase squeezing (super-Poissonian photon statistics). usage: python3 agn_sqcat.py 0.38,0.39,0.40"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "depolarizing-threshold"))
import numpy as np
from scipy.optimize import minimize
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply

from gausschan import EmergentChannel, charge_blocks, coherent_info
from limit_codes import emergent_TG

NBIG, NIN, D = 500, 140, 320
a = diags(np.sqrt(np.arange(1, NBIG)), 1).tocsc()
ad = a.T.tocsc()
vac = np.zeros(NBIG); vac[0] = 1.0


def squeezed_coherent(alpha, r):
    psi = expm_multiply((0.5 * r * (a @ a - ad @ ad)).tocsc(), vac)   # S(r)|0>, r real
    psi = expm_multiply((alpha * (ad - a)).tocsc(), psi)             # D(alpha), alpha real
    return psi


def code(alpha, r, w):
    w = min(max(w, 1e-6), 1 - 1e-6)
    psi = squeezed_coherent(alpha, r)[: NIN + 1]
    n = np.arange(NIN + 1)
    c0 = np.where(n % 3 == 0, psi, 0.0); c1 = np.where(n % 3 == 1, psi, 0.0)
    c0 /= np.linalg.norm(c0); c1 /= np.linalg.norm(c1)
    return np.array([np.sqrt(w) * c0, np.sqrt(1 - w) * c1])


if __name__ == "__main__":
    bl = charge_blocks([0, 1], 3, D)
    for s2 in [float(x) for x in sys.argv[1].split(",")]:
        T, G = emergent_TG("agn", s2)
        ch = EmergentChannel(0.5, NIN, D, T=T, G=G)
        f = lambda x: -coherent_info(ch, code(x[0], x[1], x[2]), blocks=bl)
        best = None
        for a0 in [3.0, 4.0, 5.0]:
            for r0 in [-0.8, -0.4, 0.4]:
                res = minimize(f, [a0, r0, 0.44], method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-13, "maxiter": 600})
                if best is None or res.fun < best.fun:
                    best = res
        al, r, w = best.x
        print(f"{s2}: best squeezed cat Ic = {-best.fun:.4e} at alpha = {al:.4f}, r = {r:.4f}, w0 = {w:.4f}", flush=True)
        np.save(f"runs/agn_sqcat_{s2}.npy", code(al, r, w))
