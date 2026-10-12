"""AGN channel: optimise a rank-r input with no residue constraint, from a residue code plus noise.
Checks whether breaking the Z_k symmetry of the residue codes increases I_c."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "depolarizing-threshold"))
import argparse
import json
import time

import numpy as np

from gausschan import EmergentChannel, coherent_info
from limit_codes import emergent_TG
from optimize import optimise, supports

ap = argparse.ArgumentParser()
ap.add_argument("--init", required=True)
ap.add_argument("--sigma2", type=float, required=True)
ap.add_argument("--rank", type=int, default=2)
ap.add_argument("--noise", type=float, default=0.02)
ap.add_argument("--Nin", type=int, default=140)
ap.add_argument("--D", type=int, default=320)
ap.add_argument("--maxiter", type=int, default=600)
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--tag", required=True)
a = ap.parse_args()

rng = np.random.default_rng(a.seed)
C0 = np.load(a.init)
C = np.zeros((a.rank, a.Nin + 1))
C[: C0.shape[0], : min(C0.shape[1], a.Nin + 1)] = C0[:, : a.Nin + 1]
C = C / np.linalg.norm(C)
C += a.noise * rng.normal(size=C.shape) * np.exp(-((np.arange(a.Nin + 1) - 14) ** 2) / 200)
C /= np.linalg.norm(C)
sup = supports([0] * a.rank, 1, a.Nin)
T, G = emergent_TG("agn", a.sigma2)
ch = EmergentChannel(0.5, a.Nin, a.D, T=T, G=G)
t = time.time()
Ic0 = coherent_info(ch, C)
for _ in range(3):
    Ic, C = optimise(ch, sup, C, maxiter=a.maxiter, verbose=False, scale=1.0 / max(abs(Ic0), 1e-10))
    Ic0 = Ic
np.save(f"runs/{a.tag}.npy", C)
print(json.dumps({"tag": a.tag, "sigma2": a.sigma2, "rank": a.rank, "Ic": Ic, "weights": np.sum(C**2, 1).round(4).tolist(),
                  "sec": round(time.time() - t)}), flush=True)
