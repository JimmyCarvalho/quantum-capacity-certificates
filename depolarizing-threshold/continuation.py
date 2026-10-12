"""Track the optimal code of a given residue family as eta decreases; report max Ic(eta)."""
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse
import json
import time

import numpy as np

from gausschan import EmergentChannel, charge_blocks, coherent_info
from optimize import gaussian_init, optimise, supports

ap = argparse.ArgumentParser()
ap.add_argument("--k", type=int, default=3)
ap.add_argument("--res", type=str, default="0,1")
ap.add_argument("--etas", type=str, default="0.76,0.75,0.745,0.74,0.735,0.732,0.7294452")
ap.add_argument("--Nin", type=int, default=200)
ap.add_argument("--D", type=int, default=300)
ap.add_argument("--mean", type=float, default=20.0)
ap.add_argument("--width", type=float, default=2.0)
ap.add_argument("--q", type=float, default=0.45)
ap.add_argument("--maxiter", type=int, default=600)
ap.add_argument("--init", type=str, default=None)
ap.add_argument("--tag", type=str, required=True)
a = ap.parse_args()
residues = [int(r) for r in a.res.split(",")]
sup = supports(residues, a.k, a.Nin)
blocks = charge_blocks(residues, a.k, a.D)
if a.init:
    C = np.load(a.init)
    C = C[:, : a.Nin + 1] if C.shape[1] > a.Nin else np.pad(C, ((0, 0), (0, a.Nin + 1 - C.shape[1])))
else:
    r = len(residues)
    C = gaussian_init(sup, a.Nin, a.mean, a.width, weights=[a.q] + [(1 - a.q) / (r - 1)] * (r - 1))
n = np.arange(a.Nin + 1)
for eta in [float(e) for e in a.etas.split(",")]:
    t = time.time()
    ch = EmergentChannel(eta, a.Nin, a.D)
    for _ in range(2):
        Ic0 = coherent_info(ch, C, blocks=blocks)
        Ic, C = optimise(ch, sup, C, maxiter=a.maxiter, verbose=False,
                         scale=1.0 / max(abs(Ic0), 1e-12), blocks=blocks)
    np.save(f"runs/{a.tag}_eta{eta}.npy", C)
    print(json.dumps({"eta": eta, "p": 0.75 * (1 - eta), "Ic": Ic, "weights": np.sum(C**2, 1).round(6).tolist(),
                      "mean_n": round(float(np.sum(n * C**2)), 3), "tail": float(np.sum(C[:, -20:] ** 2)),
                      "sec": round(time.time() - t)}), flush=True)
