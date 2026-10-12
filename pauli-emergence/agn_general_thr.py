"""AGN (or other emergent) channel: continuation in the noise parameter (sigma^2 for agn, p otherwise) for symmetry-free rank-2 inputs (no residue constraint).
Saves the optimised state at every sigma^2 where I_c > thr; stops at the first failure."""
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
ap.add_argument("--channel", default="agn")
ap.add_argument("--s0", type=float, required=True)
ap.add_argument("--step", type=float, default=0.001)
ap.add_argument("--smax", type=float, default=0.45)
ap.add_argument("--thr", type=float, default=2e-8)
ap.add_argument("--Nin", type=int, default=200)
ap.add_argument("--D", type=int, default=400)
ap.add_argument("--maxiter", type=int, default=600)
ap.add_argument("--noise", type=float, default=0.0)
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--tag", required=True)
a = ap.parse_args()

C0 = np.load(a.init)
r = C0.shape[0]
C = np.zeros((r, a.Nin + 1))
C[:, : min(C0.shape[1], a.Nin + 1)] = C0[:, : a.Nin + 1]
C /= np.linalg.norm(C)
if a.noise > 0:  # break the residue symmetry of the starting point
    # envelope of the state itself (neighbours included, so off-residue levels get weight too)
    env = np.abs(C).max(axis=0)
    env = np.maximum(env, np.maximum(np.roll(env, 1), np.roll(env, -1)))
    C += a.noise * np.random.default_rng(a.seed).normal(size=C.shape) * env
    C /= np.linalg.norm(C)
sup = supports([0] * r, 1, a.Nin)
s2 = a.s0
while s2 <= a.smax + 1e-12:
    t = time.time()
    T, G = emergent_TG(a.channel, s2)
    ch = EmergentChannel(0.5, a.Nin, a.D, T=T, G=G)
    Ic0 = coherent_info(ch, C)
    Ic, C2 = Ic0, C
    for _ in range(3):
        Ic, C2 = optimise(ch, sup, C2, maxiter=a.maxiter, verbose=False, scale=1.0 / max(abs(Ic), 1e-10))
    n = np.arange(a.Nin + 1)
    print(json.dumps({"sigma2": round(s2, 7), "Ic": Ic, "Ic_start": Ic0, "weights": np.sum(C2**2, 1).round(4).tolist(),
                      "mean_n": round(float(np.sum(n * C2**2)), 2), "tail": float(np.sum(C2[:, -20:] ** 2)),
                      "sec": round(time.time() - t)}), flush=True)
    if Ic <= a.thr:
        break
    C = C2
    np.save(f"runs/{a.tag}_s{round(s2, 7)}.npy", C)
    s2 = round(s2 + a.step, 7)
