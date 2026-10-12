"""Gaussian-limit (emergent channel) threshold of residue codes: continuation in p then bisection."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "depolarizing-threshold"))
import argparse
import json
import time

import numpy as np

from gausschan import EmergentChannel, charge_blocks, coherent_info
from optimize import gaussian_init, optimise, supports
from limit_codes import emergent_TG

ap = argparse.ArgumentParser()
ap.add_argument("--channel", default="bb84")
ap.add_argument("--k", type=int, default=4)
ap.add_argument("--res", default="0,2")
ap.add_argument("--Nin", type=int, default=200)
ap.add_argument("--D", type=int, default=360)
ap.add_argument("--p0", type=float, required=True)
ap.add_argument("--step", type=float, default=0.0005)
ap.add_argument("--tol", type=float, default=0.0001)
ap.add_argument("--init", default=None)
ap.add_argument("--mean", type=float, default=30.0)
ap.add_argument("--width", type=float, default=5.0)
ap.add_argument("--maxiter", type=int, default=500)
ap.add_argument("--thr", type=float, default=1e-6)
ap.add_argument("--tag", required=True)
a = ap.parse_args()
res = [int(x) for x in a.res.split(",")]
sup = supports(res, a.k, a.Nin)
blocks = charge_blocks(res, a.k, a.D)
r = len(res)
C = gaussian_init(sup, a.Nin, a.mean, a.width, weights=[1.0 / r] * r)
if a.init:
    C0 = np.load(a.init)
    C = np.zeros((r, a.Nin + 1))
    m = min(C0.shape[1], a.Nin + 1)
    C[:, :m] = C0[:, :m]
n = np.arange(a.Nin + 1)


def run(p, C):
    t = time.time()
    T, G = emergent_TG(a.channel, p)
    ch = EmergentChannel(0.5, a.Nin, a.D, T=T, G=G)
    for _ in range(2):
        Ic0 = coherent_info(ch, C, blocks=blocks)
        Ic, C = optimise(ch, sup, C, maxiter=a.maxiter, verbose=False, scale=1.0 / max(abs(Ic0), 1e-10), blocks=blocks)
    print(json.dumps({"p": p, "Ic": Ic, "weights": np.sum(C**2, 1).round(4).tolist(),
                      "mean_n": round(float(np.sum(n * C**2)), 2), "tail": float(np.sum(C[:, -10:] ** 2)),
                      "sec": round(time.time() - t)}), flush=True)
    return Ic, C


p, lo, hi, best = a.p0, None, None, None
while hi is None:
    Ic, C2 = run(p, C)
    if Ic > a.thr and min(np.sum(C2**2, 1)) > 0.05:
        lo, C, best = p, C2, Ic
        np.save(f"runs/{a.tag}_best.npy", C)
        p = round(p + a.step, 7)
    else:
        hi = p
        if lo is None:
            raise SystemExit("not positive at p0")
while hi - lo > a.tol:
    mid = round((lo + hi) / 2, 7)
    Ic, C2 = run(mid, C)
    if Ic > a.thr and min(np.sum(C2**2, 1)) > 0.05:
        lo, C, best = mid, C2, Ic
        np.save(f"runs/{a.tag}_best.npy", C)
    else:
        hi = mid
print(json.dumps({"RESULT": a.tag, "threshold_lo": lo, "first_fail": hi, "Ic_at_lo": best}), flush=True)
