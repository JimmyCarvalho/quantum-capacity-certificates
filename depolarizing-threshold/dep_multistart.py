"""Multi-start search at the depolarizing record point in double precision (fast, noise ~3e-14 bits):
mod-3 codes from Gaussian envelopes with random mean, width, weight and a random smooth (skew) deformation,
each optimised by L-BFGS (optimize.optimise). The best final states are re-evaluated with hp_ic.
usage: python3 dep_multistart.py eta n_starts seed tag [Nin D maxiter]"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import json
import sys
import time

import numpy as np

from gausschan import EmergentChannel, charge_blocks, coherent_info
from optimize import optimise, supports

if __name__ == "__main__":
    eta, n_starts, seed, tag = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    Nin = int(sys.argv[5]) if len(sys.argv) > 5 else 200
    D = int(sys.argv[6]) if len(sys.argv) > 6 else 320
    maxiter = int(sys.argv[7]) if len(sys.argv) > 7 else 400
    rng = np.random.default_rng(seed)
    ch = EmergentChannel(eta, Nin, D)
    sup = supports([0, 1], 3, Nin)
    bl = charge_blocks([0, 1], 3, D)
    n = np.arange(Nin + 1)
    for s in range(n_starts):
        mu, sig, q = rng.uniform(40, 75), rng.uniform(2.8, 5.0), rng.uniform(0.40, 0.55)
        amp = rng.uniform(0.0, 0.3)
        C = np.zeros((2, Nin + 1))
        for i, sp in enumerate(sup):
            x = (sp - mu) / sig
            # smooth perturbation of the log-envelope (pointwise noise costs ~1e-6 bits on this plateau)
            env = np.exp(-x ** 2 / 4 + amp * (rng.normal() * x / 3 + rng.normal() * x ** 3 / 27))
            C[i, sp] = np.sqrt(q if i == 0 else 1 - q) * env / np.linalg.norm(env)
        t = time.time()
        Ic0 = coherent_info(ch, C, blocks=bl)
        Ic, C = Ic0, C
        for _ in range(3):
            Ic, C = optimise(ch, sup, C, maxiter=maxiter, verbose=False, scale=1.0 / max(abs(Ic), 1e-14), blocks=bl)
        p = np.sum(C ** 2, axis=0)
        np.save(f"runs/{tag}_s{seed}_{s}.npy", C)
        print(json.dumps({"tag": tag, "start": s, "mu": round(mu, 2), "sig": round(sig, 2), "q": round(q, 3),
                          "pert": round(amp, 3), "Ic0": Ic0, "Ic": Ic, "nbar": round(float(np.sum(n * p)), 2),
                          "sec": round(time.time() - t)}), flush=True)
