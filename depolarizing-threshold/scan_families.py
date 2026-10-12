"""Compare residue-class code families: continuation from eta=0.75 down to a target eta, several starts each."""
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
ap.add_argument("--families", type=str, required=True, help="e.g. '3:0,1;4:0,2'")
ap.add_argument("--etas", type=str, default="0.75,0.74,0.735")
ap.add_argument("--means", type=str, default="30,50")
ap.add_argument("--width", type=float, default=2.5)
ap.add_argument("--Nin", type=int, default=200)
ap.add_argument("--D", type=int, default=300)
ap.add_argument("--maxiter", type=int, default=500)
ap.add_argument("--tag", type=str, required=True)
a = ap.parse_args()
etas = [float(e) for e in a.etas.split(",")]
chans = {}
for fam in a.families.split(";"):
    k, res = fam.split(":")
    k = int(k)
    residues = [int(r) for r in res.split(",")]
    sup = supports(residues, k, a.Nin)
    blocks = charge_blocks(residues, k, a.D)
    r = len(residues)
    for mean in [float(m) for m in a.means.split(",")]:
        C = gaussian_init(sup, a.Nin, mean, a.width, weights=[1.0 / r] * r)
        for eta in etas:
            t = time.time()
            if eta not in chans:
                chans[eta] = EmergentChannel(eta, a.Nin, a.D)
            ch = chans[eta]
            for _ in range(2):
                Ic0 = coherent_info(ch, C, blocks=blocks)
                Ic, C = optimise(ch, sup, C, maxiter=a.maxiter, verbose=False,
                                 scale=1.0 / max(abs(Ic0), 1e-12), blocks=blocks)
            name = f"{a.tag}_k{k}_r{'-'.join(map(str, residues))}_m{int(mean)}_eta{eta}"
            np.save(f"runs/{name}.npy", C)
            n = np.arange(a.Nin + 1)
            print(json.dumps({"k": k, "res": residues, "mean0": mean, "eta": eta, "Ic": Ic,
                              "weights": np.sum(C**2, 1).round(5).tolist(),
                              "mean_n": round(float(np.sum(n * C**2)), 2), "sec": round(time.time() - t)}), flush=True)
