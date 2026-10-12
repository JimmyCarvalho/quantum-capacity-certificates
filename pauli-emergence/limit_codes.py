"""Residue-class codes for the emergent Gaussian channels of Pauli channels (pole on y), Fock basis.
Reuses the depolarizing-threshold machinery (phase-insensitive channel A_G o L_T, output truncation + flag)."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "depolarizing-threshold"))
import argparse
import json
import time

import numpy as np

from gausschan import EmergentChannel, charge_blocks, coherent_info
from optimize import gaussian_init, optimise, supports


def emergent_TG(channel, p):
    """(T, G) with A_G o L_T equal to the emergent channel for the pole on y (phase-insensitive)."""
    if channel == "bb84":
        lp, lz = 1 - 2 * p, (1 - 2 * p) ** 2
    elif channel == "2pauli":
        lp, lz = 1 - 2 * p, 1 - 4 * p
    elif channel == "depol":  # per-Pauli p, pole arbitrary
        lp, lz = 1 - 4 * p, 1 - 4 * p
    elif channel == "agn":  # additive Gaussian noise of variance p = sigma^2 (vacuum variance 1/2)
        return 1 / (1 + p), 1 + p
    elif channel.startswith("th"):  # thermal attenuator, nu = float(channel[2:]) thermal photons, loss p = 1 - eta
        nu, eta = float(channel[2:]), 1 - p
        G = (1 - eta) * (nu + 0.5) + (1 + eta) / 2
        return eta / G, G
    else:
        raise ValueError(channel)
    tau = lp**2 / lz
    Y = (1 - lp**2) / (2 * lz)
    G = (2 * Y + tau + 1) / 2
    return tau / G, G


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--channel", default="bb84")
    ap.add_argument("--ps", default="0.116,0.118,0.12,0.122")
    ap.add_argument("--k", type=int, default=2)
    ap.add_argument("--res", default="0,1")
    ap.add_argument("--Nin", type=int, default=120)
    ap.add_argument("--D", type=int, default=220)
    ap.add_argument("--mean", type=float, default=3.0)
    ap.add_argument("--width", type=float, default=2.0)
    ap.add_argument("--maxiter", type=int, default=400)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--init", default=None)
    a = ap.parse_args()
    residues = [int(x) for x in a.res.split(",")]
    sup = supports(residues, a.k, a.Nin)
    blocks = charge_blocks(residues, a.k, a.D)
    r = len(residues)
    C = gaussian_init(sup, a.Nin, a.mean, a.width, weights=[1.0 / r] * r)
    if a.init:
        C0 = np.load(a.init)
        C = np.zeros((r, a.Nin + 1))
        m = min(C0.shape[1], a.Nin + 1)
        C[:, :m] = C0[:, :m]
    n = np.arange(a.Nin + 1)
    for p in [float(x) for x in a.ps.split(",")]:
        t = time.time()
        T, G = emergent_TG(a.channel, p)
        ch = EmergentChannel(0.5, a.Nin, a.D, T=T, G=G)
        for _ in range(2):
            Ic0 = coherent_info(ch, C, blocks=blocks)
            Ic, C = optimise(ch, sup, C, maxiter=a.maxiter, verbose=False,
                             scale=1.0 / max(abs(Ic0), 1e-12), blocks=blocks)
        np.save(f"runs/{a.tag}_p{p}.npy", C)
        print(json.dumps({"channel": a.channel, "p": p, "k": a.k, "res": residues, "Ic": Ic,
                          "weights": np.sum(C**2, 1).round(4).tolist(), "mean_n": round(float(np.sum(n * C**2)), 3),
                          "tail": float(np.sum(C[:, -10:] ** 2)), "T": T, "G": G, "sec": round(time.time() - t)}), flush=True)
