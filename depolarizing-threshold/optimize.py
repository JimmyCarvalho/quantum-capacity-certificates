"""Maximise the coherent information of the emergent Gaussian channel over structured codewords.

Codeword i is supported on Fock levels n = residues[i] (mod k), n <= Nin, with real amplitudes.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse
import json
import time

import numpy as np
from scipy.optimize import minimize

from gausschan import EmergentChannel, charge_blocks, coherent_info


def supports(residues, k, Nin):
    return [np.arange(r, Nin + 1, k) for r in residues]


def pack(C, sup):
    return np.concatenate([C[i, s] for i, s in enumerate(sup)])


def unpack(x, sup, Nin):
    C = np.zeros((len(sup), Nin + 1))
    o = 0
    for i, s in enumerate(sup):
        C[i, s] = x[o:o + len(s)]
        o += len(s)
    return C


def gaussian_init(sup, Nin, mean, width, weights=None, rng=None):
    C = np.zeros((len(sup), Nin + 1))
    for i, s in enumerate(sup):
        amp = np.exp(-((s - mean) ** 2) / (4 * width**2))
        if rng is not None:
            amp = amp * (1 + 0.05 * rng.normal(size=amp.shape))
        w = 1.0 if weights is None else weights[i]
        C[i, s] = np.sqrt(w) * amp / np.linalg.norm(amp)
    return C


def optimise(ch, sup, C0, maxiter=500, verbose=True, scale=1e4, blocks=None):
    Nin = ch.Nin
    t0 = time.time()
    it = [0]

    def f(x):
        Ic, g = coherent_info(ch, unpack(x, sup, Nin), grad=True, blocks=blocks)
        it[0] += 1
        if verbose and it[0] % 25 == 0:
            print(f"  it {it[0]:4d}  Ic={Ic:+.12e}  ({time.time() - t0:.0f}s)", flush=True)
        return -scale * Ic, -scale * pack(g, sup)

    res = minimize(f, pack(C0, sup), jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "maxcor": 30, "ftol": 1e-16, "gtol": 1e-12})
    C = unpack(res.x, sup, Nin)
    C /= np.linalg.norm(C)
    return coherent_info(ch, C, blocks=blocks), C


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--eta", type=float, default=0.7294452)
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--res", type=str, default="0,1")
    ap.add_argument("--Nin", type=int, default=220)
    ap.add_argument("--D", type=int, default=320)
    ap.add_argument("--mean", type=float, default=56.5)
    ap.add_argument("--width", type=float, default=8.0)
    ap.add_argument("--maxiter", type=int, default=400)
    ap.add_argument("--q", type=float, default=0.45)
    ap.add_argument("--init", type=str, default=None)
    ap.add_argument("--out", type=str, default=None)
    a = ap.parse_args()
    residues = [int(r) for r in a.res.split(",")]
    sup = supports(residues, a.k, a.Nin)
    ch = EmergentChannel(a.eta, a.Nin, a.D)
    if a.init:
        C0 = np.load(a.init)
        C0 = C0[:, : a.Nin + 1] if C0.shape[1] > a.Nin else np.pad(C0, ((0, 0), (0, a.Nin + 1 - C0.shape[1])))
    else:
        C0 = gaussian_init(sup, a.Nin, a.mean, a.width, weights=[a.q] + [(1 - a.q) / (len(sup) - 1)] * (len(sup) - 1))
    print("initial Ic", coherent_info(ch, C0))
    blocks = charge_blocks(residues, a.k, a.D)
    Ic, C = optimise(ch, sup, C0, maxiter=a.maxiter, blocks=blocks)
    n = np.arange(a.Nin + 1)
    w = np.sum(C**2, axis=1)
    print(json.dumps({"eta": a.eta, "p": 0.75 * (1 - a.eta), "k": a.k, "res": residues, "Ic": Ic,
                      "weights": w.tolist(), "mean_n": float(np.sum(n * C**2))}))
    if a.out:
        np.save(a.out, C)
