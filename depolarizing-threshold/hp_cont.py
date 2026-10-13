"""Continuation in eta with high-precision L-BFGS (hp_ic.ic_and_grad) for mod-3 residue codes of the
depolarizing emergent channel. At each eta: optimise from the previous state; if I_c > 0 keep it and step
down, otherwise halve the step. Saves runs/{tag}_eta{eta}.npy for every positive point.
usage: python3 hp_cont.py init.npy eta0 step min_step Nin D iters tag [procs] [grow]
With procs > 1 the six charge blocks are evaluated in parallel (hp_par)."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import json
import sys
import time
from fractions import Fraction

import flint
import numpy as np
from scipy.optimize import minimize

import hp_ic
from certify import ArbChannel
from optimize import pack, supports, unpack

flint.ctx.prec = 128


def optimise(icg, C0, sup, Nin, iters):
    Ic0, _ = icg(C0)
    s = 1.0 / max(abs(Ic0), 1e-40)
    best = [Ic0, C0]

    def f(x):
        C = unpack(x, sup, Nin)
        Ic, g = icg(C)
        if Ic > best[0]:
            best[0], best[1] = Ic, C / np.linalg.norm(C)
        return -s * Ic, -s * pack(g, sup)

    minimize(f, pack(C0, sup), jac=True, method="L-BFGS-B",
             options={"maxiter": iters, "maxcor": 30, "ftol": 1e-15, "gtol": 1e-14})
    return best


if __name__ == "__main__":
    path, eta, step, min_step = sys.argv[1], Fraction(sys.argv[2]), Fraction(sys.argv[3]), Fraction(sys.argv[4])
    Nin, D, iters, tag = int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), sys.argv[8]
    procs = int(sys.argv[9]) if len(sys.argv) > 9 else 1
    grow = float(sys.argv[10]) if len(sys.argv) > 10 else 1.0
    max_step = step
    par = None
    if procs > 1:
        import hp_par
        par = hp_par.ParIc(procs=procs)
    sup = supports([0, 1], 3, Nin)
    S = np.load(path)
    C = np.zeros((2, Nin + 1))
    C[:, : min(Nin + 1, S.shape[1])] = S[:, : Nin + 1]
    C = unpack(pack(C, sup), sup, Nin)
    C /= np.linalg.norm(C)
    good = None
    t0 = time.time()
    while step >= min_step:
        if par is None:
            ch = ArbChannel(eta, Nin, D)
            icg = lambda X: hp_ic.ic_and_grad(ch, X)
        else:
            icg = lambda X, e=eta: par(str(e), Nin, D, X)
        Ic, Cn = optimise(icg, C, sup, Nin, iters)
        n = np.arange(Nin + 1)
        p = np.sum(Cn ** 2, axis=0)
        rec = {"tag": tag, "eta": str(eta), "eta_f": float(eta), "Ic": Ic, "nbar": float(np.sum(n * p)),
               "w": np.sum(Cn ** 2, axis=1).round(4).tolist(), "sec": round(time.time() - t0)}
        print(json.dumps(rec), flush=True)
        if Ic > 0:
            good = (eta, Ic)
            C = Cn
            np.save(f"runs/{tag}_eta{float(eta):.7f}.npy", Cn)
            step = min(step * grow, max_step)  # grow back after a success (grow = 1: fixed/halving only)
            eta -= step
        else:
            eta += step
            step /= 2
            eta -= step
    print(json.dumps({"RESULT": tag, "last_positive_eta": str(good[0]) if good else None,
                      "Ic": good[1] if good else None}), flush=True)
