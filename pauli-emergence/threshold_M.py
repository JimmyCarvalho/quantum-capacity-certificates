"""Finite-M threshold search: continuation in p (step up) then bisection, warm starts, parity codes."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import warnings; warnings.filterwarnings("ignore")
import argparse, json, time
import numpy as np
from finite_opt import Transfer, optimise, channel_params, parity_init, coherent_info

ap = argparse.ArgumentParser()
ap.add_argument("--channel", default="2pauli")
ap.add_argument("--M", type=int, required=True)
ap.add_argument("--K", type=int, required=True)
ap.add_argument("--p0", type=float, default=0.117)
ap.add_argument("--step", type=float, default=0.0005)
ap.add_argument("--tol", type=float, default=0.00005)
ap.add_argument("--mean", type=float, default=4.0)
ap.add_argument("--width", type=float, default=3.0)
ap.add_argument("--init", default=None)
ap.add_argument("--maxiter", type=int, default=250)
ap.add_argument("--kmod", type=int, default=2)
ap.add_argument("--res", default="0,1")
ap.add_argument("--thr", type=float, default=1e-6)
ap.add_argument("--tag", required=True)
ap.add_argument("--complex", action="store_true")
ap.add_argument("--noise", type=float, default=0.0)
a = ap.parse_args()
res = [int(x) for x in a.res.split(",")]
if a.kmod <= 1:
    res = None
if a.init:
    C = np.load(a.init)
    C = C[:, : a.K + 1] if C.shape[1] > a.K else np.pad(C, ((0, 0), (0, a.K + 1 - C.shape[1])))
else:
    C = parity_init(a.K, a.mean, a.width, np.random.default_rng(0))
    if a.kmod != 2:
        k = np.arange(a.K + 1)
        env = np.exp(-((k - a.mean) ** 2) / (4 * a.width**2))
        C = np.array([np.where(k % a.kmod == r_, env, 0.0) for r_ in (res or [0, 0])])
        C = C / np.linalg.norm(C)
if a.noise:
    rng = np.random.default_rng(7)
    C = C + a.noise * np.abs(C).max() * rng.normal(size=C.shape)
if a.complex:
    C = C.astype(complex)
    if a.noise:
        C = C + 1j * a.noise * np.abs(C).max() * np.random.default_rng(8).normal(size=C.shape)
best = None

def run(p, C):
    t = time.time()
    tr = Transfer(a.M, a.K, *channel_params(a.channel, p))
    Ic, C2 = optimise(tr, C, maxiter=a.maxiter, residues=res, kmod=a.kmod, real=not a.complex)
    k = np.arange(a.K + 1)
    print(json.dumps({"channel": a.channel, "M": a.M, "K": a.K, "p": p, "Ic": Ic,
                      "weights": np.round(np.sum(np.abs(C2)**2, 1), 4).tolist(), "mean_k": round(float(np.sum(k * np.abs(C2)**2)), 2),
                      "tail": float(np.sum(np.abs(C2[:, -5:]) ** 2)), "sec": round(time.time() - t)}), flush=True)
    return Ic, C2

p = a.p0
lo, hi = None, None
while hi is None:
    Ic, C2 = run(p, C)
    if Ic > a.thr:
        lo, C, best = p, C2, (p, Ic)
        np.save(f"runs/{a.tag}_best.npy", C)
        p = round(p + a.step, 7)
    else:
        hi = p
        if lo is None:
            raise SystemExit("not positive at p0")
while hi - lo > a.tol:
    mid = round((lo + hi) / 2, 7)
    Ic, C2 = run(mid, C)
    if Ic > a.thr:
        lo, C, best = mid, C2, (mid, Ic)
        np.save(f"runs/{a.tag}_best.npy", C)
    else:
        hi = mid
print(json.dumps({"RESULT": a.tag, "M": a.M, "K": a.K, "threshold_lo": lo, "first_fail": hi, "Ic_at_lo": best[1]}), flush=True)
