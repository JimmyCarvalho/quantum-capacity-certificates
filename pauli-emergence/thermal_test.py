"""Exact Ic of a truncated thermal input on Sym^M (Dicke excitations around the pole) through N^{(x)M},
compared with the emergent Gaussian channel prediction."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse
import json
import time

import numpy as np
from scipy.sparse.linalg import expm_multiply

from dicke_lindblad import PauliSym


def g(x):
    return 0.0 if x <= 0 else (x + 1) * np.log2(x + 1) - x * np.log2(x)


def gaussian_ic(Nbar, X, Y):
    """Ic of a thermal input (mean Nbar) through a single-mode Gaussian channel with diagonal X, Y."""
    a = Nbar + 0.5
    c = np.sqrt(a * a - 0.25)
    V = np.block([[a * np.eye(2), c * np.diag([1, -1]) @ X.T], [X @ np.diag([1, -1]) * c, X @ X.T * a + Y]])
    Om = np.kron(np.eye(2), np.array([[0, 1], [-1, 0]]))
    nus = np.sort(np.abs(np.linalg.eigvals(1j * Om @ V)))[::2]
    VB = V[2:, 2:]
    nB = np.sqrt(np.linalg.det(VB))
    return g(nB - 0.5) - sum(g(n - 0.5) for n in nus)


def entropy_bits(w):
    w = w[w > 1e-300]
    return -np.sum(w * np.log2(w))


def thermal_ic(ch, p):
    """p: probabilities of |k>, k=0..K (diagonal input)."""
    K = len(p) - 1
    sp = ch.space
    pairs = [(k, l) for k in range(K + 1) for l in range(k, K + 1)]
    V = np.zeros((sp.size, len(pairs)), dtype=complex)
    d0 = sp.dims[0]
    for i, (k, l) in enumerate(pairs):
        V[k * d0 + l, i] = np.sqrt(p[k] * p[l])
    W = expm_multiply(ch.L, V)
    out = {pr: sp.blocks(W[:, i]) for i, pr in enumerate(pairs)}
    Ic = 0.0
    for bi, d in enumerate(sp.dims):
        Bdiag = np.real(sum(np.diag(out[(k, k)][bi]) for k in range(K + 1)))
        SB = entropy_bits(Bdiag)
        SRB = 0.0
        for c in range(-K, d):
            idx = [(k, a) for k in range(K + 1) for a in [k + c] if 0 <= a < d]
            n = len(idx)
            if n == 0:
                continue
            A = np.zeros((n, n), dtype=complex)
            for u, (k, a) in enumerate(idx):
                for v, (l, b) in enumerate(idx):
                    if k <= l:
                        A[u, v] = out[(k, l)][bi][a, b]
                    else:
                        A[u, v] = np.conj(out[(l, k)][bi][b, a])
            SRB += entropy_bits(np.linalg.eigvalsh(0.5 * (A + A.conj().T)))
        Ic += SB - SRB
    return Ic


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--channel", default="bb84")
    ap.add_argument("--q", type=float, default=0.115)
    ap.add_argument("--Ms", default="10,20,40")
    ap.add_argument("--Nbar", type=float, default=1.0)
    ap.add_argument("--K", type=int, default=24)
    a = ap.parse_args()
    if a.channel == "bb84":
        s = 1 - 2 * a.q
        lam_perp, lam_par = s, s * s
    elif a.channel == "2pauli":
        lam_perp, lam_par = 1 - 2 * a.q, 1 - 4 * a.q
    else:
        raise SystemExit("unknown channel")
    X = np.eye(2) * lam_perp / np.sqrt(lam_par)
    Y = np.eye(2) * (1 - lam_perp**2) / (2 * lam_par)
    print(json.dumps({"channel": a.channel, "q": a.q, "Nbar": a.Nbar,
                      "gaussian_limit_Ic": gaussian_ic(a.Nbar, X, Y)}), flush=True)
    for M in [int(m) for m in a.Ms.split(",")]:
        K = min(a.K, M)
        x = a.Nbar / (1 + a.Nbar)
        p = (1 - x) * x ** np.arange(K + 1)
        p /= p.sum()
        t = time.time()
        ch = PauliSym(M, lam_perp, lam_par)
        Ic = thermal_ic(ch, p)
        print(json.dumps({"M": M, "K": K, "Ic": Ic, "sec": round(time.time() - t, 1)}), flush=True)
