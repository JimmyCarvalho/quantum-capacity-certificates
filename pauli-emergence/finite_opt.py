"""Exact optimisation of the coherent information of rank-r inputs on Sym^M (support: k <= K excitations
around the pole of the phase-covariance axis) through N^{(x)M}, N a z-phase-covariant unital Pauli channel.

The channel is linear, so the outputs of all |k><l| (k, l <= K) are computed once (band by band) and every
input state is evaluated by contraction. Everything is finite dimensional: no truncation of the output."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import warnings
warnings.filterwarnings("ignore")
import argparse
import json
import time

import numpy as np
from scipy.optimize import minimize

from bandsym import BandChannel

LN2 = np.log(2.0)


def channel_params(name, p):
    if name == "bb84":
        return 1 - 2 * p, (1 - 2 * p) ** 2
    if name == "2pauli":
        return 1 - 2 * p, 1 - 4 * p
    if name == "depol":  # per-Pauli p
        return 1 - 4 * p, 1 - 4 * p
    if name.startswith("pr"):  # (p_X, p_Y, p_Z) = p (1, 1, r)/(2 + r), total error p, pole on z
        r = float(name[2:])
        return 1 - 2 * p * (1 + r) / (2 + r), 1 - 4 * p / (2 + r)
    raise ValueError(name)


class Transfer:
    """Per-sector linear maps X (on span{|0..K>}) -> Y_j (on the spin-j block)."""

    def __init__(self, M, K, lam_perp, lam_par):
        self.M, self.K = M, K
        ch = BandChannel(M, lam_perp, lam_par)
        self.blocks = []  # list over j of (dim_j, {delta: matrix (n_j^delta, K+1-delta)})
        outs = {}
        for d in range(K + 1):
            cols = np.eye(K + 1 - d)
            outs[d] = np.real(ch.evolve(d, cols.astype(complex)))
            bands = ch.band(d)
            outs[d] = (bands, outs[d])
        b0 = outs[0][0]
        for j in b0.js:
            dimj = int(round(2 * j + 1))
            per = {}
            for d in range(K + 1):
                band, arr = outs[d]
                if j not in band.jpos:
                    continue
                jp = band.jpos[j]
                o, n = band.off[jp], band.lens[jp]
                per[d] = arr[o:o + n, :]
            self.blocks.append((dimj, per))

    def apply(self, bi, X):
        dimj, per = self.blocks[bi]
        Y = np.zeros((dimj, dimj), dtype=X.dtype)
        K = self.K
        for d, A in per.items():
            n = A.shape[0]
            u = np.diagonal(X, offset=d)[: K + 1 - d]
            Y[np.arange(n), np.arange(n) + d] = A @ u
            if d:
                w = np.diagonal(X, offset=-d)[: K + 1 - d]
                Y[np.arange(n) + d, np.arange(n)] = A @ w
        return Y

    def adjoint(self, bi, Z):
        dimj, per = self.blocks[bi]
        K = self.K
        X = np.zeros((K + 1, K + 1), dtype=Z.dtype)
        for d, A in per.items():
            n = A.shape[0]
            zu = Z[np.arange(n), np.arange(n) + d]
            X[np.arange(K + 1 - d), np.arange(K + 1 - d) + d] = A.T @ zu
            if d:
                zl = Z[np.arange(n) + d, np.arange(n)]
                X[np.arange(K + 1 - d) + d, np.arange(K + 1 - d)] = A.T @ zl
        return X


def _ent_log(Mat):
    w, V = np.linalg.eigh(0.5 * (Mat + Mat.conj().T))
    w = np.clip(w, 1e-300, None)
    return -np.sum(w * np.log(w)), (V * np.log(w)) @ V.conj().T


def _blocked_ent_log(Mat, groups):
    if groups is None:
        return _ent_log(Mat)
    L = np.zeros_like(Mat)
    S = 0.0
    for gidx in groups:
        Sg, Lg = _ent_log(Mat[np.ix_(gidx, gidx)])
        S += Sg
        L[np.ix_(gidx, gidx)] = Lg
    return S, L


def residue_groups(dimj, r, residues, kmod):
    """Charge blocks when codeword a lives on k = residues[a] (mod kmod): B_j and RB_j commute with the
    rotation exp(2 pi i (x - a_R)/kmod), x = excitations from the top of the spin-j block."""
    if kmod is None or kmod <= 1:
        return None, None
    x = np.arange(dimj)
    gB = [x[x % kmod == c] for c in range(kmod)]
    full = np.concatenate([a * dimj + x for a in range(r)])
    ch = np.concatenate([(x - residues[a]) % kmod for a in range(r)])
    gRB = [full[ch == c] for c in range(kmod)]
    return [g for g in gB if len(g)], [g for g in gRB if len(g)]


def coherent_info(tr, C, grad=False, residues=None, kmod=None):
    nrm = np.sqrt(np.sum(np.abs(C) ** 2))
    Phi = C / nrm
    r = Phi.shape[0]
    Ic = 0.0
    g = np.zeros_like(Phi)
    for bi, (dimj, _) in enumerate(tr.blocks):
        Y = [[tr.apply(bi, np.outer(Phi[a], Phi[b].conj())) for b in range(r)] for a in range(r)]
        B = sum(Y[a][a] for a in range(r))
        RB = np.block(Y)
        if np.trace(B).real < 1e-300:
            continue
        gB, gRB = residue_groups(dimj, r, residues, kmod) if residues is not None else (None, None)
        SB, LB = _blocked_ent_log(B, gB)
        SRB, LRB = _blocked_ent_log(RB, gRB)
        Ic += (SB - SRB) / LN2
        if grad:
            ZB = tr.adjoint(bi, LB)
            for a in range(r):
                g[a] -= ZB @ Phi[a]
                for b in range(r):
                    g[a] += tr.adjoint(bi, LRB[a * dimj:(a + 1) * dimj, b * dimj:(b + 1) * dimj]) @ Phi[b]
    if not grad:
        return Ic
    g = 2.0 * g / LN2
    gC = (g - np.real(np.vdot(Phi, g)) * Phi) / nrm
    return Ic, gC


def optimise(tr, C0, maxiter=300, real=True, residues=None, kmod=None):
    shape = C0.shape

    def pack(C):
        return C.real.ravel() if real else np.concatenate([C.real.ravel(), C.imag.ravel()])

    def unpack(x):
        if real:
            return x.reshape(shape)
        h = x.size // 2
        return (x[:h] + 1j * x[h:]).reshape(shape)

    scale = [1.0]
    mask = np.ones(shape)
    if residues is not None:
        kk = np.arange(shape[1])
        mask = np.array([(kk % kmod == residues[a]).astype(float) for a in range(shape[0])])
        C0 = C0 * mask

    def f(x):
        Ic, g = coherent_info(tr, unpack(x), grad=True, residues=residues, kmod=kmod)
        if residues is not None:
            g = g * mask
        gg = g.real.ravel() if real else np.concatenate([g.real.ravel(), g.imag.ravel()])
        return -scale[0] * Ic, -scale[0] * gg

    Ic0 = coherent_info(tr, C0, residues=residues, kmod=kmod)
    scale[0] = 1.0 / max(abs(Ic0), 1e-8)
    res = minimize(f, pack(C0), jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "maxcor": 30, "ftol": 1e-15, "gtol": 1e-11})
    C = unpack(res.x)
    C = C / np.linalg.norm(C)
    return coherent_info(tr, C, residues=residues, kmod=kmod), C


def parity_init(K, mean, width, rng=None):
    k = np.arange(K + 1)
    env = np.exp(-((k - mean) ** 2) / (4 * width**2))
    C = np.zeros((2, K + 1))
    C[0, k % 2 == 0] = env[k % 2 == 0]
    C[1, k % 2 == 1] = env[k % 2 == 1]
    if rng is not None:
        C = C + 0.02 * rng.normal(size=C.shape) * env
    return C / np.linalg.norm(C)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--channel", default="2pauli")
    ap.add_argument("--M", type=int, default=30)
    ap.add_argument("--K", type=int, default=30)
    ap.add_argument("--ps", default="0.116,0.117,0.118")
    ap.add_argument("--mean", type=float, default=2.0)
    ap.add_argument("--width", type=float, default=2.0)
    ap.add_argument("--maxiter", type=int, default=300)
    ap.add_argument("--init", default=None)
    ap.add_argument("--complex", action="store_true")
    ap.add_argument("--tag", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(0)
    C = np.load(a.init) if a.init else parity_init(a.K, a.mean, a.width, rng)
    if a.complex:
        C = C.astype(complex)
    for p in [float(x) for x in a.ps.split(",")]:
        t = time.time()
        tr = Transfer(a.M, a.K, *channel_params(a.channel, p))
        t1 = time.time() - t
        Ic, C = optimise(tr, C, maxiter=a.maxiter, real=not a.complex)
        np.save(f"runs/{a.tag}_M{a.M}_p{p}.npy", C)
        k = np.arange(a.K + 1)
        print(json.dumps({"channel": a.channel, "M": a.M, "K": a.K, "p": p, "Ic": Ic,
                          "weights": np.round(np.sum(np.abs(C) ** 2, 1), 4).tolist(),
                          "mean_k": round(float(np.sum(k * np.abs(C) ** 2)), 3),
                          "t_transfer": round(t1, 1), "sec": round(time.time() - t)}), flush=True)
