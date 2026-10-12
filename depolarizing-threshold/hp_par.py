"""hp_ic.ic_and_grad split over the 2k charge blocks (B and RB for each class c), one process per block.
Each worker rebuilds the channel (cached per eta) and the output bands, diagonalises its block at the
working precision and returns its entropy and its share of the gradient as decimal strings, so nothing
below double precision is lost in transit. Same value and gradient as hp_ic.ic_and_grad."""
import multiprocessing as mp
import os
from fractions import Fraction

import flint
import numpy as np
from flint import arb, arb_mat

import hp_ic
from certify import ArbChannel, exact

_CH = {}
DIGITS = 45


def _channel(eta, Nin, D):
    key = (eta, Nin, D)
    if key not in _CH:
        _CH.clear()
        _CH[key] = ArbChannel(Fraction(eta), Nin, D)
    return _CH[key]


def _block_job(args):
    eta, Nin, D, C, residues, k, kind, c, prec, tol, floor = args
    flint.ctx.prec = prec
    floor = arb(floor)
    ch = _channel(eta, Nin, D)
    code = hp_ic.HPCode(ch, C, residues, k)
    r = len(residues)
    Phi = [[exact(C[i, m]) for m in range(Nin + 1)] for i in range(r)]
    Z = sum((x * x for row in Phi for x in row), arb(0)).sqrt()
    Phi = [[x / Z for x in row] for row in Phi]
    if kind == "B":
        ids = [(None, a) for a in range(D + 1) if a % k == c]
        A = arb_mat([[sum((code.Y(i, i, a, b) for i in range(r)), arb(0)) for (_, b) in ids] for (_, a) in ids])
    else:
        ids = [(i, a) for i in range(r) for a in range(D + 1) if (a - residues[i]) % k == c]
        A = arb_mat([[code.Y(i, j, a, b) for (j, b) in ids] for (i, a) in ids])
    n = A.nrows()
    V = hp_ic.refine_eigh(A, tol=tol)
    dmat = V.transpose() * A * V
    ent = arb(0)
    logs = []
    for u in range(n):
        lam = arb(dmat[u, u].mid())
        if lam > 0:
            ent += hp_ic._h(lam)
        logs.append((lam if lam > floor else arb(floor)).log())
    Lm = V * arb_mat([[logs[u] if u == v else arb(0) for v in range(n)] for u in range(n)]) * V.transpose()
    pos = {key: u for u, key in enumerate(ids)}
    zero = arb(0)

    def get_factory(i, j):
        if kind == "B":
            return lambda a, b: (Lm[pos[None, a], pos[None, b]] if (None, a) in pos and (None, b) in pos else zero)
        return lambda a, b: (Lm[pos[i, a], pos[j, b]] if (i, a) in pos and (j, b) in pos else zero)

    def adj_apply(get, i, j):
        out = [arb(0)] * (Nin + 1)
        for d in range(-Nin, Nin + 1):
            if (d - (residues[j] - residues[i])) % k:
                continue
            ad = abs(d)
            if d >= 0:
                z = [get(a, a + d) for a in range(D + 1 - d)]
            else:
                z = [get(a + ad, a) for a in range(D + 1 - ad)]
            if all(x == 0 for x in z):
                continue
            x = hp_ic._apply_adjoint_band(ch, ad, z)
            for m in range(Nin + 1 - ad):
                if d >= 0:
                    out[m] += x[m] * Phi[j][m + d]
                else:
                    out[m + ad] += x[m] * Phi[j][m]
        return out

    g = [[arb(0)] * (Nin + 1) for _ in range(r)]
    for i in range(r):
        if kind == "B":
            gi = adj_apply(get_factory(None, None), i, i)
            g[i] = [-x for x in gi]
        else:
            for j in range(r):
                gj = adj_apply(get_factory(i, j), i, j)
                g[i] = [x + y for x, y in zip(g[i], gj)]
    return ent.mid().str(DIGITS, radius=False), [[x.mid().str(DIGITS, radius=False) for x in row] for row in g]


class ParIc:
    def __init__(self, procs=None, prec=128, tol=1e-26, floor=1e-60):
        self.pool = mp.get_context("fork").Pool(procs or os.cpu_count())
        self.prec, self.tol, self.floor = prec, tol, floor

    def __call__(self, eta, Nin, D, C, residues=(0, 1), k=3):
        flint.ctx.prec = self.prec
        C = np.asarray(C, dtype=float)[:, : Nin + 1]
        jobs = [(str(eta), Nin, D, C, tuple(residues), k, kind, c, self.prec, self.tol, self.floor)
                for kind in ("RB", "B") for c in range(k)]
        res = self.pool.map(_block_job, jobs, chunksize=1)
        ch = _channel(str(eta), Nin, D)
        code = hp_ic.HPCode(ch, C, residues, k)
        r = len(residues)
        SB, SRB = arb(0), arb(0)
        G = [[arb(0)] * (Nin + 1) for _ in range(r)]
        for (kind, _), (ent, g) in zip([(j[6], j[7]) for j in jobs], res):
            if kind == "B":
                SB += arb(ent)
            else:
                SRB += arb(ent)
            for i in range(r):
                G[i] = [x + arb(y) for x, y in zip(G[i], g[i])]
        fB = sum(code.flags, arb(0))
        SB += hp_ic._h(arb(fB.mid()))
        for f in code.flags:
            SRB += hp_ic._h(arb(f.mid()))
        ln2 = arb(2).log()
        Ic = (SB - SRB) / ln2
        Phi = [[exact(C[i, m]) for m in range(Nin + 1)] for i in range(r)]
        Z = sum((x * x for row in Phi for x in row), arb(0)).sqrt()
        Phi = [[x / Z for x in row] for row in Phi]
        keptm = hp_ic._apply_adjoint_band(ch, 0, [arb(1)] * (D + 1))
        lnfB = arb(fB.mid()).log()
        g = np.zeros_like(C)
        for i in range(r):
            lnfi = arb(code.flags[i].mid()).log()
            gi = [x + (lnfi - lnfB) * (1 - keptm[m]) * Phi[i][m] for m, x in enumerate(G[i])]
            g[i] = [2.0 * float((x / ln2).mid()) for x in gi]
        Phif = np.array([[float(x.mid()) for x in row] for row in Phi])
        gC = (g - np.sum(Phif * g) * Phif) / float(Z.mid())
        return float(Ic.mid()), gC
