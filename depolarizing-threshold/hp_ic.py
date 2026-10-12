"""High-precision coherent information of a residue code through the truncated-with-flag emergent channel.

Same quantity as gausschan.coherent_info (and the lower bound certified by certify.py), but every channel
coefficient and output entry is an Arb ball at `prec` bits and the block spectra come from Arb's
eigensolver at that precision. Double precision has a noise floor of ~3e-14 bits here, which hides the
near-classical plateau where the depolarizing record codes live (I_c ~ 1e-12 at ~56 photons).

The value returned is a float approximation (not a certificate); certify_hp() turns it into an interval.
"""
import sys
from fractions import Fraction

import flint
import numpy as np
from flint import acb_mat, arb, arb_mat

from certify import ArbChannel, exact, h_range


def _h(x):
    return arb(0) if x == 0 else -x * x.log()


class HPCode:
    """Builds the B and RB charge blocks of rho = sum_i |phi_i><phi_i| (rows of C, normalised jointly)."""

    def __init__(self, ch, C, residues, k):
        self.ch, self.residues, self.k = ch, list(residues), k
        Nin, D = ch.Nin, ch.D
        C = np.asarray(C)[:, : Nin + 1]
        n = np.arange(Nin + 1)
        for i, res in enumerate(self.residues):
            if np.any(C[i][n % k != res % k] != 0):
                raise ValueError("codeword %d is not supported on residue %d mod %d" % (i, res, k))
        if len({x % k for x in self.residues}) != len(self.residues):
            raise ValueError("residues must be distinct mod k")
        r = len(self.residues)
        Cb = [[exact(C[i, m]) for m in range(Nin + 1)] for i in range(r)]
        Z = sum((c * c for row in Cb for c in row), arb(0))
        self.bands = {}
        for i in range(r):
            for j in range(i, r):
                for d in range(-Nin, Nin + 1):
                    if (d + self.residues[i] - self.residues[j]) % k:
                        continue  # band d of phi_i phi_j^T vanishes unless d = res_j - res_i (mod k)
                    ad = abs(d)
                    if d >= 0:
                        x = [Cb[i][m] * Cb[j][m + d] for m in range(Nin + 1 - d)]
                    else:
                        x = [Cb[i][m + ad] * Cb[j][m] for m in range(Nin + 1 - ad)]
                    y = ch.apply_band(ad, x)
                    self.bands[i, j, d] = [v / Z for v in y]
        self.flags = [sum((Cb[i][m] ** 2 for m in range(Nin + 1)), arb(0)) / Z
                      - sum(self.bands[i, i, 0], arb(0)) for i in range(r)]

    def Y(self, i, j, a, b):
        """<a| N(phi_i phi_j^T) |b> (0 outside the stored bands)."""
        if i > j:
            return self.Y(j, i, b, a)
        d = b - a
        y = self.bands.get((i, j, d))
        if y is None:
            return arb(0)
        return y[min(a, b)]

    def blocks(self):
        D, k, r = self.ch.D, self.k, len(self.residues)
        out = []
        for c in range(k):
            idx = [a for a in range(D + 1) if a % k == c]
            B = arb_mat([[sum((self.Y(i, i, a, b) for i in range(r)), arb(0)) for b in idx] for a in idx])
            idxRB = [(i, a) for i in range(r) for a in range(D + 1) if (a - self.residues[i]) % k == c]
            RB = arb_mat([[self.Y(i, j, a, b) for (j, b) in idxRB] for (i, a) in idxRB])
            out.append((B, RB))
        return out


def spectrum_entropy(A):
    """-Tr A ln A from the approximate eigenvalues of A at the working precision (negative ones -> 0)."""
    ev = acb_mat(A).eig(algorithm="approx")
    s = arb(0)
    for e in ev:
        lam = arb(e.real.mid())
        if lam > 0:
            s += _h(lam)
    return s


def ic(ch, C, residues=(0, 1), k=3):
    """I_c in bits (Arb value, accurate far below double precision; not a rigorous enclosure)."""
    code = HPCode(ch, C, residues, k)
    SB, SRB = arb(0), arb(0)
    for B, RB in code.blocks():
        SB += spectrum_entropy(B)
        SRB += spectrum_entropy(RB)
    fB = sum(code.flags, arb(0))
    SB += _h(arb(fB.mid())) if fB.mid() > 0 else arb(0)
    for f in code.flags:
        SRB += _h(arb(f.mid())) if f.mid() > 0 else arb(0)
    return (SB - SRB) / arb(2).log()


def _to_arb_mat(X):
    """float64 array -> exact Arb matrix (every float64 is a dyadic rational)."""
    return arb_mat(np.asarray(X, dtype=float).tolist())


def _mid(A):
    n, m = A.nrows(), A.ncols()
    return np.array([[float(A[i, j].mid()) for j in range(m)] for i in range(n)])


def _rotation(M, d):
    """Orthogonal Q (double precision) that nearly diagonalises the nearly diagonal symmetric M.
    Strongly coupled indices (|M_ij| > 0.01 |d_i - d_j|) form clusters, each diagonalised exactly after
    scaling to its own size; the remaining weak couplings get first-order Jacobi angles. A simultaneous
    first-order update alone does not converge when near-degenerate pairs share an index."""
    from scipy.sparse.csgraph import connected_components
    n = len(d)
    off = M - np.diag(d)
    gap = np.abs(d[:, None] - d[None, :])
    strong = np.abs(off) > 0.01 * gap
    np.fill_diagonal(strong, False)
    ncomp, lab = connected_components(strong.astype(np.int8), directed=False)
    Q = np.eye(n)
    for c in range(ncomp):
        idx = np.nonzero(lab == c)[0]
        if len(idx) < 2:
            continue
        sub = M[np.ix_(idx, idx)]
        sc = np.abs(sub).max()
        _, W = np.linalg.eigh(sub / sc)
        Q[np.ix_(idx, idx)] = W
    with np.errstate(divide="ignore", invalid="ignore"):
        X = off / (d[None, :] - d[:, None])  # X[j, i] = M_ji / (d_i - d_j)
    X[strong | ~np.isfinite(X)] = 0.0
    same = lab[:, None] == lab[None, :]
    X[same] = 0.0
    X = 0.5 * (X - X.T)
    return (np.eye(n) + X + 0.5 * (X @ X)) @ Q


def refine_eigh(A, iters=8, tol=1e-34):
    """Real, nearly orthonormal eigenbasis V of a symmetric Arb matrix A at the working precision.
    Start from double-precision eigh; each round re-orthonormalises V (Newton-Schulz polar step in Arb)
    and rotates it by the 2x2 Jacobi angles of the off-diagonal part of V^T A V (angles in double
    precision, applied exactly), until max|offdiag| <= tol * max|diag|. Off-diagonal entries go
    1e-17 -> 1e-34 -> working precision. The basis is only approximate: callers bound the error
    (Weyl/Ostrowski), nothing relies on it."""
    n = A.nrows()
    mid = _mid(A)
    _, V0 = np.linalg.eigh(0.5 * (mid + mid.T))
    V = _to_arb_mat(V0)
    I3 = arb_mat(n, n)
    for i in range(n):
        I3[i, i] = arb(3)

    def exact_mid(W):  # drop the radii: V only has to be some exactly known matrix
        return arb_mat([[W[i, j].mid() for j in range(n)] for i in range(n)])

    for _ in range(iters):
        for _ in range(2):
            V = exact_mid(V * (I3 - V.transpose() * V) * arb(0.5))
        M = _mid(V.transpose() * A * V)
        d = np.diag(M)
        off = np.abs(M - np.diag(d)).max()
        if off <= tol * np.abs(d).max():
            return V
        Q = _rotation(M, d)
        V = exact_mid(V * _to_arb_mat(Q))
    for _ in range(2):
        V = exact_mid(V * (I3 - V.transpose() * V) * arb(0.5))
    return V


def entropy_enclosure(A):
    """Rigorous [lo, hi] of -Tr A ln A for a real symmetric PSD Arb matrix A with trace <= 1.
    With V from refine_eigh (an exact matrix of dyadics), A' = V^T A V = diag(d) + E.
    Weyl: the sorted eigenvalues of A' are within ||E||_F of the sorted d. Ostrowski: lambda_i(A') =
    theta_i lambda_i(A) with theta_i in [1 - delta, 1 + delta], delta >= ||V^T V - I||_F. A is PSD."""
    n = A.nrows()
    V = refine_eigh(A)
    Ap = V.transpose() * A * V
    G = V.transpose() * V
    eps2, del2 = arb(0), arb(0)
    d = []
    for i in range(n):
        for j in range(n):
            if i == j:
                d.append(Ap[i, i])
                del2 += abs(G[i, i] - 1).upper() ** 2
            else:
                eps2 += abs(Ap[i, j]).upper() ** 2
                del2 += abs(G[i, j]).upper() ** 2
    eps, delta = arb(eps2.sqrt().upper()), arb(del2.sqrt().upper())
    if not (eps.is_finite() and delta.is_finite()) or not delta < 0.5:
        raise ValueError("eigenbasis too inaccurate for the Weyl/Ostrowski enclosure")
    lo, hi = arb(0), arb(0)
    for di in d:
        a = (di.lower() - eps) / (1 + delta)
        b = (di.upper() + eps) / (1 - delta)
        a = arb(0) if a < 0 else arb(a.lower())
        b = arb(b.upper())
        hl, hh = h_range(a.union(b))
        lo += hl
        hi += hh
    return lo.lower(), hi.upper(), float(eps)


def certify_hp(ch, C, residues=(0, 1), k=3):
    """Rigorous lower/upper bound on I_c (bits) for the truncated-with-flag channel (a lower bound on I_c of
    the untruncated channel by data processing). Returns (lo, hi, largest Weyl radius)."""
    code = HPCode(ch, C, residues, k)
    lo_tot, hi_tot, rad = arb(0), arb(0), 0.0
    for B, RB in code.blocks():
        lo, hi, r1 = entropy_enclosure(B)
        lo_tot += lo
        hi_tot += hi
        lo, hi, r2 = entropy_enclosure(RB)
        lo_tot -= hi
        hi_tot -= lo
        rad = max(rad, r1, r2)
    fB = sum(code.flags, arb(0))
    hlo, hhi = h_range(fB)
    lo_tot += hlo
    hi_tot += hhi
    for f in code.flags:  # distinct residues: the R (x) flag block is diagonal
        hlo, hhi = h_range(f)
        lo_tot -= hhi
        hi_tot -= hlo
    ln2 = arb(2).log()
    return (lo_tot / ln2).lower(), (hi_tot / ln2).upper(), rad


def _apply_adjoint_band(ch, d, z):
    """Band d of N^dagger(Z) from band d of Z (list of Arb of length D + 1 - d): W_d^T z."""
    if not hasattr(ch, "_lossT"):
        ch._lossT = [L.transpose() for L in ch.loss]
        ch._ampT = [A.transpose() for A in ch.amp]
    v = arb_mat([[x] for x in z])
    y = ch._lossT[d] * (ch._ampT[d] * v)
    return [y[i, 0] for i in range(y.nrows())]


def ic_and_grad(ch, C, residues=(0, 1), k=3, iters=8, floor=1e-60, tol=1e-26):
    """I_c (float, bits) and its gradient (float64 array shaped like C, before normalisation), both from
    working-precision spectra and matrix logarithms (refine_eigh). Same formula as gausschan.coherent_info:
      g_i = N^+(L_RB,ii - L_B) phi_i + sum_{j != i} N^+(L_RB,ij) phi_j + (ln f_i - ln f_B)(1 - N^+(1)) phi_i."""
    Nin, D = ch.Nin, ch.D
    C = np.asarray(C, dtype=float)[:, : Nin + 1]
    code = HPCode(ch, C, residues, k)
    r = len(code.residues)
    Phi = [[exact(C[i, m]) for m in range(Nin + 1)] for i in range(r)]
    Z = sum((x * x for row in Phi for x in row), arb(0)).sqrt()
    Phi = [[x / Z for x in row] for row in Phi]
    SB, SRB = arb(0), arb(0)
    # L[(i, j)][(a, b)] = log-matrix entry coupling (i, a) and (j, b); LB[(a, b)] Bob's log
    LRB, LB = {}, {}
    for c in range(k):
        idx = [a for a in range(D + 1) if a % k == c]
        idxRB = [(i, a) for i in range(r) for a in range(D + 1) if (a - code.residues[i]) % k == c]
        for kind, ids in (("B", idx), ("RB", idxRB)):
            if kind == "B":
                A = arb_mat([[sum((code.Y(i, i, a, b) for i in range(r)), arb(0)) for b in ids] for a in ids])
            else:
                A = arb_mat([[code.Y(i, j, a, b) for (j, b) in ids] for (i, a) in ids])
            n = A.nrows()
            V = refine_eigh(A, iters=iters, tol=tol)
            dmat = V.transpose() * A * V
            ent = arb(0)
            logs = []
            for u in range(n):
                lam = arb(dmat[u, u].mid())
                if lam > 0:
                    ent += _h(lam)
                logs.append((lam if lam > floor else arb(floor)).log())
            Lm = V * arb_mat([[logs[u] if u == v else arb(0) for v in range(n)] for u in range(n)]) * V.transpose()
            if kind == "B":
                SB += ent
                for u, a in enumerate(ids):
                    for v, b in enumerate(ids):
                        LB[a, b] = Lm[u, v]
            else:
                SRB += ent
                for u, (i, a) in enumerate(ids):
                    for v, (j, b) in enumerate(ids):
                        LRB[i, j, a, b] = Lm[u, v]
    fB = sum(code.flags, arb(0))
    SB += _h(arb(fB.mid()))
    for f in code.flags:
        SRB += _h(arb(f.mid()))
    ln2 = arb(2).log()
    Ic = (SB - SRB) / ln2

    def adj_apply(get, i, j):
        """(N^+(Z) phi_j) where Z[a, b] = get(a, b) is supported on b - a = res_j - res_i (mod k)."""
        out = [arb(0)] * (Nin + 1)
        for d in range(-D, D + 1):
            if (d - (code.residues[j] - code.residues[i])) % k:
                continue
            ad = abs(d)
            if ad > Nin:
                continue
            if d >= 0:
                z = [get(a, a + d) for a in range(D + 1 - d)]
            else:
                z = [get(a + ad, a) for a in range(D + 1 - ad)]
            x = _apply_adjoint_band(ch, ad, z)  # band d of N^+(Z), input levels
            for m in range(Nin + 1 - ad):
                if d >= 0:
                    out[m] += x[m] * Phi[j][m + d]
                else:
                    out[m + ad] += x[m] * Phi[j][m]
        return out

    zero = arb(0)
    # (1 - N^+(1)) is diagonal: 1 - (mass of N(|m><m|) kept below the cut)
    one_vec = [arb(1)] * (D + 1)
    keptm = _apply_adjoint_band(ch, 0, one_vec)
    g = np.zeros_like(C)
    lnfB = arb(fB.mid()).log()
    for i in range(r):
        gi = adj_apply(lambda a, b: LRB.get((i, i, a, b), zero) - LB.get((a, b), zero), i, i)
        for j in range(r):
            if j != i:
                gj = adj_apply(lambda a, b: LRB.get((i, j, a, b), zero), i, j)
                gi = [x + y for x, y in zip(gi, gj)]
        lnfi = arb(code.flags[i].mid()).log()
        gi = [x + (lnfi - lnfB) * (1 - keptm[m]) * Phi[i][m] for m, x in enumerate(gi)]
        g[i] = [2.0 * float((x / ln2).mid()) for x in gi]
    Phif = np.array([[float(x.mid()) for x in row] for row in Phi])
    nrm = float(Z.mid())
    gC = (g - np.sum(Phif * g) * Phif) / nrm
    return float(Ic.mid()), gC


if __name__ == "__main__":
    # usage: python3 hp_ic.py state.npy eta Nin D [prec] [--certify]
    path, eta, Nin, D = sys.argv[1], Fraction(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    flint.ctx.prec = int(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[5].isdigit() else 128
    ch = ArbChannel(eta, Nin, D)
    C = np.load(path)
    print("Ic (high precision) =", ic(ch, C).str(12), "bits")
    if "--certify" in sys.argv:
        lo, hi, rad = certify_hp(ch, C)
        print("Ic in [", lo.str(15, radius=False), ",", hi.str(15, radius=False), "] bits; largest Weyl radius %.1e" % rad)
        print("CERTIFIED POSITIVE" if lo > 0 else "NOT CERTIFIED")
