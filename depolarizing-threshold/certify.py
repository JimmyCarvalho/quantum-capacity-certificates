"""Rigorous lower bound on Ic(rho, T_D o G_eta) with ball arithmetic (python-flint / Arb).

* eta is an exact rational; T = 2 eta^2/(1+eta), G = (1+eta)/(2 eta) are exact rationals.
* The input amplitudes are the given float64 numbers, read as exact dyadic rationals (no rounding).
* Every channel coefficient and every output matrix entry is an Arb ball, so the output blocks
  B and RB are known up to rigorous radii.
* Eigenvalues: with an approximate eigenbasis V (float64, exact as a matrix of dyadics), set
  A' = V^T A V = diag(d) + E. Weyl: |lambda_i(A') - d_(i)| <= ||E||_F (d sorted). Ostrowski:
  lambda_i(A') = theta_i lambda_i(A), theta_i in [1-delta, 1+delta], delta >= ||V^T V - I||_F.
  A is a density block, so lambda_i(A) >= 0. This gives an interval per eigenvalue, hence an
  interval for the von Neumann entropy (h(x) = -x ln x is monotone on [0,1/e] and [1/e,1]).
* T_D collapses all output levels > D into one flag state, a channel; by data processing the
  certified number is a lower bound on Ic(rho, G_eta).
"""
import sys
from fractions import Fraction

import flint
import numpy as np
from flint import arb, arb_mat

flint.ctx.prec = 192


def exact(x):
    """float64 -> exact Arb (a float64 is a dyadic rational)."""
    fr = Fraction(float(x))
    return arb(flint.fmpq(fr.numerator, fr.denominator))


def binom_table(nmax):
    """Exact integer binomials as Arb (exact up to 192 bits of mantissa; radii tracked)."""
    from math import comb
    return [[arb(comb(n, k)) for k in range(n + 1)] for n in range(nmax + 1)]


class ArbChannel:
    def __init__(self, eta, Nin, D, T=None, G=None):
        q = lambda fr: arb(flint.fmpq(fr.numerator, fr.denominator))
        if T is None:
            self.eta = Fraction(eta)
            e = self.eta
            T = 2 * e * e / (1 + e)
            G = (1 + e) / (2 * e)
        else:
            T, G = Fraction(T), Fraction(G)
        self.T, self.G = T, G
        self.Nin, self.D = Nin, D
        sT, s1T = q(T).sqrt(), q(1 - T).sqrt()
        r = q((G - 1) / G).sqrt()
        isG = q(1 / G).sqrt()  # G^{-1/2}
        pT = [arb(1)]
        for _ in range(2 * Nin + 2):
            pT.append(pT[-1] * sT)  # sqrt(T)^j
        p1T = [arb(1)]
        for _ in range(2 * Nin + 2):
            p1T.append(p1T[-1] * s1T)
        pr = [arb(1)]
        for _ in range(2 * D + 2):
            pr.append(pr[-1] * r)
        piG = [arb(1)]
        for _ in range(2 * D + 4):
            piG.append(piG[-1] * isG)  # G^{-j/2}
        self.binom = binom_table(D + Nin + 2)
        sq = {}

        def sb(n, k):
            key = (n, k)
            if key not in sq:
                sq[key] = self.binom[n][k].sqrt()
            return sq[key]

        self.loss, self.amp = [], []
        for d in range(Nin + 1):
            n = Nin - d + 1
            L = arb_mat(n, n)
            for m in range(n):
                for k in range(m + 1):
                    # sqrt(C(m,k) C(m+d,k)) T^{m + d/2 - k} (1-T)^k
                    L[m - k, m] = sb(m, k) * sb(m + d, k) * pT[2 * m + d - 2 * k] * p1T[2 * k]
            nout = D - d + 1
            A = arb_mat(nout, n)
            for j in range(n):
                for k in range(nout - j):
                    # sqrt(C(j+k,k) C(j+d+k,k)) G^{-(j + d/2 + 1)} ((G-1)/G)^{k/2 * 2}
                    A[j + k, j] = sb(j + k, k) * sb(j + d + k, k) * piG[2 * j + d + 2] * pr[2 * k]
            self.loss.append(L)
            self.amp.append(A)

    def apply_band(self, d, x):
        """x: list of Arb (band d of the input operator) -> list of Arb (band d of P G(X) P)."""
        v = arb_mat([[xi] for xi in x])
        y = self.amp[d] * (self.loss[d] * v)
        return [y[i, 0] for i in range(y.nrows())]


def entropy_interval(A, dim):
    """Rigorous [lo, hi] of -Tr A ln A for a PSD Arb matrix A (dim x dim). Returns (lo, hi) as Arb."""
    mid = np.array([[float(A[i, j].mid()) for j in range(dim)] for i in range(dim)])
    mid = 0.5 * (mid + mid.T)
    _, V = np.linalg.eigh(mid)
    Vb = arb_mat([[exact(V[i, j]) for j in range(dim)] for i in range(dim)])
    Ap = Vb.transpose() * A * Vb
    R = Vb.transpose() * Vb
    eps2 = arb(0)
    del2 = arb(0)
    d = []
    for i in range(dim):
        for j in range(dim):
            if i == j:
                d.append(Ap[i, i])
                del2 += abs(R[i, i] - 1).upper() ** 2
            else:
                # |x|.upper()^2 is a valid upper bound of x^2 also when the ball contains 0
                # (x ** 2 of such a ball is NaN in python-flint 0.9)
                eps2 += abs(Ap[i, j]).upper() ** 2
                del2 += abs(R[i, j]).upper() ** 2
    eps = eps2.sqrt().upper()
    delta = del2.sqrt().upper()
    if not (eps.is_finite() and delta.is_finite()) or not delta < 0.5:
        raise ValueError("eigenbasis too inaccurate for the Weyl/Ostrowski enclosure")
    one = arb(1)
    lo_tot, hi_tot = arb(0), arb(0)
    inv_e = (-one).exp()
    for di in d:
        lo_l = (di.lower() - eps) / (1 + delta)
        if lo_l < 0:
            lo_l = arb(0)
        hi_l = (di.upper() + eps) / (1 - delta)
        lo_l, hi_l = arb(lo_l.lower()), arb(hi_l.upper())
        # h(x) = -x ln x on [lo_l, hi_l]
        def h(x):
            return arb(0) if x == 0 else -x * x.log()
        hlo, hhi = h(lo_l), h(hi_l)
        if not (hlo.is_finite() and hhi.is_finite()):
            raise ValueError("non-finite entropy enclosure")
        mn = arb(min(hlo.lower(), hhi.lower()))
        if hi_l <= inv_e or lo_l >= inv_e:
            mx = arb(max(hlo.upper(), hhi.upper()))
        else:
            mx = inv_e
        lo_tot += mn
        hi_tot += mx
    return lo_tot.lower(), hi_tot.upper(), float(eps), float(delta)


def h_range(x):
    """[min, max] of h(t) = -t ln t over the part of the ball x inside [0, 1] (x encloses a probability).
    Uses the exact endpoints, so a ball that touches 0 (flag mass below the working precision) is fine."""
    def h(t):
        return arb(0) if t == 0 else -t * t.log()
    a, b = arb(x.lower()), arb(x.upper())
    if a < 0:
        a = arb(0)
    if b > 1:
        b = arb(1)
    inv_e = arb(1) / arb(1).exp()
    ha, hb = h(a), h(b)
    if b < inv_e:
        return ha.lower(), hb.upper()
    if a > inv_e:
        return hb.lower(), ha.upper()
    return min(ha.lower(), hb.lower()), inv_e.upper()


def certify(C, residues, k, eta, Nin, D, verbose=True, T=None, G=None):
    r = len(residues)
    if k is not None and k > 1:
        # the charge-block split (and the diagonal flag block) are exact only for genuine residue codes
        n = np.arange(Nin + 1)
        for i in range(r):
            if np.any(np.asarray(C)[i, : Nin + 1][n % k != residues[i] % k] != 0):
                raise ValueError(f"codeword {i} has weight outside residue {residues[i]} mod {k}; use k = 1")
    ch = ArbChannel(eta, Nin, D, T=T, G=G)
    Cb = [[exact(C[i, n]) for n in range(Nin + 1)] for i in range(r)]
    Z = sum((c * c for row in Cb for c in row), arb(0))
    # Y[i][j]: dict band-wise of P G(phi_i phi_j^T) P / Z
    n1 = D + 1
    Y = {}
    for i in range(r):
        for j in range(r):
            M = [[arb(0)] * n1 for _ in range(n1)]
            for d in range(-Nin, Nin + 1):
                ad = abs(d)
                if d >= 0:
                    x = [Cb[i][m] * Cb[j][m + d] for m in range(Nin + 1 - d)]
                else:
                    x = [Cb[i][m + ad] * Cb[j][m] for m in range(Nin + 1 - ad)]
                if all(xx == 0 for xx in x):
                    continue
                y = ch.apply_band(ad, x)
                for a, ya in enumerate(y):
                    if d >= 0:
                        M[a][a + d] = ya / Z
                    else:
                        M[a + ad][a] = ya / Z
            Y[i, j] = M
    # flags
    f = [arb(1) * sum((Cb[i][m] ** 2 for m in range(Nin + 1)), arb(0)) / Z - sum((Y[i, i][a][a] for a in range(n1)), arb(0))
         for i in range(r)]
    fB = sum(f, arb(0))
    # blocks
    nB = np.arange(n1)
    S_B_lo, S_B_hi = arb(0), arb(0)
    S_RB_lo, S_RB_hi = arb(0), arb(0)
    report = []
    for c in range(k):
        idx = [a for a in nB if a % k == c]
        Bc = arb_mat(len(idx), len(idx))
        for u, a in enumerate(idx):
            for v, b in enumerate(idx):
                Bc[u, v] = sum((Y[i, i][a][b] for i in range(r)), arb(0))
        lo, hi, eps, delta = entropy_interval(Bc, len(idx))
        S_B_lo += lo
        S_B_hi += hi
        report.append(("B", c, len(idx), eps, delta, float(hi - lo)))
        idxRB = [(i, a) for i in range(r) for a in nB if (a - residues[i]) % k == c]
        Rc = arb_mat(len(idxRB), len(idxRB))
        for u, (i, a) in enumerate(idxRB):
            for v, (j, b) in enumerate(idxRB):
                Rc[u, v] = Y[i, j][a][b]
        lo, hi, eps, delta = entropy_interval(Rc, len(idxRB))
        S_RB_lo += lo
        S_RB_hi += hi
        report.append(("RB", c, len(idxRB), eps, delta, float(hi - lo)))
    hlo, hhi = h_range(fB)
    S_B_lo += hlo
    S_B_hi += hhi
    if k is not None and k > 1 and len({x % k for x in residues}) == r:
        # distinct residues: N(phi_i phi_j^T), i != j, has no diagonal, so the R (x) flag block is diagonal
        for fi in f:
            hlo, hhi = h_range(fi)
            S_RB_lo += hlo
            S_RB_hi += hhi
    else:
        # general input: the R (x) flag block is an r x r PSD matrix F with trace fB, and
        # S(F) = h(fB) + fB * H(spectrum / fB) with 0 <= H <= ln r
        S_RB_lo += hlo
        S_RB_hi += hhi + arb(fB.upper()) * arb(r).log()
    ln2 = arb(2).log()
    Ic_lo = (S_B_lo - S_RB_hi) / ln2
    Ic_hi = (S_B_hi - S_RB_lo) / ln2
    if verbose:
        for rep in report:
            print("  block %s c=%d dim=%d eps=%.2e delta=%.2e width(S)=%.2e" % rep)
        print("  flag mass:", [float(x.mid()) for x in f])
    return Ic_lo.lower(), Ic_hi.upper()


if __name__ == "__main__":
    path, res, k, eta, Nin, D = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], int(sys.argv[5]), int(sys.argv[6])
    C = np.load(path)[:, : Nin + 1]
    lo, hi = certify(C, [int(x) for x in res.split(",")], k, eta, Nin, D)
    print("eta =", eta, " Ic in [", lo.str(20, radius=False), ",", hi.str(20, radius=False), "] bits")
    print("CERTIFIED POSITIVE" if lo > 0 else "NOT CERTIFIED")
