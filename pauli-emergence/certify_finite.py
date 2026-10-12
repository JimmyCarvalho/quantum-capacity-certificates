"""Rigorous lower bound (Arb ball arithmetic) on Ic(rho, N^{(x)M}) for a rank-r input on Sym^M, N a
z-phase-covariant unital Pauli channel with rational Bloch factors (lam_perp, lam_par).

N^{(x)M} on permutation-invariant operators equals exp(L), L = sum_i (g D[s-_i] + g D[s+_i] + gD/4 D[sz_i])
with g = -ln(lam_par)/2, gD = -2 ln(lam_perp) + ln(lam_par) (the per-qubit map exp(L_i) is exactly N, for any
sign of gD). L conserves delta = m - m'; each band is exponentiated with arb_mat.exp, which returns rigorous
enclosures. The input lives in the j = M/2 block. Entropies are bounded with the Weyl/Ostrowski argument of
depolarizing-threshold/certify.py. Nothing here is truncated: Sym^M and all output sectors are finite."""
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "depolarizing-threshold"))
import flint
import numpy as np
from flint import arb, arb_mat, fmpq

from bandsym import Band
from certify import entropy_interval, exact

flint.ctx.prec = 160


def q(fr):
    fr = Fraction(fr)
    return arb(fmpq(fr.numerator, fr.denominator))


def sq(fr):
    fr = Fraction(fr)
    return arb(0) if fr == 0 else q(fr).sqrt()


def band_generator_arb(band, g, gD):
    """Exact (ball) generator restricted to one band; same formulas as bandsym.Band.generator."""
    M = band.M
    N = Fraction(M)
    n = band.size
    L = arb_mat(n, n)
    two = Fraction(2)
    for r in range(n):
        j = Fraction(band.j[r]).limit_denominator(2)
        m = Fraction(band.m[r]).limit_denominator(2)
        m1 = Fraction(band.m1[r]).limit_denominator(2)
        if j <= 0:
            deph = gD * q(N / 4)
        else:
            deph = gD / 2 * q(N / 2 - m * m1 * (N / 2 + 1) / j / (j + 1))
        L[r, r] = -(g * q(N) + deph)
        terms = []
        # emission-type sources (m+1, m1+1)
        A, B = m + 1, m1 + 1
        J = j
        if J > 0:
            terms.append((J, A, B, g / 2, (J + A) * (J - A + 1) * (J + B) * (J - B + 1), (N / 2 + 1) / (J * (J + 1))))
        J = j + 1
        terms.append((J, A, B, g / 2, (J + A) * (J + A - 1) * (J + B) * (J + B - 1), (N / 2 + J + 1) / (J * (2 * J + 1))))
        J = j - 1
        if J + 1 > 0 and J >= 0:
            terms.append((J, A, B, g / 2, (J - A + 1) * (J - A + 2) * (J - B + 1) * (J - B + 2), (N / 2 - J) / ((J + 1) * (2 * J + 1))))
        # dephasing sources (m, m1)
        J = j + 1
        terms.append((J, m, m1, gD / 2, (J * J - m * m) * (J * J - m1 * m1), (N / 2 + J + 1) / (J * (2 * J + 1))))
        J = j - 1
        if J >= 0:
            terms.append((J, m, m1, gD / 2, ((J + 1) ** 2 - m * m) * ((J + 1) ** 2 - m1 * m1), (N / 2 - J) / ((J + 1) * (2 * J + 1))))
        # pumping sources (m-1, m1-1)
        A, B = m - 1, m1 - 1
        J = j + 1
        terms.append((J, A, B, g / 2, (J - A - 1) * (J - A) * (J - B - 1) * (J - B), (N / 2 + J + 1) / (J * (2 * J + 1))))
        J = j
        if J > 0:
            terms.append((J, A, B, g / 2, (J + A + 1) * (J - A) * (J + B + 1) * (J - B), (N / 2 + 1) / (J * (J + 1))))
        J = j - 1
        if J >= 0:
            terms.append((J, A, B, g / 2, (J + A + 1) * (J + A + 2) * (J + B + 1) * (J + B + 2), (N / 2 - J) / ((J + 1) * (2 * J + 1))))
        for (J, A, B, rate, rad, fac) in terms:
            c = int(band.index(float(J), float(A), float(B)))
            if c < 0 or rad <= 0 or fac == 0:
                continue
            L[r, c] += rate * sq(rad) * q(fac)
    return L


def certify_finite(C, M, lam_perp, lam_par, residues=None, kmod=None, verbose=True):
    lam_perp, lam_par = Fraction(lam_perp), Fraction(lam_par)
    g = -q(lam_par).log() / 2
    gD = -2 * q(lam_perp).log() + q(lam_par).log()
    r, K1 = C.shape
    K = K1 - 1
    Cb = [[exact(C[a, k]) for k in range(K1)] for a in range(r)]
    Z = sum((x * x for row in Cb for x in row), arb(0))
    outs = {}
    t0 = time.time()
    for d in range(K + 1):
        band = Band(M, d)
        L = band_generator_arb(band, g, gD)
        E = L.exp()
        outs[d] = (band, E)
        if verbose:
            print(f"  band {d}: size {band.size}, exp done ({time.time() - t0:.0f}s)", flush=True)
    js = Band(M, 0).js
    Ic_lo, Ic_hi = arb(0), arb(0)
    for j in js:
        dimj = int(round(2 * j + 1))
        Y = {}
        for a in range(r):
            for b in range(r):
                Mab = [[arb(0)] * dimj for _ in range(dimj)]
                for d in range(K + 1):
                    band, E = outs[d]
                    if j not in band.jpos:
                        continue
                    jp = band.jpos[j]
                    o, n = band.off[jp], band.lens[jp]
                    for k in range(K + 1 - d):
                        up = Cb[a][k] * Cb[b][k + d]
                        lo = Cb[a][k + d] * Cb[b][k]
                        if up == 0 and lo == 0:
                            continue
                        for x in range(n):
                            e = E[o + x, k]
                            if up != 0:
                                Mab[x][x + d] += up * e
                            if d and lo != 0:
                                Mab[x + d][x] += lo * e
                Y[a, b] = Mab
        # blocks with residue charge structure
        xs = np.arange(dimj)
        if residues is None:
            gB, gRB = [xs], [np.arange(r * dimj)]
        else:
            gB = [xs[xs % kmod == c] for c in range(kmod)]
            full = np.concatenate([a * dimj + xs for a in range(r)])
            ch = np.concatenate([(xs - residues[a]) % kmod for a in range(r)])
            gRB = [full[ch == c] for c in range(kmod)]
        for gidx in gB:
            if len(gidx) == 0:
                continue
            A_ = arb_mat(len(gidx), len(gidx))
            for u, x in enumerate(gidx):
                for v, y in enumerate(gidx):
                    A_[u, v] = sum((Y[a, a][x][y] for a in range(r)), arb(0)) / Z
            lo_, hi_, _, _ = entropy_interval(A_, len(gidx))
            Ic_lo += lo_
            Ic_hi += hi_
        for gidx in gRB:
            if len(gidx) == 0:
                continue
            A_ = arb_mat(len(gidx), len(gidx))
            for u, X in enumerate(gidx):
                a, x = divmod(int(X), dimj)
                for v, Yv in enumerate(gidx):
                    b, y = divmod(int(Yv), dimj)
                    A_[u, v] = Y[a, b][x][y] / Z
            lo_, hi_, _, _ = entropy_interval(A_, len(gidx))
            Ic_lo -= hi_
            Ic_hi -= lo_
    ln2 = arb(2).log()
    return (Ic_lo / ln2).lower(), (Ic_hi / ln2).upper()



def outward(lo, hi):
    """Float bounds rounded outward: nextafter of the nearest float is always on the safe side."""
    import math
    return math.nextafter(float(lo), -math.inf), math.nextafter(float(hi), math.inf)



def exact_lambdas(channel, p):
    """Exact (lambda_perp, lambda_par) with the pole on the phase-covariance axis; unknown names fail."""
    if channel == "bb84":
        return 1 - 2 * p, (1 - 2 * p) ** 2
    if channel == "2pauli":
        return 1 - 2 * p, 1 - 4 * p
    if channel == "depol":
        return 1 - 4 * p, 1 - 4 * p
    if channel.startswith("pr"):
        r = Fraction(channel[2:])
        return 1 - 2 * p * (1 + r) / (2 + r), 1 - 4 * p / (2 + r)
    raise ValueError("unknown channel " + channel)


if __name__ == "__main__":
    channel, p, path, M = sys.argv[1], Fraction(sys.argv[2]), sys.argv[3], int(sys.argv[4])
    kmod = int(sys.argv[5]) if len(sys.argv) > 5 else None
    res = [int(x) for x in sys.argv[6].split(",")] if len(sys.argv) > 6 else None
    lp, lz = exact_lambdas(channel, p)
    import hashlib
    print("argv:", " ".join(sys.argv))
    print("input sha256:", hashlib.sha256(open(path, "rb").read()).hexdigest())
    print("python-flint", flint.__version__, "prec", flint.ctx.prec)
    C = np.load(path)
    nz = np.nonzero(np.abs(C).sum(0) > 0)[0]
    C = C[:, : nz.max() + 1]
    t = time.time()
    lo, hi = certify_finite(C, M, lp, lz, residues=res, kmod=kmod)
    olo, ohi = outward(lo, hi)
    print(f"{channel} p={p} M={M} K={C.shape[1]-1}: Ic in [{olo!r}, {ohi!r}] bits (outward-rounded floats; {time.time()-t:.0f}s)")
    print("CERTIFIED POSITIVE" if lo > 0 else "NOT CERTIFIED")
