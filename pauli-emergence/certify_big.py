"""Scalable rigorous certificate for Ic(rho, N^{(x)M}) (rank-2 residue codes on Sym^M), for large M.

Instead of exponentiating each band in ball arithmetic (cost ~ n^3), the band ODE v' = L v is integrated in
float64 with an explicit, rigorous error budget:

1. Generator. Every entry L_rc is evaluated in Arb from the exact formulas (rates -ln(lam)/2 included). The
   float matrix used is Lf = round(mid(L)); |Lf - L| <= rad + ulp. Hence ||L - Lf||_1 (max column sum) is a
   rigorous bound, as are ||Lf||_1 and the logarithmic norm mu1(Lf) = max_c (Lf_cc + sum_{r!=c} |Lf_rc|),
   which gives ||exp(t Lf)||_1 <= exp(t * max(mu1, 0)).
2. Integration. exp(Lf) = exp(h Lf)^s with h = 2^-k, h ||Lf||_1 <= 1/4. Each step is a degree-m Taylor sum in
   float64. Per step and per column we bound in 1-norm: truncation, the rounding of every sparse product
   (gamma_q |Lf||t|, q = nonzeros per row + 1) propagated through the remaining Taylor recursion, and the
   rounding of the summation. Step errors are propagated to t = 1 with exp((1 - t) mu1+).
3. Model error. ||exp(L) - exp(Lf)||_1 <= exp(max(mu1(L), mu1(Lf), 0)) ||L - Lf||_1 (Duhamel).
4. Assembly. The float blocks B_j, RB_j of the (unnormalised) code are exact dyadic matrices sigma; their
   entropies are enclosed with depolarizing-threshold/certify.entropy_interval. The true blocks differ from
   sigma by at most eps in total absolute entry sum (all budget terms above, weighted by |c_a(k) c_b(l)|, plus
   the rounding of the assembly), hence in trace norm. With x = clip(eig(sigma), 0), ||eig(true) - x||_1 <=
   eps' = eps + (negative eigenvalue mass of sigma), and |sum h(eig) - sum h(x)| <= eps' ln(D/eps') for
   eps' <= 1/e (h(t) = -t ln t, D = number of eigenvalues), giving the final interval.
Every float bound is computed from a handful of float operations; a factor 2 safety margin is applied to all
budget totals, far above their own rounding error (< 1e-12 relative)."""
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "depolarizing-threshold"))
import flint
import numpy as np
from flint import arb, arb_mat
from scipy.sparse import coo_matrix

from bandsym import Band
from certify import exact
from certify_finite import q, sq

flint.ctx.prec = 128
U = 2.0 ** -53
LD = os.environ.get("CERT_LD") == "1"   # integrate in x87 extended precision (u = 2^-64)
if LD:
    # the budget below assumes a 64-bit mantissa and correctly rounded string parsing; long double is plain
    # double on some platforms (MSVC, macOS arm64) or can run with 53-bit precision control
    _one, _ulp = np.longdouble(1), np.longdouble(2.0 ** -63)
    if (np.finfo(np.longdouble).nmant != 63 or (_one + _ulp) - _one != _ulp
            or np.longdouble("1.0000000000000000001") == _one):
        raise RuntimeError("CERT_LD=1 needs x87 extended precision (64-bit mantissa) for np.longdouble")
UP = float(np.finfo(np.longdouble).eps) / 2 if LD else U
DT = np.longdouble if LD else np.float64


def band_entries_arb(band, g, gD):
    """Same formulas as certify_finite.band_generator_arb, returned as sparse (rows, cols, arb values)."""
    N = Fraction(band.M)
    rows, cols, vals = [], [], []
    for r in range(band.size):
        j = Fraction(band.j[r]).limit_denominator(2)
        m = Fraction(band.m[r]).limit_denominator(2)
        m1 = Fraction(band.m1[r]).limit_denominator(2)
        deph = gD * q(N / 4) if j <= 0 else gD / 2 * q(N / 2 - m * m1 * (N / 2 + 1) / j / (j + 1))
        rows.append(r); cols.append(r); vals.append(-(g * q(N) + deph))
        terms = []
        A, B = m + 1, m1 + 1
        if j > 0:
            terms.append((j, A, B, g / 2, (j + A) * (j - A + 1) * (j + B) * (j - B + 1), (N / 2 + 1) / (j * (j + 1))))
        J = j + 1
        terms.append((J, A, B, g / 2, (J + A) * (J + A - 1) * (J + B) * (J + B - 1), (N / 2 + J + 1) / (J * (2 * J + 1))))
        J = j - 1
        if J >= 0:
            terms.append((J, A, B, g / 2, (J - A + 1) * (J - A + 2) * (J - B + 1) * (J - B + 2), (N / 2 - J) / ((J + 1) * (2 * J + 1))))
        J = j + 1
        terms.append((J, m, m1, gD / 2, (J * J - m * m) * (J * J - m1 * m1), (N / 2 + J + 1) / (J * (2 * J + 1))))
        J = j - 1
        if J >= 0:
            terms.append((J, m, m1, gD / 2, ((J + 1) ** 2 - m * m) * ((J + 1) ** 2 - m1 * m1), (N / 2 - J) / ((J + 1) * (2 * J + 1))))
        A, B = m - 1, m1 - 1
        J = j + 1
        terms.append((J, A, B, g / 2, (J - A - 1) * (J - A) * (J - B - 1) * (J - B), (N / 2 + J + 1) / (J * (2 * J + 1))))
        if j > 0:
            terms.append((j, A, B, g / 2, (j + A + 1) * (j - A) * (j + B + 1) * (j - B), (N / 2 + 1) / (j * (j + 1))))
        J = j - 1
        if J >= 0:
            terms.append((J, A, B, g / 2, (J + A + 1) * (J + A + 2) * (J + B + 1) * (J + B + 2), (N / 2 - J) / ((J + 1) * (2 * J + 1))))
        for (JJ, AA, BB, rate, rad, fac) in terms:
            c = int(band.index(float(JJ), float(AA), float(BB)))
            if c < 0 or rad <= 0 or fac == 0:
                continue
            rows.append(r); cols.append(c); vals.append(rate * sq(rad) * q(fac))
    return np.array(rows), np.array(cols), vals


def up(x):
    """Upper bound of a nonnegative arb as float."""
    return float(arb(x).upper()) * (1 + 4 * U)


def band_operator(band, g, gD):
    rows, cols, vals = band_entries_arb(band, g, gD)
    n = band.size
    if LD:
        Lf = np.array([np.longdouble(v.mid().str(30, radius=False)) for v in vals], dtype=np.longdouble)
        # decimal string with 30 digits then rounding to 64-bit mantissa: relative error <= UP + 1e-29
        dl = np.array([float(v.rad()) for v in vals]) * (1 + 4 * U) + np.abs(Lf).astype(float) * (UP + 1e-29) * (1 + 4 * U)
    else:
        Lf = np.array([float(v.mid()) for v in vals])
        # |Lf - L| <= rad + |Lf| u (rounding of mid to float)
        dl = np.array([float(v.rad()) for v in vals]) * (1 + 4 * U) + np.abs(Lf) * U
    # column sums
    dL1 = np.zeros(n)
    np.add.at(dL1, cols, dl)
    absL = np.zeros(n)
    np.add.at(absL, cols, np.abs(Lf).astype(float))
    diag = np.zeros(n)
    dmask = rows == cols
    np.add.at(diag, cols[dmask], Lf[dmask].astype(float))
    mu_f = np.max(absL - np.abs(diag) + diag)          # mu1(Lf)
    normL = np.max(absL)
    dLmax = np.max(dL1)
    mu_true = mu_f + dLmax                               # mu1(L) <= mu1(Lf) + ||L - Lf||_1
    nnz_row = np.bincount(rows, minlength=n).max()
    Lcsr = coo_matrix((Lf, (rows, cols)), shape=(n, n)).tocsr()
    Labs = coo_matrix((np.abs(Lf), (rows, cols)), shape=(n, n)).tocsr()
    return Lcsr, Labs, dict(mu=mu_f, mu_true=mu_true, normL=normL, dL=dLmax, q=int(nnz_row) + 1)


def propagate(L, info, V, m=16):
    """Float integration of exp(Lf) V with a rigorous per-column 1-norm error bound (see module docstring).
    Per step, exact terms T_n = (h Lf)^n v / n!, computed t_n; e_n = t_n - T_n obeys
    ||e_n|| <= (h|Lf|/n) ||e_{n-1}|| + ||delta_n||, ||delta_n|| <= (h/n)(gamma_q + 4u) |Lf| ||t_{n-1}||."""
    normL, qn = info["normL"], info["q"]
    mu = max(info["mu"] + 1e-12 * normL, 0.0)
    k = 0
    while normL * 2.0 ** (-k) > 0.25:
        k += 1
    h = 2.0 ** (-k)
    s = 2 ** k
    gq = qn * UP / (1 - qn * UP)
    gm = (m + 1) * UP / (1 - (m + 1) * UP)
    hn = h * normL
    t = 1.0
    for n in range(1, m + 1):
        t *= hn / n
    trunc = 0.0
    for n in range(m + 1, m + 80):
        t *= hn / n
        trunc += t
    trunc *= 2.0  # tail of a convergent series dominated by its first terms (hn <= 1/4)
    grow = np.exp(h * mu)
    v = V.astype(DT)
    err = np.zeros(V.shape[1])
    for _ in range(s):
        vnorm = np.abs(v).sum(axis=0).astype(float)
        tn = v
        acc = v.copy()
        e = np.zeros(V.shape[1])
        sumE = np.zeros(V.shape[1])
        sumT = vnorm.copy()
        for n in range(1, m + 1):
            tprev = np.abs(tn).sum(axis=0).astype(float)
            # h/n formed in the working precision: a float64 quotient would carry a 2^-53 error that the
            # 4u term does not cover when u = 2^-64
            tn = (L @ tn) * (DT(h) / DT(n))
            e = (hn / n) * e + (h / n) * (gq + 4 * UP) * normL * tprev
            sumE += e
            sumT += np.abs(tn).sum(axis=0).astype(float)
            acc += tn
        err = err * grow + sumE + gm * sumT + trunc * vnorm
        v = acc
    # final rounding to float64 for the assembly: |fl(v) - v| <= U |v|
    vf = v.astype(np.float64)
    err = err + U * np.abs(vf).sum(axis=0) * (1 + 4 * U)
    return vf, err


def _band_job(job):
    M, d, K, lam_perp, lam_par = job
    g = -q(lam_par).log() / 2
    gD = -2 * q(lam_perp).log() + q(lam_par).log()
    band = Band(M, d)
    Lf, Labs, info = band_operator(band, g, gD)
    if LD:
        Lf = Lf.astype(np.longdouble)
    V = np.zeros((band.size, K + 1 - d))
    V[np.arange(K + 1 - d), np.arange(K + 1 - d)] = 1.0
    W, e_int = propagate(Lf, info, V)
    e_model = np.exp(max(info["mu_true"], info["mu"], 0.0)) * info["dL"]
    return d, W, 2.0 * (e_int + e_model), info, band.size  # safety factor 2


def certify_big(C, M, lam_perp, lam_par, residues, kmod, verbose=True):
    from certify import entropy_interval
    lam_perp, lam_par = Fraction(lam_perp), Fraction(lam_par)
    g = -q(lam_par).log() / 2
    gD = -2 * q(lam_perp).log() + q(lam_par).log()
    r, K1 = C.shape
    K = K1 - 1
    # the entropies are computed sector by sector of the residue charge; that equals the true value only if
    # the code really lives on its residues
    if not (np.all(np.isfinite(C)) and 0 <= K <= M and len(residues) == r
            and all(0 <= x < kmod for x in residues)):
        raise ValueError("invalid code: need finite C, K <= M, one residue in [0, kmod) per codeword")
    ks_all = np.arange(K + 1)
    for a in range(r):
        if np.any(C[a, ks_all % kmod != residues[a]] != 0):
            raise ValueError(f"codeword {a} has weight outside residue {residues[a]} mod {kmod}")
    outs, errs = {}, {}
    t0 = time.time()
    procs = int(os.environ.get("CERT_PROCS", "1"))
    jobs = [(M, d, K, lam_perp, lam_par) for d in range(K + 1)]
    if procs > 1:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(procs) as pool:
            results = pool.imap_unordered(_band_job, jobs)
            for d, W, err, info, n in results:
                outs[d] = (Band(M, d), W)
                errs[d] = err
                if verbose:
                    print(f"  band {d}: n={n} |L|={info['normL']:.1f} mu={info['mu']:.2e} dL={info['dL']:.1e} "
                          f"err<= {err.max():.2e} ({time.time() - t0:.0f}s)", flush=True)
    else:
        for job in jobs:
            d, W, err, info, n = _band_job(job)
            outs[d] = (Band(M, d), W)
            errs[d] = err
            if verbose:
                print(f"  band {d}: n={n} |L|={info['normL']:.1f} mu={info['mu']:.2e} dL={info['dL']:.1e} "
                      f"err<= {err.max():.2e} ({time.time() - t0:.0f}s)", flush=True)
    js = Band(M, 0).js
    Cn = np.asarray(C, dtype=float)
    Z = sum((exact(x) ** 2 for x in Cn.ravel()), arb(0))
    # error budget (absolute entry sums, unnormalised) for B and for RB
    epsB, epsRB = 0.0, 0.0
    for d in range(K + 1):
        ks = np.arange(K + 1 - d)
        for a in range(r):
            for b in range(r):
                w = np.abs(Cn[a, ks] * Cn[b, ks + d]) * errs[d]
                wl = np.abs(Cn[a, ks + d] * Cn[b, ks]) * errs[d] if d else 0.0
                tot = w.sum() + (np.sum(wl) if d else 0.0)
                epsRB += tot
                if a == b:
                    epsB += tot
    # rounding of the float assembly: every entry is a sum of <= K+1 products
    asm = 0.0
    for d in range(K + 1):
        band, W = outs[d]
        colnorm = np.abs(W).sum(axis=0)
        ks = np.arange(K + 1 - d)
        for a in range(r):
            for b in range(r):
                asm += ((K + 3) * U) * (np.abs(Cn[a, ks] * Cn[b, ks + d]) * colnorm).sum() * 2
    epsB = 2.0 * (epsB + asm)
    epsRB = 2.0 * (epsRB + asm)
    S_B = [arb(0), arb(0)]
    S_RB = [arb(0), arb(0)]
    negB, negRB, DB, DRB = 0.0, 0.0, 0, 0
    for jpos, j in enumerate(js):
        dimj = int(round(2 * j + 1))
        Y = {}
        for a in range(r):
            for b in range(r):
                Mab = np.zeros((dimj, dimj))
                for d in range(K + 1):
                    band, W = outs[d]
                    if j not in band.jpos:
                        continue
                    jp = band.jpos[j]
                    o, n = band.off[jp], band.lens[jp]
                    blk = W[o:o + n, :]
                    ks = np.arange(K + 1 - d)
                    u = Cn[a, ks] * Cn[b, ks + d]
                    Mab[np.arange(n), np.arange(n) + d] = blk @ u
                    if d:
                        wv = Cn[a, ks + d] * Cn[b, ks]
                        Mab[np.arange(n) + d, np.arange(n)] = blk @ wv
                Y[a, b] = Mab
        xs = np.arange(dimj)
        gB = [xs[xs % kmod == c] for c in range(kmod)]
        full = np.concatenate([a * dimj + xs for a in range(r)])
        chg = np.concatenate([(xs - residues[a]) % kmod for a in range(r)])
        gRB = [full[chg == c] for c in range(kmod)]
        Bm = sum(Y[a, a] for a in range(r))
        RBm = np.block([[Y[a, b] for b in range(r)] for a in range(r)])
        for gi in gB:
            if len(gi) == 0:
                continue
            sub = 0.5 * (Bm[np.ix_(gi, gi)] + Bm[np.ix_(gi, gi)].T)
            A_ = arb_mat([[exact(x) for x in row] for row in sub]) / Z
            lo, hi, neg = entropy_interval_neg(A_, len(gi))
            S_B[0] += lo; S_B[1] += hi; negB += neg; DB += len(gi)
        for gi in gRB:
            if len(gi) == 0:
                continue
            sub = 0.5 * (RBm[np.ix_(gi, gi)] + RBm[np.ix_(gi, gi)].T)
            A_ = arb_mat([[exact(x) for x in row] for row in sub]) / Z
            lo, hi, neg = entropy_interval_neg(A_, len(gi))
            S_RB[0] += lo; S_RB[1] += hi; negRB += neg; DRB += len(gi)
    Zl = float(Z.lower())
    eB = epsB / Zl * (1 + 1e-12) + negB
    eRB = epsRB / Zl * (1 + 1e-12) + negRB

    def corr(e, D):
        if not (np.isfinite(e) and 0 <= e < 1 / np.e):
            raise ValueError(f"error budget {e!r} outside the range of the continuity bound")
        return arb(e) * (arb(D) / arb(e)).log() if e > 0 else arb(0)

    cB, cRB = corr(eB, DB), corr(eRB, DRB)
    ln2 = arb(2).log()
    lo = (S_B[0] - cB.upper() - S_RB[1] - cRB.upper()) / ln2
    hi = (S_B[1] + cB.upper() - S_RB[0] + cRB.upper()) / ln2
    if verbose:
        print(f"  budget: eps_B={eB:.2e} eps_RB={eRB:.2e} (negative mass {negB:.1e}, {negRB:.1e}); "
              f"corrections {float(cB.mid()):.2e}, {float(cRB.mid()):.2e} nats")
    return lo.lower(), hi.upper()


def entropy_interval_neg(A, dim):
    """certify.entropy_interval plus an upper bound on the negative eigenvalue mass of A (A need not be PSD)."""
    from certify import entropy_interval as ei
    mid = np.array([[float(A[i, j].mid()) for j in range(dim)] for i in range(dim)])
    mid = 0.5 * (mid + mid.T)
    _, Vv = np.linalg.eigh(mid)
    Vb = arb_mat([[exact(Vv[i, j]) for j in range(dim)] for i in range(dim)])
    Ap = Vb.transpose() * A * Vb
    R = Vb.transpose() * Vb
    eps2, del2 = arb(0), arb(0)
    diag = []
    for i in range(dim):
        for j in range(dim):
            if i == j:
                diag.append(Ap[i, i]); del2 += abs(R[i, i] - 1).upper() ** 2
            else:
                eps2 += abs(Ap[i, j]).upper() ** 2; del2 += abs(R[i, j]).upper() ** 2
    eps = eps2.sqrt().upper(); delta = del2.sqrt().upper()
    neg = 0.0
    for di in diag:
        lo_l = (di.lower() - eps) / (1 - delta)  # most negative possible eigenvalue
        if lo_l < 0:
            neg += -float(arb(lo_l).lower())
    lo, hi, _, _ = ei(A, dim)
    return lo, hi, neg * (1 + 1e-12)



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
    channel, p, path, M, kmod, res = sys.argv[1], Fraction(sys.argv[2]), sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), [int(x) for x in sys.argv[6].split(",")]
    lp, lz = exact_lambdas(channel, p)
    import hashlib
    print("argv:", " ".join(sys.argv))
    print("input sha256:", hashlib.sha256(open(path, "rb").read()).hexdigest())
    print("python-flint", flint.__version__, "prec", flint.ctx.prec, "| extended-precision integration (CERT_LD):", LD,
          "| step scalar in working precision: True | input checks: True")
    C = np.load(path)
    t = time.time()
    lo, hi = certify_big(C, M, lp, lz, res, kmod)
    olo, ohi = outward(lo, hi)
    print(f"{channel} p={p} M={M} K={C.shape[1]-1}: Ic in [{olo!r}, {ohi!r}] bits (outward-rounded floats; {time.time()-t:.0f}s)")
    print("CERTIFIED POSITIVE" if lo > 0 else "NOT CERTIFIED")
