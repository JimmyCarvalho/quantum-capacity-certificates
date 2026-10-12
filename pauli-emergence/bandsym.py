"""Band-restricted exact evolution for z-phase-covariant unital Pauli channels on permutation-invariant
operators of M qubits. The generator conserves delta = m - m', so each band is an invariant subspace of size
~M^2/4 (instead of ~M^3/6 for the full Dicke space). Same coefficients as dicke_lindblad (float64)."""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply


class Band:
    def __init__(self, M, delta):
        self.M, self.delta = M, delta
        jmax = M / 2
        js = [jmax - i for i in range(int(np.floor(jmax)) + 1)]
        self.js = [j for j in js if 2 * j + 1 - delta > 0]
        self.lens = [int(round(2 * j + 1 - delta)) for j in self.js]
        self.off = np.cumsum([0] + self.lens)[:-1]
        self.size = int(sum(self.lens))
        jj, mm = [], []
        for j, n in zip(self.js, self.lens):
            jj.append(np.full(n, j))
            mm.append(j - np.arange(n))  # m from j downwards; m1 = m - delta >= -j
        self.j = np.concatenate(jj)
        self.m = np.concatenate(mm)
        self.m1 = self.m - delta
        self.jpos = {j: i for i, j in enumerate(self.js)}

    def index(self, j, m, m1):
        j = np.asarray(j, float); m = np.asarray(m, float); m1 = np.asarray(m1, float)
        ok = (j <= self.M / 2 + 1e-9) & (j >= 0) & (np.abs(m) <= j + 1e-9) & (np.abs(m1) <= j + 1e-9) & (2 * j + 1 - self.delta > 1e-9)
        ij = np.rint(self.M / 2 - j).astype(int)
        ij = np.where(ok, ij, 0)
        # js list may have dropped small j at the end only, so position == ij when ok
        ok = ok & (ij < len(self.js))
        off = self.off[np.clip(ij, 0, len(self.off) - 1)]
        idx = off + np.rint(j - m).astype(int)
        return np.where(ok, idx, -1)

    def generator(self, yE, yP, yD):
        N = float(self.M)
        j, m, m1 = self.j, self.m, self.m1
        r = np.arange(self.size)
        rows, cols, vals = [r], [r], []
        with np.errstate(divide="ignore", invalid="ignore"):
            jp = np.where(j > 0, j, 1.0)
            deph = np.where(j <= 0, yD * N / 4, yD / 2 * (N / 2 - m * m1 * (N / 2 + 1) / jp / (jp + 1)))
            vals.append(-(yE / 2 * (N + m + m1) + yP / 2 * (N - m - m1) + deph))

            def add(J, A, B, coef):
                c = self.index(J, A, B)
                ok = (c >= 0) & np.isfinite(coef) & (coef != 0)
                rows.append(r[ok]); cols.append(c[ok]); vals.append(coef[ok])

            A, B = m + 1, m1 + 1
            J = j; Jp = np.where(J > 0, J, 1.0)
            add(J, A, B, np.where(J > 0, yE / 2 * np.sqrt(np.clip((J + A) * (J - A + 1) * (J + B) * (J - B + 1), 0, None)) * (N / 2 + 1) / (Jp * (Jp + 1)), 0.0))
            J = j + 1
            add(J, A, B, yE / 2 * np.sqrt(np.clip((J + A) * (J + A - 1) * (J + B) * (J + B - 1), 0, None)) * (N / 2 + J + 1) / (J * (2 * J + 1)))
            J = j - 1
            add(J, A, B, np.where(J + 1 > 0, yE / 2 * np.sqrt(np.clip((J - A + 1) * (J - A + 2) * (J - B + 1) * (J - B + 2), 0, None)) * (N / 2 - J) / ((J + 1) * (2 * J + 1)), 0.0))
            J = j + 1
            add(J, m, m1, yD / 2 * np.sqrt(np.clip((J**2 - m**2) * (J**2 - m1**2), 0, None)) * (N / 2 + J + 1) / (J * (2 * J + 1)))
            J = j - 1
            add(J, m, m1, yD / 2 * np.sqrt(np.clip(((J + 1) ** 2 - m**2) * ((J + 1) ** 2 - m1**2), 0, None)) * (N / 2 - J) / ((J + 1) * (2 * J + 1)))
            A, B = m - 1, m1 - 1
            J = j + 1
            add(J, A, B, yP / 2 * np.sqrt(np.clip((J - A - 1) * (J - A) * (J - B - 1) * (J - B), 0, None)) * (N / 2 + J + 1) / (J * (2 * J + 1)))
            J = j; Jp = np.where(J > 0, J, 1.0)
            add(J, A, B, np.where(J > 0, yP / 2 * np.sqrt(np.clip((J + A + 1) * (J - A) * (J + B + 1) * (J - B), 0, None)) * (N / 2 + 1) / (Jp * (Jp + 1)), 0.0))
            J = j - 1
            add(J, A, B, yP / 2 * np.sqrt(np.clip((J + A + 1) * (J + A + 2) * (J + B + 1) * (J + B + 2), 0, None)) * (N / 2 - J) / ((J + 1) * (2 * J + 1)))
        R = np.concatenate(rows); C = np.concatenate(cols); V = np.concatenate(vals)
        return coo_matrix((V, (R, C)), shape=(self.size, self.size)).tocsr()


def rates(lam_perp, lam_par):
    g = -np.log(lam_par) / 2.0
    return g, g, -2.0 * np.log(lam_perp) + np.log(lam_par)


class BandChannel:
    """Exact N^{(x)M} on operators |k><l| of Sym^M (k excitations), via band decomposition."""

    def __init__(self, M, lam_perp, lam_par):
        self.M = M
        self.rates = rates(lam_perp, lam_par)
        self.bands = {}

    def band(self, delta):
        if delta not in self.bands:
            b = Band(self.M, delta)
            b.L = b.generator(*self.rates)
            self.bands[delta] = b
        return self.bands[delta]

    def evolve(self, delta, coeffs):
        """coeffs: dict k -> amplitude for input operators |k><k+delta| (delta >= 0), several columns allowed:
        coeffs is an array (K+1-delta, ncols) of weights of |k><k+delta|. Returns array (band.size, ncols)."""
        b = self.band(delta)
        V = np.zeros((b.size, coeffs.shape[1]), dtype=complex)
        n0 = b.lens[0]
        V[: coeffs.shape[0], :] = coeffs  # input block j=M/2: m = M/2 - k at position k
        assert coeffs.shape[0] <= n0
        return expm_multiply(b.L, V)


def _ent(w):
    w = np.where(w > 1e-300, w, 1.0)
    return -np.sum(np.where(w < 1.0, w * np.log2(w), 0.0) * (w > 0))


def ent_bits(w):
    w = np.asarray(w).ravel()
    w = w[w > 1e-300]
    return -np.sum(w * np.log2(w))


def diag_input_ic(ch, p, return_parts=False):
    """Exact Ic of rho = sum_k p_k |k><k| (k = excitations, k <= K <= M) through N^{(x)M}."""
    K = len(p) - 1
    p = np.asarray(p, float)
    outs = {}
    for d in range(K + 1):
        coeff = np.zeros((K + 1 - d, K + 1 - d), dtype=complex)
        ks = np.arange(K + 1 - d)
        coeff[ks, ks] = np.sqrt(p[ks] * p[ks + d])
        outs[d] = ch.evolve(d, coeff)
    Ic = 0.0
    bands = {d: ch.band(d) for d in range(K + 1)}
    jlist = bands[0].js
    for jpos, j in enumerate(jlist):
        dimj = int(round(2 * j + 1))
        b0 = bands[0]
        o0 = b0.off[jpos]
        Bdiag = np.real(outs[0][o0:o0 + dimj, :].sum(axis=1))
        SB = ent_bits(Bdiag)
        # RB sectors c = a - k, a in [0, dimj), k in [0, K]
        cs = np.arange(-K, dimj)
        A = np.zeros((len(cs), K + 1, K + 1), dtype=complex)
        for d in range(K + 1):
            b = bands[d]
            if j not in b.jpos:
                continue
            jp = b.jpos[j]
            n = b.lens[jp]
            blk = outs[d][b.off[jp]:b.off[jp] + n, :]  # (a, k): value of <a| N(|k><k+d|) |a+d>
            for k in range(K + 1 - d):
                a = cs + k
                ok = (a >= 0) & (a < n)
                A[ok, k, k + d] = blk[a[ok], k]
                if d:
                    A[ok, k + d, k] = np.conj(blk[a[ok], k])
        w = np.linalg.eigvalsh(A)
        SRB = ent_bits(w)
        Ic += SB - SRB
    return Ic
