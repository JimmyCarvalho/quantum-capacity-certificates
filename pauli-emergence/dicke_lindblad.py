"""Float64 generator of local emission / pumping / dephasing on permutation-invariant operators of M qubits,
stored compactly as the blocks rho_j[m, m'] (Dicke basis, total weight of the d_j copies).

Coefficients follow Shammah et al., PRA 98, 063815 (2018) as implemented in qutip.piqs (_piqs.pyx), which
evaluates them in float32; here they are float64 (and exact formulas, usable with Arb in certify code).
A local Pauli channel that is phase covariant about z and unital, with transverse/longitudinal Bloch factors
(lam_perp, lam_par), equals exp(L) with emission = pumping = g = -ln(lam_par)/2 and dephasing
gD = -2 ln(lam_perp) + ln(lam_par) (gD < 0 is a valid linear map; the composite stays CPTP).
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply


class DickeSpace:
    def __init__(self, M):
        self.M = M
        jmax = M / 2
        self.js = [jmax - i for i in range(int(np.floor(jmax)) + 1)]  # largest j first
        self.dims = [int(round(2 * j + 1)) for j in self.js]
        self.off = np.cumsum([0] + [d * d for d in self.dims])[:-1]
        self.size = int(sum(d * d for d in self.dims))
        jj, mm, m1 = [], [], []
        for j, d in zip(self.js, self.dims):
            m = j - np.arange(d)
            mm.append(np.repeat(m, d))
            m1.append(np.tile(m, d))
            jj.append(np.full(d * d, j))
        self.j = np.concatenate(jj)
        self.m = np.concatenate(mm)
        self.m1 = np.concatenate(m1)
        self.jindex = {j: i for i, j in enumerate(self.js)}

    def index(self, j, m, m1):
        """Vectorised index of (j, m, m1); -1 where it does not exist."""
        j = np.asarray(j, float)
        m = np.asarray(m, float)
        m1 = np.asarray(m1, float)
        ok = (j >= (self.M % 2) / 2 - 1e-9) & (j <= self.M / 2 + 1e-9) & (np.abs(m) <= j + 1e-9) & (np.abs(m1) <= j + 1e-9)
        i_j = np.where(ok, np.rint(self.M / 2 - j), 0).astype(int)
        d = np.rint(2 * j + 1).astype(int)
        off = self.off[np.clip(i_j, 0, len(self.off) - 1)]
        idx = off + np.rint(j - m).astype(int) * d + np.rint(j - m1).astype(int)
        return np.where(ok, idx, -1)

    def generator(self, yE, yP, yD):
        N = float(self.M)
        j, m, m1 = self.j, self.m, self.m1
        rows, cols, vals = [], [], []
        r = np.arange(self.size)
        with np.errstate(divide="ignore", invalid="ignore"):
            jp = np.where(j > 0, j, 1.0)
            deph = np.where(j <= 0, yD * N / 4, yD / 2 * (N / 2 - m * m1 * (N / 2 + 1) / jp / (jp + 1)))
            g1 = -(yE / 2 * (N + m + m1) + yP / 2 * (N - m - m1) + deph)
            rows.append(r); cols.append(r); vals.append(g1)

            def add(jj, ma, m1a, coef):
                c = self.index(jj, ma, m1a)
                ok = (c >= 0) & np.isfinite(coef) & (coef != 0)
                rows.append(r[ok]); cols.append(c[ok]); vals.append(coef[ok])

            # 2: (j, m+1, m1+1)
            J, A, B = j, m + 1, m1 + 1
            Jp = np.where(J > 0, J, 1.0)
            add(J, A, B, np.where(J > 0, yE / 2 * np.sqrt(np.clip((J + A) * (J - A + 1) * (J + B) * (J - B + 1), 0, None)) * (N / 2 + 1) / (Jp * (Jp + 1)), 0.0))
            # 3: (j+1, m+1, m1+1)
            J = j + 1
            add(J, A, B, yE / 2 * np.sqrt(np.clip((J + A) * (J + A - 1) * (J + B) * (J + B - 1), 0, None)) * (N / 2 + J + 1) / (J * (2 * J + 1)))
            # 4: (j-1, m+1, m1+1)
            J = j - 1
            add(J, A, B, np.where(J + 1 > 0, yE / 2 * np.sqrt(np.clip((J - A + 1) * (J - A + 2) * (J - B + 1) * (J - B + 2), 0, None)) * (N / 2 - J) / ((J + 1) * (2 * J + 1)), 0.0))
            # 5: (j+1, m, m1)
            J = j + 1
            add(J, m, m1, yD / 2 * np.sqrt(np.clip((J**2 - m**2) * (J**2 - m1**2), 0, None)) * (N / 2 + J + 1) / (J * (2 * J + 1)))
            # 6: (j-1, m, m1)
            J = j - 1
            add(J, m, m1, yD / 2 * np.sqrt(np.clip(((J + 1) ** 2 - m**2) * ((J + 1) ** 2 - m1**2), 0, None)) * (N / 2 - J) / ((J + 1) * (2 * J + 1)))
            # 7: (j+1, m-1, m1-1)
            J, A, B = j + 1, m - 1, m1 - 1
            add(J, A, B, yP / 2 * np.sqrt(np.clip((J - A - 1) * (J - A) * (J - B - 1) * (J - B), 0, None)) * (N / 2 + J + 1) / (J * (2 * J + 1)))
            # 8: (j, m-1, m1-1)
            J = j
            Jp = np.where(J > 0, J, 1.0)
            add(J, A, B, np.where(J > 0, yP / 2 * np.sqrt(np.clip((J + A + 1) * (J - A) * (J + B + 1) * (J - B), 0, None)) * (N / 2 + 1) / (Jp * (Jp + 1)), 0.0))
            # 9: (j-1, m-1, m1-1)
            J = j - 1
            add(J, A, B, yP / 2 * np.sqrt(np.clip((J + A + 1) * (J + A + 2) * (J + B + 1) * (J + B + 2), 0, None)) * (N / 2 - J) / ((J + 1) * (2 * J + 1)))
        R = np.concatenate(rows); C = np.concatenate(cols); V = np.concatenate(vals)
        return coo_matrix((V, (R, C)), shape=(self.size, self.size)).tocsr()

    def from_sym(self, A):
        """Operator on Sym^M in the excitation basis |k> = |j=M/2, m=M/2-k> -> compact vector."""
        v = np.zeros(self.size, dtype=complex)
        d = self.dims[0]
        v[: d * d] = np.asarray(A, dtype=complex).reshape(-1)
        return v

    def blocks(self, v):
        return [v[o:o + d * d].reshape(d, d) for o, d in zip(self.off, self.dims)]


class PauliSym:
    """N^{(x)M} on permutation-invariant operators, N unital and phase covariant about z."""

    def __init__(self, M, lam_perp, lam_par):
        self.space = DickeSpace(M)
        g = -np.log(lam_par) / 2.0
        gD = -2.0 * np.log(lam_perp) + np.log(lam_par)
        self.L = self.space.generator(g, g, gD)

    def apply(self, ops):
        V = np.stack([self.space.from_sym(o) for o in ops], axis=1)
        W = expm_multiply(self.L, V)
        return [self.space.blocks(W[:, i]) for i in range(W.shape[1])]


def ent(A):
    w = np.linalg.eigvalsh(0.5 * (A + A.conj().T))
    w = w[w > 1e-300]
    return -np.sum(w * np.log2(w))


def coherent_info(ch, vecs):
    r = len(vecs)
    pairs = [(a, b) for a in range(r) for b in range(r)]
    outs = ch.apply([np.outer(vecs[a], np.conj(vecs[b])) for a, b in pairs])
    Y = dict(zip(pairs, outs))
    Ic = 0.0
    for bi in range(len(ch.space.js)):
        B = sum(Y[(a, a)][bi] for a in range(r))
        RB = np.block([[Y[(a, b)][bi] for b in range(r)] for a in range(r)])
        Ic += ent(B) - ent(RB)
    return Ic
