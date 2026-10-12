"""Emergent Gaussian channel of the qubit depolarizing channel (Ahmed et al., arXiv:2609.39747).

G_eta = A_G o L_T with T = 2 eta^2 / (1 + eta), G = (1 + eta) / (2 eta): a thermal attenuator
with transmissivity eta and environment photon number 1/(2 eta).

The channel is phase covariant, so it maps the band |m><m+d| into the same band d.
For each band we store W[d] with  G(|m><m+d|) = sum_a W[d][a, m] |a><a+d|  (output truncated to a+d <= D).
Output levels above D are collapsed into one flag state |perp> (a channel post-processing, so by data
processing the truncated coherent information is a lower bound on the true one).
"""
import numpy as np
from scipy.special import gammaln

LN2 = np.log(2.0)


def _lchoose(n, k):
    return gammaln(n + 1.0) - gammaln(k + 1.0) - gammaln(n - k + 1.0)


def loss_band(T, Nin, d):
    """Matrix L[d][m', m]: pure loss on band d, input m = 0..Nin-d, output m' = m - k."""
    n = Nin - d + 1
    m = np.arange(n)[None, :]
    mp = np.arange(n)[:, None]
    k = m - mp
    mask = k >= 0
    kk = np.where(mask, k, 0)
    with np.errstate(divide="ignore"):
        lt = np.log(T)
        l1t = np.log1p(-T) if T < 1 else -np.inf
    val = 0.5 * (_lchoose(m, kk) + _lchoose(m + d, kk)) + (m + 0.5 * d - kk) * lt + kk * l1t
    out = np.where(mask, np.exp(np.where(mask, val, -np.inf)), 0.0)
    return out


def amp_band(G, Nin, D, d):
    """Matrix A[d][a, j]: quantum-limited amplifier on band d, input j = 0..Nin-d, output a = j + k <= D - d."""
    nout = D - d + 1
    nin = Nin - d + 1
    j = np.arange(nin)[None, :]
    a = np.arange(nout)[:, None]
    k = a - j
    mask = k >= 0
    kk = np.where(mask, k, 0)
    lg = np.log(G)
    lr = np.log((G - 1.0) / G) if G > 1 else -np.inf
    val = 0.5 * (_lchoose(j + kk, kk) + _lchoose(j + d + kk, kk)) - (j + 0.5 * d + 1.0) * lg + kk * lr
    return np.where(mask, np.exp(np.where(mask, val, -np.inf)), 0.0)


class EmergentChannel:
    def __init__(self, eta, Nin, D, T=None, G=None):
        assert D >= Nin
        self.eta = eta
        self.T = 2 * eta**2 / (1 + eta) if T is None else T
        self.G = (1 + eta) / (2 * eta) if G is None else G
        self.Nin, self.D = Nin, D
        self.W = [amp_band(self.G, Nin, D, d) @ loss_band(self.T, Nin, d) for d in range(Nin + 1)]
        # trace of the kept (un-flagged) part of the output of |m><m|
        self.kept = self.W[0].sum(axis=0)

    def apply(self, X):
        """P G(X) P for an operator X on levels 0..Nin; returns (D+1)x(D+1)."""
        D, Nin = self.D, self.Nin
        Y = np.zeros((D + 1, D + 1), dtype=X.dtype)
        idx_out = np.arange(D + 1)
        idx_in = np.arange(Nin + 1)
        for d in range(Nin + 1):
            Wd = self.W[d]
            xu = X[idx_in[: Nin + 1 - d], idx_in[: Nin + 1 - d] + d]
            if xu.any():  # residue codes leave most bands exactly zero
                Y[idx_out[: D + 1 - d], idx_out[: D + 1 - d] + d] = Wd @ xu
            if d:
                xl = X[idx_in[: Nin + 1 - d] + d, idx_in[: Nin + 1 - d]]
                if xl.any():
                    Y[idx_out[: D + 1 - d] + d, idx_out[: D + 1 - d]] = Wd @ xl
        return Y

    def adjoint(self, Z):
        """Adjoint of X -> P G(X) P; Z is (D+1)x(D+1), returns (Nin+1)x(Nin+1)."""
        D, Nin = self.D, self.Nin
        X = np.zeros((Nin + 1, Nin + 1), dtype=Z.dtype)
        io = np.arange(D + 1)
        ii = np.arange(Nin + 1)
        for d in range(Nin + 1):
            Wd = self.W[d]
            zu = Z[io[: D + 1 - d], io[: D + 1 - d] + d]
            if zu.any():
                X[ii[: Nin + 1 - d], ii[: Nin + 1 - d] + d] = Wd.T @ zu
            if d:
                zl = Z[io[: D + 1 - d] + d, io[: D + 1 - d]]
                if zl.any():
                    X[ii[: Nin + 1 - d] + d, ii[: Nin + 1 - d]] = Wd.T @ zl
        return X


def _entropy_and_log(M, floor=1e-300):
    w, V = np.linalg.eigh(M)
    w = np.clip(w, floor, None)
    S = -np.sum(w * np.log(w))
    L = (V * np.log(w)) @ V.conj().T
    return S, L, w


def _block_entropy_and_log(M, blocks):
    """Entropy and matrix log of M, which is block diagonal on the index sets `blocks` (None = dense)."""
    if blocks is None:
        S, L, _ = _entropy_and_log(M)
        return S, L
    L = np.zeros_like(M)
    S = 0.0
    for b in blocks:
        Sb, Lb, _ = _entropy_and_log(M[np.ix_(b, b)])
        S += Sb
        L[np.ix_(b, b)] = Lb
    return S, L


def charge_blocks(residues, k, D):
    """Index sets on which B and RB are block diagonal when codeword i lives on n = residues[i] (mod k).
    Phase covariance of the channel makes the output commute with exp(2 pi i (n_B - a_R) / k)."""
    if k is None or k <= 1:
        return None
    n = np.arange(D + 1)
    bB = [n[n % k == c] for c in range(k)]
    r = len(residues)
    full = np.concatenate([i * (D + 1) + n for i in range(r)])
    charge = np.concatenate([(n - residues[i]) % k for i in range(r)])
    bRB = [full[charge == c] for c in range(k)]
    return bB, bRB


def coherent_info(ch, C, grad=False, blocks=None):
    """Coherent information (bits) of rho = sum_i |phi_i><phi_i| (phi_i = rows of C, normalised jointly)
    through the truncated-with-flag channel. Optionally the gradient w.r.t. C (Euclidean, before normalisation).
    `blocks` = charge_blocks(...) exploits the residue structure of the codewords."""
    nrm = np.sqrt(np.sum(np.abs(C) ** 2))
    Phi = C / nrm
    r, D = Phi.shape[0], ch.D
    Y = [[None] * r for _ in range(r)]
    flag = np.zeros((r, r), dtype=Phi.dtype)
    for i in range(r):
        for j in range(i, r):
            X = np.outer(Phi[i], Phi[j].conj())
            Y[i][j] = ch.apply(X)
            flag[i, j] = np.trace(X) - np.trace(Y[i][j])
            if j != i:
                Y[j][i] = Y[i][j].conj().T
                flag[j, i] = np.conj(flag[i, j])
    B = sum(Y[i][i] for i in range(r))
    fB = np.real(np.trace(flag))
    RB = np.block(Y)
    bB, bRB = (None, None) if blocks is None else blocks
    SB_main, LB = _block_entropy_and_log(B, bB)
    SB = SB_main - (fB * np.log(fB) if fB > 0 else 0.0)
    SRB_main, LRB = _block_entropy_and_log(RB, bRB)
    SF, LF, _ = _entropy_and_log(flag)
    SRB = SRB_main + SF
    Ic = (SB - SRB) / LN2
    if not grad:
        return Ic
    I = np.eye(D + 1)
    one_minus_kept = np.eye(ch.Nin + 1) - ch.adjoint(I)
    lfB = np.log(fB) if fB > 0 else 0.0
    ZB = ch.adjoint(LB + I) + (lfB + 1.0) * one_minus_kept
    g = np.zeros_like(Phi)
    n1 = D + 1
    for i in range(r):
        g[i] -= ZB @ Phi[i]
        for j in range(r):
            Lij = LRB[i * n1:(i + 1) * n1, j * n1:(j + 1) * n1] + (I if i == j else 0.0)
            Zij = ch.adjoint(Lij) + (LF[i, j] + (1.0 if i == j else 0.0)) * one_minus_kept
            g[i] += Zij @ Phi[j]
    # d(-Tr Y log Y) = -Tr[(log Y + 1) dY]; dIc/dPhi_i^* = -ZB phi_i + sum_j Z_ij phi_j ; real params get factor 2
    g = 2.0 * g / LN2
    # chain rule through normalisation Phi = C / |C|
    gC = (g - np.real(np.vdot(Phi, g)) * Phi) / nrm
    return Ic, gC
