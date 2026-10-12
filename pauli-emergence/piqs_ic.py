"""Exact coherent information of symmetric-subspace inputs through N^{(x)M}, N a z-phase-covariant unital
Pauli channel (lam_perp, lam_par), via qutip.piqs (permutation-invariant Dicke-basis evolution)."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from qutip import piqs
from scipy.sparse.linalg import expm_multiply

def rates(lam_perp, lam_par):
    g = -np.log(lam_par) / 2.0
    gD = -2.0 * np.log(lam_perp) + np.log(lam_par)
    return g, gD

class SymChannel:
    def __init__(self, M, lam_perp, lam_par):
        g, gD = rates(lam_perp, lam_par)
        self.M = M
        # qutip's Cython stores the rates as float32: build unit-rate generators and scale in float64
        LE = piqs.Dicke(M, emission=1.0).liouvillian().data.as_scipy()
        LP = piqs.Dicke(M, pumping=1.0).liouvillian().data.as_scipy()
        LD = piqs.Dicke(M, dephasing=1.0).liouvillian().data.as_scipy()
        self.L = (g * LE + g * LP + gD * LD).tocsr()
        self.nds = piqs.num_dicke_states(M)
        self.js = piqs.j_vals(M)[::-1]  # piqs orders blocks from largest j
        self.deg = [piqs.state_degeneracy(M, j) for j in self.js]
        # sanity: block layout
        self.blocks_start = np.cumsum([0] + [int(round(2 * j + 1)) for j in self.js])[:-1]
        # the block entries form an invariant subspace of the Liouvillian: restrict to it
        idx = []
        for st, j in zip(self.blocks_start, self.js):
            n = int(round(2 * j + 1))
            for c in range(n):
                for r in range(n):
                    idx.append((st + c) * self.nds + (st + r))  # column-stacked index of entry (r, c)
        self.idx = np.array(sorted(idx))
        Lfull = self.L
        self.Lsub = Lfull[self.idx][:, self.idx].tocsr()
        leak = abs(Lfull[:, self.idx]).sum() - abs(self.Lsub).sum()
        assert leak < 1e-9 * abs(self.Lsub).sum(), "block subspace not invariant"

    def op_from_sym(self, A):
        """A: (M+1)x(M+1) operator on Sym^M in basis |k excitations> -> Dicke-basis nds x nds matrix."""
        R = np.zeros((self.nds, self.nds), dtype=complex)
        R[: self.M + 1, : self.M + 1] = A
        return R

    def apply(self, ops):
        """ops: list of nds x nds matrices -> list of outputs (column-stacked superoperator convention)."""
        V = np.stack([o.reshape(-1, order="F")[self.idx] for o in ops], axis=1)
        W = expm_multiply(self.Lsub, V)
        outs = []
        for i in range(W.shape[1]):
            full = np.zeros(self.nds * self.nds, dtype=complex)
            full[self.idx] = W[:, i]
            outs.append(full.reshape(self.nds, self.nds, order="F"))
        return outs

    def blocks(self, R):
        out = []
        for s, j in zip(self.blocks_start, self.js):
            n = int(round(2 * j + 1))
            out.append(R[s:s + n, s:s + n])
        return out

def ent(A):
    w = np.linalg.eigvalsh(0.5 * (A + A.conj().T))
    w = w[w > 1e-300]
    return -np.sum(w * np.log2(w))

def coherent_info(ch, vecs):
    """rho = sum_a |v_a><v_a| (v_a unnormalised, total trace 1), vectors on Sym^M (k-excitation basis)."""
    r = len(vecs)
    pairs = [(a, b) for a in range(r) for b in range(r)]
    outs = ch.apply([ch.op_from_sym(np.outer(vecs[a], vecs[b].conj())) for a, b in pairs])
    Y = {pr: ch.blocks(o) for pr, o in zip(pairs, outs)}
    Ic = 0.0
    for bi, (j, dj) in enumerate(zip(ch.js, ch.deg)):
        B = sum(Y[(a, a)][bi] for a in range(r))
        RB = np.block([[Y[(a, b)][bi] for b in range(r)] for a in range(r)])
        # piqs stores the total weight of the d_j copies; the log d_j terms cancel between B and RB
        Ic += ent(B) - ent(RB)
    return Ic

if __name__ == "__main__":
    import sys
    sys.path.insert(0, ".")
    from piqs_check import pauli_channel_bruteforce
    from math import comb
    rng = np.random.default_rng(1)
    q = 0.11
    s = 1 - 2 * q
    probs_y_frame = [(1 - q) ** 2, q * (1 - q), q * (1 - q), q * q]
    for M in [4, 5]:
        ch = SymChannel(M, s, s * s)
        vecs = [rng.normal(size=M + 1) + 1j * rng.normal(size=M + 1) for _ in range(2)]
        nrm = np.sqrt(sum(np.vdot(v, v).real for v in vecs)); vecs = [v / nrm for v in vecs]
        Ic_piqs = coherent_info(ch, vecs)
        # brute force
        dim = 2 ** M
        D = np.zeros((dim, M + 1))
        for x in range(dim):
            k = bin(x).count("1"); D[x, k] = 1 / np.sqrt(comb(M, k))
        V = [D @ v for v in vecs]
        rhoRA = np.zeros((2 * dim, 2 * dim), dtype=complex)
        psi = np.concatenate(V)  # |0>_R V0 + |1>_R V1
        outRB = np.zeros((2 * dim, 2 * dim), dtype=complex)
        for a in range(2):
            for b in range(2):
                outRB[a*dim:(a+1)*dim, b*dim:(b+1)*dim] = pauli_channel_bruteforce(np.outer(V[a], V[b].conj()), probs_y_frame, M)
        outB = outRB[:dim, :dim] + outRB[dim:, dim:]
        Ic_bf = ent(outB) - ent(outRB)
        print(f"M={M}: piqs Ic={Ic_piqs:.12f}  brute force Ic={Ic_bf:.12f}")
