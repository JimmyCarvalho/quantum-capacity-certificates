"""Initial states for the high-photon search: the moderate-photon mod-3 branch state, shifted up by s
excitations per residue step (n -> n + 3s, shape kept) or stretched in the residue index m (n = 3m + r,
m -> lam * m, linear interpolation of the amplitudes, signs included).
usage: python3 stretch_init.py in.npy out.npy shift|stretch value Nin"""
import sys

import numpy as np

if __name__ == "__main__":
    src, out, mode, val, Nin = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), int(sys.argv[5])
    S = np.load(src)
    C = np.zeros((2, Nin + 1))
    for r in (0, 1):
        m_src = np.arange(r, S.shape[1], 3)
        a = S[r, m_src]
        m_new = np.arange(r, Nin + 1, 3)
        idx_new = (m_new - r) // 3
        if mode == "shift":
            j = idx_new - int(val)
            ok = (j >= 0) & (j < len(a))
            C[r, m_new[ok]] = a[j[ok]]
        else:
            x = idx_new / val
            C[r, m_new] = np.interp(x, np.arange(len(a)), a, right=0.0)
        C[r] *= np.sqrt(np.sum(S[r] ** 2)) / np.linalg.norm(C[r])
    n = np.arange(Nin + 1)
    print("mean photon number %.2f" % (np.sum(n * C ** 2) / np.sum(C ** 2)))
    np.save(out, C)
