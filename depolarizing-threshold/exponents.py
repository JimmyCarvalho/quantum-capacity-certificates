"""Overlap-decay exponents of coherent 'blobs' seen by Bob (B) and by the environment (E,E') of the
emergent thermal attenuator (eta, nu=1/(2 eta)). Convention: vacuum quadrature variance 1/2, a = (x+ip)/sqrt2.
For equal-covariance Gaussian states displaced by d, in normal modes (Williamson) the Chernoff overlap is
Q = exp(-sum_i |d'_i|^2 (sqrt(N_i+1)-sqrt(N_i))^2) and the fidelity |<.|.>|-type exponent uses
F = exp(-sum_i |d'_i|^2 / (2 N_i + 1)) (Uhlmann fidelity of displaced thermal states, squared form)."""
import numpy as np
from scipy.linalg import sqrtm
from scipy.optimize import brentq

Om1 = np.array([[0, 1], [-1, 0]])

def williamson(V):
    n = V.shape[0] // 2
    Om = np.kron(np.eye(n), Om1)
    w = np.linalg.eigvals(1j * Om @ V)
    nus = np.sort(np.abs(w.real))[::2]
    # symplectic S with V = S diag(nu) S^T : use standard construction
    Vh = np.real(sqrtm(V)); Vhi = np.linalg.inv(Vh)
    A = Vhi @ Om @ Vhi  # antisymmetric
    # real Schur form of antisymmetric A
    from scipy.linalg import schur
    T, O = schur(A, output='real')
    # reorder blocks to [[0,a],[-a,0]] with a>0
    P = np.zeros_like(O)
    for i in range(n):
        blk = T[2*i:2*i+2, 2*i:2*i+2]
        if blk[0, 1] < 0:
            P[:, 2*i] = O[:, 2*i+1]; P[:, 2*i+1] = O[:, 2*i]
        else:
            P[:, 2*i] = O[:, 2*i]; P[:, 2*i+1] = O[:, 2*i+1]
    a = np.array([abs(T[2*i, 2*i+1]) for i in range(n)])
    nu = 1/a
    S = Vh @ P @ np.diag(np.repeat(1/np.sqrt(nu), 2))
    assert np.allclose(S @ np.diag(np.repeat(nu, 2)) @ S.T, V), 'Williamson V'
    assert np.allclose(S @ Om @ S.T, Om) or np.allclose(S @ Om @ S.T, -Om), 'symplectic'
    return nu, S  # V = S diag(nu) S^T, S symplectic

def exponents(eta, nu_env=None):
    nu = 1/(2*eta) if nu_env is None else nu_env
    v = nu + 0.5; c = np.sqrt(v*v - 0.25)
    NB = (1 - eta) * nu
    # displacement of unit |Delta alpha| along x: Delta x = sqrt2
    dx = np.sqrt(2)
    xiB_ch = eta * (np.sqrt(NB+1) - np.sqrt(NB))**2
    xiB_F = eta / (2*NB + 1)
    # E E' covariance
    Z = np.diag([1, -1]); I = np.eye(2)
    x = (1-eta)/2 + eta*v; y = np.sqrt(eta)*c
    V = np.block([[x*I, y*Z], [y*Z, v*I]])
    nus, S = williamson(V)
    out = {}
    for name, ang in [("x", 0.0), ("p", np.pi/2)]:
        d = np.array([-np.sqrt(1-eta)*dx*np.cos(ang), -np.sqrt(1-eta)*dx*np.sin(ang), 0, 0])
        dp = np.linalg.solve(S, d)
        ch = 0; F = 0
        for i in range(2):
            Ni = nus[i] - 0.5
            amp2 = (dp[2*i]**2 + dp[2*i+1]**2)/2  # |d'|^2 as complex amplitude
            ch += amp2 * (np.sqrt(Ni+1) - np.sqrt(max(Ni, 0)))**2
            F += amp2 / (2*Ni + 1)
        out[name] = (ch, F, nus)
    return xiB_ch, out["x"][0], out["p"][0], xiB_F, out["x"][1], nus

for eta in [0.70, 1/np.sqrt(2), 0.72, 0.7294452, 0.74, 0.76]:
    b, ex, ep, bF, eF, nus = exponents(eta)
    print(f"eta={eta:.6f}  Chernoff: B={b:.6f} E(x)={ex:.6f} E(p)={ep:.6f} | 'F': B={bF:.6f} E={eF:.6f} | E-modes nu={np.round(nus,4)}")
f = lambda e: exponents(e)[0] - exponents(e)[1]
print("Chernoff crossing eta =", brentq(f, 0.55, 0.95))
g = lambda e: exponents(e)[3] - exponents(e)[4]
print("F crossing eta =", brentq(g, 0.55, 0.95))
