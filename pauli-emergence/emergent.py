"""Emergent single-mode Gaussian channel of a qubit Pauli channel restricted to the symmetric subspace
near the pole of axis c (Holstein-Primakoff limit), from first and second moments:
    X = diag(lam_a, lam_b) / sqrt(lam_c),   Y = diag(1 - lam_a^2, 1 - lam_b^2) / (2 lam_c)
(vacuum quadrature variance 1/2). For the depolarizing channel this is the thermal attenuator of
Ahmed et al. (arXiv:2609.39747): transmissivity eta, environment 1/(2 eta)."""
import numpy as np
from scipy.optimize import brentq

def lambdas(pI, pX, pY, pZ):
    return {"x": 1 - 2 * (pY + pZ), "y": 1 - 2 * (pX + pZ), "z": 1 - 2 * (pX + pY)}

def emergent(probs, pole):
    lam = lambdas(*probs)
    a, b = [k for k in "xyz" if k != pole]
    X = np.diag([lam[a], lam[b]]) / np.sqrt(lam[pole])
    Y = np.diag([1 - lam[a] ** 2, 1 - lam[b] ** 2]) / (2 * lam[pole])
    return X, Y

def g(x):
    return 0.0 if x <= 0 else (x + 1) * np.log2(x + 1) - x * np.log2(x)

def canonical(X, Y):
    tau = np.linalg.det(X)
    dY = np.sqrt(np.linalg.det(Y))
    if abs(tau - 1) < 1e-14:
        return tau, dY, None
    N = dY / abs(1 - tau) - 0.5
    return tau, dY, N

def gaussian_Q(X, Y):
    """Infinite-energy Gaussian coherent information (Holevo-Werner, canonical form)."""
    tau, dY, N = canonical(X, Y)
    if N is None:
        return np.log2(1 / (np.e * dY))
    return np.log2(tau / abs(1 - tau)) - g(N)

def antidegradable(X, Y):
    """Phase-insensitive canonical forms: attenuator tau<1 is AD iff tau <= (N+1/2)/(N+1)
    (Lami-Khatri-Adesso-Wilde 2019); additive noise AD iff sigma^2 >= 1/2."""
    tau, dY, N = canonical(X, Y)
    if N is None:
        return dY >= 0.5
    if tau < 1:
        return tau <= (N + 0.5) / (N + 1)
    return None  # amplifier: handled separately

def threshold(fam, pole, f, lo, hi):
    return brentq(lambda p: f(*emergent(fam(p), pole)), lo, hi, xtol=1e-14)

dep = lambda p: (1 - 3 * p, p, p, p)            # per-Pauli p ; total error 3p
indep = lambda q: ((1 - q) ** 2, q * (1 - q), q * q, q * (1 - q))
twop = lambda p: (1 - 2 * p, p, 0.0, p)

if __name__ == "__main__":
    print("check depolarizing: eta=0.7294452 ->", emergent(dep(0.75 * (1 - 0.7294452) / 3), "z"))
    for name, fam, poles, hashing in [("depolarizing (per-Pauli p)", dep, "z", None),
                                       ("independent X-Z / BB84 (q)", indep, "yz", None),
                                       ("2-Pauli (p)", twop, "yz", None)]:
        for pole in poles:
            try:
                tG = threshold(fam, pole, gaussian_Q, 1e-4, 0.2)
            except ValueError:
                tG = float("nan")
            print(f"{name:28s} pole={pole}: Gaussian-input threshold {tG:.6f}", end="")
            X, Y = emergent(fam(tG if tG == tG else 0.1), pole)
            print(f"   canonical(tau, sqrtdetY, N) at threshold = {np.round(canonical(X, Y)[0],5), np.round(canonical(X,Y)[1],5), canonical(X,Y)[2]}")
    # hashing thresholds for reference
    H = lambda ps: -sum(x * np.log2(x) for x in ps if x > 0)
    for name, fam in [("depolarizing", dep), ("indep X-Z", indep), ("2-Pauli", twop)]:
        print(f"hashing {name}: {brentq(lambda p: 1 - H(fam(p)), 1e-4, 0.2):.6f}")
    # ceilings from antidegradability of the emergent channel
    print("depolarizing ceiling (eta=1/sqrt2): per-Pauli p =", (1 - 1 / np.sqrt(2)) / 4, " total =", 0.75 * (1 - 1 / np.sqrt(2)))
    print("indep X-Z, y-pole ceiling sigma^2=1/2: q =", (1 - 1 / np.sqrt(2)) / 2)
