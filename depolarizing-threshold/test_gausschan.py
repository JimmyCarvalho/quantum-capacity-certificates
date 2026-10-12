"""Checks of gausschan against exact Gaussian formulas and finite differences."""
import numpy as np
from gausschan import EmergentChannel, coherent_info

def g(x):
    return 0.0 if x <= 0 else (x + 1) * np.log2(x + 1) - x * np.log2(x)

def gaussian_ic_thermal(eta, nE, n):
    """Ic of thermal input (mean n) through thermal attenuator (eta, nE). Units: vacuum variance 1/2."""
    # R-A two-mode squeezed vacuum with local variance a = n + 1/2, correlation c
    a = n + 0.5; c = np.sqrt(a * a - 0.25)
    b = eta * a + (1 - eta) * (nE + 0.5)
    cc = np.sqrt(eta) * c
    # RB covariance [[a I, cc Z],[cc Z, b I]]; symplectic eigenvalues
    Delta = a * a + b * b - 2 * cc * cc
    det = (a * b - cc * cc) ** 2
    nu_p = np.sqrt((Delta + np.sqrt(Delta**2 - 4 * det)) / 2)
    nu_m = np.sqrt((Delta - np.sqrt(Delta**2 - 4 * det)) / 2)
    return g(b - 0.5) - g(nu_p - 0.5) - g(nu_m - 0.5)

eta = 0.75
ch = EmergentChannel(eta, Nin=120, D=260)
nE = 1 / (2 * eta)
# 1. vacuum -> thermal with mean (1-eta) nE
Y = ch.apply(np.diag(np.r_[1.0, np.zeros(120)]))
mean = np.sum(np.arange(261) * np.diag(Y))
print("vacuum out mean", mean, "expected", (1 - eta) * nE, "trace", np.trace(Y))
# 2. Gaussian Ic of thermal input: purification of thermal state = TMSV; rank = full but truncated
for n in [0.5, 2.0, 5.0]:
    lam = n / (n + 1)
    p = (1 - lam) * lam ** np.arange(121)
    C = np.diag(np.sqrt(p))  # rows = unnormalised phi_i = sqrt(p_i)|i>
    Ic = coherent_info(ch, C)
    print(f"thermal n={n}: Fock Ic={Ic:.10f}  Gaussian formula={gaussian_ic_thermal(eta, nE, n):.10f}")
# 3. gradient vs finite differences
rng = np.random.default_rng(0)
ch2 = EmergentChannel(0.8, Nin=30, D=80)
C = rng.normal(size=(2, 31)) * np.exp(-0.1 * np.arange(31))
Ic, gC = coherent_info(ch2, C, grad=True)
for (i, k) in [(0, 3), (1, 7), (0, 20)]:
    e = np.zeros_like(C); e[i, k] = 1e-6
    fd = (coherent_info(ch2, C + e) - coherent_info(ch2, C - e)) / 2e-6
    print("grad", i, k, gC[i, k], "fd", fd)
