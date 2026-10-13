"""Thermal attenuator E_{eta,nu}: where Q > 0 is certified (values from CERTIFICATES.md)."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

S2_AGN = 0.4111875  # Q > 0 certified for the additive Gaussian noise channel up to this variance

# (nu, smallest certified eta) from direct certificates
CERTIFIED = [(0.25, 0.626), (0.5, 0.6965), (1.0, 0.78), (2.0, 0.8587), (3.0, 0.8955), (5.0, 0.932375)]
MELE = (1.0, 0.7841)


def g(n):
    n = np.asarray(n, dtype=float)
    return (n + 1) * np.log2(n + 1) - n * np.log2(n)


def eta_gauss(n):  # thermal inputs give Q_G > 0 iff eta > eta_gauss
    r = 2 ** g(n)
    return r / (1 + r)


def eta_ad(n):  # antidegradable (Q = 0) for eta <= eta_ad
    return (n + 0.5) / (n + 1)


def eta_agn(n):  # corollary of the AGN certificate: Q > 0 for eta >= eta_agn
    return (n + 1) / (n + 1 + S2_AGN)


def frac(n, eta):  # fraction of the open window (eta_ad, eta_gauss) certified
    return (eta_gauss(n) - eta) / (eta_gauss(n) - eta_ad(n))


nu = np.linspace(0.02, 10, 500)
fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, 4.2))

a.fill_between(nu, 0.5, eta_ad(nu), color="0.85", label="antidegradable (Q = 0)")
a.plot(nu, eta_gauss(nu), "k-", lw=1.2, label="thermal (Gaussian) inputs work above")
a.plot(nu, eta_agn(nu), color="C0", ls="--", lw=1.2, label="Q > 0 above (AGN corollary, this work)")
a.plot(*zip(*CERTIFIED), "o", color="C3", label="Q > 0 certified directly (this work)")
a.plot(*MELE, "s", mfc="none", color="C2", ms=7, label="Mele et al. 2026 (certified)")
a.set_xlabel("thermal photons in the environment, ν")
a.set_ylabel("transmissivity η")
a.set_xlim(0, 8)
a.set_ylim(0.55, 1.0)
a.legend(fontsize=7, loc="lower right")
a.set_title("Thermal attenuator: certified Q > 0")

fa = np.clip(frac(nu, eta_agn(nu)), 0, None)
b.plot(nu, 100 * fa, color="C0", ls="--", lw=1.2, label="AGN corollary")
b.plot([x for x, _ in CERTIFIED], [100 * frac(x, e) for x, e in CERTIFIED], "o-", color="C3", label="direct certificates")
b.plot(MELE[0], 100 * frac(*MELE), "s", mfc="none", color="C2", ms=7, label="Mele et al. 2026")
b.set_xlabel("thermal photons in the environment, ν")
b.set_ylabel(r"% of the open window $(\eta_{AD}, \eta_G)$ certified")
b.set_xlim(0, 10)
b.set_ylim(0, 60)
b.legend(fontsize=7.5, loc="upper right")
b.set_title("How much of the unknown region is closed")

fig.tight_layout()
fig.savefig("attenuator_map.png", dpi=150)
