"""Rigorous (Arb ball arithmetic) lower bound on Ic(rho, T_D o G) for the emergent channel of a Pauli channel."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "depolarizing-threshold"))
from fractions import Fraction
import time
import numpy as np
from certify import certify


def exact_TG(channel, p):
    p = Fraction(p)
    if channel == "bb84":
        lp, lz = 1 - 2 * p, (1 - 2 * p) ** 2
    elif channel == "2pauli":
        lp, lz = 1 - 2 * p, 1 - 4 * p
    elif channel == "depol":
        lp = lz = 1 - 4 * p
    elif channel == "agn":  # additive Gaussian noise, p = sigma^2: A_G o L_T with G = 1 + sigma^2, T = 1/G
        return 1 / (1 + p), 1 + p
    elif channel.startswith("th"):  # thermal attenuator with nu thermal photons, p = 1 - eta (exact rationals)
        nu, eta = Fraction(channel[2:]), 1 - p
        G = (1 - eta) * (nu + Fraction(1, 2)) + (1 + eta) / 2
        return eta / G, G
    else:
        raise ValueError(channel)
    tau = lp * lp / lz
    Y = (1 - lp * lp) / (2 * lz)
    G = (2 * Y + tau + 1) / 2
    return tau / G, G


if __name__ == "__main__":
    channel, p, path, res, k, Nin, D = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
    T, G = exact_TG(channel, Fraction(p))
    C = np.load(path)[:, : Nin + 1]
    t = time.time()
    lo, hi = certify(C, [int(x) for x in res.split(",")], k, None, Nin, D, T=T, G=G)
    print(f"{channel} p={p} T={T} G={G}")
    print("Ic in [", lo.str(15, radius=False), ",", hi.str(15, radius=False), "] bits;", round(time.time() - t), "s")
    print("CERTIFIED POSITIVE" if lo > 0 else "NOT CERTIFIED")
