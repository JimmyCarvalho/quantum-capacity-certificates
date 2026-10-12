"""Depolarizing record point, high precision: mod-3 codes (residues 0, 1) with a Gaussian envelope
|c(n)|^2 ~ exp(-(n - mu)^2 / (2 sig^2)), positive amplitudes, weights q, 1 - q.
usage: python3 dep_hp_envelope.py eta mus sigs q"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import json
import sys
from fractions import Fraction
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "depolarizing-threshold"))
import flint
import numpy as np

import hp_ic
from certify import ArbChannel

NIN, D = 220, 320
flint.ctx.prec = 128


def code(mu, sig, q):
    n = np.arange(NIN + 1)
    env = np.exp(-((n - mu) ** 2) / (4 * sig ** 2))
    cs = []
    for x in (0, 1):
        c = np.where(n % 3 == x, env, 0.0)
        cs.append(c / np.linalg.norm(c))
    return np.array([np.sqrt(q) * cs[0], np.sqrt(1 - q) * cs[1]])


if __name__ == "__main__":
    eta = Fraction(sys.argv[1])
    ch = ArbChannel(eta, NIN, D)
    q = float(sys.argv[4])
    for mu in [float(x) for x in sys.argv[2].split(",")]:
        for sig in [float(x) for x in sys.argv[3].split(",")]:
            v = float(hp_ic.ic(ch, code(mu, sig, q)).mid())
            print(json.dumps({"eta": str(eta), "mu": mu, "sig": sig, "q": q, "Ic": v}), flush=True)
