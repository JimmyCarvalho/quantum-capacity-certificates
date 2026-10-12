"""Three-component cat inputs on the AGN channel (usage: python3 agn_cat3.py 0.36,0.38)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "depolarizing-threshold"))
import numpy as np
from math import lgamma
from scipy.optimize import minimize_scalar, minimize
from gausschan import EmergentChannel, charge_blocks, coherent_info
from limit_codes import emergent_TG
Nin, D = 140, 320
n = np.arange(Nin + 1)
def cat(alpha, w):
    # three-component cat codewords: support n = 0 and 1 (mod 3), Poisson amplitudes alpha^n / sqrt(n!)
    logamp = n * np.log(alpha) - 0.5 * np.array([lgamma(k + 1) for k in n])
    amp = np.exp(logamp - logamp.max())
    c0 = np.where(n % 3 == 0, amp, 0.0); c1 = np.where(n % 3 == 1, amp, 0.0)
    c0 /= np.linalg.norm(c0); c1 /= np.linalg.norm(c1)
    return np.array([np.sqrt(w) * c0, np.sqrt(1 - w) * c1])
bl = charge_blocks([0, 1], 3, D)
for s2 in [float(x) for x in sys.argv[1].split(",")]:
    T, G = emergent_TG("agn", s2); ch = EmergentChannel(0.5, Nin, D, T=T, G=G)
    f = lambda x: -coherent_info(ch, cat(x[0], x[1]), blocks=bl)
    best = None
    for a0 in [2.5, 3.5, 4.5, 5.5]:
        r = minimize(f, [a0, 0.45], method="Nelder-Mead", bounds=[(0.5, 9), (0.05, 0.95)], options={"xatol": 1e-5, "fatol": 1e-12, "maxiter": 400})
        if best is None or r.fun < best.fun: best = r
    print(s2, "best cat Ic = %.4e at alpha = %.4f (mean n = %.2f), w0 = %.4f" % (-best.fun, best.x[0], best.x[0]**2, best.x[1]), flush=True)
