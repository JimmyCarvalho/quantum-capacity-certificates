import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys, json
import numpy as np
from scipy.optimize import minimize
from gausschan import EmergentChannel, charge_blocks, coherent_info
from optimize import supports

eta = float(sys.argv[1]); k = int(sys.argv[2]); residues = [int(x) for x in sys.argv[3].split(",")]
Nin, D = 200, 300
ch = EmergentChannel(eta, Nin, D)
sup = supports(residues, k, Nin)
blocks = charge_blocks(residues, k, D)

def comb(params):
    mu, sig = params[0], abs(params[1])
    qs = np.r_[params[2:], 0.0]
    w = np.exp(qs) / np.exp(qs).sum()
    C = np.zeros((len(residues), Nin + 1))
    for i, s in enumerate(sup):
        a = np.exp(-((s - mu) ** 2) / (4 * sig**2))
        C[i, s] = np.sqrt(w[i]) * a / np.linalg.norm(a)
    return C

best = None
for mu0 in [20, 40, 60, 80]:
    for sig0 in [2.0, 4.0, 7.0]:
        x0 = np.r_[mu0, sig0, np.zeros(len(residues) - 1)]
        r = minimize(lambda x: -1e6 * coherent_info(ch, comb(x), blocks=blocks), x0, method="Nelder-Mead",
                     options={"maxiter": 400, "xatol": 1e-4, "fatol": 1e-9})
        if best is None or r.fun < best.fun:
            best = r
        print(json.dumps({"mu0": mu0, "sig0": sig0, "x": np.round(r.x, 4).tolist(), "Ic": -r.fun / 1e6}), flush=True)
print("BEST", json.dumps({"eta": eta, "k": k, "res": residues, "x": best.x.tolist(), "Ic": -best.fun / 1e6}))
np.save(f"runs/comb_k{k}_r{'-'.join(map(str,residues))}_eta{eta}.npy", comb(best.x))
