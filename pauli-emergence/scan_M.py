import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import warnings; warnings.filterwarnings("ignore")
import sys, json, time
import numpy as np
from bandsym import BandChannel, diag_input_ic
ch_name, q, Nbar, K = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
Ms = [int(m) for m in sys.argv[5].split(",")]
if ch_name == "bb84":
    lp, lz = 1 - 2 * q, (1 - 2 * q) ** 2
else:
    lp, lz = 1 - 2 * q, 1 - 4 * q
x = Nbar / (1 + Nbar)
p = x ** np.arange(K + 1); p /= p.sum()
for M in Ms:
    t = time.time()
    Ic = diag_input_ic(BandChannel(M, lp, lz), p)
    print(json.dumps({"ch": ch_name, "q": q, "Nbar": Nbar, "K": K, "M": M, "Ic": Ic, "sec": round(time.time() - t, 1)}), flush=True)
