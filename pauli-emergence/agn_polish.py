"""Multi-round L-BFGS polish of a mod-3 residue code on the AGN channel (double precision).
usage: python3 agn_polish.py sigma2 init.npy out.npy Nin D maxiter_per_round [rounds=6]"""
import os, sys, json
os.environ["OMP_NUM_THREADS"]="1"
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "depolarizing-threshold"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from gausschan import EmergentChannel, charge_blocks, coherent_info
from optimize import optimise, supports
from limit_codes import emergent_TG
s2, path, out, Nin, D, iters = float(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
rounds = int(sys.argv[7]) if len(sys.argv) > 7 else 6
T, G = emergent_TG("agn", s2)
ch = EmergentChannel(0.5, Nin, D, T=T, G=G)
bl = charge_blocks([0,1], 3, D); sup = supports([0,1], 3, Nin)
C = np.load(path)[:, :Nin+1]
Ic = coherent_info(ch, C, blocks=bl); print(json.dumps({"start": Ic}), flush=True)
for r in range(rounds):
    Ic, C = optimise(ch, sup, C, maxiter=iters, verbose=False, scale=1.0/max(abs(Ic),1e-12), blocks=bl)
    np.save(out, C); print(json.dumps({"round": r, "Ic": Ic}), flush=True)
