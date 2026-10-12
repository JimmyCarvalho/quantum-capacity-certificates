import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys, json
import numpy as np
from gausschan import EmergentChannel, charge_blocks, coherent_info
from optimize import supports, gaussian_init
k=3; res=[0,1]
Nin, D = 160, 260
etas=[0.70,0.71,0.72,0.7294452,0.74,0.75]
chs={e:EmergentChannel(e,Nin,D) for e in etas}
sup=supports(res,k,Nin); bl=charge_blocks(res,k,D)
for mu in [6,10,15,20,25,30,40]:
    for sig in [1.5,2.5,4.0]:
        row=[]
        for e in etas:
            best=max(coherent_info(chs[e],gaussian_init(sup,Nin,mu,sig,weights=[q,1-q]),blocks=bl) for q in [0.4,0.45,0.5,0.55])
            row.append(best)
        print(f"mu={mu:3d} sig={sig:3.1f} "+" ".join(f"{x:+.2e}" for x in row), flush=True)
