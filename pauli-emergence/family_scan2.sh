#!/bin/bash
# mod-3 codes for r = 0.5 and 0.25 (M = 60), then certification of each threshold state.
cd "$(dirname "$0")"
for r in 0.5 0.25; do
  tag=tpr${r}_k3_M60
  p0=$(python3 -c "
import numpy as np
from scipy.optimize import brentq
H=lambda ps: -sum(x*np.log2(x) for x in ps if x>0)
f=lambda p: (1-p, p/(2+$r), p/(2+$r), p*$r/(2+$r))
print(round(brentq(lambda p: 1-H(f(p)),1e-4,0.49)-0.002,4))")
  python3 -u threshold_M.py --channel pr$r --M 60 --K 60 --p0 $p0 --step 0.001 --tol 0.0001 --kmod 3 --res 0,1 --mean 6 --width 4 --tag $tag > runs/$tag.log 2>&1
  p=$(python3 -c "
import json
for l in open('runs/$tag.log'):
    if 'RESULT' in l: d=json.loads(l); print(d['threshold_lo'])")
  frac=$(python3 -c "from fractions import Fraction; print(Fraction('$p'))")
  cp runs/${tag}_best.npy runs/state_pr${r}_k3_M60_thr.npy
  python3 -u certify_big.py pr$r $frac runs/state_pr${r}_k3_M60_thr.npy 60 3 0,1 > runs/certB_pr${r}_k3_M60.log 2>&1
done
