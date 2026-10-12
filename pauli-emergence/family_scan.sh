#!/bin/bash
# Certified-threshold map for (p_X, p_Y, p_Z) = p (1, 1, r)/(2 + r), M = 60, between 2-Pauli (r=0) and depolarizing (r=1).
cd "$(dirname "$0")"
for spec in "0.25 2 0,1" "0.5 2 0,1" "0.75 2 0,1" "0.75 3 0,1"; do
  set -- $spec
  r=$1; kmod=$2; res=$3
  p0=$(python3 -c "
import numpy as np
from scipy.optimize import brentq
H=lambda ps: -sum(x*np.log2(x) for x in ps if x>0)
f=lambda p: (1-p, p/(2+$r), p/(2+$r), p*$r/(2+$r))
print(round(brentq(lambda p: 1-H(f(p)),1e-4,0.49)-0.002,4))")
  python3 -u threshold_M.py --channel pr$r --M 60 --K 60 --p0 $p0 --step 0.001 --tol 0.0001 --kmod $kmod --res $res --mean 6 --width 4 --tag tpr${r}_k${kmod}_M60 > runs/tpr${r}_k${kmod}_M60.log 2>&1
done
