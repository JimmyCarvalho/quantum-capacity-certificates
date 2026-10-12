#!/bin/bash
# Thermal attenuator E_{eta,nu} for several nu: threshold in the loss p = 1 - eta for mod-3 {0,1} codes,
# started just inside the Gaussian region (eta_G = 2^g/(1+2^g)). Runs after th1_scan*.sh finish.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
while ps -eo args | grep -q "[t]h1_scan"; do sleep 60; done
run() {
  nu=$1
  p0=$(python3 -c "
import math
nu=$nu; g=(nu+1)*math.log2(nu+1)-nu*math.log2(nu); r=2**g
print(round(1-r/(1+r)-0.002, 4))")
  python3 -u threshold_limit.py --channel th$nu --k 3 --res 0,1 --Nin 120 --D 300 --p0 $p0 --step 0.002 --tol 0.0002 \
    --maxiter 400 --mean 12 --width 4 --tag thnu_${nu} > runs/thnu_${nu}.log 2>&1
}
(run 0.5; run 2; run 5) &
(run 0.25; run 3) &
wait
