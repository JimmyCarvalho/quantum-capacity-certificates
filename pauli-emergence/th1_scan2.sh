#!/bin/bash
# nu = 1 thermal attenuator: symmetry-free continuation from the mod-3 and mod-4 states (p = 1 - eta).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 -u agn_general_thr.py --channel th1 --init runs/th1_k3_0-1_m4_best.npy --s0 0.217 --step 0.0005 \
  --Nin 120 --D 300 --maxiter 600 --noise 0.005 --seed 1 --thr 1e-9 --tag th1gen_k3 > runs/th1gen_k3.log 2>&1 &
python3 -u threshold_limit.py --channel th1 --k 3 --res 0,1 --Nin 160 --D 360 --p0 0.2165 --step 0.0005 --tol 0.0001 \
  --maxiter 500 --init runs/th1_k3_0-1_m4_best.npy --tag th1_k3_big > runs/th1_k3_big.log 2>&1 &
wait
