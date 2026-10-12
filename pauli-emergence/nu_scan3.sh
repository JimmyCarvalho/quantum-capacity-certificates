#!/bin/bash
# Thermal attenuator, third pass: nu = 1/4 and nu = 5 were pressing against N_in = 160; larger cutoff.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 -u agn_general_thr.py --channel th0.25 --init runs/th0.25_best.npy --s0 0.3745 --step 0.0005 \
  --Nin 260 --D 520 --maxiter 600 --noise 0.003 --seed 3 --thr 1e-9 --tag thgen_0.25_big > runs/thgen_0.25_big.log 2>&1 &
python3 -u threshold_limit.py --channel th5 --k 3 --res 0,1 --Nin 260 --D 520 --p0 0.0674 --step 0.0003 --tol 0.0001 \
  --maxiter 500 --init runs/th5_best.npy --tag thnu_5_big > runs/thnu_5_big.log 2>&1 &
wait
