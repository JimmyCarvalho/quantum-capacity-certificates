#!/bin/bash
# Depolarizing channel, Gaussian limit (thermal attenuator of Ahmed et al.): continue the M = 240 mod-3 code
# in the limit, then break the residue symmetry (as for the AGN channel). Per-Pauli p; record is 0.0676387.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 -u threshold_limit.py --channel depol --k 3 --res 0,1 --Nin 300 --D 600 --p0 0.0655 --step 0.0005 \
  --tol 0.0001 --maxiter 300 --init runs/state_dep_M240_thr.npy --tag limdep_k3 > runs/limdep_k3.log 2>&1
python3 -u agn_general_thr.py --channel depol --init runs/limdep_k3_best.npy --s0 0.066 --step 0.0002 \
  --Nin 200 --D 400 --noise 0.02 --thr 1e-9 --tag limdep_gen > runs/limdep_gen.log 2>&1
