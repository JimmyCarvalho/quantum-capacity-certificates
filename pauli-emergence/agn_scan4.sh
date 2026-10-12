#!/bin/bash
# AGN channel beyond rank-2 residue codes: rank 3 (mod 3 {0,1,2}) and symmetry-free rank-2/3 inputs.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 -u threshold_limit.py --channel agn --k 3 --res 0,1,2 --Nin 300 --D 600 --p0 0.39 --step 0.004 --tol 0.0005 --maxiter 300 --mean 15 --width 5 --tag agnL_k3_0-1-2 > runs/agnL_k3_0-1-2.log 2>&1
python3 -u agn_general.py --init runs/agn_k3_sigma2_0.405.npy --sigma2 0.405 --rank 2 --Nin 200 --D 400 --noise 0.02 --tag agn_gen_r2_0.405 > runs/agn_gen_r2_0.405.log 2>&1
python3 -u agn_general.py --init runs/agn_k3_sigma2_0.405.npy --sigma2 0.405 --rank 3 --Nin 200 --D 400 --noise 0.02 --tag agn_gen_r3_0.405 > runs/agn_gen_r3_0.405.log 2>&1
