#!/bin/bash
# AGN channel beyond sigma^2 = 0.407: symmetry-free continuation with a larger cutoff (lane A) and from the
# high-energy mod-4 {0,2} family with its symmetry broken (lane B).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 -u agn_general_thr.py --channel agn --init runs/agn_gen_sigma2_0.407.npy --s0 0.4075 --step 0.0005 \
  --Nin 300 --D 600 --maxiter 400 --noise 0.005 --seed 1 --thr 1e-9 --tag agnA > runs/agnA.log 2>&1 &
python3 -u agn_general_thr.py --channel agn --init runs/agnL_k4_0-2_best.npy --s0 0.4 --step 0.001 \
  --Nin 300 --D 600 --maxiter 400 --noise 0.01 --seed 2 --thr 1e-9 --tag agnB > runs/agnB.log 2>&1 &
wait
