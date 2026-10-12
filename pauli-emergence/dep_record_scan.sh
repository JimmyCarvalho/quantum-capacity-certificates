#!/bin/bash
# Depolarizing record point (eta = 0.7294452, per-Pauli p = 0.0676387; Ahmed et al. report Ic ~ 1e-12 in float).
# Symmetry-free rank-2 optimisation from the mod-3 plateau state (mean 56.5 photons), small perturbations.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for seed in 1 2; do
  python3 -u agn_general_thr.py --channel depol --init ../depolarizing-threshold/runs/mod3_01_w3.5.npy \
    --s0 0.0676387 --step 0.0001 --Nin 220 --D 400 --maxiter 800 --noise 0.003 --seed $seed --thr 1e-12 \
    --tag deprec_s$seed > runs/deprec_s$seed.log 2>&1 &
done
wait
