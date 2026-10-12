#!/bin/bash
# Thermal attenuator map, second pass: finer continuation where the first pass collapsed (nu = 3, 5),
# larger cutoff for nu = 0.25, and symmetry-free continuation for every nu (p = 1 - eta).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
gen() {  # nu init s0 step tag
  python3 -u agn_general_thr.py --channel th$1 --init $2 --s0 $3 --step $4 --Nin 160 --D 360 --maxiter 600 \
    --noise 0.005 --seed 1 --thr 1e-9 --tag $5 > runs/$5.log 2>&1
}
fine() {  # nu init p0 step tag Nin D
  python3 -u threshold_limit.py --channel th$1 --k 3 --res 0,1 --Nin $6 --D $7 --p0 $3 --step $4 --tol 0.0001 \
    --maxiter 500 --init $2 --tag $5 > runs/$5.log 2>&1
}
(fine 0.25 runs/thnu_0.25_best.npy 0.370 0.001 thnu_0.25_big 220 440; gen 0.25 runs/thnu_0.25_big_best.npy 0.371 0.0005 thgen_0.25) &
(fine 3 runs/thnu_3_best.npy 0.1013 0.0005 thnu_3_fine 160 360; gen 3 runs/thnu_3_fine_best.npy 0.1015 0.0005 thgen_3; gen 2 runs/thnu_2_best.npy 0.1398 0.0005 thgen_2) &
(gen 0.5 runs/thnu_0.5_best.npy 0.302 0.0005 thgen_0.5; fine 5 runs/thnu_5_best.npy 0.0626 0.0003 thnu_5_fine 160 360) &
wait
