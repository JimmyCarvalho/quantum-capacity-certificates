#!/bin/bash
# Thermal attenuator with nu = 1 thermal photon (Mele et al., arXiv:2607.27449: Q > 0 certified down to
# eta = 0.7841, i.e. loss p = 0.2159; Gaussian inputs stop at eta = 0.8; antidegradable for eta <= 0.75).
# Residue families at low energy; p = 1 - eta increases with the noise.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
O="--channel th1 --Nin 80 --D 240 --p0 0.205 --step 0.004 --tol 0.0005 --maxiter 400"
lane() {
  for spec in "$@"; do
    IFS=: read k res mean width <<< "$spec"
    tag=th1_k${k}_$(echo $res | tr ',' '-')_m${mean}
    python3 -u threshold_limit.py $O --k $k --res $res --mean $mean --width $width --tag $tag > runs/$tag.log 2>&1
  done
}
lane "2:0,1:3:2" "3:0,1:4:2" "4:0,2:6:3" "3:0,1:10:4" &
lane "2:0,1:8:3" "4:0,1:5:3" "3:0,1,2:5:2" "5:0,2:8:3" &
wait
