#!/bin/bash
# AGN channel: residue-code families with a large Fock cutoff (Nin = 300, D = 600), two lanes.
# (The Nin = 140 scan truncated the high-energy families: mod-4 {0,2} wanted ~60-75 photons.)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
O="--Nin 300 --D 600 --tol 0.0005 --maxiter 300"
lane() {
  for spec in "$@"; do
    IFS=: read k res p0 init <<< "$spec"
    tag=agnL_k${k}_$(echo $res | tr ',' '-')
    if [ -n "$init" ]; then ini="--init $init"; else ini="--mean 40 --width 10"; fi
    python3 -u threshold_limit.py --channel agn --k $k --res $res --p0 $p0 --step 0.004 $O $ini --tag $tag > runs/$tag.log 2>&1
  done
}
lane "4:0,2:0.394:runs/agn_k4_0-2_N140_best.npy" "4:0,1:0.37:" "6:0,3:0.37:" &
lane "3:0,1:0.4:runs/agn_k3_01_best.npy" "5:0,2:0.37:" "5:0,1:0.37:" "2:0,1:0.37:" &
wait
