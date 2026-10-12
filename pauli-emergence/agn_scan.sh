#!/bin/bash
# Threshold (in sigma^2) of residue-code families on the additive Gaussian noise channel (Gaussian limit of BB84).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
python3 -u threshold_limit.py --channel agn --k 3 --res 0,1 --Nin 140 --D 320 --p0 0.4 --step 0.002 --tol 0.0002 --init runs/agn_k3_sigma2_0.4.npy --tag agn_k3_01 > runs/agn_k3_01.log 2>&1
for spec in 4:0,2 4:0,1 5:0,1 5:0,2 2:0,1 6:0,3; do
  k=${spec%%:*}; res=${spec#*:}; tag=agn_k${k}_$(echo $res | tr ',' '-')
  python3 -u threshold_limit.py --channel agn --k $k --res $res --Nin 140 --D 320 --p0 0.37 --step 0.004 --tol 0.0005 --mean 12 --width 4 --tag $tag > runs/$tag.log 2>&1
done
