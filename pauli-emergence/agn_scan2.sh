#!/bin/bash
# After agn_scan.sh: rank-3 mod-3 codes and symmetry-free inputs on the AGN channel.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
while pgrep -f "agn_scan.s[h]" > /dev/null; do sleep 30; done
python3 -u threshold_limit.py --channel agn --k 3 --res 0,1,2 --Nin 140 --D 320 --p0 0.37 --step 0.004 --tol 0.0005 --mean 14 --width 4 --tag agn_k3_0-1-2 > runs/agn_k3_0-1-2.log 2>&1
python3 -u agn_general.py --init runs/agn_k3_01_best.npy --sigma2 0.404 --rank 2 --tag agn_gen_r2_0.404 > runs/agn_gen_r2_0.404.log 2>&1
python3 -u agn_general.py --init runs/agn_k3_01_best.npy --sigma2 0.404 --rank 3 --tag agn_gen_r3_0.404 > runs/agn_gen_r3_0.404.log 2>&1
