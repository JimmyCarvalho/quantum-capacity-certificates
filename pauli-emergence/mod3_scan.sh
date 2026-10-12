#!/bin/bash
# mod-3 codes for 2-Pauli and BB84 at M = 60 (compare with parity: 0.1210938 and 0.1217188), and r=0.5 retry.
cd "$(dirname "$0")"
while pgrep -f "family_scan2.s[h]" > /dev/null; do sleep 30; done
python3 -u threshold_M.py --channel 2pauli --M 60 --K 60 --p0 0.117 --step 0.0005 --tol 0.00005 --kmod 3 --res 0,1 --mean 6 --width 4 --tag t2p3_M60 > runs/t2p3_M60.log 2>&1
python3 -u threshold_M.py --channel bb84 --M 60 --K 60 --p0 0.117 --step 0.0005 --tol 0.00005 --kmod 3 --res 0,1 --mean 6 --width 4 --tag tbb3_M60 > runs/tbb3_M60.log 2>&1
python3 -u threshold_M.py --channel pr0.5 --M 60 --K 60 --p0 0.196 --step 0.0005 --tol 0.0001 --kmod 3 --res 0,1 --mean 6 --width 4 --tag tpr0.5_k3b_M60 > runs/tpr0.5_k3b_M60.log 2>&1
