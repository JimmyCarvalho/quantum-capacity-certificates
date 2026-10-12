#!/bin/bash
cd "$(dirname "$0")"
until grep -q RESULT runs/tbb3_M60.log 2>/dev/null; do sleep 30; done
nohup python3 -u threshold_M.py --channel 2pauli --M 120 --K 90 --p0 0.1225 --step 0.0005 --tol 0.00005 --kmod 3 --res 0,1 --init runs/t2p3_M60_best.npy --tag t2p3_M120 > runs/t2p3_M120.log 2>&1 &
nohup python3 -u threshold_M.py --channel bb84 --M 120 --K 90 --p0 0.1235 --step 0.0005 --tol 0.00005 --kmod 3 --res 0,1 --init runs/tbb3_M60_best.npy --tag tbb3_M120 > runs/tbb3_M120.log 2>&1 &
wait
