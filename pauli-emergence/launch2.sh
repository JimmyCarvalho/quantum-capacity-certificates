#!/bin/bash
cd "$(dirname "$0")"
nohup python3 -u threshold_M.py --channel depol --M 120 --K 90 --p0 0.065 --step 0.0005 --kmod 3 --res 0,1 --init runs/cert_dep_M60_thr.npy --tag tdep_M120 > runs/tdep_M120.log 2>&1 &
nohup python3 -u threshold_M.py --channel depol --M 90 --K 80 --p0 0.065 --step 0.0005 --kmod 3 --res 0,1 --init runs/cert_dep_M60_thr.npy --tag tdep_M90 > runs/tdep_M90.log 2>&1 &
nohup python3 -u threshold_M.py --channel bb84 --M 240 --K 110 --p0 0.1235 --step 0.0005 --init runs/state_bb_M120_thr.npy --tag tbb_M240 > runs/tbb_M240.log 2>&1 &
