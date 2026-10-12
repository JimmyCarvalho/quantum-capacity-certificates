#!/bin/bash
cd "$(dirname "$0")"
nohup python3 -u threshold_M.py --channel 2pauli --M 240 --K 110 --p0 0.1232 --step 0.0005 --init runs/state_2p_M120_thr.npy --tag t2p_M240 > runs/t2p_M240.log 2>&1 &
nohup python3 -u threshold_M.py --channel depol --M 240 --K 110 --p0 0.06565 --step 0.00025 --kmod 3 --res 0,1 --init runs/state_dep_M120_thr.npy --tag tdep_M240 > runs/tdep_M240.log 2>&1 &
