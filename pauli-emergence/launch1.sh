#!/bin/bash
cd "$(dirname "$0")"
nohup python3 -u limit_codes.py --channel 2pauli --k 4 --res 0,2 --ps 0.124,0.125,0.126,0.127,0.128,0.13 --Nin 220 --D 380 --init runs/2p_k4_p0.124.npy --tag 2p_k4b > runs/2p_k4b.log 2>&1 &
nohup python3 -u limit_codes.py --channel bb84 --k 4 --res 0,2 --ps 0.114,0.117,0.12,0.121,0.122,0.123 --Nin 180 --D 320 --mean 20 --width 4 --tag bb84_k4 > runs/bb84_k4.log 2>&1 &
nohup python3 -u threshold_M.py --channel 2pauli --M 60 --K 60 --p0 0.117 --tag t2p_M60 > runs/t2p_M60.log 2>&1 &
nohup python3 -u threshold_M.py --channel bb84 --M 60 --K 60 --p0 0.117 --tag tbb_M60 > runs/tbb_M60.log 2>&1 &
