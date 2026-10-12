#!/bin/bash
# After the M = 180 2-Pauli certificate, certify the M = 240 2-Pauli state at p = 251/2000 (extended precision).
cd "$(dirname "$0")"
until grep -q "CERTIF" runs/certB_2p3_M180thr_LD.log; do sleep 30; done
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 CERT_LD=1 CERT_PROCS=4 python3 -u certify_big.py 2pauli 251/2000 runs/state_2p3_M240_p0.1255.npy 240 3 0,1 > runs/certB_2p3_M240_LD.log 2>&1
