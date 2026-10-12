#!/bin/bash
# Re-run every certify_big certificate of CERTIFICATES.md with the final code (one code version for all),
# then the new 2-Pauli M = 240 point. Bands in parallel; one BLAS thread per process.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 CERT_PROCS=4
run() { log=$1; shift; python3 -u certify_big.py "$@" > runs/$log.log 2>&1; }
run certB_2p_M60      2pauli 605469/5000000 runs/cert_2p_M60_p0.1210938.npy 60 2 0,1
run certB_bb_M60      bb84 304297/2500000 runs/cert_bb_M60_thr.npy 60 2 0,1
run certB_dep_M60     depol 652187/10000000 runs/cert_dep_M60_thr.npy 60 3 0,1
run certB_pr0.25_k3_M60 pr0.25 17827/80000 runs/state_pr0.25_k3_M60_thr.npy 60 3 0,1
run certB_pr0.5_k3b_M60 pr0.5 3359/16000 runs/state_pr0.5_k3b_M60_thr.npy 60 3 0,1
run certB_pr0.75_k3_M60 pr0.75 16127/80000 runs/state_pr0.75_k3_M60_thr.npy 60 3 0,1
run certB_bb_M120     bb84 619531/5000000 runs/state_bb_M120_thr.npy 120 2 0,1
run certB_2p_M120     2pauli 308047/2500000 runs/state_2p_M120_thr.npy 120 2 0,1
run certB_dep_M120    depol 1051/16000 runs/state_dep_M120_thr.npy 120 3 0,1
run certB_bb3_M120    bb84 156211/1250000 runs/state_bb3_M120_thr.npy 120 3 0,1
run certB_2p3_M120    2pauli 154961/1250000 runs/state_2p3_M120_thr.npy 120 3 0,1
CERT_LD=1 run certB_2p3_M240_LD 2pauli 251/2000 runs/state_2p3_M240_p0.1255.npy 240 3 0,1
run certB_dep_M240    depol 41207/625000 runs/state_dep_M240_thr.npy 240 3 0,1
run certB_bb3_M240    bb84 63/500 runs/state_bb3_M240_q0.126.npy 240 3 0,1
CERT_LD=1 run certB_2p3_M180thr_LD 2pauli 2003/16000 runs/state_2p3_M180_thr.npy 180 3 0,1
echo ALL DONE > runs/recertify_all.done
