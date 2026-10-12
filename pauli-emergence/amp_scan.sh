#!/bin/bash
# Thermal amplifier (emergent 2-Pauli channel): symmetry-free rank-2 continuation in p from the mod-3 state.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 -u agn_general_thr.py --channel 2pauli --init runs/amp_2p3_p0.1275.npy --s0 0.1275 --step 0.00025 \
  --Nin 200 --D 400 --noise 0.02 --tag amp_gen_thr > runs/amp_gen_thr.log 2>&1
