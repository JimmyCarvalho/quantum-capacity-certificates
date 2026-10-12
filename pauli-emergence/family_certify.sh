#!/bin/bash
# Wait for each family threshold search and certify its best state with certify_big.py.
cd "$(dirname "$0")"
for spec in "0.25 2 0,1" "0.5 2 0,1" "0.75 2 0,1" "0.75 3 0,1"; do
  set -- $spec
  r=$1; kmod=$2; res=$3
  tag=tpr${r}_k${kmod}_M60
  until grep -q RESULT runs/$tag.log 2>/dev/null; do sleep 30; done
  p=$(python3 -c "
import json
for l in open('runs/$tag.log'):
    if 'RESULT' in l: d=json.loads(l); print(d['threshold_lo'])")
  frac=$(python3 -c "from fractions import Fraction; print(Fraction('$p'))")
  cp runs/${tag}_best.npy runs/state_pr${r}_k${kmod}_M60_thr.npy
  python3 -u certify_big.py pr$r $frac runs/state_pr${r}_k${kmod}_M60_thr.npy 60 $kmod $res > runs/certB_pr${r}_k${kmod}_M60.log 2>&1
done
