#!/bin/bash
# Other residue-code families for BB84 at M = 60 (mod-3 {0,1} gives 0.1227813, parity 0.1217188).
cd "$(dirname "$0")"
part=$1
if [ "$part" = "a" ]; then specs="4:0,1 4:0,2 3:0,1,2"; else specs="5:0,1 5:0,2 4:0,1,2"; fi
for spec in $specs; do
  k=${spec%%:*}; res=${spec#*:}
  tag=tbbk${k}_$(echo $res | tr ',' '-')_M60
  python3 -u threshold_M.py --channel bb84 --M 60 --K 60 --p0 0.118 --step 0.001 --tol 0.0001 --kmod $k --res $res --mean 6 --width 4 --tag $tag > runs/$tag.log 2>&1
done
