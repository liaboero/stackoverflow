for c in "3 16" "3 17" "3 18" "3 19" "2 18"; do set -- $c
  while [ "$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits | head -1)" -gt 30 ]; do sleep 60; done
  AILE_MOY=4 AILE_NREL=16 AILE_DECOL=1 target/release/examples/aile_infinie.exe 16 $1 6e6 $2 15 - 60 > decrochage/sep_r$1_$2.log 2>&1
  grep -E 'décollement|relevés sur' decrochage/sep_r$1_$2.log; sleep 10
done
