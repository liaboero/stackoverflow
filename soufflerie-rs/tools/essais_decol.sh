#!/bin/sh
# Essais de la loi de paroi et de la viscosité turbulente sur le décrochage
# (aile infinie NACA 0012, Re = 6·10⁶). Usage (depuis soufflerie-rs) :
#   tools/essais_decol.sh <dossier> <mailles 64|128> "<α…>" "<essai…>"
# Une ligne CSV par (essai, α) dans <dossier>/essais.csv ; AILE_NUCAP passe.
set -u
OUT=$1; L=$2; ALPHAS=$3; ESSAIS=$4
case $L in 64) R=2 ;; 128) R=3 ;; *) echo "mailles : 64 ou 128"; exit 1 ;; esac
mkdir -p "$OUT"
X=target/release/examples/aile_infinie.exe
[ -f "$OUT/essais.csv" ] || echo "essai,L,alpha,Cz_boite,sd,Cx_boite,Cz_surface,Cm,x_decol_ext,x_decol_ext_normale,duree_s" > "$OUT/essais.csv"
for e in $ESSAIS; do
  for a in $ALPHAS; do
    m=1; [ "$(awk "BEGIN{print ($a > 14)}")" = 1 ] && m=4
    t=$((12 + m - 1))
    log="$OUT/${e}_${L}_$a.log"
    while [ "$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits | head -1)" -gt 30 ]; do sleep 60; done
    AILE_MOY=$m AILE_NREL=$((4 * m)) AILE_DECOL=1 "$X" 16 $R 6e6 "$a" "$t" - 60 12 "$e" > "$log" 2>&1
    s=$(grep 'relevés sur' "$log")
    czb=$(echo "$s" | sed -n 's/.*bilan (marge [0-9.]*) Cz \([-0-9.]*\) ± \([0-9.]*\), Cx.*/\1/p')
    sdb=$(echo "$s" | sed -n 's/.*bilan (marge [0-9.]*) Cz \([-0-9.]*\) ± \([0-9.]*\), Cx.*/\2/p')
    cxb=$(echo "$s" | sed -n 's/.*bilan (marge [0-9.]*) Cz [-0-9.]* ± [0-9.]*, Cx \([-0-9.]*\).*/\1/p')
    czs=$(echo "$s" | sed -n 's/.*surface Cz \([-0-9.]*\) ±.*/\1/p')
    cm=$(echo "$s" | sed -n 's/.*Cm \([-0-9.]*\) ± \([0-9.]*\) .*/\1/p')
    xe=$(grep 'décollement x/c (vitesse' "$log" | sed -n 's/.*extrados \([0-9.]*\) ;.*/\1/p')
    xeh=$(grep 'au-dessus de la paroi' "$log" | sed -n 's/.*extrados \([0-9.]*\) ·.*/\1/p')
    du=$(sed -n 's/.* · \([0-9]*\) s$/\1/p' "$log")
    echo "$e,$L,$a,$czb,$sdb,$cxb,$czs,$cm,$xe,$xeh,$du" | tee -a "$OUT/essais.csv"
    sleep 10
  done
done
