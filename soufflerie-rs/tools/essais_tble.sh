#!/bin/sh
# Loi de paroi TBLE (SOUFFLERIE_TBLE) sur l'aile infinie NACA 0012,
# Re = 6·10⁶. Usage (depuis soufflerie-rs) :
#   tools/essais_tble.sh <dossier> <mailles 64|128> "<α…>" "<modèle…>"
# modèle : base (Spalding), eq, gp, gpc. Une ligne CSV par calcul dans
# <dossier>/tble.csv (mêmes colonnes que tools/essais_decol.sh).
set -u
OUT=$1; L=$2; ALPHAS=$3; MODES=$4
case $L in 64) R=2 ;; 128) R=3 ;; *) echo "mailles : 64 ou 128"; exit 1 ;; esac
mkdir -p "$OUT"
X=target/release/examples/aile_infinie.exe
[ -f "$OUT/tble.csv" ] || echo "modele,L,alpha,Cz_boite,sd,Cx_boite,Cz_surface,Cm,x_decol_ext,x_decol_ext_normale,duree_s" > "$OUT/tble.csv"
for e in $MODES; do
  for a in $ALPHAS; do
    m=1; [ "$(awk "BEGIN{print ($a > 14)}")" = 1 ] && m=4
    t=$((12 + m - 1))
    log="$OUT/${e}_${L}_$a.log"
    if [ "$e" = base ]; then unset SOUFFLERIE_TBLE; else export SOUFFLERIE_TBLE=$e; fi
    AILE_MOY=$m AILE_NREL=$((4 * m)) AILE_DECOL=1 "$X" 16 $R 6e6 "$a" "$t" - 60 12 > "$log" 2>&1
    s=$(grep 'relevés sur' "$log")
    czb=$(echo "$s" | sed -n 's/.*bilan (marge [0-9.]*) Cz \([-0-9.]*\) ± \([0-9.]*\), Cx.*/\1/p')
    sdb=$(echo "$s" | sed -n 's/.*bilan (marge [0-9.]*) Cz \([-0-9.]*\) ± \([0-9.]*\), Cx.*/\2/p')
    cxb=$(echo "$s" | sed -n 's/.*bilan (marge [0-9.]*) Cz [-0-9.]* ± [0-9.]*, Cx \([-0-9.]*\).*/\1/p')
    czs=$(echo "$s" | sed -n 's/.*surface Cz \([-0-9.]*\) ±.*/\1/p')
    cm=$(echo "$s" | sed -n 's/.*Cm \([-0-9.]*\) ± \([0-9.]*\) .*/\1/p')
    xe=$(grep 'décollement x/c (vitesse' "$log" | sed -n 's/.*extrados \([0-9.]*\) ;.*/\1/p')
    xeh=$(grep 'au-dessus de la paroi' "$log" | sed -n 's/.*extrados \([0-9.]*\) ·.*/\1/p')
    du=$(sed -n 's/.* · \([0-9]*\) s$/\1/p' "$log")
    echo "$e,$L,$a,$czb,$sdb,$cxb,$czs,$cm,$xe,$xeh,$du" | tee -a "$OUT/tble.csv"
  done
done
