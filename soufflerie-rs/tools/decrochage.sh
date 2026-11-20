#!/bin/sh
# Balayage d'incidence du NACA 0012 jusqu'au décrochage, Re = 6·10⁶ :
# 2D (panneaux + couche limite, couplage Newton) puis aile d'envergure
# infinie 3D (`aile_infinie`, 64 puis 128 mailles par corde, ny = 60).
# Usage (depuis soufflerie-rs) :
#   tools/decrochage.sh [dossier] [incidences 3D] [mailles : 64, 128 ou « 64 128 »]
# Défauts : decrochage, « 0 2 4 6 8 10 12 14 15 16 17 18 19 20 », « 64 128 ».
# T : durée en cordes (12) ; MOY : durée de moyenne en cordes (1 ; 4 au-delà
# de 14°, l'écoulement décollé pouvant être instationnaire) ; PAUSE entre les
# calculs (s, 20) ; SANS2D=1 : 3D seulement. Écrit <dossier>/2d.txt, <dossier>/3d_<L>.csv et les
# journaux <dossier>/3d_<L>_<α>.log.
# GPU partagé : vérifier avec nvidia-smi qu'aucun autre calcul ne tourne.
set -u
OUT=${1:-decrochage}
ALPHAS=${2:-"0 2 4 6 8 10 12 14 15 16 17 18 19 20"}
MAILLES=${3:-"64 128"}
T=${T:-12}
PAUSE=${PAUSE:-20}
mkdir -p "$OUT"
cargo build --release -p lbm --example decrochage2d --example aile_infinie 2>/dev/null
X=target/release/examples
X2=
[ -x "$X/aile_infinie.exe" ] && X2=.exe

[ "${SANS2D:-}" = "" ] && "$X/decrochage2d$X2" 0012 6e6 0:20:0.5 | tee "$OUT/2d.txt"

for L in $MAILLES; do
  case $L in 64) R=2 ;; 128) R=3 ;; *) echo "mailles : 64 ou 128"; exit 1 ;; esac
  echo "alpha,Cz_boite,Cx_boite,Cz_surface,Cm_surface,sd_Cz_boite,sd_Cm,x_decol_ext,x_decol_int,x_decol_ext_normale,x_decol_int_normale,duree_s" > "$OUT/3d_$L.csv"
  for a in $ALPHAS; do
    m=${MOY:-1}
    [ "${MOY:-}" = "" ] && [ "$(awk "BEGIN{print ($a > 14)}")" = 1 ] && m=4
    t=$(awk "BEGIN{print $T + $m - 1}")
    log="$OUT/3d_${L}_$a.log"
    # GPU partagé : attendre qu'il soit libre (utilisation sous 30 %).
    while [ "$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits | head -1)" -gt 30 ]; do sleep 60; done
    AILE_MOY=$m AILE_NREL=$((4 * m)) AILE_DECOL=1 "$X/aile_infinie$X2" 16 $R 6e6 "$a" "$t" - 60 > "$log" 2>&1
    s=$(grep 'relevés sur' "$log")
    czb=$(echo "$s" | sed -n 's/.*bilan (marge [0-9.]*) Cz \([-0-9.]*\) ± \([0-9.]*\), Cx \([-0-9.]*\).*/\1/p')
    sdb=$(echo "$s" | sed -n 's/.*bilan (marge [0-9.]*) Cz \([-0-9.]*\) ± \([0-9.]*\), Cx.*/\2/p')
    cxb=$(echo "$s" | sed -n 's/.*bilan (marge [0-9.]*) Cz [-0-9.]* ± [0-9.]*, Cx \([-0-9.]*\).*/\1/p')
    czs=$(echo "$s" | sed -n 's/.*surface Cz \([-0-9.]*\) ±.*/\1/p')
    cm=$(echo "$s" | sed -n 's/.*Cm \([-0-9.]*\) ± \([0-9.]*\) .*/\1/p')
    sdm=$(echo "$s" | sed -n 's/.*Cm \([-0-9.]*\) ± \([0-9.]*\) .*/\2/p')
    d=$(grep 'décollement x/c (vitesse' "$log")
    dh=$(grep 'au-dessus de la paroi' "$log")
    xeh=$(echo "$dh" | sed -n 's/.*extrados \([0-9.]*\) ·.*/\1/p')
    xih=$(echo "$dh" | sed -n 's/.*intrados \([0-9.]*\)$/\1/p')
    xe=$(echo "$d" | sed -n 's/.*extrados \([0-9.]*\) ;.*/\1/p')
    xi=$(echo "$d" | sed -n 's/.*intrados \([0-9.]*\) ;.*/\1/p')
    du=$(sed -n 's/.* · \([0-9]*\) s$/\1/p' "$log")
    echo "$a,$czb,$cxb,$czs,$cm,$sdb,$sdm,$xe,$xi,$xeh,$xih,$du" | tee -a "$OUT/3d_$L.csv"
    sleep "$PAUSE"
  done
done
