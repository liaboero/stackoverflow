#!/bin/sh
# Validation de la soufflerie 3D : polaire Cx(Re) de la sphère, et aile
# rectangulaire NACA 0012 d'allongement 4 en 3D contre la surface portante.
# Usage : tools/validation3d.sh [dossier de sortie]   (depuis soufflerie-rs)
# Pause entre les calculs : le GPU d'un portable se bride vite en température.
set -u
E=${SOUFFLERIE:-target/release/soufflerie}
[ -x "$E" ] || E=target/release/soufflerie.exe
OUT=${1:-validation3d}
mkdir -p "$OUT"
PAUSE=${PAUSE:-20}

echo "Re,L,Cx,Cz,Cy,St,Cx_reference" > "$OUT/sphere.csv"
for re in 20 50 100 200 300 500 1000; do
  for l in 12 16; do
    r=$($E run3d --shape sphere --re $re --l $l --t 50 --avg 20 2>/dev/null)
    cx=$(echo "$r" | sed -n 's/^Cx = \([-0-9.]*\).*/\1/p')
    cz=$(echo "$r" | sed -n 's/.*Cz = \([-0-9.]*\).*/\1/p' | head -1)
    cy=$(echo "$r" | sed -n 's/.*Cy = \([-0-9.]*\).*/\1/p' | head -1)
    st=$(echo "$r" | sed -n 's/.*St = \([-0-9.]*\).*/\1/p' | head -1)
    ref=$(echo "$r" | sed -n 's/.*Cx ≈ \([-0-9.]*\).*/\1/p')
    echo "$re,$l,$cx,$cz,$cy,$st,$ref" | tee -a "$OUT/sphere.csv"
    sleep "$PAUSE"
  done
done

echo "alpha_deg,Re,L,Cx,Cz" > "$OUT/aile3d.csv"
for a in 0 4 8 12; do
  for re in 200 1000; do
    r=$($E run3d --shape wing --alpha $a --re $re --l 16 --t 30 --avg 10 2>/dev/null)
    cx=$(echo "$r" | sed -n 's/^Cx = \([-0-9.]*\).*/\1/p')
    cz=$(echo "$r" | sed -n 's/.*Cz = \([-0-9.]*\).*/\1/p' | head -1)
    echo "$a,$re,16,$cx,$cz" | tee -a "$OUT/aile3d.csv"
    sleep "$PAUSE"
  done
done

# Surface portante (polaires de section au plus bas Reynolds disponible).
$E wing --naca 0012 --ar 4 --taper 1 --re 5e4 --alpha 0:12:4 > "$OUT/aile_portante.csv" 2>/dev/null
echo "Terminé : $OUT"
