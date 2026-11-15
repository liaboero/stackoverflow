#!/bin/sh
# Compare la soufflerie 3D à VSPAERO et à la mesure sur l'aile NACA TR-1208.
# Usage (depuis soufflerie-rs) :
#   OPENVSP=/chemin/OpenVSP-3.56.0-win64 tools/comparaison_vspaero.sh [dossier] [L] [raffinement]
# OPENVSP : dossier d'OpenVSP (vspscript, vspaero) ; L : mailles par corde
# moyenne de la veine (8) ; raffinement : grilles fines (2, soit 4 L près de
# l'aile). Écrit <dossier>/comparaison.csv et affiche le tableau.
set -eu
VSP=${OPENVSP:?"OPENVSP : dossier d'OpenVSP (https://openvsp.org/download.php)"}
E=${SOUFFLERIE:-target/release/soufflerie}
[ -x "$E" ] || E=target/release/soufflerie.exe
OUT=${1:-comparaison}
L=${2:-8}
R=${3:-2}
mkdir -p "$OUT"
cp tools/vsp/tr1208.vspscript tools/vsp/NACA63A012.af "$OUT/"
S="$VSP/vspscript"
[ -x "$S" ] || S="$VSP/vspscript.exe"
( cd "$OUT" && "$S" -script tr1208.vspscript > vspaero.log 2>&1 )
grep -E '^ *vlm ' "$OUT/vspaero.log" | awk '{print $2","$3","$4}' > "$OUT/vspaero.csv"
: > "$OUT/soufflerie.csv"
for a in 0 2 4 6 8; do
  r=$("$E" run3d --stl "$OUT/TR1208.stl" --unite in --lref 16.672 --sref 2024.35 --marge 0.4 \
      --l "$L" --raff "$R" --re 4e6 --alpha "$a" --t 20 --avg 8 2>/dev/null | grep '^Cx')
  cx=$(echo "$r" | sed -n 's/^Cx = \([-0-9.]*\).*/\1/p')
  cz=$(echo "$r" | sed -n 's/.*Cz = \([-0-9.]*\).*/\1/p')
  echo "$a,$cz,$cx" >> "$OUT/soufflerie.csv"
done
# Mesure : NACA TR-1208, figure 2, lue sur le graphique (±0,01).
printf '0,0.000\n2,0.136\n4,0.257\n6,0.378\n8,0.499\n' > "$OUT/mesure.csv"
echo "alpha,CL_mesure,CL_vspaero,CL_soufflerie,CD_vspaero,CD_soufflerie" > "$OUT/comparaison.csv"
awk -F, 'FILENAME ~ /mesure/ {m[$1+0]=$2; next}
         FILENAME ~ /vspaero/ {vc[$1+0]=$2; vd[$1+0]=$3; next}
         {printf "%s,%s,%.4f,%s,%.4f,%s\n", $1, m[$1+0], vc[$1+0], $2, vd[$1+0], $3}' \
    "$OUT/mesure.csv" "$OUT/vspaero.csv" "$OUT/soufflerie.csv" >> "$OUT/comparaison.csv"
column -s, -t "$OUT/comparaison.csv" 2>/dev/null || cat "$OUT/comparaison.csv"
