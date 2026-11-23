#!/bin/sh
# Banc de vitesse : deux exécutables de `banc3d` (avant, après) lancés en
# alternance sur les cas de référence, médiane des durées de calcul, et
# écart relatif maximal des efforts et sommes de contrôle du champ.
#   banc/vitesse.sh avant.exe apres.exe [n] [cas…]
# Vérifier que le GPU est libre (nvidia-smi) avant de lancer.
A=$1; B=$2; N=${3:-3}
shift 3 2>/dev/null
CAS=${*:-"ahmed:1 ahmedr:1 aile:1 ailer:0.5 ar4c3:0.25"}
T=${TMPDIR:-/tmp}/vitesse.$$
mkdir -p "$T"
printf "%-8s %10s %10s %8s %10s %10s %7s  %s\n" cas "ms/pas av." "ms/pas ap." "gain" "MLUPS av." "MLUPS ap." "pas" "écart max (efforts, champ)"
for ct in $CAS; do
  c=${ct%%:*}; t=${ct#*:}
  for i in $(seq "$N"); do
    for e in A B; do
      [ $e = A ] && x=$A || x=$B
      "$x" "$c" --t "$t" --sortie "$T/$c.$e.txt" 2>/dev/null | grep MLUPS >> "$T/$c.$e.log"
    done
  done
  for e in A B; do
    # Durée de calcul (s), mises à jour par seconde, pas.
    awk '{for (k = 1; k <= NF; k++) { if ($k == "calcul") c = $(k + 1); if ($k == "MLUPS)") m = substr($(k - 1), 2); if ($k == "pas)") p = substr($(k - 1), 2) } print c, m, p}' "$T/$c.$e.log" | sort -n | awk -v n="$N" 'NR == int((n + 1) / 2)' > "$T/$c.$e.med"
  done
  read ca ma pa < "$T/$c.A.med"
  read cb mb pb < "$T/$c.B.med"
  d=$(paste "$T/$c.A.txt" "$T/$c.B.txt" | awk '{a = $2; b = $4; r = (a != 0) ? (b - a) / a : b - a; if (r < 0) r = -r; if ($1 ~ /^champ/) { if (r > dc) dc = r } else if (r > de) de = r } END {printf "%.1e, %.1e", de, dc}')
  awk -v c="$c" -v ca="$ca" -v cb="$cb" -v ma="$ma" -v mb="$mb" -v p="$pa" -v d="$d" 'BEGIN {printf "%-8s %10.2f %10.2f %7.2fx %10.0f %10.0f %7d  %s\n", c, 1000 * ca / p, 1000 * cb / p, ca / cb, ma, mb, p, d}'
done
rm -rf "$T"
