#!/bin/sh
# Comparaison A/B de shaders (dossiers SOUFFLERIE_WGSL), lancements alternés
# pour amortir les variations du GPU (partage, température).
#   banc/ab.sh cas passe dossierA dossierB [n]
# Rend la médiane de la durée de la passe (ms) pour chaque dossier.
B=${B:-target/release/examples/banc3d.exe}
CAS=$1; PASSE=$2; A=$3; BB=$4; N=${5:-5}
for i in $(seq "$N"); do
  for d in "$A" "$BB"; do
    t=$(SOUFFLERIE_WGSL=$d $B "$CAS" --t 0.05 --profil 1 --reps 100 2>/dev/null | grep "$PASSE" | head -1 | awk '{print $(NF-3)}')
    echo "$d $t"
  done
done | awk '{v[$1]=v[$1]" "$2} END {for (d in v) {n=split(v[d],a," "); asort(a); printf "%s médiane %s ms (%s)\n", d, a[int((n+1)/2)], v[d]}}'
