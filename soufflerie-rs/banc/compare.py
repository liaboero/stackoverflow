"""Compare les résultats du banc (banc3d --sortie) à une ou plusieurs références.

  python banc/compare.py nouveau.txt ref1.txt [ref2.txt ...]

Pour chaque grandeur : valeur, moyenne et dispersion (max − min) des
références, et écart à la moyenne rapporté à cette dispersion (avec une
seule référence, ou des références identiques : écart relatif seul).
"""
import sys


def lire(p):
    d = {}
    for ligne in open(p, encoding="utf-8"):
        m = ligne.split()
        if len(m) == 2:
            d[m[0]] = float(m[1])
    return d


nouveau = lire(sys.argv[1])
refs = [lire(p) for p in sys.argv[2:]]
pire = 0.0
for k, v in nouveau.items():
    rv = [r[k] for r in refs if k in r]
    if not rv:
        continue
    moy = sum(rv) / len(rv)
    disp = max(rv) - min(rv)
    rel = (v - moy) / abs(moy) if moy != 0 else v - moy
    if disp > 0:
        x = (v - moy) / disp
        pire = max(pire, abs(x))
        print(f"{k:18} {v:15.8f}  réf. {moy:15.8f} ± {disp / 2:.2e}  écart {rel:+.2e}  ({x:+.2f} dispersion)")
    else:
        print(f"{k:18} {v:15.8f}  réf. {moy:15.8f}            écart {rel:+.2e}{'' if v == moy else '  DIFFÉRENT'}")
if any(max(r[k] for r in refs) != min(r[k] for r in refs) for k in nouveau if all(k in r for r in refs)):
    print(f"plus grand écart : {pire:.2f} fois la dispersion des références")
