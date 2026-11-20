# Évolutions à faire

Travaux repérés mais non faits, avec ce qu'on sait déjà.

## Cylindre 3D à Re = 3900 (Parnaudeau et al. 2008)

But : St ≈ 0,21, Cd ≈ 0,98–1,0, longueur de recirculation ≈ 1,5 D, à
résolution légère, puis le bruit rayonné (FW-H 3D).

Constats (branche `claude/valid-acoustique-192776`, exemple
`lbm/examples/cylindre3d.rs`, veine 3D existante 16 × 4 × 8 D, sous-maille
de Smagorinsky, sans loi de paroi) :

- diverge avant 12 D/U à D = 16 et à D = 24 (τ₀ = 0,5012 / 0,5018) ;
- à Re = 1000, D = 16, le calcul tient mais donne Cd ≈ 2,0 (attendu ≈ 1) ;
- `periodic_span` n'est jamais mis à vrai dans `gpu3d.rs`, alors que la
  maquette se présente comme « envergure périodique 4D » ;
- blocage D/h = 1/8 (Parnaudeau : domaine 20 D × πD × 20 D).

À faire : stabilité de la paroi 3D à τ ≈ 0,5 sans loi de paroi ; envergure
réellement périodique ; veine plus haute ou correction de blocage fondée
(Allen–Vincenti) ; vérifier la surface de référence des efforts ; puis
tableau St / Cd / Cl′ dans le README. Aucun facteur de réglage.

## Aéroacoustique : suites

- FW-H 3D temporelle (formulation 1A de Farassat) pour le cas ci-dessus ;
  le module `lbm/src/fwh.rs` ne traite que le 2D (forme fréquentielle).
- Sondes et couches absorbantes sur quatre côtés en 3D (`gpu3d.rs`) ; elles
  n'existent qu'en 2D (`gpu.rs`).
- Rayonnement de traînée (2f₀) : les surfaces FW-H traversées par le
  sillage le surestiment de 1,8 à 3,2 fois ; essayer la moyenne de plusieurs
  « bouchons » de sortie (Shur et al. 2005) plutôt qu'une surface ouverte.
- Coût des couches absorbantes sur quatre côtés : non mesuré séparément.
