# Liste de travail

Établie le 25 novembre 2026, après la propagation en place. État de
chaque point tenu à jour ici ; le détail des travaux faits va dans le
README (section du sujet, et « Coût de calcul des modifications »).

## Précision (écarts restants avec les logiciels commerciaux)

| # | Travail | État |
| --- | --- | --- |
| 1 | Décollement : loi de paroi avec gradient de pression relevée dans la couche interne, raffinement local de paroi (une dizaine de mailles dans la couche limite) là où le gradient est adverse (bord de fuite, lunette) | en cours : 16° Cz 1,74 → 1,69 (couche de loi) ; couche ×32 possible en mémoire mais diverge (couches très fines à forte incidence instables, à instruire) |
| 2 | Frottement 2 à 3 fois trop fort en aval du recollement de la marche (ν_t de Spalart–Allmaras) | en partie : couche de la loi de paroi cohérente (défaut), −20 % sur la marche, aile et plaque améliorées ; reste ≈ 3 fois la mesure ; k-ε à faire |
| 3 | Corps d'Ahmed : 35° en régime faux (pas de décollement franc à l'arête) ; 25° : Cx 1,2 à 1,3 fois la mesure (culot) | à faire |
| 4 | Arêtes vives : juger `SOUFFLERIE_ARETE_CONVEXE` sur Ahmed 25° à L = 136 avant de l'activer par défaut | à faire |
| 5 | Surfels : perte à l'arrêt et au bord d'attaque, voxels du bord de fuite non découpés | à faire |
| 6 | `--raff-auto` : bilan faux | fait : faces du bilan hors des nœuds recopiés (marge repoussée) ; bilan indépendant de la marge à 0,5 % près (`ar4b`) |
| 7 | Bord d'attaque en tronçons (TR-1208) : la mémoire tient, le raccord non | à faire |
| 8 | Branches de sauvegarde `claude/bord-de-fuite-dcd548` (diagnostics te3d, `--decalage-z`) et `claude/valid-decrochage-0f9767` (journaux) : exploiter ou clore | fait : `te3d` et `run3d --decalage-z` repris ; essais de paroi 16384/32768/65536 écartés (bit 16384 réattribué depuis, question du bord de fuite close par son calage) ; journaux bruts du décrochage rangés dans `validation/decrochage/journaux` (résultats déjà dans le README). Les deux branches peuvent être supprimées |

## Vitesse

| # | Travail | Gain attendu |
| --- | --- | --- |
| 9 | Propagation en place, veine boîte : écritures positives chez la voisine mieux groupées | les 11 % perdus |
| 10 | Propagation en place : interpolation des fantômes (2,3 → 3,0 ms) | ≈ 2 % du pas |
| 11 | Collision pariétale (39 % de la bande passante) : listes de liens compactées | ≈ 5 % du pas |
| 12 | Vitesses en demi-précision | ≈ 8 %, à juger sur la précision |

## Validation

| # | Travail |
| --- | --- |
| 13 | Ahmed et aile à plus haute résolution (≈ 66 M mailles possibles), contre les mesures |
| 14 | Cylindre 3D à Re = 3900 et suites aéroacoustiques (voir [EVOLUTIONS.md](EVOLUTIONS.md)) |

## Interface

| # | Travail | État |
| --- | --- | --- |
| 15 | Vérifier à l'écran la barre de progression de la mise en place | fait : vérifiée par capture (préréglage Validation) ; avancement pondéré par la durée estimée de chaque grille, étape aussi au centre de la vue |
| 16 | Rendu directement dans la fenêtre (sans relecture de l'image) | plus long que prévu : la fenêtre tourne sous DirectX 12 (choisi pour la robustesse des pilotes) et le calcul sous Vulkan (DirectX 12 y est 100 fois plus lent) ; deux périphériques ne partagent pas de texture. Il faudrait passer la fenêtre sous Vulkan avec le périphérique du calcul (repli DirectX 12 et relecture sinon), puis rendre dans une texture enregistrée auprès d'egui. Gain : surtout l'image en pleine résolution des grandes fenêtres (plafond actuel 900 pixels de large) |
| 17 | Version Android | abandonnée |
