# Liste de travail

Établie le 25 novembre 2026, après la propagation en place. État de
chaque point tenu à jour ici ; le détail des travaux faits va dans le
README (section du sujet, et « Coût de calcul des modifications »).

## Précision (écarts restants avec les logiciels commerciaux)

| # | Travail | État |
| --- | --- | --- |
| 1 | Décollement : loi de paroi avec gradient de pression relevée dans la couche interne, raffinement local de paroi (une dizaine de mailles dans la couche limite) là où le gradient est adverse (bord de fuite, lunette) | en cours : divergence des couches fines au pic d'aspiration (Mach du réseau 0,57 à U = 0,1) ; l'écart de Cz entre U = 0,1 et 0,05 n'est pas la compressibilité mais une erreur d'ordre 3 du schéma en 1/U près de la paroi (README, « Vitesse du réseau ») : U = 0,05 n'est pas plus juste, vitesse gardée à 0,1 ; source identifiée : couche limite sous-résolue (≈ 1 maille à L = 64), l'écart baisse avec le raffinement de paroi (22 % → 5 % au ×2, f32) ; fp16 avec raffinement rend le calcul instationnaire (θ ±15 % selon l'instant ; seul le tout-f32 est stationnaire, précision par niveau insuffisante) : moyenner ou passer en f32 (TRT à Λ fixe : borné à τ⁻ ≤ 2, θ indépendant de U mais trop grand ; Λ = 1/16 juste à U = 0,05, instable à 0,1 ; essai `SOUFFLERIE_LAMBDA`) ; décollement toujours absent |
| 2 | Frottement 2 à 3 fois trop fort en aval du recollement de la marche (ν_t de Spalart–Allmaras) | en cours : couche de la loi de paroi cohérente (défaut, SA et k-ε) ; source trouvée : mode pair-impair (damier en x, période 2 pas) des premiers rangs quand la paroi sous loi est à q ≈ ½ (rebond pur), qui porte −u′w′ ≈ 0,009 q vers la paroi ; la marche avait ses parois à q = ½. Paroi à q = 0,45 : frottement de la loi 0,0011 / 0,0019 à 10 / 20 h (mesure ≈ 0,001 / 0,002), couche amont juste (Re_θ 4 400 contre 5 000) ; restent l'échange 0,0005 sous la loi et le recollement (5,6 / 9,7 h). Essai de correction (paroi virtuelle à q = 1, `SOUFFLERIE_PAROI_Q1=1`) : éteint le mode, mais change les cas sans mode (plaque +12 % de frottement, Ahmed Cx boîte +24 %) ; non retenu, correction à trouver |
| 3 | Corps d'Ahmed : 35° en régime faux (pas de décollement franc à l'arête) ; 25° : Cx 1,2 à 1,3 fois la mesure (culot) | en partie : 35° en bon régime (Cx 0,571 → 0,402, mesure 0,26) ; reste le culot (C_B 0,15 contre 0,09) et l'avant (C_K 0,03 contre 0,016) |
| 4 | Arêtes vives : juger `SOUFFLERIE_ARETE_CONVEXE` sur Ahmed 25° à L = 136 avant de l'activer par défaut | fait : activée par défaut (arêtes de moins de 90° seulement) ; arêtes vives à l'échelle de la loi (pieds exclus) ; Ahmed L = 136 : 25° Cx 0,387, 35° Cx 0,402 (lunette décollée, profils 8 fois mieux suivis) |
| 5 | Surfels : perte à l'arrêt et au bord d'attaque, voxels du bord de fuite non découpés | à faire |
| 6 | `--raff-auto` : bilan faux | fait : faces du bilan hors des nœuds recopiés (marge repoussée) ; bilan indépendant de la marge à 0,5 % près (`ar4b`) |
| 7 | Bord d'attaque en tronçons (TR-1208) : la mémoire tient, le raccord non | fait : le raccord tient (Cz +5 % au lieu de −23 %), après les corrections du raccord qui coupe la paroi |
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
