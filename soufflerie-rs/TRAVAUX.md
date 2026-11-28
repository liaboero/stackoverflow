# Liste de travail

Établie le 25 novembre 2026, après la propagation en place. État de
chaque point tenu à jour ici ; le détail des travaux faits va dans le
README (section du sujet, et « Coût de calcul des modifications »).

## Précision (écarts restants avec les logiciels commerciaux)

| # | Travail | État |
| --- | --- | --- |
| 1 | Décollement : loi de paroi avec gradient de pression relevée dans la couche interne, raffinement local de paroi (une dizaine de mailles dans la couche limite) là où le gradient est adverse (bord de fuite, lunette) | en grande partie : les divergences à forte incidence n'étaient pas le Mach du réseau mais un nœud pariétal du pic d'aspiration emballé par sa paroi, le terme de paroi mobile de Bouzidi étant pris à ρ_w = 1 ; corrigé par défaut (ρ du nœud, `SOUFFLERIE_RHO_PAROI=0` pour l'ancien) : aile 64 mailles, Cz max 1,61 à 17° (Ladson 1,61–1,63), puis 1,52 / 1,43 à 18 / 19°, 20° ne diverge plus ; couches fines stables (×32 à 17° : Cz 1,624, la portance ne croît plus quand on affine la paroi). Prix : traînée en attaché +10 à 16 % (4° : Cx 0,0181 contre 0,0156 ; Cz 0,452), plaque échange / loi 0,00454 / 0,00401 (0,00417 / 0,00404), mode pair-impair de la marche au deuxième rang 0,12 (0,026) ; Ahmed inchangé. Reste : aucun reflux à la paroi (le glissement de Spalding garde le signe de u_τ), chute après le maximum plus douce que la mesure (comme SA), creux de 15° à 64 mailles ; mesure de l'aile à 128 mailles et raffinée à faire |
| 2 | Frottement 2 à 3 fois trop fort en aval du recollement de la marche (ν_t de Spalart–Allmaras) | en cours : source trouvée (mode pair-impair des premiers rangs, paroi sous loi à q ≈ ½, glissement lissé = paroi adhérente pour les fluctuations) et corrigée par défaut : glissement moyen relaxé + fluctuation donnée par la loi (u_n′ − y_n g′) ; marche q = ½ : Re_θ 4 200 (mesure 5 000), x_r 6,30 h au premier nœud (6,26), Cf de la loi 0,0009 / 0,0018 à 10 / 20 h (≈ 0,001 / 0,002) ; plaque et aile inchangées. Écart échange / loi en aval : bilan du premier rang (inertie de la maille pariétale qui accélère, moins la pression), cohérent, la loi d'équilibre l'ignore ; plancher indispensable (plaque, aile), TBLE gp/gpc moins bonne (mode au rang 2). Reste : relèvement du frottement après la bulle trop lent (x_r frottement 10,7 h contre 6,26) ; glissement asservi à l'échange essayé et écarté (juste sur la marche, x_r 7 h ; l'échange d'une maille n'est pas une mesure locale du frottement sur paroi courbe : aile Cz −20 %) |
| 3 | Corps d'Ahmed : 35° en régime faux (pas de décollement franc à l'arête) ; 25° : Cx 1,2 à 1,3 fois la mesure (culot) | en partie : 35° en bon régime (Cx 0,571 → 0,402, mesure 0,26). Soubassement freiné trouvé et corrigé : le sol mobile (rebond à q = ½, τ ≈ ½) portait le mode pair-impair, premier rang en retour à côté du tapis, ν̃ démesuré, soubassement à mi-pression totale ; sol glissant par défaut (limite sans frottement du tapis) : jet de sortie 42,9 m/s (42,8 mesurés), Cx 0,391 → 0,335 à 25°, 0,402 → 0,347 à 35° (L = 136), Cx boîte 0,323 (mesure 0,285). Reste : avant C_K −0,013 contre 0,016 (tapis idéal contre plancher fixe de l'essai : plancher fixe sous loi de paroi à faire) ; culot C_B 0,14 contre 0,09, source dans le sillage proche (tourbillon bas absent, descente de l'écoulement de lunette deux à trois fois trop faible), pas dans le soubassement |
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
