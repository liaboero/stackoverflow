# Soufflerie de poche

Petite soufflerie numérique 2D pour mobile, basée sur la méthode de
Boltzmann sur réseau (Lattice Boltzmann, D2Q9). Elle calcule en direct les
efforts aérodynamiques sur des formes simples : **Cx** (traînée), **Cz**
(portance), **Cm** (tangage) et le nombre de **Strouhal**.

C'est une application web progressive (PWA) : HTML/CSS/JavaScript sans
dépendance ni étape de compilation. Elle s'installe sur l'écran d'accueil
d'un téléphone (iOS ou Android) et fonctionne ensuite hors ligne.

Deux moteurs de calcul partagent le même modèle : un solveur **GPU**
(WebGL2, choisi automatiquement quand le navigateur le permet) et un solveur
**CPU** en JavaScript, qui sert de repli et de référence pour les tests.

## Fonctions

- **Maquettes** : cylindre (avec rotation, effet Magnus), carré, ellipse
  (épaisseur réglable), plaque plane, profil NACA à 4 chiffres, triangle.
  Incidence α réglable (positive = cabrée).
- **Écoulement** : Reynolds de 10 à 20 000, résolutions de 12 à 32 mailles
  par longueur de référence sur CPU, jusqu'à 48 et 64 mailles sur GPU.
- **Visualisation** : vorticité, vitesse ou coefficient de pression Cp,
  traceurs (lignes d'émission), force résultante moyenne et ses composantes,
  trois cadrages (maquette, sillage, veine entière). Palettes adaptées au
  thème clair ou sombre.
- **Mesures** : moyennes glissantes calculées sur un nombre entier de
  périodes quand l'écoulement est périodique, amplitude des oscillations,
  Strouhal, état de convergence (transitoire, stationnaire, périodique).
- **Échelle physique et vitesse** : longueur de référence en mètres et
  fluide (air, eau) → vitesse amont V = Re·ν/L, durée réelle simulée et
  efforts en N/m d'envergure ; vitesse de simulation imposée (temps réel,
  ×0,1, ×0,5, ×2, ×10) ou maximale.
- **Polaires** : enregistrement de points, balayage automatique en incidence
  (ou en vitesse de rotation pour le cylindre), tracé Cx, Cz = f(α) ou
  polaire d'Eiffel Cz(Cx), copie du tableau vers un tableur.

## Utilisation sur téléphone

1. Héberger le dossier `soufflerie-lbm/` en HTTPS (voir « Déploiement »).
2. Ouvrir l'adresse sur le téléphone, puis :
   - iOS (Safari) : Partager → « Sur l'écran d'accueil » ;
   - Android (Chrome) : menu ⋮ → « Installer l'application ».

Variante sans serveur : `npm run build` produit
`dist/soufflerie-de-poche.html`, un fichier unique (CSS et JS intégrés) à
ouvrir directement dans le navigateur du téléphone.

En local :

```bash
cd soufflerie-lbm
npm run serve        # http://localhost:8000
```

## Méthode numérique

| Élément | Choix |
| --- | --- |
| Réseau | D2Q9, Δx = Δt = 1, c_s² = 1/3 |
| Collision | Régularisée récursive à l'ordre 3 (Malaspinas 2015) par défaut, BGK en option ; viscosité de sous-maille de Smagorinsky (C_s = 0,15) quand τ < 0,51, amortie près de la paroi (type Van Driest) pour ne pas épaissir la couche limite |
| Paroi | Rebond à mi-chemin interpolé de Bouzidi (forme exacte, pas l'escalier), y compris pour les parois plus minces qu'une maille (bord de fuite, plaque) ; correction de paroi mobile pour la rotation |
| Efforts | Échange de quantité de mouvement (Ladd), forme de Wen et al. (2014) pour les parois mobiles |
| Entrée | Vitesse imposée, partie hors équilibre extrapolée (Guo et al. 2002) |
| Sortie | Pression imposée, précédée d'une zone éponge à viscosité croissante |
| Haut et bas | Glissement + bandes absorbantes (relaxation vers l'écoulement amont) qui évacuent les ondes acoustiques |
| Démarrage | Rampe de vitesse sur 5 t* avec accélération uniforme du fluide (pas de choc de démarrage) |
| Domaine | 24 L × 12 L, maquette à 8 L de l'entrée ; vitesse réseau U = 0,1 (Ma ≈ 0,17) |

Coefficients (par unité d'envergure, ρ∞ = 1) :
Cx = Fx / (½ρU²L), Cz = Fz / (½ρU²L), Cm = M / (½ρU²L²), positif
cabreur, réduit au point indiqué (centre, ou quart de corde pour les
profils). Strouhal St = f·L/U, mesuré sur Cz.

## Validation

`npm test` (Node ≥ 18, ~2 min, moteur CPU, collision régularisée) vérifie :

| Cas | Calcul | Référence (milieu infini) |
| --- | --- | --- |
| Écoulement uniforme sans obstacle | écart < 1·10⁻⁶ | 0 |
| Cylindre, Re = 20 | Cx = 2,37 | 2,05 – 2,09 |
| Cylindre, Re = 100 | Cx = 1,50 ; St = 0,176 ; Cz' = 0,375 | 1,33 – 1,35 ; 0,164 – 0,166 ; 0,33 |
| NACA 0012, α = 0° | Cz = 0,000 | 0 |
| NACA 0012, α = ±4°, Re = 500 | Cz = ±0,207 (antisymétrie exacte) | ≈ 0,22 (pente ≈ π/rad) |
| NACA 0012, α = 0°, Re = 1000 | Cx = 0,145 / 0,123 / 0,115 (corde de 16 / 24 / 32 mailles) | 0,119 (Di Ilio et al. 2018 ; XFOIL) |
| NACA 0012, α = 0°, Re = 2000, corde de 24 mailles | Cx = 0,096 | 0,084 |
| Cylindre tournant ΩR/U = 1, Re = 100 | Cz = 2,64 | ≈ 2,5 (Mittal & Kumar 2003) |
| Cylindre, Re = 5000 | stable | — |

Les écarts sur le cylindre viennent surtout du **confinement** du domaine
(8 % de blocage, entrée à 8 diamètres) : dans un domaine de 40 L × 24 L, le
même solveur donne Cx = 1,43 et St = 0,167 à Re = 100. La convergence en
maillage est bonne grâce à la paroi interpolée (Re = 20 : Cx = 2,750 / 2,728 /
2,715 pour 16 / 24 / 32 mailles par diamètre dans l'ancien domaine 16 L × 8 L).

## Collision régularisée

Le schéma BGK relaxe toutes les populations vers l'équilibre au même taux.
Près de τ = 0,5 (Reynolds élevé), ses modes non hydrodynamiques sont à peine
amortis et le calcul devient bruité puis instable. La collision régularisée
récursive ne conserve que ρ, u et le tenseur hors équilibre
a⁽¹⁾ = Π − ρc_s²I − ρuu. Elle reconstruit les populations par projection
d'Hermite jusqu'à l'ordre 3, en calculant a⁽¹⁾ d'ordre 3 par récurrence. Les
modes parasites sont filtrés à chaque pas. Elle ne coûte pas plus cher que
BGK : l'équilibre n'a pas à être évalué population par population.

| Cas | BGK | Régularisée |
| --- | --- | --- |
| Cylindre, Re = 100, 16 mailles | Cx = 1,490 ; St = 0,176 | Cx = 1,496 ; St = 0,176 |
| NACA 0012, Re = 1000, 24 mailles | Cx = 0,125 | Cx = 0,124 |
| NACA 0012, α = 4°, Re = 500, 12 mailles | Cz = 0,183 | Cz = 0,207 (réf. ≈ 0,22) |
| Cylindre, Re = 5000, 12 mailles, sans LES | — | stable |
| NACA 0012, Re = 10 000, 32 mailles | Cz′ ≈ 0,13 parasite | sans oscillation parasite |

## Moteur GPU

Le solveur WebGL2 (`js/lbm-gpu.js`) reprend exactement le modèle du CPU :
même collision, mêmes parois de Bouzidi (géométrie construite sur le CPU par
`js/geometry.js`, puis envoyée en textures), mêmes conditions aux limites.
Les 9 populations sont stockées en ping-pong dans deux textures RGBA32F et
une R32F (36 octets par maille) ; un pas de temps est un seul tracé avec
trois cibles de rendu. Les efforts sont calculés à la demande par une passe
sur les mailles de paroi, réduite sur le GPU puis relue de façon asynchrone
(tampon de pixels et barrière), sans jamais attendre le GPU. Le champ
affiché et les traceurs sont calculés et dessinés sur le GPU
(`js/render-gpu.js`). La maquette et la force sont dessinées sur un calque
2D superposé.

Optimisations du moteur :

- 36 au lieu de 48 octets par maille : ρ et u ne sont plus écrits à chaque
  pas mais recalculés une fois par image pour l'affichage ;
- textures liées une fois pour toutes à des unités réservées, une copie du
  programme par sens du ping-pong et uniformes mis en cache : 3 à 4 appels
  WebGL par pas au lieu d'environ 36, ce qui compte sur téléphone ;
- relecture asynchrone des efforts et réglage du nombre de pas par image
  sur le retard du GPU (au plus une image de retard) ;
- 16 échantillons d'efforts par t* au lieu de 40.

`tests/gpu-bench.html` mesure la vitesse du moteur sur l'appareil qui
l'ouvre. Avec le moteur logiciel SwiftShader, seul disponible dans mon
environnement de test : +31 % à 16 mailles, +18 % à 32 mailles. Je n'ai
pas pu mesurer de GPU matériel.

Seule différence : les colonnes d'entrée et de sortie sont extrapolées
depuis le pas précédent. `tests/gpu-compare.html` compare les deux moteurs
(à ouvrir via un serveur local). Mesuré avec le moteur logiciel SwiftShader,
l'écart sur Cx et Cz reste sous 1 % : cylindre à Re = 100 et 5000, cylindre
tournant, NACA 0012 à 8°, plaque à 30° en BGK.

Le moteur exige WebGL2 et l'extension `EXT_color_buffer_float` (navigateurs
récents, iOS ≥ 15, la plupart des Android). Sinon, ou si le GPU refuse la
grille, l'application repasse sur le CPU. Le choix se fait dans « Réglages
numériques ». Mémoire GPU : ≈ 115 Mo pour la résolution Max (1536 × 768).

## Profils : portance, oscillations et bord de fuite

**Bord de fuite recalé.** La pointe d'un bord de fuite aigu est plus fine
qu'une maille ; c'est pourtant elle qui fixe la circulation (condition de
Kutta). Sur un réseau cartésien, sa position effective change d'une maille
à l'autre quand l'incidence varie, ce qui revient à une erreur d'incidence
de ±1 à 2° à 16 mailles. Sans correction, Cz(α) variait en dents de scie
(NACA 0012, Re = 1000, 16 mailles : Cz = 0,003 / 0,050 / 0,170 à 1 / 2 /
3°). Le profil est donc décalé de moins d'une demi-maille pour que sa pointe
garde toujours la même position par rapport au réseau : à l'abscisse d'un
nœud, à mi-hauteur entre deux nœuds, comme à α = 0. Parmi les positions
essayées, c'est celle qui donne la portance la plus régulière. Les efforts ne dépendent pas de
cette translation. Résultat : Cz = 0,047 / 0,083 / 0,101 à 16 mailles, et
une pente constante de 0,046 ± 0,003 par degré à 32 mailles
(référence ≈ 0,055).

**Efforts moyennés sur plusieurs pas.** Avec la collision régularisée,
certaines configurations de bord de fuite excitent un mode parasite qui
alterne d'un pas de temps au suivant. Il reste local et ne perturbe pas
l'écoulement moyen, mais il s'ajoute aux efforts instantanés (Cx alternant
entre 0,24 et 0,06, par exemple). Un échantillonnage régulier le replie en
un biais : jusqu'à Cx = −0,28 au lieu de 0,135 sur un cas. Les efforts sont
désormais moyennés sur tous les pas (CPU) ou sur deux pas consécutifs (GPU),
ce qui annule exactement cette alternance.

**Oscillations physiques.** Au-delà d'une incidence critique, le sillage
d'un profil à bas Reynolds devient instationnaire. Des tourbillons se
détachent du bord de fuite, avec un Strouhal (rapporté à la corde) de
l'ordre de 1. Calculé pour le NACA 0012 :

| Re | stationnaire jusqu'à | périodique dès |
| --- | --- | --- |
| 1000 (24 mailles) | 6° | ≈ 8 à 10° (Cz′ = 0,002 à 8° ; littérature : 7 à 8°) |
| 2000 (24 mailles) | 6° (Cz′ = 0,005) | ≈ 8° (Cz′ = 0,035, St ≈ 0,9) |
| 5000 (32 mailles) | 4° | ≈ 6° (Cz′ = 0,016, St ≈ 1,3) |

L'angle critique baisse quand Re augmente, comme l'observe Kurtulus (2022)
entre Re = 1000 et 4000. Avec le moteur GPU, qui rend accessibles les
Reynolds de quelques milliers, les profils sont donc souvent oscillants dès
5 à 6° : c'est la physique 2D de ces régimes. Les oscillations d'amplitude
négligeable (Cz′ < 0,01) sont affichées « quasi stationnaire ».

**Paroi en marches.** Sur un profil mince (2 mailles d'épaisseur à
16 mailles de corde), le rebond simple perd la symétrie, même à α = 0. Ce
mode ne sert qu'à comparer, sur des corps épais.

## Traînée et décrochage

**Traînée.** Aux Reynolds accessibles (10 à 5000), la traînée est dominée
par le frottement d'une couche limite laminaire épaisse : Cx ≈ 0,12 pour un
NACA 0012 à Re = 1000, contre ≈ 0,006 à Re = 10⁶. Ce n'est pas un défaut du
calcul : à corde de 24 ou 32 mailles, l'écart à la référence fine est de
quelques pourcents. En résolution rapide (12 mailles), la traînée d'un
profil est surestimée d'environ 20 à 30 %.

**Décrochage.** Rien n'est modélisé : le décollement, la recirculation et le
lâcher tourbillonnaire sortent directement de la résolution des équations.
À Re = 1000, le calcul retrouve le décollement de bord de fuite et le début
du lâcher tourbillonnaire vers α ≈ 8° (7 à 8° dans la littérature), puis un
plafonnement de Cz avec forte hausse de traînée. Il s'agit du décrochage
progressif des bas Reynolds, pas du décrochage d'une aile réelle
(Re ~ 10⁶ – 10⁷), qui dépend de la transition et de la turbulence de couche
limite, absentes d'un calcul 2D aussi grossier. La portance des profils est
sensible à la résolution du bord de fuite (± 15 %) et reste en dessous de la
référence fine (pente ≈ π par radian à Re = 1000).

## Jusqu'à quel Reynolds ?

Deux limites, la plus basse l'emporte.

**Résolution.** Erreur sur la traînée du NACA 0012 à α = 0 (références :
Di Ilio et al. 2018 et XFOIL, Cx = 0,119 / 0,084 / 0,053 / 0,037 – 0,040
à Re = 1000 / 2000 / 5000 / 10 000). Jusqu'à 32 mailles : collision BGK ;
48 mailles : collision régularisée.

| L (mailles) | Re = 1000 | Re = 2000 | Re = 5000 | Re = 10 000 |
| --- | --- | --- | --- | --- |
| 12 (Rapide) | +54 % | +96 % | — | — |
| 16 (Standard) | +22 % | +42 % | +99 % | — |
| 24 (Fine) | +5 % | +14 % | +35 % | — |
| 32 (Très fine) | −4 % | −4 % | +1 % | +19 à +29 % |
| 48 (Ultra, GPU) | — | — | +2 % | +5 à +13 % |
| 64 (Max, GPU) | — | — | — | −6 à −13 % |

Règle retenue : traînée à ± 10 % environ jusqu'à Re ≈ 350 (Rapide), 650
(Standard), 1450 (Fine), 5000 (Très fine) et 10 000 (Ultra et Max). À
Re = 10 000, les deux résolutions GPU encadrent les références, qui
diffèrent elles-mêmes de 8 %. Re = 15 000 et 20 000 restent accessibles mais
sans référence pour les valider ; l'application les signale comme hors du
domaine fiable.

**Coût.** Le coût par temps convectif croît comme L³ (L² mailles,
L pas de temps par t*) : chaque palier de résolution coûte 2 à 4 fois plus.
Moteur CPU, mesuré dans Chromium sur un processeur de serveur (un seul
cœur) : 2,7 / 1,1 / 0,28 / 0,15 t*/s pour Rapide / Standard / Fine / Très
fine à Re = 1000. Les résolutions Ultra (664 000 mailles) et Max
(1,2 million) ne sont proposées qu'avec le moteur GPU. Je n'ai pas pu
mesurer de GPU matériel ; la ligne de performance de l'application donne la
vitesse réelle.

Avec BGK, Re = 10 000 restait inexploitable : la traînée oscillait de façon
parasite à 32 mailles, et le sillage du bord de fuite devenait
asymétrique à 48 mailles. La collision régularisée donne un écoulement
propre et symétrique dans les deux cas.

Vers Re ≈ 10⁴, un profil réel présente des bulbes de décollement laminaire
et une transition, phénomènes 3D absents d'un calcul 2D ; l'application le
signale.

**Physique 2D.** Le sillage d'un cylindre devient tridimensionnel vers
Re ≈ 190 (mode A, Williamson 1996 ; Barkley & Henderson 1996). Au-delà, un
calcul 2D surestime Cx, St et les fluctuations de portance (cylindre à
Re = 1000 : Cx ≈ 1,4 en 2D contre ≈ 1 mesuré). Pour les corps non profilés
(cylindre, carré, triangle, plaque ou ellipse très inclinées), rester sous
Re ≈ 200 pour des valeurs quantitatives. Un profil peu incliné reste
laminaire et 2D bien plus loin : c'est alors la résolution qui limite.

## Vers Re = 10⁷ : feuille de route

**Pourquoi c'est un changement de nature.** À Re = 10⁷, la couche limite
d'un profil est turbulente et mince : δ ≈ 0,37 c·Re^(−1/5) ≈ 0,015 c au bord
de fuite, et la sous-couche visqueuse (y⁺ = 1) est à y ≈ 2·10⁻⁶ c. La
résoudre directement demanderait de l'ordre d'un million de mailles par
corde. En outre, la turbulence est tridimensionnelle et un calcul 2D ne peut
pas la résoudre. Comme les codes industriels pour les polaires de profils,
il faut modéliser la turbulence en moyenne (RANS) et ne pas résoudre la
sous-couche (loi de paroi). Côté LBM, τ − ½ = 3ν ≈ 3·10⁻⁵ en unités réseau :
seule la viscosité turbulente du modèle rend alors le calcul stable et
physique.

| Étape | Contenu | Reynolds visé |
| --- | --- | --- |
| A | Modèle de Spalart–Allmaras résolu par différences finies sur le réseau (une passe GPU de plus par pas), ν_t injecté dans τ comme la LES actuelle ; loi de paroi de Spalding (u_τ tiré de la vitesse au premier nœud, imposé par un rebond modifié) | 10⁵ – 10⁶, qualitatif sur grille uniforme |
| B | Raffinement local emboîté, facteur 2 par niveau, 5 à 6 niveaux | 10⁶ – 10⁷, quantitatif |
| C | Transition (déclenchement imposé à x/c donné, puis modèle de transition), passage à WebGPU (calcul, réductions, grilles multiniveaux plus simples) | polaires de profils lisses, ×1,5 à 3 en vitesse |

Ordres de grandeur à Re = 10⁷ (u_τ/U ≈ 0,04) : un premier nœud à y⁺ ≈ 500,
encore dans la zone logarithmique (δ⁺ ≈ 6000), demande ≈ 800 mailles par
corde près de la paroi, soit une douzaine de mailles dans l'épaisseur de
couche limite. En grille uniforme sur 24 c × 12 c, cela ferait 180 millions
de mailles (13 Go). En grille emboîtée, avec une boîte fine de 1,2 c × 0,4 c
et 5 niveaux plus grossiers, on tombe à ≈ 2 millions de mailles (≈ 130 Mo).
Un point de polaire, soit 30 t* et 240 000 pas du niveau fin, coûterait
≈ 1,5·10¹¹ mises à jour : quelques minutes sur un GPU de portable à
1 GLUPS, plutôt une demi-heure sur téléphone.

Références : modèle SA dans une LBM multidomaine (Pellerin, Leclaire et
Reggio, *Computers & Mathematics with Applications*, 2015) ; loi de paroi
analytique couplée au SA sur grilles cartésiennes adaptatives (Ponsin et
Lozano, 2026). Cas de validation visé : NACA 0012 à Re = 6·10⁶, M = 0,15, de
la *Turbulence Modeling Resource* de la NASA (SA : Cd ≈ 0,0081 à 0° ; codes
publiés : Cz ≈ 1,09 à 10° et ≈ 1,54 à 15°).

## Limites

- Calcul **2D** et laminaire : pas de transition ni de turbulence 3D.
- Écoulement faiblement compressible (Ma ≈ 0,17 par défaut, réglable).
- Moteur CPU : un seul cœur, ≈ 15 à 30 millions de mailles mises à jour par
  seconde. Moteur GPU : typiquement 10 à 100 fois plus, selon la carte
  graphique ; la ligne de performance de l'application affiche la valeur
  réelle.

## Structure

```
soufflerie-lbm/
├── index.html            interface
├── css/app.css           styles (thèmes clair et sombre)
├── js/geometry.js        géométrie sur le réseau (commune aux deux moteurs)
├── js/lbm.js             solveur CPU D2Q9 (sans DOM, testable sous Node)
├── js/lbm-gpu.js         solveur GPU WebGL2 (mêmes équations)
├── js/shapes.js          géométries et placement en incidence
├── js/sim.js             veine d'essai : domaine, démarrage, coefficients
├── js/stats.js           moyennes, périodes, Strouhal
├── js/render.js          rendu CPU du champ et des traceurs ; maquette, force
├── js/render-gpu.js      rendu GPU du champ et des traceurs
├── js/chart.js           graphiques (historique, polaire)
├── js/app.js             logique de l'interface
├── sw.js, manifest.webmanifest, icons/   PWA
├── tests/validate.mjs    validation physique (moteur CPU)
├── tests/gpu-compare.html  comparaison GPU / CPU dans le navigateur
├── tests/gpu-bench.html    vitesse du moteur GPU sur l'appareil
└── tools/build.mjs       version « fichier unique »
```

## Déploiement

Le workflow `.github/workflows/soufflerie-pages.yml` lance les tests à chaque
pull request touchant ce dossier et publie l'application sur GitHub Pages à
chaque push sur `master`. Il faut au préalable activer Pages dans les
réglages du dépôt : *Settings → Pages → Build and deployment → Source :
GitHub Actions*. L'application est alors servie à
`https://<utilisateur>.github.io/<dépôt>/`.

Tout autre hébergement statique en HTTPS convient (Netlify, Cloudflare
Pages…) ; le service worker exige HTTPS, sauf sur `localhost`.

## Pistes

- Raffinement local de maillage près de la maquette, qui résoudrait mieux
  le bord de fuite sans payer une grille fine partout.
- Correction de blocage, mode dimensionnel (corde en m, vitesse en m/s,
  air ou eau → efforts en N/m).
- Dessin libre de la maquette au doigt.
