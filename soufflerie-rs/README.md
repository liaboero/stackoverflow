# Soufflerie de poche — version native (Rust)

Fork natif de l'application web [`soufflerie-lbm`](../soufflerie-lbm) : même
soufflerie numérique 2D Lattice Boltzmann (D2Q9), aile finie et soufflerie
3D (D3Q19 sur GPU), en Rust, pour ordinateur (Linux, Windows, macOS).
Interface graphique et ligne de commande.

![Interface : NACA 2412 à 6°, Re = 1000, fronts de traceurs](capture.png)

## Pourquoi une version native ?

Le moteur processeur de l'application web tourne sur un seul cœur, en
JavaScript. Ici, la même boucle est compilée, vectorisée (SSE, AVX2,
AVX-512) et répartie sur tous les cœurs. Mesures sur le même cas (NACA 0012,
α = 4°, Re = 1000, 32 mailles par corde, grille 768 × 384), sur un serveur à
4 cœurs :

| Moteur | Millions de mailles par seconde | Rapport |
| --- | --- | --- |
| Application web, CPU (JavaScript, 1 cœur) | 16,5 | 1 |
| Rust, 1 fil | 78 | ×4,7 |
| Rust, 4 fils | 236 | ×14 |
| Rust, 1 fil, `-C target-cpu=native` (AVX-512) | 121 | ×7,3 |
| Rust, 4 fils, `-C target-cpu=native` | 309 | ×19 |

Sur un portable à 8 cœurs, compter de l'ordre de 400 à 800 millions de
mailles par seconde. La polaire (`soufflerie polar`) calcule en plus les
incidences en parallèle : 5 points (NACA 0012, Re = 1000, 16 mailles) en
11 s.

**Moteur GPU natif** (wgpu : Vulkan, Metal, DirectX 12) : mêmes équations, en
« compute shaders ». Il donne les mêmes efforts que le moteur processeur à
10⁻⁵ près (test `gpu`). Sur GPU, c'est la bande passante mémoire qui limite,
comme dans le navigateur : le gain attendu par rapport au moteur WebGL2 de
l'application web est modeste (×1,5 à 2, sans les limites mémoire du
navigateur). Je n'ai pu le mesurer que sur un GPU logiciel (SwiftShader),
sans signification pour la vitesse.

## Ce qui est porté

- D2Q9, collision régularisée récursive (Malaspinas 2015) ou BGK ;
  sous-maille de Smagorinsky amortie près de la paroi (Van Driest) ; entrée
  et sortie de Guo ; bandes absorbantes ; zone éponge ; mise en vitesse
  progressive ; paroi interpolée de Bouzidi, y compris les parois plus
  minces qu'une maille et la paroi mobile (cylindre tournant) ; efforts par
  échange de quantité de mouvement ; moyennes sur un nombre entier de
  périodes et Strouhal.
- Maquettes : cylindre (rotation), carré, ellipse, plaque, profil NACA à 4
  chiffres, triangle ; bord de fuite recalé sur le réseau.
- **Turbulence de Spalart–Allmaras et loi de paroi de Spalding** (Re de
  5·10⁴ à 10⁷), sur les deux moteurs : viscosité turbulente ajoutée au temps
  de relaxation, glissement de paroi qui transmet le frottement de la loi,
  ν̃ = κu_τy au premier nœud, sous-maille conservée loin de la paroi.
- **Méthode des panneaux et couche limite intégrale** (tourbillons
  linéaires, Thwaites, transition eᴺ, Head, Squire–Young, décrochage
  approché de Kirchhoff, correction de compressibilité) : référence Cz, Cx,
  Cm en écoulement attaché, et **couplage hybride** : en transition libre,
  la paroi du calcul LBM reste laminaire en amont de la transition prédite.
- **Distributions pariétales** Cp et Cf du calcul LBM, moyennées dans le
  temps, comparées aux panneaux (profils) ou à l'écoulement potentiel
  (cylindre).
- Interface : veine zoomable (molette, glisser, double-clic pour recadrer),
  vorticité, vitesse, pression ou viscosité turbulente, traceurs en filets
  ou en **fronts** (lignes de temps, couleurs alternées), force résultante,
  coefficients, référence « panneaux + couche limite », historique et
  distributions pariétales, copie CSV. Reynolds de 10 à 10⁷, résolutions
  jusqu'à 192 mailles par longueur (13 millions de mailles).

![NACA 2412 à 6°, Re = 10⁶, Spalart–Allmaras et transition libre : −Cp à la paroi, points LBM et tirets des panneaux](capture-sa.png)

Validation (`cargo test --release`) : mêmes cas que l'application web, mêmes
valeurs à 10⁻³ près.

| Cas | Rust | Application web |
| --- | --- | --- |
| Cylindre, Re = 20, L = 16 | Cx = 2,373 | 2,373 |
| Cylindre, Re = 100, L = 16 | Cx = 1,496, St = 0,176, Cz' = 0,376 | 1,496, 0,176, 0,375 |
| NACA 0012, α = ±4°, Re = 500, L = 12 | Cz = ±0,1972 | ±0,1972 |
| NACA 0012, α = 0, Re = 1000, L = 16 | Cx = 0,1454 | 0,1454 |
| Cylindre tournant ΩR/U = 1, Re = 100 | Cz = 2,642 | 2,642 |
| SA, NACA 0012, Re = 10⁶, α = 4°, transition libre, L = 16 | Cz = 0,490, Cx = 0,0536 | 0,491, 0,0536 |
| SA, NACA 0012, Re = 10⁷, α = 8°, turbulent, L = 16 | Cz = 0,931 | 0,934 |
| SA, NACA 0012, Re = 10⁶, L = 32 : α = 0 turbulent / libre, α = 4° libre | Cx = 0,0126 / 0,0091, Cz = 0,505 | 0,0125 / 0,0091, 0,505 |
| Panneaux, NACA 0012 non visqueux, α = 4° | Cz = 0,4797 | 0,4797 |
| Panneaux + couche limite, Re = 6·10⁶ turbulent : α = 0, α = 10° | Cx = 0,00816 ; Cz = 1,168, Cx = 0,01285 | identiques |
| Cylindre, Re = 40 : Cp au point d'arrêt | 1,145 | (≈ 1,14 publié) |
| Aile finie A = 6, profil 2π, α = 4° | CL = 0,29376, e = 0,9955 | identiques |
| Aile NACA 2412, A = 8, λ = 0,5, Re = 10⁶, α = 4° / 12° / 16° | CL = 0,5056 / 1,1236 / 1,3659 | identiques |
| 3D, sphère Re = 100, D = 8 / 12 mailles | Cx = 1,259 / 1,175 | 1,26 / 1,18 |
| 3D, sphère Re = 100, D = 12, simple / demi-précision | Cx = 1,178 / 1,175 | écart 0,3 % |

Le moteur GPU donne les mêmes efforts que le processeur avec Spalart–Allmaras
(écart ≤ 7·10⁻⁵ sur un profil à Re = 10⁶ et un cylindre tournant à 10⁵).
Spalart–Allmaras coûte environ 45 % de la vitesse sur processeur (131 au lieu
de 239 millions de mailles par seconde sur 4 cœurs), un peu moins sur GPU.

Comme dans l'application web, la portance et le décollement sont
qualitatifs à ces résolutions, et la traînée LBM d'un profil mince reste
approximative : la référence « panneaux + couche limite » est la bonne
valeur de traînée en écoulement attaché. À Re = 10⁶, la traînée LBM d'un
NACA 0012 à 4° vaut 0,006 à L = 32 et 0,020 à L = 48 (0,008 attendu) : elle
dépend de la finesse du bord d'attaque, pas d'un défaut de convergence.

**Traînée retenue.** Elle ne converge pas non plus en affinant la grille
(NACA 2412, Re = 10⁶, α = 4°, panneaux : 0,0078) :

| L (mailles par corde) | 32 | 48 | 64 | 96 |
| --- | --- | --- | --- | --- |
| Cx LBM | 0,022 | 0,018 | 0,021 | 0,016 |

Le bilan de quantité de mouvement sur un volume de contrôle donne les mêmes
valeurs que l'échange aux parois : c'est l'écoulement calculé lui-même qui
porte cette traînée. Le bord d'attaque tient dans une maille, la loi de paroi
impose le frottement, et quelques dixièmes de degré d'erreur sur la direction
de la résultante basculent de la portance (Cz ≈ 1) dans une traînée de
l'ordre de 0,01. Pour un profil attaché sous loi de paroi, l'interface
affiche donc le Cx LBM comme non fiable et retient la traînée « panneaux +
couche limite » (lecture, Cz/Cx, polaire Cz(Cx), colonne `Cx_retenu` du
CSV). Dès que les panneaux prédisent un décollement, la traînée de pression
domine et le Cx LBM est de nouveau retenu.

**Ce qui porte l'excès.** En séparant frottement (∫Cf dx) et pression, la
loi de paroi donne le bon frottement (0,0050 contre 0,0053 pour la couche
limite intégrale) : l'excès est une traînée de pression, qui décroît quand
le bord d'attaque, le bord de fuite et la couche limite sont mieux résolus
(0,017 à L = 64, 0,007 à 128, 0,005 à 192 ; panneaux 0,0025).

**Raffinement local (GPU).** Plutôt que d'affiner toute la veine (coût en
L³), `--raff 1|2` (interface : « Près de la maquette ») ajoute autour de la
maquette une ou deux grilles imbriquées, chacune deux fois plus fine en
espace et en temps que la précédente (emprise : 0,9 corde en amont, 1,9 en
aval, ±0,75 ; puis 0,65, 1,5, ±0,58). Le couplage passe par les moments
(ρ, u, partie hors équilibre de Π remise à l'échelle de chaque grille) :
interpolation cubique le long du bord de la grille fine et linéaire en
temps, recopie filtrée sur 3 × 3 nœuds vers la grille parente (Lagrava
et al., 2012). ν̃ suit la même voie. Efforts et distributions pariétales
viennent de la grille la plus fine, l'affichage de la veine. NACA 2412,
Re = 10⁶, α = 4°, Cx moyen une fois le tourbillon de démarrage sorti :

| Près du corps | Grille | Cx | Calcul (34 t*) |
| --- | --- | --- | --- |
| 128 | L = 128 uniforme | 0,0121 | ≈ 3 min |
| 128 | L = 64 + ×2 | 0,0115 | 23 s |
| 128 | L = 32 + ×4 | 0,0115 | 8 s |
| 192 | L = 48 + ×4 | 0,0088 | 21 s |
| 256 | L = 128 + ×2 | 0,0108 | 241 s |
| 256 | L = 64 + ×4 | 0,0106 | 50 s |

À résolution égale près du corps, les grilles imbriquées retrouvent la
grille uniforme à quelques pour cent, 5 à 20 fois plus vite. La traînée
reste toutefois au-dessus des panneaux (0,0078) et ne décroît plus
au-delà de 192 : elle dépend de la viscosité de sous-maille près de la
paroi, nécessaire à la stabilité à τ ≈ 0,5 (une longueur d'amortissement
pariétal de 4 mailles au lieu de 2 la ramène à 0,008, mais le calcul
devient instable). Trois niveaux (×8) divergent de même dans la couche
limite laminaire, et ne sont pas proposés.

Le moteur GPU 2D choisit désormais Vulkan (ou Metal) avant DirectX 12, comme
la soufflerie 3D : sous DirectX 12 (pilote NVIDIA), il divergeait dès le
premier pas, ce qui touchait l'interface, dont la fenêtre est en DirectX 12.
`SOUFFLERIE_GPU2D=dx12|vulkan|gl` impose un pilote.

**Blocage de la veine.** La veine 2D ne fait que 12 longueurs de haut :
l'écoulement accélère autour de la maquette et de son sillage, ce qui
surestimait la traînée d'un cylindre de 14 % à Re = 20 et de 8 % à Re = 100,
quelle que soit la résolution. Les coefficients (interface, `run`, `polar`)
sont désormais corrigés par la méthode d'Allen et Vincenti (vitesse effective
U(1 + ε), ε = blocage de volume + blocage de sillage) ; `--paroi 0` ou la case
de l'interface rendent les valeurs brutes :

| Cylindre, L = 32 | brut | corrigé | référence (veine infinie) |
| --- | --- | --- | --- |
| Re = 20, Cx | 2,347 | 2,077 | 2,05 à 2,09 |
| Re = 100, Cx | 1,442 | 1,331 | 1,33 à 1,35 |
| Re = 100, St | 0,175 | 0,170 | 0,164 à 0,167 |

La soufflerie 3D n'a pas encore cette correction.

**Lecture des efforts sans attente.** Le moteur GPU relève les efforts par
copies asynchrones au lieu de vider la file de calcul à chaque relevé : un
calcul de cylindre à L = 32 passe de 750 à 1 630 millions de mailles par
seconde (GTX 1070), à résultats identiques.

L'aile finie (onglet « Aile finie », commande `wing`) reprend la surface
portante de Weissinger non linéaire de l'application web : polaires de
section par la méthode des panneaux au Reynolds local (emplanture et bout,
interpolées en log Re), décambrage, traînée induite dans le plan de Trefftz ;
flèche, dièdre, effilement et vrillage. Une polaire complète prend 0,3 s.

La soufflerie 3D (onglet « Soufflerie 3D », commande `run3d`) reprend le
solveur D3Q19 WebGPU de l'application web : le WGSL est engendré à
l'identique (propagation et collision fusionnées, collision régularisée,
populations en demi-précision, rebond de Bouzidi, efforts accumulés sur le
GPU) ; sphère, cube, cylindre, aile NACA 0012 d'allongement 4, de 8 à
32 mailles par longueur de référence. Sur une GTX 1070 non bridée : environ
850 millions de mailles par seconde (sphère, 0,75 M mailles), soit 30 t* en
5 s. Le calcul 3D préfère Vulkan (ou Metal) : sous DirectX 12, ce gros shader
est environ 100 fois plus lent ; `SOUFFLERIE_GPU3D=dx12|vulkan|gl` impose un
pilote. Les soumissions sont courtes (≈ 3·10⁷ mises à jour de maille) pour
ne pas déclencher la réinitialisation du pilote Windows (2 s). Affichage :
vue 3D orbitale calculée sur le GPU par lancer de rayons dans le champ
(maquette ombrée, plan de coupe en vitesse ou vorticité, isosurfaces du
critère Q translucides teintées par la vitesse, cadre de la veine), comme
la version web ; Q est annulé à moins de deux mailles du corps, où la
couche limite mal résolue donne des valeurs parasites. Une coupe « à plat »
(vitesse, vorticité, Q) reste disponible.

L'interface 2D a aussi la polaire automatique (balayage d'incidence ou de
rotation, stabilisation puis mesure adaptatives, onglet « Polaire », copie
CSV) et l'échelle physique (air ou eau, longueur de référence : vitesse,
durées, efforts par mètre d'envergure, alerte de Mach réel). Le moteur GPU
2D peut stocker ses populations en demi-précision (`--half 1`, ou
« Précision GPU » dans l'interface) : deux fois moins de mémoire, écarts de
0,1 à 0,2 % sur les efforts moyens (cylindre Re = 100, NACA 0012 Re = 10⁶
avec Spalart–Allmaras), 3 % sur l'amplitude de Cz'.

## Loi de paroi en 3D

Comme en 2D, la soufflerie 3D passe à Spalart–Allmaras avec loi de paroi dès
Re = 5·10⁴ (`--turb auto|sa|off`) : ν_t = ν̃ f_v1 s'ajoute à la viscosité
de la collision, ν̃ est transporté par une seconde passe (advection décentrée
amont, diffusion, production et destruction de Spalart–Allmaras, ν̃ = κ u_τ y
au premier nœud), et la paroi glisse : chaque maille pariétale tire u_τ de la
loi de Spalding à partir de sa vitesse tangentielle et de sa distance à la
paroi, puis le rebond de Bouzidi la traite en paroi mobile à la vitesse de
glissement (efforts avec le terme correspondant). Distance à la paroi :
exacte à moins de 3 mailles du corps (triangles ou formes analytiques),
propagée au-delà par une transformée de chanfrein ; normales par gradient de
la distance exacte. La couche limite est supposée turbulente dès l'avant
(pas de transition en 3D).

| L = 16 mailles, Re = 10⁶ | Sous-maille seule | Spalart–Allmaras + loi de paroi | Référence |
| --- | --- | --- | --- |
| Sphère, Cx | 0,594 | 0,254 | ≈ 0,2 (après la crise de traînée) |
| Aile NACA 0012, allongement 4, α = 8°, Cz | 0,153 | 0,494 | ≈ 0,56 (ligne portante) |

À Re = 10⁵, la sphère réelle a une couche limite laminaire (Cx ≈ 0,47) :
Spalart–Allmaras, turbulent partout, donne alors le régime d'après la crise
(0,28). Test `loi_de_paroi_3d`.

A320 (STL) en approche : 70 m/s au niveau de la mer, α = 8°, Re = 1,7·10⁸,
L = 48 mailles (0,75 m), coefficients rapportés à la surface alaire
(122,6 m²) : CL = 0,22 et CD = 0,54 sans modèle de paroi, CL = 0,72 et
CD = 0,31 avec (ordre de grandeur réel en lisse : CL ≈ 0,9, CD de quelques
centièmes, davantage train sorti). La traînée reste trop forte à cette
résolution, et la coupe de vorticité garde un bruit en damier à τ ≈ 0,5.
Avec la collision régularisée d'ordre 3, les bandes absorbantes et une
grille fine ×2 (L = 32 + ×2, soit 64 mailles près de l'avion, 3,5 min) :
CL = 0,63, CD = 0,27, champ sans bruit.

## Stabilité et raffinement local en 3D

**Collision et bords.** La collision D3Q19 est régularisée récursivement
(ordre 3, Malaspinas 2015 : moments d'ordre 3 tirés de ρ, u et de la partie
hors équilibre d'ordre 2), comme en 2D. Aux Reynolds de maille élevés
(sous-maille active), des bandes absorbantes rappellent l'écoulement vers
l'amont près des bords latéraux périodiques, de l'entrée et dans la zone
éponge de sortie : les ondes de pression de la mise en vitesse ne rebondissent
plus dans la veine. Sphère à Re = 10⁶ (D = 16, loi de paroi) : vorticité
résiduelle en amont 2,63 → 0,004 U/L, fluctuation rapide de Cx
0,014 → 0,0006, Cx 0,254 → 0,214. Sans loi de paroi, la sous-maille est
amortie près de la paroi (Van Driest, longueur 8 L/√Re bornée à [1, 4]
mailles, comme en 2D) ; avec elle, elle reste entière : c'est elle qui
stabilise les premiers nœuds à τ ≈ 0,5.

**Raffinement local** (`run3d --raff 1|2`, interface : « Près de la
maquette »). Une ou deux grilles imbriquées autour de la maquette, deux fois
plus fines en espace et en temps (emprise : 0,6 L en amont, 1,5 L en aval,
0,6 L de côté, puis 0,35 / 0,7 / 0,4 L), couplées comme en 2D : moments
(ρ, u, Π hors équilibre remis à l'échelle, ν̃) interpolés en cubique sur
l'enveloppe de mailles fantômes et en linéaire dans le temps, recopie vers
la veine sous la zone raffinée ; efforts de la grille la plus fine.

| Cas | Veine uniforme | Raffinée | Calcul |
| --- | --- | --- | --- |
| Sphère, Re = 100 : D = 16 contre D = 8 + ×2 | 1,1472 | 1,1504 (+0,3 %) | 31 s contre 3 s |
| Aile, Re = 300, α = 8° : L = 32 contre 16 + ×2 (même calage) | Cx 0,2620, Cz 0,2628 | 0,2618, 0,2634 | 74 s contre 30 s |
| Aile, Re = 10⁶ : idem | Cx 0,0858, Cz 0,5600 | 0,0858, 0,5444 | — |

« Même calage » : la maquette d'une grille fine tombe à une demi-maille près
d'un nœud là où la veine la place à un quart ; à L = 32, l'aile (4 mailles
d'épaisseur) y est sensible (Re = 10⁶ : Cx 0,057 ou 0,086 selon le calage),
ce qu'atténue justement le raffinement (×4 : Cx 0,056, Cz 0,600).

**Transition** (`run3d --xtr 0.3`, interface : curseur « transition (x/c) »
quand la loi de paroi est active). Couche limite laminaire du bord d'attaque
jusqu'à la fraction xtr de la corde locale, mesurée dans chaque tranche
d'envergure comme l'étendue en x du corps : frottement de Blasius
(0,664/√Re_x) à la paroi et pas de production de Spalart–Allmaras autour,
comme la zone laminaire du 2D. Aile d'allongement 4 à 8°, Re = 10⁶,
L = 16 : Cx 0,1163 turbulente, 0,1107 laminaire partout (écart de frottement
attendu ≈ 0,006, test `transition_3d`) ; sur une grille plus grossière, le
glissement laminaire réduit aussi la traînée de pression. Le décollement
laminaire des corps épais n'est pas capturé (sphère à Re = 10⁵ : 0,25
turbulente, 0,21 laminaire, au lieu de 0,47) : pour ce régime sous-critique,
couper la loi de paroi (« Paroi : sans », `--turb off`) donne une couche
limite épaisse qui décolle tôt (sphère à Re = 10⁶ : Cx = 0,46).

## Comparaison avec VSPAERO et la mesure

VSPAERO 7.2.2 (OpenVSP 3.56, NASA, méthodes des tourbillons et des
panneaux, traînée de profil estimée) sur deux ailes publiques, contre la
soufflerie 3D avec loi de paroi et raffinement local. Données complètes :
`validation/comparaison_vspaero.csv`. FlightStream (commercial) n'a pas pu
être utilisé.

**Aile NACA TR-1208** (Schneider, 1951 : flèche 45°, allongement 8,02,
effilement 0,45, NACA 63A012, Re = 4·10⁶ sur la corde moyenne ; géométrie
du tutoriel OpenVSP exportée en STL, `--unite in --lref 16.672
--sref 2024.35 --marge 0.4`, L = 8 + ×4, soit 32 mailles par corde moyenne) :

| α | Mesure (figure 2, ±0,01) | VSPAERO VLM | Soufflerie |
| --- | --- | --- | --- |
| 2° | CL 0,136 | 0,132 | 0,161 |
| 4° | CL 0,257, CD ≈ 0,011 | 0,266, CD 0,015 | 0,314, CD 0,042 |
| 6° | CL 0,378 | 0,399 | 0,460 |
| 8° | CL 0,499, CD ≈ 0,021 | 0,531, CD 0,026 | 0,575, CD 0,052 |

**Aile rectangulaire NACA 0012, allongement 4**, Re = 10⁶ (pas de mesure ;
ligne portante d'Helmbold : CL ≈ 0,27 à 4°) :

| α | VSPAERO VLM | VSPAERO panneaux | Soufflerie L = 16 + ×2 | + ×4 |
| --- | --- | --- | --- | --- |
| 4° | 0,251 (CD 0,021) | 0,265 | 0,306 (CD 0,059) | 0,306 (0,032) |
| 8° | 0,500 (0,039) | 0,529 | 0,545 (0,086) | 0,600 (0,056) |
| 12° | 0,747 (0,070) | 0,787 | 0,958 (0,143) | 0,839 (0,100) |

Ce que cela montre :

- VSPAERO suit la mesure de TR-1208 à 5 % près sur la pente de portance
  (0,066 contre 0,060 par degré, le surcroît habituel d'une méthode non
  visqueuse).
- La soufflerie porte trop : +20 à 30 % sur la pente, au-dessus même de la
  méthode non visqueuse. Ce n'est pas le confinement de la veine (même aile
  en STL : Cz = 0,295 dans une veine serrée, 0,306 dans une veine 2,5 fois
  plus large) ; c'est la loi de paroi : à y⁺ de l'ordre de 10³, le
  glissement de Spalding devient presque total et le bord de fuite tient mal
  la condition de Kutta (déjà +7 % sur le NACA 2412 en 2D). Piste : une loi
  de paroi qui impose la contrainte pariétale sans laisser glisser la paroi
  au bord de fuite.
- La traînée est 2 à 4 fois trop forte et dépend de la taille de la veine
  (0,042 ou 0,062 selon la marge) : en 3D, seule la portance se compare
  aujourd'hui, et la traînée reste qualitative.

## Maquettes STL (3D)

Dans l'interface, onglet Soufflerie 3D : « Charger un STL… » (ou glisser le
fichier sur la fenêtre, ou `soufflerie --stl fichier.stl` au lancement). Les
mêmes réglages qu'en ligne de commande y figurent (axes, unité, incidence,
dérapage, références), avec les dimensions, le maître-couple et le blocage ;
« Échelle physique » fixe le Reynolds par la vitesse et donne les efforts
en newtons.

En ligne de commande :
`soufflerie run3d --stl fichier.stl --l 32 --re 5000` étudie une maquette
maillée (STL binaire ou ASCII, volume fermé). Mise en place :

| Option | Rôle | Défaut |
| --- | --- | --- |
| `--unite mm\|cm\|m\|in` | unité des coordonnées du fichier | mm |
| `--ecoulement +x` | axe du fichier dans le sens de l'écoulement (de l'avant vers l'arrière) | +x |
| `--haut +z` | axe du fichier vers le haut | +z |
| `--alpha 6` | incidence (°), avant vers le haut | 0 |
| `--beta 10` | dérapage (°), vent venant de la droite de la maquette | 0 |
| `--lref 250` | longueur de référence (unités du fichier) : `--l` mailles | longueur hors tout selon l'écoulement |
| `--sref 4.2e4` | surface de référence (unités du fichier²) | maître-couple calculé |
| `--vitesse 30` | vitesse (m/s) : fixe Re sur la longueur de référence, efforts en newtons | — |
| `--fluide air\|eau` | fluide de l'échelle physique | air |

La veine s'ajuste à la maquette : 3,5 M en amont, 7,5 M en aval et 2,5 M de
part et d'autre (M : la plus grande de L et des dimensions de la maquette),
au moins 6 L de large et de haut ; le blocage (maître-couple sur section de
veine) est affiché, avec un avertissement au-delà de 5 %. L'appartenance
d'une maille vient de la parité des traversées d'un rayon selon x, les
fractions de lien de Bouzidi d'une bisection, et l'affichage d'une distance
signée aux triangles. Une icosphère de 5 120 triangles donne la même
voxelisation que la sphère analytique et une traînée à 1 % près (test
`sphere_stl`) ; un même fuseau décrit avec d'autres axes (écoulement selon
−y, haut selon +x) donne des coefficients identiques.

**Parois minces.** Les liens du réseau qui traversent la surface du maillage
reçoivent une fraction de Bouzidi par intersection exacte avec les triangles,
y compris entre deux mailles fluides : une aile, une dérive ou une plaque plus
minces qu'une maille restent des parois (moins précises qu'à une résolution
qui les résout). Une plaque carrée face au vent de 0,08 maille d'épaisseur,
sans aucune maille solide, donne la traînée d'une plaque de 2 mailles
(Cx = 1,829 contre 1,831, test `paroi_mince_stl`) ; sur un A320 à 8°
(L = 32 mailles, Re = 2000), la portance passe de Cz = 0,10 à 0,24. À
l'affichage, les maillages sont épaissis de 0,6 pas de la grille de distance
pour que ces parois restent visibles.

**Voir le maillage.** Vue 3D, « Maquette » : *lisse* (forme exacte),
*grille* (forme exacte, grille du calcul tracée dessus), ou *mailles du
calcul* (ce que voit le solveur : mailles solides en cubes, et en orange les
mailles qui ne portent que des liens de paroi mince). En ligne de commande :
`--maquette lisse|grille|mailles` avec `--image`.

**Maillages imparfaits et gros fichiers.** À la lecture, les sommets sont
soudés (à 10⁻⁷ de la diagonale près), les triangles dégénérés et en double
retirés, et les arêtes comptées : trous (arêtes ouvertes), arêtes partagées
par plus de deux faces, orientations incohérentes. Un maillage ouvert est
signalé et l'appartenance y est votée à la majorité sur trois rayons (x, y,
z) au lieu de la parité d'un seul : une sphère privée d'un triangle sur 500
donne la même traînée que la sphère fermée (1,1470 contre 1,1472, D = 16
mailles). Quand plus de 20 % du volume du maillage disparaît à la
voxelisation (parties plus fines qu'une maille), un avertissement invite à
augmenter la résolution. Une sphère de 327 680 triangles se prépare en 4 s
(D = 16 mailles) ; dans l'interface, l'aperçu de mise en place d'un gros
maillage attend que le réglage soit stable.

## Validation 3D

Résultats bruts dans `validation/` ; `tools/validation3d.sh` refait les
calculs (GTX 1070, une pause entre deux calculs pour limiter l'échauffement).

**Sphère, traînée selon Re** (50 t*, moyenne sur les 20 derniers ; référence
de Schiller–Naumann et Clift) :

| Re | Cx, D = 12 mailles | Cx, D = 16 mailles | Référence | Écart à 16 mailles |
| --- | --- | --- | --- | --- |
| 20 | 2,848 | 2,838 | 2,610 | +8,7 % |
| 50 | 1,665 | 1,648 | 1,539 | +7,1 % |
| 100 | 1,168 | 1,145 | 1,094 | +4,7 % |
| 200 | 0,856 | 0,823 | 0,810 | +1,7 % |
| 300 | 0,739 | 0,695 | 0,691 | +0,5 % |
| 500 | 0,671 | 0,617 | 0,576 | +7,2 % |
| 1000 | 0,580 | 0,533 | 0,466 | +14 % |

(Avant la collision régularisée d'ordre 3 et l'amortissement pariétal de la
sous-maille : +18 % à Re = 500 et +32 % à Re = 1000.)

L'écart diminue partout quand la grille s'affine. Aux faibles Reynolds, il
vient surtout du confinement (veine de 6 D, bords périodiques), qui pèse
plus quand la zone d'influence visqueuse s'étend loin. Au-delà de Re = 300,
la couche limite n'a plus que 2 à 3 mailles à D = 16, et la sous-maille
(active dès τ < 0,51) ajoute de la traînée : il faut D = 24 à 32 pour ces
Reynolds. Au-delà de Re ≈ 210, la portance latérale apparaît comme attendu
(sillage oblique) ; le lâcher périodique à Re = 300 (St ≈ 0,13) n'est pas
encore détecté en 50 t*.

**Aile rectangulaire NACA 0012 d'allongement 4** : portance 3D (corde
16 mailles) comparée à la correction d'allongement de Helmbold appliquée à
la pente 2D du même profil, au même Reynolds, calculée par la soufflerie 2D
(`polar`, L = 16) :

| α | Re = 200, 3D | Helmbold + 2D | Re = 1000, 3D | Helmbold + 2D |
| --- | --- | --- | --- | --- |
| 4° | 0,103 | 0,149 | 0,111 | 0,101 |
| 8° | 0,249 | 0,289 | 0,214 | 0,180 |
| 12° | 0,383 | 0,415 | 0,330 | 0,350 |

À Re = 1000, l'aile 3D porte 10 à 19 % de plus que la correction
d'allongement (elle portait 1 à 7 % de moins avant l'amortissement pariétal
de la sous-maille, qui épaississait la couche limite) ; à 12°, le profil 2D
est déjà instationnaire et la comparaison ne tient plus. À Re = 200, l'aile
3D porte 10 à 30 % de moins que la correction d'allongement : la couche
limite épaisse et les tourbillons marginaux visqueux sortent du cadre de la
ligne portante. La surface portante alimentée par les panneaux n'est valable
qu'à partir de Re = 5·10⁴ (CL = 0,253 à 4°, deux fois et demie plus que
l'aile 3D à Re = 1000) : l'écart mesure l'effet du Reynolds, pas une erreur.
À incidence nulle, l'aile 3D garde un Cz de −0,03 à Re = 1000 (−0,007 à
Re = 200) alors qu'il devrait être nul : à examiner (décalage d'un quart de
maille de la maquette, ou sillage instationnaire).

## Installation

Installer Rust (<https://rustup.rs>), puis, dans ce dossier :

```sh
cargo run --release                     # interface graphique
cargo run --release -- --shape naca --naca 2412 --alpha 6 --re 1e6 --l 24 --graph cp   # interface réglée
cargo run --release -- run --shape naca --naca 2412 --alpha 6 --re 1000 --l 32
cargo run --release -- run --shape naca --naca 0012 --alpha 4 --re 1e6 --l 32 --turb sa
cargo run --release -- polar --shape naca --naca 0012 --re 2000 --l 24 --alpha -4:12:2 > polaire.csv
cargo run --release -- panel --naca 2412 --re 3e6 --alpha -4:18:1 > panneaux.csv   # instantané
cargo run --release -- wing --naca 2412 --ar 8 --taper 0.5 --re 1e6 --alpha -6:20:1 > aile.csv
cargo run --release -- wing --alpha 12 --span 1      # répartition en envergure
cargo run --release -- --tab wing       # interface, onglet aile finie
cargo run --release -- run3d --shape sphere --re 300 --l 16 --t 40 --coupe sillage.ppm --image vue.ppm
cargo run --release -- --tab 3d         # interface, onglet soufflerie 3D
cargo run --release -- bench --l 32     # vitesse des moteurs
cargo test --release                    # validation
```

Pour le processeur de la machine (AVX2, AVX-512) :
`RUSTFLAGS="-C target-cpu=native" cargo run --release`.

Options de calcul : `--shape cylinder|square|ellipse|plate|naca|triangle`,
`--naca 2412`, `--alpha 4`, `--re 1000`, `--l 32` (mailles par longueur),
`--t 40` (durée en temps convectifs), `--avg 15` (fenêtre de moyenne),
`--engine cpu|gpu`, `--threads 0` (tous les cœurs), `--spin 1`,
`--thickness 0.4`, `--collision rr|bgk`, `--les auto|on|off`, `--u 0.1`,
`--turb auto|sa|off` (auto : Spalart–Allmaras dès Re = 5·10⁴),
`--transition free|turb` (libre : position prédite par la couche limite
intégrale, paroi laminaire en amont ; turb : couche limite turbulente dès le
bord d'attaque). `panel` accepte aussi `--mach 0.5` (Prandtl–Glauert,
Kármán–Tsien, détection de la zone supersonique locale). Sans sous-commande,
les mêmes options règlent l'interface au démarrage, avec `--field
vorticity|speed|pressure|eddy` et `--graph history|cp|cf`.

`polar` ajoute les colonnes de la méthode des panneaux (Re ≥ 2·10⁴) ;
`panel` donne Cz, Cx, Cm (décrochage de Kirchhoff compris), le Cz attaché,
les transitions, le décollement et le Cp minimal, toutes incidences en
quelques millisecondes.

`wing` : `--naca`, `--ar` (allongement), `--taper` (effilement), `--sweep`
(flèche du quart de corde, °), `--dihedral`, `--twist` (vrillage en bout, °),
`--re` (rapporté à la corde moyenne), `--transition free|turb`, `--strips 40`
(bandes par demi-aile) ; CL, CD, CDi, CDp, Cm, Oswald e et finesse par
incidence, ou avec `--span 1` le Cl local, le Cl max et la charge le long de
l'envergure.

`run3d` : `--shape sphere|cube|cylinder|wing`, `--re` (10 à 2000 environ),
`--l 12` (mailles par longueur de référence), `--alpha` (cube, aile), `--t 40`
et `--avg 15` (temps convectifs), `--fp32 1` (simple précision), `--coupe
fichier.ppm` (vorticité ωy dans le plan de symétrie), `--image fichier.ppm`
(vue 3D 960 × 540 avec isosurfaces ; `--az 62 --el 24 --dist 6 --q 0.4`). Affiche Cx, Cz, Cy
(latéral), le Strouhal si l'écoulement est périodique, et la référence de
Schiller–Naumann pour la sphère.

Windows avec la chaîne GNU de Rust : si la compilation réclame `dlltool`,
mettre dans le PATH un dossier contenant `dlltool.exe` et `as.exe` de MinGW
(binutils), sans `gcc`, pour que Rust garde son propre éditeur de liens.

L'interface tourne dans un processus enfant surveillé : si un pilote
graphique la fait planter (souvent sans aucun message), le code d'erreur est
expliqué et un autre rendu est essayé (Windows : wgpu sur DirectX 12, puis
OpenGL, puis wgpu avec tous les pilotes ; ailleurs : wgpu puis OpenGL).
`--rendu glow|dx12|vulkan|wgpu` en impose un.

Linux : l'interface demande `libxkbcommon-x11` et un pilote Vulkan ou
OpenGL. Sans interface :
`cargo build --release --no-default-features --features gpu`.

## Structure

```
soufflerie-rs/
├── lbm/                bibliothèque de calcul
│   ├── src/d2q9.rs     réseau, équilibre, collision (RR, BGK, sous-maille)
│   ├── src/cpu.rs      moteur processeur (rayon, boucles vectorisées)
│   ├── src/gpu.rs      moteur GPU (wgpu) ; d2q9.wgsl, sa.wgsl, macros.wgsl
│   ├── src/sa.rs       Spalart–Allmaras, loi de paroi de Spalding
│   ├── src/panel.rs    méthode des panneaux, couche limite intégrale, Kirchhoff
│   ├── src/wall.rs     distributions pariétales Cp, Cf
│   ├── src/wing.rs     aile finie : surface portante non linéaire, Trefftz, polaires de section
│   ├── src/gpu3d.rs    soufflerie 3D : D3Q19 sur GPU (WGSL engendré), voxelisation, coupes, rendu 3D
│   ├── src/bodies3d.rs maquettes 3D (sphère, cube, cylindre, aile)
│   ├── src/geometry.rs mailles solides, liens de Bouzidi, distance et normale à la paroi
│   ├── src/shapes.rs   maquettes
│   ├── src/tunnel.rs   veine d'essai, mise en vitesse, turbulence, transition, coefficients
│   ├── src/stats.rs    moyennes, périodes, Strouhal
│   └── tests/          validation, panneaux, aile finie, 3D, comparaison GPU/CPU
├── validation/         résultats de validation 3D (CSV)
├── tools/validation3d.sh  calculs de validation 3D
└── app/                application : interface (egui) et commandes
    ├── src/main.rs     ligne de commande (run, polar, panel, wing, run3d, bench)
    ├── src/gui.rs      interface graphique, onglets
    ├── src/wing_tab.rs onglet aile finie (vue 3D, envergure, polaire)
    ├── src/tab3d.rs    onglet soufflerie 3D (fil de calcul GPU, vue 3D, coupe, efforts)
    └── src/worker.rs   fil de calcul de l'interface (image, traceurs)
```

## Pistes

- Soufflerie 3D : rendu directement dans la fenêtre (texture partagée avec
  egui) plutôt que par relecture de l'image, pour les grandes fenêtres.
- Stockage en demi-précision et propagation en place sur GPU (deux à quatre
  fois moins de mémoire), raffinement local de maillage.
- Version mobile : egui fonctionne sur Android ; il faudrait l'empaqueter
  (cargo-apk / xbuild). Sur téléphone, l'application web reste la plus
  simple à installer.
