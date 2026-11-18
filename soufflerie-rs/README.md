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

## Suites de la comparaison

**Outils ajoutés.**

- `run3d --envergure fichier.csv` : répartition de la portance en
  envergure (c·cl/c_ref par tranche, efforts de la grille la plus fine sur
  deux temps convectifs de plus). Maquette maillée : l'effort de chaque
  triangle est réparti sur son étendue en envergure selon la densité d'aire
  du triangle (`span_distribution`). Posé au centre du triangle, il faisait
  apparaître des creux et des bosses parasites : sur l'aile AR 4, maillée
  par stations de 0,087 m, chaque tranche de 1/32 m recevait une ou deux
  rangées de triangles (aire par tranche de 1 à 1,5). La répartition est
  maintenant lisse (écart à la moyenne glissante < 0,002 hors des saumons).
- Moment de tangage Cm des maquettes STL (et avec `--xref`), autour du
  point situé à xref longueur de référence de l'avant (quart de corde par
  défaut), positif à cabrer.
- `--sol <garde au sol>` : sol mobile (plancher qui défile à la vitesse
  amont, plafond glissant) pour les véhicules ; pas encore avec le
  raffinement local.
- `--marge` : veine resserrée autour d'une maquette très allongée.
- Raffinement 3D jusqu'à ×8 (`--raff 3`, interface « ×8 ») : aile
  d'allongement 4 à 4°, L = 8 + ×8 → Cz 0,298, Cx 0,032, contre L = 16 + ×4
  (même résolution près de l'aile) → 0,306, 0,032. À L = 16 + ×8 (128
  mailles par corde), le calcul dépasse 30 min sur une GTX 1070.
- `tools/comparaison_vspaero.sh` : construit l'aile de TR-1208 dans
  OpenVSP (`tools/vsp/`), lance VSPAERO, exporte le STL, lance la soufflerie
  à 0–8° et écrit le tableau mesure / VSPAERO / soufflerie.

**Répartition en envergure** (TR-1208, α = 4,7°,
`validation/tr1208_envergure.csv`) : l'excès de portance de la soufflerie
est réparti sur toute la mi-envergure (c·cl/(CL·c̄) ≈ 1,2 entre 25 et 50 %
de la demi-envergure, contre 1,1 mesuré et 1,07 pour VSPAERO), pas
concentré au saumon ni à l'emplanture. Les efforts par maille sont bruités
par l'escalier des bords : il faut des tranches d'environ 0,13 demi-envergure.

**Loi de paroi : ce qui a été essayé** (aile d'allongement 4, α = 4°,
référence ≈ 0,26 à 0,27) :

| Variante | Cz |
| --- | --- |
| Glissement de Spalding à la première maille (actuel) | 0,306 |
| Sans glissement dans les 2 ou 4 dernières mailles (bord de fuite) | 0,299 / 0,288 |
| Vitesse échantillonnée à 2 ou 3 mailles de la paroi | 0,391 / 0,345 |
| Sans glissement du tout (Spalart–Allmaras gardé) | 0,147 |

Le glissement encadre la bonne valeur mais en fait trop, et la cause n'est
ni le bord de fuite seul ni une rétroaction de la première maille. Une
paroi adhérente avec une force tangentielle qui rend au fluide l'excédent
de frottement du rebond diverge (la composante tangentielle des liens d'une
paroi en escalier mêle pression et frottement) ; plafonnée, elle donne
Cz 0,205 et une traînée de 0,17 : un réglage, pas un modèle. Aucune de ces
variantes n'est gardée.

Les efforts eux-mêmes sont justes : un bilan de quantité de mouvement sur
une boîte autour de l'aile (flux de ρuu + p à travers les faces, champ de
la veine) donne Cz 0,314 à 0,321 contre 0,308 par l'échange aux parois à
Re = 10⁶, et 0,212 à 0,216 contre 0,213 à Re = 1000 (`Gpu3D::macros`).
L'excès de portance est donc dans l'écoulement calculé ; il existe aussi
sans loi de paroi (Re = 1000, sous-maille seule : +19 % sur la correction
d'allongement), ce qui oriente vers la résolution du bord de fuite et du
sillage près de l'aile plutôt que vers la seule loi de paroi.

**Aile infinie : le 3D contre le 2D.** Une aile NACA 0012 qui occupe
toute la largeur périodique de la veine 3D (veine haute de 12 cordes,
L = 32, Re = 10⁶, α = 4°) porte Cz 0,547 contre 0,505 pour la soufflerie
2D au même maillage (panneaux et couche limite : 0,471). Avec la
sous-maille amortie près de la paroi comme en 2D, elle porte exactement
0,505. L'aile d'allongement 4 porte alors 0,321, soit 0,63 fois l'aile
infinie, le rapport de la ligne portante d'Helmbold (0,62) : le passage
du 2D au 3D est juste, et l'excès par rapport à VSPAERO et à la mesure
vient du profil lui-même, qui porte déjà 7 % de trop en 2D avec la loi de
paroi (épaisseur de couche limite et décambrage visqueux mal rendus par
le glissement). Cet amortissement n'est pas gardé par défaut : sans
plancher, il diverge (aile finie, L = 32, 8°) ; avec un plancher, il
augmente la portance de l'aile finie (0,621 → 0,655). Prochaine étape :
corriger la portance du profil en 2D (comparaison panneaux + couche limite
déjà en place), le 3D suivra. Dans ce montage, la traînée de l'aile
infinie (0,041 à 0°, contre 0,009 en 2D, une fois les bandes latérales
retirées) reste à comprendre.

**Traînée et taille de la veine.** Sur l'aile d'allongement 4, la traînée
est établie dès t* = 5 et ne bouge plus jusqu'à t* = 40, mais vaut 0,042
dans une veine serrée et 0,062 dans une veine 2,5 fois plus large : un biais
de domaine stable, à l'opposé d'un effet de blocage (la veine plus large
traîne plus). Expliqué depuis : c'est la place du bord de fuite dans la
maille (voir « Traînée 3D » ci-dessous).

**Sphère à travers la crise de traînée** (D = 16 + ×2) :

| Re | Loi de paroi | Sans loi de paroi | Achenbach (mesure) |
| --- | --- | --- | --- |
| 10⁵ | 0,151 | 0,334 | ≈ 0,47 |
| 3·10⁵ | 0,136 | 0,336 | ≈ 0,1 (crise) |
| 10⁶ | 0,123 | 0,334 | ≈ 0,13 |
| 3·10⁶ | 0,114 | 0,339 | ≈ 0,18 |

La loi de paroi donne le bon ordre de grandeur au-delà de la crise, mais
ni la crise ni la remontée ; sans elle, la sphère reste sous-critique.

**Corps d'Ahmed à 25° sur sol mobile** (60 m/s, Re = 4,3·10⁶, garde au sol
50 mm, sans pieds ; mesure : Cx = 0,285) : Cx 0,78 à L = 32 (garde au sol
d'une maille et demie : soubassement bouché, Cz = 1,5), 0,60 à L = 64
(Cz = 1,0, lâcher périodique). Hors de portée sans raffinement près du sol.
Une sphère à un diamètre du sol garde une portance faible (Cz = 0,03) :
la condition de sol elle-même est saine.

**Transition sur TR-1208** (α = 4°, laminaire jusqu'à 30 % de la corde) :
Cz 0,317 au lieu de 0,314, Cx 0,0396 au lieu de 0,0415 : effet faible,
l'excès de portance ne vient pas de là.

**Moment de tangage.** Le point par défaut (quart de la longueur de
référence depuis l'avant de la maquette) n'est pas le quart de la corde
moyenne d'une aile en flèche : pour comparer à une mesure, donner `--xref`
(TR-1208 à 4° : Cm = −0,58 autour de ce point par défaut, très en avant du
foyer).

**Non fait.** Ailes de Sivells : les géométries et mesures sont dans des
rapports NACA scannés, non extraits.

## Traînée 3D : placement et intégration sur la surface

**Placement dans la maille.** Sur l'aile d'allongement 4 (STL, L = 32 sans
raffinement, Re = 10⁶, α = 4°), déplacer la maquette d'une fraction de
maille en x fait varier Cz de 0,21 à 0,44 et Cx de 0,032 à 0,059, avec un
écoulement établi et stationnaire dans chaque cas. C'est l'erratique
« biais de domaine » noté plus haut : changer la longueur de la veine
déplaçait en fait la maquette dans la maille. L'effet est périodique d'une
maille, existe aussi à Re = 2000, sans loi de paroi, en simple précision,
sans parois minces et sur l'aile analytique : il vient de la place du bord
de fuite par rapport aux nœuds. Quand il tombe juste avant un nœud, deux
rangées de mailles pariétales disparaissent en envergure et la portance
chute d'un quart : à 32 mailles par corde, une marche d'une maille au
bord de fuite agit comme un volet Gurney. La portance est stable (plateau)
quand le bord de fuite est entre 0,4 et 0,8 maille après un nœud, et le
raffinement ×2 converge vers ce plateau. Le point le plus en aval des
maquettes (et le bord de fuite de l'aile analytique) est donc calé à
0,7 maille après un nœud, fraction qui reste dans le plateau sur les grilles
fines (×2 : 0,4 ; ×4 : 0,8 ; ×8 : 0,6). Même aile, veine allongée de 0,6 à
2 fois vers l'amont : Cz 0,300 à 0,303 au lieu de 0,21 à 0,44.
`MeshSetup::offset` (`run3d --decalage dx`) décale la maquette pour
mesurer cette sensibilité. Une aile en flèche ou effilée, dont le bord de
fuite coupe les nœuds à toutes les fractions, y est peu sensible.

**Efforts par intégration sur la surface.** L'échange de quantité de
mouvement aux parois surestime nettement la traînée avec la loi de paroi.
Sur l'aile AR 4 à 0° (L = 16 + ×2), il donne Cx = 0,064 alors que la
surface donne 0,018 ; relevé maille par maille, l'excès n'est pas réparti
le long de la corde mais concentré sur trois rangées en envergure, aux
marches de l'escalier, là où la surface coupe un plan de la grille sous un
angle rasant. Les liens presque parallèles à la paroi y renvoient la
pression et la quantité de mouvement longitudinale dans leur propre
direction : une marche d'une maille agit comme une petite rugosité, que
l'écoulement ressent réellement (le bilan de quantité de mouvement sur une
boîte donne la même traînée). L'excès décroît avec le raffinement (0,064,
0,041, 0,034 à ×2, ×4, ×8) sans disparaître. Pour les maquettes STL, les
efforts sont donc intégrés sur les triangles du maillage
(`Gpu3D::surface_forces`, grille la plus fine) : −p_p n dA, plus τ_p dA avec
τ_p = ρ u_τ² de la loi de Spalding sur la vitesse tangentielle relevée à
une maille le long de la normale. `run3d` en fait la moyenne sur la fenêtre
de moyenne (un relevé par demi-temps convectif), en tire aussi Cm et la
répartition en envergure (`--envergure`), et affiche l'échange pour
comparaison ; l'interface trace ces coefficients pour les maquettes STL.
Les formes analytiques gardent l'échange.

**Problèmes connus de la méthode et remèdes retenus.** La littérature sur
les lois de paroi en LBM sur grille cartésienne signale :
- *Le décalage de la zone logarithmique* (Kawai et Larsson, *Phys. Fluids*
  24, 2012) : au premier nœud, l'erreur numérique et de sous-maille fausse
  la vitesse qui alimente la loi de paroi. Remède retenu : la loi de paroi
  lit la vitesse une maille plus loin le long de la normale (interpolation
  trilinéaire, repli sur le nœud si un voisin est solide), comme le variant
  « ESG » d'Asmuth et al. (*Phys. Fluids* 33, 2021). Cz baisse de 4 à 5 %,
  vers la mesure.
- *La pression lue hors de la paroi* : relevée à une maille, la forte
  dépression du bord d'attaque est atténuée et sa succion manque. Remède
  retenu : p_p extrapolée linéairement depuis les relevés à une et deux
  mailles. La traînée baisse de 10 à 20 %.
- *L'incohérence entre la vitesse de glissement et la viscosité au premier
  nœud* : la vitesse de glissement suppose ν + κ u_τ y, alors que la maille
  y ajoute Spalart–Allmaras et Smagorinsky. Calculer le glissement avec la
  viscosité effective du nœud diverge (rétroaction entre glissement,
  cisaillement et sous-maille).
- *L'échange de quantité de mouvement inverse* (Asmuth et al. 2021), qui
  choisit la vitesse de paroi pour que l'effort tangentiel cédé par la
  maille vaille τ_p A (A = 6 Σ w_i |c_i·n|, aire exacte en moyenne pour un
  plan quelconque) : essayé, il ne réduit pas l'excès des marches, qui
  vient de la pression statique renvoyée par les liens rasants, et le
  compenser aussi exige des vitesses de paroi instables.
- *La reconstruction du premier nœud* (Malaspinas et Sagaut, *J. Comput.
  Phys.* 275, 2014 ; Degrigny et al., *Comput. Fluids* 227, 2021), qui
  impose vitesse et contrainte du modèle au nœud : essayée avec le rebond
  de Bouzidi gardé pour les liens, elle n'est pas conservative ici (le
  rebond retire à chaque pas la quantité de mouvement que le nœud remet) et
  fausse le champ de pression (Cx de pression négatif).

| Cas | Échange Cz / Cx | Surface Cz / Cx (pression + frottement) | Référence |
| --- | --- | --- | --- |
| Aile AR 4, 4°, L = 16 + ×2 | 0,309 / 0,062 | 0,263 / 0,030 (0,023 + 0,007) | VSPAERO 0,251 à 0,265 / 0,021 à 0,022 |
| Aile AR 4, 4°, L = 16 + ×4 | 0,314 / 0,044 | 0,289 / 0,030 (0,023 + 0,008) | id. |
| TR-1208, 4°, L = 8 + ×4 | 0,314 / 0,041 | 0,261 / 0,021 (0,016 + 0,006) | mesure 0,257 / 0,011 ; VSPAERO 0,266 / 0,015 |
| TR-1208, 8°, L = 8 + ×4 | 0,581 / 0,053 | 0,489 / 0,053 (0,048 + 0,005) | mesure 0,499 / 0,021 ; VSPAERO 0,531 / 0,026 |

Ces valeurs (pression relevée à une maille, loi de paroi alimentée au
premier nœud) sont celles de la première version. Avec la loi de paroi
alimentée une maille plus loin et la pression extrapolée à la paroi :

| Cas | Échange Cz / Cx | Surface Cz / Cx (pression + frottement) | Référence |
| --- | --- | --- | --- |
| Aile AR 4, 4°, L = 16 + ×4 | 0,308 / 0,048 | 0,288 / 0,024 (0,017 + 0,007) | VSPAERO 0,251 à 0,265 / 0,021 à 0,022 |
| TR-1208, 4°, L = 8 + ×4 | 0,304 / 0,045 | 0,275 / 0,018 (0,013 + 0,005) | mesure 0,257 / 0,011 ; VSPAERO 0,266 / 0,015 |
| TR-1208, 8°, L = 8 + ×4 | 0,551 / 0,056 | 0,506 / 0,047 (0,042 + 0,005) | mesure 0,499 / 0,021 ; VSPAERO 0,531 / 0,026 |
| TR-1208, 8°, L = 8 + ×8 | 0,593 / 0,043 | 0,570 / 0,041 (0,035 + 0,006) | id. |

La traînée baisse de 10 à 20 % et rejoint VSPAERO à 4° sur l'aile AR 4 ; la
portance reste à quelques pour cent de VSPAERO à ×4. Le frottement (0,005
à 0,008) a l'ordre de grandeur d'une plaque plane turbulente. Restent deux
écarts, qui sont dans l'écoulement calculé et non dans la mesure des
efforts (l'échange les porte aussi) :
- *Traînée de pression à 8°*, encore double de la mesure. Elle décroît
  avec le raffinement (0,042 à ×4, 0,035 à ×8) : le rayon du bord
  d'attaque ne couvre qu'une à deux mailles, et la succion qu'il porte
  (∝ CL²) n'est pas résolue.
- *Portance à 8° et ×8*, 14 % au-dessus de la mesure et au-dessus de
  VSPAERO non visqueux, et qui croît avec le raffinement. Cause trouvée et
  corrigée depuis : voir « Loi de paroi 3D : couche pariétale cohérente ».

Le moment de tangage et la répartition en envergure des maquettes STL
viennent aussi de l'intégration sur la surface (effort de chaque triangle).
Pour les formes analytiques, ils viennent de l'échange aux parois, dont les
efforts par maille retirent la part de la pression de référence (ρ₀ c_s²,
nulle en somme mais pas maille par maille), ce qui supprime les pics
parasites aux saumons.

## Loi de paroi 3D : couche pariétale cohérente

**Le défaut, vu en 2D.** NACA 0012, Re = 6·10⁶, α = 10°, couche limite
turbulente (référence NASA TMR : mesure de Ladson CL ≈ 1,08, Spalart–Allmaras
sur grille fine CL = 1,09, CD = 0,0123). Avec la loi de paroi d'origine, la
soufflerie 2D donne CL = 1,12 à 64 mailles par corde, 1,29 à 128 (au-dessus
de la solution non visqueuse, 1,19), puis un écoulement décollé et
instationnaire à 256 (CL ≈ 0,95, CD = 0,045) : la solution s'éloigne de la
référence quand la grille s'affine. Le Cp pariétal montre pourquoi : un pic
de pression suivi d'une dépression, à intervalles réguliers, là où la paroi
traverse une rangée de nœuds sous un angle rasant. Ce sont les marches de
l'échange de quantité de mouvement relevées plus haut en 3D ; sur
l'intrados, elles agissent comme de petits volets Gurney. Elles
disparaissent en laminaire (paroi adhérente) : elles viennent de la loi de
paroi.

**La cause.** Le glissement u_s = u_n − u_τ² y_n / (ν + κ u_τ y_n) vaut
u_n − u_τ/κ à peu près quelle que soit la distance y_n du nœud à la paroi :
deux nœuds pariétaux voisins, l'un à 0,05 maille et l'autre à 1 maille de
la paroi, en déduisent des vitesses de paroi qui diffèrent de tout l'écart
de vitesse résolu entre ces deux hauteurs (≈ u_τ ln(20)/κ). La vitesse
extrapolée à la paroi saute donc à chaque traversée d'une rangée de nœuds.

**Le remède retenu** (`gpu3d.rs`, shader de propagation) :
- la loi de paroi est alimentée par la vitesse relevée à une hauteur fixe
  h = 2 mailles au-dessus du pied de la normale (Kawai et Larsson 2012),
  avec repli sur le nœud si le relevé touche une maille solide ou ne tombe
  pas à la hauteur h de la paroi la plus proche (bords minces, coins) ;
- la couche [0, h] a la viscosité turbulente uniforme κ u_τ h : profil
  linéaire de pente τ_p / (ν + κ u_τ h), et chaque nœud glisse à
  u_n − pente · y_n. La vitesse extrapolée à la paroi ne dépend plus de la
  place du nœud dans la maille, et la contrainte transmise est τ_p partout ;
  le premier nœud reçoit la même viscosité (collision) et ν̃ = κ u_τ h
  (Spalart–Allmaras) ;
- le frottement intégré sur la surface est relevé à la même hauteur h.

**Stabilité.** Près de la paroi, ce glissement devient presque libre, et
une paroi qui suit la vitesse de son propre nœud n'a pas de rappel : le
calcul diverge à 64 mailles par corde. Trois mesures, sans effet sur les
efforts en régime établi :
- La vitesse de glissement est relaxée sur une dizaine de pas
  (u_w ← u_w + 0,1 (u_s − u_w)). En régime établi, sans effet : 0,1 ou
  0,02 donnent les mêmes efforts à 0,3 % près.
- Le glissement est borné par celui que donne aussi la vitesse relevée à h
  (profil monotone ; d'abord u_s ≤ u_h − pente · h, désormais
  u_s ≤ min(u_n, u_h) − pente · y_n, voir « Couche limite de la plaque
  plane »). Les
  deux estimations coïncident quand le profil est établi ; sans la borne,
  un nœud du bord d'attaque de TR-1208 à 8° s'emballait lentement (vitesse
  maximale stable à 1,37 U jusqu'à t* = 6,5, puis multipliée par trois en
  3 t*).
- Aux arêtes (saumon plat, bord de fuite mince), la normale d'une maille
  n'est celle d'aucune des faces que touchent ses liens, et sa vitesse
  « tangentielle » traverse ces faces : la paroi soufflait. Chaque lien
  garde désormais la normale de la face qu'il touche (triangle du maillage,
  ou gradient de la distance juste devant la paroi pour les formes
  analytiques), et la vitesse de paroi est projetée dans le plan de cette
  face.

- Aux arêtes vives (faces touchées par les liens d'une maille à plus de
  60° l'une de l'autre : saumon plat, bord de fuite), l'écoulement décolle et
  la loi de paroi, faite pour une couche limite attachée le long d'une face,
  ne s'applique pas : la paroi y reste adhérente, avec la viscosité de
  Spalart–Allmaras et de la sous-maille. Sans cela, la vitesse de
  contournement de l'arête (singulière en non visqueux) croissait jusqu'à
  la divergence (aile AR 4 à 8°, coin bord d'attaque / saumon). Les facettes
  d'une surface lisse font quelques degrés entre elles, une arête vive 90° à
  165° : tout seuil entre 30° et 80° classe les mêmes mailles.

**Essayés et écartés.**
- Faire glisser la paroi à la vitesse relevée à h (u_h − u_τ/κ) : stable,
  mais la paroi entraîne alors le fluide plus vite que le nœud ne va
  (Cz 0,35 à 0,37 contre 0,30 en glissement libre sur l'aile AR 4).
- ν̃ du premier nœud non imposé (Spalart–Allmaras seul) : la turbulence ne
  s'établit pas sur ces grilles, la couche limite décolle. Plafonné par la
  valeur de Clauser 0,0168 U_e δ* : sans effet sur la portance.

**Contrôle non visqueux.** Paroi en glissement libre (vitesse du nœud), aile
analytique AR 4 à 4°, L = 16 + ×2, échange aux parois : Cz = 0,299 ; la
nouvelle loi de paroi donne 0,301 (+ ×2) et 0,310 (+ ×4) : elle n'ajoute
plus de portance au glissement libre (l'ancienne en ajoutait à mesure que la
grille s'affinait).

**Résultats** (efforts intégrés sur la surface ; « avant » : loi de paroi
de la section précédente) :

| Cas | Avant Cz / Cx | Couche pariétale cohérente Cz / Cx | Référence |
| --- | --- | --- | --- |
| Aile AR 4, 4°, L = 16 + ×2 | 0,263 / 0,030 | 0,235 / 0,023 | VSPAERO 0,251 à 0,265 / 0,021 à 0,022 |
| Aile AR 4, 4°, L = 16 + ×4 | 0,288 / 0,024 | 0,262 / 0,022 | id. |
| Aile AR 4, 8°, L = 16 + ×4 | 0,600 / 0,056 (aile analytique, échange) | 0,528 / 0,040 | VSPAERO 0,501 à 0,529 / 0,039 à 0,045 |
| TR-1208, 4°, L = 8 + ×4 | 0,275 / 0,018 | 0,258 / 0,015 | mesure 0,257 / 0,011 ; VSPAERO 0,266 / 0,015 |
| TR-1208, 4°, L = 8 + ×8 | — | 0,268 / 0,015 | id. |
| TR-1208, 8°, L = 8 + ×4 | 0,506 / 0,047 | 0,493 / 0,044 | mesure 0,499 / 0,021 ; VSPAERO 0,531 / 0,026 |
| TR-1208, 8°, L = 8 + ×8 | 0,570 / 0,041 | 0,526 / 0,036 | id. |

La portance de TR-1208 est à 1 % de la mesure à ×4 et à 4 à 5 % à ×8, au
niveau de VSPAERO non visqueux ; celle de l'aile AR 4 tombe dans la
fourchette VSPAERO (tourbillons / panneaux). Elle croît encore un peu avec
le raffinement (+4 à 7 % de ×4 à ×8, contre +13 % avant) : à 64 mailles par
corde, la couche limite du bord de fuite ne couvre qu'une à deux mailles, et
son décambrage visqueux (5 à 6 % ici : mesure contre VSPAERO) n'est pas
encore rendu. À ×2 (32 mailles par corde), le bord de fuite mince, couvert
de mailles d'arête adhérentes, décambre au contraire un peu trop (0,235).

La traînée est juste à 4° (0,015 : VSPAERO, et la mesure à 0,004 près, la
couche limite étant supposée turbulente dès le bord d'attaque) mais reste
trop forte à 8° (0,036 à ×8 contre 0,021 à 0,026). Le relevé de Cx par
tranche de corde (TR-1208, 8°, ×4) la situe : sur les 20 premiers pour cent,
la succion de l'extrados ne donne qu'une poussée de −0,015 (Cp ≈ −0,3 juste
derrière le bord d'attaque à mi-envergure, au lieu d'un pic de −2 à −3),
et l'arrière de l'extrados, qu'elle devrait compenser, ajoute +0,026. Le
rayon du bord d'attaque du NACA 63A012 (≈ 1 % de la corde) tient dans moins
d'une maille à ×8 : c'est une limite de résolution, qui décroît en affinant
(traînée de pression hors traînée induite : 0,028 à ×4, 0,019 à ×8), et non
plus un défaut de la loi de paroi.

**Correction.** Les portances de ce tableau sont sous-estimées de 5 à 10 % :
la pression pariétale y était relevée à une ou deux mailles de la paroi,
trop loin pour la dépression du bord d'attaque. Le bilan de quantité de
mouvement sur une boîte donne une portance plus forte de 6 à 12 % que les
références ; voir la section suivante.

La soufflerie 2D (processeur et GPU) garde pour l'instant la loi de paroi
d'origine, identique à celle de l'application web : le défaut décrit
ci-dessus y reste au-delà d'environ 64 mailles par corde. Un essai sur le
GPU 2D (u_τ et glissement tirés de la vitesse relevée à h) ramène le
NACA 0012 à 10° de CL = 1,29 à 1,12 (128 mailles par corde) et supprime le
décollement à 256 (CL = 1,13 à 1,17), mais c'est la variante qui entraîne
le fluide en 3D ; la formulation retenue en 3D (glissement depuis le nœud,
relaxé, arêtes adhérentes) reste à porter sur les deux moteurs 2D.

## Bord d'attaque, pression pariétale et bilans de quantité de mouvement

Cinq pistes pour la traînée trop forte à 8° (TR-1208), examinées une à une.

**Arêtes vives au bord d'attaque : non.** Une maille du bord d'attaque dont
les liens touchent des faces à plus de 60° l'une de l'autre perd la loi de
paroi (arête vive) ; avec un rayon de moins d'une maille, on pouvait le
craindre. Relevé : 1 maille sur 46 à ×4, aucune sur 57 à ×8 dans les
premiers 0,5 % de corde. Les mailles d'arête sont au bord de fuite.

**Demi-précision : non.** Le même calcul en simple précision (TR-1208, 8°,
×4) donne les mêmes Cp à 0,006 près et la même traînée par tranche de corde
à 0,0001 près.

**Pression pariétale : relevée trop loin.** Le champ contient la dépression
du bord d'attaque, mais près de la paroi : à 2 % de corde sur l'extrados,
Cp = −0,66 à un quart de maille et −0,27 à une maille (×4). L'intégration
sur la surface relevait la pression là où les huit nœuds voisins sont
fluides, soit à une à trois mailles. Elle extrapole désormais depuis ½ et
1 maille, interpolée sur les seuls nœuds fluides du côté de la paroi (pas
ceux de l'autre face d'une paroi mince). Le résultat dépend du choix des
hauteurs, et c'est la portance qui bouge le plus :

| TR-1208, 8°, Cz / Cx | ×4 | ×8 |
| --- | --- | --- |
| relevés à ½ et 1 maille (retenu) | 0,551 / 0,044 | 0,575 / 0,039 |
| ¾ et 1½ | 0,519 / 0,038 | 0,548 / 0,033 |
| 1 et 2 | 0,496 / 0,042 | 0,528 / 0,034 |
| ancien relevé (voisinage entièrement fluide) | 0,493 / 0,044 | 0,526 / 0,036 |
| bilan sur une boîte (ci-dessous) | 0,545 / 0,049 | 0,554 / — |

Essayée et écartée : la correction de courbure de l'équilibre normal
(dp/dn = ρ u_t² κ, κ courbure du maillage dans la direction de
l'écoulement, p_p = p(½) − ρ u_t² ln(1 + κ/2)). Elle ramène la portance à
0,510 quand le bilan sur une boîte donne 0,545 : le modèle d'écoulement
concentrique ne vaut pas quand le rayon tient dans une maille.

**Bilan sur une boîte : la portance de référence du calcul.**
`Gpu3D::box_forces_vol` intègre −∮ (ρ u u·n + p n) dA sur une boîte autour
de la maquette, dans la grille la plus fine qui l'entoure, indépendamment
du traitement des parois. La portance ne dépend ni de la taille de la boîte
(0,5446 à 0,5457 de 0,05 à 0,25 longueur de marge) ni du temps. La traînée
croît linéairement avec le volume de la boîte (0,058 à 0,064) : la veine a
un léger gradient de pression longitudinal qui pousse aussi le fluide de la
boîte. `run3d` affiche ce bilan pour deux boîtes (0,1 et 0,15 longueur) et
la traînée extrapolée à volume nul (0,048 contre 0,044 sur la surface à
×4, 8°).

**Raffinement du bord d'attaque (`--ba 1`).** Une grille deux fois plus
fine encore, limitée aux faces tournées vers l'amont (n_x < −0,5), avec 10 %
de longueur de référence vers l'aval. Les efforts de la surface viennent,
triangle par triangle, de la grille la plus fine qui l'entoure ; l'échange
aux parois combine les deux grilles. Sur TR-1208, 8°, ×4 : la traînée de
pression passe de 0,0387 à 0,0332, celle du raffinement complet ×8 (0,0335) :
la baisse de ×4 à ×8 vient du bord d'attaque. La portance (boîte) passe de
0,545 à 0,549. Une seule boîte alignée sur les axes couvre toutefois le
bord d'attaque d'une aile en flèche sur presque toute sa corde : 9,4 M
mailles à ×8 (la grille ×4 en a 2,5 M), 54 M à ×16, au-delà de la mémoire du
GPU ; `--ba` refuse au-delà de quatre fois la grille parente.
*Correction :* une grille partielle qui coupe la maquette fausse
l'écoulement (section « Grilles des bords » plus bas) ; l'accord de
TR-1208 et la conclusion sur la traînée du bord d'attaque ne tiennent pas.

**Traînée en champ lointain (`--sillage 1`).** `Gpu3D::wake_drag` décompose
la traînée dans un plan transversal (Betz, Maskell) : perte de pression
totale dans le sillage visqueux (ν̃ > 10 ν), perte hors du sillage
(dissipation numérique de l'écoulement extérieur), énergie de l'écoulement
transversal (traînée induite). La pression totale de référence est la
médiane du même plan hors du sillage : les ondes de pression de la veine,
uniformes sur la section, s'annulent. Le bilan complet entre deux plans
n'est pas exploitable : les bandes absorbantes des bords latéraux, qui
ramènent l'écoulement vers l'amont, exercent leurs propres efforts (la
portance du plan ne vaut que 0,28 à 0,37 au lieu de 0,55). TR-1208, 8° :
traînée induite 0,011 à ×4, 0,015 à 0,019 à ×8 (CL²/(πA) = 0,012 pour un
allongement 8) ; profil (perte dans le sillage) 0,054 à ×4, 0,034 à 0,037 à
×8 ; hors sillage ±0,005, sans signe constant. Deux vérifications : la
même intégration sur la seule grille de la veine donne les mêmes
composantes à 0,005 près (le résultat ne vient pas du raccord des grilles) ;
avec la pression totale de référence prise dans le plan amont, la perte
« hors sillage » allait de 0 à +0,028 selon le plan et le calcul, et le
bilan complet de 0,03 à 0,14 d'un relevé à l'autre (ondes de pression de la
veine), d'où la référence prise dans le plan lui-même.

**Ce que disent les bilans.** Avec les trois mesures de la portance qui
s'accordent (boîte, échange aux parois, surface), le calcul porte 6 à 12 %
de plus que les références :

| Cas | Surface Cz / Cx | Échange Cz | Boîte Cz / Cx (volume nul) | Référence |
| --- | --- | --- | --- | --- |
| Aile AR 4, 4°, L = 16 + ×2 | 0,288 / 0,025 | 0,271 | 0,268 / 0,076 | VSPAERO 0,251 à 0,265 / 0,021 à 0,022 |
| Aile AR 4, 4°, L = 16 + ×4 | 0,274 / 0,033 | 0,276 | 0,272 / 0,044 | id. |
| Aile AR 4, 8°, L = 16 + ×4 | 0,545 / 0,046 | 0,569 | 0,563 / 0,059 | VSPAERO 0,501 à 0,529 / 0,039 à 0,045 |
| TR-1208, 4°, L = 8 + ×4 | 0,289 / 0,016 | 0,290 | 0,286 / 0,044 | mesure 0,257 / 0,011 ; VSPAERO 0,266 / 0,015 |
| TR-1208, 4°, L = 8 + ×8 | 0,299 / 0,018 | 0,291 | — | id. |
| TR-1208, 8°, L = 8 + ×4 | 0,551 / 0,044 | 0,550 | 0,545 / 0,049 | mesure 0,499 / 0,021 ; VSPAERO 0,531 / 0,026 |
| TR-1208, 8°, L = 8 + ×8 | 0,574 / 0,039 | 0,560 | 0,554 / — | id. |

(« — » : boîte trop grande pour la grille ×8 dans ces calculs ; `run3d`
utilise maintenant 0,1 et 0,15 longueur.) La portance dépasse VSPAERO non
visqueux : le décambrage visqueux de la couche limite (5 à 6 % ici) manque,
et un peu plus. La traînée reste mal déterminée : la surface et l'échange
donnent 0,016 à 0,046, la boîte 0,044 à 0,076, au-dessus des références dans
tous les cas. Hypothèse pour la suite : la loi de paroi, qui fait glisser
l'écoulement au premier nœud, rend une couche limite trop mince pour
décambrer le profil ; la résolution du bord d'attaque n'explique qu'une
partie de la traînée.

### D'où vient l'excès de portance

Trois pistes écartées, une confirmée.

- *Interférence de la veine : non.* Bords latéraux périodiques : la veine
  simule une rangée d'ailes, voisines à 2 envergures environ avec
  `--marge 0.4`. Une somme sur le réseau des tourbillons en fer à cheval
  images prévoit +4 % de portance sur l'aile AR 4 si rien ne les amortit.
  Calculée (AR 4, 4°, ×2) avec `--marge` 0,4, 0,8 et 1,2 (2, 13 et 39 M
  mailles), la portance du bilan vaut 0,2682, 0,2683 et 0,2693 : les bandes
  absorbantes effacent les images.
- *Écoulement amont : non.* Vitesse et masse volumique moyennes sur les
  sections de la veine à 0,2 % de l'entrée près (`freestream3d`). Le
  nombre de Mach du réseau (u = 0,1, Ma = 0,17) ajoute environ 1 % à la
  portance d'une aile d'allongement 4 (Prandtl–Glauert).
- *Résolution : non.* La portance croît encore de ×2 à ×4 (0,268, 0,272)
  et de ×4 à ×8 (TR-1208 8° : 0,545, 0,554).
- *Décambrage visqueux absent : oui.* En glissement libre, l'aile AR 4
  porte autant qu'avec la loi de paroi (0,299 et 0,301, échange, ×2) :
  la couche limite du calcul ne retire rien à la portance. Dans
  l'écoulement réel, elle lui retire 6 à 10 % (NACA 0012, Re = 6·10⁶,
  10° : Ladson 1,08 contre 1,19 non visqueux). À ces résolutions, la couche
  limite tient presque entière sous la hauteur de raccord de la loi de
  paroi (2 mailles, 3 % de la corde à ×8) : son épaisseur de déplacement,
  qui grandit vers le bord de fuite de l'extrados, n'existe pas dans le
  calcul. La loi de paroi locale ne peut pas la rendre : à hauteur fixe,
  son défaut de débit (∝ u_τ h) diminue vers le bord de fuite au lieu de
  croître. *Correction :* les profils de vitesse relevés depuis montrent
  une couche déficitaire résolue au-dessus de WM_H, de déplacement plus
  fort que celui d'une couche limite réelle (section suivante) ; cette
  explication ne tient pas.

Il faudrait un couplage visqueux–non visqueux : épaisseur de déplacement
tirée d'une couche limite intégrale le long de la surface, imposée par
transpiration (vitesse normale d(u_e δ*)/ds aux parois). La référence 2D
elle-même ne le fait pas assez : la méthode des panneaux + couche limite
(couplage par profil épaissi, refermé sur la position moyenne des deux
faces au bord de fuite) ne retire que 1 % de portance (`decamber2d` :
1,182 contre 1,193 à Re = 6·10⁶, 10°, 20 itérations ; le test garde
1,168 à une itération, valeur de l'application web). Essayé : couplage par
sources équivalentes sur le profil et le sillage (Lighthill, Lock), condition
de Kutta par égalité des vitesses au bord de fuite. Stable seulement avec la
vitesse de la couche limite figée au-delà de 95 % de corde (sinon la
décélération vers le bord de fuite fermé du contour fait diverger δ*), il
ne retire alors que 1,2 % (4°) et 1,7 % (10°) : la croissance de δ* sur les
derniers pour cent de l'extrados, qui fait l'essentiel du décambrage, est
perdue. Un couplage simultané (Newton, comme XFOIL) était nécessaire : il
est décrit ci-dessous et retire bien 8,5 % de portance à 10°.

### Couplage simultané (Newton) en 2D

`PanelOptions { coupling: Coupling::Newton, .. }` (`lbm/src/vii.rs`),
`soufflerie panel --couplage newton` ; par défaut, le couplage faible reste
celui de l'application web (`PanelOptions::free` et `tripped` inchangés).
Méthode de XFOIL (Drela & Giles 1987, Drela 1989) :

- *Écoulement potentiel* : mêmes panneaux tourbillons (Kuethe & Chow,
  80 panneaux par face), plus des sources de transpiration σ = d(ue δ*)/ds
  sur le profil et sur un sillage d'une corde, tracé le long de la ligne de
  courant non visqueuse issue du bord de fuite (bissectrice sur les deux
  premiers centièmes). La condition de glissement est écrite du côté
  intérieur de la nappe : l'intérieur reste au repos et la vitesse normale
  extérieure vaut σ. ue = ue_inv + D·m, m = ue δ*, D matrice d'influence
  (sources → tourbillons par le glissement et la condition de Kutta →
  vitesse aux stations). Sources à densité linéaire, continues aux nœuds,
  valeurs nodales par différences compactes de m entre stations voisines
  (des sources uniformes par panneau, vues aux points de contrôle, ne voient
  pas un m en dents de scie : il apparaissait dans les zones décollées
  laminaires). Vérifié sur une bosse de δ* : la transpiration rend la
  portance du corps épaissi géométriquement à 2 % près de l'écart (essai
  `transpiration_contre_corps_epaissi`).
- *Couche limite* : intégrales de quantité de mouvement et d'énergie
  cinétique, différences logarithmiques de XFOIL (BLDIF), décentrement selon
  Hk. Fermetures de Drela & Giles (laminaire : H*, Cf, dissipation) et de
  Drela (turbulent : H* , frottement de Swafford, dissipation
  Cf/2·Us + Cτ(1 − Us), équation de retard de Green pour √Cτ, K_C = 5,6,
  G–β A = 6,7, B = 0,75 ; frottement et dissipation bornés par leurs valeurs
  laminaires aux faibles Re_θ, comme XFOIL). Sillage : Cf = 0, dissipation
  doublée, constante de retard 0,9. Point d'arrêt : similitude de Hiemenz.
- *Transition* : enveloppe eᴺ de Drela & Giles (N_crit = 9) avec la rampe
  de XFOIL autour du Re_θ critique ; point de transition interpolé dans
  l'intervalle, N₂ implicite (TRCHEK2), √Cτ initial
  1,8 e^(−3,3/(Hk−1)) √Cτ_éq ; transition imposée par `x_trip`.
- *Newton* : inconnues θ, m et N ou √Cτ à chaque station du profil et du
  sillage (600 inconnues), jacobienne par différences finies locales
  (vérifiée contre une différentiation du système complet), y compris la
  sensibilité de l'abscisse curviligne au déplacement du point d'arrêt.
  Comme XFOIL, ue est une variable dont l'écart à ue_inv + D·m est résorbé
  par le pas de Newton ; pas limités (variations relatives entre −50 % et
  +150 %). Départ : marche directe à ue non visqueuse, inverse (Hk imposé)
  au décollement. Rien n'est figé au bord de fuite. Si Newton échoue, il
  repart de 0° par pas de 2° au plus. Convergence en 5 à 25 itérations,
  0,3 à 2 s par incidence.
- *Traînée* : Squire–Young au bout du sillage.

Validation (`decamber2d`, Cz/Cz_inv : rapport à la portance non visqueuse
des mêmes panneaux) :

| Cas | α | Cz | Cz/Cz_inv | Cx | Référence |
| --- | --- | --- | --- | --- | --- |
| NACA 0012, Re = 6·10⁶, déclenchée | 0° | 0 | — | 0,00826 | TMR (SA) 0,00809 |
| | 4° | 0,4462 | 0,930 | 0,00879 | |
| | 10° | 1,0924 | 0,915 | 0,01222 | TMR 1,0909 / 0,01231 ; Ladson ≈ 1,08 |
| | 15° | 1,5526 | 0,874 | 0,01958 | TMR 1,5461 / 0,02124 |
| NACA 0012, Re = 10⁶, libre | 0° | 0 | — | 0,00523 | transition 0,70 / 0,70 |
| | 2° | 0,2074 | 0,864 | 0,00564 | transition 0,49 / 0,87 |
| | 4° | 0,4251 | 0,886 | 0,00709 | transition 0,27 / 0,96 |
| | 6° | 0,6768 | 0,942 | 0,00943 | transition 0,09 / 0,99 |
| | 8° | 0,9273 | 0,969 | 0,01187 | transition 0,04 / 1,00 |
| NACA 2412, Re = 3·10⁶, libre | −2° | 0,0186 | — | 0,00569 | |
| | 0° | 0,2361 | 0,931 | 0,00531 | |
| | 2° | 0,4504 | 0,912 | 0,00493 | |
| | 4° | 0,6639 | 0,906 | 0,00548 | |
| | 8° | 1,1154 | 0,923 | 0,00959 | |
| | 10° | 1,3074 | 0,906 | 0,01155 | |

À Re = 6·10⁶, la perte de portance visée est atteinte : 8,5 % à 10° (TMR et
Ladson : 8,6 à 9,5 % sous la valeur non visqueuse 1,193), à 0,1 % de la
portance TMR ; 12,6 % à 15°. La traînée est à 2 % de TMR à 0° et 10°, 8 %
en dessous à 15°. Le résultat ne dépend ni du nombre de panneaux (rapport
0,916, 0,915 et 0,914 avec 60, 80 et 120 panneaux par face) ni de la
longueur du sillage (½, 1 et 2 cordes : 1,0855 à 1,0858 avant le passage
aux sources linéaires). Pour le NACA 2412, α₀ ≈ −2,2° et la pente
0,107 par degré entre 0 et 4° ; Cm ≈ −0,05. Abbott & von Doenhoff donnent
α₀ ≈ −2° et une pente voisine (valeurs lues sur les courbes, à revérifier
avant d'en faire un critère) ; la traînée minimale calculée (0,0049) est
plus basse que leur 0,006 environ, la transition libre calculée tombant
plus en arrière que sur la maquette.

Le couplage faible, lui, ne retire que 1 à 2 % ; il reste l'option par
défaut (valeurs de l'application web).

Limites :

- *Transition libre aux faibles Re* : à Re = 10⁶, la perte relative de
  portance diminue avec l'incidence (14 % à 2°, 3 % à 8°). Elle ne dépend pas
  du maillage (13,5 et 13,6 % à 2° avec 60 et 120 panneaux) : c'est la
  réponse du modèle à une transition très dissymétrique (extrados turbulent
  bien avant l'intrados, qui décolle en laminaire près du bord de fuite).
  Faute de polaire XFOIL de référence sous la main, elle n'est pas validée.
- *Convergence* : Newton échoue (Cz et Cx vides, `NaN`) quand la transition
  tombe sur les toutes dernières stations du bord de fuite (NACA 2412,
  Re = 3·10⁶, 6°) et au-delà de 10° environ à Re = 10⁶ et de 12° pour le
  2412 (bulbe laminaire au bord d'attaque). La continuation en incidence
  ne suffit pas toujours. Le point de transition ne se déplace que d'une
  station par itération.
- Pas de bord de fuite épais (profils refermés), pas de compressibilité
  dans la couche limite, pas de décrochage au-delà de ce que la couche
  limite intégrale décrit (la correction de Kirchhoff n'est pas appliquée
  au couplage simultané).

### Décambrage par transpiration : essayé, et ce qu'il révèle

Essai du couplage proposé ci-dessus (`--transpiration 1`, `lbm/src/ibl3d.rs`,
`Gpu3D::enable_transpiration` et `couple_transpiration`). Résultat : il ne
corrige pas l'excès de portance, et la mesure faite en chemin contredit
l'hypothèse du « déplacement manquant ».

**Méthode.** Tous les demi-temps convectifs, à partir de t* = 6 :
- vitesse extérieure u_e de chaque triangle du maillage, relevée sur la
  grille entière la plus fine à la hauteur de raccord WM_H le long de la
  normale (plus loin si un nœud voisin est solide), projetée dans le plan ;
- couche limite intégrale marchée sur le graphe des triangles voisins, en
  volumes finis décentrés amont sur les équations en lignes de courant sans
  écoulement transversal (∇·(u_e² θ ê) + δ* u_e ê·∇u_e = τ_p/ρ, et
  ∇·(u_e θ H₁ ê) = u_e C_E), ce qui compte la divergence des lignes de
  courant ; Thwaites en laminaire, Head et Ludwieg–Tillmann en turbulent
  (fermetures copiées de `panel.rs`) ; turbulente dès le triangle qui suit
  la source, sauf `--xtr` ; ordre topologique (un triangle après tous ses
  amonts, boucles de flux rompues au plus tôt atteint) ; sans amont, la même
  équation avec Thwaites donne Hiemenz au point d'arrêt (θ² ≈ 0,09 ν/a) et
  Blasius au bord d'une plaque ; H plafonné à 2,4 (décollement signalé) ;
- v_n = ∇·(u_e δ* ê) par les mêmes flux (Lighthill 1958), le déplacement
  quittant la surface par les arêtes vives convergentes (bord de fuite),
  moyenné par l'aire sur chaque triangle et ses voisins (même débit, sans
  valeur démesurée sur un triangle effilé) ;
- à chaque maille pariétale, v_n du triangle le plus proche dont la face
  la regarde, sous-relaxée (ω = 0,5), ajoutée selon la normale de la face
  à la vitesse de la paroi mobile de Bouzidi (débit injecté, efforts
  d'échange cohérents). Grille entière la plus fine et grilles du bord
  d'attaque.

Tests unitaires de la couche limite (`ibl3d`) : plaque laminaire contre
Blasius (θ à 4 %), plaque turbulente contre la loi en puissance 1/7 (à
15 %), point d'arrêt contre Hiemenz, débit injecté égal au débit qui quitte
la surface.

**Résultats** (bilan sur une boîte pour Cz, ligne `Cx =` de `run3d` pour Cx ;
sans transpiration : section précédente) :

| Cas | Sans | Avec transpiration | Référence |
| --- | --- | --- | --- |
| AR 4, 4°, L = 16 + ×2 | 0,268 | 0,235 / 0,024, stable | VSPAERO 0,251 à 0,265 ; ligne portante 0,251 |
| AR 4, 4°, L = 16 + ×4 | 0,272 | divergence à t* = 9,8 | id. |
| AR 4, 8°, L = 16 + ×4 | 0,563 | divergence à t* = 7,3 | VSPAERO 0,501 à 0,529 ; ligne portante 0,498 |
| TR-1208, 4°, L = 8 + ×4 | 0,286 | divergence à t* = 6,3 ; avec la moyenne par voisinage : 0,321 / 0,017 | mesure 0,257 / 0,011 |
| TR-1208, 8°, L = 8 + ×4 | 0,545 | divergence à t* = 6,3 (avant la moyenne par voisinage) | mesure 0,499 / 0,021 |

À ×2, l'aile AR 4 perd 12 %, un peu sous les références ; à ×4, le calcul
diverge ; sur TR-1208, la portance monte au lieu de baisser, avec des v_n
locales absurdes (jusqu'à 450 u∞ avant la moyenne d'un triangle à l'autre,
δ* de l'intrados au bord de fuite à 5 % de la corde). Les ×8 n'ont pas été
lancés.

**Pourquoi : le champ contient déjà une couche déficitaire épaisse.**
Profils de vitesse tangentielle le long de la normale (aile AR 4, 4°,
mi-envergure, exemple `ibl3d`) :

| Extrados, u_t/u∞ à h mailles | h = 1 | 2 | 3 | 4 | 6 | 10 |
| --- | --- | --- | --- | --- | --- | --- |
| ×2, x/c = 0,90 | 0,76 | 0,87 | 0,97 | 1,01 | 1,01 | 1,02 |
| ×2, x/c = 0,95 | 0,73 | 0,81 | 0,93 | 1,00 | 1,00 | 1,01 |
| ×4, x/c = 0,90 | 0,68 | 0,76 | 0,87 | 0,95 | 1,02 | 1,02 |
| ×4, x/c = 0,96 | 0,59 | 0,69 | 0,81 | 0,90 | 0,99 | 1,00 |

Au bord de fuite, la vitesse n'atteint l'extérieur qu'à 4 mailles (×2) ou
6 (×4) de la paroi, soit 9 à 12 % de la corde, quand la couche limite réelle
en fait 3 à 4 %. Hors de la couche [0, WM_H] de la loi de paroi, le calcul
résout donc lui-même une couche déficitaire (Spalart–Allmaras sur ces
grilles, mailles d'arête adhérentes au bord de fuite). Son épaisseur de
déplacement ∫(1 − u_t/u_e) dh, jusqu'au premier maximum de u_t au-dessus de
WM_H (vitesse à la paroi extrapolée linéairement depuis 1 et 2 mailles),
comparée à celle de la couche limite intégrale marchée avec u_e pris à ce
maximum (aile AR 4, 4°, ×2, mi-envergure, derniers 5 % de corde) :

| δ*/c au bord de fuite | extrados | intrados | différence |
| --- | --- | --- | --- |
| couche limite intégrale (u_e au bord de la couche résolue) | 0,005 | 0,0035 | 0,0015 |
| déjà dans le champ du calcul | 0,022 à 0,026 | 0,012 à 0,014 | 0,010 |

Le calcul décambre donc déjà, par le déplacement, plus qu'une couche limite
réelle (quatre fois l'épaisseur sur l'extrados, six fois la dissymétrie).
L'hypothèse de la section précédente (couche limite sous la hauteur de
raccord, donc sans déplacement) ne tient pas au bord de fuite. Deux
conséquences :
- Le couplage essayé (u_e à WM_H, δ* entier) relève u_e au fond de cette
  couche : 0,80 u∞ au bord de fuite de l'extrados à ×2 au lieu d'environ
  1,0 au-dessus, 0,69 à ×4. Le gradient adverse, trop fort, épaissit la
  couche limite intégrale (δ*/c = 0,013 au lieu de 0,005), d'autant plus
  que la grille est fine (décollement sur 4 % de la surface à ×4, 11 % à 8°)
  ; le déplacement s'ajoute en outre à celui du champ. La baisse de 12 % à
  ×2 vient de ce double compte, pas d'un décambrage physique.
- Le couplage cohérent (vitesse extérieure équivalente : u_e au bord de la
  couche résolue, et seul le déplacement qui manque au champ,
  δ*_intégrale − δ*_calcul) donne une transpiration surtout négative
  (aspiration), qui augmenterait la portance. Essayé dans l'exemple `ibl3d`
  (mode 2) : le bord de la couche résolue, premier maximum du profil, saute
  d'un triangle à l'autre, et sa dérivée donne |v_n| jusqu'à 3 u∞ ; la
  portance oscille entre 0,19 et 0,24. Non retenu dans `run3d`.

**Ce qui reste.** L'excès de portance de 6 à 12 % ne vient pas d'un
déplacement manquant : le champ en a déjà trop. Pistes : ce qui fixe la
circulation au bord de fuite (mailles d'arête adhérentes, glissement de la
loi de paroi juste en amont, épaisseur du bord de fuite voxelisé), et
pourquoi une couche déficitaire aussi épaisse décambre si peu (en
glissement libre, la portance est la même : vérifier si la couche existe
aussi sans loi de paroi, auquel cas elle est numérique). La couche limite
intégrale et sa mesure restent utiles pour ces diagnostics.
`--transpiration 1` reste disponible pour reproduire l'essai, mais ne doit
pas servir aux calculs.

### Bord de fuite : une erreur de placement d'ordre un

L'excès de portance existe déjà en glissement libre (aile AR 4, 4°,
×2 : 0,299 contre 0,265 pour VSPAERO en panneaux épais) : il n'est pas
visqueux. Il dépend de la place du bord de fuite dans la maille.
`run3d --decalage d` décale la maquette de d mailles de la veine ; la
fraction f du bord de fuite dans la maille de la grille fine en découle
(`te_snap` la cale à 0,7 sur la veine, soit 0,4 à ×2, 0,8 à ×4, 0,6 à ×8).
Portance du bilan sur boîte, aile AR 4, 4°, L = 16 :

| f (×2) | 0,15 | 0,28 | 0,40 | 0,53 | 0,65 | 0,78 | 0,90 | 0,03 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| loi de paroi, ×2 | 0,252 | 0,269 | 0,268 | 0,277 | 0,283 | 0,291 | 0,294 | 0,279 |

| f (×4) | 0,80 | 0,05 | 0,30 | 0,55 |
| --- | --- | --- | --- | --- |
| loi de paroi, ×4 | 0,273 | 0,275 | 0,275 | 0,260 |

En glissement libre (×2), même dent de scie : 0,287 et 0,299 pour f = 0,4,
0,322 et 0,323 pour f = 0,9. La portance croît pendant que le bord de fuite
avance vers le nœud suivant, puis retombe quand il le passe : la condition
de Kutta discrète reste attachée aux nœuds pendant que la géométrie se
déplace. L'écart (±8 % à ×2, ±3 % à ×4) est d'ordre un en maille ; la
moyenne sur le placement passe de 0,276 à 0,271, et son extrapolation
(Richardson, ordre un) donne 0,266, la valeur de VSPAERO en panneaux
épais. L'« excès de portance » des sections précédentes est donc pour
l'essentiel cette erreur de discrétisation au bord de fuite, plus ou moins
forte selon où tombe le bord de fuite (le calage de `te_snap` place la
grille ×4 près du haut de la dent). La traînée du bilan sur boîte varie
plus encore avec le placement (0,032 à 0,070 à ×4). Le frottement retire
7 à 9 % de portance à placement égal (0,287 → 0,268 ; 0,322 → 0,294) : le
décambrage visqueux existe bien.

### Grilles des bords (`--ba 1`, `--bf 1`) : essai

Pour réduire cette erreur sans tout raffiner, les grilles partielles
suivent maintenant les bords : plusieurs boîtes sœurs au même niveau
(`Refinement::parent`, `depth`, `core`), qui suivent en envergure le bord
d'attaque (faces tournées vers l'amont) et le bord de fuite (arêtes vives
saillantes orientées en envergure, `geom3d::feature_edges`) ; sur
TR-1208, 9 boîtes par bord au lieu d'une boîte qui couvrait presque toute
l'aile. Chaque boîte ne compte, pour les efforts, que dans son domaine
propre (pas de double compte entre sœurs). Une boîte qui entoure toute
l'aile (diagnostic `SOUFFLERIE_BOITE_ENTIERE=1`) redonne le calcul ×4
(0,2733 contre 0,2734) : le mécanisme fonctionne.

Mais une boîte qui coupe la maquette fausse l'écoulement : aile AR 4, 4°,
×2, portance 0,2685 sans boîte, 0,189 avec la boîte du bord d'attaque,
0,174 avec celle du bord de fuite ; avant cette refonte déjà, la boîte
unique du bord d'attaque donnait 0,318. Le Cp à mi-envergure montre une
succion plus faible dans toute la boîte, jusqu'à 6 mailles de la paroi, et
un saut là où la face de la boîte coupe l'extrados : la circulation
change. Deux défauts du couplage près d'une paroi ont été corrigés (sans
effet sur les grilles qui entourent toute la maquette) :
- remplissage des mailles fantômes : les populations tirées d'un nœud
  solide de la grille parente sont remplacées par le rebond du nœud
  lui-même, et l'interpolation (linéaire près d'une paroi) ne lit plus
  les nœuds solides ;
- recopie vers la grille parente : une maille pariétale de la grille fine
  rend ses moments avec leur partie hors équilibre (la contrainte), et non
  un équilibre.

Le défaut reste (0,189 avec ces corrections). Cause probable : le couplage
par les moments n'est pas conservatif là où une paroi traverse la face
d'une grille (couche limite et loi de paroi à deux résolutions de part et
d'autre). Il faudrait un raccord conservatif (volumique, à la Rohde) ou
des boîtes qui ne coupent pas la maquette. Jusque-là, `--ba` et `--bf`
affichent un avertissement et ne doivent pas servir aux calculs ;
l'incertitude de placement (±3 % à ×4 sur la portance) est à garder en
tête, ou à moyenner en calculant plusieurs `--decalage`.

### Couche pariétale fine (`--bande D`)

Plutôt que des boîtes, une grille deux fois plus fine sur une couche
d'épaisseur D (mailles fines) autour de toute la surface, comme les
régions de résolution décalées de la paroi des codes Boltzmann
commerciaux. Les faces de raccord avec la grille parente sont parallèles à
la paroi, dans le fluide : elles ne coupent jamais la maquette (le défaut
des grilles des bords). Première étape : la grille fine reste une boîte
qui entoure la maquette, mais les mailles fluides à plus de D de la paroi
sont des mailles fantômes, remplies depuis la grille parente et non
calculées. Seules les fantômes voisines d'une maille active sont remplies
(`band_fill`), et seules les mailles actives à plus de 2 MARGIN + 1 de
toute fantôme sont recopiées vers la grille parente (`band_deep` ; la
recopie suit désormais une liste de nœuds, boîte ou couche).

Aile AR 4, 4°, L = 16, grille ×2 et couche ×4, contre la grille ×4
complète (même placement du bord de fuite) :

| D (mailles ×4) | 8 | 12 | 16 | ×4 complet |
| --- | --- | --- | --- | --- |
| Cz (boîte) | 0,257 | 0,2775 | 0,2732 | 0,2734 |
| Cx (boîte, volume nul) | 0,062 | 0,0421 | 0,0422 | 0,0421 |
| Cz (échange) | 0,292 | 0,280 | 0,279 | — |
| durée | 1,0 min | 1,2 min | 1,4 min | ≈ 10 min |

À 8 mailles, la couche est trop mince : le raccord tombe dans la couche
limite et la zone de la loi de paroi, et l'échange et la boîte se
contredisent. À 16 mailles (un quart de corde à ×4), la couche redonne le
×4 complet, sept fois plus vite. TR-1208, 8°, grille ×4 et couche ×8
(D = 16) : Cz 0,558 (boîte), 0,564 (échange), Cx 0,0388 (surface), contre
0,5545, 0,561 et 0,0393 pour le ×8 complet ; 4,6 min au lieu d'environ
30.

Erreur de placement du bord de fuite avec la couche ×8 (aile AR 4, grille
×4 et couche D = 16) :

| f (×8) | 0,60 | 0,85 | 0,10 | 0,35 | moyenne |
| --- | --- | --- | --- | --- | --- |
| Cz (boîte) | 0,2685 | 0,2767 | 0,2677 | 0,2626 | 0,269 |

La moyenne sur le placement converge (0,276 à ×2, 0,271 à ×4, 0,269 à ×8,
vers environ 0,267 ; VSPAERO en panneaux épais : 0,265), et l'écart diminue
(±8 %, ±3 %, ±2,6 % sur quatre placements seulement).

**Stockage creux par blocs.** La couche est rangée par blocs de 8³
mailles : seuls les blocs qui contiennent une maille active ou une
fantôme à remplir sont gardés (`SparseMap`), plus un bloc vide de mailles
fantômes qui reçoit toute lecture hors des blocs. Les shaders trouvent le
rang d'une maille par une table des blocs (`id_of`) et la position d'un
rang par les coordonnées du bloc (`pos_of`) ; le raccord des grilles lit la
position des fantômes dans une liste et la grille fine par la même table ;
les relevés sur le processeur (surface, échange, sillage) passent par la
même correspondance. Mêmes résultats que la boîte : aile AR 4, Cz 0,2733
(boîte : 0,2732), 1,34 M mailles gardées au lieu de 2,88 M ; TR-1208, 8°,
couche ×8, Cz 0,5582 (0,5583), 2,6 M mailles au lieu de 16 M, 3,6 min.

La couche ×16 devient possible (8,3 M mailles ; la boîte en ferait plus
de 50 M) :

| TR-1208, 8° | Cz boîte | Cz échange | Cx surface |
| --- | --- | --- | --- |
| ×4 | 0,545 | 0,550 | 0,044 |
| ×8 complet | 0,5545 | 0,561 | 0,039 |
| ×4 + couche ×8 | 0,558 | 0,564 | 0,039 |
| ×8 + couche ×16 (26 min) | 0,562 | 0,568 | 0,0375 |
| mesure / VSPAERO | 0,499 / 0,531 | | 0,021 / 0,026 |

La traînée continue de baisser en affinant ; la portance converge vers
environ 0,565, 13 % au-dessus de la mesure et 6 % au-dessus de VSPAERO
non visqueux. Sur TR-1208, l'excès de portance ne se résorbe donc pas
avec la résolution près de la paroi, contrairement à l'aile AR 4 où la
moyenne sur le placement rejoint VSPAERO.

**Couches imbriquées (`--bande 16,16,…`, `--bande auto`).** Chaque couche
est deux fois plus fine que la précédente et n'existe qu'à moins de son
épaisseur (en ses propres mailles) de la paroi ; une couche est la grille
parente de la suivante (remplissage et recopie lisent la couche parente
par sa table de blocs). La mise en place est elle aussi creuse : les blocs
sont choisis par la distance de leur centre au maillage, puis classés
maille par maille (solide, active à moins de D de la paroi, fantôme), et
seuls sont gardés ceux qui contiennent une maille active ou une fantôme à
remplir ; mémoire et temps de mise en place suivent la surface. Les
efforts d'une couche comptent hors des couches plus fines (échange aux
parois sans la part de la pression de référence, qui ne s'annule que sur
une surface fermée). `--bande auto` choisit le nombre de couches pour
mettre au moins 2 mailles dans le plus petit rayon de bord d'attaque des
coupes (au plus 4 couches de 16 mailles) et prévient si le bord d'attaque
est à facettes (coin de plus de 60° : le rayon mesuré est celui des
facettes, `--lisser 1` le rend). Au-delà de 40 M mailles en tout, la mise
en place s'arrête.

Aile AR 4, 4°, grille ×2 :

| couches | mailles des couches | Cz boîte | Cz échange | Cz surface | Cx boîte | durée |
| --- | --- | --- | --- | --- | --- | --- |
| aucune | — | 0,2684 | 0,2715 | — | 0,073 | 0,4 min |
| ×4 (16) | 1,4 M | 0,2714 | 0,2783 | — | 0,043 | 1,8 min |
| ×4, ×8 (16, 16) | 7,0 M | 0,2694 | 0,2692 | — | 0,046 | 10,5 min |
| ×4, ×8, ×16 (`auto`) | 21,2 M | 0,2797 | 0,2820 | 0,2712 | 0,033 | 55 min |

Avec deux couches, le calcul rejoint la grille ×4 complète plus une couche
×8 (0,2685). La troisième couche n'est pas encore sûre : les trois mesures
de la portance s'écartent (0,271 à 0,282) et la tendance s'inverse. Seize
mailles ×16 ne font que 6 % de la corde, moins que la couche de 8 mailles
×4 (12 %) déjà trop mince : l'épaisseur nécessaire est sans doute physique
(couche limite, zone de la loi de paroi, écoulement du bord de fuite), pas
un nombre de mailles. À vérifier avant de s'y fier.

`--bande-corde f` impose à chaque couche une épaisseur d'au moins f fois
la longueur de référence (en ses mailles). Premier essai, f = 0,125 :
épaisseurs 16, 16 et 32 mailles, 34 M mailles, 100 min ; écoulement faux
dès t* = 5 (Cz 0,03 à 0,06). La couche ×16 de 32 mailles atteignait toute
l'épaisseur de sa couche parente (16 mailles ×8) : ses fantômes étaient
remplies depuis des mailles que la couche parente ne calcule pas.
`nested_bands` impose désormais D_parent ≥ D/2 + 2 MARGIN + 3 (16, 16, 32
devient 19, 23, 32) ; l'essai à trois couches de 16 mailles respectait
déjà cette règle, son incohérence a une autre cause.

Essai refait avec la règle (19, 23 et 32 mailles, 36 M mailles, 86 min) :

| Aile AR 4, 4°, ×2 + couches ×4, ×8, ×16 | Cz boîte | Cz échange | Cz surface | Cx boîte | Cx surface |
| --- | --- | --- | --- | --- | --- |
| 16, 16, 16 mailles (×16 : 6 % de corde) | 0,2797 | 0,2820 | 0,2712 | 0,033 | 0,030 |
| 19, 23, 32 mailles (×16 : 12,5 % de corde) | 0,2724 | 0,2809 | 0,2702 | 0,023 | 0,030 |

Avec la couche la plus fine épaisse de 12,5 % de corde, la boîte et la
surface s'accordent à 0,8 % près (3 % avec 16 mailles) ; l'échange reste
3 % au-dessus. La portance (0,270 à 0,272, un seul placement du bord de
fuite) prolonge la suite 0,276, 0,271, 0,269 des moyennes à ×2, ×4, ×8 ;
VSPAERO en panneaux épais donne 0,265, la ligne portante à polaires
visqueuses 0,251. La traînée de surface descend à 0,030 (0,042 à ×4 ;
VSPAERO 0,021 à 0,022). Règle retenue : `--bande auto --bande-corde 0.125`.

Le rangement des triangles pour la distance au maillage utilisait des
cases d'une maille quelle que soit la grille : sur une couche ×16, chaque
recherche parcourait des dizaines de milliers de cases. Les cases ont
maintenant la longueur médiane des arêtes (au moins une maille).

### Normales lisses et relevés multiples sur la surface

Le maillage d'une maquette ne décrit qu'une surface plane par facettes ;
deux corrections la rapprochent de la surface réelle sans rien changer à
la grille :
- **normales lisses** : chaque coin de triangle reçoit la moyenne des
  normales des faces voisines de son groupe lisse, sans traverser d'arête
  vive (`geom3d::corner_normals`, seuil de 60°), et la normale en un point
  du triangle est l'interpolation de ses coins (Phong). Elle sert au point
  où un lien coupe la paroi (normale de la face touchée : plan de la
  vitesse de glissement, critère d'arête vive des mailles) et aux relevés
  de la surface. Sur une sphère facettée, l'écart à la normale radiale
  passe de 1 − cos = 1,5·10⁻² à zéro (test `normale_lisse_sphere`) ;
- **relevés multiples** : un triangle plus grand qu'une maille de la grille
  où il est relevé est découpé en m² sous-triangles (m ≤ 8), chacun relevé
  en son centre le long de sa normale lisse, au lieu d'un seul relevé au
  centre le long de la normale plate.

Le bilan sur boîte, indépendant de la surface, sert d'arbitre :

| Cas | Cz surface / boîte avant | après | Cx surface avant → après |
| --- | --- | --- | --- |
| TR-1208, 8°, L = 8 + ×4 | 0,5510 / 0,5446 | 0,5394 / 0,5394 | 0,0441 → 0,0387 |
| TR-1208, 8°, ×4, STL lissé (`--lisser 1`) | — | 0,5556 / 0,5550 | 0,0431 |
| Aile AR 4, 4°, L = 16 + ×2 | 0,2897 / 0,2685 | 0,2859 / 0,2675 | 0,0244 → 0,0164 |

Sur TR-1208, la surface et la boîte coïncident désormais (écart 0 à 0,1 %,
contre 1,2 %) ; l'échange aux parois de l'aile AR 4 rejoint lui aussi la
boîte (0,2679). Les normales lisses des liens changent un peu l'écoulement
(portance de la boîte −1 % sur TR-1208) : la vitesse de glissement suit la
surface réelle et non la facette. Sur l'aile AR 4 à ×2 (32 mailles par
corde), la surface reste 7 % au-dessus de la boîte : le bord de fuite y
tient dans une à deux mailles. Le STL lissé de TR-1208 porte 3 % de plus
(rayon du bord d'attaque rendu).

### Programmes de diagnostic

Les essais décrits dans ce document se refont avec les programmes de
`lbm/examples` (`cargo run --release -p lbm --example <nom> -- <arguments>`,
arguments en tête de chaque fichier) :

| Programme | Ce qu'il mesure |
| --- | --- |
| `cp3d` | Cp pariétal d'une tranche d'envergure, traînée par tranche de corde, profils de pression le long de la normale ; dernier argument 0 : simple précision |
| `edge3d` | mailles d'arête vive (loi de paroi coupée) par tranche de corde |
| `pwh3d` | portance et traînée de surface selon les hauteurs de relevé de la pression, contre le bilan sur une boîte |
| `box3d` | bilan sur des boîtes de tailles croissantes, traînée extrapolée à volume nul |
| `wake3d` | traînée en champ lointain dans cinq plans aval, relevé par relevé |
| `bilan3d` | bilan sur une boîte dans le champ de la veine (aile analytique) |
| `surface3d`, `trainee3d` | efforts intégrés sur la surface ; gradient de pression le long de la veine |
| `blow`, `probe` | où et quand un calcul 3D (STL) ou raffiné diverge |
| `te_map`, `voxel_check`, `wing_vox`, `wing_offset`, `link_sum` | voxelisation : bord de fuite, volume conservé, décalage d'une demi-maille, poussée de la pression de référence |
| `refine3d_check`, `sa3d_check`, `transition3d`, `noise3d`, `time_mesh` | raffinement contre grille uniforme, Spalart–Allmaras 3D, transition imposée, bruit, coût de mise en place |
| `cp2d`, `split` | 2D : Cp et Cf pariétaux, part frottement / pression, contre les panneaux |
| `freestream3d` | vitesse et masse volumique réellement vues par la maquette, le long de la veine et par section |
| `ibl3d` | couche limite intégrale sur le maillage (θ, H, δ*, v_n par triangle, CSV), profils de vitesse le long de la normale, couplages de transpiration (mode 1 : u_e à WM_H ; mode 2 : déplacement manquant au champ) |
| `decamber2d` | décambrage visqueux 2D : couplage faible (1 et 20 itérations) contre couplage simultané (Newton), traînée, transition ; `decamber2d 0012 6e6 0 0:15:1`, `decamber2d 2412 3e6 libre -2:12:2` |
| `geom3d_check`, `geom_edges` | géométrie d'une aile maillée : rayon du bord d'attaque, facettes, parois minces, mailles d'arête vive loin des arêtes du maillage ; angles dièdres |
| `edge_boxes`, `te_map` | boîtes des grilles des bords (étendue, domaine propre) ; mailles autour du bord de fuite selon `--decalage` |
| `plaque3d` | couche limite turbulente d'une plaque plane d'envergure périodique contre la théorie (δ*, θ, H, Cf, profil de Coles), perte de pression totale le long de x, variantes de diagnostic de la loi de paroi (`Gpu3D::set_essai`) |
| `banc3d` | banc de mesure de la vitesse 3D (voir « Vitesse de calcul 3D ») : durées de mise en place, de calcul et de relevé, profil par passe, résultats à temps fixé comparés à une référence |

L'essai d'aile d'envergure infinie (section « Aile infinie ») reposait sur
une option de maquette retirée depuis ; son programme n'est pas conservé.

## Couche limite de la plaque plane : la couche lente est numérique

Le couplage par transpiration avait montré que le champ des ailes porte au
bord de fuite une couche lente de 4 à 6 mailles, d'épaisseur de déplacement
environ 4 fois celle d'une vraie couche limite. Cas canonique pour en
trouver l'origine : `plaque3d`, plaque d'épaisseur ≈ L/16 à bords d'attaque
et de fuite demi-elliptiques 4:1, alignée sur l'écoulement, d'envergure
périodique (elle traverse les bords y de la veine, où `Gpu3D::periodic_span`
retire les bandes absorbantes latérales) : écoulement plan, couche limite
turbulente dès le bord d'attaque. Profils u(y) de l'extrados moyennés sur
l'envergure, rapportés à l'écoulement extérieur (droite ajustée entre 0,1 et
0,2 L et prolongée jusqu'à la paroi : la plaque épaisse a sa propre
survitesse potentielle, qui décroît en s'éloignant de la paroi) ; sous le
premier nœud, profil de Spalding du u_τ de la loi de paroi. Épaisseur ajustée
pour que le premier nœud soit à la même distance (« calage », 0,3 maille
ici) des deux faces : sinon les deux faces n'ont pas la même couche lente et
la plaque symétrique porte (±10 % de vitesse extérieure).

**Constat** (borne de glissement d'origine, voir plus bas ; x = 0,5 L,
Re_L = 10⁶ sauf mention ; δ* et θ en mailles de la
veine ; théorie : δ* = 0,046 x Re_x^−0,2, θ = 0,036 x Re_x^−0,2,
Cf = 0,026 Re_x^−1/7, H ≈ 1,3) :

| Cas | δ* | θ | H | Cf (loi de paroi) | δ* / théorie |
| --- | --- | --- | --- | --- | --- |
| théorie, L = 64 | 0,107 | 0,083 | 1,29 | 0,00399 | 1 |
| L = 32 | 0,319 | 0,254 | 1,26 | 0,00307 | 6,0 |
| L = 64 | 0,271 | 0,212 | 1,28 | 0,00359 | 2,5 |
| L = 128 | 0,381 | 0,297 | 1,28 | 0,00392 | 1,8 |
| L = 32 + ×2 | 0,174 | 0,136 | 1,28 | 0,00347 | 3,3 |
| L = 64, Re = 10⁷ | 0,226 | 0,190 | 1,19 | 0,00242 | 3,4 |
| L = 64, calage 0,7 | 0,366 | 0,286 | 1,28 | 0,00344 | 3,4 |

(théorie en mailles : × L/64.) La couche lente est numérique :
- elle couvre toujours 3 à 4 mailles de la grille la plus fine (u/U_e =
  0,87, 0,92, 0,97, 1,00 aux nœuds 0,3, 1,3, 2,3, 3,3 à L = 64), quand le
  profil de Coles atteint U_e à 0,65 maille ; δ* reste de 0,2 à 0,4 maille
  quels que soient L et Re : δ*/L ne converge qu'en 1/L ;
- elle dépend de la place de la paroi dans la maille (calage 0,3 ou 0,7 :
  δ* de 0,27 à 0,37) ;
- la forme du profil est juste (H ≈ 1,27), et le frottement aussi : Cf de la
  loi de paroi à 10 % de la théorie, frottement transmis par la paroi
  (échange de quantité de mouvement) égal à ρu_τ² sur la partie plane, et θ
  croît vers l'aval au rythme Cf/2. L'excès de θ (0,13 à 0,2 maille) est
  acquis là où la loi de paroi commence, au bord d'attaque.

**Où naît l'excès.** La perte de pression totale de la colonne au-dessus de
la plaque (∫ (1 − p₀/p₀∞) dz) passe de 0,03 à 0,26 maille sur le nez, entre
2 et 4 mailles du bord d'attaque, puis ne croît plus que du frottement ; en
théorie, θ vaut 0,02 maille à cet endroit. Les efforts des mailles pariétales
du nez sont pourtant ceux du glissement libre : la perte est une dissipation
dans le fluide, pas une traînée de paroi. Le nez en glissement libre et la
loi de paroi seulement sur la partie plane : la perte apparaît au raccord, là
où commence la loi de paroi.

**Causes écartées**, un facteur à la fois (`set_essai`, L = 64, θ à
mi-plaque 0,21 à 0,22 dans tous les cas sauf mention) :
- viscosité turbulente de Spalart–Allmaras hors du premier nœud retirée :
  couche plus mince (2,5 mailles) mais même θ ;
- sous-maille de Smagorinsky : amortie près de la paroi (Van Driest), −5 à
  −10 % ; constante divisée par 2 ou par 4, sans effet ; sans sous-maille,
  le calcul diverge ;
- viscosité du premier nœud (κ u_τ h ou celle de Spalart–Allmaras), hauteur
  de raccord h = 1,5 ou 3 mailles, relaxation du glissement, demi-précision :
  sans effet ;
- forme du nez : plus émoussé, pire (nez demi-circulaire : θ = 0,30) ; plus
  effilé (10:1), sans gain ;
- viscosité du premier nœud réduite (κ u_τ y_n, ou plafond de Clauser
  0,0168 U δ*) : la paroi transmet alors 2,5 à 4,5 fois ρu_τ² ; avec le
  rebond, le premier nœud doit garder une viscosité de l'ordre de κ u_τ h ;
- glissement libre (témoin non visqueux) : il reste un défaut de 2 à 3 % sur
  3 mailles (δ* ≈ 0,05 à 0,1 maille), perte numérique au bord d'attaque.

**Cause retenue : la borne de monotonie du glissement.** Elle imposait
u_s ≤ u_h − pente · h : le profil linéaire de la couche pariétale, prolongé
jusqu'à h, ne devait pas dépasser la vitesse relevée à h. C'est supposer que
la couche pariétale atteint h. Quand la couche limite est plus mince que h
(au bord d'attaque, et partout aux résolutions courantes : 0,65 maille ici
contre h = 2), u_h ≈ u_n et la borne retranche au glissement pente · (h − y_n)
≈ 2 u_τ, soit 9 % de U ; sur une paroi convexe, où la vitesse potentielle
décroît en s'éloignant de la paroi (u_n > u_h), davantage. Le premier nœud
est freiné sans frottement correspondant ; la quantité de mouvement perdue
reste dans la couche lente. La borne porte désormais sur le profil à la
hauteur du nœud : u_s ≤ min(u_n, u_h) − pente · y_n. Elle arrête toujours un
nœud qui s'emballe (u_n ≫ u_h), raison pour laquelle elle avait été
introduite.

| x = 0,5 L | Borne à h : δ* / θ | Borne au nœud : δ* / θ | Théorie |
| --- | --- | --- | --- |
| L = 32 | 0,319 / 0,254 | 0,314 / 0,250 | 0,053 / 0,042 |
| L = 64 | 0,271 / 0,212 | 0,245 / 0,192 | 0,107 / 0,083 |
| L = 128 | 0,381 / 0,297 | 0,304 / 0,232 | 0,213 / 0,167 |
| L = 64, Re = 10⁷ | 0,226 / 0,190 | 0,182 / 0,151 | 0,067 / 0,053 |

L'écart à la théorie passe de 1,8 à 1,4 fois à L = 128 ; la perte de
pression totale du nez baisse d'un tiers.

**Ce qui reste.** Le premier nœud doit garder la viscosité κ u_τ h pour que
le rebond transmette ρu_τ², et cette viscosité fait passer le frottement aux
nœuds 2 à 4, au-dessus de la vraie couche limite : la couche lente garde
3 mailles. La corriger demanderait d'imposer la contrainte pariétale
directement (populations reconstruites au premier nœud) plutôt que par un
glissement et une viscosité.

**Damier au calage ½.** Paroi à mi-chemin entre deux rangées de nœuds
(q = 0,5 à ±0,05 près) : un mode pair–impair de période 2 mailles s'installe
au premier nœud (δ* de 0,8 à 1,3 maille d'une maille à l'autre), en
glissement libre comme avec la loi de paroi. C'est le rebond à mi-chemin, où
l'interpolation de Bouzidi n'amortit plus rien, à τ ≈ 0,5. Sur une paroi
courbe, seules quelques mailles ont q ≈ 0,5 ; non corrigé.

**Effet sur les ailes** (même code, borne à h contre borne au nœud ; bilan
sur une boîte, traînée extrapolée à volume nul) :

| Cas | Borne à h : Cz / Cx | Borne au nœud : Cz / Cx | Référence |
| --- | --- | --- | --- |
| Aile AR 4, 4°, L = 16 + ×4 | 0,272 / 0,044 | 0,274 / 0,042 | VSPAERO 0,251 à 0,265 / 0,021 à 0,022 |
| TR-1208, 4°, L = 8 + ×4 | 0,287 / 0,044 | 0,288 / 0,040 | mesure 0,257 / 0,011 |
| TR-1208, 8°, L = 8 + ×4 | 0,544 / 0,049 | 0,548 / 0,045 | mesure 0,499 / 0,021 |

La traînée du bilan baisse de 5 à 9 % (échange aux parois : 2 à 9 %), la
portance ne bouge pas (+0,5 %), et TR-1208 à 8° reste stable. L'excès de
portance ne vient donc pas de la couche lente : en l'amincissant, on ne le
réduit pas. La part de la traînée qu'elle porte (perte de quantité de
mouvement au bord d'attaque) n'en explique qu'une petite fraction.

## Maquettes STL (3D)

Dans l'interface, onglet Soufflerie 3D : « Charger un STL… » (ou glisser le
fichier sur la fenêtre, ou `soufflerie --stl fichier.stl` au lancement). Les
mêmes réglages qu'en ligne de commande y figurent (axes, unité, incidence,
dérapage, références), avec les dimensions, le maître-couple et le blocage ;
« Échelle physique » fixe le Reynolds par la vitesse et donne les efforts
en newtons.
« Couches pariétales » (avec au moins une grille ×2 « Près de la
maquette ») ajoute 1 à 3 couches de 16 mailles, chacune deux fois plus
fine, ou « auto » comme `--bande auto` : les mêmes couches creuses et les
mêmes optimisations qu'en ligne de commande ; la ligne d'information donne
le nombre de couches, leurs mailles et la vitesse sur toutes les grilles.

En ligne de commande :
`soufflerie run3d --stl fichier.stl --l 32 --re 5000` étudie une maquette
maillée (STL binaire ou ASCII, volume fermé, ou fichier `.sgeo` qui porte
aussi sa mise en place). `soufflerie geom --stl fichier.stl` en fait le
contrôle géométrique. Mise en place :

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
| `--lisser 1` | facettes raffinées sur une surface lisse, arêtes vives gardées (voir « Prétraitement géométrique ») | 0 |
| `--aretes maillage\|mailles` | loi de paroi coupée aux arêtes vives du maillage, ou par le critère par maille | mailles |

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

## Prétraitement géométrique des maquettes STL

Examen des défauts géométriques de la soufflerie 3D (escalier, place du bord
de fuite, parois minces, rayon du bord d'attaque sous la maille, facettes
grossières, arêtes vives détectées maille par maille, boîte du bord
d'attaque alignée sur les axes) : ce que l'on mesure, ce qu'un prétraitement
y change. Module `lbm/src/geom3d.rs`, commande `soufflerie geom`.

**Contrôle de résolution (`soufflerie geom`).** Pour une aile (envergure
selon y), coupe en tranches d'envergure dans la grille la plus fine : corde,
épaisseur, rayon du bord d'attaque (ajusté sur y² = 2Rs + c s^1,5 entre 0,2
et 2 % de corde, sur le polygone des facettes, c'est-à-dire ce que voit le
solveur ; NACA 0012 finement facetté : 1,588 pour 1,587), longueur et plus
grand angle des facettes du bord d'attaque, épaisseur du culot, longueur où
l'épaisseur passe sous 1 et 2 mailles ; arêtes vives du maillage ; coût du
raffinement du bord d'attaque. Sur les deux maquettes de validation, grille
×4 :

| | Aile AR 4 (L = 16, ×4) | TR-1208 (L = 8, ×4) |
| --- | --- | --- |
| Triangles | 3 984 | 2 984 |
| Corde (mailles) | 64 | 22 à 40 |
| Rayon du bord d'attaque exact | 1,02 maille | ≈ 0,2 à 0,4 |
| Rayon du polygone des facettes | 0,69 | 0,03 à 0,05 |
| Facettes du bord d'attaque | 0,4 à 0,5 maille, 48° entre elles | 0,2 à 0,3 maille, 59° (78° en coupe selon l'écoulement) |
| Épaisseur < 1 maille / < 2 mailles | 3,5 / 7,5 dernières mailles | 4 / 8 dernières mailles |
| Culot du bord de fuite | fermé | 0,01 à 0,02 maille |
| Arêtes vives (> 60°) | bord de fuite (165°), saumons (90°) | culot du bord de fuite (81°), saumons (67° à 95°) |

Les maillages exportés d'OpenVSP décrivent le bord d'attaque par des
facettes de la taille de son rayon : le polygone inscrit tronque le rayon de
32 % sur l'aile AR 4 et en fait presque un coin sur TR-1208. La paroi mince
du bord de fuite occupe toujours 3,5 à 4 mailles (angle du dièdre de fuite),
une part de corde qui diminue de moitié à chaque raffinement : rien à
prétraiter. Le culot de TR-1208 (0,02 maille) et la place du bord de fuite
(déjà calée par `te_snap`, et sans objet pour un bord de fuite en flèche)
non plus.

**Raffinement lisse des facettes (`--lisser 1`).** Chaque face est découpée
en n² sous-faces placées sur la surface de Phong qui passe par les sommets du
fichier et leurs normales (moyennes des faces voisines, pondérées par
l'angle, sans traverser d'arête vive) ; facteur ½, qui rend un arc de cercle
entre deux sommets à normales radiales (le ¾ de Boubekeur et Alexa,
choisi pour l'aspect, gonfle les arêtes de 50 %). Les arêtes vives restent
droites, le maillage reste fermé (points d'arête partagés), n est choisi
pour que les sous-facettes fassent moins de 10° entre elles (5 sur l'aile
AR 4, 6 sur TR-1208, 100 000 triangles). Icosaèdre subdivisé (80 faces) :
volume à −2,4 % de la sphère au lieu de −12,7 %. Rayon du bord d'attaque vu
par le solveur : 0,97 maille au lieu de 0,69 sur l'aile AR 4 (exact 1,02),
0,16 à 0,25 au lieu de 0,04 sur TR-1208 ; angle entre facettes 19° et 36°.
Effet sur les efforts (grille ×4) :

| Cas | Boîte Cz : brut / lissé | Surface Cz / Cx : brut | lissé | Référence |
| --- | --- | --- | --- | --- |
| Aile AR 4, 4° | 0,2720 / 0,2853 | 0,2735 / 0,0330 | 0,3012 / 0,0218 | VSPAERO 0,251 à 0,265 / 0,021 à 0,022 |
| Aile AR 4, 8° | 0,5630 / 0,5619 | 0,5452 / 0,0455 | 0,5658 / 0,0517 | VSPAERO 0,501 à 0,529 / 0,039 à 0,045 |
| TR-1208, 4° | 0,2862 / 0,2889 | 0,2888 / 0,0155 | 0,2877 / 0,0149 | mesure 0,257 / 0,011 |
| TR-1208, 8° | 0,5444 / 0,5551 | 0,5509 / 0,0441 | 0,5551 / 0,0438 | mesure 0,499 / 0,021 |

Décaler la maquette d'un quart de maille change la portance de la boîte de
0,2 % (aile AR 4, 4° : 0,2715 contre 0,2720 ; lissée, 0,2849 contre 0,2853) : l'effet du lissage
n'est pas un effet de placement. La géométrie lissée est plus juste, mais
elle ne rapproche pas des références : la portance monte de 0 à 5 % (la
dépression du bord d'attaque arrondi est mieux rendue), la traînée de
surface baisse fortement à 4° sur l'aile AR 4 (succion du bord d'attaque
retrouvée, au niveau de VSPAERO) mais monte à 8°, et bouge de 1 à 3 % sur
TR-1208. Les facettes grossières ne sont donc pas la cause de l'excès de
portance ; elles en masquaient plutôt une petite part. Le lissage reste une
option (`--lisser 1`, case « Lisser les facettes » de l'interface), à
réserver aux fichiers dont `soufflerie geom` montre un bord d'attaque
nettement plus anguleux que le profil, et à réévaluer une fois le
décambrage visqueux rendu.

**Arêtes vives tirées du maillage (`--aretes maillage`).** Le solveur coupe
la loi de paroi aux mailles dont les liens touchent des faces à plus de 60°
l'une de l'autre. Comparé aux arêtes vives du maillage (angle dièdre > 60°,
bords ouverts), ce critère par maille est juste sur l'aile AR 4 : les 97
mailles d'arête de la grille ×4 sont toutes à moins de 1,5 maille d'une
arête du maillage. Sur TR-1208, 14 des 30 sont à plus de 4 mailles de toute
arête : ce sont des mailles du bord d'attaque, plus fin qu'une maille, dont
les liens touchent les deux faces ; avec le lissage, 36 sur 65. L'option ne
coupe la loi de paroi qu'à portée de lien (√2 maille) d'une arête du
maillage. Effet : nul (TR-1208, boîte Cz 0,2864 contre 0,2862 à 4°, 0,5458
contre 0,5444 à 8°, traînée inchangée ; lissage et arêtes ensemble :
0,2888 et 0,5574, comme le lissage seul). Pas de divergence non plus au bord
d'attaque. Gardée en option ; le critère par maille reste le défaut.

**Raffinement du bord d'attaque suivant la flèche : non fait, chiffré.**
Une boîte unique alignée sur les axes couvre presque toute une aile en
flèche. Découpée en tronçons d'envergure qui suivent le bord d'attaque, la
grille fine coûterait sur TR-1208 :

| Grille du bord d'attaque | Une boîte | 4 tronçons | 8 tronçons | 16 tronçons |
| --- | --- | --- | --- | --- |
| ×8 (sur ×4) | 6,2 M mailles | 2,5 M | 1,2 M | 0,6 M |
| ×16 (sur ×8) | 33,9 M | 11,9 M | 4,7 M | 2,0 M |

Huit tronçons rendraient `--raff 3 --ba 1` possible sur TR-1208 (4,7 M
mailles au lieu de 34 M, au-delà de la mémoire du GPU). Mais les grilles
fines forment une chaîne (chaque grille n'a qu'une grille fille) : des
grilles sœurs demandent de reprendre le couplage (`refine3d.rs`) et la
répartition des efforts entre grilles (`gpu3d.rs`). À faire séparément ;
c'est la seule piste géométrique qui attaque la résolution du bord
d'attaque, première cause de la traînée de pression excessive à 8°.

**Écartés sans calcul.** Champs de distance signée et de normales
précalculés : la distance exacte aux triangles (rangés par cases) sert déjà
à la voxelisation, aux fractions de Bouzidi et aux normales, et rien ne
s'en trouverait plus juste. Courbure par triangle pour corriger la pression
pariétale : essayée plus haut (équilibre normal), elle ne vaut pas quand le
rayon tient dans une maille. Remaillage à l'échelle de la grille (facettes
d'une maille) : sans objet, les fractions de Bouzidi sont exactes quelle
que soit la taille des triangles ; seule compte la forme, d'où le lissage.

**Format `.sgeo`.** Un STL ne porte ni unité, ni axes, ni références, ni
incidence : chaque calcul les redemande (`--unite in --lref 16.672 --sref
2024.35`), et une erreur d'unité ou de surface passe inaperçue. Il n'a pas
non plus de connectivité (sommets répétés en simple précision) ni de moyen
de dire quelles arêtes sont vives. Le format `.sgeo` (`lbm/src/sgeo.rs`) est
un fichier texte : en-tête `soufflerie-geometrie 1`, puis `nom`, `unite`,
`ecoulement`, `haut`, `lref`, `sref`, `alpha`, `beta`, `sol` (tous
facultatifs), les sommets en double précision, les faces (orientées vers
l'extérieur) et les arêtes vives (paires de sommets, calculées à l'export
d'après l'angle dièdre et modifiables à la main). La ligne de commande et
l'interface le lisent comme un STL (`--stl fichier.sgeo`, « Charger un
STL… », glisser-déposer) ; ses valeurs servent de défaut, les options les
remplacent ; avec `--aretes maillage`, ses arêtes vives remplacent celles de
l'angle dièdre. Export :

```
soufflerie geom --stl TR1208.stl --unite in --lref 16.672 --sref 2024.35 --alpha 4 --exporter TR1208.sgeo
soufflerie run3d --stl TR1208.sgeo --l 8 --raff 2 --re 4e6 --marge 0.4
```

`--exporter fichier.stl` écrit les triangles nettoyés (lissés avec
`--lisser 1`). Tests : aller-retour au bit près, défauts signalés avec leur
ligne (`sgeo::tests`), même voxelisation qu'avec le STL et arêtes placées
avec la maquette, y compris dans les grilles fines (`bodies3d::tests::
pretraitement_stl`), raffinement lisse fermé et arêtes vives conservées
(`geom3d::tests`).

**Bilan.** Ce qui aide : le contrôle de résolution, qui chiffre ce que voit
le solveur (rayon du bord d'attaque en mailles, facettes, parois minces) et
annonce un bord d'attaque sous la maille, et le format `.sgeo`, qui fixe la
mise en place. Ce qui rend la géométrie plus juste sans rapprocher des
références : le lissage des facettes. Ce qui ne change rien de mesurable :
les arêtes vives tirées du maillage. Ce qui aiderait le plus, et reste à
faire : un raffinement du bord d'attaque en tronçons qui suivent la flèche.

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

## Vitesse de calcul 3D

Le calcul 3D a été accéléré de deux à trois fois sans rien changer au
modèle, au schéma, à la résolution ni à la précision des nombres (la
demi-précision reste la demi-précision).

**Banc de mesure.** `lbm/examples/banc3d.rs` reprend la boucle de
`soufflerie run3d` sur quatre petits cas : `ar4` (aile AR 4 du STL
`banc/wing_ar4.stl`, 4°, Re = 10⁶, L = 16, grille ×2), `ar4b` (idem plus
une couche pariétale ×4 de 16 mailles), `ar4n` (sans raffinement) et
`sphere` (Re = 1000, L = 16, sous-maille sans loi de paroi). Il donne les
durées de mise en place, de calcul et de relevé, les mises à jour de
maille par seconde (Σ n_k 2^niveau par pas de la veine), et à temps fixé
les efforts (échange, surface, boîtes) et des sommes de contrôle du champ
de la grille la plus fine ; `--champs f.bin --champs-ref g.bin` compare
le champ maille par maille, `--profil 1` mesure chaque passe isolément.
`banc/compare.py` rapporte les écarts à la dispersion de plusieurs
références (`banc/ref` : version d'origine, t* = 6), `banc/ab.sh` compare
deux jeux de shaders en lancements alternés (`SOUFFLERIE_WGSL=dossier`
écrit les shaders engendrés dans le dossier, ou les relit s'ils y sont :
essais sans recompiler).

```bash
cargo run --release -p lbm --example banc3d -- ar4b --t 6 --profil 1
```

**Gains** (GTX 1070 seule sur le GPU, t* = 6, 960 pas de la veine ;
ancienne et nouvelle version lancées en alternance) :

| cas | mailles (mises à jour par pas de la veine) | calcul avant | après | gain | mise en place avant | après |
| --- | --- | --- | --- | --- | --- | --- |
| `sphere` | 1,77 M (1,77 M) | 3,18 s | 1,01 s | ×3,1 | 0,75 s | 0,75 s |
| `ar4` | 2,91 M (3,68 M) | 7,65 s | 3,35 s | ×2,3 | 1,3 s | 1,2 s |
| `ar4b` | 4,27 M (9,10 M) | 29,4 s | 8,6 s | ×3,4 | 4,2 s | 2,2 s |
| `ar4r2` (×2 puis ×4) | 5,52 M (14,1 M) | 35 s | 15,2 s | ×2,3 | 2,5 s | 2,5 s |
| `ar4c3` (×2, couches ×4, ×8, ×16), 40 pas | 24,1 M (292 M) | 43,0 s | 11,1 s | ×3,9 | 18 s | 9,2 s |

`ar4c3` reproduit `--bande auto` sur l'aile AR 4 (trois couches de 16
mailles) : c'est le cas qui ressemble le plus aux gros calculs, et celui
qui gagne le plus (la couche la plus fine y fait les deux tiers du pas).
Mises à jour de maille par seconde, toutes grilles : 60 à 530 M avant,
890 à 1 690 M après. Sur un calcul de validation (t* = 20) : `ar4b`
108 s → 31 s, `ar4` ≈ 25 s → 11,6 s (avant les dernières optimisations).

Par optimisation (profil par passe, `ar4` ; durées par lancement) :

| optimisation | passe | avant | après | résultat |
| --- | --- | --- | --- | --- |
| tirage des raccords en deux passes : moments de chaque nœud parent une seule fois (un nœud sert à jusqu'à 64 mailles fantômes), puis interpolation | tirage | 3,9 ms | 0,6 ms | identique au bit |
| index des voisins sans modulo (rangées voisines calculées une fois, début de rangée précalculé) | collision | 11,4 ms | 8,8 ms | identique au bit |
| collision en deux passes : mailles sans lien de paroi, sans le code des liens (recopié dans les 18 directions, il réservait des registres à toutes les mailles), puis mailles pariétales sur leur liste | collision | 8,8 ms | ≈ 6 ms (−24 % en alterné) | identique au bit |
| termes de coefficient nul retirés de la collision (0 · x n'est pas simplifiable en IEEE) | collision | | −3 % | identique au bit |
| voisins de ν̃ sans modulo | ν̃ | 1,2 ms | 0,9 ms | identique au bit |
| shaders de calcul sans les contrôles d'exécution de wgpu (bornes, divisions protégées, compteurs de boucle ; `SOUFFLERIE_CONTROLES=1` les rétablit) | toutes | | faible | identique au bit |
| loi de paroi dans une passe à part (voir ci-dessous) | collision | | faible | reproductible |
| collision par paires de directions opposées : parties paire et impaire calculées une fois pour f_i et f_−i | collision | 1,70 ms | 1,54 ms | arrondi (voir ci-dessous) |
| mailles des grilles parentes recopiées de la grille fine, et dont les six voisins le sont aussi : collision et ν̃ sautés (`SOUFFLERIE_RECOUVERT=0` pour comparer) | toutes, grilles parentes | | 28 % de la grille ×2 sautés sur `ar4r2` (−2 % au total) | identique au bit |
| collision spécialisée par grille : paramètres fixes (SA, sous-maille, bords, diagnostics, stockage) en constantes du shader, pipelines recompilés s'ils changent | collision | | −2 à −7 % selon le cas | identique au bit |
| tirage des raccords et remplissage du premier demi-pas de la grille fine en une passe | raccords | | −1,5 % | identique au bit |
| couches creuses : le cœur de la collision ne lance que les mailles à calculer (un tiers des mailles stockées sont fantômes ou solides) | collision, couche ×16 | 7,2 ms | 6,9 ms | identique au bit |
| collision pariétale : fractions q et codes des normales des 18 liens lus d'un coup au lieu d'un par direction derrière ses branches | collision pariétale, couche ×16 | 2,6 ms | 1,7 ms | identique au bit |
| couches creuses : collision du cœur sans lecture du type des mailles de la liste (toutes fluides sans paroi) | collision, couche ×16 | 6,85 ms | 6,64 ms | identique au bit |
| mise en place des couches : liens de paroi cherchés à moins de 1,5 maille de la surface seulement (un lien de √2 ne peut pas la couper plus loin), fantômes à remplir cherchés en parallèle | mise en place, couche ×16 | 3,9 s | 1,6 s | identique au bit |
| masque d'un bit des mailles fluides ordinaires : collision des grilles boîtes et transport de ν̃ ne lisent le type complet (4 octets) que pour les autres | collision, ν̃ | 1,37 ms, 0,338 ms | 1,35 ms, 0,320 ms | identique au bit |
| mise en place : distance de chaque maille bornée par celle du centre de son bloc (au-delà du plafond sans recherche, recherche arrêtée plus tôt ; test `distance_bornee_exacte`) | mise en place, couche ×16 | 2,7 s | 2,4 s | identique au bit |
| distance à la surface par coques croissantes (même minimum que le parcours complet, test `distance_par_coques_exacte`) | mise en place des couches | 4,3 s | 1,4 s | identique au bit |

**Précision.** Toutes les optimisations sauf deux rendent des champs et des
efforts identiques au bit (comparaison maille par maille des champs de
`ar4`, `ar4b` et `sphere`). L'exception corrige un défaut : la collision
appliquait la loi de paroi en lisant la vitesse des nœuds voisins (relevé à
la hauteur WM_H) pendant que la même passe la réécrivait, tantôt celle du
pas précédent, tantôt celle du pas courant selon l'ordre d'exécution des
groupes. Les calculs avec loi de paroi n'étaient pas reproductibles :
d'un lancement à l'autre, 10⁻⁵ à 10⁻³ d'écart relatif sur les efforts à
t* = 6 (`banc/ref`). La loi de paroi est maintenant une passe sur les
mailles pariétales, avant la collision, sur le champ du pas précédent :
le calcul est reproductible au bit. La seconde, la collision par paires,
calcule la même formule dans un autre ordre : l'arrondi de chaque
population change d'environ 10⁻⁷, ce qui suffit pour que les champs
instantanés divergent ensuite comme d'un lancement à l'autre de l'ancienne
version. Les efforts moyens restent à moins d'une largeur de dispersion
de la moyenne de l'ancienne version (nouvelle version complète,
`banc/long/new2_*`) ; Cz échange et Cz boîte dépassent le plus grand des
six tirages de 3 et 5·10⁻⁵ (10⁻⁴ et 2·10⁻⁴ en relatif) :

| `ar4`, t* = 20, moyenne sur 8 | ancienne version, 6 calculs (min – max) | nouvelle |
| --- | --- | --- |
| Cx échange | 0,05712 – 0,05714 | 0,05713 |
| Cz échange | 0,27212 – 0,27223 | 0,27226 |
| Cx surface | 0,02435 – 0,02435 | 0,02435 |
| Cz surface | 0,28962 – 0,28979 | 0,28969 |
| Cx frottement | 0,00749 – 0,00749 | 0,00749 |
| Cx boîte | 0,07193 – 0,07342 | 0,07332 |
| Cz boîte | 0,26892 – 0,26903 | 0,26908 |

Même constat sur `ar4b` (t* = 20, trois calculs de l'ancienne version,
`banc/long/orig_ar4b_*`) : moyennes à moins d'une dispersion ; seule une
valeur instantanée (Cx au dernier relevé) en sort, de 6·10⁻⁴ en relatif.

**Essais écartés.**

- Bandes absorbantes : `absorb` coûte environ 15 % de la collision, mais
  une sortie anticipée à l'intérieur de la veine n'y change rien (le coût
  vient des 19 branches `sig > 0` intercalées, qui occupent des
  registres) ; regrouper ces 19 branches en une seule change l'arrondi
  (le compilateur du pilote contracte autrement les multiplications-
  additions) sans gain mesurable : écarté.
- Vecteurs AVX2 pour le moteur 2D sur processeur (`target-cpu=native`) :
  +60 % sur un fil, rien sur 12 fils (limité par la mémoire) ; écarté
  (binaire non portable).
- Taille des groupes de travail de la collision (`SOUFFLERIE_WG`, 32 à
  256), GPU libre : 64 à 256 équivalents, 32 plus lent ; 64 gardé.
- Parcours de la veine par tranches de quelques rangées y sur toute la
  hauteur z (voisins en z relus pendant qu'ils sont dans le cache) :
  ≤ 5 %, dans le bruit ; écarté.
- Mailles fantômes des raccords triées par parité de position (pochoirs de
  même taille dans un groupe) : deux fois plus lent (accès dispersés) ;
  écarté.
- Mêmes fantômes rangés par bloc de 8³ puis par parité : plus lent
  aussi. L'interpolation des fantômes (10 % du pas de `ar4c3`) est limitée
  par ses lectures, environ 800 octets par fantôme ; seul un regroupement
  des fantômes d'une même maille parente en mémoire partagée les
  réduirait (de moitié environ), refonte lourde non faite.
- Divisions entières de `pos_of` remplacées par une estimation flottante
  corrigée (quotient exact) : identique au bit, aucun gain.
- Écriture des seules vitesses (sans ρ) dans le champ à chaque pas :
  aucun gain, une écriture partielle de 16 octets coûte autant qu'une
  complète ; il faudrait des tableaux séparés (≈ 1,5 %).
- Voisins dans le même bloc de 8³ sans passer par la table des blocs
  (collision et recopie des couches) : plus lent (branche), la table est
  dans le cache.
- Normales des liens décodées une fois à la mise en place : plus lent
  (lectures dépendantes au lieu de calcul).
- Populations de la maille préchargées dans la collision pariétale : plus
  lent.
- Transport de ν̃ et collision limités aux mailles lues sur les couches
  (sans solides ni fantômes jamais lus) : aucun gain pour ν̃, limité par
  la mémoire des mailles actives.
- Champ de vitesses en plans séparés (u_x, u_y, u_z, puis ρ écrit au
  dernier pas seulement) : ν̃ −6 %, mais collision des couches +5 % (trois
  flux d'écriture au lieu d'un) ; écarté.
- Mailles pariétales lancées par motif de liens (groupes homogènes) :
  collision pariétale 1,5 fois plus lente (localité perdue) ; écarté.
- Boucles de 19 des raccords déroulées à la génération : à peine plus
  rapide, et plus identique au bit (le compilateur replie autrement) ;
  écarté.
- Distance à la surface par listes de triangles par bloc de 8³ : plus
  lent (listes trop longues) ; écarté.
- Conversions demi-précision ↔ simple précision par opérations entières
  sur les bits (au lieu de l'unité de conversion, lente sur Pascal) :
  identiques au bit (le matériel arrondit au plus proche pair), mais 3 %
  plus lentes ; écarté.
- Accès robuste aux tampons de Vulkan (`robustBufferAccess`, que wgpu
  active sans option pour le couper) : copie de wgpu-hal sans lui, aucune
  différence mesurable (< 0,3 %) ; rien à changer.
- Barrières entre lancements (wgpu en place une à chaque passe qui écrit)
  et encodage des commandes : moins de 1 % du pas.
- Listes exactes de triangles par bloc pour la distance (d_centre + deux
  demi-diagonales) : couche ×16 plus rapide, couches ×4 et ×8 plus
  lentes ; écarté.

**Ce qui reste.** GPU libre, la collision de la veine prend 1,40 ms ;
sans calcul ni tirage des voisins (mêmes lectures et écritures), elle
prendrait 1,11 ms (≈ 200 Go/s, 80 % de la bande passante de la carte).
Sur `ar4c3`, le pas se répartit en collision de la couche ×16 (cœur et
mailles pariétales, 54 %), transport de ν̃ (14 %), raccord ×8 → ×16
(19 %, dont l'interpolation des fantômes 10 %) et le reste. Pistes non
faites : interpolation des fantômes par maille parente en mémoire
partagée (≈ 5 %) ; transport de ν̃ dans la collision du pas suivant
(≈ 5 %, mais un second champ de vitesses, 16 octets de plus par maille,
et un décalage d'un pas à reporter dans les raccords) ; vitesses rangées
en tableaux séparés (≈ 2 %).


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
Kármán–Tsien, détection de la zone supersonique locale) et `--couplage
newton` (couplage visqueux–non visqueux simultané, voir « Couplage
simultané (Newton) en 2D » ; par défaut `faible`). Sans sous-commande,
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
│   ├── src/ibl3d.rs    couche limite intégrale sur un maillage triangulé, vitesse de transpiration
│   ├── src/bodies3d.rs maquettes 3D (sphère, cube, cylindre, aile)
│   ├── src/stl.rs      lecture STL, nettoyage, appartenance, distance à la surface
│   ├── src/geom3d.rs   prétraitement géométrique : arêtes vives, raffinement lisse, contrôle par tranche
│   ├── src/sgeo.rs     format de géométrie .sgeo (maillage indexé, mise en place, arêtes vives)
│   ├── src/geometry.rs mailles solides, liens de Bouzidi, distance et normale à la paroi
│   ├── src/shapes.rs   maquettes
│   ├── src/tunnel.rs   veine d'essai, mise en vitesse, turbulence, transition, coefficients
│   ├── src/stats.rs    moyennes, périodes, Strouhal
│   ├── examples/       programmes de diagnostic (voir « Programmes de diagnostic »)
│   └── tests/          validation, panneaux, aile finie, 3D, comparaison GPU/CPU
├── banc/               banc de mesure de la vitesse 3D : maquette, références, comparaisons
├── validation/         résultats de validation 3D (CSV)
├── tools/validation3d.sh  calculs de validation 3D
└── app/                application : interface (egui) et commandes
    ├── src/main.rs     ligne de commande (run, polar, panel, wing, run3d, geom, bench)
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
