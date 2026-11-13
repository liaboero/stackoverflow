# Soufflerie de poche — version native (Rust)

Fork natif de l'application web [`soufflerie-lbm`](../soufflerie-lbm) : même
soufflerie numérique 2D Lattice Boltzmann (D2Q9), en Rust, pour ordinateur
(Linux, Windows, macOS). Interface graphique et ligne de commande.

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

Le moteur GPU donne les mêmes efforts que le processeur avec Spalart–Allmaras
(écart ≤ 7·10⁻⁵ sur un profil à Re = 10⁶ et un cylindre tournant à 10⁵).
Spalart–Allmaras coûte environ 45 % de la vitesse sur processeur (131 au lieu
de 239 millions de mailles par seconde sur 4 cœurs), un peu moins sur GPU.

Comme dans l'application web, la portance et le décollement sont
qualitatifs à ces résolutions, et la traînée LBM d'un profil mince reste
approximative : la référence « panneaux + couche limite » est la bonne
valeur de traînée en écoulement attaché.

Pas encore porté : aile finie, soufflerie 3D, balayage automatique dans
l'interface (les commandes `polar` et `panel` le remplacent), unités
physiques (fluide, corde, vitesse) et Mach réel dans l'interface.

## Installation

Installer Rust (<https://rustup.rs>), puis, dans ce dossier :

```sh
cargo run --release                     # interface graphique
cargo run --release -- --shape naca --naca 2412 --alpha 6 --re 1e6 --l 24 --graph cp   # interface réglée
cargo run --release -- run --shape naca --naca 2412 --alpha 6 --re 1000 --l 32
cargo run --release -- run --shape naca --naca 0012 --alpha 4 --re 1e6 --l 32 --turb sa
cargo run --release -- polar --shape naca --naca 0012 --re 2000 --l 24 --alpha -4:12:2 > polaire.csv
cargo run --release -- panel --naca 2412 --re 3e6 --alpha -4:18:1 > panneaux.csv   # instantané
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

Linux : l'interface demande `libxkbcommon-x11` et un pilote Vulkan ou
OpenGL (rendu wgpu, repli OpenGL automatique ; `SOUFFLERIE_RENDU=glow`
force OpenGL). Sans interface :
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
│   ├── src/geometry.rs mailles solides, liens de Bouzidi, distance et normale à la paroi
│   ├── src/shapes.rs   maquettes
│   ├── src/tunnel.rs   veine d'essai, mise en vitesse, turbulence, transition, coefficients
│   ├── src/stats.rs    moyennes, périodes, Strouhal
│   └── tests/          validation, panneaux, comparaison GPU/CPU
└── app/                application : interface (egui) et commandes
    ├── src/main.rs     ligne de commande (run, polar, panel, bench)
    ├── src/gui.rs      interface graphique
    └── src/worker.rs   fil de calcul de l'interface (image, traceurs)
```

## Pistes

- Porter l'aile finie (la méthode des panneaux, déjà portée, fournit les
  polaires de section) et la 3D (le solveur WGSL 3D de l'application web se
  transpose presque tel quel).
- Stockage en demi-précision et propagation en place sur GPU (deux à quatre
  fois moins de mémoire), raffinement local de maillage.
- Version mobile : egui fonctionne sur Android ; il faudrait l'empaqueter
  (cargo-apk / xbuild). Sur téléphone, l'application web reste la plus
  simple à installer.
