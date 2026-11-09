# Soufflerie de poche

Petite soufflerie numérique 2D pour mobile, basée sur la méthode de
Boltzmann sur réseau (Lattice Boltzmann, D2Q9). Elle calcule en direct les
efforts aérodynamiques sur des formes simples : **Cx** (traînée), **Cz**
(portance), **Cm** (tangage) et le nombre de **Strouhal**.

C'est une application web progressive (PWA) : HTML/CSS/JavaScript sans
dépendance ni étape de compilation. Elle s'installe sur l'écran d'accueil
d'un téléphone (iOS ou Android) et fonctionne ensuite hors ligne.

## Fonctions

- **Maquettes** : cylindre (avec rotation, effet Magnus), carré, ellipse
  (épaisseur réglable), plaque plane, profil NACA à 4 chiffres, triangle.
  Incidence α réglable (positive = cabrée).
- **Écoulement** : Reynolds de 10 à 5000, trois résolutions (12, 16 ou
  24 mailles par longueur de référence).
- **Visualisation** : vorticité, vitesse ou coefficient de pression Cp,
  traceurs (lignes d'émission), force résultante moyenne et ses composantes,
  trois cadrages (maquette, sillage, veine entière). Palettes adaptées au
  thème clair ou sombre.
- **Mesures** : moyennes glissantes calculées sur un nombre entier de
  périodes quand l'écoulement est périodique, amplitude des oscillations,
  Strouhal, état de convergence (transitoire, stationnaire, périodique).
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
| Collision | BGK ; viscosité de sous-maille de Smagorinsky (C_s = 0,15) quand τ < 0,51 |
| Paroi | Rebond à mi-chemin interpolé de Bouzidi (forme exacte, pas l'escalier), correction de paroi mobile pour la rotation |
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

`npm test` (Node ≥ 18, ~1 min 20 s) vérifie :

| Cas | Calcul | Référence (milieu infini) |
| --- | --- | --- |
| Écoulement uniforme sans obstacle | écart < 1·10⁻⁶ | 0 |
| Cylindre, Re = 20 | Cx = 2,38 | 2,05 – 2,09 |
| Cylindre, Re = 100 | Cx = 1,49 ; St = 0,176 ; Cz' = 0,37 | 1,33 – 1,35 ; 0,164 – 0,166 ; 0,33 |
| NACA 0012, α = 0° | Cz = 0,000 | 0 |
| NACA 0012, α = ±4°, Re = 500 | Cz = ±0,172 (antisymétrie exacte) | — |
| Cylindre tournant ΩR/U = 1, Re = 100 | Cz = 2,59 | ≈ 2,5 (Mittal & Kumar 2003) |
| Cylindre, Re = 5000 | stable (LES) | — |

Les écarts sur le cylindre viennent surtout du **confinement** du domaine
(8 % de blocage, entrée à 8 diamètres) : dans un domaine de 40 L × 24 L, le
même solveur donne Cx = 1,43 et St = 0,167 à Re = 100. La convergence en
maillage est bonne grâce à la paroi interpolée (Re = 20 : Cx = 2,750 / 2,728 /
2,715 pour 16 / 24 / 32 mailles par diamètre dans l'ancien domaine 16 L × 8 L).

## Limites

- Calcul **2D** : pas de transition, pas de turbulence 3D. Au-delà de
  Re ≈ 1000 le résultat est qualitatif.
- Résolution modeste : au-delà de Re ≈ 300 la couche limite fait moins
  d'une maille ; l'application l'indique.
- Écoulement faiblement compressible (Ma ≈ 0,17 par défaut, réglable).
- Calcul sur un seul cœur du processeur : ≈ 20 à 30 millions de mailles mises
  à jour par seconde sur un téléphone récent, soit 1 à 3 temps convectifs
  par seconde en résolution standard.

## Structure

```
soufflerie-lbm/
├── index.html            interface
├── css/app.css           styles (thèmes clair et sombre)
├── js/lbm.js             solveur D2Q9 (sans DOM, testable sous Node)
├── js/shapes.js          géométries et placement en incidence
├── js/sim.js             veine d'essai : domaine, démarrage, coefficients
├── js/stats.js           moyennes, périodes, Strouhal
├── js/render.js          rendu du champ, traceurs, maquette, force
├── js/chart.js           graphiques (historique, polaire)
├── js/app.js             logique de l'interface
├── sw.js, manifest.webmanifest, icons/   PWA
├── tests/validate.mjs    validation physique
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

- Solveur sur GPU (WebGL2/WebGPU) pour des grilles 10 à 50 fois plus fines.
- Collision MRT ou régularisée pour les Reynolds élevés.
- Correction de blocage, mode dimensionnel (corde en m, vitesse en m/s,
  air ou eau → efforts en N/m).
- Dessin libre de la maquette au doigt.
