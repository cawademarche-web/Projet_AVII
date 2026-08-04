# ACTU-F502 — Checklist de contrôle du code 

À garder ouvert pendant que Claude Code tourne. Dès qu'une valeur ne matche pas un signal attendu ci-dessous : stop, on regarde le code avant de continuer.

> **Note de correspondance avec le brief officiel** : la numérotation de cette
> checklist est antérieure au brief. Ici « A1 » = **A.1** du brief (taux MLE)
> et « A2 » = **A.3** du brief (Lee-Carter). La discussion **A.2** du brief
> (espérance de vie périodique, âge médian/IQR, expansion-rectangularisation)
> n'a pas de section dédiée dans cette checklist.

## A1 — Taux MLE (µ̂ₓ = Dₓₜ/ETRₓₜ)
- µ̂ₓ globalement croissant avec l'âge, quasi linéaire en échelle log au-delà
  de ~40 ans (Gompertz). Quelques décrochages locaux = bruit d'échantillonnage,
  attendus, et argument en faveur du lissage Lee-Carter.
- Précision des IC : largeur relative = 2·1,96·√(µ̂/ETR)/µ̂ = 2·1,96/√D.
  Elle ne dépend QUE du nombre de décès, jamais de l'exposition — ce n'est donc
  pas un test d'exposition. Profil attendu : courbe en U, large aux âges jeunes
  (peu de décès), MINIMALE autour du mode des décès (~85–90 ans), large à
  nouveau aux âges extrêmes (peu de survivants).
- 🚩 Largeur relative ≠ 2·1,96/√D → erreur de formule dans le calcul de l'IC.
- TEST D'EXPOSITION : contrôle croisé de µ̂ contre Mx_1x1.txt publié
  par la HMD. Il valide la chaîne de mise en forme — appariement des colonnes
  Male/Female, recodage 110+, ordre de remplissage de la matrice, alignement
  âge/année. Une inversion H/F ou un décalage d'indexation produirait un écart
  d'ordre 1.
- 🚩 Écart max > 1e-4 sur les cellules à ETR ≥ 1000 → défaut de mise en forme.

## A2 — Lee-Carter (ln µₓ(t) = αₓ + βₓ·κₜ)
- κₜ décroît ~linéairement dans le temps.
- βₓ majoritairement positifs.
- Σβₓ=1 et Σκₜ=0 **vérifiées numériquement dans le code**, pas juste supposées.
- Heatmap des résidus : pas de bandes diagonales marquées (sinon = effet cohorte non capturé — à *documenter* dans le rapport, pas à corriger en douce).
- 🚩 Contraintes non vérifiées, ou κ ne décroît pas → fit cassé.

## B — Simulation de trajectoires (5000 × κ)
- Les 5000 trajectoires simulent **κ seul**, paramètres (α, β) fixés au fit central — **pas** un re-bootstrap.
- L'éventail des trajectoires s'élargit en cône avec l'horizon (variance croît avec h).
- 🚩 La fonction B ressemble à la fonction bootstrap de A → confusion incertitude de paramètres / incertitude de projection, ruine la comparaison demandée en rapport.
- 🚩 IC de A.3 plus étroits que ceux de B = bug certain
## C — VAP et primes
- VAP[B] < VAP[A] toujours (la rente temporaire ampute la queue longue, la plus incertaine).
- VAP décroît quand t augmente ; VAP[A] plus sensible à t que VAP[B] (effet duration).
- Ordre de grandeur : rente viagère à 65 ans à t=3% ≈ 13–16× la rente annuelle.
- 🚩 VAP[B] > VAP[A], ou VAP insensible à t → logique de calcul cassée.

## D — Solvabilité / SCR (le plus risqué)
- **Trois taux, trois rôles stricts, jamais interchangeables :**
  - t = 3% → déjà utilisé en amont (prime, section C). N'apparaît **plus** dans la boucle D.
  - r = 2% → actualise BE à chaque date, et actualise BOF₂₀₂₃ dans la formule du SCR.
  - i = 4% → fait croître les **actifs uniquement**, une fois par an.
- Toutes les valeurs calculées en t⁺ (après paiement primes/prestations de l'année).
- BOF₂₀₂₂ **vraisemblablement négatif** (r < t ⇒ BE actualisé plus généreusement que la prime n'a été tarifée ⇒ BE₂₀₂₂ > A₂₀₂₂). C'est un résultat attendu, pas un bug.
- SCR[A] > SCR[B] (portefeuille A porte plus de risque de longévité).
- Sensibilité i=1% doit dégrader nettement le tableau (spread i−r devient négatif).
- 🚩 **BOF₂₀₂₂ positif** → vérifier en premier que `calcule_BE()` utilise bien r=2% et pas t=3% ou i=4% par confusion de variable. C'est le bug le plus probable, pas un problème de tarification.
- 🚩 La boucle D re-simule ses propres trajectoires au lieu de réutiliser celles de B → gaspillage et incohérence entre sections.

## Réflexe à chaque section, à mettre dans le prompt initial de Claude Code
- Seed fixé dès la première ligne.
- Une fonction = un concept actuariel (`vap_rente()`, `simule_kappa()`, `calcule_BE()` — pas un script monolithique).
- Chaque bloc commenté avec la notation du cours (ex. `# ceci est ₖpₓ`, `# tirage ε ~ N(0,σ²) de la RWD`).
