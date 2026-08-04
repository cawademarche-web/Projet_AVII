# ACTU-F502 — Checklist de contrôle du code 

À garder ouvert pendant que Claude Code tourne. Dès qu'une valeur ne matche pas un signal attendu ci-dessous : stop, on regarde le code avant de continuer.

## A1 — Taux MLE (µ̂ₓ = Dₓₜ/ETRₓₜ)
- µ̂ₓ croissant avec l'âge, quasi linéaire en échelle log au-delà de ~40 ans (Gompertz).
- IC 95% étroits aux âges à forte exposition (30–70 ans), larges aux grands âges.
- 🚩 IC plus large à 50 ans qu'à 90 ans → exposition mal calculée.

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

## C — VAP et primes
- VAP[B] < VAP[A] toujours (la rente temporaire ampute la queue longue, la plus incertaine).
- VAP décroît quand t augmente ; VAP[A] plus sensible à t que VAP[B] (effet duration).
- Ordre de grandeur : rente viagère à 65 ans à t=3% ≈ 15–18× la rente annuelle.
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
