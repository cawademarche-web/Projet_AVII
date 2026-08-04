# PIPELINE — Projet AVII (ACTU-F502)

Données : HMD Autriche, période 1947–2023, âges 0–110+ (Deaths_1x1, Exposures_1x1).
Numérotation alignée sur le brief officiel : A.1 (MLE), A.2 (discussion), A.3 (Lee-Carter), B, C, D.

| Section | Script | Méthode | Outputs (chiffres clés) | TODO / à justifier |
|---------|--------|---------|-------------------------|--------------------|
| Préparation | `scripts/01_donnees_taux.R` | Lecture HMD, matrices D_{x,t} et ETR_{x,t} H/F | `resultats/01_donnees.rds` — 4 matrices 111 × 77 (âges 0–110 × années 1947–2023), 0 NA | — |
| A.1 (taux MLE + IC 95 %) | `scripts/01_donnees_taux.R` | µ̂_x(t) = D_{x,t}/ETR_{x,t} (MLE Poisson) ; Var(µ̂) = µ̂/ETR ⇒ IC 95 % = µ̂ ± 1,96·√(µ̂/ETR) | `resultats/01_taux_mle.rds` (6 matrices) ; µ̂(2022) : 40/65/90 ans = 0,00138 / 0,01418 / 0,21308 (H) et 0,00052 / 0,00712 / 0,16410 (F) ; contrôle croisé vs Mx publié — écart max 0,1176471 (H et F) sur toutes cellules finies, 1,713452e-06 (H) et 2,035366e-06 (F) sur ETR ≥ 1000 ; largeur relative d'IC en U (âges 25/50/70/88/100) : 0,68238 / 0,29890 / 0,13681 / 0,10414 / 0,45880 (H) et 1,23961 / 0,39598 / 0,17324 / 0,08984 / 0,20805 (F), minimum à 88 ans ; identité largeur relative = 2·1,96/√D vérifiée TRUE (H et F) ; figures `fig_A1_taux_2022_ic.png`, `fig_A1_taux_evolution.png` | - A.1 illustrée en profils périodiques (2022 + années de comparaison) —
  décision PROVISOIRE, confirmation demandée à l'assistant. La dimension
  cohorte apparaît avec les projections (A.3.v, B.1). Ne pas solder tant
  que la réponse n'est pas reçue. - Table de survie tronquée à ℓ₁₀₂ (µ̂ défini jusqu'à 101 pour les deux
  sexes : ETR = 0 dès l'âge 102 chez les femmes en 1949-1951). Impact
  chiffré par comparaison 100 / 101 / 102 / 103 — à mentionner dans le rapport. - Borne de troncature uniforme à 102 sur toutes les années, alors que la
  contrainte ETR = 0 à l'âge 102 ne concerne que 1949-1951 chez les femmes ;
  une borne annuelle serait plus fine mais l'écart est de ~3 centièmes d'année
  sur e₆₅ — à mentionner.|
| A.2 (discussion : espérance de vie, âge médian/IQR, expansion-rectangularisation) | `scripts/01_donnees_taux.R` (indicateurs) | Table de survie périodique ℓ_x = exp(−Σ_{k<x} µ_k) ; e_x ≈ Σ_{k≥1} ₖp_x + ½ (curtate + ½) ; quantiles de l'âge au décès par interpolation de ℓ_x | `resultats/01_indicateurs.rds` ; 1947 → 2023 : e₀ 58,87 → 79,45 (H) et 63,96 → 84,23 (F) ; e₆₅ 11,97 → 18,40 (H) et 13,51 → 21,57 (F) ; médiane 67,85 → 82,43 (H) et 72,95 → 87,08 (F) ; IQR 27,94 → 15,79 (H) et 22,88 → 12,70 (F) ; sensibilité troncature 2023 (bornes 100→103) : e₆₅ 18,388 → 18,398 (H), 21,545 → 21,578 (F) ; 4 figures `fig_A2_courbes_survie.png`, `fig_A2_courbe_deces.png`, `fig_A2_esperance_vie.png`, `fig_A2_mediane_iqr.png` | Rédaction de la discussion (expansion / rectangularisation) = rapport |
| A.3 (Lee-Carter : fit, résidus, projection, bootstrap N=5000) | `scripts/02_leecarter_stmomo.R` | — | — | — |
| B (5000 trajectoires κ, IC, comparaison avec bootstrap) | `scripts/03_simulation_kappa.R` | — | — | — |
| C (VAP rentes [A]/[B], primes) | `scripts/04_pricing_vap.R` | — | — | — |
| D (BOF, SCR 2022/2023, sensibilité) | `scripts/05_solvabilite_scr.R` | — | — | — |

## TODO globaux

- Plage d'âges de calibration LC : borne haute à fixer bien en dessous de 110
  (groupe ouvert) — à justifier (A.3.i).
- Période de calibration : les années 1947–1950 (immédiat après-guerre)
  sont-elles à exclure ? — à justifier (A.3.i).
- Lecture de « les cohortes concernées à partir de 2022 » en A.1 et B.1 :
  profil périodique 2022 vs diagonale de cohorte — question posée à
  l'assistant, réponse en attente.
- Bornes basses d'IC ≤ 0 aux grands âges (10 âges en 2022 chez les hommes,
  9 chez les femmes) : attendu (approximation normale non contrainte à ℝ⁺),
  masquées en échelle log — à mentionner dans le rapport.
- Choix des années de comparaison 1947 / 1972 / 1997 / 2022 (pas de 25 ans)
  — à justifier dans le rapport.
- µ̂ n'est pas strictement croissant en âge sur 40–95 ans (5 décrochages
  locaux en 2022, H comme F) : bruit d'échantillonnage des taux bruts ; la
  quasi-linéarité en échelle log tient (R² de log µ̂ ~ âge = 0,9966 H /
  0,9935 F) — à mentionner comme motivation du lissage Lee-Carter.
