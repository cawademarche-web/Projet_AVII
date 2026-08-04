# PIPELINE — Projet AVII (ACTU-F502)

Données : HMD Autriche, période 1947–2023, âges 0–110+ (Deaths_1x1, Exposures_1x1).
Numérotation alignée sur le brief officiel : A.1 (MLE), A.2 (discussion), A.3 (Lee-Carter), B, C, D.

| Section | Script | Méthode | Outputs (chiffres clés) | TODO / à justifier |
|---------|--------|---------|-------------------------|--------------------|
| Préparation | `scripts/01_donnees_taux.R` | Lecture HMD, matrices D_{x,t} et ETR_{x,t} H/F | `resultats/01_donnees.rds` — 4 matrices 111 × 77 (âges 0–110 × années 1947–2023), 0 NA | — |
| A.1 (taux MLE + IC 95 %) | `scripts/01_donnees_taux.R` | — | — | Encore à coder dans ce script |
| A.2 (discussion : espérance de vie, âge médian/IQR, expansion-rectangularisation) | `scripts/01_donnees_taux.R` (indicateurs) | — | — | — |
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
