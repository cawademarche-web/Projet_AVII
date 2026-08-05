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
| A.3 (Lee-Carter : fit, résidus, projection, bootstrap N=5000) | `scripts/02_leecarter_stmomo.R` | LC Poisson via StMoMo `fit(lc(link="log"), Dxt=, Ext=, ages=, years=)` — expositions **centrales** (validées par le contrôle croisé d'A.1 contre le Mx publié) ; poids w_{x,t} = 0 sur les cellules ETR = 0 ; projection `forecast(h=50, kt.method="mrwd", jumpchoice="fit")` ; bootstrap `type="semiparametric", deathType="observed"` | voir les sous-sections ci-dessous | voir TODO globaux |
| B (5000 trajectoires κ, IC, comparaison avec bootstrap) | `scripts/03_simulation_kappa.R` | `simulate(LCfit, nsim = 5000, h = 50, kt.method = "mrwd", jumpchoice = "fit")` — α, β, dérive d̂ et volatilité σ̂ **gelés au fit central**, seul l'aléa futur de κ est simulé : **UNE** source d'incertitude contre trois en A.3.vi. Diagonale de cohorte µ_{65+k}(2022+k) extraite par trajectoire | `resultats/03_diagonales_cohorte.rds` (4 matrices 37 × 5000, 0 NA) ; e₆₅ de cohorte sous B : 19,762 [18,823 ; 20,660] (H) et 23,149 [22,123 ; 24,162] (F) ; largeur relative **9,30 % (H)** et **8,81 % (F)** contre 9,83 % / 9,12 % en A.3.vi ⇒ **rapport A.3/B = 1,058 (H) et 1,036 (F)** ; contrôles RWD à h = 50 : sd(κ) = 6,2391 vs σ̂·√50 = 6,3100 (H), 7,1962 vs 7,0641 (F) ; 5 figures `fig_B1_taux_cohorte.png`, `fig_B2_trajectoires_kappa.png`, `fig_B2_ic_taux.png`, `fig_B3_comparaison_ic.png`, `fig_B3_e65_comparaison.png` | voir sous-section « Section B — détail » |
| C (VAP rentes [A]/[B], primes) | `scripts/04_pricing_vap.R` | **100 % custom, base R, aucun `library()`.** Chaîne µ → ₖp₆₅ → VAP : `probas_survie_trajectoire()` = `apply(M, 2, function(mu) cumprod(exp(-mu)))` — cumprod **colonne par colonne**, force constante par morceaux, même convention qu'A.1 et qu'`esperance_vie_cohorte()` des scripts 02/03 ; `vap_rente(S, taux, n)` = `colSums(v^(1:n) · S[1:n, ])`, **terme échu**, somme de k = 1 à n (n = 37 pour [A], borne ℓ₁₀₂ ; n = 15 pour [B]). Un seul taux dans le script : `taux_tarif = 3 %` | `resultats/04_vap.rds` (3,09 Mo) ; **VAP[A] = 13,7899 (H) et 15,6995 (F)**, VAP[B] = 10,4086 (H) et 11,1011 (F) ; variances 0,061949 / 0,067253 ([A]) et 0,003615 / 0,001602 ([B]) ; primes uniques = ces moyennes (rente annuelle de 1) ; sensibilité 1 % → 5 % : **−33,48 % [A] contre −23,92 % [B]** (H) et −35,96 % contre −24,48 % (F) ; sensibilité COVID à 3 % : **+3,48 % / +2,52 %** sur [A], +1,60 % / +0,89 % sur [B] ; 2 figures `fig_C_hist_vap.png`, `fig_C_sensibilite_taux.png` | voir sous-section « Section C — détail » |
| D (BOF, SCR 2022/2023, sensibilité) | `scripts/05_solvabilite_scr.R` | — | — | — |

## Section A.3 — détail (chiffres recopiés de la console du run `n_boot = 5000`)

### Expositions nulles et matrice de poids
`ETR = 0` par plafond d'âge, années 1947-2023 : **âges 0-100 → 0 cellule (H et F)** ;
**âges 0-102 → 0 (H) et 3 (F)** ; âges 0-105 → 32 (H) et 26 (F).
Les 3 cellules féminines sont toutes à **l'âge 102, années 1949, 1950, 1951**.
Avec le plafond retenu (102), `w_{x,t}` est donc active pour les femmes seulement.
Ordres de grandeur aux grands âges en 2022 : H, âge 101 → D = 48, ETR = 70,1, µ̂ = 0,685 ;
âge 102 → D = 21, ETR = 30,6, µ̂ = 0,686. F, âge 101 → D = 232, ETR = 429,96, µ̂ = 0,540 ;
âge 102 → D = 123, ETR = 221,64, µ̂ = 0,555.
NB : `fit()` zéro-pondère déjà automatiquement les ETR ≤ 0 ; `w_{x,t}` est construite
explicitement pour être comptée et documentée. 12 avis émis au bloc 2a = 2 warnings
(exposition non positive + valeur manquante) × les 6 variantes féminines démarrant en 1947.

### A.3.i — grille de 18 variantes par sexe (36 fits, **tous convergés**)

⚠️ **d̂ n'est PAS comparable d'une plage d'âges à l'autre.** Sous la contrainte
d'identification Σβ_x = 1, l'échelle de κ dépend du nombre d'âges de la plage : 103 âges
pour 0:102 contre 43 pour 60:102. Le rapport ~3,6 entre les d̂ de ces deux plages est un
**pur artefact d'identification**, pas une différence de dynamique. La quantité invariante
est **β₆₅ · d̂**, taux annuel d'amélioration de ln µ₆₅, et c'est elle qui pilote la VAP
d'une rente à 65 ans. Les comparaisons **à plage d'âges fixée** (2019 vs 2023, H vs F)
restent valides sur d̂.

| sexe | ages | période | R²(κ~t) | d̂ | β₆₅ | **β₆₅·d̂** | σ̂(Δκ) | ρ₁(Δκ) | déviance | BIC | RMSE commune |
|---|---|---|---|---|---|---|---|---|---|---|---|
| h | 0:102 | 1947-2023 | 0,9530 | −1,9739 | 0,007786 | −0,01537 | 2,8022 | −0,1424 | 25305 | 82342 | 0,13314 |
| h | 0:102 | 1947-2019 | 0,9448 | −2,1205 | 0,007708 | −0,01634 | 2,3916 | −0,2140 | 23804 | 78043 | 0,13216 |
| h | 0:102 | 1970-2023 | 0,9851 | −2,1934 | 0,008273 | −0,01815 | 2,6353 | −0,1462 | 9337 | 49447 | 0,13335 |
| h | 0:102 | 1970-2019 | 0,9926 | −2,4101 | 0,008246 | −0,01987 | 2,0463 | −0,3361 | 8209 | 45515 | 0,13253 |
| h | 0:102 | 1980-2023 | 0,9769 | −2,2418 | 0,008227 | −0,01844 | 2,7366 | −0,1402 | 6458 | 39141 | 0,13098 |
| h | 0:102 | 1980-2019 | 0,9945 | −2,5124 | 0,008180 | −0,02055 | 2,0512 | −0,3836 | 5500 | 35375 | 0,13029 |
| h | 50:102 | 1947-2023 | 0,9094 | −0,5600 | 0,024264 | −0,01359 | 1,3701 | −0,2960 | 8514 | 41847 | 0,13275 |
| h | 50:102 | 1947-2019 | 0,8945 | −0,6037 | 0,024168 | −0,01459 | 1,2855 | −0,3409 | 7809 | 39355 | 0,13183 |
| h | 50:102 | 1970-2023 | 0,9824 | −0,7974 | 0,022954 | −0,01830 | 1,1213 | −0,2135 | 5148 | 28973 | 0,13326 |
| h | 50:102 | 1970-2019 | 0,9901 | −0,8810 | 0,022775 | −0,02007 | 0,8983 | −0,3525 | 4395 | 26430 | 0,13244 |
| h | 50:102 | 1980-2023 | 0,9746 | −0,8326 | 0,022135 | −0,01843 | 1,1731 | −0,2275 | 3783 | 23367 | 0,13086 |
| h | 50:102 | 1980-2019 | 0,9934 | −0,9418 | 0,021808 | −0,02054 | 0,9010 | −0,4307 | 3146 | 20940 | 0,13018 |
| h | 60:102 | 1947-2023 | 0,8994 | −0,4077 | 0,032192 | −0,01312 | 1,1249 | −0,3299 | 6015 | 33095 | 0,13254 |
| h | 60:102 | 1947-2019 | 0,8836 | −0,4444 | 0,031743 | −0,01411 | 1,0705 | −0,3754 | 5578 | 31176 | 0,13166 |
| **h** | **60:102** | **1970-2023** | **0,9799** | **−0,6023** | **0,030263** | **−0,01823** | **0,8924** | **−0,2189** | **3864** | **23321** | **0,13312** |
| h | 60:102 | 1970-2019 | 0,9902 | −0,6732 | 0,029640 | −0,01995 | 0,7296 | −0,3585 | 3414 | 21388 | 0,13229 |
| h | 60:102 | 1980-2023 | 0,9676 | −0,6097 | 0,029897 | −0,01823 | 0,9243 | −0,2255 | 2862 | 18865 | 0,13075 |
| h | 60:102 | 1980-2019 | 0,9908 | −0,6993 | 0,029008 | −0,02029 | 0,7294 | −0,4160 | 2465 | 16984 | 0,13007 |
| f | 0:102 | 1947-2023 | 0,9826 | −2,2679 | 0,007843 | −0,01779 | 3,3054 | −0,3488 | 16883 | 72074 | 0,07401 |
| f | 0:102 | 1947-2019 | 0,9819 | −2,4241 | 0,007881 | −0,01910 | 3,0074 | −0,4292 | 15687 | 68208 | 0,07363 |
| f | 0:102 | 1970-2023 | 0,9825 | −2,0481 | 0,008498 | −0,01740 | 2,8116 | −0,3216 | 7954 | 46577 | 0,06999 |
| f | 0:102 | 1970-2019 | 0,9941 | −2,2272 | 0,008649 | −0,01926 | 2,3891 | −0,4625 | 6944 | 42891 | 0,06958 |
| f | 0:102 | 1980-2023 | 0,9687 | −2,0158 | 0,007664 | −0,01545 | 2,9455 | −0,3395 | 6059 | 37532 | 0,06854 |
| f | 0:102 | 1980-2019 | 0,9906 | −2,2334 | 0,007802 | −0,01742 | 2,4617 | −0,5082 | 5126 | 33920 | 0,06798 |
| f | 50:102 | 1947-2023 | 0,9735 | −0,7222 | 0,022965 | −0,01659 | 1,3994 | −0,4330 | 7602 | 41569 | 0,07423 |
| f | 50:102 | 1947-2019 | 0,9712 | −0,7696 | 0,023276 | −0,01791 | 1,3166 | −0,4991 | 6793 | 38977 | 0,07414 |
| f | 50:102 | 1970-2023 | 0,9815 | −0,8247 | 0,021073 | −0,01738 | 1,2133 | −0,3380 | 4807 | 29095 | 0,07003 |
| f | 50:102 | 1970-2019 | 0,9932 | −0,9006 | 0,021390 | −0,01926 | 1,0456 | −0,4733 | 4059 | 26562 | 0,06965 |
| f | 50:102 | 1980-2023 | 0,9675 | −0,8134 | 0,018957 | −0,01542 | 1,2520 | −0,3551 | 3769 | 23702 | 0,06853 |
| f | 50:102 | 1980-2019 | 0,9903 | −0,9059 | 0,019218 | −0,01741 | 1,0528 | −0,5262 | 3082 | 21229 | 0,06799 |
| f | 60:102 | 1947-2023 | 0,9708 | −0,5559 | 0,029448 | −0,01637 | 1,1494 | −0,4405 | 6215 | 34361 | 0,07435 |
| f | 60:102 | 1947-2019 | 0,9685 | −0,5976 | 0,029664 | −0,01773 | 1,0849 | −0,5088 | 5515 | 32158 | 0,07433 |
| **f** | **60:102** | **1970-2023** | **0,9796** | **−0,6490** | **0,026609** | **−0,01727** | **0,9990** | **−0,3398** | **3989** | **24251** | **0,06991** |
| f | 60:102 | 1970-2019 | 0,9927 | −0,7183 | 0,026707 | −0,01918 | 0,8663 | −0,4798 | 3410 | 22168 | 0,06956 |
| f | 60:102 | 1980-2023 | 0,9643 | −0,6371 | 0,024012 | −0,01530 | 1,0255 | −0,3578 | 3067 | 19722 | 0,06836 |
| f | 60:102 | 1980-2019 | 0,9896 | −0,7225 | 0,023956 | −0,01731 | 0,8697 | −0,5351 | 2555 | 17704 | 0,06786 |

En **gras** : le cas de base retenu.

Lectures à porter au rapport :
- **Vérification de l'artefact d'échelle.** Période 1970-2023, les trois plages d'âges :
  d̂ = −2,1934 / −0,7974 / −0,6023 (hommes), soit un rapport 3,6 ; mais
  β₆₅·d̂ = −0,01815 / −0,01830 / −0,01823, **amplitude relative 0,9 %**. Chez les femmes :
  −0,01740 / −0,01738 / −0,01727, **0,8 %**. Les trois plages décrivent donc la **même**
  dynamique d'amélioration à 65 ans : ~1,8 %/an (H), ~1,7 %/an (F).
- **Mais cette neutralité dépend de la période** — amplitude de β₆₅·d̂ **entre** les trois
  plages d'âges :

  | début de période | hommes | femmes |
  |---|---|---|
  | 1947 (fin 2023) | **16,0 %** (−0,01537 / −0,01359 / −0,01312) | **8,4 %** (−0,01779 / −0,01659 / −0,01637) |
  | 1970 (fin 2023) | 0,9 % | 0,8 % |
  | 1980 (fin 2023) | 1,2 % | 1,0 % |

  La plage d'âges n'est neutre **qu'une fois l'après-guerre exclu**. Sur 1947-2023, elle
  déplace β₆₅·d̂ de 16 % chez les hommes — ce qui n'est pas négligeable pour une VAP. C'est
  un **argument supplémentaire en faveur d'un début en 1970** : quand les trois plages
  cessent de diverger, c'est que la structure β_x est devenue stable dans le temps, ce que
  Lee-Carter suppose précisément.
- **Déviance et BIC ne sont pas comparables entre variantes** : le jeu de données change
  (nobs va de 2 322 à 7 931 cellules). Ils ne discriminent rien ici.
- La **RMSE sur fenêtre commune** (60-102 × 1980-2019) est comparable, mais ne discrimine
  quasiment pas non plus : 0,130 à 0,133 (H), 0,068 à 0,074 (F), ~2 % d'écart entre les
  extrêmes. Le choix ne se joue **pas** sur la qualité d'ajustement in-sample.
- Le critère qui discrimine est la **linéarité de κ** : R² passe de 0,8994 (60:102,
  1947-2023) à 0,9799 (60:102, 1970-2023). Exclure l'après-guerre immédiat change
  réellement quelque chose ; restreindre encore à 1980 ne gagne plus rien (0,9676).
- **ρ₁(Δκ) systématiquement négatif** (−0,14 à −0,54) : ce n'est **pas** une violation de
  l'hypothèse de projection, c'est la signature attendue de l'erreur d'estimation sur κ̂ —
  voir TODO globaux.

### A.3.i — effet des années COVID (comparaison à plage d'âges constante, donc valide sur d̂)
Écart relatif de d̂ entre fin 2019 et fin 2023, par plage × début : **+6,15 % à +12,82 %**.
Le signe est toujours positif : inclure 2020-2021 **réduit systématiquement l'amplitude de
la dérive** (surmortalité ⇒ κ remonte ⇒ tendance baissière atténuée). L'effet croît quand
la période est courte (1980 : +9,7 à +12,8 %) et quand la plage d'âges est haute.

Pour le cas de base (60:102, début 1970), sur la quantité qui pilote réellement la VAP :

| | β₆₅·d̂ fin 2019 | β₆₅·d̂ fin 2023 | écart | (écart sur d̂ seul) |
|---|---|---|---|---|
| Hommes | −0,01995 | −0,01823 | **+8,64 %** | +10,52 % |
| Femmes | −0,01918 | −0,01727 | **+9,98 %** | +9,64 % |

**Arbitrage tranché : cas de base = 1970-2023** (données complètes). Inclure 2020-2021
réduit |β₆₅·d̂| de ~9 %, donc **réduit les VAP** : c'est le sens **anti-prudentiel** pour
un assureur de rentes, dont le risque est de sous-estimer l'amélioration future de la
mortalité. La variante **1970-2019 devient une sensibilité prudentielle**, à propager en
sections C et D avec son effet chiffré sur VAP et SCR (voir TODO section C).

### A.3.i — profils séquentiels à rebours (Denuit-Goderniaux / HR09 §3.10)
`figures/fig_A3_profil_kappa_sequentiel.png` — panneau de droite tracé en **β₆₅·d̂**, pas
en d̂ : superposer trois plages d'âges sur l'échelle de κ inviterait au contre-sens.

| série | R² max | atteint en t₀ | amplitude de β₆₅·d̂ sur t₀ = 1947-1990 | sur tout le profil |
|---|---|---|---|---|
| h 0:102 | 0,9839 | 1970 | 31,8 % (−0,01986 à −0,01445) | 31,4 % |
| h 50:102 | 0,9827 | 1969 | 49,1 % (−0,02138 à −0,01272) | 48,6 % |
| h 60:102 | 0,9804 | 1969 | 50,9 % (−0,02093 à −0,01211) | 50,6 % |
| f 0:102 | 0,9843 | 1967 | 19,3 % (−0,01965 à −0,01620) | 38,4 % |
| f 50:102 | 0,9831 | 1967 | 25,3 % (−0,02039 à −0,01580) | 41,4 % |
| f 60:102 | 0,9813 | 1967 | 26,4 % (−0,02029 à −0,01551) | 44,0 % |

- Le **maximum de R² tombe en t₀ = 1967-1970** pour les six séries : argument chiffré en
  faveur d'un début de calibration autour de 1970, indépendant du sexe et de la plage.
- **Recalcul sur la quantité invariante effectué** : les amplitudes relatives sont
  **inchangées** par rapport au calcul sur d̂ seul. C'est attendu — β₆₅ est une constante à
  l'intérieur d'un profil donné, donc la multiplication ne modifie pas une amplitude
  relative. L'écart entre sexes n'était donc **pas** un artefact d'échelle : il porte sur
  des séries comparées **à plage d'âges fixée**. Le constat subsiste après vérification :
  chez les hommes, restreindre aux t₀ ≤ 1990 ne réduit pas l'amplitude (50,9 % vs 50,6 %
  pour 60:102), alors que chez les femmes elle tombe de 44,0 % à 26,4 %.

### A.3.ii — paramètres du fit de travail (60:102 × 1970-2023)
Contraintes d'identification **vérifiées numériquement** : Σβ_x = 1 (H **et** F) ;
Σκ_t = −2,053913e−15 (H) et −1,94289e−15 (F).
Pente de la régression κ_t ~ t : **−0,6645 (H)**, **−0,7117 (F)** — κ décroît, conforme au
signal de la checklist. **β_x < 0 : 0 âge** chez les deux sexes — également conforme.
Figure `fig_A3_parametres.png`.

### A.3.iii — MLE (A.1) vs Lee-Carter
RMSE sur les log-taux : **0,1353 (H)** et **0,0768 (F)** sur la fenêtre de calibration
complète ; **0,0767 (H)** et **0,0614 (F)** sur la seule année 2022.
Cellules exclues (µ̂ = 0 ⇒ ln µ̂ = −∞) : 2 (H), 0 (F).
Figures `fig_A3_mle_vs_lc.png`, `fig_A3_mle_vs_lc_age65.png`.
*Lecture retenue* : profil périodique 2022 + série temporelle à l'âge 65. La diagonale de
cohorte n'a que **deux points observés** (65 ans en 2022, 66 ans en 2023) : une comparaison
cohorte sur données historiques est matériellement impossible.

### A.3.iv — résidus
Résidus NA (cellules w = 0) : **0** chez les deux sexes. Moyenne −0,0141 (H) / −0,0171 (F),
écart-type 0,9699 (H et F).
Écart-type des moyennes de résidus standardisées par √n, par groupe (**attendu ~1 en
l'absence de structure**) :

| regroupement | hommes | femmes |
|---|---|---|
| par âge x | 0,14 | 0,08 |
| par année t | 0,58 | 0,83 |
| **par cohorte c = t − x** | **2,98** | **2,76** |

🚩 **Effet de cohorte avéré, chiffré et non capturé par Lee-Carter.** Le modèle
α_x + β_x·κ_t ne contient aucun terme en t − x ; les regroupements par âge et par année
sont attendus **en dessous** de 1 (α_x et κ_t sont précisément estimés pour annuler ces
marges), et c'est bien le cas. Seule la cohorte est libre, et elle sort à ~3× le niveau
d'absence de structure. Les heatmaps `fig_A3_residus_H.png` / `fig_A3_residus_F.png` le
montrent visuellement. **Documenté, pas corrigé** : le projet ne demande que Lee-Carter.
Sens du biais indiqué par HR09 : en ajoutant un effet cohorte (APC), la mortalité projetée
devient plus **basse** que sous LC — donc les VAP de rentes seraient probablement
sous-estimées ici.

### A.3.v — projection centrale
`forecast(h = 50, kt.method = "mrwd", jumpchoice = "fit")`, horizon 2024-2073.
Dérive estimée de la marche aléatoire : **d̂ = −0,6023 (H)** et **−0,6490 (F)**.
Amorce à l'âge 65 : µ ajusté 2023 = 0,013424 (H) / 0,006706 (F) ; µ observé 2023 =
0,013327 (H) / 0,007090 (F) ; écart **−0,72 % (H)** et **+5,73 % (F)**.
⚠️ L'écart est presque huit fois plus grand chez les femmes : le choix `jumpchoice` n'est
pas symétrique entre les deux cohortes. Verrouillé à `"fit"` partout (A.3.v, A.3.vi **et
section B**) — voir TODO global.
Figures `fig_A3_projection_centrale.png`, `fig_A3_projection_cohorte.png`.

### A.3.vi — bootstrap semi-paramétrique, N = 5000
`bootstrap(type = "semiparametric", deathType = "observed")` puis
`simulate(nsim = 1, h = 50, kt.method = "mrwd", jumpchoice = "fit")` → 1 × 5000 = 5000 chemins.
Le tirage est D*_{x,t} ~ Poisson(ETR_{x,t}·µ̂_x(t)) où µ̂ est le **taux MLE brut de A.1** :
ETR·µ̂ = D, donc le bootstrap tire autour des **décès observés** (BDVK05 §4.2, `deathType = "observed"` ;
la variante `"fitted"` est celle de Renshaw-Haberman 2008 et n'est pas retenue).

**Coût réel : 5139,4 s = 85,7 min** pour 5000 réplications × 2 sexes (run de validation à
200 : 181,9 s). **Réplications à paramètres NA : 0 (H) et 0 (F).**

**Test de convention de clôture** (µ̂ **bruts** de A.1, x = 65…101, année 2023) :
e₆₅ = **18,396 (H)** et **21,572 (F)** — reproduit les 18,40 / 21,57 d'A.1. Clôture cohérente
entre les scripts 01 et 02. À titre d'information (**pas un test**), la même formule sur les
taux **ajustés** LC donne 18,455 (H) et 21,638 (F) : le lissage relève e₆₅ périodique de
~0,06 an.

**Espérance de vie de COHORTE à 65 ans** (cohorte des 65 ans en 2022, lecture diagonale
µ_{65+k}(2022+k), k = 0…36) :

| | prévision ponctuelle | moyenne bootstrap | IC 95 % | largeur relative |
|---|---|---|---|---|
| Hommes | 19,778 | 19,766 | [18,795 ; 20,738] | 9,83 % |
| Femmes | 23,156 | 23,138 | [22,068 ; 24,179] | 9,12 % |

Écart cohorte / périodique 2023 : +1,38 an (H) et +1,58 an (F) — c'est l'apport des gains
de mortalité futurs projetés. Figure `fig_A3_hist_e65_cohorte.png` (distributions unimodales,
quasi symétriques). IC sur les log-taux projetés : `fig_A3_bootstrap_ic.png` (cône
s'élargissant avec l'horizon, conforme).

**Les quatre diagnostics de contrôle**, confrontés aux étalons de `docs/BOOTSTRAP_NOTES.md` §9 :

| # | Obtenu | Étalon | Verdict |
|---|---|---|---|
| i — biais moyenne/ponctuel | **−0,06 % (H)**, **−0,07 % (F)** | 0,1 à 1,0 % (BDVK05) | ✅ L'étalon est un **plafond** (« au-delà du pourcent, cherche le bug »), pas une cible. Un biais dix fois plus petit est un bon résultat. À n = 200 on avait ±0,15 % : le passage à 5000 a absorbé l'erreur de Monte-Carlo, ce qui **valide le protocole en deux temps**. |
| ii — largeur relative IC 95 % de e₆₅ | **9,83 % (H)**, **9,12 % (F)** | quelques % à ~15 % (BDVK05, mais à **90 %** et sur a₆₅) | ✅ dans la fourchette |
| iii — asymétrie H/F | rapport **1,08** ; σ̂(Δκ) = **0,8924 (H)** vs **0,9990 (F)** | hommes systématiquement plus larges (BDVK05) | ✅ sens conforme, et l'ampleur s'explique par nos propres estimations : BDVK05 attribue l'écart aux IC plus larges sur la projection des κ masculins, or **ici c'est la volatilité féminine qui est supérieure** (0,9990 > 0,8924). Le ratio proche de 1 est donc cohérent avec les données autrichiennes, il n'a pas à ressembler au cas belge. |
| iv — dispersion paramétrique | voir tableau ci-dessous | β̂ nettement plus dispersé | ✅ confirmé, **mais seulement après normalisation** |

Diagnostic (iv), écart-type bootstrap moyen — **brut** puis **normalisé par l'amplitude
(max − min) du paramètre au fit central**, α, β et κ n'ayant pas la même échelle :

| | α̂ brut | α̂ norm. | β̂ brut | β̂ norm. | κ̂ brut | κ̂ norm. |
|---|---|---|---|---|---|---|
| Hommes | 0,01015 | 0,0026 | 0,00098 | **0,0306** | 0,21490 | 0,0066 |
| Femmes | 0,00699 | 0,0016 | 0,00062 | **0,0202** | 0,19189 | 0,0055 |

⚠️ **En écart-type brut, β̂ est le paramètre le MOINS dispersé** (0,00098 contre 0,215 pour κ̂) :
lu ainsi, le signal attendu semblerait inversé. C'est un artefact d'échelle — β̂ vit sur
[0,004 ; 0,035] quand κ̂ couvre ~33 unités. **Normalisé par l'amplitude, β̂ est 4,6× plus
dispersé que κ̂ et 12× plus qu'α̂** (hommes), ce qui confirme la hiérarchie annoncée par la
vignette StMoMo §8.
Dispersion de β̂ par âge (écart-type bootstrap) : 0,00054 (65 ans) → 0,00046 (85) →
0,00467 (101) → **0,00637 (102)**, soit ×12 entre 65 et 102 ans chez les hommes ; 0,00057 →
0,00032 → 0,00226 → 0,00298 chez les femmes. Conforme à l'attendu « plus dispersé aux âges
extrêmes, où l'exposition est faible ». Figure `fig_A3_bootstrap_parametres.png`.

**Constat lié — l'enveloppe bootstrap de β̂ traverse zéro aux grands âges.** Le fit central
donne β_x > 0 **partout** (minimum 0,00190 chez les hommes, 0,00374 chez les femmes), mais
le quantile 2,5 % bootstrap passe sous zéro :

| | âges concernés | nombre d'âges | fraction max de réplications à β̂ₓ < 0 |
|---|---|---|---|
| Hommes | **98 à 102** | 5 | **24,8 %** à l'âge 101 |
| Femmes | **101 à 102** | 2 | **7,5 %** à l'âge 102 |

Une fraction des réplications projette donc une mortalité **croissante** à ces âges (β_x < 0
avec κ décroissant). Sans effet matériel sur la VAP — un rentier de 65 ans a une
probabilité de l'ordre de 0,2 % d'atteindre 101 ans — mais à documenter : c'est la
conséquence directe de la faible exposition aux grands âges, déjà acceptée lors du choix
du plafond 102.

**Trois sources d'incertitude en A.3, une seule en B** (BDVK05 §4.2, *spécification gelée,
paramètres relâchés*) :

| | α, β | dérive & σ du RWD | aléa futur de κ |
|---|---|---|---|
| `simulate(LCfit, nsim, h)` — **section B** | fixés | **fixés** | simulé |
| `simulate(LCboot, h)` — **section A.3** | ré-estimés par réplication | **ré-estimés par réplication** | simulé |

La formulation « A.3 = incertitude sur (α, β) » est **fausse** et ne doit apparaître nulle part.

### Fichiers produits par le script 02

| Fichier | Taille | Versionné |
|---|---|---|
| `resultats/02_variantes.rds` | 0,01 Mo | oui |
| `resultats/02_fit_lc.rds` | 1,60 Mo | oui |
| `resultats/02_forecast_lc.rds` | 1,68 Mo | oui |
| `resultats/02_boot_ic.rds` | 0,19 Mo | oui — c'est ce que B.3 comparera |
| `resultats/LCboot_5000.rds` | 14,09 Mo | oui |
| `resultats/LCsimPU_5000.rds` | **339,87 Mo** | **non** — `.gitignore` : `resultats/LCsimPU_*.rds` |

Les artefacts du run de validation (`LCboot_200.rds`, `LCsimPU_200.rds`) restent sur disque
mais ne sont pas versionnés : ils sont reproductibles en ~3 min.

## Décisions de conception pour le script 03 (section B) — arrêtées avant codage

### Sensibilité prudentielle COVID : portée
La sensibilité 1970-2019 exige un **second jeu de 5 000 trajectoires**, issu du fit
`60:102 × 1970-2019` — **fit ré-ajusté dans le script 03** : les objets `fitStMoMo` ne sont
pas dans `resultats/02_variantes.rds`, qui ne contient que les tableaux de synthèse
(`grille`, `covid`, `profils`). Voir la sous-section « Section B — détail » pour le
contrôle de reproduction.

| | Décision | Raison |
|---|---|---|
| (a) | **Pas de second bootstrap.** | La sensibilité porte sur la **VAP centrale**, pas sur les intervalles. Re-bootstrapper coûterait 86 min pour une information qui n'est pas demandée. |
| (b) | Le script 03 produit **deux jeux de trajectoires** : cas de base (1970-2023) et variante (1970-2019). | Les deux alimentent C ; `simulate(LCfit, nsim = 5000, h = 50)` coûte quelques secondes, contrairement au bootstrap. |
| (c) | **Portée limitée à la section C** (VAP et primes). La section D **n'est pas dupliquée**. | La sensibilité exigée par le brief en D porte sur le **rendement des actifs** (`rdt_actifs` = 1 % au lieu de 4 %), pas sur la calibration de mortalité. |

### Interface officielle vers C et D : la diagonale de cohorte
`resultats/LCsimPU_5000.rds` pèse **340 Mo** et la section B produira un array de taille
comparable. **C et D n'ont pas besoin de l'array complet** : la seule quantité qu'ils
utilisent est la **diagonale de cohorte** µ_{65+k}(2022+k), k = 0…36 (âges 65 à 101,
années 2022 à 2058), par trajectoire et par sexe.

- **Format** : deux matrices **37 × 5000** (âges en lignes, **trajectoires en colonnes**,
  conformément à CLAUDE.md), une par sexe. Environ **3 Mo**, donc **versionnable**.
- k = 0…36 couvre exactement la convention de clôture d'A.1 (table jusqu'à ℓ₁₀₂, donc µ
  mobilisé jusqu'à l'âge 101) : ₖp₆₅ pour k = 1…37 pour la rente [A], k = 1…15 pour la
  rente [B], et µ₆₅(2022) pour les décès de la première année en section D.
- **Déclaré comme l'interface officielle** vers les sections C et D. L'array complet reste
  sur disque, **hors versionnement**, pour les seules figures d'IC de la section B.
- Même contrainte qu'en A.3 : `jumpchoice = "fit"` et `kt.method = "mrwd"`, sinon la
  comparaison B.3 avec les intervalles bootstrap devient ininterprétable.

## Section B — détail (chiffres recopiés de la console du run de `scripts/03_simulation_kappa.R`)

### Ce qui distingue B de A.3.vi
`simulate(LCfit, ...)` est appelé sur le **fit**, jamais sur un objet bootstrap : le script
ne contient **aucun** appel à `bootstrap()`. Preuve numérique que α, β et κ sont bien gelés
jusqu'en 2023 : `max|fitted − fit central| = 0` (H **et** F), et l'étendue de µ₆₅(2022) et
de µ₆₆(2023) sur les 5000 trajectoires vaut **0** et **0**.

### B.1 — diagonale de cohorte (`fig_B1_taux_cohorte.png`)
| | 65 ans (2022, ajusté) | 75 ans (2032) | 85 ans (2042) | 101 ans (2058) |
|---|---|---|---|---|
| Hommes | 0,014179 | 0,026586 | 0,078172 | 0,510320 |
| Femmes | 0,007104 | 0,014754 | 0,054585 | 0,448534 |

Taux **observés** d'A.1 sur les deux seuls points historiques de la diagonale :
µ₆₅(2022) = 0,014177 (H) et 0,007119 (F) ; µ₆₆(2023) = 0,014339 (H) et 0,006991 (F).
C'est ici que la lecture cohorte devient **matériellement possible** — 37 points, dont 35
projetés — alors qu'A.3.iii n'en avait que 2. Argument versé au TODO ouvert
« périodique vs diagonale », qui **reste ouvert**.

### B.2 — contrôles de la marche aléatoire
κ₂₀₂₃ = −14,7349 (H) et −15,2863 (F) ; d̂ = −0,6023 / −0,6490 ; σ̂ = 0,8924 / 0,9990.

| | h | moyenne obtenue | attendue (κ₂₀₂₃ + h·d̂) | sd obtenue | attendue (σ̂·√h) |
|---|---|---|---|---|---|
| H | 10 (2033) | −20,7433 | −20,7583 | 2,7872 | 2,8219 |
| H | 25 (2048) | −29,6854 | −29,7935 | 4,4481 | 4,4618 |
| H | 50 (2073) | −44,6048 | −44,8520 | 6,2391 | 6,3100 |
| F | 10 (2033) | −21,7906 | −21,7762 | 3,1485 | 3,1592 |
| F | 25 (2048) | −31,5305 | −31,5111 | 5,0874 | 4,9951 |
| F | 50 (2073) | −47,8787 | −47,7358 | 7,1962 | 7,0641 |

**Le cône s'élargit** : largeur de l'IC 95 % sur ln µ₆₅(t) en 2024 / 2048 / 2073 =
0,1057 / 0,5283 / 0,7422 (H) et 0,1050 / 0,5262 / 0,7455 (F) — `croissante : TRUE` pour
les deux sexes. Confirmé visuellement par `fig_B2_trajectoires_kappa.png`, où l'enveloppe
analytique κ₂₀₂₃ + h·d̂ ± 1,96·σ̂·√h épouse l'éventail simulé.

**Innovations standardisées** z = (Δκ* − d̂)/σ̂, n = 250 000 tirages par sexe :
moyenne **+0,00554 (t = +2,77)** et **sd = 0,99900** chez les hommes ; **−0,00286
(t = −1,43)** et **sd = 1,00226** chez les femmes. C'est ce contrôle qui tranche l'écart
de −44,6048 contre −44,8520 : la spécification du RWD est intacte (c'est sd(z) = 1 qui la
teste, et il est exact), l'écart sur la moyenne est du bruit de Monte-Carlo — voir TODO.

### B.3 — comparaison avec le bootstrap d'A.3.vi
| | moyenne | IC 95 % | largeur relative |
|---|---|---|---|
| **Hommes — B (κ seul)** | 19,762 | [18,823 ; 20,660] | **9,30 %** |
| Hommes — A.3.vi (bootstrap) | 19,766 | [18,795 ; 20,738] | 9,83 % |
| **Femmes — B (κ seul)** | 23,149 | [22,123 ; 24,162] | **8,81 %** |
| Femmes — A.3.vi (bootstrap) | 23,138 | [22,068 ; 24,179] | 9,12 % |

**Oracle satisfait** : A.3 est strictement plus large, rapport **1,058 (H)** et **1,036 (F)**
— très loin du facteur 2 qui aurait imposé un diagnostic. Sur les log-taux en 2058, le
rapport A.3/B vaut 1,06 / 1,06 / 1,05 aux âges 65 / 75 / 85 (H) et 1,05 / 1,05 / 1,04 (F).
Figures `fig_B3_comparaison_ic.png` (l'enveloppe A.3 contient l'enveloppe B partout, la
différence n'est qu'une frange) et `fig_B3_e65_comparaison.png`.

### B.3 — POURQUOI l'élargissement est de second ordre (le point central de la section)
Décomposition de Var(κ₂₀₇₃) sous A.3, à partir des paramètres bootstrap **déjà calculés**
par le script 02 (aucun re-bootstrap), avec
κ*_{2023+h} = κ*₂₀₂₃ + h·d̂* + σ̂*·Σz* :

| | amorce | dérive | covariance | projection | sd A.3 | sd B | rapport |
|---|---|---|---|---|---|---|---|
| Hommes | 0,052 | 0,142 | 0,146 | **43,477** | 6,619 | 6,310 | **1,049** |
| Femmes | 0,044 | 0,085 | 0,098 | **53,476** | 7,328 | 7,064 | **1,037** |

Le rapport **prédit analytiquement** (1,049 / 1,037) reproduit le rapport **mesuré** sur
e₆₅ (1,058 / 1,036) : le mécanisme est compris, pas seulement constaté.

⚠️ **Le canal dominant n'est PAS la dérive** (0,142 + 0,146 = 0,29 sur 43,8, soit 0,7 %)
mais **σ̂\* > σ̂** : σ̂* moyen = 0,93156 contre σ̂ = 0,89237 (**+4,4 %**) chez les hommes,
1,03326 contre 0,99902 (**+3,4 %**) chez les femmes. Le bruit d'estimation sur κ̂ ajoute un
MA(1) qui gonfle sd(Δκ*) — **même mécanisme que le ρ₁(Δκ) négatif** documenté en A.3.i.

**Ce que NI A.3 NI B ne capturent** : `sd(d̂*)` bootstrap = **0,00755** contre
σ̂/√n = **0,12258** (H, n = 53 accroissements), soit un **facteur 16** ; 0,00583 contre
0,13723 chez les femmes, **facteur 24**. Le bootstrap de Poisson **conditionne sur le
chemin κ̂ observé** et ne perturbe que les décès : il ne peut pas produire une autre
réalisation de la marche aléatoire, donc pas l'incertitude d'échantillonnage **temporelle**
de la dérive.

Chiffrage de ce qui manque, au sens de Lee & Carter (1992, annexe B), où
**Var[κ_{T+h}] = h·σ̂² + h²·σ̂²/n** *(évalué sur console à part, pas produit par le script 03
— il s'agit d'un élément de discussion du rapport, pas d'un résultat de la section)* :

| | h·σ̂² | h²·σ̂²/n | sd totale | sd sous B | élargissement |
|---|---|---|---|---|---|
| Hommes (h = 50, n = 53) | 39,816 | 37,562 | 8,7965 | 6,3100 | **+39,4 %** |
| Femmes (h = 50, n = 53) | 49,902 | 47,077 | 9,8478 | 7,0641 | **+39,4 %** |

Les deux termes sont du **même ordre** : à h = 50 l'incertitude sur la dérive pèse presque
autant que l'aléa de projection lui-même. L'élargissement vaut en forme fermée
**√(1 + h/n) = 1,39406**, indépendant de σ̂ — d'où un chiffre **identique pour les deux
sexes** : c'est un résultat structurel (rapport horizon / longueur de calibration), pas une
coïncidence.

⚠️ **Portée exacte de ce +39,4 %, à ne pas outrepasser.** Il s'applique à **sd(κ₂₀₇₃)**,
c'est-à-dire à l'écart-type de κ **au seul horizon h = 50**. e₆₅ et la VAP ne sont pas des
fonctions de κ₂₀₇₃ : elles intègrent **tous** les horizons k = 1…37, dont les horizons
courts, où le facteur √(1 + k/n) est nettement plus faible (à k = 5, il vaut 1,05). L'effet
sur leurs intervalles serait donc **moindre**, et **différent entre les deux sexes** — la
neutralité en σ̂ de la forme fermée ne survit pas à l'agrégation sur les horizons, puisque
la pondération par ₖp₆₅ diffère entre H et F. Même ordre de grandeur, pas égalité : tout
chiffre précis sur e₆₅ ou la VAP demanderait un calcul dédié, non effectué ici.

**Seconde réserve** : Var[κ_{T+h}] = h·σ̂² + h²·σ̂²/n est l'**approximation** de Lee &
Carter (1992, annexe B). Elle traite σ̂ comme connue et néglige donc l'incertitude sur σ̂
lui-même. C'est acceptable ici (n = 53 accroissements), mais à **nommer comme une
approximation**, pas à présenter comme la variance exacte.

**Conclusion à écrire au rapport** : ni le bootstrap de Poisson ni la simulation de κ seul
ne capturent l'incertitude d'échantillonnage temporel de la dérive. C'est la **limite
commune aux deux jeux d'intervalles**. Le faible écart entre A.3 et B (+3 à +6 %) ne dit
donc pas que l'incertitude d'estimation est négligeable — il dit que **la source
d'incertitude d'estimation qui compte n'est pas celle que le bootstrap échantillonne**.

### Sensibilité prudentielle COVID (alimente la section C)
Le fit `60:102 × 1970-2019` est **ré-ajusté** dans le script 03 (absent de
`02_variantes.rds`, qui ne contient que les tableaux de synthèse). Contrôle de reproduction
contre `02_variantes.rds$grille` :

| | d̂ obtenu | σ̂ obtenu | β₆₅ obtenu | écarts au tableau du script 02 |
|---|---|---|---|---|
| Hommes | −0,6731840 | 0,7295808 | 0,0296397 | 4,7e-10 / 3,4e-10 / 2,3e-12 |
| Femmes | −0,7182583 | 0,8662930 | 0,0267074 | 4,0e-11 / 1,0e-08 / 3,0e-11 |

`tous < 1e-6 : TRUE`. Pas de matrice de poids dans ce refit : la fenêtre
60:102 × 1970-2019 ne contient **aucune** cellule ETR = 0 (vérifié : 0 cellule, ETR minimum
1,49 chez les hommes et 5,24 chez les femmes), donc w_{x,t} y vaut identiquement 1.

**Test RNG imprimé** : `identical(.Random.seed) = FALSE` — `fit()` **déplace** le flux
(gnm tire ses valeurs de départ). D'où le `set.seed(1234)` avant **chaque** couple de
`simulate()`, cas de base inclus.

e₆₅ de cohorte sous la variante : **20,621** (cas de base 19,762) chez les hommes,
**23,906** (23,149) chez les femmes ; largeurs relatives 9,05 % et 8,67 %.
Amorce µ₆₅(2022) médiane : 0,014179 → **0,012477** (H) et 0,007104 → **0,006187** (F) — la
variante ne se contente pas d'une dérive plus forte, elle part aussi d'un niveau plus bas
(voir TODO).

### Fichiers produits par le script 03
| Fichier | Taille | Versionné |
|---|---|---|
| `resultats/03_diagonales_cohorte.rds` | 5,26 Mo | **oui — interface officielle vers C et D** |
| `resultats/03_ic_trajectoires.rds` | 3,89 Mo | oui (quantiles, κ, e₆₅ et largeurs des deux scénarios) |
| `resultats/LCsimK_h_5000.rds` | 80,67 Mo | **non** — `.gitignore` : `resultats/LCsimK_*.rds` |
| `resultats/LCsimK_f_5000.rds` | 80,69 Mo | **non** |

`K` = incertitude de κ **seul** (section B) ; `PU` = *parameter uncertainty* (A.3.vi).
Les sections C et D ne chargeront **que** `03_diagonales_cohorte.rds`.

## Section C — détail (chiffres recopiés de la console du run de `04_pricing_vap.R`)

### Ce qui distingue C des sections amont
**Aucun appel à StMoMo, aucun `library()`, aucun `set.seed`.** L'absence de graine est un
**signal**, pas un oubli : le script ne tire rien, il consomme les 5000 trajectoires du
script 03. Deux `readRDS` seulement — `03_diagonales_cohorte.rds` (tarification) et
`03_ic_trajectoires.rds` (**bloc oracle uniquement**, aucune quantité tarifaire n'en dérive).

Convention verrouillée : rentes **à terme échu**, a₆₅ = Σ_{k≥1} vᵏ·ₖp₆₅. Le premier µ
consommé est µ₆₅(2022) ; il produit ₁p₆₅ = survie jusqu'à 66 ans, qui porte le paiement
k = 1 versé **fin 2023**. Aucun terme k = 0, ni dans `S` (que `cumprod` ne produit pas) ni
dans `v^(1:n)` (qui démarre à v¹).

### C.1 / C.2 — VAP à 3 %, moyenne et variance sur les 5000 trajectoires
| | moyenne | variance | écart-type | quantiles 2,5 / 50 / 97,5 % |
|---|---|---|---|---|
| **[A] Hommes** | **13,7899** | 0,061949 | 0,248895 | 13,2877 / 13,7939 / 14,2566 |
| **[A] Femmes** | **15,6995** | 0,067253 | 0,259332 | 15,1843 / 15,7053 / 16,1995 |
| [B] Hommes | 10,4086 | 0,003615 | 0,060125 | 10,2838 / 10,4107 / 10,5193 |
| [B] Femmes | 11,1011 | 0,001602 | 0,040024 | 11,0196 / 11,1018 / 11,1750 |

⚠️ **La variance de [B] s'inverse entre sexes** : 0,003615 (H) > 0,001602 (F), alors que [A]
ordonne F > H. Ce n'est pas un bug — explication au TODO dédié.

### C.3 — VAP moyenne selon le taux technique
| taux | [A] H | [A] F | [B] H | [B] F |
|---|---|---|---|---|
| 0 % *(contrôle, pas un scénario)* | 19,2622 | 22,6487 | 12,9161 | 13,8542 |
| 1 % | 17,1143 | 19,8904 | 11,9896 | 12,8355 |
| 2 % | 15,3124 | 17,6060 | 11,1576 | 11,9221 |
| **3 % (brief)** | **13,7899** | **15,6995** | **10,4086** | **11,1011** |
| 4 % | 12,4945 | 14,0964 | 9,7327 | 10,3614 |
| 5 % | 11,3851 | 12,7388 | 9,1213 | 9,6934 |

**Effet duration** (variation 1 % → 5 %) : **−33,48 % sur [A] contre −23,92 % sur [B]** chez
les hommes, −35,96 % contre −24,48 % chez les femmes.

### C.4 — primes uniques par le principe d'équivalence, t = 3 %
Prime = E[VAP] sur les 5000 trajectoires, pour une **rente annuelle de 1** : **13,7899 (H)**
et **15,6995 (F)** pour [A] ; 10,4086 et 11,1011 pour [B]. L'agrégation portefeuille
(500 H + 500 F) est une hypothèse de la **section D** et n'est pas faite ici.

### Sensibilité prudentielle COVID (calibration 1970-2019, 3 % seulement)
| | base | COVID | écart | sd(base) | sd(COVID − base) | rapport |
|---|---|---|---|---|---|---|
| [A] Hommes | 13,7899 | 14,2693 | **+3,48 %** | 0,24890 | 0,03646 | **0,15** |
| [A] Femmes | 15,6995 | 16,0948 | **+2,52 %** | 0,25933 | 0,03725 | **0,14** |
| [B] Hommes | 10,4086 | 10,5749 | +1,60 % | 0,06012 | 0,01959 | 0,33 |
| [B] Femmes | 11,1011 | 11,2001 | +0,89 % | 0,04002 | 0,01332 | 0,33 |

L'effondrement de `sd(ΔVAP)` (rapport 0,14 à 0,33) est la **preuve chiffrée de
l'appariement** par nombres aléatoires communs : l'écart mesuré est un effet de
**calibration pur**, sans bruit de Monte-Carlo différentiel. La variante prudentielle
renchérit les rentes, davantage sur la viagère que sur la temporaire.

### Oracles — tous satisfaits
| Oracle | Résultat |
|---|---|
| Alignement `S["1", ] = exp(−µ["65", ])` sur les 5000 colonnes | écart **0** (H et F) — preuve que `cumprod` n'a franchi aucune frontière de colonne |
| Intégrité `S` ∈ ]0,1], décroissante en k, 0 NA | TRUE / TRUE / 0 ; rownames `"1".."37"`, `attr age_atteint = 66:102` |
| 1 — VAP[B] < VAP[A] **par trajectoire** | 0 violation / 5000, H et F |
| 2 — bande recalibrée à terme échu | 13,7899 ∈ [13 ; 14] · 15,6995 ∈ [14 ; 16] |
| 3a — ancrage inter-scripts à taux nul | `max\|VAP₀ − (e₆₅ − ½)\|` = **0** (tolérance 1e-12) ; e₆₅ reconstitué **19,7622 / 23,1487** contre 19,762 / 23,149 publiés en B.3 |
| 3b — boucle scalaire naïve vs matriciel | écart max **2,13e-14** sur les 8 cas (2 sexes × 2 produits × taux 0 et 3 %) ; trajectoire 1 : 14,140819 des deux côtés |
| 4 — monotonie en t sur {0…5} % | TRUE pour les 4 séries |
| 5 — duration [A] > [B] | −33,48 / −23,92 (H) · −35,96 / −24,48 (F) |
| 6 — ordre H/F, moyenne **et** quantiles | F > H partout, [A] et [B] |
| Dernier terme inclus v³⁷·₃₇p₆₅ | 0,00231 (H) = **0,017 %** de VAP[A] · 0,00745 (F) = **0,047 %** |

⚠️ **Portée de l'oracle 3a, à savoir énoncer.** L'écart vaut exactement 0 parce que le script
03 a produit ses e₆₅ par `apply(diagonale, 2, esperance_vie_cohorte)` sur la matrice même que
C recharge : mêmes données, même formule, même ordre d'opérations flottantes. L'oracle teste
donc la **cohérence de convention entre 03 et 04**, et ne pourrait pas détecter une erreur
d'alignement **partagée**. La preuve d'alignement absolu vient d'ailleurs : du contrôle
`S["1", ] = exp(−µ["65", ])`, et des e₆₅ de B déjà validés contre A.1. C'est aussi la raison
d'être de l'oracle 3b, dont la boucle scalaire ne partage **aucune** primitive avec
`vap_rente()`.

### Fichiers produits par le script 04
| Fichier | Taille | Contenu |
|---|---|---|
| `resultats/04_vap.rds` | 3,09 Mo | `$base` = {vap (4 × 5000), survies (2 × 37 × 5000), primes, sensibilite} · `$covid` = {vap, ecart_relatif} |

Séparation **structurelle** base / covid : la section D charge `$base` et rien d'autre.
`$base$sensibilite` ne contient que des **moyennes** — aucune VAP par trajectoire à un autre
taux que 3 % n'est sauvée, le BE de D s'actualisant à `taux_actu_be = 2 %` depuis
`$base$survies`.

## TODO globaux

- **A.3.i — calibration retenue : `60:102 × 1970-2023`** (cas de base). Justifications en
  main : plafond 102 aligné sur `age_max_table` d'A.1 ; maximum du profil R²(t₀) en
  1967-1970 pour les six séries ; déviance, BIC et RMSE ne discriminent pas ; et le choix
  de plage d'âges devient quasi neutre sur ce qui pilote la VAP **une fois 1970 retenu**
  (β₆₅·d̂ varie de 0,8-0,9 % entre les trois plages sur 1970-2023, contre **16,0 % (H)** et
  **8,4 % (F)** sur 1947-2023). Cette neutralité est elle-même un argument pour 1970 : elle
  signale que la structure β_x s'est stabilisée, hypothèse centrale de Lee-Carter. Reste à
  **rédiger** la justification.
- **A.3.i — d̂ n'est comparable qu'à plage d'âges fixée.** Σβ_x = 1 fixe l'échelle de κ en
  fonction du nombre d'âges : comparer les d̂ de 0:102 et 60:102 est un contre-sens. Toujours
  raisonner sur **β₆₅·d̂**. Vaut aussi pour toute comparaison future en sections B, C et D.
- **A.3.i — instabilité de β₆₅·d̂ selon t₀, asymétrique entre sexes** (constat **vérifié sur
  la quantité invariante**, il ne s'agit pas d'un artefact d'échelle : la comparaison est
  faite à plage d'âges fixée). Chez les hommes, restreindre aux t₀ ≤ 1990 ne réduit pas
  l'amplitude (50,9 % vs 50,6 %) ; chez les femmes elle tombe de 44,0 % à 26,4 %. Le choix
  de t₀ pèse donc plus lourd pour les hommes — à mentionner, sans en tirer de correction.
- **A.3.i — ρ₁(Δκ) négatif : argument de défense, pas une limite.** On n'observe pas κ_t
  mais son estimateur κ̂_t = κ_t + η_t. Donc Δκ̂_t = Δκ_t + η_t − η_{t−1} : même si les vrais
  Δκ sont i.i.d., l'erreur d'estimation ajoute un **MA(1) à coefficient négatif**, dont
  l'autocorrélation d'ordre 1 est bornée par **−0,5**. Nos valeurs (−0,14 à −0,54) sont
  donc la **signature attendue** du bruit d'estimation, et non une violation de
  l'hypothèse de marche aléatoire. **Corollaire à écrire au rapport** : σ̂ = sd(Δκ̂) est
  gonflé par ce bruit, donc les intervalles de projection sont plutôt **conservateurs** —
  ce qui va dans le sens de la prudence pour le SCR de la section D.
- **A.3.ii — β̂ₓ change de signe dans le bootstrap aux grands âges.** Fit central positif
  partout, mais l'enveloppe bootstrap traverse zéro aux âges 98-102 (H, jusqu'à 24,8 % des
  réplications à l'âge 101) et 101-102 (F, jusqu'à 7,5 %) : une fraction des réplications
  y projette une mortalité croissante. Sans effet matériel sur la VAP, à documenter.
- **A.3.v/vi + B — `jumpchoice`, arbitrage GLOBAL : APPLIQUÉ ET VÉRIFIÉ.** `"fit"` est
  employé dans les **quatre** appels : `forecast` (A.3.v), `simulate(LCboot)` (A.3.vi), et
  **les deux couples de `simulate()` de la section B** (cas de base et variante 1970-2019).
  La comparaison d'IC de B.3 est donc bien à amorce identique — elle ne mélange pas
  incertitude et décalage de départ. **Ce qui reste ouvert** : la justification de `"fit"`
  **contre** `"actual"` comme choix global, à documenter au rapport. Écart ajusté/observé à
  65 ans en 2023 : −0,72 % (H) mais **+5,73 % (F)** — le choix n'est pas symétrique entre
  les deux cohortes.
- **A.3.iv — effet de cohorte avéré et chiffré** : écart-type des moyennes de résidus
  standardisées = 2,98 (H) / 2,76 (F) par cohorte, contre 0,14/0,08 par âge et 0,58/0,83
  par année. Lee-Carter n'a aucun terme en t − x. Documenté, **pas corrigé** (le projet ne
  demande que Lee-Carter). Sens du biais d'après HR09 : un modèle APC projetterait une
  mortalité plus basse ⇒ les VAP calculées ici sont probablement sous-estimées.
- **A.3.vi — terminologie à neutraliser au rapport** : « bootstrap semi-paramétrique au sens
  de StMoMo, désigné *Poisson bootstrap* par Brouhns, Denuit & Van Keilegom (2005) ». La
  vraie alternative discutée par le papier est le bootstrap **paramétrique** (BDV 2002b,
  tirage dans la normale multivariée asymptotique), **pas** le bootstrap résiduel.
- **A.3.vi — les quatre diagnostics sont conformes** ; deux méritent une phrase
  d'explication au rapport plutôt qu'un simple « conforme ».
  **(i) Biais de −0,06 % (H) / −0,07 % (F)** : l'étalon BDVK05 de 0,1-1,0 % est un
  **plafond** (« au-delà du pourcent, cherche le bug »), pas une cible ; un biais dix fois
  plus petit est un bon résultat. À n = 200 on mesurait ±0,15 % : le passage à 5000 a
  absorbé l'erreur de Monte-Carlo, ce qui **valide le protocole en deux temps**.
  **(iii) Asymétrie H/F de 1,08** : BDVK05 attribue l'écart aux IC plus larges sur la
  projection des κ masculins ; **ici c'est la volatilité féminine qui est supérieure**
  (σ̂(Δκ) = 0,9990 pour les femmes contre 0,8924 pour les hommes). Un ratio proche de 1 est
  donc **cohérent avec les données autrichiennes** — il n'a aucune raison de reproduire le
  cas belge.
- **B — année 1 sans risque systématique (à énoncer en défense, pour la section D).**
  Sous `jumpchoice = "fit"`, α, β et κ sont fixés jusqu'en 2023 : les taux 2022-2023 sont
  **identiques sur les 5000 trajectoires** (vérifié : `max|fitted − fit central| = 0`,
  étendue de µ₆₅(2022) et µ₆₆(2023) = 0). La mortalité de l'année 1 (2022 → 2023) ne porte
  donc **que** le risque diversifiable de la binomiale. Le risque de longévité systématique
  entre dans SCR₂₀₂₂ par la **réévaluation de BE₂₀₂₃**, pas par le nombre de décès.
- **B.3 — la limite commune aux deux jeux d'intervalles, quantifiée.** Ni le bootstrap de
  Poisson (qui conditionne sur le chemin κ̂ observé) ni la simulation de κ seul ne
  capturent l'incertitude d'échantillonnage **temporel** de la dérive :
  `sd(d̂*) = 0,00755` contre `σ̂/√n = 0,12258` (facteur 16 chez les hommes, 24 chez les
  femmes). Avec le terme de Lee & Carter (1992, annexe B), les deux jeux d'IC seraient
  **~39 % plus larges** à h = 50 — élargissement en forme fermée √(1 + h/n), identique
  pour les deux sexes. À rédiger : le faible écart A.3/B (+3 à +6 %) ne signifie pas que
  l'incertitude d'estimation est négligeable, mais que **celle qui compte n'est pas celle
  que le bootstrap échantillonne**.
- **B.2 — bruit de Monte-Carlo de la graine 1234, assumé et non corrigé.** La moyenne des
  250 000 innovations standardisées vaut **+0,00554 (t = 2,77)** chez les hommes, d'où
  l'écart de la moyenne de κ₂₀₇₃ (−44,6048 au lieu de −44,8520). La spécification du RWD
  est intacte — `sd(z) = 0,99900`, et c'est cette statistique-là qui la teste. Changer de
  graine pour effacer l'écart serait du **seed-hacking** : on ne le fait pas. Sous nombres
  aléatoires communs, ce bruit **s'annule** dans la sensibilité COVID de la section C,
  puisque les deux scénarios consomment le même flux.
- **B — nombres aléatoires communs base / variante COVID.** Les deux jeux de trajectoires
  sont **positivement corrélés par construction** : c'est le but (la variance de la
  *différence* de VAP s'effondre). Corollaire à savoir énoncer : ces 10 000 trajectoires
  ne forment **pas** un échantillon de taille 10 000 et ne peuvent jamais être fusionnées.
  L'appariement se fait par **indice de pas**, pas par année civile (base : pas 1 = 2024
  depuis κ₂₀₂₃ ; variante : pas 1 = 2020 depuis κ₂₀₁₉), soit un décalage de 4 indices ;
  les sommes cumulées partagent malgré tout ~33 termes sur 36.
- **B — la variante 1970-2019 projette ses années 2020-2023 au lieu de les ajuster.** Son
  amorce diffère donc du cas de base : µ₆₅(2022) médian passe de 0,014179 à 0,012477 (H)
  et de 0,007104 à 0,006187 (F). La sensibilité VAP de la section C mêlera l'effet
  « dérive plus forte » et l'effet « amorce sans surmortalité COVID ». À **énoncer** au
  rapport ; **ne pas** corriger par un recalage d'amorce (méthode non demandée).
- **B — `02_variantes.rds` ne contient pas les fits de la grille**, seulement `grille`,
  `covid` et `profils`. Le script 03 ré-ajuste donc la variante 1970-2019. Argument de
  défense à énoncer correctement : ce n'est **pas** que `fit()` ne tirerait rien — gnm tire
  ses valeurs de départ, et le test `identical(.Random.seed) = FALSE` le prouve — c'est que
  **l'optimum est unique** à contrainte d'identification près et que l'IRLS y converge quel
  que soit le départ, d'où une reproduction à la **tolérance de convergence** (écarts
  4,7e-10 à 1,0e-08), pas au bit près. La correction propre (sauver `fits_grille` dans le
  script 02) a été écartée pour ne pas re-bootstrapper 86 min sur un script déjà validé.
- **C.2 — ce que la variance des VAP mesure, et surtout ce qu'elle NE mesure PAS.** Les 5000
  VAP ne diffèrent que par la trajectoire de κ : c'est l'incertitude **systématique de
  longévité**, **non diversifiable** — 1000 assurés ne la réduisent pas. Le risque
  **idiosyncratique** des têtes individuelles n'y est pas ; il n'apparaîtra qu'en section D,
  par la binomiale des décès de l'année 1. Distinction centrale de la Q4 d'examen, et c'est
  ce risque non diversifiable qui fonde le SCR.
- **C.2 — l'inversion de la variance sur [B], à expliquer et non à corriger.**
  var(VAP[B]) = **0,003615 (H) contre 0,001602 (F)**, soit l'ordre **inverse** de [A]
  (0,061949 H contre 0,067253 F). Explication : la VAP[B] féminine est presque **saturée** —
  11,1011 sur **11,938**, l'annuité certaine 15 ans à 3 % (Σ_{k=1}^{15} 1,03^{−k}), soit
  **93,0 %** contre 87,2 % chez les hommes. Les ₖp₆₅ féminins sont si proches de 1 sur
  15 ans qu'ils laissent peu de prise aux chocs de κ : la sensibilité de exp(−µ) est
  proportionnelle au niveau de µ, et l'avantage de mortalité féminin **réduit** donc aussi la
  dispersion. Cohérent avec le rapport d'appariement COVID, **0,33 sur [B] contre 0,14-0,15
  sur [A]**. Question de défense probable.
- **C.2 — autres sources d'incertitude** (à lister au rapport, non chiffrées en C) :
  idiosyncratique (binomiale des décès) ; **modèle** (effet cohorte non capturé, cf. A.3.iv —
  d'après HR09 un modèle APC projetterait une mortalité plus basse, donc les VAP d'ici sont
  probablement **sous-estimées**) ; **estimation de la dérive** (cf. B.3 : ni A.3 ni B ne
  capturent l'incertitude d'échantillonnage temporelle de d̂, facteur 16 à 24 sur sd(d̂)) ;
  taux d'intérêt (traité en scénarios par C.3, pas en stochastique) ; frais.
- **C.2/B.3 — ne pas transporter le +39,4 % sur la VAP.** Ce chiffre porte sur sd(κ₂₀₇₃), au
  seul horizon h = 50. La VAP intègre k = 1…37, dont des horizons courts où √(1+k/n) est bien
  plus faible (1,05 à k = 5). Aucun chiffre précis sur la VAP ne s'en déduit.
- **C.3 — le taux technique n'est pas stochastique.** C.3 est une étude de **scénarios
  déterministes**, pas une modélisation du risque de taux. Le taux 0 % de la grille est un
  contrôle (oracle 3), pas un scénario tarifaire.
- **C.4 — principe d'équivalence, convention retenue.** Prime unique = E[VAP] sur les 5000
  trajectoires, pour une rente annuelle de 1 : c'est l'espérance **sous mortalité
  stochastique**. Prime **pure** — ni chargement, ni marge de prudence. À opposer au rapport
  à un chargement pour risque systématique, qui serait un **quantile** et non une moyenne.
- **C — clôture à ℓ₁₀₂ : lecture correcte du dernier terme.** [A] est bornée à 37 termes. La
  quantité imprimée v³⁷·₃₇p₆₅ (0,00231 H = 0,017 % de VAP[A] ; 0,00745 F = 0,047 %) est le
  **dernier terme INCLUS**, **pas** l'erreur de troncature. La queue omise est
  Σ_{k≥38} vᵏ·ₖp₆₅, que la clôture annule par convention (ℓ₁₀₃ = 0) ; le dernier terme en
  fournit un **majorant d'ordre de grandeur** — l'espérance de vie résiduelle à 102 ans valant
  ≈ 1,5-2 ans, la queue pèse grossièrement 1,5 à 2 fois ce terme, donc reste négligeable.
  **Ne jamais écrire « erreur de troncature = v³⁷·₃₇p₆₅ ».**
- **C — figure de rapport non produite par le script** : distribution des différences ΔVAP
  **par trajectoire** entre base et COVID. C'est là que les nombres aléatoires communs se
  voient (variance de la différence effondrée). Vecteurs disponibles dans
  `resultats/04_vap.rds$covid$vap`, sans relancer le script 04.
- **Section D — précision path-wise (à savoir énoncer, rien à coder).** `simulate()` de
  StMoMo produit de **vraies marches aléatoires** (algorithme 2 de HR09 §4.7 :
  κ*_{t+j} = κ_t + j·θ̂ + σ̂·Σᵢ z*ᵢ, z*ᵢ i.i.d.), et non des droites déviées par un tirage
  unique (algorithme 1 — marginalement identique, mais trajectoires lisses et ordonnées, et
  IC systématiquement plus larges). Le SCR = VaR₉₉,₅ % de −ΔBOF dépend de la queue du
  **chemin complet** (décès de l'année 1 → état 2023 → BE₂₀₂₃), pas de la marginale à
  l'horizon final : un éventail lisse fausserait la dépendance temporelle, donc le SCR.
- **Limites à énoncer en défense** (`docs/BOOTSTRAP_NOTES.md` §11) : expositions traitées
  comme fixes ; hypothèse de Poisson non testée formellement ; risque de modèle non couvert ;
  extrapolation pure du passé ; N = 200 insuffisant pour une VaR 99,5 % (d'où le run à 5000).
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
