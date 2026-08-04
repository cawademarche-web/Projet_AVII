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
| B (5000 trajectoires κ, IC, comparaison avec bootstrap) | `scripts/03_simulation_kappa.R` | — | — | — |
| C (VAP rentes [A]/[B], primes) | `scripts/04_pricing_vap.R` | — | — | **Sensibilité prudentielle COVID** : rejouer VAP et primes sur la calibration 1970-2019 (β₆₅·d̂ = −0,01995 H / −0,01918 F contre −0,01823 / −0,01727 en cas de base) et chiffrer l'écart ; à propager en D sur le SCR |
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
  β₆₅·d̂ = −0,01815 / −0,01830 / −0,01823, **amplitude relative 0,86 %**. Chez les femmes :
  −0,01740 / −0,01738 / −0,01727, **0,78 %**. Les trois plages décrivent donc la **même**
  dynamique d'amélioration à 65 ans : ~1,8 %/an chez les hommes, ~1,7 %/an chez les femmes.
  Le choix de plage d'âges ne change pratiquement rien à ce qui pilote la VAP.
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

## TODO globaux

- **A.3.i — calibration retenue : `60:102 × 1970-2023`** (cas de base). Justifications en
  main : plafond 102 aligné sur `age_max_table` d'A.1 ; maximum du profil R²(t₀) en
  1967-1970 pour les six séries ; déviance, BIC et RMSE ne discriminent pas ; et le choix
  de plage d'âges est de toute façon quasi neutre sur ce qui pilote la VAP (β₆₅·d̂ varie de
  moins de 1 % entre les trois plages). Reste à **rédiger** la justification.
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
- **A.3.v/vi + B — `jumpchoice`, arbitrage GLOBAL.** Verrouillé à `"fit"` dans les trois
  appels (`forecast` de A.3.v, `simulate(LCboot)` de A.3.vi, et **la section B devra faire
  de même**). Si A.3 et B partaient d'amorces différentes, la comparaison d'IC demandée en
  B.3 mélangerait incertitude et décalage de départ, et deviendrait ininterprétable. Écart
  ajusté/observé à 65 ans en 2023 : −0,72 % (H) mais **+5,73 % (F)** — le choix n'est pas
  symétrique entre les deux cohortes. À appliquer partout ou nulle part.
- **A.3.iv — effet de cohorte avéré et chiffré** : écart-type des moyennes de résidus
  standardisées = 2,98 (H) / 2,76 (F) par cohorte, contre 0,14/0,08 par âge et 0,58/0,83
  par année. Lee-Carter n'a aucun terme en t − x. Documenté, **pas corrigé** (le projet ne
  demande que Lee-Carter). Sens du biais d'après HR09 : un modèle APC projetterait une
  mortalité plus basse ⇒ les VAP calculées ici sont probablement sous-estimées.
- **A.3.vi — terminologie à neutraliser au rapport** : « bootstrap semi-paramétrique au sens
  de StMoMo, désigné *Poisson bootstrap* par Brouhns, Denuit & Van Keilegom (2005) ». La
  vraie alternative discutée par le papier est le bootstrap **paramétrique** (BDV 2002b,
  tirage dans la normale multivariée asymptotique), **pas** le bootstrap résiduel.
- **A.3.vi — diagnostics divergeant des étalons**, à commenter plutôt qu'à taire : biais
  (i) sous l'étalon BDVK05 (−0,06 % / −0,07 % contre 0,1-1,0 %) ; asymétrie H/F (iii)
  conforme en sens mais très atténuée (1,08 contre ≈ 4,0).
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
