# StMoMo — note de référence pour ce projet

Distillé de la vignette StMoMo, du papier JSS (Villegas, Millossovich & Kaishev) et du
manuel de référence du package. **À lire avant de coder les sections A2 et B.**
En cas de doute sur un argument : `?fit.StMoMo`, `?simulate.fitStMoMo`,
`?bootstrap.fitStMoMo` dans R — la doc installée fait foi sur cette note.

---

## 1. Chaîne d'appels pour NOS sections

### Préparation des données
```r
library(demography); library(StMoMo)
# Depuis les fichiers HMD téléchargés :
donnees <- read.demogdata("donnees/Mx_1x1.txt", "donnees/Exposures_1x1.txt",
                          type = "mortality", label = "<PAYS>")
dat_h <- StMoMoData(donnees, series = "male")
dat_f <- StMoMoData(donnees, series = "female")
```
`StMoMoData` produit un objet contenant `Dxt`, `Ext`, `ages`, `years`, `type`.
Pour un Lee-Carter Poisson, `type` doit être `"central"` (expositions centrales,
ETR). Si les données sont en expositions initiales, convertir avec
`central2initial()` / `initial2central()`.

### Section A2 — fit Lee-Carter Poisson
```r
LCfit_h <- fit(lc(link = "log"), data = dat_h,
               ages.fit = <plage>, years.fit = <periode>)
```
- `lc(link = "log")` = Lee-Carter **Poisson** (log-bilinéaire), celui du cours.
  `link = "logit"` donnerait la version binomiale — pas ce qu'on veut.
- Contraintes appliquées par StMoMo : Σβₓ = 1 et Σκₜ = 0. **À vérifier
  numériquement** : `sum(LCfit_h$bx)`, `sum(LCfit_h$kt)`.
- Composants utiles : `$ax`, `$bx`, `$kt`, `$Dxt`, `$Ext`, `$ages`, `$years`.
- Diagnostics : `plot(LCfit_h)` (α, β, κ) ;
  `plot(residuals(LCfit_h), type = "colourmap")` pour la heatmap des résidus.

### Section A2 — projection centrale + bootstrap (incertitude de PARAMÈTRES)
```r
LCfor  <- forecast(LCfit_h, h = <horizon>)                    # RWD par défaut (mrwd)
LCboot <- bootstrap(LCfit_h, nBoot = 5000, type = "semiparametric")
LCsimPU <- simulate(LCboot, h = <horizon>)   # trajectoires AVEC erreur d'estimation
```
- `forecast()` : `kt.method = "mrwd"` par défaut = random walk with drift, exactement
  la spécification du cours. (`"iarima"` permettrait un ARIMA quelconque — ne pas
  s'en servir sans le justifier.)
- Bootstrap semi-paramétrique = ré-échantillonnage des décès selon la loi de Poisson
  supposée (Brouhns, Denuit & Van Keilegom 2005). L'alternative `type = "residual"`
  ré-échantillonne les résidus de déviance. **Semi-paramétrique = le choix par défaut
  et celui à retenir**, mais le noter comme choix justifiable dans le rapport.
- ⚠️ **COÛT : la vignette indique ~2 h pour nBoot = 5000.** Stratégie imposée :
  d'abord `nBoot = 200` pour valider le pipeline de bout en bout, puis relance à 5000
  en arrière-plan avec `saveRDS()` immédiat du résultat. Ne jamais re-bootstrapper
  ensuite : on recharge le RDS.

### Section B — trajectoires de κ SEUL (incertitude de PROJECTION)
```r
set.seed(1234)
LCsim <- simulate(LCfit_h, nsim = 5000, h = <horizon>)   # α, β FIXÉS au fit central
```
**C'est LA distinction du projet, et elle est native dans StMoMo :**

| Appel | Ce qui varie | Section |
|---|---|---|
| `simulate(LCboot, h = ...)` | κ futur **+** α, β ré-estimés **+ dérive et σ du RWD ré-estimées** | A2 (paramètres) |
| `simulate(LCfit, nsim, h = ...)` | κ futur seul, α/β **et** dérive/σ fixés au fit central | **B** (projection) |

⚠️ **Correction** : la ligne A2 ci-dessus portait « κ futur + α, β ré-estimés », ce qui
est incomplet. `simulate(LCboot, ...)` ré-estime AUSSI la dérive et la volatilité de la
marche aléatoire sur la série κ*ₜ de chaque réplication — BDVK05 §4.2 : *spécification
gelée, paramètres relâchés*. La section A.3 capture donc **trois** sources d'incertitude
(α/β, paramètres de la série temporelle, aléa futur de κ) contre **une** en section B.
Détail et références : `docs/BOOTSTRAP_NOTES.md` §3.

La vignette compare explicitement les deux jeux d'IC à 95 % sur un même graphe
(taux à 40, 60, 80 ans) — c'est le graphe à reproduire pour la comparaison du rapport.

Sorties de `simulate()` (classe `simStMoMo`) : `$rates` (array 3D âge × année × sim),
`$kt.s` (trajectoires simulées de κ), `$years`. Pour des quantiles :
```r
q025 <- apply(LCsim$rates, c(1, 2), quantile, probs = 0.025)
```

### Ce que StMoMo NE fait PAS
Aucune fonction de tarification, de VAP, de best estimate ni de SCR.
**Sections C et D : 100 % code custom**, à partir de `LCsim$kt.s` (ou `$rates`)
sauvegardé en RDS par la section B.

---

## 2. Pièges à éviter

- **`ages.fit` / `years.fit`** : ce sont NOS choix de calibration, à noter en TODO
  à justifier (la vignette utilise 55–89 / 1961–2011 pour l'Angleterre-Galles,
  0–89 / 1985–2008 pour la Nouvelle-Zélande — ce sont des exemples, pas une norme).
- **`jumpchoice`** dans `simulate()` : `"fit"` (défaut) part des taux ajustés de la
  dernière année, `"actual"` des taux observés. Choix à connaître si le jury demande
  d'où partent les projections.
- **Ne pas confondre `forecast()` et `simulate()`** : le premier donne la projection
  centrale déterministe (+ IC analytiques), le second des trajectoires stochastiques.
  Les sections C et D ont besoin des trajectoires, pas de la centrale.
- **Le seed** : `simulate()` accepte un argument `seed`, mais on fixe `set.seed(1234)`
  en tête de script pour tout le pipeline (cf. CLAUDE.md).
- **Fonctions modèles disponibles** si besoin de comparaison : `lc()`, `cbd()`,
  `apc()`, `rh()`, `m6()`, `m7()`, `m8()`. Le projet ne demande que Lee-Carter —
  ne pas en ajouter spontanément.
- **Ne pas inventer d'arguments.** Si une option n'apparaît ni dans cette note ni
  dans `?fonction`, elle n'existe probablement pas : vérifier avant d'écrire.
