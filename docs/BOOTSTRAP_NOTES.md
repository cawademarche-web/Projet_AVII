# A.3.vi — Bootstrap semi-paramétrique du Lee-Carter Poisson
### Note théorique/technique — v2, adossée aux papiers originaux

**Sources primaires utilisées :**
- **BDVK05** — Brouhns, N., Denuit, M. & Van Keilegom, I. (2005), « Bootstrapping the Poisson log-bilinear model for mortality forecasting », *Scandinavian Actuarial Journal* 2005(3), 212–224.
- **BDV02a** — Brouhns, N., Denuit, M. & Vermunt, J.K. (2002), « A Poisson log-bilinear regression approach to the construction of projected lifetables », *Insurance: Mathematics and Economics* 31, 373–393.
- **HR09** — Haberman, S. & Renshaw, A. (2009), « On age-period-cohort parametric mortality rate projections », *Insurance: Mathematics and Economics* 45, 255–270.
- Documentation StMoMo (vignette, papier JSS, `?bootstrap.fitStMoMo`).

**Ce qui change par rapport à la v1 :** une correction de terminologie (§2), la mécanique exacte de la boucle bootstrap y compris le ré-ajustement ARIMA (§3), l'identification de la *vraie* alternative discutée par BDVK05 (§4, ce n'était pas le bootstrap résiduel), un résultat de littérature qui prédit ce que la comparaison A3 vs B doit donner (§6), et une précision path-wise critique pour la section D (§7).

---

## 0. Contrat de la section

- **Input** : `LCfit_h`, objet `fitStMoMo` du Lee-Carter Poisson ajusté (section A2).
- **Appel** : `bootstrap(LCfit_h, nBoot = ..., type = "semiparametric")` puis `simulate(LCboot, h = ...)`.
- **Ce que ça mesure** : l'incertitude d'**estimation** des paramètres (α̂ₓ, β̂ₓ, κ̂ₜ), propagée jusqu'aux quantités actuarielles — *plus* l'incertitude de projection, puisque `simulate()` reprojette chaque réplication (cf. §3).
- Cette section produit le second jeu d'intervalles que le rapport doit comparer à ceux de la section B (`simulate(LCfit, nsim, h)`, projection seule).

---

## 1. Pourquoi un bootstrap ici (intuition avant formalisme)

Le Lee-Carter original s'estime par SVD des log-taux observés, en minimisant Σ(ln μ̃ₓ(t) − αₓ − βₓκₜ)². BDV02a (§3.2, §4.1) rappelle le défaut central : la SVD suppose des erreurs **homoscédastiques**, c'est-à-dire implicitement normales et de variance constante — hypothèse irréaliste, puisque le log du taux brut est bien plus variable aux grands âges, où les décès sont peu nombreux en valeur absolue. Le nombre de décès étant une variable de comptage, l'hypothèse de Poisson est la spécification naturelle (BDV02a cite Brillinger 1986 sur ce point).

D'où la reformulation :

> Dₓₜ ~ Poisson(ETRₓₜ · μₓ(t)), avec μₓ(t) = exp(αₓ + βₓκₜ), sous Σκₜ = 0 et Σβₓ = 1 [BDVK05 (3.3) ; BDV02a (4.1)]

estimée par maximisation de la log-vraisemblance

> L(α, β, κ) = Σₓ Σₜ { Dₓₜ(αₓ + βₓκₜ) − ETRₓₜ·exp(αₓ + βₓκₜ) } + cte [BDVK05 §3.2 ; BDV02a §4.2]

**C'est ce modèle-là que `fit(lc(link = "log"))` ajuste dans StMoMo** — pas la SVD.

**Trois points de défense qui viennent directement de BDV02a §4.2 :**

1. Le terme bilinéaire βₓκₜ empêche l'usage d'une routine GLM standard : BDV02a résout par une méthode de Newton uni-dimensionnelle (Goodman 1979), qui met à jour un bloc de paramètres à la fois (α, puis κ, puis β) en gelant les autres. Savoir *pourquoi* ce n'est pas un GLM ordinaire est le genre de question que pose un jury.
2. **Les contraintes ne sont pas dans la vraisemblance, elles sont ré-imposées à chaque itération** : Σκₜ = 0 par centrage après mise à jour des κ, Σβₓ = 1 par rescaling. Ce sont des contraintes d'identification, pas des restrictions de modèle — c'est précisément pourquoi le simple centrage/rescaling suffit (BDV02a le dit explicitement).
3. **Pas de seconde étape de recalage des κ.** Dans le Lee-Carter SVD classique, on ré-estime les κ̂ₜ pour reproduire le total de décès observés Σₓ Dₓₜ (BDV02a éq. 3.4). En Poisson, l'erreur porte directement sur le nombre de décès, donc « il n'y a pas besoin d'une estimation de seconde étape » (BDV02a §4.2, fin). Si le jury demande « et le second-stage adjustment de Lee-Carter ? », la réponse est : sans objet ici, et voici pourquoi.

Une fois le modèle probabiliste posé, on peut ré-échantillonner. Mais l'intractabilité analytique demeure : BDVK05 §4.1 est explicite — il faut combiner **deux sources d'incertitude de nature différente** (erreur d'échantillonnage sur α, β, κ et erreur de prévision de l'ARIMA), et les quantités d'intérêt (eₓ(t), aₓ(t)) sont des fonctions **non linéaires compliquées** de ces paramètres. Pas de formule fermée → bootstrap.

**Message pour la défense** : le bootstrap n'est pas un supplément optionnel. C'est la réponse standard à une intractabilité analytique explicitement documentée dans la littérature de référence.

---

## 2. ⚠️ Terminologie : « semi-paramétrique » est le mot de StMoMo, pas celui de BDVK05

**À savoir avant la défense.** L'article original n'emploie **jamais** le terme *semiparametric*. BDVK05 §1 annonce comparer « une procédure de bootstrap **non paramétrique** à la procédure de bootstrap **paramétrique** proposée par Brouhns et al. [7] », et la section 4.2 s'intitule **« Poisson bootstrap »**.

L'étiquette *semiparametric* vient de la documentation StMoMo (vignette §8, `?bootstrap.fitStMoMo`), qui l'utilise pour distinguer cette approche du *residual bootstrap*.

Pourquoi les deux appellations se défendent :
- **« Non paramétrique » (BDVK05)** : par contraste avec le bootstrap paramétrique de Brouhns-Denuit-Vermunt (2002b), qui tire directement dans la loi normale multivariée asymptotique des estimateurs du maximum de vraisemblance. Ici on ne postule aucune loi *sur les paramètres*.
- **« Semi-paramétrique » (StMoMo)** : par contraste avec le bootstrap résiduel, qui ré-échantillonne empiriquement. Ici on impose bien une famille paramétrique (Poisson) *au tirage des pseudo-données*.

**À faire dans le rapport** : écrire « bootstrap semi-paramétrique au sens de StMoMo, désigné *Poisson bootstrap* par Brouhns, Denuit & Van Keilegom (2005) ». Une phrase, et le point est neutralisé. Sans elle, un membre du jury qui a l'article en tête peut légitimement demander « où voyez-vous *semiparametric* dans BDVK05 ? » — et il aura raison.

---

## 3. Mécanique exacte, telle que spécifiée dans BDVK05 §4.2

Point de départ : les observations (ETRₓₜ, Dₓₜ) et le fit central α̂ₓ, β̂ₓ, κ̂ₜ.

Pour n = 1, …, N :

**Étape 1 — tirage des pseudo-décès.** On génère D*ₓₜ ~ Poisson(moyenne = ETRₓₜ · μ̂ₓ(t)) où μ̂ₓ(t) est l'estimateur **non contraint** du maximum de vraisemblance, c'est-à-dire le taux brut μ̂ₓ(t) = Dₓₜ / ETRₓₜ (BDVK05 éq. 2.5). Donc :

> **ETRₓₜ · μ̂ₓ(t) = Dₓₜ** — la moyenne du tirage Poisson **est le nombre de décès observé**.

BDVK05 le formule d'ailleurs directement ainsi : les pseudo-décès sont obtenus « en appliquant un bruit de Poisson aux nombres de décès observés ».

**Trois conséquences opérationnelles :**
- Cela correspond exactement à `deathType = "observed"` dans StMoMo, le défaut. La variante `"fitted"` (tirage autour des D̂ₓₜ ajustés) est celle de Renshaw & Haberman (2008), **pas** celle de BDVK05. Garder `"observed"`.
- Le μ̂ₓ(t) de l'étape 1 est **littéralement le taux MLE brut calculé en section A1** de ton projet. C'est un lien inter-sections à faire explicitement dans le rapport : A1 ne sert pas qu'à produire un graphe, il fournit l'objet autour duquel le bootstrap de A3 tire.
- Les ETRₓₜ sont **fixes** à chaque réplication — jamais ré-échantillonnées (elles apparaissent inchangées dans le couple (ETRₓₜ, D*ₓₜ) chez BDVK05). Hypothèse simplificatrice à nommer en limite (§9).

**Étape 2 — ré-estimation du modèle.** Les αₓ, βₓ et κₜ sont ré-estimés sur (ETRₓₜ, D*ₓₜ) par la même procédure MLE, avec ré-imposition des contraintes.

**Étape 3 — reprojection ARIMA. ⚠️ Point que la v1 de cette note omettait.** BDVK05 §4.2 est explicite : « les κₜ sont ensuite projetés sur la base du modèle ARIMA **ré-estimé** ». Et immédiatement la précision qui compte :

> On **ne sélectionne pas** un nouveau modèle ARIMA — on garde la spécification retenue sur les données originales (ARIMA(0,1,0) hommes, ARIMA(0,1,1) femmes chez eux). En revanche, **les paramètres de ce modèle sont ré-estimés sur les données bootstrappées**.

Autrement dit : **spécification gelée, paramètres relâchés**. Chaque réplication a sa propre dérive et sa propre volatilité estimées.

C'est exactement ce que fait `simulate(LCboot, h = ...)` avec `kt.method = "mrwd"` (le défaut) : la spécification random walk with drift est fixée par toi, mais la dérive et σ sont ré-estimées sur la série κ*ₜ de chaque réplication.

**Ce que ça implique pour la comparaison A3 / B — et c'est le cœur du rapport :**

| | α, β | dérive & σ du RWD | aléa futur de κ |
|---|---|---|---|
| `simulate(LCfit, nsim, h)` — **section B** | fixés au fit central | **fixés** au fit central | simulé |
| `simulate(LCboot, h)` — **section A3** | ré-estimés par réplication | **ré-estimés** par réplication | simulé |

La section A3 capture donc **trois** sources, pas deux : incertitude sur (α, β), incertitude sur les paramètres de la série temporelle, et aléa de projection. La v1 de cette note ne mentionnait que la première. Si le jury demande « qu'est-ce que votre bootstrap ajoute exactement par rapport à la section B ? », c'est cette ligne-là qu'il faut pouvoir énoncer.

**Étape 4** — calcul de la quantité d'intérêt sur chaque réplication → distribution empirique → quantiles.

---

## 4. La vraie alternative de BDVK05 : le bootstrap **paramétrique**, pas le résiduel

Correction importante par rapport à la v1. Le papier construit sa comparaison principale non pas contre le bootstrap résiduel, mais contre le **bootstrap paramétrique** de Brouhns, Denuit & Vermunt (2002b, *Bulletin of the Swiss Association of Actuaries*), décrit en BDVK05 §4.1 :

> « Brouhns et al. [7] ont tiré directement dans la loi normale multivariée approchée des estimateurs du maximum de vraisemblance α̂, β̂, κ̂. Nous proposons ici une approche alternative. »

**Il y a donc trois familles, pas deux :**

| Approche | Ce qu'on tire | Statut chez BDVK05 | Dans StMoMo |
|---|---|---|---|
| **Paramétrique** (BDV 2002b) | Directement les paramètres, dans leur loi normale multivariée asymptotique (matrice d'information de Fisher) | Le point de comparaison principal du papier | ❌ non implémenté |
| **Poisson / semi-paramétrique** (BDVK05) | Les pseudo-décès D*ₓₜ ~ Poisson(Dₓₜ), puis ré-estimation complète | La contribution du papier | ✅ `type = "semiparametric"` |
| **Résiduel** (Koissi et al. 2006) | Les résidus de déviance, avec remise, puis inversion vers des D*ₓₜ | Mentionné en une phrase, Remarque 4.1 fin | ✅ `type = "residual"` |

**Ce que ça change pour l'argumentaire du rapport** — et c'est un gain net :

Le bootstrap paramétrique repose sur la **normalité asymptotique de l'estimateur du MLE**, exactement la machinerie de la delta-method de ton Q2 d'examen. L'argument devient donc beaucoup plus fort et beaucoup plus « cours » que ce que j'avais écrit en v1 :

> Le bootstrap Poisson évite d'avoir à invoquer la normalité asymptotique des estimateurs, et donc à supposer que l'approximation asymptotique est valable à taille d'échantillon finie — hypothèse d'autant plus discutable aux grands âges, où les effectifs de décès sont faibles et où c'est précisément le comportement de la queue qui pilote la valeur d'une rente viagère.

Et BDVK05 fournit la validation empirique de ce choix : sur les données belges, les deux procédures donnent des intervalles **comparables** et des histogrammes de forme similaire. Les auteurs en tirent une conclusion méthodologique explicite et réutilisable :

> Cette concordance entre bootstrap paramétrique et bootstrap Poisson peut être considérée comme **une indication que le modèle log-bilinéaire Poisson ajuste raisonnablement bien les données** (BDVK05 §4.3).

Autrement dit, l'accord entre les deux méthodes est lui-même un **test de spécification informel**. À citer si tu veux muscler la justification.

**Une quatrième variante, la Remarque 4.1 de BDVK05** (mentionne-la seulement si tu veux montrer que tu as lu le papier en entier) : bootstrapper **par génération** plutôt que par cellule, en tirant les pseudo-décès le long d'une diagonale (ETRₓₜ, Dₓₜ), (ETRₓ₊₁,ₜ₊₁, Dₓ₊₁,ₜ₊₁), … dans une loi multinomiale d'exposant D = Σₖ Dₓ₊ₖ,ₜ₊ₖ et de paramètres Dₓₜ/D, Dₓ₊₁,ₜ₊₁/D, … BDVK05 note que c'est très proche du bootstrap Poisson, puisque **la loi conditionnelle de comptages de Poisson indépendants sachant leur somme est multinomiale** — la différence étant que le bootstrap Poisson ne garantit pas la conservation du total D. C'est une remarque de fond utile si le jury demande « et si on voulait préserver le nombre total de décès ? ».

---

## 5. Justification du choix retenu (structure pour le rapport)

Quatre arguments, du plus fort au plus circonstanciel. À reformuler avec tes mots.

1. **Cohérence avec le mécanisme d'estimation.** Le fit de A2 est un MLE sous hypothèse Poisson ; le bootstrap ré-échantillonne selon ce même mécanisme générateur. C'est le bootstrap naturel du modèle, pas une approximation externe.
2. **Pas de recours à la normalité asymptotique.** Contrairement au bootstrap paramétrique de BDV (2002b), on n'a pas à supposer que la loi normale multivariée approche correctement la loi de (α̂, β̂, κ̂) à taille finie — hypothèse fragile aux âges élevés où les décès sont rares (§4).
3. **Propagation complète et documentée des sources d'incertitude**, y compris la ré-estimation des paramètres de la série temporelle à chaque réplication (BDVK05 §4.2), ce que ne fait pas la simulation de la section B.
4. **Diagnostic empirique en main.** Si la heatmap des résidus de A2 ne montre pas de structure marquée, l'hypothèse Poisson n'est pas contredite. ⚠️ **Vérifie ce diagnostic avant d'écrire cet argument** : BDVK05 §3.2 fonde exactement le même raisonnement, en observant que « l'absence de structure à la plupart des âges (sauf aux tout premiers, généralement inutilisés en calcul actuariel) soutient le modèle ». Note bien la réserve sur les jeunes âges — si ta plage `ages.fit` démarre bas, attends-toi à des résidus dégradés en bas de grille, et dis-le plutôt que de le cacher.

Le résiduel reste mentionné comme alternative légitime en cas de mauvaise spécification soupçonnée (sur-dispersion, structure résiduelle diagonale non captée) : son tirage i.i.d. des résidus détruit toute dépendance âge/période réelle, ce qui en fait un mauvais choix précisément quand un effet cohorte est visible. À formuler comme limite reconnue, pas comme argument central.

**Formule de déviance à citer si on te demande la définition exacte des résidus** (BDVK05 §3.2) :

> rₓₜ = signe(Dₓₜ − D̂ₓₜ) · √( 2 [ Dₓₜ·ln(Dₓₜ/D̂ₓₜ) − (Dₓₜ − D̂ₓₜ) ] ),  avec D̂ₓₜ = ETRₓₜ·exp(α̂ₓ + β̂ₓκ̂ₜ)

*(BDVK05 imprime cette expression sans le facteur 2 sous la racine ; c'est une coquille manifeste de l'article — la déviance de Poisson est D = 2Σ[d·ln(d/d̂) − (d − d̂)]. Ne pas recopier la version de l'article telle quelle.)*

---

## 6. 🎯 Ce que la comparaison A3 vs B **doit** donner — résultat de littérature à connaître

C'est l'apport le plus utile des papiers uploadés, et ça te donne un contrôle qualitatif sur ta section B.

La chaîne de références est nette et convergente :

- **Lee & Carter (1992, Appendice B)** : pour les espérances de vie, il est raisonnable de se restreindre aux erreurs de prévision de l'indice de mortalité et d'ignorer celles d'ajustement de la matrice — y compris à court terme. BDV02a §5.4 s'appuie explicitement là-dessus pour ne fonder ses intervalles que sur la variabilité de l'indice κ.
- **BDVK05 §1** rappelle ce même constat : dans Brouhns et al. [6], les intervalles étaient obtenus en ignorant toutes les erreurs sauf celle de prévision de l'indice, ces erreurs dominant les autres pour les rentes et les espérances de vie résiduelles.
- **HR09 §2.9** : la résolution retenue « implique la simulation de l'erreur de prévision de l'ARIMA(p,1,q) appliqué à l'effet période principal, tout en ignorant l'erreur d'ajustement du modèle : cette dernière composante se révèle d'effet négligeable en comparaison, résultat globalement en accord avec les conclusions originales de Lee et Carter (1992) ».
- **HR09 §5 (Summary)** confirme : il existe des éléments (Denuit et al., 2009) en faveur du constat original selon lequel l'erreur de prévision de la série temporelle de l'indice période **domine** la source d'erreur d'estimation paramétrique dans la structure LC.

**Conséquence directe pour ton rapport :**

> Les intervalles de A3 (avec incertitude de paramètres) doivent être **modérément plus larges** que ceux de B (projection seule). Pas spectaculairement.

- ✅ Résultat attendu : élargissement perceptible mais de second ordre. Tu l'écris, et tu le rattaches à Lee & Carter (1992, App. B), BDVK05 et HR09 — ça transforme un simple graphe en discussion étayée.
- 🚩 Si tes intervalles A3 sont **massivement** plus larges (facteur 2 ou plus) : suspecte d'abord un bug, pas une découverte. Vérifie que la section B ne fixe pas trop de choses, ou que le bootstrap n'accumule pas les tirages.
- 🚩 Si les intervalles A3 sont **plus étroits** que ceux de B : bug certain. A3 englobe strictement les sources de B.

**Nuance à garder en réserve** (BDVK05 §8 de la vignette StMoMo, et §3 du papier JSS) : l'incertitude de paramètres pèse davantage sur les **petites populations**. C'est pourquoi la vignette bascule sur la Nouvelle-Zélande (4,3 M d'habitants) plutôt que l'Angleterre-Galles (54,8 M) pour illustrer le bootstrap. Si ton pays assigné est petit, attends-toi à un élargissement plus net, et dis-le — c'est un point d'interprétation qui rapporte.

---

## 7. ⚠️ Précision path-wise — critique pour ta section D (SCR)

HR09 §4.7–4.9 compare deux algorithmes de simulation qui ont **les mêmes lois marginales** mais produisent des **trajectoires de formes différentes**. Sous RWD pur (φ = 0, dérive θ̂, volatilité σ̂) :

- **Algorithme 1** : κ*ₜₙ₊ⱼ = κₜₙ + j·θ̂ + √j·σ̂·z*, avec **un seul** z* ~ N(0,1) par trajectoire.
- **Algorithme 2** : κ*ₜₙ₊ⱼ = κₜₙ + j·θ̂ + σ̂·Σᵢ₌₁ʲ z*ᵢ, avec des z*ᵢ ~ N(0,1) **i.i.d.** — la vraie marche aléatoire.

HR09 montre que ces deux expressions sont marginalement identiques (puisque Σⱼ z*ᵢ ~ N(0, j)). Mais la Figure 15 est parlante : sous l'algorithme 1, les 40 trajectoires forment un éventail lisse et ordonné — chaque trajectoire est une droite déviée ; sous l'algorithme 2, on obtient de vraies marches aléatoires qui se croisent et serpentent. HR09 §4.8 note aussi que l'algorithme 1 produit systématiquement des intervalles de prédiction plus larges.

**Pourquoi ça t'importe :** ta section D calcule un **SCR = VaR₉₉,₅% de −ΔBOF** sur les trajectoires. Une VaR à 99,5 % dépend de la queue de la distribution **du chemin complet** (décès année 1, état 2023, BE₂₀₂₃ …), pas seulement de la marginale à l'horizon final. Un éventail lisse type algorithme 1 déformerait la dépendance temporelle et fausserait le SCR.

**Bonne nouvelle :** `simulate()` de StMoMo simule de vraies trajectoires de marche aléatoire (logique algorithme 2). Tu n'as rien à corriger. Mais **sache le dire** — « pourquoi vos 5000 trajectoires et pas simplement 5000 tirages de la loi de κ à l'horizon final ? » est exactement le type de question de défense où cette distinction sépare celui qui a compris de celui qui a lancé une fonction.

**Une source d'incertitude que même A3 ne capture pas** (HR09 §4.7, Algorithme 2 étape 1) : on pourrait aussi simuler l'erreur d'estimation des paramètres de la série temporelle depuis leurs lois a posteriori marginales — (n−1)σ̂²/σ² ~ χ²ₙ₋₁ et θ ~ N(θ̂, σ̂²/(n−1)) sous a priori vague. HR09 §4.9 le teste et conclut à un effet **marginal à 65 ans et négligeable à 70 ans et au-delà**. À citer en limite : « une provision supplémentaire pour l'erreur d'estimation de la dérive aurait pu être ajoutée ; la littérature (Haberman & Renshaw 2009) la documente comme d'effet négligeable aux âges de retraite. » Note toutefois que ton bootstrap la capture déjà partiellement via la ré-estimation du RWD à chaque réplication (§3).

---

## 8. Contrat technique pour Claude Code

```r
# Section A2 (déjà fait) : LCfit_h <- fit(lc(link = "log"), data = dat_h, ...)

set.seed(1234)  # en tête de script (cf. CLAUDE.md). bootstrap() n'a pas d'argument
                 # seed propre : set.seed() gouverne les tirages Poisson internes.

# 1) VALIDATION du pipeline — nBoot réduit, étape non négociable
LCboot_test <- bootstrap(LCfit_h, nBoot = 200,
                          type      = "semiparametric",  # = Poisson bootstrap, BDVK05
                          deathType = "observed")        # = BDVK05 §4.2 au sens strict
saveRDS(LCboot_test, "resultats/LCboot_test_200.rds")

# Valider AUSSI l'aval avant de lancer le run long :
LCsimPU_test <- simulate(LCboot_test, h = <horizon>)   # kt.method = "mrwd" par défaut

# 2) Run de production, en arrière-plan (~2 h d'après la vignette pour nBoot = 5000)
LCboot <- bootstrap(LCfit_h, nBoot = 5000,
                     type = "semiparametric", deathType = "observed")
saveRDS(LCboot, "resultats/LCboot_5000.rds")

LCsimPU <- simulate(LCboot, h = <horizon>)
saveRDS(LCsimPU, "resultats/LCsimPU_5000.rds")

# 3) Ne JAMAIS re-bootstrapper ensuite : recharger le RDS.
```

Vigilances :
- Écrire `type` **et** `deathType` explicitement, même si ce sont les défauts — ça documente le choix dans le code et protège d'un changement de défaut en amont.
- Ne pas fusionner cette fonction avec celle de la section B. Deux fonctions distinctes (cf. le 🚩 du checklist jour 2).
- L'objet qui alimente la comparaison du rapport est `LCsimPU`, pas `LCboot`.
- BDVK05 tourne à **N = 10 000** réplications ; la vignette StMoMo à 5 000. Ton `nBoot = 5000` est dans la norme de la littérature — c'est un chiffre justifiable, à mentionner en une ligne dans le rapport.

---

## 9. Diagnostics de contrôle après bootstrap

**Contrôle 1 — biais faible.** La moyenne des N réplications d'une quantité d'intérêt doit être proche de la prévision ponctuelle. BDVK05 §4.3 donne les ordres de grandeur sur données belges (N = 10 000, e₆₅(2000) et a₆₅(2000)) :

| | moyenne bootstrap | prévision ponctuelle | écart |
|---|---|---|---|
| e₆₅ hommes | 16,04 | 16,01 | 0,2 % |
| e₆₅ femmes | 20,60 | 20,40 | 1,0 % |
| a₆₅ hommes | 10,71 | 10,69 | 0,2 % |
| a₆₅ femmes | 13,01 | 13,02 | 0,1 % |

Si ton écart moyenne/ponctuel dépasse largement le pourcent, cherche le bug.

**Contrôle 2 — largeur relative des intervalles.** BDVK05 obtient sur intervalles à 90 % : a₆₅ hommes [9,86 ; 11,55] et femmes [12,76 ; 13,28], soit des largeurs relatives de **15,7 % (hommes) et 3,9 % (femmes)**, qualifiées de « raisonnablement précises compte tenu de la projection à 50 ans qu'elles requièrent ». C'est ton étalon d'ordre de grandeur — quelques pourcents à ~15 % sur une VAP de rente. Une largeur de 50 % ou de 0,5 % signale un problème. Attention : leurs intervalles sont à **90 %**, les tiens à 95 % seront mécaniquement plus larges.

**Contrôle 3 — asymétrie hommes/femmes.** BDVK05 obtient systématiquement des intervalles **plus larges pour les hommes**, et l'attribue aux intervalles plus larges sur la projection des κ masculins (§4.3, renvoi à leur Figure 3). Si tu obtiens l'inverse sur ton pays, ce n'est pas nécessairement faux, mais c'est à examiner et à commenter.

**Contrôle 4 — hiérarchie des dispersions paramétriques.** La vignette StMoMo (§8) constate que l'incertitude sur αₓ et κₜ est modeste, tandis que celle sur β⁽¹⁾ₓ est nettement plus marquée. Dispersion de β̂ₓ plus large aux âges extrêmes (faible exposition). Si tu observes l'inverse, creuse.

**Contrôle 5 — cohérence des ordres de grandeur avec ta section C.** BDV02a calcule ses rentes à i = 4 % et obtient a₆₅ ≈ 10,7 (hommes) et 13,2 (femmes) sur données belges. Ton taux de tarification est de 3 %, donc tes VAP seront **plus élevées** (actualisation plus faible). Ton heuristique de contrôle existante (rente viagère à 65 ans à t = 3 % ≈ 15–18× la rente annuelle) reste cohérente avec cet ancrage.

**Contrôle 6 — qualité d'ajustement globale, en amont.** BDV02a §5.1 utilise la statistique du rapport de vraisemblance L² = 2ΣₓΣₜ Dₓₜ·ln(Dₓₜ/D̂ₓₜ) contre le modèle saturé, et rapporte que le Lee-Carter réduit L² de **91,4 % (femmes) et 85,7 % (hommes)** par rapport au modèle à taux constants dans le temps. Si tu veux un chiffre synthétique d'ajustement à mettre dans le rapport avant la partie bootstrap, c'est celui-là — comparable, et tu as le benchmark belge.

---

## 10. Position dans le pipeline

```
A1 — taux MLE bruts μ̂ₓ(t) = Dₓₜ/ETRₓₜ
   │        └──────────────────┐  (c'est la moyenne du tirage Poisson, BDVK05 §4.2)
A2 — fit(lc(link="log"))       │
   │  → LCfit_h                │
   ├── forecast(LCfit_h, h)    │        → projection centrale déterministe
   │                           ▼
   ├── bootstrap(LCfit_h, nBoot=5000, "semiparametric")  → LCboot   [A.3.vi]
   │        └── simulate(LCboot, h)  → LCsimPU  (α,β + dérive RWD + aléa κ)
   │
   └── simulate(LCfit_h, nsim=5000, h) → LCsim  [section B] (aléa κ seul)
                │
                └──→ sections C (VAP, primes) et D (BE, BOF, SCR) — code custom
```

Le graphe de référence du rapport : comparaison des largeurs d'IC de `LCsimPU` et `LCsim` sur un même graphique (taux à 40/60/80 ans par exemple), avec la discussion du §6 en commentaire.

---

## 11. Limites à énoncer en défense

1. **Expositions traitées comme fixes.** Le bootstrap ne propage aucune incertitude sur ETRₓₜ (BDVK05 conserve (ETRₓₜ, D*ₓₜ)). Note connexe : BDV02a §2.2 rappelle que sous hypothèse de force de mortalité constante par morceaux, l'exposition elle-même est reconstruite par ETRₓₜ = Lₓₜ · qₓ(t)/ln(1 − qₓ(t)) quand on ne dispose que des effectifs au 1ᵉʳ janvier — donc elle est déjà un objet estimé, traité ici comme connu.
2. **Hypothèse Poisson non testée formellement.** Diagnostic visuel seulement (heatmap des résidus). BDVK05 procède exactement de même (§3.2). Une piste de renforcement citable : l'accord entre bootstrap Poisson et bootstrap paramétrique est interprété par BDVK05 §4.3 comme un indice d'ajustement correct.
3. **Risque de modèle non couvert.** Le bootstrap quantifie l'incertitude d'estimation *à l'intérieur* de la famille Lee-Carter Poisson, pas le risque que Lee-Carter soit le mauvais modèle. HR09 est précisément une démonstration de ce point : en ajoutant un effet cohorte modulé par l'âge (APC), les prévisions de mortalité deviennent **plus favorables** (mortalité plus basse) que sous LC pour l'Angleterre-Galles, et les espérances de vie prédites à 65 ans passent par exemple de 19,16 (LC) à 19,81 (APC) pour l'époque 2007. Ton projet ne demande que Lee-Carter — mais si ta heatmap de résidus A2 montre des bandes diagonales, tu peux dire précisément dans quel sens ton modèle est probablement biaisé, ce qui vaut bien mieux qu'un « le modèle a des limites » générique.
4. **Extrapolation pure du passé.** BDV02a §3.4 le formule bien : la méthode n'incorpore aucune hypothèse sur les progrès médicaux ou les changements environnementaux, et ne peut donc anticiper ni une percée thérapeutique ni une épidémie. Ils citent aussi la réponse standard à cette critique (Wilmoth 2000) : elle n'est valable que dans la mesure où les mécanismes sous-jacents seraient compris assez précisément pour offrir une méthode de prédiction alternative légitime — ce qui n'est pas le cas. **Sache citer les deux côtés**, c'est le genre de nuance qui distingue une limite récitée d'une limite comprise.
5. **N = 200 est insuffisant pour des quantiles extrêmes** (VaR 99,5 % en section D) — d'où le passage obligé à 5000 avant tout calcul de SCR.

---

## 12. Ce que les papiers uploadés apportent **au-delà** de la section A.3.vi

Trois pistes rentables ailleurs dans ton travail :

- **BDV02a §6 — antisélection.** Modèle relationnel de type Brass reliant la mortalité des rentiers à celle de la population générale : ln μᴿᴬₓ = ϑ₁ + ϑ₂ ln μᴺᴵˢₓ, estimé sur données belges. L'impact sur les primes pures de rentes atteint **~15 % pour les femmes et ~26 % pour les hommes**. Ton projet tarifie sur mortalité de population générale ; si le jury demande « et la mortalité de vos assurés ? », ce chiffre est une réponse chiffrée et sourcée. C'est aussi du matériau direct pour la partie « perspective marché de l'assurance » du cours.
- **HR09 §2.5 — identifiabilité APC.** Le papier démontre visuellement (Fig. 1) qu'en ajustant log μₓₜ = αₓ + κₜ + ιₜ₋ₓ en un seul passage, l'**ordre de spécification des effets dans la formule change complètement les patterns** de κ et ι, seul αₓ restant stable. C'est le problème d'identifiabilité APC en image, causé par l'identité année de naissance = année calendaire − âge. **Ton Q1 d'examen porte sur l'identifiabilité APC** — cette figure est le meilleur support d'intuition que tu aies dans ce corpus.
- **HR09 §3.10 & 4.2 — choix de la période de calibration.** Ton `years.fit` est un TODO à justifier. HR09 documente la procédure : suivre le profil du R² d'une régression linéaire sur κₜ construite séquentiellement à rebours (méthode de Denuit & Goderniaux 2005) pour repérer le point de départ de la non-linéarité, puis ré-estimer αₓ et βₓ sur la période réduite via le GLM Poisson log μₓₜ = αₓ + βₓ(t − t̄). Ils notent que κₜ est nettement plus linéaire sous APC que sous LC, la courbure résiduelle du κ de LC étant précisément ce que le terme cohorte absorbe. Voilà une justification méthodologique référencée pour ton choix de période, au lieu d'un choix implicite.
