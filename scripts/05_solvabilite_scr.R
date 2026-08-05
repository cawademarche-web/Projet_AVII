set.seed(1234)
# ==============================================================
# Script 05 — Section D : fonds propres (BOF) et capital de solvabilité (SCR)
# Projet ACTU-F502 — Assurance Vie II
#
#   BOF_t = A_t − BE_t                       (fonds propres = actifs − provision)
#   BE_t  = Σ_sexe (assurés en vie) × a_{65+m} ,  actualisé à taux_actu_be
#   SCR_t = VaR_{99,5 %}( BOF_t − BOF_{t+1}/(1+taux_actu_be) )        (brief, D)
#
# D.1  bloc 4 — A₂₀₂₂, BE₂₀₂₂, BOF₂₀₂₂ du portefeuille [A]
# D.2  bloc 4 — état 2023 par trajectoire, BOF₂₀₂₃, SCR₂₀₂₂
# D.3  bloc 4 — idem portefeuille [B], comparaison
# D.4  bloc 5 — SCR₂₀₂₃ en lecture « état central »
# D.5  bloc 6 — sensibilité rdt_actifs = 1 %
#      bloc 7 — diagnostics de structure   · bloc 8 — oracles
#      bloc 9 — limites chiffrées          · blocs 10-11 — figures, sauvegardes
#
# Code 100 % CUSTOM : StMoMo n'a aucune fonction de solvabilité (décision
# méthodologique verrouillée, CLAUDE.md). Base R uniquement, aucun library().
#
# set.seed EN LIGNE 1, et c'est un SIGNAL : ce script TIRE. La binomiale des
# décès est la SEULE source d'aléa nouvelle de la section D — les 5000
# trajectoires de κ arrivent toutes faites du script 03 via le script 04, et
# re-simuler κ ici serait une violation du contrat d'interface.
#
# CONVENTION DE DATE (t⁺, après flux de l'année) : le paiement k tombe EN DATE
# 2022+k, c'est-à-dire à la clôture de la k-ième année d'assurance. ₁p₆₅ =
# exp(−µ₆₅(2022)) porte sur l'ANNÉE CIVILE 2022 : les survivants de la 1ʳᵉ année
# touchent leur rente en date 2023. Ne jamais écrire « fin de l'année 2022+k »,
# ce serait décalé d'un an.
# ==============================================================

# ==============================================================
# Bloc 1 — Les deux taux, et le chargement (contrat d'interface)
# ==============================================================
# DEUX taux dans ce script, et deux seuls. Le troisième — le taux technique de
# 3 % de la section C — N'EST PAS INSTANCIÉ ICI : les primes arrivent toutes
# faites du RDS, déjà tarifées. La valeur numérique du taux technique n'apparaît
# nulle part dans ce fichier, et aucune variable ne s'appelle r, i ou t : ces
# deux greps sont le contrôle de confusion de taux le plus rapide de la section.
taux_actu_be <- 0.02   # r du brief : actualise le BE à CHAQUE date, ET BOF_{t+1}
                       # dans la formule du SCR. Le seul taux de calcule_BE().
rdt_actifs   <- 0.04   # i du brief : fait croître les ACTIFS uniquement, une
                       # fois par an. Sensibilité D.5 : 1 %.
rdt_grille   <- c(0.04, 0.01)
n_tetes      <- 500    # 500 hommes ET 500 femmes, dans chaque portefeuille

# Entrée UNIQUE : la branche $base du script 04. La branche $covid en est la
# voisine et n'est JAMAIS lue (portée COVID bornée à la section C, décision
# actée). 03_diagonales_cohorte.rds n'est pas rechargé non plus : q₆₅, q₆₆ et
# toutes les annuités se dérivent de S.
vap     <- readRDS("resultats/04_vap.rds")
S       <- vap$base$survies   # list(h, f) : 37 × 5000, S[k, j] = ₖp₆₅ de la
                              # trajectoire j, rownames "1".."37"
primes  <- vap$base$primes    # list(A, B), primes unitaires H et F (rente de 1)
sens    <- vap$base$sensibilite   # bloc ORACLE 2 uniquement
n_sim   <- vap$n_sim          # 5000
n_A     <- vap$n_A            # 37 termes pour la viagère (clôture ℓ₁₀₂)
duree_B <- vap$duree_B        # 15 termes pour la temporaire

# ==============================================================
# Bloc 2 — Fonctions (une fonction = un concept actuariel)
# ==============================================================
# Le taux est passé en ARGUMENT et jamais lu comme variable globale, alors même
# qu'il vaut toujours taux_actu_be : c'est le seul taux du passif. Écart assumé
# à la règle « paramètre jamais varié → constante », parce que le risque n°1 de
# la section est la confusion de taux et qu'un taux visible au SITE D'APPEL rend
# le 🚩 de la checklist vérifiable d'un coup d'œil.

# annuite_conditionnelle : a_{65+m} à TERME ÉCHU, conditionnelle à la survie
# jusqu'à l'âge 65+m :
#   a_{65+m} = Σ_{k=1}^{n} vᵏ · ₖp_{65+m} ,  avec  ₖp_{65+m} = ₖ₊ₘp₆₅ / ₘp₆₅
# Le passage par le quotient est ce qui rend la ré-évaluation du BE possible :
# à la date 2023 l'assuré A SURVÉCU sa première année, on renormalise donc par
# ₁p₆₅. Indexation par NOM de ligne (règle CLAUDE.md), jamais par position.
annuite_conditionnelle <- function(S, m, n, taux) {
  p_m    <- if (m == 0) 1 else S[as.character(m), ]     # ₘp₆₅ (5000 valeurs, ou 1)
  lignes <- as.character((m + 1):(m + n))               # bornes explicites
  # unname : S porte les étiquettes de simulation en colonnes, qui n'ont aucun
  # sens actuariel ici et se propageraient jusque dans les BOF sauvegardés.
  unname(colSums((1 / (1 + taux))^(1:n) * S[lignes, , drop = FALSE]) / p_m)
}

# n_reste : nombre de paiements RESTANT après la date 2022+m.
# ⚠️ PIÈGE D'INDICE DE LA SECTION : à la date 2023 (m = 1) la prestation k = 1 de
# la temporaire [B] est DÉJÀ versée — il en reste 14, pas 15 ; à la date 2024
# (m = 2), 13. Sur [A] une borne trop haute ferait planter R (m+n > 37) ; sur
# [B] elle renverrait SILENCIEUSEMENT les lignes "2".."16" et un a₆₆ faux de
# +4,8 %. Rien ne protège la temporaire : c'est l'ORACLE 3 qui la garde.
n_reste <- function(produit, m) (if (produit == "A") n_A else duree_B) - m

# calcule_BE : BE_t = Σ_sexe (assurés en vie en t) × a_{65+m}.
# LE SEUL TAUX qui entre ici est taux_actu_be = 2 %, par les annuités qu'on lui
# passe — jamais le 3 % de la tarification, jamais rdt_actifs. Le 🚩 de la
# checklist (« BOF₂₀₂₂ positif ») se vérifie sur cette ligne et sur elle seule.
calcule_BE <- function(n_h, a_h, n_f, a_f) n_h * a_h + n_f * a_f

# perte_annuelle : −ΔBOF actualisé, la variable dont le SCR est le quantile.
# BOF_{t+1} est ramené en t au taux de MARCHÉ taux_actu_be, pas au rendement
# des actifs : c'est une actualisation, pas une capitalisation.
perte_annuelle <- function(bof_t, bof_t1, taux) bof_t - bof_t1 / (1 + taux)

# calcule_SCR : VaR_{99,5 %} de la perte, en convention STATISTIQUE D'ORDRE
# (type = 1 = inverse de la fonction de répartition empirique, VaR_α(X) =
# inf{x : F(x) ≥ α}). C'est la définition du cours, et le SCR est alors une
# perte RÉELLEMENT OBSERVÉE dans l'échantillon (la 26ᵉ pire sur 5000), pas une
# valeur interpolée entre deux scénarios. Écart avec type = 7 : ORACLE 10.
calcule_SCR <- function(bof_t, bof_t1, taux)
  quantile(perte_annuelle(bof_t, bof_t1, taux), 0.995, type = 1, names = FALSE)

nom_sexe <- function(sexe) if (sexe == "h") "Hommes" else "Femmes"
nom_rdt  <- function(rdt)  paste0("rdt", round(100 * rdt))

# ==============================================================
# Bloc 3 — Structure de l'aléa : q₆₅, q₆₆, et les QUATRE tirages
# ==============================================================
# q₆₅ = 1 − ₁p₆₅ et q₆₆ = 1 − ₂p₆₅/₁p₆₅ : probabilités de décès des deux
# premières années, calculées PAR TRAJECTOIRE. Elles se révèlent dégénérées
# (ORACLES 1 et 1bis) — c'est un résultat de jumpchoice = "fit", pas une
# hypothèse — mais on les calcule en vecteurs, ce qui rend le code juste
# indépendamment de cette dégénérescence.
q65 <- lapply(S, function(M) unname(1 - M["1", ]))
q66 <- lapply(S, function(M) unname(1 - M["2", ] / M["1", ]))

cat("=== Bloc 3 — Structure de l'alea de l'annee 1 =================\n")
cat("Les taux 2022-2023 sont AJUSTES (le fit LC va jusqu'en 2023), donc\n")
cat("identiques sur les 5000 trajectoires : la mortalite REALISEE des deux\n")
cat("premieres annees ne porte que du risque IDIOSYNCRATIQUE (binomial).\n")
cat("Le risque SYSTEMATIQUE entre par la REEVALUATION du BE : a66(j) pour\n")
cat("SCR_2022, a67(j) pour SCR_2023.\n")
for (sexe in c("h", "f"))
  cat("  ", nom_sexe(sexe), " : q65 = ", sprintf("%.6f", q65[[sexe]][1]),
      " (etendue ", format(diff(range(q65[[sexe]])), digits = 3), ")",
      " | q66 = ", sprintf("%.6f", q66[[sexe]][1]),
      " (etendue ", format(diff(range(q66[[sexe]])), digits = 3), ")\n", sep = "")

# ---- PROTOCOLE RNG (contraignant) ----------------------------
# QUATRE rbinom, et quatre seulement, tirés ICI — avant toute boucle sur le
# produit et avant toute boucle sur le rendement. Trois raisons :
#   1. [A] et [B] sont deux produits vendus à la MÊME population, de MÊME
#      mortalité : mêmes tirages ⇒ l'écart SCR[A] − SCR[B] commenté en D.3 est
#      un effet PRODUIT pur, débarrassé du bruit de Monte-Carlo (± 14 sur [A],
#      cf. bloc 9). Même protocole d'appariement que la sensibilité COVID de C.
#   2. Les deux scénarios de rendement partagent les mêmes tirages ⇒ l'identité
#      de translation de l'ORACLE 6 est EXACTE, et non approchée à ± 14.
#   3. Deux rbinom SÉPARÉS pour H et F : un rbinom(n_sim, 500, c(qh, qf))
#      recyclerait les deux probabilités en alternance, silencieusement.
deces_2023 <- list(h = rbinom(n_sim, n_tetes, q65$h),
                   f = rbinom(n_sim, n_tetes, q65$f))
survivants_2023 <- lapply(deces_2023, function(d) n_tetes - d)

# État CENTRAL de 2023 pour le SCR₂₀₂₃ (D.4) : survivants ATTENDUS sous q₆₅.
# 500·(1−q₆₅) n'est pas entier alors que la binomiale de l'année suivante exige
# un effectif entier ; on arrondit au plus proche, seul traitement qui laisse un
# ÉTAT DE PORTEFEUILLE RÉEL (on ne peut pas avoir 492,96 assurés en vie).
# L'effet de cet arrondi est chiffré en forme fermée par l'ORACLE 8.
n_central <- c(h = round(n_tetes * (1 - q65$h[1])),
               f = round(n_tetes * (1 - q65$f[1])))
deces_2024 <- list(h = rbinom(n_sim, n_central[["h"]], q66$h),
                   f = rbinom(n_sim, n_central[["f"]], q66$f))
survivants_2024 <- list(h = n_central[["h"]] - deces_2024$h,
                        f = n_central[["f"]] - deces_2024$f)

cat("Etat central 2023 : H ", sprintf("%.2f", n_tetes * (1 - q65$h[1])),
    " -> ", n_central[["h"]], " | F ", sprintf("%.2f", n_tetes * (1 - q65$f[1])),
    " -> ", n_central[["f"]], "  (arrondis, cf. ORACLE 8)\n", sep = "")

# ---- Annuités, calculées UNE SEULE FOIS -----------------------
# Hors de toute boucle sur le rendement : le contrôle « BE inchangé sous 1 % »
# (ORACLE 7) devient ainsi STRUCTUREL, garanti par la forme du code et non par
# un test — même logique que la séparation base/covid du script 04.
#   m = 0 → a₆₅ (date 2022)  ·  m = 1 → a₆₆ (date 2023)  ·  m = 2 → a₆₇ (2024)
annuites <- list(A = list(), B = list())
for (produit in c("A", "B"))
  for (m in 0:2)
    annuites[[produit]][[as.character(m)]] <-
      lapply(S, function(M) annuite_conditionnelle(M, m, n_reste(produit, m),
                                                   taux_actu_be))

cat("Annuites a ", 100 * taux_actu_be, " %, terme echu (a65 / a66 / a67) :\n", sep = "")
for (produit in c("A", "B"))
  for (sexe in c("h", "f")) {
    moyennes <- sapply(as.character(0:2), function(m)
                       mean(annuites[[produit]][[m]][[sexe]]))
    termes   <- sapply(0:2, function(m) n_reste(produit, m))
    cat("  [", produit, "] ", nom_sexe(sexe), " : ",
        paste(sprintf("%.4f", moyennes), collapse = " / "),
        "   (", paste(termes, collapse = " / "), " termes)\n", sep = "")
  }

# ==============================================================
# Bloc 4 — D.1 à D.3 : BOF₂₀₂₂, BOF₂₀₂₃ et SCR₂₀₂₂
# ==============================================================
# scenario_solvabilite : la projection en t⁺ du portefeuille sur deux exercices,
# pour un produit et un rendement d'actifs donnés. Décès, survivants et annuités
# sont fixés en amont : changer de rendement ne consomme pas un seul nombre
# aléatoire et ne recalcule pas un seul BE.
scenario_solvabilite <- function(produit, rdt) {
  a <- annuites[[produit]]
  N_h <- survivants_2023$h ; N_f <- survivants_2023$f

  # --- D.1 : date 2022⁺, après encaissement des primes ---------
  # BE = ESPÉRANCE sur les 5000 trajectoires (définition du Best Estimate) :
  # BE₂₀₂₂ est donc un SCALAIRE, l'incertitude de trajectoire n'entre qu'en 2023.
  A_2022   <- n_tetes * sum(primes[[produit]])
  BE_2022  <- calcule_BE(n_tetes, mean(a[["0"]]$h), n_tetes, mean(a[["0"]]$f))
  BOF_2022 <- A_2022 - BE_2022

  # --- D.2 : date 2023⁺, par trajectoire -----------------------
  # Les actifs croissent une année pleine, PUIS la prestation k = 1 (rente de 1
  # par survivant, terme échu) est versée. Les décédés de l'année ne touchent
  # rien : c'est la définition même du terme échu.
  A_2023   <- A_2022 * (1 + rdt) - (N_h + N_f)
  BE_2023  <- calcule_BE(N_h, a[["1"]]$h, N_f, a[["1"]]$f)
  BOF_2023 <- A_2023 - BE_2023
  SCR_2022 <- calcule_SCR(BOF_2022, BOF_2023, taux_actu_be)

  # --- D.4 : état central 2023, puis date 2024⁺ ----------------
  # Lecture « état central » : l'état 2023 est l'état ATTENDU (q₆₅ étant commune
  # aux 5000 trajectoires, il n'y en a qu'un), et l'on distribue BOF₂₀₂₄ depuis
  # cet état. Un seul SCR₂₀₂₃, comparable à SCR₂₀₂₂, sans simulations imbriquées.
  A_2023c   <- A_2022 * (1 + rdt) - sum(n_central)
  BE_2023c  <- calcule_BE(n_central[["h"]], mean(a[["1"]]$h),
                          n_central[["f"]], mean(a[["1"]]$f))
  BOF_2023c <- A_2023c - BE_2023c
  Np_h <- survivants_2024$h ; Np_f <- survivants_2024$f
  A_2024   <- A_2023c * (1 + rdt) - (Np_h + Np_f)   # 2ᵉ année de rendement
  BE_2024  <- calcule_BE(Np_h, a[["2"]]$h, Np_f, a[["2"]]$f)
  BOF_2024 <- A_2024 - BE_2024
  SCR_2023 <- calcule_SCR(BOF_2023c, BOF_2024, taux_actu_be)

  list(A_2022 = A_2022, BE_2022 = BE_2022, BOF_2022 = BOF_2022,
       A_2023 = A_2023, BE_2023 = BE_2023, BOF_2023 = BOF_2023, SCR_2022 = SCR_2022,
       A_2023c = A_2023c, BE_2023c = BE_2023c, BOF_2023c = BOF_2023c,
       BE_2024 = BE_2024, BOF_2024 = BOF_2024, SCR_2023 = SCR_2023)
}

res <- list()
for (rdt in rdt_grille) {
  res[[nom_rdt(rdt)]] <- list()
  for (produit in c("A", "B"))
    res[[nom_rdt(rdt)]][[produit]] <- scenario_solvabilite(produit, rdt)
}
base <- res[[nom_rdt(rdt_actifs)]]   # scénario de référence, i = 4 %

cat("\n=== Bloc 4 — D.1 a D.3 : BOF et SCR_2022 (rdt actifs = 4 %) ===\n")
for (produit in c("A", "B")) {
  pf <- base[[produit]]
  cat("  [", produit, "] A_2022 = ", sprintf("%.2f", pf$A_2022),
      " | BE_2022 = ", sprintf("%.2f", pf$BE_2022),
      " | BOF_2022 = ", sprintf("%.2f", pf$BOF_2022), "\n", sep = "")
  cat("        BOF_2023 : moyenne = ", sprintf("%.2f", mean(pf$BOF_2023)),
      " | ecart-type = ", sprintf("%.2f", sd(pf$BOF_2023)),
      " | quantile 0.5 % = ",
      sprintf("%.2f", quantile(pf$BOF_2023, 0.005, type = 1, names = FALSE)),
      "\n", sep = "")
  cat("        SCR_2022 = ", sprintf("%.2f", pf$SCR_2022), "\n", sep = "")
}
# BOF₂₀₂₂ < 0 est le résultat ATTENDU, pas un bug : la prime a été tarifée à 3 %
# alors que le BE s'actualise à 2 %, donc BE₂₀₂₂ > A₂₀₂₂ mécaniquement.
cat("  D.3 — comparaison : SCR_2022[A] - SCR_2022[B] = ",
    sprintf("%.2f", base$A$SCR_2022 - base$B$SCR_2022), "\n", sep = "")

# ==============================================================
# Bloc 5 — D.4 : SCR₂₀₂₃ en lecture « état central »
# ==============================================================
cat("\n=== Bloc 5 — D.4 : SCR_2023 (etat central, rdt actifs = 4 %) ==\n")
for (produit in c("A", "B")) {
  pf <- base[[produit]]
  cat("  [", produit, "] A_2023c = ", sprintf("%.2f", pf$A_2023c),
      " | BE_2023c = ", sprintf("%.2f", pf$BE_2023c),
      " | BOF_2023c = ", sprintf("%.2f", pf$BOF_2023c), "\n", sep = "")
  cat("        BOF_2024 : moyenne = ", sprintf("%.2f", mean(pf$BOF_2024)),
      " | ecart-type = ", sprintf("%.2f", sd(pf$BOF_2024)), "\n", sep = "")
  cat("        SCR_2023 = ", sprintf("%.2f", pf$SCR_2023),
      " | SCR_2023 - SCR_2022 = ", sprintf("%+.2f", pf$SCR_2023 - pf$SCR_2022),
      "\n", sep = "")
}
cat("  /!\\ Cet ecart SCR_2023 - SCR_2022 doit etre lu avec l'erreur de\n")
cat("      Monte-Carlo du quantile (bloc 9.1) avant toute interpretation.\n")

# ==============================================================
# Bloc 6 — D.5 : sensibilité rdt_actifs = 1 %, et tableau des 8 SCR
# ==============================================================
# Seuls les ACTIFS changent. Les BE sont inchangés par construction (annuités
# calculées hors de la boucle) — l'ORACLE 7 l'imprime quand même.
cat("\n=== Bloc 6 — D.5 : les 8 SCR (2 portefeuilles x 2 dates x 2 rdt) ===\n")
cat("  Le SCR BRUT est la valeur de reference. La colonne « plancher » est la\n")
cat("  lecture prudentielle : sous Solvabilite II le capital requis est\n")
cat("  plancherise a 0 — un VaR negatif n'est PAS un besoin de capital negatif.\n")
cat("  rdt   portef.   date     SCR brut     plancher max(SCR, 0)\n")
for (rdt in rdt_grille)
  for (produit in c("A", "B"))
    for (date in c("2022", "2023")) {
      scr <- res[[nom_rdt(rdt)]][[produit]][[paste0("SCR_", date)]]
      cat(sprintf("  %2.0f %%    [%s]      %s   %10.2f   %10.2f\n",
                  100 * rdt, produit, date, scr, max(scr, 0)))
    }

# ==============================================================
# Bloc 7 — Diagnostics de structure
# ==============================================================
cat("\n=== Bloc 7 — Diagnostics ======================================\n")

# ---- 7.1 Décomposition de la variance de BOF₂₀₂₃ ---------------
# Loi de la variance totale, en CONDITIONNANT SUR LA TRAJECTOIRE :
#   var(BOF₂₀₂₃) = E[var(BOF | traj)]  +  var(E[BOF | traj])
#                  \___ idiosyncratique ___/  \___ systematique ___/
# En réécrivant BOF₂₀₂₃ = A₂₀₂₂(1+i) − N_h·ä₆₆ʰ − N_f·ä₆₆ᶠ, avec
# ä₆₆ = 1 + a₆₆ : la prestation de 1 versée en date 2023 PLUS la valeur du
# reliquat. Aucun double comptage — A₂₀₂₃ retranche le paiement, BE₂₀₂₃ ne le
# contient pas.
cat("\n7.1 — Variance de BOF_2023, par canal (rdt actifs = 4 %) :\n")
decomposition <- data.frame(empirique = c(NA, NA), idiosyncratique = c(NA, NA),
                            systematique = c(NA, NA), row.names = c("A", "B"))
for (produit in c("A", "B")) {
  adot_h  <- 1 + annuites[[produit]][["1"]]$h
  adot_f  <- 1 + annuites[[produit]][["1"]]$f
  # var(N) = n·q·(1−q) : les deux sexes sont tirés SÉPARÉMENT, donc les
  # variances CONDITIONNELLES s'additionnent.
  idio <- mean(adot_h^2 * n_tetes * q65$h * (1 - q65$h) +
               adot_f^2 * n_tetes * q65$f * (1 - q65$f))
  # var EMPIRIQUE de la somme : elle embarque déjà la covariance H/F des
  # trajectoires. Ne JAMAIS l'écrire en forme fermée comme somme de deux
  # variances, ce serait supposer une indépendance au lieu de la mesurer.
  syst <- var(-n_tetes * (1 - q65$h) * adot_h - n_tetes * (1 - q65$f) * adot_f)
  emp  <- var(base[[produit]]$BOF_2023)
  decomposition[produit, ] <- c(emp, idio, syst)
  cat("  [", produit, "] var empirique = ", sprintf("%.1f", emp),
      " | idiosyncratique = ", sprintf("%.1f", idio),
      " (", sprintf("%.1f %%", 100 * idio / (idio + syst)), ")",
      " | systematique = ", sprintf("%.1f", syst),
      " (", sprintf("%.1f %%", 100 * syst / (idio + syst)), ")\n", sep = "")
  cat("        correlation H/F des annuites a66 : ",
      sprintf("%+.4f", cor(adot_h, adot_f)), "\n", sep = "")
}

# ---- 7.2 Gain de spread contre queue de distribution -----------
# Identité EXACTE : SCR = queue − gain, avec
#   gain  = −E[perte] = A_t·(i − r)/(1 + r)   en FORME FERMÉE
#   queue = VaR_{99,5 %}(perte) − E[perte]
# La forme fermée tient parce que la récurrence a₆₅ = v·₁p₆₅·(1 + a₆₆) fait que
# (prestation + BE₂₀₂₃)/(1+r) a pour espérance BE₂₀₂₂ EXACTEMENT : le bruit de
# trajectoire se compense entre les deux BE, et le seul gain espéré est le
# spread i − r sur les actifs. C'est l'explication du SIGNE du SCR.
cat("\n7.2 — Signe du SCR : gain de spread contre queue a 99,5 % :\n")
cat("  rdt   portef. date    gain empirique  gain forme fermee   queue      SCR\n")
spread <- list()
for (rdt in rdt_grille)
  for (produit in c("A", "B"))
    for (date in c("2022", "2023")) {
      pf     <- res[[nom_rdt(rdt)]][[produit]]
      bof_t  <- if (date == "2022") pf$BOF_2022 else pf$BOF_2023c
      bof_t1 <- if (date == "2022") pf$BOF_2023 else pf$BOF_2024
      A_t    <- if (date == "2022") pf$A_2022   else pf$A_2023c
      perte  <- perte_annuelle(bof_t, bof_t1, taux_actu_be)
      gain_e <- -mean(perte)
      gain_f <- A_t * (rdt - taux_actu_be) / (1 + taux_actu_be)
      scr    <- pf[[paste0("SCR_", date)]]
      spread[[paste(nom_rdt(rdt), produit, date)]] <-
        c(gain_empirique = gain_e, gain_ferme = gain_f, queue = scr + gain_e, scr = scr)
      cat(sprintf("  %2.0f %%    [%s]   %s  %13.2f  %17.2f %8.2f %8.2f\n",
                  100 * rdt, produit, date, gain_e, gain_f, scr + gain_e, scr))
    }

# ==============================================================
# Bloc 8 — Oracles
# ==============================================================
cat("\n=== Bloc 8 — Oracles ==========================================\n")
ok <- function(test) if (test) "OK" else "*** ECHEC"
v_be <- 1 / (1 + taux_actu_be)

# ---- ORACLE 1 et 1bis : q₆₅ ET q₆₆ sont communes ---------------
# Le fit LC va jusqu'en 2023 : µ₆₅(2022) ET µ₆₆(2023) sont des valeurs AJUSTÉES,
# donc identiques sur les 5000 trajectoires. La première année stochastique est
# 2024, portée par µ₆₇ — donc par S["3", ]. C'est un RÉSULTAT VÉRIFIÉ de
# jumpchoice = "fit", jamais une hypothèse.
cat("\nORACLE 1/1bis — mortalite des 2 premieres annees, deterministe :\n")
for (sexe in c("h", "f")) {
  e1 <- diff(range(S[[sexe]]["1", ])) ; e2 <- diff(range(q66[[sexe]]))
  e3 <- diff(range(S[[sexe]]["3", ]))
  cat("  ", nom_sexe(sexe), " : etendue 1p65 = ", format(e1, digits = 3),
      " | etendue q66 = ", format(e2, digits = 3), "  ", ok(e1 == 0 && e2 == 0),
      "\n            etendue 3p65 = ", format(e3, digits = 3),
      " (> 0 : 2024 est la 1ere annee stochastique)  ", ok(e3 > 0), "\n", sep = "")
}

# ---- ORACLE 2 : ANCRAGE ABSOLU sur la section C ----------------
# a₆₅ reconstruite ici doit reproduire EXACTEMENT la VAP moyenne à 2 % du
# script 04 — même chemin de calcul, donc écart nul au bit près. C'est le seul
# contrôle qui pince l'HORIZON ABSOLU des bornes de sommation, et il le fait
# sans instancier le taux technique de la section C.
cat("\nORACLE 2 — ancrage absolu de a65 sur sensibilite[taux = 2 %] du script 04 :\n")
for (produit in c("A", "B"))
  for (sexe in c("h", "f")) {
    ecart <- abs(mean(annuites[[produit]][["0"]][[sexe]]) -
                 sens[sens$taux == taux_actu_be, paste0(produit, "_", sexe)])
    cat("  [", produit, "] ", nom_sexe(sexe), " : ecart = ",
        format(ecart, digits = 3), "  ", ok(ecart == 0), "\n", sep = "")
  }

# ---- ORACLE 3 : ANCRAGE RELATIF, la récurrence d'annuité -------
#   a₆₅ = v·₁p₆₅·(1 + a₆₆)   et   a₆₆ = v·₁p₆₆·(1 + a₆₇)
# Le « 1 + » est le paiement de la date suivante, le facteur v·ₖp l'escompte et
# la survie. Sur [B], c'est CE test qui attrape une annuité conditionnelle
# laissée à 15 termes au lieu de 14 (résidu ~5e-1 au lieu de ~1e-15).
# ANGLE MORT CONNU, et c'est pourquoi l'ORACLE 2 existe : la récurrence est
# invariante à un décalage COMMUN des deux bornes (a₆₅ à 16 termes contre a₆₆ à
# 15 la satisferait aussi).
cat("\nORACLE 3 — recurrence a_{65+m} = v * 1p_{65+m} * (1 + a_{66+m}) :\n")
for (produit in c("A", "B"))
  for (sexe in c("h", "f")) {
    a <- annuites[[produit]]
    p65 <- S[[sexe]]["1", ] ; p66 <- S[[sexe]]["2", ] / S[[sexe]]["1", ]
    e_65 <- max(abs(a[["0"]][[sexe]] - v_be * p65 * (1 + a[["1"]][[sexe]])))
    e_66 <- max(abs(a[["1"]][[sexe]] - v_be * p66 * (1 + a[["2"]][[sexe]])))
    cat("  [", produit, "] ", nom_sexe(sexe), " : max|a65 - ...| = ",
        format(e_65, digits = 3), " | max|a66 - ...| = ", format(e_66, digits = 3),
        "  ", ok(e_65 < 1e-12 && e_66 < 1e-12), "\n", sep = "")
  }

# ---- ORACLE 4 : la binomiale tire bien q₆₅ ---------------------
# Tolérance en unités d'ERREUR-TYPE, pas en absolu : un seuil dur clignoterait
# pour rien. SE = sqrt(q(1−q)/n_tetes/n_sim).
cat("\nORACLE 4 — frequence empirique des deces de l'annee 1 vs q65 :\n")
for (sexe in c("h", "f")) {
  q  <- q65[[sexe]][1]
  se <- sqrt(q * (1 - q) / n_tetes / n_sim)
  ecart <- abs(mean(deces_2023[[sexe]]) / n_tetes - q)
  cat("  ", nom_sexe(sexe), " : |frequence - q65| = ", format(ecart, digits = 3),
      " | 3 SE = ", format(3 * se, digits = 3), "  ", ok(ecart < 3 * se), "\n", sep = "")
}

# ---- ORACLE 5 : hierarchie des portefeuilles -------------------
# Sur les valeurs BRUTES. [A] porte 37 termes de risque de longévité contre 15,
# donc davantage de risque systématique. Le SIGNE de SCR[B] n'est PAS un oracle :
# c'est un résultat à chiffrer (bloc 7.2 en donne le mécanisme).
cat("\nORACLE 5 — SCR[A] > SCR[B] sur les valeurs BRUTES :\n")
for (rdt in rdt_grille)
  for (date in c("2022", "2023")) {
    sA <- res[[nom_rdt(rdt)]]$A[[paste0("SCR_", date)]]
    sB <- res[[nom_rdt(rdt)]]$B[[paste0("SCR_", date)]]
    cat("  rdt ", sprintf("%2.0f %%", 100 * rdt), ", ", date, " : ",
        sprintf("%.2f", sA), " > ", sprintf("%.2f", sB), "  ", ok(sA > sB),
        "  (SCR[A] > 0 : ", ok(sA > 0), ")\n", sep = "")
  }

# ---- ORACLE 6 : la sensibilité au rendement est une TRANSLATION -
# Changer i déplace BOF_{t+1} d'une constante identique sur les 5000
# trajectoires ; les quantiles étant équivariants par translation, le SCR se
# déplace EXACTEMENT de   [A_t(i)·(i − r) − A_t(i')·(i' − r)] / (1 + r).
# La constante DIFFÈRE entre 2022 et 2023 : en 2023 l'actif A_t dépend lui-même
# de i (les actifs ont capitalisé une fois) alors que la prestation de l'année 1
# n'en dépend pas. Testé sur les SCR BRUTS uniquement — le plancher à 0 casse
# la linéarité. Cette identité n'est exacte que sous nombres aléatoires communs.
cat("\nORACLE 6 — translation exacte du SCR sous rdt 4 % -> 1 % :\n")
for (produit in c("A", "B"))
  for (date in c("2022", "2023")) {
    pf_haut <- res[[nom_rdt(rdt_grille[1])]][[produit]]
    pf_bas  <- res[[nom_rdt(rdt_grille[2])]][[produit]]
    A_haut <- if (date == "2022") pf_haut$A_2022 else pf_haut$A_2023c
    A_bas  <- if (date == "2022") pf_bas$A_2022  else pf_bas$A_2023c
    obs <- pf_bas[[paste0("SCR_", date)]] - pf_haut[[paste0("SCR_", date)]]
    att <- (A_haut * (rdt_grille[1] - taux_actu_be) -
            A_bas  * (rdt_grille[2] - taux_actu_be)) / (1 + taux_actu_be)
    cat("  [", produit, "] ", date, " : observe = ", sprintf("%.4f", obs),
        " | attendu = ", sprintf("%.4f", att), " | ecart = ",
        format(abs(obs - att), digits = 3), "  ", ok(abs(obs - att) < 1e-9),
        "\n", sep = "")
  }

# ---- ORACLE 7 : les BE ne bougent pas avec le rendement --------
# Garanti STRUCTURELLEMENT (annuités et effectifs hors de la boucle), imprimé
# quand même : c'est le contrôle explicitement demandé en D.5.
cat("\nORACLE 7 — BE identiques entre rdt 4 % et rdt 1 % :\n")
for (produit in c("A", "B")) {
  pf_haut <- res[[nom_rdt(rdt_grille[1])]][[produit]]
  pf_bas  <- res[[nom_rdt(rdt_grille[2])]][[produit]]
  e <- max(abs(pf_haut$BE_2022  - pf_bas$BE_2022),
           abs(pf_haut$BE_2023c - pf_bas$BE_2023c),
           max(abs(pf_haut$BE_2023 - pf_bas$BE_2023)),
           max(abs(pf_haut$BE_2024 - pf_bas$BE_2024)))
  cat("  [", produit, "] max|BE(4 %) - BE(1 %)| aux 3 dates = ",
      format(e, digits = 3), "  ", ok(e == 0), "\n", sep = "")
}

# ---- ORACLE 8 : cohérence entre D.2 et D.4 ---------------------
# BOF₂₀₂₃ᶜ (état central, D.4) doit différer de mean(BOF₂₀₂₃) (D.2) du SEUL
# effet d'arrondi des effectifs, calculable en forme fermée. Un écart résiduel
# de plusieurs dizaines signalerait une incohérence entre les deux blocs.
cat("\nORACLE 8 — coherence D.2 <-> D.4, l'ecart est l'arrondi des effectifs :\n")
for (produit in c("A", "B")) {
  pf   <- base[[produit]]
  adot <- sapply(c("h", "f"), function(s) mean(1 + annuites[[produit]][["1"]][[s]]))
  att  <- -sum((n_central - n_tetes * c(h = 1 - q65$h[1], f = 1 - q65$f[1])) * adot)
  obs  <- pf$BOF_2023c - mean(pf$BOF_2023)
  cat("  [", produit, "] observe = ", sprintf("%+.2f", obs),
      " | attendu (arrondi seul) = ", sprintf("%+.2f", att),
      " | residuel = ", sprintf("%+.2f", obs - att), "  (erreur MC)\n", sep = "")
}

# ---- ORACLE 8bis : l'identité de spread, aux deux dates --------
# Tolérance ±2, pas 1e-12 : le gain empirique porte l'erreur de Monte-Carlo des
# binomiales (~0,8).
cat("\nORACLE 8bis — gain empirique = A_t (i - r)/(1 + r), tolerance +/- 2 :\n")
for (cle in names(spread)) {
  s <- spread[[cle]]
  cat("  ", cle, " : ecart = ", sprintf("%+.3f", s["gain_empirique"] - s["gain_ferme"]),
      "  ", ok(abs(s["gain_empirique"] - s["gain_ferme"]) < 2), "\n", sep = "")
}

# ---- ORACLE 9 : la décomposition de variance se referme --------
# Tolérance 2 % : les deux membres sont estimés sur les MÊMES 5000 tirages, et
# l'erreur relative d'une variance empirique vaut sqrt(2/(n_sim-1)) ~ 2 %.
cat("\nORACLE 9 — somme des deux canaux vs variance empirique (tolerance 2 %) :\n")
for (produit in c("A", "B")) {
  d <- decomposition[produit, ]
  ec <- abs((d$idiosyncratique + d$systematique) / d$empirique - 1)
  cat("  [", produit, "] somme = ", sprintf("%.1f", d$idiosyncratique + d$systematique),
      " | empirique = ", sprintf("%.1f", d$empirique),
      " | ecart relatif = ", sprintf("%.2f %%", 100 * ec), "  ", ok(ec < 0.02),
      "\n", sep = "")
}

# ---- ORACLE 10 : la convention de quantile est immatérielle ----
# Imprimé UNE FOIS. À comparer à l'erreur de Monte-Carlo du quantile lui-même
# (bloc 9.1) : le choix de convention est du bruit devant elle.
perte_A <- perte_annuelle(base$A$BOF_2022, base$A$BOF_2023, taux_actu_be)
cat("\nORACLE 10 — convention de quantile, SCR_2022[A] :\n")
cat("  type = 1 (statistique d'ordre, retenue) = ",
    sprintf("%.4f", quantile(perte_A, 0.995, type = 1, names = FALSE)),
    "\n  type = 7 (defaut de R, interpolation)   = ",
    sprintf("%.4f", quantile(perte_A, 0.995, type = 7, names = FALSE)),
    "\n  ecart = ", sprintf("%.4f", quantile(perte_A, 0.995, type = 7, names = FALSE) -
                              quantile(perte_A, 0.995, type = 1, names = FALSE)),
    "\n", sep = "")

# ==============================================================
# Bloc 9 — Limites CHIFFRÉES (distinct des oracles)
# ==============================================================
# Le SCR du projet reste celui de la formule VaR du brief, calculé au bloc 4.
# Ce qui suit sont des BORNES et des marges d'erreur — jamais des scénarios
# alternatifs. Elles sont imprimées ici parce qu'elles seront citées CHIFFRÉES
# au rapport, et qu'un chiffre du rapport doit venir d'une console réelle.
cat("\n=== Bloc 9 — Limites chiffrees ================================\n")
limites <- list()

# ---- 9.1 Erreur de Monte-Carlo du quantile 99,5 % --------------
# Le q₉₉,₅ de 5000 points est la ~25ᵉ plus grande valeur : son incertitude
# d'échantillonnage n'est pas négligeable. Bootstrap du vecteur de perte.
# Ce tirage vient EN DERNIER, après toutes les quantités substantielles : il ne
# peut rien perturber en amont du flux aléatoire.
cat("\n9.1 — Erreur de Monte-Carlo du SCR (bootstrap, B = 1000) :\n")
for (produit in c("A", "B")) {
  pf <- base[[produit]]
  perte <- perte_annuelle(pf$BOF_2022, pf$BOF_2023, taux_actu_be)
  boot  <- replicate(1000, quantile(sample(perte, n_sim, replace = TRUE), 0.995,
                                    type = 1, names = FALSE))
  limites[[paste0("mc_", produit)]] <- sd(boot)
  cat("  [", produit, "] SCR_2022 = ", sprintf("%.2f", pf$SCR_2022),
      " | ecart-type bootstrap = ", sprintf("%.2f", sd(boot)),
      " -> a annoncer ", sprintf("%.0f +/- %.0f", pf$SCR_2022, 2 * sd(boot)),
      "\n", sep = "")
  cat("        SCR_2023 - SCR_2022 = ", sprintf("%+.2f", pf$SCR_2023 - pf$SCR_2022),
      if (abs(pf$SCR_2023 - pf$SCR_2022) < 2 * sd(boot))
        "  -> DANS LE BRUIT, ne pas interpreter" else
        "  -> hors bruit, interpretable", "\n", sep = "")
}

# ---- 9.2 Comonotonie H/F : borne haute --------------------------
# Les trajectoires H et F sont deux tirages SÉQUENTIELS du même flux (protocole
# RNG du script 03) : leur corrélation vaut ~0 PAR CONSTRUCTION. Le SCR encaisse
# donc un crédit de diversification longévité qui est un ARTEFACT DE SIMULATION,
# pas un fait modélisé. On en calcule la borne haute par appariement par rang.
# Portée : SCR₂₀₂₂ seulement. L'étendre à 2023 demanderait de réordonner aussi
# a₆₇ᶠ avec le même vecteur de rangs, ce qu'on ne fait pas.
cat("\n9.2 — Comonotonie H/F, borne haute sur SCR_2022 :\n")
for (produit in c("A", "B")) {
  pf <- base[[produit]]
  a66_h <- annuites[[produit]][["1"]]$h
  a66_f <- annuites[[produit]][["1"]]$f
  # ⚠️ a66_f[rank(a66_h)] serait FAUX : cela permute a66_f sans l'APPARIER (les
  # valeurs restent en ordre de simulation) et la borne sortirait égale au cas
  # de base, silencieusement vide. Il faut TRIER avant d'indexer par le rang :
  # la j-ième plus grande annuité H reçoit la j-ième plus grande annuité F.
  a66_f_como <- sort(a66_f)[rank(a66_h, ties.method = "first")]
  BOF_como <- pf$A_2023 - calcule_BE(survivants_2023$h, a66_h,
                                     survivants_2023$f, a66_f_como)
  scr_como <- calcule_SCR(pf$BOF_2022, BOF_como, taux_actu_be)
  limites[[paste0("como_", produit)]] <- scr_como
  cat("  [", produit, "] correlation a66 H/F : base = ", sprintf("%+.4f", cor(a66_h, a66_f)),
      " -> comonotone = ", sprintf("%+.4f", cor(a66_h, a66_f_como)), "\n", sep = "")
  cat("        permutation pure : |mean(a66_f_como) - mean(a66_f)| = ",
      format(abs(mean(a66_f_como) - mean(a66_f)), digits = 3),
      " -> BE_2022 et BOF_2022 invariants\n", sep = "")
  # Écart ABSOLU : une variation relative sur une base négative s'inverse et se
  # lit à contresens (sur [B] le SCR MONTE alors que le rapport baisse).
  cat("        SCR_2022 : base = ", sprintf("%.2f", pf$SCR_2022),
      " | comonotone = ", sprintf("%.2f", scr_como),
      " | ecart = ", sprintf("%+.2f", scr_como - pf$SCR_2022),
      if (pf$SCR_2022 > 0) sprintf(" (%+.1f %%)", 100 * (scr_como / pf$SCR_2022 - 1))
      else "  (pas de % : base negative)", "\n", sep = "")
}
# Limite jumelle, gratuite à énoncer : q₆₅ étant constante, le NOMBRE de décès
# est indépendant de la réévaluation du BE. Le modèle ne peut pas produire la
# dépendance réelle « moins de décès qu'attendu ⟺ améliorations futures
# meilleures qu'attendu ».

# ---- 9.3 Les deux définitions du SCR du brief -------------------
# Le brief donne AUSSI P(BOF_{t+1} > 0 | BOF_t = SCR_t) = 99,5 %, qui suppose
# A_t = BE_t + SCR_t (point fixe), alors que la formule VaR retenue utilise
# l'actif RÉEL. Les deux ne coïncident que si A_t = BE_t + SCR_t, très loin du
# cas présent. Solution du point fixe, en forme fermée :
#   SCR^pf = VaR_{99,5 %}(prestation + BE_{t+1}) / (1 + i) − BE_t
# soit l'actif nécessaire aujourd'hui pour couvrir le décaissement et le passif
# du centile 0,5 %, moins le BE déjà provisionné.
cat("\n9.3 — Lecture « point fixe » du brief, NON retenue :\n")
for (produit in c("A", "B")) {
  pf <- base[[produit]]
  sortie <- (survivants_2023$h + survivants_2023$f) + pf$BE_2023
  scr_pf <- quantile(sortie, 0.995, type = 1, names = FALSE) /
              (1 + rdt_actifs) - pf$BE_2022
  limites[[paste0("pf_", produit)]] <- scr_pf
  cat("  [", produit, "] SCR_2022 retenu (VaR) = ", sprintf("%.2f", pf$SCR_2022),
      " | lecture point fixe = ", sprintf("%.2f", scr_pf),
      " (ecart ", sprintf("%+.1f %%", 100 * (scr_pf / pf$SCR_2022 - 1)), ")\n", sep = "")
}

# ---- 9.4 BE_{t+1} path-wise (énoncé, rien à coder) --------------
# a₆₆(j) et a₆₇(j) sont calculées le long de TOUTE la trajectoire réalisée :
# l'assureur est supposé apprendre l'intégralité du futur en t+1. Le vrai BE
# serait E[· | F_{t+1}], de variance strictement plus faible. Notre SCR est donc
# CONSERVATEUR sur sa composante systématique — c'est ce qui produit les parts
# du bloc 7.1 — et c'est asymétrique avec BE₂₀₂₂, qui est bien une espérance.

# ==============================================================
# Bloc 10 — Figures
# ==============================================================
# Titres en ASCII : le device PNG sous Windows ne rend pas les indices Unicode
# (même remarque qu'aux scripts 02, 03 et 04).
png("figures/fig_D_hist_bof2023.png", width = 2000, height = 900, res = 150)
par(mfrow = c(1, 2))
for (produit in c("A", "B")) {
  pf <- base[[produit]]
  q005 <- quantile(pf$BOF_2023, 0.005, type = 1, names = FALSE)
  hist(pf$BOF_2023, breaks = 40, xlim = range(c(pf$BOF_2023, pf$BOF_2022)),
       col = adjustcolor("steelblue", alpha.f = 0.55), border = "white",
       xlab = "BOF 2023 (par trajectoire)", ylab = "Frequence",
       main = paste0("Portefeuille [", produit, "] — lecture du SCR 2022"))
  abline(v = pf$BOF_2022, lwd = 2, col = "firebrick")
  abline(v = q005, lwd = 2, lty = 2, col = "grey20")
  # Légende à droite et sur fond opaque : les deux traits verticaux tombent dans
  # la moitié gauche et la traverseraient.
  legend("topright", c("BOF 2022", "quantile 0.5 % de BOF 2023"),
         col = c("firebrick", "grey20"), lwd = 2, lty = c(1, 2), cex = 0.7,
         bg = "white", box.col = "grey70")
}
par(mfrow = c(1, 1))
dev.off()

png("figures/fig_D_scr_comparaison.png", width = 1500, height = 1000, res = 150)
# Colonne et étiquette ajoutées DANS LA MÊME itération : une matrice remplie
# d'un côté et des colnames écrits à la main de l'autre se désynchronisent sans
# rien signaler (bug constaté et corrigé ici — les colonnes 2 et 3 étaient
# interverties, « [B] i=4 % » affichait le SCR de [A] a 1 %).
mat <- NULL
etiquettes <- NULL
for (rdt in rdt_grille)
  for (produit in c("A", "B")) {
    pf <- res[[nom_rdt(rdt)]][[produit]]
    mat <- cbind(mat, c(pf$SCR_2022, pf$SCR_2023))
    etiquettes <- c(etiquettes, sprintf("[%s] i=%.0f%%", produit, 100 * rdt))
  }
dimnames(mat) <- list(c("SCR 2022", "SCR 2023"), etiquettes)
# Axe NON tronqué et zero visible : le basculement de signe sous i = 1 % doit
# sauter aux yeux, barres traversantes comprises.
etendue <- range(c(0, mat))
barplot(mat, beside = TRUE, ylim = etendue + c(-0.12, 0.12) * diff(etendue),
        col = c("steelblue", "firebrick"), ylab = "SCR brut (VaR 99,5 %)",
        main = "D.5 — SCR brut par portefeuille, date et rendement d'actifs")
abline(h = 0, lwd = 2)
legend("topleft", rownames(mat), fill = c("steelblue", "firebrick"),
       cex = 0.8, bty = "n")
dev.off()

# ==============================================================
# Bloc 11 — Sauvegardes
# ==============================================================
saveRDS(list(
  rdt_4pc = res$rdt4, rdt_1pc = res$rdt1,
  deces = list(annee_2023 = deces_2023, annee_2024 = deces_2024),
  survivants = list(annee_2023 = survivants_2023, annee_2024 = survivants_2024),
  n_central = n_central, annuites = annuites,
  decomposition = decomposition, spread = spread, limites = limites,
  taux_actu_be = taux_actu_be, rdt_grille = rdt_grille, n_tetes = n_tetes,
  n_sim = n_sim), file = "resultats/05_scr.rds")

cat("\n=== Bloc 11 — sauvegardes =====================================\n")
cat("resultats/05_scr.rds ecrit.\n")
cat("  $rdt_4pc / $rdt_1pc : par portefeuille, BOF_2022 scalaire, vecteurs\n")
cat("                        BOF_2023 et BOF_2024 (", n_sim, "), SCR_2022, SCR_2023\n")
cat("  $deces, $survivants : les 4 tirages binomiaux (nombres aleatoires communs)\n")
cat("  $decomposition, $spread, $limites : diagnostics des blocs 7 et 9\n")
cat("Figures : figures/fig_D_hist_bof2023.png, figures/fig_D_scr_comparaison.png\n")
