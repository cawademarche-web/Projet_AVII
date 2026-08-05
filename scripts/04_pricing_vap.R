# ==============================================================
# Script 04 — Section C : valeur actuarielle présente (VAP) des rentes
#             [A] viagère et [B] temporaire 15 ans, et primes uniques
# Projet ACTU-F502 — Assurance Vie II
#
#   ₖp₆₅ = exp(−Σ_{j=0}^{k−1} µ₆₅₊ⱼ)         (force constante par morceaux, A.1)
#   a₆₅  = Σ_{k≥1} vᵏ · ₖp₆₅ ,  v = 1/(1+t)   (rente à TERME ÉCHU)
#
# C.1  bloc 2 — VAP des rentes [A] et [B] par trajectoire, taux technique 3 %
# C.2  bloc 2 — moyenne, variance et écart-type sur les 5000 trajectoires
# C.3  bloc 3 — sensibilité au taux technique t = 1, 2, 3, 4, 5 %
# C.4  bloc 4 — primes uniques par le principe d'équivalence
#      bloc 5 — sensibilité prudentielle COVID (hors brief, décision actée)
#      bloc 6 — oracles
#
# Code 100 % CUSTOM : StMoMo n'a aucune fonction de tarification (décision
# méthodologique verrouillée, CLAUDE.md). Base R uniquement, aucun library().
#
# PAS de set.seed ici, et c'est un SIGNAL, pas un oubli : ce script ne tire
# rien. Les 5000 trajectoires sont lues sur disque (contrat d'interface). Un
# set.seed signalerait une re-simulation, donc une violation de ce contrat.
# ==============================================================

taux_tarif <- 0.03   # taux technique du brief (C.2, C.4). SEUL taux de ce
                     # script : ni taux_actu_be (2 %, section D) ni rdt_actifs
                     # (4 %, section D) n'y ont leur place.
duree_B    <- 15     # rente [B] : temporaire 15 ans, à terme échu

# ==============================================================
# Bloc 1 — Chargement (contrat d'interface) et fonctions
# ==============================================================
# Entrée UNIQUE de la tarification : les diagonales de cohorte du script 03,
# µ_{65+k}(2022+k) pour k = 0…36 — matrices 37 × 5000, trajectoires en
# COLONNES, lignes nommées par l'ÂGE "65".."101". Les arrays complets
# LCsimK_*.rds (80 Mo chacun) ne sont jamais chargés : la tarification n'a
# besoin que de cette diagonale.
diagonales <- readRDS("resultats/03_diagonales_cohorte.rds")
n_sim      <- diagonales$n_sim
n_A        <- length(diagonales$ages)   # 37 termes : k = 1…37, l'âge 102 est
                                        # atteint en k = 37 (clôture ℓ_102, A.1)

# probas_survie_trajectoire : ₖp₆₅ = exp(−Σ_{j=0}^{k−1} µ₆₅₊ⱼ), force de
# mortalité constante par morceaux — MÊME convention qu'en A.1 et que
# esperance_vie_cohorte() des scripts 02/03. Aucune convention nouvelle.
# apply(M, 2, ...) appelle cumprod SUR CHAQUE COLONNE séparément : la
# récurrence repart de 1 à chaque trajectoire par CONSTRUCTION. Le piège évité
# est cumprod(exp(-M)), qui lirait M comme un seul vecteur de 185 000 éléments
# en ordre colonne-majeur et courrait à travers les 5000 trajectoires sans
# jamais se réinitialiser.
probas_survie_trajectoire <- function(M) {
  S <- apply(M, 2, function(mu) cumprod(exp(-mu)))
  # M est indexée par ÂGE ; S est indexée par k, et sa ligne 1 vaut ₁p₆₅ =
  # survie jusqu'à 66 ans. Sans ce renommage, S hériterait des noms d'âge de M,
  # FAUX d'un an — et notre règle d'indexation par nom les croirait.
  #   k  →  âge atteint 65 + k  →  paiement en fin d'année 2022 + k
  rownames(S) <- as.character(1:nrow(S))
  attr(S, "age_atteint") <- 65 + 1:nrow(S)   # 66:102 ; le script 05 l'imprimera
  S                                          # en contrôle, il n'indexe pas avec
}

# vap_rente : a₆₅ à TERME ÉCHU = Σ_{k=1}^{n} vᵏ · ₖp₆₅, avec v = 1/(1+taux).
# (à échoir, ce serait ä₆₅ = a₆₅ + 1 : le terme k = 0 est absent ici, ni dans
# S — que cumprod ne produit pas — ni dans v^(1:n), qui démarre à v¹.)
# v^(1:n) est recyclé le long de chaque colonne en ordre colonne-majeur :
# alignement exact, n étant la hauteur de la sous-matrice. Aucune boucle sur
# les 5000 trajectoires.
vap_rente <- function(S, taux, n) colSums((1 / (1 + taux))^(1:n) * S[1:n, ])

# resume_vap : les statistiques exigées en C.2 (le brief demande la VARIANCE ;
# l'écart-type sert la discussion du rapport) et les quantiles de l'oracle 6.
resume_vap <- function(v) c(moyenne = mean(v), variance = var(v), sd = sd(v),
                            quantile(v, c(0.025, 0.5, 0.975)))

nom_sexe <- function(sexe) if (sexe == "h") "Hommes" else "Femmes"

# ==============================================================
# Bloc 2 — C.1 et C.2 : VAP des rentes [A] et [B] au taux technique de 3 %
# ==============================================================
# [A] rente viagère à terme échu : Σ_{k=1}^{37} vᵏ ₖp₆₅
# [B] rente temporaire 15 ans    : Σ_{k=1}^{15} vᵏ ₖp₆₅
# Les matrices de survie sont construites UNE SEULE FOIS : la sensibilité du
# bloc 3 ne recalculera que les facteurs d'actualisation vᵏ.
S_base <- lapply(diagonales$base, probas_survie_trajectoire)

vap_base <- list(
  A = lapply(S_base, vap_rente, taux = taux_tarif, n = n_A),
  B = lapply(S_base, vap_rente, taux = taux_tarif, n = duree_B))

cat("=== Bloc 2 — C.1/C.2 : VAP a terme echu, taux technique 3 % ===\n")
cat("Rente annuelle de 1, cohorte des 65 ans en 2022, N =", n_sim,
    "trajectoires.\n")
cat("[A] viagere : k = 1..", n_A, " | [B] temporaire : k = 1..", duree_B,
    "\n", sep = "")
for (produit in c("A", "B")) {
  for (sexe in c("h", "f")) {
    r <- resume_vap(vap_base[[produit]][[sexe]])
    cat("  [", produit, "] ", nom_sexe(sexe),
        " : moyenne = ", sprintf("%.4f", r["moyenne"]),
        " | variance = ", sprintf("%.6f", r["variance"]),
        " | ecart-type = ", sprintf("%.6f", r["sd"]), "\n", sep = "")
    cat("             quantiles 2.5 / 50 / 97.5 % : ",
        sprintf("%.4f", r["2.5%"]), " / ", sprintf("%.4f", r["50%"]), " / ",
        sprintf("%.4f", r["97.5%"]), "\n", sep = "")
  }
}
# Ce que cette variance mesure — et ce qu'elle ne mesure PAS. Les 5000 VAP
# diffèrent par la seule trajectoire de κ : c'est l'incertitude SYSTÉMATIQUE de
# longévité, NON DIVERSIFIABLE (1000 assurés ne la réduisent pas). Le risque
# idiosyncratique des têtes individuelles n'est pas ici : il apparaîtra en
# section D, par la binomiale des décès de l'année 1.

png("figures/fig_C_hist_vap.png", width = 2000, height = 1400, res = 150)
par(mfrow = c(2, 2))
for (produit in c("A", "B")) {
  for (sexe in c("h", "f")) {
    v <- vap_base[[produit]][[sexe]]
    # Titres en ASCII : le device PNG sous Windows ne rend pas les indices
    # Unicode (même remarque qu'aux scripts 02 et 03).
    hist(v, breaks = 40, col = adjustcolor(if (sexe == "h") "steelblue" else
                                             "firebrick", alpha.f = 0.55),
         border = "white", xlab = "VAP (rente annuelle de 1)",
         ylab = "Fréquence",
         main = paste0(nom_sexe(sexe), " — rente [", produit, "], t = 3 %"))
    abline(v = mean(v), lwd = 2, col = "grey20")
    abline(v = quantile(v, c(0.025, 0.975)), lwd = 2, lty = 2, col = "grey20")
    legend("topright", c("moyenne", "quantiles 2.5 / 97.5 %"),
           col = "grey20", lwd = 2, lty = c(1, 2), cex = 0.75)
  }
}
par(mfrow = c(1, 1))
dev.off()

# ==============================================================
# Bloc 3 — C.3 : sensibilité au taux technique
# ==============================================================
# Étude de SCÉNARIOS déterministes de taux, pas une modélisation du risque de
# taux. Le brief demande t = 1, 2, 4, 5 % ; 3 % est conservé en colonne
# centrale pour la lecture. Le taux 0 % n'a aucun sens tarifaire : il sert
# l'oracle 3 du bloc 6, où VAP[A] doit valoir e₆₅ − ½.
# S_base n'est PAS reconstruite : seuls les vᵏ changent d'un taux à l'autre.
taux_grille <- c(0, 0.01, 0.02, 0.03, 0.04, 0.05)

sensibilite <- data.frame(taux = taux_grille)
for (produit in c("A", "B")) {
  n <- if (produit == "A") n_A else duree_B
  for (sexe in c("h", "f")) {
    sensibilite[[paste0(produit, "_", sexe)]] <-
      sapply(taux_grille, function(taux) mean(vap_rente(S_base[[sexe]], taux, n)))
  }
}

cat("\n=== Bloc 3 — C.3 : VAP moyenne selon le taux technique ========\n")
cat("  taux        [A] H      [A] F      [B] H      [B] F\n")
for (i in seq_along(taux_grille)) {
  cat(sprintf("  %4.0f %%   %8.4f   %8.4f   %8.4f   %8.4f%s\n",
              100 * sensibilite$taux[i], sensibilite$A_h[i], sensibilite$A_f[i],
              sensibilite$B_h[i], sensibilite$B_f[i],
              if (sensibilite$taux[i] == 0) "   <- controle, pas un scenario" else
                if (sensibilite$taux[i] == taux_tarif) "   <- taux du brief" else ""))
}

png("figures/fig_C_sensibilite_taux.png", width = 1500, height = 1000, res = 150)
tarif <- sensibilite[sensibilite$taux > 0, ]   # 0 % exclu : contrôle, pas un taux
plot(100 * tarif$taux, tarif$A_f, type = "n", ylim = range(tarif[, -1]),
     xlab = "Taux technique t (%)", ylab = "VAP moyenne (rente annuelle de 1)",
     main = "C.3 — VAP moyenne vs taux technique")
for (produit in c("A", "B")) {
  for (sexe in c("h", "f")) {
    lines(100 * tarif$taux, tarif[[paste0(produit, "_", sexe)]], lwd = 2,
          lty = if (produit == "A") 1 else 2,
          col = if (sexe == "h") "steelblue" else "firebrick")
  }
}
abline(v = 100 * taux_tarif, lwd = 1, lty = 3, col = "grey40")
legend("topright", c("[A] Hommes", "[A] Femmes", "[B] Hommes", "[B] Femmes"),
       col = c("steelblue", "firebrick", "steelblue", "firebrick"),
       lty = c(1, 1, 2, 2), lwd = 2, cex = 0.8)
dev.off()

# ==============================================================
# Bloc 4 — C.4 : primes uniques par le principe d'équivalence
# ==============================================================
# Principe d'équivalence : prime unique = E[engagements de l'assureur], ici
# l'espérance SOUS MORTALITÉ STOCHASTIQUE, c'est-à-dire la moyenne des VAP sur
# les 5000 trajectoires de κ. Prime PURE : ni chargement, ni marge de
# prudence — une prime prudente serait un quantile, pas une moyenne.
# Rente annuelle de 1 : l'agrégation du portefeuille (500 H + 500 F) est une
# hypothèse de la SECTION D, elle ne descend pas dans un script de tarification.
primes <- lapply(vap_base, function(p) sapply(p, mean))

cat("\n=== Bloc 4 — C.4 : primes uniques (equivalence, t = 3 %) ======\n")
for (produit in c("A", "B")) {
  cat("  [", produit, "] Hommes = ", sprintf("%.4f", primes[[produit]]["h"]),
      " | Femmes = ", sprintf("%.4f", primes[[produit]]["f"]),
      "   (rente annuelle de 1)\n", sep = "")
}

# ==============================================================
# Bloc 5 — Sensibilité prudentielle COVID (hors brief)
# ==============================================================
# Mêmes calculs sur la calibration 1970-2019 (dérive plus forte, sens PRUDENT
# pour un assureur de rentes). À 3 % SEULEMENT : pas de sensibilité croisée
# COVID × taux.
# Les deux jeux de trajectoires sont issus du MÊME flux aléatoire (protocole
# RNG du script 03 : set.seed(1234) avant chaque couple, ordre H puis F), donc
# APPARIÉS par nombres aléatoires communs. L'écart mesuré ici est un effet de
# CALIBRATION PUR, sans bruit de Monte-Carlo différentiel — ce que prouve
# l'effondrement de sd(ΔVAP) imprimé ci-dessous.
S_covid <- lapply(diagonales$covid2019, probas_survie_trajectoire)
vap_covid <- list(
  A = lapply(S_covid, vap_rente, taux = taux_tarif, n = n_A),
  B = lapply(S_covid, vap_rente, taux = taux_tarif, n = duree_B))

cat("\n=== Bloc 5 — Sensibilite COVID (1970-2019), t = 3 % ===========\n")
ecart_covid <- list()
for (produit in c("A", "B")) {
  ecart_covid[[produit]] <- c(h = NA, f = NA)
  for (sexe in c("h", "f")) {
    vb <- vap_base[[produit]][[sexe]]
    vc <- vap_covid[[produit]][[sexe]]
    ecart_covid[[produit]][sexe] <- mean(vc) / mean(vb) - 1
    cat("  [", produit, "] ", nom_sexe(sexe),
        " : base = ", sprintf("%.4f", mean(vb)),
        " | COVID = ", sprintf("%.4f", mean(vc)),
        " | ecart = ", sprintf("%+.2f %%", 100 * ecart_covid[[produit]][sexe]),
        "\n", sep = "")
    cat("             appariement : sd(base) = ", sprintf("%.5f", sd(vb)),
        " | sd(COVID - base) = ", sprintf("%.5f", sd(vc - vb)),
        " (rapport ", sprintf("%.2f", sd(vc - vb) / sd(vb)), ")\n", sep = "")
  }
}
# La variante 1970-2019 mêle DEUX effets — dérive plus forte ET amorce sans
# surmortalité COVID (µ₆₅(2022) médian 0,014179 → 0,012477 chez les hommes).
# À énoncer au rapport, PAS à corriger par un recalage d'amorce.

# ==============================================================
# Bloc 6 — Oracles
# ==============================================================
cat("\n=== Bloc 6 — Oracles ==========================================\n")

# ---- Intégrité et ALIGNEMENT de la matrice de survie -----------
# Contrôle décisif de la vectorisation par colonne : la ligne 1 de S doit valoir
# exp(−µ₆₅(2022)) sur CHACUNE des 5000 colonnes. Si un cumprod avait franchi une
# frontière de colonne, S["1", j] traînerait le produit accumulé des colonnes
# précédentes pour tout j ≥ 2.
cat("\nIntegrite et alignement de S (kp65) :\n")
for (sexe in c("h", "f")) {
  S <- S_base[[sexe]]
  M <- diagonales$base[[sexe]]
  decroissante <- all(apply(S, 2, function(col) all(diff(col) < 0)))
  cat("  ", nom_sexe(sexe), " : max|S[\"1\", ] - exp(-mu[\"65\", ])| = ",
      format(max(abs(S["1", ] - exp(-M["65", ]))), digits = 3),
      " | S dans ]0,1] : ", all(S > 0 & S <= 1),
      " | decroissante en k : ", decroissante,
      " | NA : ", sum(is.na(S)), "\n", sep = "")
  cat("            rownames = \"", rownames(S)[1], "\"..\"", rownames(S)[nrow(S)],
      "\" | attr age_atteint = ", min(attr(S, "age_atteint")), ":",
      max(attr(S, "age_atteint")), "\n", sep = "")
}

# ---- ORACLE 1 : VAP[B] < VAP[A] sur CHACUNE des 5000 trajectoires
# [B] tronque la somme de [A] et tous les termes vᵏ·ₖp₆₅ sont > 0 : l'inégalité
# est structurelle AU SEIN d'une trajectoire. Une seule violation = bug
# d'indexation, pas un aléa.
cat("\nORACLE 1 — VAP[B] < VAP[A] sur chacune des", n_sim, "trajectoires :\n")
for (sexe in c("h", "f")) {
  violations <- sum(vap_base$B[[sexe]] >= vap_base$A[[sexe]])
  cat("  ", nom_sexe(sexe), " : violations = ", violations, "  ",
      if (violations == 0) "OK" else
        "*** BUG CERTAIN : indexation de la borne de sommation", "\n", sep = "")
}

# ---- ORACLE 2 : ordre de grandeur, bande recalibrée à TERME ÉCHU
# La bande 15-18 du kit initial valait à ÉCHOIR (ä = a + 1) : elle ne s'applique
# pas ici. À terme échu, on attend 13-14 (H) et 14-16 (F).
cat("\nORACLE 2 — VAP[A] moyenne a 3 %, bande recalibree a terme echu :\n")
bandes <- list(h = c(13, 14), f = c(14, 16))
for (sexe in c("h", "f")) {
  m <- mean(vap_base$A[[sexe]])
  cat("  ", nom_sexe(sexe), " : ", sprintf("%.4f", m),
      " | bande attendue [", bandes[[sexe]][1], " ; ", bandes[[sexe]][2], "]  ",
      if (m >= bandes[[sexe]][1] && m <= bandes[[sexe]][2]) "OK" else
        "*** HORS BANDE : diagnostic avant de conclure", "\n", sep = "")
}

# ---- ORACLE 3 : contrôle à taux nul, dispositif à DEUX étages
# (a) Ancrage INTER-SCRIPTS : à taux nul, VAP[A] = Σ_{k=1}^{37} ₖp₆₅ = e₆₅ − ½,
#     et le script 03 a stocké les e₆₅ PAR TRAJECTOIRE. C'est le seul contrôle
#     qui traverse une frontière de script : une réécriture de la convention en
#     C serait attrapée ici.
#     Dérogation assumée à la règle « entrée unique » : 03_ic_trajectoires.rds
#     est lu DANS CE BLOC UNIQUEMENT, aucune quantité tarifaire n'en dérive.
ic_B <- readRDS("resultats/03_ic_trajectoires.rds")

cat("\nORACLE 3a — taux nul : VAP[A] = e65 - 1/2, ancrage sur le script 03 :\n")
for (sexe in c("h", "f")) {
  vap0  <- vap_rente(S_base[[sexe]], 0, n_A)
  ecart <- max(abs(vap0 - (ic_B$e65$base[[sexe]] - 0.5)))
  cat("  ", nom_sexe(sexe), " : moyenne VAP0 = ", sprintf("%.4f", mean(vap0)),
      " | e65 reconstitue = ", sprintf("%.4f", mean(vap0) + 0.5),
      " (section B : ", sprintf("%.3f", mean(ic_B$e65$base[[sexe]])), ")\n", sep = "")
  cat("            max|VAP0 - (e65 - 1/2)| = ", format(ecart, digits = 3),
      "  ", if (ecart < 1e-12) "OK (tolerance 1e-12)" else
        "*** convention divergente entre les scripts 03 et 04", "\n", sep = "")
}
# LIMITE À ÉNONCER, PAS À CORRIGER : si les scripts 03 et 04 partageaient la
# même erreur d'alignement, cet ancrage croisé passerait quand même. La preuve
# d'alignement ABSOLU ne vient pas d'ici : elle repose sur le contrôle
# S["1", ] == exp(-µ["65", ]) ci-dessus, et sur les e₆₅ de la section B, déjà
# validés contre les taux observés d'A.1.

# (b) Chemin numérique INDÉPENDANT : la même quantité, calculée à la main.
# vap_rente_naive reconstruit la récurrence ₖp₆₅ = ₖ₋₁p₆₅ · exp(−µ₆₅₊ₖ₋₁) et
# somme terme à terme, sans aucune primitive vectorisée. Elle ne partage rien
# avec vap_rente() : un décalage d'indice dans la version matricielle ne peut
# pas se reproduire ici. (exp(-cumsum(µ)) aurait été le MÊME calcul par un
# autre chemin flottant — aussi creux que pas de test du tout.)
vap_rente_naive <- function(mu, taux, n) {
  p   <- 1        # ₀p₆₅ = 1 — jamais payé, la rente est à terme échu
  vap <- 0
  for (k in 1:n) {
    p   <- p * exp(-mu[as.character(64 + k)])   # k = 1 → µ₆₅ ; indexation par NOM
    vap <- vap + (1 / (1 + taux))^k * p         # paiement en fin d'année 2022+k
  }
  unname(vap)
}

cat("\nORACLE 3b — boucle scalaire naive vs version matricielle :\n")
for (sexe in c("h", "f")) {
  M <- diagonales$base[[sexe]]
  for (produit in c("A", "B")) {
    n <- if (produit == "A") n_A else duree_B
    for (taux in c(0, taux_tarif)) {
      naive <- apply(M, 2, vap_rente_naive, taux = taux, n = n)
      ecart <- max(abs(naive - vap_rente(S_base[[sexe]], taux, n)))
      cat("  ", nom_sexe(sexe), " [", produit, "] t = ",
          sprintf("%.0f %%", 100 * taux), " : max|naive - matriciel| = ",
          format(ecart, digits = 3), "  ",
          if (ecart < 1e-8) "OK" else "*** decalage d'indice", "\n", sep = "")
    }
  }
}
cat("  exemple lisible, trajectoire 1 (Hommes, [A], 3 %) : naive = ",
    sprintf("%.6f", vap_rente_naive(diagonales$base$h[, 1], taux_tarif, n_A)),
    " | matriciel = ", sprintf("%.6f", vap_base$A$h[1]), "\n", sep = "")

# ---- ORACLE 4 : monotonie en fonction du taux technique
# vᵏ = (1+t)^{-k} est strictement décroissant en t pour tout k ≥ 1, et les ₖp₆₅
# ne dépendent pas de t : la somme Σ vᵏ·ₖp₆₅ hérite de cette stricte décroissance.
cat("\nORACLE 4 — VAP moyenne strictement decroissante en t :\n")
for (produit in c("A", "B")) {
  for (sexe in c("h", "f")) {
    d <- diff(sensibilite[[paste0(produit, "_", sexe)]])
    cat("  [", produit, "] ", nom_sexe(sexe), " : decroissante sur 0-5 % : ",
        all(d < 0), "  ", if (all(d < 0)) "OK" else
          "*** logique d'actualisation cassee", "\n", sep = "")
  }
}

# ---- ORACLE 5 : effet DURATION, [A] plus sensible au taux que [B]
# La viagère paie jusqu'à 37 ans, la temporaire s'arrête à 15 : sa duration est
# bien plus courte, donc sa sensibilité au taux plus faible.
cat("\nORACLE 5 — variation relative de la VAP moyenne entre 1 % et 5 % :\n")
i1 <- which(sensibilite$taux == 0.01)
i5 <- which(sensibilite$taux == 0.05)
for (sexe in c("h", "f")) {
  var_A <- sensibilite[[paste0("A_", sexe)]][i5] / sensibilite[[paste0("A_", sexe)]][i1] - 1
  var_B <- sensibilite[[paste0("B_", sexe)]][i5] / sensibilite[[paste0("B_", sexe)]][i1] - 1
  cat("  ", nom_sexe(sexe), " : [A] ", sprintf("%+.2f %%", 100 * var_A),
      " | [B] ", sprintf("%+.2f %%", 100 * var_B), "  ",
      if (abs(var_A) > abs(var_B)) "OK — duration de [A] plus longue" else
        "*** effet duration inverse", "\n", sep = "")
}

# ---- ORACLE 6 : ordre H/F
# ATTENTION à la conception du test : les trajectoires H et F sont deux tirages
# SÉQUENTIELS du même flux (protocole RNG du script 03, ordre H puis F). La
# colonne i de H et la colonne i de F sont donc INDÉPENDANTES : leur
# appariement est un artefact de stockage, pas une structure du modèle. Un test
# par paire — all(vap_f > vap_h) — serait MAL CONÇU, et une violation isolée ne
# serait pas un bug. Le test porte sur des statistiques AGRÉGÉES.
# (Symétriquement, base et COVID sont eux VRAIMENT appariés : tous deux sont le
# premier simulate() après set.seed(1234) — d'où le diagnostic du bloc 5.)
cat("\nORACLE 6 — VAP(F) > VAP(H) sur la moyenne ET les quantiles :\n")
for (produit in c("A", "B")) {
  qh <- quantile(vap_base[[produit]]$h, c(0.025, 0.5, 0.975))
  qf <- quantile(vap_base[[produit]]$f, c(0.025, 0.5, 0.975))
  ok <- mean(vap_base[[produit]]$f) > mean(vap_base[[produit]]$h) && all(qf > qh)
  cat("  [", produit, "] moyennes H/F : ", sprintf("%.4f", mean(vap_base[[produit]]$h)),
      " / ", sprintf("%.4f", mean(vap_base[[produit]]$f)), "\n", sep = "")
  cat("        quantiles H : ", paste(sprintf("%.4f", qh), collapse = " / "),
      "\n        quantiles F : ", paste(sprintf("%.4f", qf), collapse = " / "),
      "  ", if (ok) "OK" else "*** ordre H/F viole", "\n", sep = "")
}

# ---- DERNIER TERME INCLUS (pas l'erreur de troncature) ---------
# v³⁷·₃₇p₆₅ est le DERNIER TERME RETENU par [A], pas l'erreur commise. La queue
# omise est Σ_{k≥38} vᵏ·ₖp₆₅, que la clôture de la table annule par convention
# (ℓ₁₀₃ = 0). Ce terme en fournit un MAJORANT D'ORDRE DE GRANDEUR : l'espérance
# de vie résiduelle à 102 ans valant ≈ 1,5-2 ans, la queue omise pèse
# grossièrement 1,5 à 2 fois ce terme. Ne jamais écrire « erreur de troncature
# = v³⁷·₃₇p₆₅ ».
cat("\nDernier terme inclus de [A] — majorant d'ordre de grandeur de la queue omise :\n")
for (sexe in c("h", "f")) {
  dernier <- mean((1 / (1 + taux_tarif))^n_A * S_base[[sexe]][as.character(n_A), ])
  cat("  ", nom_sexe(sexe), " : v^", n_A, " * ", n_A, "p65 = ",
      sprintf("%.5f", dernier), " soit ",
      sprintf("%.3f %%", 100 * dernier / mean(vap_base$A[[sexe]])),
      " de VAP[A]\n", sep = "")
}

# ==============================================================
# Bloc 7 — Sauvegardes (contrat d'interface vers D)
# ==============================================================
# Séparation STRUCTURELLE base / covid : la section D charge $base et rien
# d'autre. Aucun vecteur covid n'est le voisin d'un vecteur base — la forme de
# l'objet rend la confusion impossible, un commentaire seul n'y suffirait pas.
#
# $base$survies : les ₖp₆₅ par trajectoire, pour que le script 05 ne refasse pas
# la chaîne µ → p en recalculant BE₂₀₂₃ à taux_actu_be = 2 %.
# $base$primes  : primes pour une rente annuelle de 1 ; l'agrégation portefeuille
#                 (500 H + 500 F) est faite en SECTION D, où le portefeuille est
#                 défini.
# $base$sensibilite : MOYENNES seules. Aucune VAP par trajectoire à un autre
#                 taux que 3 % n'est sauvée — le BE de D s'actualise à 2 % et se
#                 recalcule depuis $base$survies, jamais depuis une VAP de C.
# $covid        : réservé au RAPPORT de la section C. La section D ne le charge
#                 JAMAIS (portée COVID bornée à C, décision actée).
saveRDS(list(
  base = list(vap = vap_base, survies = S_base, primes = primes,
              sensibilite = sensibilite),
  covid = list(vap = vap_covid, ecart_relatif = ecart_covid),
  taux_tarif = taux_tarif, n_sim = n_sim, duree_B = duree_B, n_A = n_A),
  file = "resultats/04_vap.rds")

cat("\n=== Bloc 7 — sauvegardes ======================================\n")
cat("resultats/04_vap.rds ecrit.\n")
cat("  $base$vap      : 4 vecteurs de", n_sim, "VAP (A/B x H/F), t = 3 %\n")
cat("  $base$survies  : 2 matrices ", paste(dim(S_base$h), collapse = " x "),
    " (kp65, rownames 1..", n_A, ")\n", sep = "")
cat("  $base$primes   : primes unitaires | $base$sensibilite :",
    nrow(sensibilite), "taux x 4 colonnes\n")
cat("  $covid$vap     : reserve au rapport de la section C\n")
cat("Figures : figures/fig_C_hist_vap.png, figures/fig_C_sensibilite_taux.png\n")
