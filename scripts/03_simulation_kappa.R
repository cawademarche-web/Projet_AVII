set.seed(1234)
# ==============================================================
# Script 03 — Section B : trajectoires projetées de κ_t, intervalles
#             de prédiction, comparaison avec le bootstrap d'A.3.vi
# Projet ACTU-F502 — Assurance Vie II
#
#   ln µ_x(t) = α_x + β_x·κ_t     avec κ_t = κ_{t-1} + d + ε_t, ε ~ N(0, σ²)
#
# B.1  bloc 3 — log-taux historiques et projetés de la cohorte (65 ans en 2022)
# B.2  blocs 1, 2, 4 — 5000 trajectoires de κ seul, IC de prédiction 95 %
# B.3  bloc 5 — comparaison avec les IC bootstrap d'A.3.vi
#      bloc 6 — sensibilité prudentielle COVID (alimente la section C)
#
# Références : docs/STMOMO_NOTES.md §1, docs/BOOTSTRAP_NOTES.md §3.
# ==============================================================
library(StMoMo)

n_sim   <- 5000   # N imposé par le brief (B.2)
horizon <- 50     # même horizon qu'en A.3.v/A.3.vi, pour que B.3 compare
                  # deux jeux d'IC définis sur les mêmes années

# ==============================================================
# Bloc 0 — Chargement (contrat d'interface)
# ==============================================================
# Le script 03 CHARGE les objets des scripts 01 et 02. Il ne refitte rien
# (sauf la variante 1970-2019, absente des RDS — voir bloc 6) et ne
# re-bootstrappe rien.
fit_lc         <- readRDS("resultats/02_fit_lc.rds")
LCfit_h        <- fit_lc$fit_h
LCfit_f        <- fit_lc$fit_f
ages_travail   <- fit_lc$ages_travail
annees_travail <- fit_lc$annees_travail

boot_ic     <- readRDS("resultats/02_boot_ic.rds")       # A.3.vi : ce que B.3 compare
forecast_lc <- readRDS("resultats/02_forecast_lc.rds")   # A.3.v  : projection centrale
donnees     <- readRDS("resultats/01_donnees.rds")       # pour le refit du bloc 6
taux_mle    <- readRDS("resultats/01_taux_mle.rds")      # A.1    : taux observés
variantes   <- readRDS("resultats/02_variantes.rds")     # grille des 36 fits (tableaux)

mu_brut <- list(h = taux_mle$mu_h, f = taux_mle$mu_f)
# Paramètres de la marche aléatoire estimés en A.3.v, repris tels quels : ce
# sont EXACTEMENT ceux que simulate() gèle (même appel interne à mrwd()).
d_chapeau     <- c(h = forecast_lc$h$kt.f$model$drift,
                   f = forecast_lc$f$kt.f$model$drift)
sigma_chapeau <- c(h = sqrt(forecast_lc$h$kt.f$model$sigma),
                   f = sqrt(forecast_lc$f$kt.f$model$sigma))
kappa_2023    <- c(h = LCfit_h$kt[1, "2023"], f = LCfit_f$kt[1, "2023"])

# ==============================================================
# Bloc 1 — LA distinction du projet : une source d'incertitude, pas trois
# ==============================================================
# Tableau repris À L'IDENTIQUE de docs/PIPELINE.md (§ A.3.vi) :
#
#                                    | α, β        | dérive & σ du RWD | aléa futur de κ
#   simulate(LCfit, nsim, h)   [B]   | fixés       | fixés             | simulé
#   simulate(LCboot, h)        [A.3] | ré-estimés  | ré-estimés        | simulé
#                                    | par réplic. | par réplication   |
#
# BDVK05 §4.2 : « spécification gelée, paramètres relâchés » — le bootstrap
# ne re-sélectionne pas de modèle ARIMA, mais il ré-estime ses paramètres sur
# chaque série κ*_t. La section A.3 capture donc TROIS sources d'incertitude,
# la section B UNE seule.
# La formulation « A.3 = incertitude sur (α, β) » est FAUSSE et ne doit
# apparaître nulle part.
#
# ---- PROTOCOLE RNG (contraignant) ----------------------------
# fit() consomme du RNG : gnm initialise ses paramètres non linéaires par des
# valeurs de départ ALÉATOIRES (test imprimé au bloc 6). Le dispositif de
# nombres aléatoires communs ne doit donc pas dépendre de la position du refit
# dans le script. Quatre règles :
#   1. set.seed(1234) immédiatement avant CHAQUE couple de simulate() — cas de
#      base inclus, pas seulement la variante COVID ;
#   2. même ordre H puis F dans les deux couples ;
#   3. aucun tirage entre les deux simulate() d'un couple (extraction de
#      diagonale, quantiles, saveRDS et rm n'en consomment pas) ;
#   4. sous-échantillon de la figure des trajectoires = les 200 PREMIÈRES
#      colonnes, jamais un sample() — les trajectoires sont i.i.d., donc les
#      200 premières valent un tirage, sans toucher au flux.

# simule_kappa : 5000 trajectoires de κ_t SEUL. α, β, la dérive d̂ et la
# volatilité σ̂ sont FIXÉS au fit central ; seul l'aléa futur ε_t est simulé.
# jumpchoice = "fit" : même amorce qu'en A.3.v et A.3.vi (VERROU GLOBAL).
# Si A.3 et B partaient d'amorces différentes, la comparaison B.3 mélangerait
# incertitude et simple décalage du point de départ.
# cellules : extrait les cellules (âge, année) d'un array âge × année ×
# trajectoire. Indexation par NOM (dimnames), jamais par position.
# Le retour NULL n'est pas un garde-fou : c'est le CAS de la variante
# 1970-2019, dont aucune année de la diagonale n'est dans $fitted.
cellules <- function(A, ages_d, annees_d) {
  garde <- as.character(annees_d) %in% dimnames(A)[[2]]
  if (!any(garde)) return(NULL)
  ij <- cbind(match(as.character(ages_d[garde]),   dimnames(A)[[1]]),
              match(as.character(annees_d[garde]), dimnames(A)[[2]]))
  apply(A, 3, function(M) M[ij])          # trajectoires en COLONNES
}

# diagonale_cohorte : lecture DIAGONALE µ_{65+k}(2022+k), k = 0…36 (âges 65 à
# 101, années 2022 à 2058) — la cohorte des assurés de 65 ans en 2022, pas un
# profil périodique. 37 valeurs = exactement ce que mobilise la convention de
# clôture d'A.1 (table jusqu'à ℓ_102, donc µ jusqu'à l'âge 101).
# Les années déjà calibrées viennent de $fitted, les suivantes de $rates :
# dans le cas de base 2022-2023 sont AJUSTÉES, dans la variante 1970-2019
# elles sont PROJETÉES. rbind(NULL, M) vaut M, donc aucun branchement.
diagonale_cohorte <- function(sim, age0 = 65, annee0 = 2022, n = 37) {
  ages_d   <- age0 + 0:(n - 1)
  annees_d <- annee0 + 0:(n - 1)
  M <- rbind(cellules(sim$fitted, ages_d, annees_d),
             cellules(sim$rates,  ages_d, annees_d))
  rownames(M) <- as.character(ages_d)
  M
}

# simule_kappa : produit les 5000 trajectoires et n'en garde que ce que la
# section B utilise. Sépare la production (~180 Mo en mémoire) de son résumé
# (quelques Mo) — l'array complet est libéré au retour de la fonction.
simule_kappa <- function(LCfit, fichier = NULL) {
  sim <- simulate(LCfit, nsim = n_sim, h = horizon,
                  kt.method = "mrwd", jumpchoice = "fit")
  if (!is.null(fichier)) saveRDS(sim, fichier)
  list(diagonale = diagonale_cohorte(sim),
       # IC de prédiction 95 % + médiane, cellule par cellule. Mêmes dimnames
       # qu'en A.3.vi, pour que B.3 confronte deux arrays de forme identique.
       quantiles = apply(sim$rates, c(1, 2), quantile, probs = c(0.025, 0.5, 0.975)),
       kappa     = sim$kt.s$sim[1, , ],                # années × trajectoires
       # Contrôle : les taux in-sample sont-ils bien ceux du fit central ?
       ecart_fitted = max(abs(sim$fitted[, , 1] - fitted(LCfit, type = "rates"))))
}

cat("=== Bloc 1 — 5000 trajectoires de κ, cas de base =============\n")
set.seed(1234)   # règle 1 du protocole RNG : re-seed avant le couple
temps <- system.time({
  base_h <- simule_kappa(LCfit_h, "resultats/LCsimK_h_5000.rds")
  base_f <- simule_kappa(LCfit_f, "resultats/LCsimK_f_5000.rds")
})
base <- list(h = base_h, f = base_f)
cat("simulate() x 2 sexes (nsim =", n_sim, ", h =", horizon, ") :",
    sprintf("%.1f", temps[["elapsed"]]), "s\n")
cat("dimensions : κ simulé", paste(dim(base_h$kappa), collapse = " x "),
    "| diagonale", paste(dim(base_h$diagonale), collapse = " x "),
    "| quantiles", paste(dim(base_h$quantiles), collapse = " x "), "\n")

# ==============================================================
# Bloc 2 — B.2 : contrôles de la mécanique de la marche aléatoire
# ==============================================================
# simulate.mrwd construit κ*_{2023+j} = κ_2023 + j·d̂ + σ̂·Σ_{i≤j} z*_i avec
# z* i.i.d. N(0,1) — ce sont de VRAIES marches aléatoires (algorithme 2 de
# HR09 §4.7), pas des droites déviées par un tirage unique. Donc :
#   E[κ*_{2023+j}] = κ_2023 + j·d̂     et     sd[κ*_{2023+j}] = σ̂·√j.
# Un écart notable signale une erreur de spécification du RWD.
cat("\n=== Bloc 2 — B.2 : contrôles du RWD ==========================\n")
for (sexe in c("h", "f")) {
  K <- base[[sexe]]$kappa
  cat("  ", if (sexe == "h") "Hommes" else "Femmes",
      ": κ_2023 =", sprintf("%.4f", kappa_2023[sexe]),
      "| d̂ =", sprintf("%.4f", d_chapeau[sexe]),
      "| σ̂ =", sprintf("%.4f", sigma_chapeau[sexe]), "\n")
  for (j in c(10, 25, 50)) {
    an <- as.character(2023 + j)
    cat("      h =", sprintf("%2d", j), "(", an, ") : moyenne =",
        sprintf("%9.4f", mean(K[an, ])),
        "[attendu", sprintf("%9.4f", kappa_2023[sexe] + j * d_chapeau[sexe]), "]",
        "| écart-type =", sprintf("%.4f", sd(K[an, ])),
        "[attendu σ̂·√h =", sprintf("%.4f", sigma_chapeau[sexe] * sqrt(j)), "]\n")
  }
  # Contrôle croisé direct avec la projection centrale d'A.3.v
  kt_central <- as.vector(forecast_lc[[sexe]]$kt.f$mean)
  cat("       moyenne des κ simulés en 2073 vs projection centrale A.3.v :",
      sprintf("%.4f", mean(K["2073", ])), "vs", sprintf("%.4f", tail(kt_central, 1)), "\n")
  # Les écarts ci-dessus sont-ils du bruit de Monte-Carlo ou une erreur de
  # spécification du RWD ? On standardise les innovations simulées :
  # z = (Δκ* − d̂)/σ̂ doit être N(0,1) sur les 5000 × 50 tirages. Un |t| de
  # quelques unités sur la MOYENNE est du bruit ; un sd(z) qui s'écarte de 1
  # serait une erreur de spécification. NB : les écarts aux trois horizons
  # sont emboîtés (mêmes trajectoires), ils ne constituent qu'UN test.
  inc <- rbind(K[1, ] - kappa_2023[sexe], diff(K))
  z   <- (inc - d_chapeau[sexe]) / sigma_chapeau[sexe]
  cat("       innovations standardisées (n =", length(z), ") : moyenne =",
      sprintf("%+.5f", mean(z)), "| t =", sprintf("%+.2f", mean(z) * sqrt(length(z))),
      "| sd =", sprintf("%.5f", sd(z)), "[attendu 1]\n")
}

# Le cône : la largeur de l'IC 95 % sur ln µ_65(t) doit CROÎTRE avec l'horizon.
cat("\nLargeur de l'IC 95 % sur ln µ_65(t) (le cône s'élargit-il ?) :\n")
for (sexe in c("h", "f")) {
  q <- base[[sexe]]$quantiles
  larg <- sapply(c("2024", "2048", "2073"),
                 function(an) log(q["97.5%", "65", an]) - log(q["2.5%", "65", an]))
  cat("  ", if (sexe == "h") "Hommes" else "Femmes", ":",
      paste(names(larg), sprintf("%.4f", larg), sep = " = ", collapse = " | "),
      "| croissante :", all(diff(larg) > 0), "\n")
}

# Contrôle destiné à la section D : sous jumpchoice = "fit", α, β et κ sont
# fixés jusqu'en 2023, donc les taux 2022-2023 sont IDENTIQUES sur les 5000
# trajectoires. La mortalité de l'année 1 ne porte aucun risque systématique.
cat("\nTaux in-sample : écart au fit central et dispersion sur les trajectoires\n")
for (sexe in c("h", "f")) {
  D <- base[[sexe]]$diagonale
  cat("  ", if (sexe == "h") "Hommes" else "Femmes",
      ": max|fitted - fit central| =", base[[sexe]]$ecart_fitted,
      "| étendue sur les", n_sim, "trajectoires : µ_65(2022) =",
      diff(range(D["65", ])), ", µ_66(2023) =", diff(range(D["66", ])), "\n")
}

# ==============================================================
# Bloc 3 — B.1 : log-taux historiques et projetés de la cohorte
# ==============================================================
# La diagonale de cohorte devient ici matériellement lisible : 37 points,
# contre 2 seulement sur données observées (A.3.iii).
cat("\n=== Bloc 3 — B.1 : diagonale de cohorte ======================\n")
diag_centrale <- forecast_lc$diagonale    # calculée en A.3.v, rechargée
for (sexe in c("h", "f")) {
  d <- diag_centrale[[sexe]]
  cat("  ", if (sexe == "h") "Hommes" else "Femmes",
      ": µ ajusté 65 ans (2022) =", sprintf("%.6f", d[1]),
      "| µ projeté 75 ans (2032) =", sprintf("%.6f", d[11]),
      "| 85 ans (2042) =", sprintf("%.6f", d[21]),
      "| 101 ans (2058) =", sprintf("%.6f", d[37]), "\n")
  cat("       µ observés (A.1) : 65 ans en 2022 =",
      sprintf("%.6f", mu_brut[[sexe]]["65", "2022"]),
      "| 66 ans en 2023 =", sprintf("%.6f", mu_brut[[sexe]]["66", "2023"]), "\n")
}

png("figures/fig_B1_taux_cohorte.png", width = 1500, height = 1000, res = 150)
plot(65:101, log(diag_centrale$h), type = "l", lwd = 2, col = "steelblue",
     ylim = range(log(unlist(diag_centrale))),
     xlab = "Âge x (année = 2022 + x - 65)",
     ylab = expression(ln~mu[x](2022 + x - 65)),
     main = "B.1 — log-taux de la cohorte des 65 ans en 2022")
lines(65:101, log(diag_centrale$f), lwd = 2, col = "firebrick")
points(65:66, log(c(mu_brut$h["65", "2022"], mu_brut$h["66", "2023"])),
       pch = 19, col = "steelblue")
points(65:66, log(c(mu_brut$f["65", "2022"], mu_brut$f["66", "2023"])),
       pch = 19, col = "firebrick")
abline(v = 66.5, lty = 3)
text(66.5, min(log(unlist(diag_centrale))), "ajusté | projeté", pos = 4, cex = 0.8)
legend("topleft", c("Hommes (Lee-Carter)", "Femmes (Lee-Carter)", "taux observés A.1"),
       col = c("steelblue", "firebrick", "black"), lwd = c(2, 2, NA),
       pch = c(NA, NA, 19), cex = 0.85)
dev.off()

# ==============================================================
# Bloc 4 — B.2 : intervalles de prédiction 95 %
# ==============================================================
couleurs_ages <- c("steelblue", "darkorange", "firebrick")
ages_sel      <- c("65", "75", "85")
annees_proj   <- 2024:2073

png("figures/fig_B2_trajectoires_kappa.png", width = 2000, height = 950, res = 150)
par(mfrow = c(1, 2))
for (sexe in c("h", "f")) {
  K <- base[[sexe]]$kappa
  # 200 PREMIÈRES trajectoires (cf. protocole RNG : pas de sample())
  matplot(annees_proj, K[, 1:200], type = "l", lty = 1, lwd = 1,
          col = adjustcolor("grey40", alpha.f = 0.15),
          xlab = "Année t", ylab = expression(kappa[t]),
          main = paste(if (sexe == "h") "Hommes" else "Femmes",
                       "— 200 trajectoires sur", n_sim))
  j <- seq_len(horizon)
  centre <- kappa_2023[sexe] + j * d_chapeau[sexe]
  lines(annees_proj, centre, lwd = 2, col = "firebrick")
  lines(annees_proj, centre + 1.96 * sigma_chapeau[sexe] * sqrt(j),
        lwd = 2, lty = 2, col = "firebrick")
  lines(annees_proj, centre - 1.96 * sigma_chapeau[sexe] * sqrt(j),
        lwd = 2, lty = 2, col = "firebrick")
  # Légende en ASCII (le device PNG sous Windows ne rend pas les indices
  # Unicode, cf. script 02) :
  legend("bottomleft", c("trajectoires simulees",
                         "centrale : kappa_2023 + h * d",
                         "IC 95 % : +/- 1.96 * sigma * sqrt(h)"),
         col = c("grey40", "firebrick", "firebrick"), lwd = c(1, 2, 2),
         lty = c(1, 1, 2), cex = 0.75)
}
par(mfrow = c(1, 1))
dev.off()

png("figures/fig_B2_ic_taux.png", width = 2000, height = 950, res = 150)
par(mfrow = c(1, 2))
for (sexe in c("h", "f")) {
  q <- base[[sexe]]$quantiles
  plot(annees_proj, log(q["50%", ages_sel[1], ]), type = "n",
       ylim = range(log(q[, ages_sel, ])), xlab = "Année t",
       ylab = expression(ln~mu[x](t)),
       main = paste(if (sexe == "h") "Hommes" else "Femmes",
                    "— IC 95 %, incertitude de κ seule"))
  for (k in seq_along(ages_sel)) {
    polygon(c(annees_proj, rev(annees_proj)),
            c(log(q["2.5%", ages_sel[k], ]), rev(log(q["97.5%", ages_sel[k], ]))),
            col = adjustcolor(couleurs_ages[k], alpha.f = 0.25), border = NA)
    lines(annees_proj, log(q["50%", ages_sel[k], ]), lwd = 2, col = couleurs_ages[k])
  }
  legend("bottomleft", paste("âge", ages_sel), col = couleurs_ages, lwd = 2, cex = 0.8)
}
par(mfrow = c(1, 1))
dev.off()

# ==============================================================
# Bloc 5 — B.3 : comparaison avec le bootstrap d'A.3.vi
# ==============================================================
# esperance_vie_cohorte : convention d'A.1, force de mortalité constante par
# morceaux. ₖp₆₅ = exp(−Σ_{j<k} µ_{65+j}), e₆₅ = Σ_{k≥1} ₖp₆₅ + ½.
esperance_vie_cohorte <- function(mu_diag) sum(cumprod(exp(-mu_diag))) + 0.5

# largeur_relative : largeur de l'IC 95 %, rapportée à la moyenne.
largeur_relative <- function(e) {
  q <- quantile(e, c(0.025, 0.975))
  as.vector((q[2] - q[1]) / mean(e))
}

cat("\n=== Bloc 5 — B.3 : B (κ seul) vs A.3.vi (bootstrap) ==========\n")
e65_B    <- lapply(base, function(b) apply(b$diagonale, 2, esperance_vie_cohorte))
e65_A3   <- boot_ic$e65_boot
larg_B   <- sapply(e65_B,  largeur_relative)
larg_A3  <- boot_ic$largeur_relative

cat("\nEspérance de vie de COHORTE à 65 ans (cohorte des 65 ans en 2022) :\n")
for (sexe in c("h", "f")) {
  qB  <- quantile(e65_B[[sexe]],  c(0.025, 0.975))
  qA3 <- quantile(e65_A3[[sexe]], c(0.025, 0.975))
  cat("  ", if (sexe == "h") "Hommes" else "Femmes", ":\n")
  cat("       B    (κ seul)   : moyenne =", sprintf("%.3f", mean(e65_B[[sexe]])),
      "| IC 95 % = [", sprintf("%.3f", qB[1]), ";", sprintf("%.3f", qB[2]), "]",
      "| largeur relative =", sprintf("%.2f %%", 100 * larg_B[sexe]), "\n")
  cat("       A.3.vi (boot.)  : moyenne =", sprintf("%.3f", mean(e65_A3[[sexe]])),
      "| IC 95 % = [", sprintf("%.3f", qA3[1]), ";", sprintf("%.3f", qA3[2]), "]",
      "| largeur relative =", sprintf("%.2f %%", 100 * larg_A3[sexe]), "\n")
}

# Comparaison au niveau des taux : rapport des largeurs d'IC sur ln µ en 2058
# (dernière année utile de la diagonale de cohorte).
cat("\nRapport des largeurs d'IC sur ln µ_x(2058), A.3 / B :\n")
for (sexe in c("h", "f")) {
  qB  <- base[[sexe]]$quantiles
  qA3 <- boot_ic$quantiles_taux[[sexe]]
  rap <- sapply(ages_sel, function(a) {
    lB  <- log(qB["97.5%", a, "2058"])  - log(qB["2.5%", a, "2058"])
    lA3 <- log(qA3["97.5%", a, "2058"]) - log(qA3["2.5%", a, "2058"])
    lA3 / lB
  })
  cat("  ", if (sexe == "h") "Hommes" else "Femmes", ":",
      paste("âge", names(rap), "=", sprintf("%.2f", rap), collapse = " | "), "\n")
}

# ---- ORACLE ---------------------------------------------------
# A.3 englobe STRICTEMENT les sources d'incertitude de B (trois contre une),
# donc ses intervalles DOIVENT être plus larges. Un élargissement perceptible
# mais de second ordre est le résultat classique : l'erreur de projection
# domine l'erreur d'estimation (Lee & Carter 1992 annexe B, BDVK05 §1,
# HR09 §2.9).
cat("\nORACLE B.3 — A.3 doit être STRICTEMENT plus large que B :\n")
for (sexe in c("h", "f")) {
  ratio <- larg_A3[sexe] / larg_B[sexe]
  verdict <- if (ratio < 1) {
    "*** BUG CERTAIN : les IC de B sont plus larges que ceux d'A.3 — stop et diagnostic"
  } else if (ratio >= 2) {
    "*** SUSPECT : facteur >= 2, chercher un bug avant de conclure"
  } else {
    "OK — élargissement de second ordre, l'erreur de projection domine l'erreur d'estimation"
  }
  cat("  ", if (sexe == "h") "Hommes" else "Femmes",
      ": largeur A.3 / largeur B =", sprintf("%.3f", ratio), "\n       ", verdict, "\n")
}

# ---- Pourquoi l'élargissement est-il de second ordre ? --------
# On décompose Var(κ*_{2023+h}) sous A.3 à partir des paramètres bootstrap
# DÉJÀ calculés par le script 02 (aucun re-bootstrap). Avec
#   κ*_{2023+h} = κ*_2023 + h·d̂* + σ̂*·Σ_{i≤h} z*_i :
#   Var ≈ Var(κ*_2023) + h²·Var(d̂*) + 2h·Cov(κ*_2023, d̂*) + E[σ̂*²]·h.
# Sous B, les trois premiers termes sont NULS et le dernier vaut σ̂²·h.
cat("\nDécomposition de Var(κ_2073) sous A.3 (h = 50) — d'où vient l'écart :\n")
LCboot <- readRDS("resultats/LCboot_5000.rds")
n_inc  <- length(annees_travail) - 1        # 53 accroissements Δκ
for (sexe in c("h", "f")) {
  kt_b <- sapply(LCboot[[sexe]]$bootParameters, function(p) as.vector(p$kt))
  d_b  <- apply(kt_b, 2, function(k) mean(diff(k)))
  s_b  <- apply(kt_b, 2, function(k) sd(diff(k)))
  kT_b <- kt_b[nrow(kt_b), ]
  v <- c(amorce      = var(kT_b),
         derive      = horizon^2 * var(d_b),
         covariance  = 2 * horizon * cov(kT_b, d_b),
         projection  = mean(s_b^2) * horizon)
  cat("  ", if (sexe == "h") "Hommes" else "Femmes", ":",
      paste(names(v), sprintf("%.3f", v), sep = " = ", collapse = " | "), "\n")
  cat("       sd(κ_2073) : A.3 =", sprintf("%.3f", sqrt(sum(v))),
      "| B =", sprintf("%.3f", sigma_chapeau[sexe] * sqrt(horizon)),
      "| rapport =", sprintf("%.3f", sqrt(sum(v)) / (sigma_chapeau[sexe] * sqrt(horizon))), "\n")
  # Le canal dominant n'est PAS la dérive : c'est σ̂* > σ̂. Le bruit
  # d'estimation sur κ̂ ajoute un MA(1) qui gonfle sd(Δκ*) — même mécanisme
  # que le ρ₁(Δκ) négatif documenté en A.3.i.
  cat("       σ̂* moyen =", sprintf("%.5f", mean(s_b)),
      "contre σ̂ central =", sprintf("%.5f", sigma_chapeau[sexe]),
      sprintf("(%+.1f %%)", 100 * (mean(s_b) / sigma_chapeau[sexe] - 1)), "\n")
  # Ce que NI A.3 NI B ne capturent : l'incertitude d'échantillonnage de la
  # dérive au sens de Lee & Carter (1992, annexe B), soit σ̂/√n. Le bootstrap
  # de Poisson conditionne sur le chemin κ̂ observé et ne perturbe que les
  # décès : il ne peut pas produire une autre réalisation de la marche.
  cat("       sd(d̂*) bootstrap =", sprintf("%.5f", sd(d_b)),
      "contre σ̂/√n =", sprintf("%.5f", sigma_chapeau[sexe] / sqrt(n_inc)),
      "— facteur", sprintf("%.0f", (sigma_chapeau[sexe] / sqrt(n_inc)) / sd(d_b)), "\n")
}

png("figures/fig_B3_comparaison_ic.png", width = 2100, height = 1400, res = 150)
par(mfrow = c(2, 3))
for (sexe in c("h", "f")) {
  qB  <- base[[sexe]]$quantiles
  qA3 <- boot_ic$quantiles_taux[[sexe]]
  for (a in ages_sel) {
    plot(annees_proj, log(qB["50%", a, ]), type = "n",
         ylim = range(log(qA3[, a, ]), log(qB[, a, ])),
         xlab = "Année t", ylab = expression(ln~mu[x](t)),
         main = paste(if (sexe == "h") "Hommes" else "Femmes", "— âge", a))
    polygon(c(annees_proj, rev(annees_proj)),
            c(log(qA3["2.5%", a, ]), rev(log(qA3["97.5%", a, ]))),
            col = adjustcolor("grey30", alpha.f = 0.30), border = NA)
    polygon(c(annees_proj, rev(annees_proj)),
            c(log(qB["2.5%", a, ]), rev(log(qB["97.5%", a, ]))),
            col = adjustcolor("steelblue", alpha.f = 0.40), border = NA)
    lines(annees_proj, log(qB["50%", a, ]), lwd = 2, col = "firebrick")
    legend("bottomleft", c("A.3.vi — bootstrap (3 sources)", "B — κ seul (1 source)"),
           fill = c(adjustcolor("grey30", alpha.f = 0.30),
                    adjustcolor("steelblue", alpha.f = 0.40)),
           border = NA, cex = 0.7)
  }
}
par(mfrow = c(1, 1))
dev.off()

png("figures/fig_B3_e65_comparaison.png", width = 2000, height = 950, res = 150)
par(mfrow = c(1, 2))
for (sexe in c("h", "f")) {
  # Titres en ASCII : les indices Unicode ne sont pas rendus par le device PNG
  # sous Windows (même remarque qu'au script 02).
  bornes <- range(e65_A3[[sexe]], e65_B[[sexe]])
  casiers <- seq(bornes[1], bornes[2], length.out = 40)
  hist(e65_A3[[sexe]], breaks = casiers, col = adjustcolor("grey30", alpha.f = 0.45),
       border = "white", xlab = expression(e[65]~"de cohorte"), ylab = "Fréquence",
       main = paste(if (sexe == "h") "Hommes" else "Femmes",
                    "— e_65 de cohorte, N =", n_sim))
  hist(e65_B[[sexe]], breaks = casiers, col = adjustcolor("steelblue", alpha.f = 0.55),
       border = "white", add = TRUE)
  abline(v = quantile(e65_A3[[sexe]], c(0.025, 0.975)), lwd = 2, lty = 2, col = "grey20")
  abline(v = quantile(e65_B[[sexe]],  c(0.025, 0.975)), lwd = 2, lty = 2, col = "steelblue")
  legend("topright", c("A.3.vi — bootstrap", "B — κ seul"),
         fill = c(adjustcolor("grey30", alpha.f = 0.45),
                  adjustcolor("steelblue", alpha.f = 0.55)), border = NA, cex = 0.8)
}
par(mfrow = c(1, 1))
dev.off()

# ==============================================================
# Bloc 6 — Sensibilité prudentielle COVID (second jeu, pour la section C)
# ==============================================================
# La variante 1970-2019 a une dérive plus forte (|β_65·d̂| supérieur de ~9 %) :
# c'est le sens PRUDENT pour un assureur de rentes. Elle alimente la
# sensibilité VAP de la section C. Pas de second bootstrap (la sensibilité
# porte sur la VAP centrale, pas sur les intervalles) ; section D non dupliquée.
#
# Le fit 60:102 × 1970-2019 n'est pas sauvegardé dans resultats/02_variantes.rds
# (ce fichier ne contient que les tableaux de synthèse), il est donc RÉ-AJUSTÉ
# ici. Ce qui garantit la reproduction n'est PAS une absence de tirage — gnm
# tire ses valeurs de départ — mais l'UNICITÉ DE L'OPTIMUM : le MLE de
# Lee-Carter est unique à contrainte d'identification près (Σβ = 1, Σκ = 0) et
# l'IRLS y converge quel que soit le point de départ. Le refit redonne donc les
# mêmes paramètres à la TOLÉRANCE DE CONVERGENCE près, pas au bit près — d'où
# le contrôle à 1e-6 ci-dessous.
# Pas de matrice de poids ici, contrairement au script 02 : la fenêtre
# 60:102 × 1970-2019 ne contient AUCUNE cellule ETR = 0 (les 3 cellules
# féminines recensées en A.3 sont à l'âge 102, années 1949-1951), donc
# w_{x,t} y vaut identiquement 1 et n'a rien à corriger.
ajuste_lc <- function(D, ETR, ages_fit, annees_fit) {
  fit(lc(link = "log"), Dxt = D, Ext = ETR,
      ages = donnees$ages, years = donnees$annees,
      ages.fit = ages_fit, years.fit = annees_fit, verbose = FALSE)
}

cat("\n=== Bloc 6 — sensibilité COVID : fit 60:102 x 1970-2019 ======\n")
graine_avant <- .Random.seed
LCfit19_h <- ajuste_lc(donnees$D_h, donnees$ETR_h, ages_travail, 1970:2019)
LCfit19_f <- ajuste_lc(donnees$D_f, donnees$ETR_f, ages_travail, 1970:2019)
cat("fit() laisse-t-il le flux RNG inchangé ? identical(.Random.seed) =",
    identical(graine_avant, .Random.seed),
    "\n  -> FALSE attendu : gnm tire ses valeurs de départ, d'où le re-seed",
    "avant CHAQUE couple de simulate().\n")

# Contrôle de reproduction contre la grille du script 02 (tolérance 1e-6)
cat("\nReproduction du fit 1970-2019 vs resultats/02_variantes.rds :\n")
grille19 <- variantes$grille[variantes$grille$ages == "60:102" &
                             variantes$grille$periode == "1970-2019", ]
for (sexe in c("h", "f")) {
  f  <- if (sexe == "h") LCfit19_h else LCfit19_f
  dk <- diff(as.vector(f$kt))
  ref <- grille19[grille19$sexe == sexe, ]
  obt <- c(drift = mean(dk), sigma = sd(dk), b65 = f$bx[f$ages == 65])
  att <- c(drift = ref$drift, sigma = ref$sigma, b65 = ref$b65)
  cat("  ", if (sexe == "h") "Hommes" else "Femmes", ":",
      paste(names(obt), sprintf("%.7f", obt), sep = " = ", collapse = " | "), "\n")
  cat("       écarts au tableau du script 02 :",
      paste(names(att), sprintf("%.1e", abs(obt - att)), sep = " ", collapse = " | "),
      "| tous < 1e-6 :", all(abs(obt - att) < 1e-6), "\n")
}

# NOMBRES ALÉATOIRES COMMUNS. En re-seedant, les deux scénarios consomment le
# même flux de normales par sexe : l'écart de VAP mesuré en section C est alors
# purement l'effet de la calibration (d̂ et σ̂ différents), sans bruit de
# Monte-Carlo (~0,3 % à N = 5000, soit l'ordre de l'effet cherché).
# NB : l'appariement se fait par INDICE DE PAS, pas par année civile — le cas
# de base part de κ_2023 (pas 1 = 2024), la variante de κ_2019 (pas 1 = 2020),
# soit un décalage de 4 indices. Les sommes cumulées partagent malgré tout ~33
# termes sur 36 : la corrélation reste massive et la réduction de variance
# fonctionne quasi pleinement.
# Corollaire à savoir énoncer : les deux jeux sont positivement corrélés par
# construction et ne pourront JAMAIS être fusionnés en un échantillon unique.
set.seed(1234)   # règle 1 du protocole RNG, même ordre H puis F
covid_h <- simule_kappa(LCfit19_h)   # array complet non sauvegardé :
covid_f <- simule_kappa(LCfit19_f)   # il ne sert à aucune figure de B
covid <- list(h = covid_h, f = covid_f)

e65_covid  <- lapply(covid, function(b) apply(b$diagonale, 2, esperance_vie_cohorte))
larg_covid <- sapply(e65_covid, largeur_relative)
cat("\ne₆₅ de cohorte sous la variante 1970-2019 (sert à la section C) :\n")
for (sexe in c("h", "f")) {
  cat("  ", if (sexe == "h") "Hommes" else "Femmes",
      ": moyenne =", sprintf("%.3f", mean(e65_covid[[sexe]])),
      "(cas de base :", sprintf("%.3f", mean(e65_B[[sexe]])), ")",
      "| largeur relative =", sprintf("%.2f %%", 100 * larg_covid[sexe]), "\n")
}

# ==============================================================
# Bloc 7 — Sauvegardes (contrat d'interface vers C et D)
# ==============================================================
# Deux niveaux : l'interface officielle (compacte, versionnée) et les arrays
# complets de simulate() (sur disque, hors versionnement).
diagonales <- list(
  base      = list(h = base_h$diagonale,  f = base_f$diagonale),
  covid2019 = list(h = covid_h$diagonale, f = covid_f$diagonale),
  ages = 65:101, annees = 2022:2058, n_sim = n_sim)
saveRDS(diagonales, file = "resultats/03_diagonales_cohorte.rds")

saveRDS(list(quantiles_taux   = list(h = base_h$quantiles, f = base_f$quantiles),
             kappa            = list(h = base_h$kappa, f = base_f$kappa),
             e65              = list(base = e65_B, covid2019 = e65_covid),
             largeur_relative = list(base = larg_B, covid2019 = larg_covid,
                                     bootstrap_A3 = larg_A3),
             n_sim = n_sim, horizon = horizon),
        file = "resultats/03_ic_trajectoires.rds")

cat("\n=== Bloc 7 — sauvegardes =====================================\n")
cat("Diagonales µ_{65+k}(2022+k) — interface officielle vers C et D :\n")
for (scen in c("base", "covid2019")) {
  for (sexe in c("h", "f")) {
    M <- diagonales[[scen]][[sexe]]
    cat("  ", scen, sexe, ":", paste(dim(M), collapse = " x "),
        "| NA :", sum(is.na(M)),
        "| µ_65(2022) médian =", sprintf("%.6f", median(M["65", ])),
        "| µ_101(2058) médian =", sprintf("%.6f", median(M["101", ])), "\n")
  }
}
fichiers <- c("resultats/03_diagonales_cohorte.rds", "resultats/03_ic_trajectoires.rds",
              "resultats/LCsimK_h_5000.rds", "resultats/LCsimK_f_5000.rds")
cat("\nFichiers produits :\n")
for (f in fichiers)
  cat("  ", f, ":", sprintf("%.2f", file.size(f) / 1024^2), "Mo\n")
