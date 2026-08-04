set.seed(1234)
# ==============================================================
# Script 02 — Lee-Carter Poisson : fit, diagnostic, projection,
#             bootstrap (section A.3 du brief, i à vi)
# Projet ACTU-F502 — Assurance Vie II
#
#   ln µ_x(t) = α_x + β_x·κ_t     sous  Σβ_x = 1  et  Σκ_t = 0
#   D_{x,t} ~ Poisson(ETR_{x,t} · µ_x(t))
#
# A.3.i   bloc 2 — plage d'âges et période de calibration
# A.3.ii  bloc 3 — paramètres estimés
# A.3.iii bloc 4 — log-taux MLE (A.1) vs log-taux Lee-Carter
# A.3.iv  bloc 5 — résidus
# A.3.v   bloc 6 — projection centrale
# A.3.vi  bloc 7 — bootstrap semi-paramétrique, IC 95 %, e_65 cohorte
#
# Références théoriques : docs/STMOMO_NOTES.md, docs/BOOTSTRAP_NOTES.md.
# ==============================================================
library(StMoMo)

n_boot <- 5000  # 200 pour valider le pipeline, puis 5000 pour la production.
                # ~2 h annoncées par la vignette à 5000 : coût mesuré au bloc 7.
horizon <- 50   # cf. justification bloc 6

# ==============================================================
# Bloc 0 — Chargement (contrat d'interface) et type d'exposition
# ==============================================================
# Le script 02 CHARGE les matrices du script 01, il ne relit jamais les
# fichiers HMD (pas de read.demogdata()).
donnees <- readRDS("resultats/01_donnees.rds")
D_h   <- donnees$D_h   ; D_f   <- donnees$D_f
ETR_h <- donnees$ETR_h ; ETR_f <- donnees$ETR_f
ages  <- donnees$ages  ; annees <- donnees$annees

# TYPE D'EXPOSITION. lc(link = "log") est le Lee-Carter POISSON du cours :
# la vraisemblance D_{x,t} ~ Poisson(ETR_{x,t}·µ_x(t)) suppose une exposition
# CENTRALE (ETR = années-personnes vécues), pas une exposition initiale.
# ATTENTION : quand on passe Dxt/Ext directement (au lieu d'un objet
# StMoMoData), StMoMo pose type = "central" d'office SANS RIEN VÉRIFIER.
# La vérification est donc à notre charge, et elle vient de A.1 : notre
# µ̂ = D/ETR y coïncide à ~2e-6 près avec le Mx publié par la HMD, qui est
# par définition le taux CENTRAL. Les Exposures_1x1 de la HMD sont bien des
# expositions centrales.

# ==============================================================
# Bloc 1 — Expositions nulles et matrice de poids w_{x,t}
# ==============================================================
# La log-vraisemblance de Poisson est indéfinie quand ETR = 0 (le terme
# D·ln µ n'a plus de sens, µ̂ = 0/0). On met ces cellules hors du fit par
# un poids nul plutôt que de rogner la plage d'âges.
matrice_poids <- function(ETR_fit) {
  w <- matrix(1, nrow = nrow(ETR_fit), ncol = ncol(ETR_fit),
              dimnames = dimnames(ETR_fit))
  w[ETR_fit <= 0] <- 0
  w
}

cat("=== Bloc 1 — expositions nulles ==============================\n")
cat("\nD et ETR aux âges 98-104 en 2022 (grands âges = faible exposition) :\n")
for (sexe in c("hommes", "femmes")) {
  D <- if (sexe == "hommes") D_h else D_f
  E <- if (sexe == "hommes") ETR_h else ETR_f
  cat(" ", sexe, ":\n")
  print(round(rbind(D   = D[as.character(98:104), "2022"],
                    ETR = E[as.character(98:104), "2022"],
                    mu  = D[as.character(98:104), "2022"] / E[as.character(98:104), "2022"]), 3))
}

# Comptage des cellules ETR = 0 par plage d'âges candidate (période complète)
cat("\nCellules ETR = 0 par plage candidate (années 1947-2023) :\n")
for (hi in c(100, 102, 105)) {
  cat("  ages 0-", hi, " : H =", sum(ETR_h[as.character(0:hi), ] <= 0),
      "| F =", sum(ETR_f[as.character(0:hi), ] <= 0), "\n")
}
# Détail des cellules concernées dans la plage retenue (plafond 102)
for (sexe in c("hommes", "femmes")) {
  E <- if (sexe == "hommes") ETR_h else ETR_f
  E <- E[as.character(0:102), ]
  idx <- which(E <= 0, arr.ind = TRUE)
  cat("  ", sexe, "— plafond 102 :", nrow(idx), "cellule(s) pondérée(s) à 0")
  if (nrow(idx) > 0) {
    cat(" | âges :", paste(unique(rownames(E)[idx[, 1]]), collapse = ", "),
        "| années :", paste(unique(colnames(E)[idx[, 2]]), collapse = ", "))
  }
  cat("\n")
}
# NB : fit() zéro-pondère DÉJÀ automatiquement les ETR <= 0 (avec un warning).
# On construit w_{x,t} explicitement pour pouvoir la compter, l'imprimer et la
# documenter — pas parce que StMoMo l'exigerait.

# ==============================================================
# Bloc 2 — A.3.i : plage d'âges et période de calibration
# ==============================================================
# ajuste_lc : un fit Lee-Carter Poisson sur une fenêtre (âges × années).
ajuste_lc <- function(D, ETR, ages_fit, annees_fit, bavard = FALSE) {
  ETR_fen <- ETR[as.character(ages_fit), as.character(annees_fit)]
  fit(lc(link = "log"), Dxt = D, Ext = ETR, ages = ages, years = annees,
      ages.fit = ages_fit, years.fit = annees_fit,
      wxt = matrice_poids(ETR_fen), verbose = bavard)
}

# diagnostic_kappa : les quatre chiffres qui décident si κ_t est projetable
# par une marche aléatoire avec dérive κ_t = κ_{t-1} + d + ε_t, ε ~ N(0, σ²).
#   R²   : linéarité de la tendance
#   d̂    : dérive estimée = moyenne des accroissements Δκ
#   σ̂    : volatilité des accroissements
#   ρ₁   : autocorrélation d'ordre 1 de Δκ.
# ATTENTION à l'interprétation de ρ₁ : on ne dispose pas de κ_t mais de son
# ESTIMATEUR κ̂_t = κ_t + η_t. Donc Δκ̂_t = Δκ_t + η_t − η_{t-1} : même si
# les vrais Δκ sont i.i.d., l'erreur d'estimation ajoute un MA(1) de
# coefficient négatif, dont l'autocorrélation d'ordre 1 est bornée par −0,5.
# Un ρ₁ négatif est donc la SIGNATURE ATTENDUE du bruit d'estimation, pas
# une violation de l'hypothèse de marche aléatoire. Corollaire utile :
# σ̂ = sd(Δκ̂) est gonflé par ce bruit, donc les intervalles de projection
# sont plutôt CONSERVATEURS — ce qui va dans le sens de la prudence pour un
# SCR (section D).
diagnostic_kappa <- function(kt, an) {
  kt <- as.vector(kt)
  dk <- diff(kt)
  c(R2    = summary(lm(kt ~ an))$r.squared,
    drift = mean(dk),
    sigma = sd(dk),
    rho1  = as.vector(acf(dk, lag.max = 1, plot = FALSE)$acf)[2])
}

# profil_kappa_sequentiel : procédure de Denuit & Goderniaux (2005), reprise
# par Haberman & Renshaw (2009) §3.10. On recule l'année de départ t0 et on
# regarde DEUX profils : R²(t0) localise le décrochage de linéarité, d̂(t0)
# dit si le choix de période a un effet MATÉRIEL sur les projections.
# C'est un diagnostic sur la série κ déjà estimée, pas un ré-ajustement.
profil_kappa_sequentiel <- function(kt, an, n_min = 20) {
  kt <- as.vector(kt)
  t0 <- an[1:(length(an) - n_min + 1)]
  prof <- sapply(t0, function(debut) {
    sel <- an >= debut
    c(diagnostic_kappa(kt[sel], an[sel])[c("R2", "drift")])
  })
  list(t0 = t0, R2 = prof["R2", ], drift = prof["drift", ])
}

# Grille de variantes : 3 plages × 3 débuts × 2 fins = 18 fits par sexe.
# La fin 2019 sert à MESURER l'effet des années COVID 2020-2021, pas à les
# exclure d'office.
plages   <- list("0:102" = 0:102, "50:102" = 50:102, "60:102" = 60:102)
debuts   <- c(1947, 1970, 1980)
fins     <- c(2023, 2019)

# Fenêtre commune de comparaison : déviance et BIC ne sont PAS comparables
# entre variantes (les jeux de données diffèrent), la RMSE sur une fenêtre
# identique pour toutes, si.
ages_comm   <- as.character(60:102)
annees_comm <- as.character(1980:2019)

taux_mle <- readRDS("resultats/01_taux_mle.rds")
mu_brut  <- list(h = taux_mle$mu_h, f = taux_mle$mu_f)
# Cellules où µ̂ = 0 (aucun décès observé) : ln µ̂ = -Inf, exclues des RMSE.
fini_comm <- lapply(mu_brut, function(m) is.finite(log(m[ages_comm, annees_comm])))
cat("\nCellules exclues des RMSE (µ̂ = 0) sur la fenêtre commune 60-102 x 1980-2019 :",
    "H =", sum(!fini_comm$h), "| F =", sum(!fini_comm$f), "\n")

cat("\n=== Bloc 2a — A.3.i : grille de 18 variantes par sexe ========\n")
grille <- NULL
fits_grille <- list()
for (sexe in c("h", "f")) {
  D <- if (sexe == "h") D_h else D_f
  E <- if (sexe == "h") ETR_h else ETR_f
  for (nom_plage in names(plages)) {
    for (debut in debuts) {
      for (fin in fins) {
        an <- debut:fin
        f  <- ajuste_lc(D, E, plages[[nom_plage]], an)
        cle <- paste(sexe, nom_plage, debut, fin, sep = "|")
        fits_grille[[cle]] <- f
        # RMSE sur les log-taux, fenêtre commune à toutes les variantes
        ecart <- log(fitted(f, type = "rates")[ages_comm, annees_comm]) -
                 log(mu_brut[[sexe]][ages_comm, annees_comm])
        # d̂ N'EST PAS COMPARABLE d'une plage d'âges à l'autre : sous la
        # contrainte d'identification Σβ_x = 1, l'échelle de κ dépend du
        # nombre d'âges (103 âges pour 0:102, 43 pour 60:102). La quantité
        # INVARIANTE est β_65·d̂ = taux annuel d'amélioration de ln µ_65,
        # et c'est elle qui pilote la VAP d'une rente à 65 ans.
        dg <- diagnostic_kappa(f$kt, an)
        grille <- rbind(grille, data.frame(
          sexe = sexe, ages = nom_plage, periode = paste0(debut, "-", fin),
          t(dg),
          b65 = f$bx[f$ages == 65],          # $bx est sans dimnames : on
          drift_65 = f$bx[f$ages == 65] * dg[["drift"]],   # sélectionne via $ages
          deviance = f$deviance,
          BIC      = -2 * f$loglik + f$npar * log(f$nobs),
          rmse_commune = sqrt(mean(ecart[fini_comm[[sexe]]]^2)),
          conv = f$conv, fail = f$fail, stringsAsFactors = FALSE))
      }
    }
  }
}
print(format(grille[, c("sexe", "ages", "periode", "R2", "drift", "b65",
                        "drift_65", "sigma", "rho1", "deviance", "BIC",
                        "rmse_commune")],
             digits = 4), row.names = FALSE)

# Contrôle de convergence : une variante non convergée est une ligne du
# tableau qui ne veut rien dire — à repérer AVANT de lire les R² et les d̂.
mauvais <- grille[!grille$conv | grille$fail != 0, ]
if (nrow(mauvais) == 0) {
  cat("\nConvergence :", nrow(grille), "fits, tous convergés (conv = TRUE, fail = 0).\n")
} else {
  cat("\n*** Variantes NON convergées :\n"); print(mauvais[, 1:3], row.names = FALSE)
}

cat("\n=== Bloc 2b — effet des années COVID sur la dérive d̂ =========\n")
# d̂ pilote la projection : si l'écart 2019 vs 2023 est important, le choix
# de la fin de période est matériel et doit être justifié dans le rapport.
covid <- NULL
for (sexe in c("h", "f")) {
  for (nom_plage in names(plages)) {
    for (debut in debuts) {
      sel <- function(fin) grille$drift[grille$sexe == sexe & grille$ages == nom_plage &
                                        grille$periode == paste0(debut, "-", fin)]
      d19 <- sel(2019); d23 <- sel(2023)
      covid <- rbind(covid, data.frame(sexe = sexe, ages = nom_plage, debut = debut,
                                       drift_2019 = d19, drift_2023 = d23,
                                       ecart_relatif = (d23 - d19) / abs(d19)))
    }
  }
}
print(format(covid, digits = 4), row.names = FALSE)

cat("\n=== Bloc 2c — profils séquentiels à rebours R²(t0) et d̂(t0) ==\n")
# Sur la série κ estimée sur la période COMPLÈTE 1947-2023.
profils <- list()
for (sexe in c("h", "f")) {
  for (nom_plage in names(plages)) {
    f <- fits_grille[[paste(sexe, nom_plage, 1947, 2023, sep = "|")]]
    p <- profil_kappa_sequentiel(f$kt, 1947:2023)
    # On travaille sur β_65·d̂, invariant d'échelle (cf. bloc 2a) : c'est le
    # taux annuel d'amélioration de ln µ_65 implicite dans la projection.
    p$drift_65 <- f$bx[f$ages == 65] * p$drift
    profils[[paste(sexe, nom_plage)]] <- p
  }
}
for (cle in names(profils)) {
  p <- profils[[cle]]
  # La queue du profil (t0 tardifs) ne conserve qu'une vingtaine d'années :
  # R² et la dérive y bougent pour une raison MÉCANIQUE (peu d'observations),
  # pas méthodologique. On donne donc aussi l'amplitude restreinte aux années
  # de départ où le diagnostic a du sens.
  st <- p$t0 <= 1990
  cat(" ", cle, ": R² max =", sprintf("%.4f", max(p$R2)), "en t0 =", p$t0[which.max(p$R2)],
      "| β_65·d̂ sur t0 = 1947-1990 : de", sprintf("%.5f", min(p$drift_65[st])), "à",
      sprintf("%.5f", max(p$drift_65[st])),
      "(amplitude", sprintf("%.1f %%", 100 * diff(range(p$drift_65[st])) / abs(mean(p$drift_65[st]))), ")",
      "| sur tout le profil :",
      sprintf("%.1f %%", 100 * diff(range(p$drift_65)) / abs(mean(p$drift_65))), "\n")
}
# HR09 note que κ est nettement plus linéaire sous un modèle APC que sous
# Lee-Carter, la courbure résiduelle du κ de LC étant ce que le terme cohorte
# absorbe — à relier à la heatmap des résidus du bloc 5. Aucun modèle APC
# n'est ajusté ici (le projet ne demande que Lee-Carter).

png("figures/fig_A3_profil_kappa_sequentiel.png", width = 2000, height = 1400, res = 150)
par(mfrow = c(2, 2))
couleurs_plages <- c("steelblue", "darkorange", "forestgreen")
for (sexe in c("h", "f")) {
  # Panneau de droite en β_65·d̂ et NON en d̂ : tracer d̂ superposerait trois
  # plages d'âges dont les échelles de κ diffèrent par construction
  # (Σβ_x = 1), ce qui inviterait au contre-sens.
  for (quoi in c("R2", "drift_65")) {
    ys <- lapply(names(plages), function(np) profils[[paste(sexe, np)]][[quoi]])
    plot(profils[[paste(sexe, names(plages)[1])]]$t0, ys[[1]], type = "n",
         ylim = range(unlist(ys)), xlab = expression(t[0]~"(année de départ)"),
         ylab = if (quoi == "R2") expression(R^2~"de"~kappa[t]~"~"~t)
                else expression(beta[65]%.%hat(d)),
         main = paste(if (sexe == "h") "Hommes" else "Femmes",
                      if (quoi == "R2") "— linéarité de κ"
                      else "— amélioration annuelle de ln µ(65)"))
    for (j in seq_along(ys)) lines(profils[[paste(sexe, names(plages)[j])]]$t0,
                                   ys[[j]], lwd = 2, col = couleurs_plages[j])
    legend("bottomleft", names(plages), col = couleurs_plages, lwd = 2, cex = 0.8)
  }
}
par(mfrow = c(1, 1))
dev.off()

# κ des trois plages d'âges sur la période complète (inspection visuelle de
# l'après-guerre, cf. A.3.i)
png("figures/fig_A3_variantes_kappa.png", width = 2000, height = 950, res = 150)
par(mfrow = c(1, 2))
for (sexe in c("h", "f")) {
  kts <- lapply(names(plages), function(np)
    as.vector(fits_grille[[paste(sexe, np, 1947, 2023, sep = "|")]]$kt))
  plot(1947:2023, kts[[1]], type = "n", ylim = range(unlist(kts)),
       xlab = "Année t", ylab = expression(kappa[t]),
       main = paste(if (sexe == "h") "Hommes" else "Femmes", "— κ, calibration 1947-2023"))
  for (j in seq_along(kts)) lines(1947:2023, kts[[j]], lwd = 2, col = couleurs_plages[j])
  legend("topright", names(plages), col = couleurs_plages, lwd = 2, cex = 0.8)
}
par(mfrow = c(1, 1))
dev.off()

# ---- Variante de travail : CAS DE BASE ------------------------
ages_travail   <- 60:102       # rentiers de 65 ans ; la mortalité jeune adulte
                               # relève d'un autre régime et un β_x commun aux
                               # deux serait une contrainte forte.
                               #  - âges 60-102 : plage adulte cohérente avec un produit de rente à 65 ans ;
                               #    plafond aligné sur age_max_table d'A.1
annees_travail <- 1970:2023    # exclut l'après-guerre immédiat ; le maximum du
                               # profil R²(t0) tombe en 1967-1970 pour les six
                               # séries. - début 1970 : maximum du profil R²(t0) en t0 = 1967-1970 sur les 6 séries,
                               # indépendamment du sexe et de la plage d'âges (Denuit-Goderniaux / HR09 §3.10)
# Fin de période TRANCHÉE sur 2023 (données complètes). Inclure 2020-2021
# réduit |d̂| d'environ 10,5 % (bloc 2b), donc réduit les VAP : c'est le sens
# ANTI-PRUDENTIEL pour un assureur de rentes, dont le risque est de
# sous-estimer l'amélioration future de la mortalité. La variante 1970-2019
# devient donc une SENSIBILITÉ PRUDENTIELLE, à propager en sections C et D
# avec son effet chiffré sur VAP et SCR (TODO section C).
cat("\nCas de base : âges", min(ages_travail), "-", max(ages_travail),
    "| années", min(annees_travail), "-", max(annees_travail),
    "| sensibilité prudentielle : 1970-2019\n")

# ==============================================================
# Bloc 3 — A.3.ii : paramètres estimés
# ==============================================================
cat("\n=== Bloc 3 — A.3.ii : paramètres =============================\n")
LCfit_h <- ajuste_lc(D_h, ETR_h, ages_travail, annees_travail, bavard = TRUE)
LCfit_f <- ajuste_lc(D_f, ETR_f, ages_travail, annees_travail, bavard = TRUE)

# Σβ_x = 1 et Σκ_t = 0 sont des contraintes d'IDENTIFICATION (le modèle
# α_x + β_x·κ_t est invariant par (α, β, κ) -> (α - βc, β/d, d(κ + c))),
# ré-imposées par centrage/rescaling à chaque itération de l'estimation.
# On les VÉRIFIE numériquement, on ne les suppose pas.
cat("\nContraintes d'identification :\n")
cat("  Σβ_x : H =", sum(LCfit_h$bx), "| F =", sum(LCfit_f$bx), "\n")
cat("  Σκ_t : H =", sum(LCfit_h$kt), "| F =", sum(LCfit_f$kt), "\n")

for (sexe in c("h", "f")) {
  f <- if (sexe == "h") LCfit_h else LCfit_f
  pente <- coef(lm(as.vector(f$kt) ~ annees_travail))[2]
  neg   <- f$ages[f$bx < 0]
  cat("  ", if (sexe == "h") "Hommes" else "Femmes",
      ": pente de κ_t =", sprintf("%.4f", pente),
      "| β_x < 0 :", length(neg), "âge(s)",
      if (length(neg)) paste0("(", paste(neg, collapse = ", "), ")") else "", "\n")
}

png("figures/fig_A3_parametres.png", width = 2100, height = 750, res = 150)
par(mfrow = c(1, 3))
plot(ages_travail, LCfit_h$ax, type = "l", lwd = 2, col = "steelblue",
     ylim = range(LCfit_h$ax, LCfit_f$ax), xlab = "Âge x", ylab = expression(alpha[x]),
     main = expression(alpha[x]~": niveau moyen de ln"~mu[x]))
lines(ages_travail, LCfit_f$ax, lwd = 2, col = "firebrick")
legend("topleft", c("Hommes", "Femmes"), col = c("steelblue", "firebrick"), lwd = 2)
plot(ages_travail, LCfit_h$bx, type = "l", lwd = 2, col = "steelblue",
     ylim = range(LCfit_h$bx, LCfit_f$bx), xlab = "Âge x", ylab = expression(beta[x]),
     main = expression(beta[x]~": sensibilité de l'âge x à"~kappa[t]))
abline(h = 0, lty = 3)
lines(ages_travail, LCfit_f$bx, lwd = 2, col = "firebrick")
plot(annees_travail, LCfit_h$kt, type = "l", lwd = 2, col = "steelblue",
     ylim = range(LCfit_h$kt, LCfit_f$kt), xlab = "Année t", ylab = expression(kappa[t]),
     main = expression(kappa[t]~": indice temporel de mortalité"))
lines(annees_travail, LCfit_f$kt, lwd = 2, col = "firebrick")
par(mfrow = c(1, 1))
dev.off()

# ==============================================================
# Bloc 4 — A.3.iii : log-taux MLE (A.1) vs log-taux Lee-Carter
# ==============================================================
# La diagonale de cohorte n'a que DEUX points observés (65 ans en 2022,
# 66 ans en 2023) : une comparaison cohorte sur données historiques est
# matériellement impossible. On compare donc en profil périodique (2022)
# et en série temporelle à l'âge 65.
cat("\n=== Bloc 4 — A.3.iii : MLE vs Lee-Carter =====================\n")
ages_c   <- as.character(ages_travail)
annees_c <- as.character(annees_travail)
ajuste <- list(h = fitted(LCfit_h, type = "rates"), f = fitted(LCfit_f, type = "rates"))

for (sexe in c("h", "f")) {
  brut <- log(mu_brut[[sexe]][ages_c, annees_c])
  lc   <- log(ajuste[[sexe]])
  fini <- is.finite(brut)
  rmse <- function(masque) sqrt(mean(((lc - brut)[masque & fini])^2))
  m_tout <- matrix(TRUE, nrow(brut), ncol(brut))
  m_2022 <- col(brut) == which(annees_c == "2022")
  cat("  ", if (sexe == "h") "Hommes" else "Femmes",
      ": RMSE(ln µ) fenêtre complète =", sprintf("%.4f", rmse(m_tout)),
      "| année 2022 =", sprintf("%.4f", rmse(m_2022)),
      "| cellules exclues (µ̂ = 0) :", sum(!fini), "\n")
}

png("figures/fig_A3_mle_vs_lc.png", width = 2000, height = 950, res = 150)
par(mfrow = c(1, 2))
for (sexe in c("h", "f")) {
  brut <- log(mu_brut[[sexe]][ages_c, "2022"])
  plot(ages_travail, brut, type = "p", pch = 20, col = "grey45",
       xlab = "Âge x", ylab = expression(ln~hat(mu)[x](2022)),
       main = paste(if (sexe == "h") "Hommes" else "Femmes", "— 2022"))
  lines(ages_travail, log(ajuste[[sexe]][, "2022"]), lwd = 2, col = "firebrick")
  legend("topleft", c("MLE brut (A.1)", "Lee-Carter ajusté"),
         col = c("grey45", "firebrick"), pch = c(20, NA), lty = c(NA, 1), lwd = c(NA, 2))
}
par(mfrow = c(1, 1))
dev.off()

png("figures/fig_A3_mle_vs_lc_age65.png", width = 2000, height = 950, res = 150)
par(mfrow = c(1, 2))
for (sexe in c("h", "f")) {
  brut <- log(mu_brut[[sexe]]["65", annees_c])
  plot(annees_travail, brut, type = "p", pch = 20, col = "grey45",
       xlab = "Année t", ylab = expression(ln~hat(mu)[65](t)),
       main = paste(if (sexe == "h") "Hommes" else "Femmes", "— âge 65"))
  lines(annees_travail, log(ajuste[[sexe]]["65", ]), lwd = 2, col = "firebrick")
  legend("topright", c("MLE brut (A.1)", "Lee-Carter ajusté"),
         col = c("grey45", "firebrick"), pch = c(20, NA), lty = c(NA, 1), lwd = c(NA, 2))
}
par(mfrow = c(1, 1))
dev.off()

# ==============================================================
# Bloc 5 — A.3.iv : résidus de déviance
# ==============================================================
# r_{x,t} = signe(D - D̂)·√(2[D·ln(D/D̂) - (D - D̂)]), D̂ = ETR·exp(α̂ + β̂κ̂).
# Des bandes DIAGONALES sur la heatmap = effet de cohorte non capturé par
# Lee-Carter : on le DOCUMENTE (TODO rapport), on ne change pas de modèle.
#
# Pour ne pas s'en tenir à une lecture à l'œil de la heatmap, on chiffre la
# structure : on regroupe les résidus par âge x, par année t, puis par
# COHORTE de naissance c = t − x, et on standardise chaque moyenne de groupe
# par √n. Sous l'hypothèse de résidus centrés indépendants, ces moyennes
# standardisées ont un écart-type de ~1. Un écart-type nettement > 1 pour le
# regroupement par cohorte, et plus élevé que pour les deux autres, signe un
# effet de cohorte que Lee-Carter (α_x + β_x·κ_t, sans terme en t − x) ne
# peut pas capturer.
# Lecture des deux autres regroupements : α_x et κ_t sont justement estimés
# pour annuler les moyennes par âge et par année, donc ces deux écarts-types
# sont attendus EN DESSOUS de 1 (marges sur-ajustées). Seule la cohorte
# n'est contrainte par aucun paramètre du modèle.
structure_residus <- function(r) {
  x <- r$ages[row(r$residuals)] ; t <- r$years[col(r$residuals)]
  grp <- list(age = x, annee = t, cohorte = t - x)
  sapply(grp, function(g) {
    m <- tapply(as.vector(r$residuals), g, function(z) mean(z) * sqrt(length(z)))
    n <- tapply(as.vector(r$residuals), g, length)
    sd(m[n >= 10])   # ~1 si pas de structure
  })
}

cat("\n=== Bloc 5 — A.3.iv : résidus ================================\n")
for (sexe in c("h", "f")) {
  f <- if (sexe == "h") LCfit_h else LCfit_f
  r <- residuals(f)
  cat("  ", if (sexe == "h") "Hommes" else "Femmes",
      ": résidus NA (cellules w = 0) =", sum(is.na(r$residuals)),
      "| moyenne =", sprintf("%.4f", mean(r$residuals, na.rm = TRUE)),
      "| écart-type =", sprintf("%.4f", sd(r$residuals, na.rm = TRUE)), "\n")
  st <- structure_residus(r)
  cat("       écart-type des moyennes standardisées par groupe",
      "(attendu ~1 si pas de structure) :",
      paste(names(st), sprintf("%.2f", st), sep = " = ", collapse = " | "), "\n")
  png(paste0("figures/fig_A3_residus_", toupper(sexe), ".png"),
      width = 1400, height = 1000, res = 150)
  plot(r, type = "colourmap",
       main = paste("Résidus de déviance —",
                    if (sexe == "h") "hommes" else "femmes"))
  dev.off()
}

# ==============================================================
# Bloc 6 — A.3.v : projection centrale
# ==============================================================
# Horizon : la cohorte de 65 ans en 2022 est suivie jusqu'à l'âge de clôture
# retenu en A.1 (age_max_table = 102, ce qui mobilise µ jusqu'à l'âge 101,
# atteint en 2058). Dernière année de calibration = 2023, donc h >= 35.
# h = 50 (jusqu'en 2073) couvre en outre les sections C et D.
# kt.method = "mrwd" : κ_t = κ_{t-1} + d + ε_t, la marche aléatoire avec
# dérive du cours (et non un ARIMA quelconque via "iarima").
# jumpchoice = "fit" : la projection part du taux AJUSTÉ de 2023, pas du
# taux observé. VERROUILLÉ et écrit explicitement ici, au bloc 7 et — à
# consigner — dans la section B : si A.3 et B partaient de points d'amorce
# différents, l'écart de largeur des IC mélangerait les sources
# d'incertitude et un simple décalage de départ.
cat("\n=== Bloc 6 — A.3.v : projection centrale =====================\n")
LCfor_h <- forecast(LCfit_h, h = horizon, kt.method = "mrwd", jumpchoice = "fit")
LCfor_f <- forecast(LCfit_f, h = horizon, kt.method = "mrwd", jumpchoice = "fit")

# Écart entre les deux amorces possibles, chiffré sans changer le verrou
cat("\nAmorce de la projection à l'âge 65 (TODO : arbitrage global fit/actual) :\n")
for (sexe in c("h", "f")) {
  cat("  ", if (sexe == "h") "Hommes" else "Femmes",
      ": µ ajusté 2023 =", sprintf("%.6f", ajuste[[sexe]]["65", "2023"]),
      "| µ observé 2023 =", sprintf("%.6f", mu_brut[[sexe]]["65", "2023"]),
      "| écart relatif =",
      sprintf("%.2f %%", 100 * (mu_brut[[sexe]]["65", "2023"] / ajuste[[sexe]]["65", "2023"] - 1)), "\n")
}
cat("\nDérive estimée d̂ de la marche aléatoire : H =",
    sprintf("%.4f", LCfor_h$kt.f$model$drift), "| F =",
    sprintf("%.4f", LCfor_f$kt.f$model$drift), "\n")

proj <- list(h = LCfor_h$rates, f = LCfor_f$rates)
png("figures/fig_A3_projection_centrale.png", width = 2000, height = 950, res = 150)
par(mfrow = c(1, 2))
for (sexe in c("h", "f")) {
  ages_sel <- c("65", "75", "85")
  y <- log(cbind(ajuste[[sexe]][ages_sel, ], proj[[sexe]][ages_sel, ]))
  matplot(c(annees_travail, as.integer(colnames(proj[[sexe]]))), t(y), type = "l",
          lty = 1, lwd = 2, col = c("steelblue", "darkorange", "firebrick"),
          xlab = "Année t", ylab = expression(ln~mu[x](t)),
          main = paste(if (sexe == "h") "Hommes" else "Femmes", "— ajusté puis projeté"))
  abline(v = 2023.5, lty = 3)
  legend("bottomleft", paste("âge", ages_sel),
         col = c("steelblue", "darkorange", "firebrick"), lwd = 2, cex = 0.8)
}
par(mfrow = c(1, 1))
dev.off()

# taux_diagonale : lecture DIAGONALE µ_{65+k}(2022+k) — la cohorte, pas un
# profil périodique. Indexation par NOM (dimnames), jamais par position.
taux_diagonale <- function(taux, age0 = 65, annee0 = 2022, n = 37) {
  taux[cbind(as.character(age0 + 0:(n - 1)), as.character(annee0 + 0:(n - 1)))]
}
# 37 valeurs = µ aux âges 65 à 101 : c'est exactement ce que mobilise la
# convention de clôture d'A.1 (table jusqu'à ℓ_102, donc µ jusqu'à l'âge 101).
n_diag <- 37
taux_complets <- list(h = cbind(ajuste$h, proj$h), f = cbind(ajuste$f, proj$f))
diag_centrale <- lapply(taux_complets, taux_diagonale, n = n_diag)

png("figures/fig_A3_projection_cohorte.png", width = 1400, height = 950, res = 150)
plot(65:101, log(diag_centrale$h), type = "l", lwd = 2, col = "steelblue",
     ylim = range(log(unlist(diag_centrale))), xlab = "Âge x (année = 2022 + x - 65)",
     ylab = expression(ln~mu[x](2022 + x - 65)),
     main = "Diagonale de cohorte — assurés de 65 ans en 2022")
lines(65:101, log(diag_centrale$f), lwd = 2, col = "firebrick")
legend("topleft", c("Hommes", "Femmes"), col = c("steelblue", "firebrick"), lwd = 2)
dev.off()

# ==============================================================
# Bloc 7 — A.3.vi : bootstrap semi-paramétrique
# ==============================================================
# BDVK05 §4.2, étape 1 : D*_{x,t} ~ Poisson(ETR_{x,t}·µ̂_x(t)), où µ̂_x(t)
# est l'estimateur MLE NON CONTRAINT, c'est-à-dire le taux brut D/ETR de la
# section A.1. Donc ETR_{x,t}·µ̂_x(t) = D_{x,t} : le tirage de Poisson est
# CENTRÉ SUR LES DÉCÈS OBSERVÉS. A.1 ne produit pas qu'un graphe — il
# fournit l'objet même autour duquel A.3.vi ré-échantillonne.
# Les ETR_{x,t} restent FIXES d'une réplication à l'autre (limite à énoncer).
#
# Terminologie : « semi-paramétrique » est le mot de StMoMo. BDVK05 dit
# « Poisson bootstrap » et oppose son approche au bootstrap PARAMÉTRIQUE de
# Brouhns-Denuit-Vermunt 2002b (tirage dans la loi normale multivariée
# asymptotique des estimateurs) — pas au bootstrap résiduel.
#
# TROIS sources d'incertitude ici, contre UNE seule en section B :
#   simulate(LCfit, nsim, h)  [section B] : α, β fixés | dérive et σ fixées | κ simulé
#   simulate(LCboot, h)       [section A.3] : α, β ré-estimés | dérive et σ
#                                             ré-estimées | κ simulé
# BDVK05 §4.2 : « spécification gelée, paramètres relâchés » — on ne
# re-sélectionne pas de modèle ARIMA, mais on ré-estime ses paramètres sur
# chaque série κ*_t.
cat("\n=== Bloc 7 — A.3.vi : bootstrap (nBoot =", n_boot, ") =============\n")

temps <- system.time({
  LCboot_h <- bootstrap(LCfit_h, nBoot = n_boot,
                        type = "semiparametric",   # = Poisson bootstrap (BDVK05)
                        deathType = "observed")    # BDVK05 §4.2 strict
  LCboot_f <- bootstrap(LCfit_f, nBoot = n_boot,
                        type = "semiparametric", deathType = "observed")
})
saveRDS(list(h = LCboot_h, f = LCboot_f),
        file = paste0("resultats/LCboot_", n_boot, ".rds"))   # IMMÉDIAT
cat("Bootstrap :", sprintf("%.1f", temps[["elapsed"]]), "s pour", n_boot,
    "réplications x 2 sexes =", sprintf("%.1f", temps[["elapsed"]] / 60), "min",
    "| coût projeté d'un run à 5000 :",
    sprintf("%.1f", temps[["elapsed"]] * 5000 / n_boot / 60), "min\n")

# Contrôle de convergence des réplications : à savoir AVANT de construire
# des quantiles dessus.
for (sexe in c("h", "f")) {
  bp <- (if (sexe == "h") LCboot_h else LCboot_f)$bootParameters
  na <- which(sapply(bp, function(p) anyNA(p$ax) || anyNA(p$bx) || anyNA(p$kt)))
  cat("  ", if (sexe == "h") "Hommes" else "Femmes", ": réplications à paramètres NA =",
      length(na), if (length(na)) paste0("(indices : ", paste(na, collapse = ", "), ")") else "", "\n")
}

LCsimPU_h <- simulate(LCboot_h, nsim = 1, h = horizon,
                      kt.method = "mrwd", jumpchoice = "fit")
LCsimPU_f <- simulate(LCboot_f, nsim = 1, h = horizon,
                      kt.method = "mrwd", jumpchoice = "fit")
saveRDS(list(h = LCsimPU_h, f = LCsimPU_f),
        file = paste0("resultats/LCsimPU_", n_boot, ".rds"))
sim <- list(h = LCsimPU_h, f = LCsimPU_f)

# ---- IC de prédiction 95 % sur les log-taux projetés ----------
quantiles_taux <- function(s) apply(s$rates, c(1, 2), quantile, probs = c(0.025, 0.5, 0.975))
q_taux <- lapply(sim, quantiles_taux)

png("figures/fig_A3_bootstrap_ic.png", width = 2000, height = 950, res = 150)
par(mfrow = c(1, 2))
annees_proj <- as.integer(dimnames(sim$h$rates)[[2]])
for (sexe in c("h", "f")) {
  q <- q_taux[[sexe]]
  ages_sel <- c("65", "75", "85")
  coul <- c("steelblue", "darkorange", "firebrick")
  plot(annees_proj, log(q[2, ages_sel[1], ]), type = "n",
       ylim = range(log(q[, ages_sel, ])), xlab = "Année t",
       ylab = expression(ln~mu[x](t)),
       main = paste(if (sexe == "h") "Hommes" else "Femmes", "— IC 95 % bootstrap"))
  for (j in seq_along(ages_sel)) {
    polygon(c(annees_proj, rev(annees_proj)),
            c(log(q[1, ages_sel[j], ]), rev(log(q[3, ages_sel[j], ]))),
            col = adjustcolor(coul[j], alpha.f = 0.25), border = NA)
    lines(annees_proj, log(q[2, ages_sel[j], ]), lwd = 2, col = coul[j])
  }
  legend("bottomleft", paste("âge", ages_sel), col = coul, lwd = 2, cex = 0.8)
}
par(mfrow = c(1, 1))
dev.off()

# ---- Espérance de vie de COHORTE à 65 ans ---------------------
# Convention d'A.1 : force de mortalité constante par morceaux,
#   ₖp₆₅ = exp(−Σ_{j<k} µ_{65+j}), e₆₅ = Σ_{k≥1} ₖp₆₅ + ½ (curtate + ½),
# table close à ℓ_102, donc µ mobilisé jusqu'à l'âge 101.
esperance_vie_cohorte <- function(mu_diag) sum(cumprod(exp(-mu_diag))) + 0.5

# Diagonale par réplication : années 2022-2023 depuis $fitted (les taux
# in-sample sont RÉ-ESTIMÉS à chaque réplication), 2024+ depuis $rates.
diagonale_simulee <- function(s) {
  idx <- function(A, ages_d, annees_d)
    cbind(match(as.character(ages_d), dimnames(A)[[1]]),
          match(as.character(annees_d), dimnames(A)[[2]]))
  i_h <- idx(s$fitted, 65:66, 2022:2023)
  i_p <- idx(s$rates,  67:101, 2024:2058)
  rbind(apply(s$fitted, 3, function(M) M[i_h]),   # trajectoires en COLONNES
        apply(s$rates,  3, function(M) M[i_p]))
}
e65_boot <- lapply(sim, function(s) apply(diagonale_simulee(s), 2, esperance_vie_cohorte))
e65_ponctuel <- sapply(diag_centrale, esperance_vie_cohorte)

# ---- Contrôle de convention de clôture ------------------------
# La référence 18,40 (H) / 21,57 (F) d'A.1 a été calculée sur les taux
# BRUTS. On lui repasse donc les µ̂_x(2023) bruts : c'est le TEST de clôture.
# Les taux ajustés Lee-Carter sont donnés à côté à titre d'INFORMATION
# (effet du lissage sur e₆₅ périodique) — ce n'est pas un test.
cat("\ne₆₅ périodique 2023, reconstruit par esperance_vie_cohorte() :\n")
for (sexe in c("h", "f")) {
  mu_brut_2023 <- mu_brut[[sexe]][as.character(65:101), "2023"]
  mu_lc_2023   <- ajuste[[sexe]][as.character(65:101), "2023"]
  cat("  ", if (sexe == "h") "Hommes" else "Femmes",
      ": µ̂ BRUTS (test de clôture, attendu", if (sexe == "h") "18,40" else "21,57", ") =",
      sprintf("%.3f", esperance_vie_cohorte(mu_brut_2023)),
      "| µ̂ ajustés LC (information A.3.iii) =",
      sprintf("%.3f", esperance_vie_cohorte(mu_lc_2023)), "\n")
}

# ---- Quatre diagnostics post-bootstrap ------------------------
cat("\nDiagnostics post-bootstrap (étalons : BDVK05 §4.3, vignette StMoMo §8) :\n")
cat("\ne₆₅ de COHORTE (cohorte des 65 ans en 2022, lecture diagonale) :\n")
larg_rel <- c()
for (sexe in c("h", "f")) {
  e <- e65_boot[[sexe]]
  q <- quantile(e, c(0.025, 0.975))
  larg_rel[sexe] <- (q[2] - q[1]) / mean(e)
  cat("  ", if (sexe == "h") "Hommes" else "Femmes",
      ": ponctuel =", sprintf("%.3f", e65_ponctuel[[sexe]]),
      "| moyenne bootstrap =", sprintf("%.3f", mean(e)),
      "| IC 95 % = [", sprintf("%.3f", q[1]), ";", sprintf("%.3f", q[2]), "]\n")
  # (i) biais : étalon BDVK05 = 0,1 % à 1,0 %
  cat("       (i)  biais relatif moyenne/ponctuel :",
      sprintf("%+.2f %%", 100 * (mean(e) - e65_ponctuel[[sexe]]) / e65_ponctuel[[sexe]]),
      "  [étalon BDVK05 : 0,1 à 1,0 %]\n")
  # (ii) largeur relative : étalon quelques % à ~15 %. ATTENTION : les
  # intervalles de BDVK05 sont à 90 % et portent sur a₆₅ ; les nôtres sont
  # à 95 % sur e₆₅, donc mécaniquement plus larges.
  cat("       (ii) largeur relative de l'IC 95 % :", sprintf("%.2f %%", 100 * larg_rel[sexe]),
      "  [étalon BDVK05 à 90 % : 3,9 % F à 15,7 % H]\n")
}
# (iii) asymétrie H/F : BDVK05 obtient systématiquement les hommes plus
# larges, et l'attribue à des IC plus larges sur la projection des κ
# masculins. On imprime donc les σ̂(Δκ) des deux sexes à côté du ratio : ce
# sont eux qui expliquent l'ampleur (ou la faiblesse) de l'écart.
sigma_kappa <- c(h = diagnostic_kappa(LCfit_h$kt, annees_travail)[["sigma"]],
                 f = diagnostic_kappa(LCfit_f$kt, annees_travail)[["sigma"]])
cat("       (iii) rapport des largeurs H/F :", sprintf("%.2f", larg_rel["h"] / larg_rel["f"]),
    "| σ̂(Δκ) : H =", sprintf("%.4f", sigma_kappa["h"]),
    "F =", sprintf("%.4f", sigma_kappa["f"]),
    if (larg_rel["h"] > larg_rel["f"]) "— hommes plus larges, conforme à BDVK05\n"
    else "— femmes plus larges, INVERSE de BDVK05\n")

# (iv) hiérarchie des dispersions paramétriques. α, β et κ n'ont pas la même
# échelle : on donne l'écart-type bootstrap brut ET normalisé par l'amplitude
# (max − min) du paramètre au fit central.
cat("\n       (iv) dispersion bootstrap des paramètres",
    "[attendu : β̂ nettement plus dispersé] :\n")
params_boot <- function(LCboot, quoi)
  sapply(LCboot$bootParameters, function(p) as.vector(p[[quoi]]))
for (sexe in c("h", "f")) {
  b <- if (sexe == "h") LCboot_h else LCboot_f
  f <- if (sexe == "h") LCfit_h else LCfit_f
  cat("        ", if (sexe == "h") "Hommes" else "Femmes", ":\n")
  for (quoi in c("ax", "bx", "kt")) {
    et  <- apply(params_boot(b, quoi), 1, sd)
    amp <- diff(range(as.vector(f[[quoi]])))
    cat("          ", quoi, ": écart-type moyen =", sprintf("%.5f", mean(et)),
        "| normalisé par l'amplitude =", sprintf("%.4f", mean(et) / amp), "\n")
  }
  et_b <- apply(params_boot(b, "bx"), 1, sd)
  names(et_b) <- ages_travail
  cat("           écart-type de β̂ aux âges 65 / 85 / 101 / 102 :",
      paste(sprintf("%.5f", et_b[as.character(c(65, 85, 101, 102))]), collapse = " / "), "\n")
  # Le fit central donne β_x > 0 partout, mais l'ENVELOPPE bootstrap peut
  # traverser 0 aux âges extrêmes : une fraction des réplications y projette
  # alors une mortalité CROISSANTE (β_x < 0 avec κ décroissant).
  q_b <- apply(params_boot(b, "bx"), 1, quantile, probs = 0.025)
  names(q_b) <- ages_travail
  zero <- ages_travail[q_b < 0]
  cat("           âges où l'enveloppe bootstrap de β̂ traverse 0 :",
      if (length(zero)) paste(range(zero), collapse = "-") else "aucun",
      "|", length(zero), "âge(s)\n")
}

png("figures/fig_A3_hist_e65_cohorte.png", width = 2000, height = 950, res = 150)
par(mfrow = c(1, 2))
for (sexe in c("h", "f")) {
  # Titre en ASCII : les indices Unicode (₆₅) ne sont pas rendus par le
  # device PNG sous Windows. La notation e_65 passe par expression() en xlab.
  hist(e65_boot[[sexe]], breaks = 30, col = "grey85", border = "white",
       xlab = expression(e[65]~"de cohorte"), ylab = "Fréquence",
       main = paste(if (sexe == "h") "Hommes" else "Femmes",
                    "— esperance de vie de cohorte a 65 ans, N =", n_boot))
  abline(v = e65_ponctuel[[sexe]], lwd = 2, col = "firebrick")
  abline(v = quantile(e65_boot[[sexe]], c(0.025, 0.975)), lwd = 2, lty = 2, col = "steelblue")
  legend("topright", c("prévision ponctuelle", "IC 95 %"),
         col = c("firebrick", "steelblue"), lty = c(1, 2), lwd = 2, cex = 0.8)
}
par(mfrow = c(1, 1))
dev.off()

png("figures/fig_A3_bootstrap_parametres.png", width = 2100, height = 1400, res = 150)
par(mfrow = c(2, 3))
for (sexe in c("h", "f")) {
  b <- if (sexe == "h") LCboot_h else LCboot_f
  f <- if (sexe == "h") LCfit_h else LCfit_f
  for (quoi in c("ax", "bx", "kt")) {
    P <- params_boot(b, quoi)
    q <- apply(P, 1, quantile, probs = c(0.025, 0.975))
    abs_x <- if (quoi == "kt") annees_travail else ages_travail
    plot(abs_x, as.vector(f[[quoi]]), type = "n", ylim = range(q),
         xlab = if (quoi == "kt") "Année t" else "Âge x", ylab = quoi,
         main = paste(if (sexe == "h") "Hommes" else "Femmes", "—", quoi))
    polygon(c(abs_x, rev(abs_x)), c(q[1, ], rev(q[2, ])),
            col = "grey85", border = NA)
    lines(abs_x, as.vector(f[[quoi]]), lwd = 2, col = "firebrick")
  }
}
par(mfrow = c(1, 1))
dev.off()

# ==============================================================
# Bloc 8 — Sauvegardes (contrat d'interface pour B, C et D)
# ==============================================================
saveRDS(list(grille = grille, covid = covid, profils = profils),
        file = "resultats/02_variantes.rds")
saveRDS(list(fit_h = LCfit_h, fit_f = LCfit_f,
             ages_travail = ages_travail, annees_travail = annees_travail),
        file = "resultats/02_fit_lc.rds")
saveRDS(list(h = LCfor_h, f = LCfor_f, diagonale = diag_centrale,
             e65_ponctuel = e65_ponctuel),
        file = "resultats/02_forecast_lc.rds")
saveRDS(list(quantiles_taux = q_taux, e65_boot = e65_boot,
             n_boot = n_boot, largeur_relative = larg_rel),
        file = "resultats/02_boot_ic.rds")
cat("\nSauvés : resultats/02_variantes.rds, 02_fit_lc.rds, 02_forecast_lc.rds,",
    "02_boot_ic.rds, LCboot_", n_boot, ".rds, LCsimPU_", n_boot, ".rds\n", sep = "")
