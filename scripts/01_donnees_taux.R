# ==============================================================
# Script 01 — Données HMD, taux bruts MLE et indicateurs (A.1, A.2)
# Projet ACTU-F502 — Assurance Vie II
#
# Contenu :
#   1. Lecture de Deaths_1x1.txt et Exposures_1x1.txt (HMD Autriche,
#      1947-2023) -> matrices D_{x,t} et ETR_{x,t}, hommes / femmes.
#   2. Section A.1 : taux bruts µ̂_x(t) = D_{x,t}/ETR_{x,t} (MLE du
#      modèle de Poisson) avec IC à 95 %, figures en échelle log.
#   3. Section A.2 : tables de survie périodiques -> e_0, e_65,
#      médiane et quartiles de l'âge au décès (matière pour la
#      discussion expansion / rectangularisation du rapport).
# ==============================================================

# ---- Lecture des fichiers HMD --------------------------------
# Format HMD : 2 lignes d'en-tête puis colonnes Year / Age / Female / Male / Total.
# na.strings = "." : convention HMD pour les valeurs manquantes.
deces <- read.table("donnees/Deaths_1x1.txt",    skip = 2, header = TRUE, na.strings = ".")
expos <- read.table("donnees/Exposures_1x1.txt", skip = 2, header = TRUE, na.strings = ".")

# L'âge "110+" (groupe ouvert) est recodé en 110 pour avoir un âge entier
deces$Age <- as.integer(sub("\\+", "", deces$Age))
expos$Age <- as.integer(sub("\\+", "", expos$Age))

ages   <- sort(unique(deces$Age))   # 0, 1, ..., 110
annees <- sort(unique(deces$Year))  # 1947, ..., 2023

# ---- Mise en forme en matrices âges x années -----------------
# hmd_en_matrice : passe du format long HMD (une ligne par couple (t, x))
# à une matrice avec les âges x en lignes et les années t en colonnes.
hmd_en_matrice <- function(df, sexe) {
  df <- df[order(df$Year, df$Age), ]  # garantit le remplissage colonne par colonne
  matrix(df[[sexe]], nrow = length(ages), ncol = length(annees),
         dimnames = list(ages, annees))
}

# Matrices de décès D_{x,t} et d'exposition ETR_{x,t}, hommes (h) / femmes (f)
D_h   <- hmd_en_matrice(deces, "Male")
D_f   <- hmd_en_matrice(deces, "Female")
ETR_h <- hmd_en_matrice(expos, "Male")
ETR_f <- hmd_en_matrice(expos, "Female")

# ---- Contrôles -----------------------------------------------
cat("Dimensions (âges x années) :\n")
cat("  D_h   :", dim(D_h),   "\n")
cat("  D_f   :", dim(D_f),   "\n")
cat("  ETR_h :", dim(ETR_h), "\n")
cat("  ETR_f :", dim(ETR_f), "\n")
# matrix() recycle silencieusement si le panneau (t, x) est incomplet :
# on vérifie que chaque fichier a bien une ligne par couple (t, x)
cat("Lignes lues : deces =", nrow(deces), "| expos =", nrow(expos),
    "| attendu (ages x annees) :", length(ages) * length(annees), "\n")

# La HMD répartit les décès d'âge inconnu entre les âges → décès non entiers
cat("Décès non entiers : D_h =", sum(D_h %% 1 != 0),
    "| D_f =", sum(D_f %% 1 != 0), "\n")

# Rappel : "110+" est un groupe ouvert recodé en 110 ; µ̂_110 = D_110+/ETR_110+
# est donc un taux agrégé sur les âges >= 110, pas un taux à âge unique.

cat("\nExtrait D_h (décès hommes), âges 60-70 x années 2018-2023 :\n")
print(D_h[as.character(60:70), as.character(2018:2023)])

cat("\nExtrait ETR_h (exposition hommes), âges 60-70 x années 2018-2023 :\n")
print(round(ETR_h[as.character(60:70), as.character(2018:2023)], 1))

cat("\nValeurs manquantes : D =", sum(is.na(D_h)) + sum(is.na(D_f)),
    "| ETR =", sum(is.na(ETR_h)) + sum(is.na(ETR_f)), "\n")

# ---- Sauvegarde ----------------------------------------------
# Contrat d'interface : le script 02 CHARGE ces objets, il ne recalcule rien.
saveRDS(list(D_h = D_h, D_f = D_f, ETR_h = ETR_h, ETR_f = ETR_f,
             ages = ages, annees = annees),
        file = "resultats/01_donnees.rds")
cat("\nSauvé : resultats/01_donnees.rds\n")

# ==============================================================
# Section A.1 — Taux bruts par maximum de vraisemblance
# ==============================================================
# Modèle du cours : D_{x,t} ~ Poisson(µ_x(t) · ETR_{x,t})
#   =>  µ̂_x(t) = D_{x,t} / ETR_{x,t}   (estimateur du maximum de vraisemblance)
# Aux âges 103+ certaines cellules ont D = ETR = 0 (plus personne en vie) :
# 0/0 donne NaN, laissé tel quel. Aucun cas D > 0 avec ETR = 0, donc jamais Inf.
mu_h <- D_h / ETR_h
mu_f <- D_f / ETR_f

# Var(µ̂_x(t)) = µ̂/ETR (information de Fisher du modèle de Poisson)
#   =>  IC à 95 % : µ̂ ± 1.96·√(µ̂/ETR)
ic_bas_h  <- mu_h - 1.96 * sqrt(mu_h / ETR_h)
ic_haut_h <- mu_h + 1.96 * sqrt(mu_h / ETR_h)
ic_bas_f  <- mu_f - 1.96 * sqrt(mu_f / ETR_f)
ic_haut_f <- mu_f + 1.96 * sqrt(mu_f / ETR_f)

# ---- Contrôles A.1 -------------------------------------------
# Contrôle croisé : la HMD publie elle-même Mx = D/ETR (arrondi à 6 décimales) ;
# l'écart max avec nos µ̂ doit être de l'ordre de l'arrondi.
mx_hmd <- read.table("donnees/Mx_1x1.txt", skip = 2, header = TRUE, na.strings = ".")
mx_hmd$Age <- as.integer(sub("\\+", "", mx_hmd$Age))
Mx_h <- hmd_en_matrice(mx_hmd, "Male")
Mx_f <- hmd_en_matrice(mx_hmd, "Female")
# C'est LE test de la mise en forme (appariement Male/Female, recodage 110+,
# ordre de remplissage, alignement âge/année) : une inversion ou un décalage
# donnerait un écart d'ordre 1. On imprime deux seuils :
#  - toutes cellules finies : l'écart y est dominé par les âges extrêmes, où le
#    fichier HMD publie ETR arrondi à 2 décimales (ex. D = 1, ETR = 0.17 ->
#    µ̂ = 5.88 contre Mx = 6). Ce n'est pas une erreur de calcul.
#  - cellules ETR >= 1000 : il ne reste que l'arrondi de Mx à 6 décimales (~1e-6),
#    c'est le seuil sur lequel porte le 🚩 de la checklist (> 1e-4).
ec_h <- abs(mu_h - Mx_h)
ec_f <- abs(mu_f - Mx_f)
cat("\nContrôle croisé µ̂ vs Mx publié (HMD) — écart max :\n")
cat("  toutes cellules finies : H =", max(ec_h[is.finite(ec_h)]),
    "| F =", max(ec_f[is.finite(ec_f)]), "\n")
cat("  cellules ETR >= 1000   : H =", max(ec_h[ETR_h >= 1000]),
    "| F =", max(ec_f[ETR_f >= 1000]), "\n")

# Contrôle checklist : µ̂ croissant avec l'âge
sel_ages <- as.character(c(40, 65, 90))
cat("\nµ̂(2022) aux âges 40 / 65 / 90 (attendu : croissant) :\n")
print(round(rbind(hommes = mu_h[sel_ages, "2022"],
                  femmes = mu_f[sel_ages, "2022"]), 5))

# Contrôle checklist : précision des IC en 2022.
# - largeur absolue = 2·1.96·√(µ̂/ETR) : croît mécaniquement avec l'âge
#   (µ̂(90) >> µ̂(50)) -> imprimée pour information, PAS un test.
# - largeur relative = largeur/µ̂ = 2·1.96/√D : ne dépend QUE du nombre de décès
#   (l'exposition se simplifie) -> profil attendu en U, minimal au mode des
#   décès (~85-90 ans), large aux âges jeunes et aux âges extrêmes.
sel_ic <- as.character(c(25, 50, 70, 88, 100))
larg_abs_h <- (ic_haut_h - ic_bas_h)[sel_ic, "2022"]
larg_abs_f <- (ic_haut_f - ic_bas_f)[sel_ic, "2022"]
larg_rel_h <- larg_abs_h / mu_h[sel_ic, "2022"]
larg_rel_f <- larg_abs_f / mu_f[sel_ic, "2022"]
cat("\nLargeurs d'IC 95 % en 2022 (profil en U de la largeur relative) :\n")
print(round(rbind(absolue_h = larg_abs_h, relative_h = larg_rel_h,
                  D_h = D_h[sel_ic, "2022"],
                  absolue_f = larg_abs_f, relative_f = larg_rel_f,
                  D_f = D_f[sel_ic, "2022"]), 5))

# 🚩 checklist : la largeur relative doit valoir exactement 2·1.96/√D
cat("Largeur relative == 2·1.96/√D : H =",
    isTRUE(all.equal(unname(larg_rel_h), unname(2 * 1.96 / sqrt(D_h[sel_ic, "2022"])))),
    "| F =",
    isTRUE(all.equal(unname(larg_rel_f), unname(2 * 1.96 / sqrt(D_f[sel_ic, "2022"])))), "\n")

saveRDS(list(mu_h = mu_h, mu_f = mu_f,
             ic_bas_h = ic_bas_h, ic_haut_h = ic_haut_h,
             ic_bas_f = ic_bas_f, ic_haut_f = ic_haut_f),
        file = "resultats/01_taux_mle.rds")
cat("\nSauvé : resultats/01_taux_mle.rds\n")

# ---- Figures A.1 (échelle log) -------------------------------
# En échelle log, les valeurs <= 0 sont masquées en NA : µ̂ = 0 quand D = 0
# (jeunes âges de certaines années) et bornes basses d'IC <= 0 aux grands âges
# (approximation normale non contrainte à R+) — attendu, pas un bug.
masque_log <- function(v) { v[!is.na(v) & v <= 0] <- NA; v }

y_ic <- c(ic_bas_h[, "2022"], ic_haut_h[, "2022"],
          ic_bas_f[, "2022"], ic_haut_f[, "2022"],
          mu_h[, "2022"], mu_f[, "2022"])
y_ic <- y_ic[is.finite(y_ic) & y_ic > 0]

png("figures/fig_A1_taux_2022_ic.png", width = 1400, height = 900, res = 150)
plot(ages, masque_log(mu_h[, "2022"]), log = "y", type = "l", lwd = 2,
     col = "steelblue", ylim = range(y_ic), xlab = "Âge x",
     ylab = expression(hat(mu)[x](2022)),
     main = "A.1 — Taux bruts MLE 2022 et IC 95 % (Autriche)")
lines(ages, masque_log(mu_f[, "2022"]), col = "firebrick", lwd = 2)
lines(ages, masque_log(ic_bas_h[, "2022"]),  col = "steelblue", lty = 2)
lines(ages, masque_log(ic_haut_h[, "2022"]), col = "steelblue", lty = 2)
lines(ages, masque_log(ic_bas_f[, "2022"]),  col = "firebrick", lty = 2)
lines(ages, masque_log(ic_haut_f[, "2022"]), col = "firebrick", lty = 2)
legend("topleft", c("Hommes", "Femmes", "IC 95 %"),
       col = c("steelblue", "firebrick", "grey40"), lty = c(1, 1, 2), lwd = c(2, 2, 1))
dev.off()

# Évolution du profil par âge sur 75 ans (pas de 25 ans)
ans_comp <- c("1947", "1972", "1997", "2022")
couleurs_ans <- c("grey75", "grey55", "grey30", "red")
y_evo <- c(mu_h[, ans_comp], mu_f[, ans_comp])
y_evo <- y_evo[is.finite(y_evo) & y_evo > 0]

png("figures/fig_A1_taux_evolution.png", width = 2000, height = 950, res = 150)
par(mfrow = c(1, 2))
matplot(ages, masque_log(mu_h[, ans_comp]), log = "y", type = "l", lty = 1, lwd = 2,
        col = couleurs_ans, ylim = range(y_evo), xlab = "Âge x",
        ylab = expression(hat(mu)[x](t)), main = "Hommes")
legend("topleft", ans_comp, col = couleurs_ans, lty = 1, lwd = 2)
matplot(ages, masque_log(mu_f[, ans_comp]), log = "y", type = "l", lty = 1, lwd = 2,
        col = couleurs_ans, ylim = range(y_evo), xlab = "Âge x",
        ylab = expression(hat(mu)[x](t)), main = "Femmes")
legend("topleft", ans_comp, col = couleurs_ans, lty = 1, lwd = 2)
par(mfrow = c(1, 1))
dev.off()
cat("Figures A.1 sauvées : fig_A1_taux_2022_ic.png, fig_A1_taux_evolution.png\n")

# ==============================================================
# Section A.2 — Indicateurs de survie périodiques (matière rapport)
# ==============================================================
# ℓ_x calculée jusqu'à l'âge 102 ; µ̂ défini jusqu'à 101 pour TOUTES les années
# et les deux sexes (ETR = 0 dès l'âge 102 chez les femmes en 1949-1951, dès
# l'âge 103 chez les hommes) ; clôture implicite : les survivants décèdent à 102.
# NB : 102 est le dernier âge de la TABLE, pas le dernier âge disposant d'un
# µ̂ : ℓ_102 = exp(−Σ_{k<102} µ_k) ne mobilise que µ_0, ..., µ_101 — le code
# n'indexe donc jamais mu au-delà de l'âge 101 (sinon NaN).
# Justification de la troncature : en survie absolue, ℓ_x < 1e-4 au-delà de
# 100 -> impact négligeable sur e_0 ; MAIS e_65 repose sur la survie
# CONDITIONNELLE ℓ_x/ℓ_65 (survivre de 65 à 100 ans ~ quelques % chez les
# femmes) -> sensibilité à la borne chiffrée plus bas.
age_max_table <- 102

# table_survie : ℓ_x = exp(−Σ_{k<x} µ_k) = Π_{k<x} e^{−µ_k}
# (force de mortalité constante par morceaux), x = 0, ..., age_max ; ℓ_0 = 1.
table_survie <- function(mu_annee, age_max = age_max_table) {
  lx <- c(1, cumprod(exp(-mu_annee[as.character(0:(age_max - 1))])))
  names(lx) <- 0:age_max
  lx
}

# esperance_vie : espérance de vie complète ≈ curtate + ½ :
#   e_x ≈ Σ_{k≥1} kp_x + ½   (hypothèse de décès uniformes dans l'année).
# Alternative exacte sous force constante par morceaux (non retenue ici) :
# contribution de l'année y à e_x : ℓ_y·(1 − e^{−µ_y})/µ_y, divisée par ℓ_x.
esperance_vie <- function(lx, age = 0) {
  px <- lx[as.integer(names(lx)) >= age] / lx[as.character(age)]  # kp_age, k = 0, 1, ...
  sum(px[-1]) + 0.5
}

# age_quantile : âge auquel la survie ℓ_x traverse le niveau p (interpolation
# linéaire entre âges entiers) -> quantiles de l'âge au décès :
# p = 0.50 médiane ; p = 0.75 -> Q1 ; p = 0.25 -> Q3.
age_quantile <- function(lx, p) {
  pos <- which(lx <= p)[1]
  x <- as.integer(names(lx))
  unname(x[pos - 1] + (lx[pos - 1] - p) / (lx[pos - 1] - lx[pos]))
}

# Tables de survie de toutes les années (ℓ_x en lignes, années en colonnes)
lx_h <- sapply(colnames(mu_h), function(an) table_survie(mu_h[, an]))
lx_f <- sapply(colnames(mu_f), function(an) table_survie(mu_f[, an]))

# Indicateurs par année : e_0, e_65, médiane et quartiles de l'âge au décès
e0_h  <- apply(lx_h, 2, esperance_vie)
e0_f  <- apply(lx_f, 2, esperance_vie)
e65_h <- apply(lx_h, 2, esperance_vie, age = 65)
e65_f <- apply(lx_f, 2, esperance_vie, age = 65)
mediane_h <- apply(lx_h, 2, age_quantile, p = 0.50)
mediane_f <- apply(lx_f, 2, age_quantile, p = 0.50)
q1_h <- apply(lx_h, 2, age_quantile, p = 0.75)   # 25 % de la cohorte fictive décédée
q1_f <- apply(lx_f, 2, age_quantile, p = 0.75)
q3_h <- apply(lx_h, 2, age_quantile, p = 0.25)   # 75 % décédée
q3_f <- apply(lx_f, 2, age_quantile, p = 0.25)

# ---- Contrôles A.2 -------------------------------------------
# Attendu : e_0, e_65 et médiane en hausse (expansion), IQR en baisse
# (compression / rectangularisation).
ans_ctrl <- c("1947", "2023")
cat("\nIndicateurs périodiques 1947 vs 2023 :\n")
print(round(rbind(e0_h  = e0_h[ans_ctrl],  e0_f  = e0_f[ans_ctrl],
                  e65_h = e65_h[ans_ctrl], e65_f = e65_f[ans_ctrl],
                  mediane_h = mediane_h[ans_ctrl], mediane_f = mediane_f[ans_ctrl],
                  iqr_h = (q3_h - q1_h)[ans_ctrl], iqr_f = (q3_f - q1_f)[ans_ctrl]), 2))

# Sensibilité au choix de la borne de table (e_65 est le plus exposé :
# survie conditionnelle) — année 2023, où µ̂ est défini jusqu'à 102 pour les
# deux sexes, ce qui permet de tester une borne jusqu'à 103 :
cat("\nSensibilité à la troncature de la table (année 2023) :\n")
for (b in c(100, 101, 102, 103)) {
  lx_hb <- table_survie(mu_h[, "2023"], age_max = b)
  lx_fb <- table_survie(mu_f[, "2023"], age_max = b)
  cat("  borne", b, ": e0 H =", sprintf("%.3f", esperance_vie(lx_hb)),
      "F =", sprintf("%.3f", esperance_vie(lx_fb)),
      "| e65 H =", sprintf("%.3f", esperance_vie(lx_hb, 65)),
      "F =", sprintf("%.3f", esperance_vie(lx_fb, 65)), "\n")
}

# ---- Figures A.2 ---------------------------------------------
# Rectangularisation : ℓ_x se rapproche d'un rectangle au fil des années
png("figures/fig_A2_courbes_survie.png", width = 2000, height = 950, res = 150)
par(mfrow = c(1, 2))
matplot(0:age_max_table, lx_h[, ans_comp], type = "l", lty = 1, lwd = 2,
        col = couleurs_ans, xlab = "Âge x", ylab = expression(l[x]), main = "Hommes")
legend("bottomleft", ans_comp, col = couleurs_ans, lty = 1, lwd = 2)
matplot(0:age_max_table, lx_f[, ans_comp], type = "l", lty = 1, lwd = 2,
        col = couleurs_ans, xlab = "Âge x", ylab = expression(l[x]), main = "Femmes")
legend("bottomleft", ans_comp, col = couleurs_ans, lty = 1, lwd = 2)
par(mfrow = c(1, 1))
dev.off()

# d_x = ℓ_x − ℓ_{x+1} : le mode adulte se déplace vers la droite (expansion)
# et le pic se resserre (compression)
dx_h <- apply(lx_h[, ans_comp], 2, function(l) -diff(l))
dx_f <- apply(lx_f[, ans_comp], 2, function(l) -diff(l))
png("figures/fig_A2_courbe_deces.png", width = 2000, height = 950, res = 150)
par(mfrow = c(1, 2))
matplot(0:(age_max_table - 1), dx_h, type = "l", lty = 1, lwd = 2,
        col = couleurs_ans, xlab = "Âge x", ylab = expression(d[x]), main = "Hommes")
legend("topleft", ans_comp, col = couleurs_ans, lty = 1, lwd = 2)
matplot(0:(age_max_table - 1), dx_f, type = "l", lty = 1, lwd = 2,
        col = couleurs_ans, xlab = "Âge x", ylab = expression(d[x]), main = "Femmes")
legend("topleft", ans_comp, col = couleurs_ans, lty = 1, lwd = 2)
par(mfrow = c(1, 1))
dev.off()

png("figures/fig_A2_esperance_vie.png", width = 2000, height = 950, res = 150)
par(mfrow = c(1, 2))
plot(annees, e0_h, type = "l", lwd = 2, col = "steelblue",
     ylim = range(e0_h, e0_f), xlab = "Année t", ylab = expression(e[0](t)),
     main = "Espérance de vie à la naissance")
lines(annees, e0_f, lwd = 2, col = "firebrick")
legend("bottomright", c("Hommes", "Femmes"), col = c("steelblue", "firebrick"), lwd = 2)
plot(annees, e65_h, type = "l", lwd = 2, col = "steelblue",
     ylim = range(e65_h, e65_f), xlab = "Année t", ylab = expression(e[65](t)),
     main = "Espérance de vie à 65 ans")
lines(annees, e65_f, lwd = 2, col = "firebrick")
legend("bottomright", c("Hommes", "Femmes"), col = c("steelblue", "firebrick"), lwd = 2)
par(mfrow = c(1, 1))
dev.off()

# Médiane et bande interquartile [Q1, Q3] de l'âge au décès
ylim_q <- range(q1_h, q3_h, q1_f, q3_f)
png("figures/fig_A2_mediane_iqr.png", width = 2000, height = 950, res = 150)
par(mfrow = c(1, 2))
plot(annees, mediane_h, type = "n", ylim = ylim_q, xlab = "Année t",
     ylab = "Âge au décès", main = "Hommes")
polygon(c(annees, rev(annees)), c(q1_h, rev(q3_h)), col = "grey88", border = NA)
lines(annees, mediane_h, lwd = 2)
lines(annees, q1_h, lty = 2, col = "grey40")
lines(annees, q3_h, lty = 2, col = "grey40")
legend("bottomright", c("Médiane", "Q1-Q3"), col = c("black", "grey40"),
       lty = c(1, 2), lwd = c(2, 1))
plot(annees, mediane_f, type = "n", ylim = ylim_q, xlab = "Année t",
     ylab = "Âge au décès", main = "Femmes")
polygon(c(annees, rev(annees)), c(q1_f, rev(q3_f)), col = "grey88", border = NA)
lines(annees, mediane_f, lwd = 2)
lines(annees, q1_f, lty = 2, col = "grey40")
lines(annees, q3_f, lty = 2, col = "grey40")
legend("bottomright", c("Médiane", "Q1-Q3"), col = c("black", "grey40"),
       lty = c(1, 2), lwd = c(2, 1))
par(mfrow = c(1, 1))
dev.off()
cat("Figures A.2 sauvées : fig_A2_courbes_survie.png, fig_A2_courbe_deces.png,",
    "fig_A2_esperance_vie.png, fig_A2_mediane_iqr.png\n")

saveRDS(list(annees = annees,
             e0_h = e0_h, e0_f = e0_f, e65_h = e65_h, e65_f = e65_f,
             mediane_h = mediane_h, mediane_f = mediane_f,
             q1_h = q1_h, q1_f = q1_f, q3_h = q3_h, q3_f = q3_f),
        file = "resultats/01_indicateurs.rds")
cat("Sauvé : resultats/01_indicateurs.rds\n")
