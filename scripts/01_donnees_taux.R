# ==============================================================
# Script 01 — Chargement des données HMD (Autriche, 1947-2023)
# Projet ACTU-F502 — Assurance Vie II
#
# Cette étape : lire Deaths_1x1.txt et Exposures_1x1.txt (HMD) et
# construire les matrices D_{x,t} (décès observés) et ETR_{x,t}
# (exposition au risque), hommes et femmes séparés.
# Les taux bruts µ̂_x(t) = D_{x,t} / ETR_{x,t} et leurs IC à 95 %
# (section A1) seront ajoutés plus loin dans ce même script.
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
