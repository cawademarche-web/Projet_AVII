# Projet ACTU-F502 — Assurance-Vie II (ULB)

Projet de cours de Master en actuariat : estimer et projeter la mortalité de
deux cohortes (H/F, 65 ans en 2022) pour tarifer des rentes et calculer un SCR.
Consignes complètes : `docs/brief_projet.md`. Signaux de contrôle : `docs/checklist.md`.

## Contexte de qualité — LIRE EN PREMIER
Ce code sera défendu oralement par un étudiant qui doit pouvoir expliquer
chaque ligne. Le style cible est "très bon travail d'étudiant", PAS
"production". Concrètement :
- Scripts R linéaires, exécutables de haut en bas. Pas de package maison,
  pas de classes R6/S4, pas de config externe, pas de CLI.
- Pas de programmation défensive : pas de tryCatch, pas de validation
  d'arguments, pas de gestion d'erreurs au-delà du strict nécessaire.
- Une fonction = un concept actuariel, nommée en français :
  `vap_rente()`, `simule_kappa()`, `calcule_BE()`, `calcule_SCR()`...
  Fonctions courtes (< ~25 lignes). Si une fonction n'est appelée qu'une
  fois et ne porte pas un concept du cours, inliner.
- Commentaires en français, reliant le code à la notation du cours :
  "# ln µ_x(t) = α_x + β_x·κ_t (Lee-Carter)".
- Vectoriser par matrices (trajectoires en colonnes) plutôt que boucles
  imbriquées, mais rien de plus sophistiqué.
- Packages autorisés : base R, StMoMo, demography (+ ggplot2 si utile,
  sinon graphics de base). RIEN d'autre sans demander.

## Décision méthodologique (verrouillée)
- Sections A et B : StMoMo (fit, forecast, bootstrap).
- Sections C et D : code 100 % custom (StMoMo n'a pas de fonctions
  pricing/solvabilité).

## Les trois taux — NE JAMAIS confondre
Noms de variables imposés, interdiction d'utiliser `r`, `i` ou `t` seuls :
- `taux_tarif  <- 0.03`  # tarification : VAP, primes (sections C). Ne
                          # réapparaît PLUS dans la boucle SCR.
- `taux_actu_be <- 0.02` # actualise le BE à chaque date ET le BOF_2023
                          # dans la formule du SCR (section D).
- `rdt_actifs  <- 0.04`  # fait croître les ACTIFS uniquement, une fois
                          # par an (section D). Sensibilité : 0.01.

## Convention de rente (verrouillée)
Les rentes du projet sont à terme échu (paiement en fin de période).
Donc `vap_rente()` somme à partir de k = 1, jamais k = 0 :
a_65 = Σ_{k≥1} v^k · ₖp₆₅. La rente [B] est temporaire 15 ans à terme
échu : somme de k = 1 à 15. Ne jamais écrire « rente immédiate » dans
le code ou les commentaires, c'est ambigu.

## Indexation des matrices
Les matrices âges × années ont des dimnames en caractères, ligne 1 =
âge 0. Toujours indexer par nom (`M[as.character(65), ]`), jamais par
position (`M[65, ]` renvoie l'âge 64). Cette règle vaut pour tous les
scripts, en particulier la diagonale de cohorte 65 → 66 → 67 des
sections B, C et D.

## Contrat d'interface entre sections
Chaque script sauve ses objets clés via saveRDS() dans `resultats/` ;
le script suivant les CHARGE, il ne recalcule jamais.
En particulier : la section D réutilise les 5000 trajectoires de κ
simulées en section B (`resultats/trajectoires_kappa.rds`). Interdit de
re-simuler.
- 01 → taux bruts, expositions      · 02 → fit LC, bootstrap
- 03 → trajectoires κ (5000)        · 04 → VAP et primes par trajectoire
- 05 → BOF, SCR

## Structure du repo
scripts/01_donnees_taux.R · 02_leecarter_stmomo.R · 03_simulation_kappa.R ·
04_pricing_vap.R · 05_solvabilite_scr.R
resultats/ (RDS) · figures/ (PNG, nommés fig_<section>_<contenu>.png) ·
docs/ (brief, checklist, PIPELINE.md) · donnees/ (HMD)

## Reproductibilité
`set.seed(1234)` en tête de chaque script qui simule. Toute quantité
aléatoire doit être identique d'un run à l'autre.

## Fin de chaque section (obligatoire)
1. Exécuter réellement le script (`Rscript scripts/0X_....R`) et montrer
   la sortie console.
2. Imprimer les valeurs de contrôle demandées par `docs/checklist.md`
   pour la section.
3. Mettre à jour `docs/PIPELINE.md` : section → script → méthode →
   outputs (chiffres clés) → TODO/choix à justifier.
   RÈGLE STRICTE : aucune valeur numérique n'entre dans PIPELINE.md
   autrement que **recopiée d'une sortie console réelle** du run qui vient
   d'être exécuté. Interdiction d'anticiper une valeur, de la reprendre d'un
   plan, d'un run antérieur ou d'une estimation. Si un chiffre attendu diverge
   du chiffre obtenu, c'est une information à remonter, pas à lisser.
4. Lister les TODO "à justifier dans le rapport" plutôt que de trancher
   silencieusement (ex. plage d'âges de calibration).
5. Commit avec message clair.

## Compact Instructions
En cas de résumé de conversation : préserver les valeurs numériques de
contrôle obtenues, les chemins des RDS produits, et la liste des TODO.

## Notes de référence
- `docs/STMOMO_NOTES.md` — chaîne d'appels StMoMo (quoi appeler, avec quels arguments).
- `docs/BOOTSTRAP_NOTES.md` — fondements du bootstrap semi-paramétrique
  (Brouhns-Denuit-Van Keilegom 2005 ; Brouhns-Denuit-Vermunt 2002 ;
  Haberman-Renshaw 2009). Sert à : justifier le choix méthodologique,
  fixer les arguments de bootstrap(), et fournir des ORDRES DE GRANDEUR
  de contrôle post-bootstrap.
  PÉRIMÈTRE STRICT : cette note documente des alternatives (bootstrap
  paramétrique, bootstrap résiduel, bootstrap par génération, modèle
  relationnel de Brass, modèles APC). AUCUNE ne doit être implémentée.
  Le projet ne demande que Lee-Carter + bootstrap semi-paramétrique.
  Ces alternatives sont matière de RAPPORT et de DÉFENSE, pas de code.