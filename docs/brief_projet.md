# ACTUF502 – Assurance Vie II

**Antoine DELWARDE**

---

## Consignes

- **Pondération**
  - Le travail compte pour 50% de la note finale
  - Pondération rapport/présentation : rapport 50%, présentation 50%
- **Deadline** : 24/04/2026 en session de juin et 21/08/2026 en session d'août, à télécharger sur l'UV
- **Présentation** : Le travail sera défendu en session par le groupe – dates à définir
- **Langage de programmation et package** : Vous pouvez utiliser le package R StMoMo pour l'estimation, projection et diagnostic (résidus, bootstrapping) ou coder vous-mêmes. Le développement personnel du code est avantagé comme suit :
  - Développer le code soi-même permet d'atteindre la cote de 20/20 pour le rapport écrit
  - L'utilisation du package StMoMo plafonne la cote du rapport écrit à 16/20
  - Dans le cas d'une combinaison du package StMoMo et d'un code personnel, la cote maximale du rapport écrit est adaptée proportionnellement
- **Documents à remettre** : rapport écrit en version .doc/.pdf (à faire en Word ou en LaTeX), ainsi que le code et données associés. Pas de fichier "zip". Ces documents doivent être téléchargés séparément sur l'UV.
- **Longueur du rapport écrit** : maximum 30 pages, y inclus les graphiques. Le code n'est pas compris dans ce maximum.

---

## Sujet

On se propose d'estimer et de projeter la mortalité de deux cohortes d'assuré.e.s (hommes et femmes) afin de calculer la valeur actuelle probable d'un produit de rente viagère. Plus précisément, les cohortes concernent des assuré.e.s qui ont contracté en 2022 à l'âge de 65 ans.

### Données

Les données de mortalité générale sont disponibles dans la Human Mortality Database (http://www.mortality.org/) et doivent être téléchargées par le groupe. À chaque groupe est attribué un pays différent.

---

## A. Analyse de la mortalité de la population générale

1. **Estimation** : En utilisant la procédure vue en cours pour estimer les taux de mortalité par maximum de vraisemblance, calculer les taux de mortalité en fonction de l'âge pour les cohortes concernées. Illustrer ces taux de mortalité en fonction de l'âge pour les cohortes concernées à partir de 2022. Tracer également les intervalles de confiance au seuil de 95% pour ces cohortes.

2. **Discussion** : Discuter, pour chacune des cohortes, les $\hat{\mu}_x(t)$ estimés, notamment :
   - i. Tendances associées à l'espérance de vie périodique
   - ii. Tendances associées à l'âge médian et les interquartiles (IQR)
   - iii. Phénomène d'expansion et rectangularisation

3. **Estimation, projection et diagnostic** : Estimer les paramètres d'un modèle de Lee-Carter à partir des données historiques téléchargées. Prendre bien soin de :
   - i. Commenter et justifier le choix de la plage d'âges et de la période choisies pour la calibration
   - ii. Commenter les résultats obtenus en affichant les paramètres estimés
   - iii. Comparer les (log-)taux de mortalité estimés par maximum de vraisemblance (question A.1) avec les (log-)taux de mortalité estimés du modèle de Lee-Carter pour les cohortes concernées
   - iv. Afficher les résidus et discuter les résultats
   - v. Projeter les (log-)taux de mortalité
   - vi. Calculer des intervalles de prédiction au seuil de 95% par bootstrapping, avec un échantillon de taille $N = 5000$. Justifier le choix de la méthode. Afficher l'histogramme des espérances de vie cohorte à l'âge de 65 ans pour les cohortes concernées.

---

## B. Analyse de la mortalité des assurés du portefeuille

Maintenant on s'intéresse exclusivement aux assuré.e.s du portefeuille :

1. Afficher les (log-)taux de mortalité historiques et projetés pour les cohortes concernées à partir de 2022.

2. Simuler $N = 5000$ trajectoires projetées des taux de mortalité futurs pour les cohortes concernées. Calculer des intervalles de prédiction au seuil de 95%. Ici on s'intéresse seulement à l'incertitude liée à la projection du paramètre temporel $\kappa_t$ (et pas à l'incertitude obtenue dans la méthodologie du bootstrapping).

3. Comparer les intervalles de prédiction obtenus avec bootstrapping (section A.3.vi) avec ceux basés exclusivement sur la projection des $\kappa_t$ (section B.2).

---

## C. Tarification

1. Calculer la Valeur Actuarielle Présente (VAP) de contrats souscrits en 2022 :
   - i. Dans le cas de rentes viagères à terme échu (rente **[A]**)
   - ii. Dans le cas de rentes viagères à terme échu temporaires d'une durée de 15 ans (rente **[B]**)

2. Discuter les résultats avec un taux technique $t = 3\%$. Donner la valeur moyenne obtenue et sa variance pour les deux contrats. Quelles sont les autres sources d'incertitude ? Discuter.

3. Étudier la variation de la VAP en fonction du taux technique pour $t = 1\%, 2\%, 4\%, 5\%$.

4. Déterminer la prime des 2 types de rentes (A et B) par le principe d'équivalence, avec un taux technique $t = 3\%$.

---

## D. Solvabilité

Supposons les portefeuilles suivants :

- **Portefeuille [A]** : rentes [A], 500 femmes et 500 hommes
- **Portefeuille [B]** : rentes [B], 500 femmes et 500 hommes

Le Solvency Capital Requirement ($SCR$) est défini comme le percentile 0.5% de la variation de fonds propres ($BOF$) sur un horizon d'un an. En d'autres termes, l'assureur doit détenir suffisamment de fonds propres pour couvrir ses pertes potentielles avec une probabilité de 99.5% :

$$P(BOF_{t+1} > 0 \mid BOF_t = SCR_t) = 99.5\%$$

avec :

- $BOF_t = A_t - BE_t$ les fonds propres au temps $t$
- $A_t$ les actifs au temps $t$
- $BE_t$ la réserve mathématique au temps $t$ (c'est-à-dire les cashflows restants, actualisés au taux de marché $r = 2\%$)

L'actif initial $A_{2022}$ correspond donc aux primes uniques calculées à la section C.4. L'assureur investit ces primes dans un fonds dont le rendement annuel est $i = 4\%$, sans risque.

Dans cet exercice, on calcule le $SCR$ de la façon suivante :

$$SCR_t = VaR_{99.5\%}(-\Delta BOF) = VaR_{99.5\%}\left(BOF_t - \frac{BOF_{t+1}}{1+r}\right)$$

**Attention** : toutes les valeurs doivent être calculées en $t^+$, c'est-à-dire après paiement des primes et prestations qui ont lieu en $t$ :

1. Pour le portefeuille [A], calculer $A_{2022}$ et $BE_{2022}$. En déduire $BOF_{2022}$.

2. Pour le portefeuille [A], et pour chacune des $N = 5000$ trajectoires définies en section B.2, définir l'état du portefeuille en $t = 2023$ : combien d'assuré.e.s sont encore en vie, combien d'assuré.e.s sont décédés ? En déduire $A_{2023}$, $BE_{2023}$ et $BOF_{2023}$ pour chaque trajectoire. En déduire $SCR_{2022}$.

3. Recommencer les étapes 1 et 2 pour le portefeuille [B]. Comparer les résultats.

4. Calculer $SCR_{2023}$ pour les portefeuilles [A] et [B]. Comparer $SCR_{2022}$ et $SCR_{2023}$ pour les 2 portefeuilles.

5. Que deviendraient les résultats 1, 2, 3, 4 si l'assureur n'obtient qu'un rendement sans risque de $i = 1\%$ sur ses actifs ? Discuter.
