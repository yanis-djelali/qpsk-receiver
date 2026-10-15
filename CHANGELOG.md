# Historique

## Mars 2026 : projet TS218

Scripts d'exploration et versions successives de la chaîne, conservés tels quels dans `original/` :

- `estimdeltaf.m`, `deltaf.m`, `estim2deltaf.m` : estimation de la dérive de fréquence, essais de gains de PLL ;
- `affichage1_Ts.m`, `h.m`, `trouverhh.m`, `estimeroffset.m` : rythme symbole, filtre adapté, choix de l'offset ;
- `detectiontrame.m`, `pipeline.m`, `pipeline2.m` : chaîne complète avec décision par angle et test des 8 combinaisons.

La correction de la décision (rotation de π/4 et lecture des signes) a été trouvée en fin de projet et présentée dans `docs/presentation.pdf`.

## Octobre 2026 : mise en forme du dépôt

- `recepteur_qpsk.m` : un seul script reprenant la chaîne de la présentation, avec la décision corrigée et le choix automatique de l'ambiguïté de rotation (critère de régularité de l'image) ;
- conversion de bits en pixels sans `bi2de`, pour ne pas dépendre de la Communications Toolbox ;
- `verification_python/` : réplique Python, exécutée sur le vrai signal ;
- README, figures, présentation.
