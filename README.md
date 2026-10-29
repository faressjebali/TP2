# TP2 TNS — Débruitage par filtre de Wiener

```
TP2/
├── main.pdf                compte rendu final (version montrée à la prof)
├── main.tex                préambule, page de garde, sommaire + \input des parties
├── sections/               texte du rapport, un fichier par partie
│   ├── partie1.tex         Wiener et Widrow
│   ├── partie2.tex         débruitage d'une harmonique
│   ├── partie3.tex         débruitage d'un signal audio
│   ├── partie4.tex         séparation de deux voix
│   ├── conclusion.tex
│   └── references.tex
├── figures/                figures utilisées par le rapport
│   ├── logo.jpg, logo_carthage.jpg
│   ├── partie2/            fig21 … fig36
│   ├── partie3/            p3_*.png
│   └── partie4/            p4_*.png
├── TP2_Matlab/             code MATLAB (même contenu que sur MATLAB Online)
│   ├── Partie2/            Question_a_b.m … Question_f.m
│   ├── Partie3/            Question_a.m … Question_f.m, debruitage_sig.m, monson.wav
│   └── Partie4/            Q1_global.m … Q7_export.m, fonctions, voix OSR, pesq_synthese.py
├── documents/              sujet du TP, sujet de la partie 4, cours
├── latexmkrc               réglages de compilation
└── build/                  fichiers intermédiaires LaTeX (générés, ignorés par git)
```

## Workflow

1. Sur MATLAB Online, se placer dans le dossier de la partie (`cd Partie3`) puis lancer le script.
   Chaque script enregistre ses figures dans `figs_partieN/` à côté de lui.
2. Copier les figures voulues dans `figures/partieN/` (ce sont elles que le rapport utilise).
3. Modifier le texte dans `sections/partieN.tex`.
4. Compiler à la racine : `latexmk` → `main.pdf` (les fichiers intermédiaires vont dans `build/`).

## PESQ (partie 4, synthèse)

Le PESQ MATLAB de Loizou (P-code) ne tourne plus dans MATLAB, d'où un passage par Python :

1. MATLAB : `Q7_export.m` → écrit `pesq_wav/*.wav` et `pesq_wav/gains.csv`.
2. Télécharger le dossier `Partie4/` (avec `pesq_wav/`), puis :
   ```
   pip install numpy scipy matplotlib pesq
   python pesq_synthese.py
   ```
   → affiche le tableau gain / PESQ et crée `figs_partie4/p4_q7_synthese.png`.

Sur ce PC, `pesq` est déjà installé (compilé avec le gcc de Strawberry Perl) :
il suffit de lancer `python pesq_synthese.py` dans `TP2_Matlab/Partie4/`.

