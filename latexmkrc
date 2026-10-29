# Compilation : latexmk  (ou latexmk -pdf main.tex)
# Les fichiers intermediaires (.aux, .log, .toc...) vont dans build/,
# main.pdf reste a la racine.
$pdf_mode = 1;
$aux_dir = 'build';
@default_files = ('main.tex');
