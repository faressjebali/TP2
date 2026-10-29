"""Partie 4 - Synthese : PESQ des estimees et figure p4_q7_synthese.png

Le PESQ MATLAB de Loizou (P-code d'avant 2007) ne s'execute plus dans MATLAB.
On utilise donc le package Python `pesq` (code C de reference UIT-T P.862).

Etapes :
  1. MATLAB : lancer Q7_export.m  -> pesq_wav/*.wav et pesq_wav/gains.csv
  2. Python : pip install numpy scipy matplotlib pesq
              python pesq_synthese.py
     (sous Windows, `pesq` doit etre compile : il faut Visual Studio Build Tools ;
      sinon, executer ce script sur Google Colab ou Linux/macOS)
"""
from pathlib import Path
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.io import wavfile
from pesq import pesq

DOSSIER = Path(__file__).resolve().parent
WAV = DOSSIER / "pesq_wav"
FIGS = DOSSIER / "figs_partie4"

# Nom des fichiers (Q7_export.m) -> etiquette de la figure
METHODES = ["melange", "global", "trames", "trames_ola", "dictionnaire", "ar_iteratif"]
ETIQUETTES = ["Mélange", "Global", "Trames", "Trames\n+ OLA",
              "Dictionnaire\n($N_c=8$)", "AR itératif\n($K=1$)"]
RD1_ESTIME = [False, False, False, False, True, True]   # r_d1 estime (rouge) ou connu (bleu)


def lire(fichier):
    """Lit un .wav (16 bits) et le renvoie en flottants dans [-1, 1]."""
    fs, s = wavfile.read(fichier)
    if s.dtype == np.int16:
        s = s / 32768.0
    return fs, s.astype(np.float64)


def lire_gains(fichier):
    with open(fichier, newline="") as f:
        return {ligne["methode"]: float(ligne["gain_dB"]) for ligne in csv.DictReader(f)}


def main():
    fs, ref_h = lire(WAV / "ref_homme.wav")
    _, ref_f = lire(WAV / "ref_femme.wav")
    gains = lire_gains(WAV / "gains.csv")

    G, Q_h, Q_f = [], [], []
    print(f"{'methode':<13} {'gain (dB)':>9} {'PESQ homme':>11} {'PESQ femme':>11}")
    for m in METHODES:
        _, h = lire(WAV / f"homme_{m}.wav")
        _, f = lire(WAV / f"femme_{m}.wav")
        G.append(gains[m])
        Q_h.append(pesq(fs, ref_h, h, "nb"))   # bande etroite : Fs = 8 kHz
        Q_f.append(pesq(fs, ref_f, f, "nb"))
        print(f"{m:<13} {G[-1]:9.2f} {Q_h[-1]:11.2f} {Q_f[-1]:11.2f}")

    # ----- Figure de synthese -----
    x = np.arange(len(METHODES))
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True)

    couleurs = ["#c0392b" if e else "#1f4e9c" for e in RD1_ESTIME]
    barres = ax1.bar(x, G, width=0.8, color=couleurs)
    ax1.bar_label(barres, fmt="%.2f", padding=3, fontsize=9)
    ax1.set_ylabel("gain en RSB (dB)")
    ax1.set_ylim(0, 1.15 * max(G))
    ax1.grid(axis="y", alpha=0.4)
    ax1.set_title(r"Synthèse des méthodes ($P=50$, $L=40$ ms) — "
                  r"bleu : $r_{d_1}$ connu, rouge : $r_{d_1}$ estimé")

    l = 0.38
    bh = ax2.bar(x - l / 2, Q_h, width=l, color="#1f77b4", label="homme")
    bf = ax2.bar(x + l / 2, Q_f, width=l, color="#e377c2", label="femme")
    ax2.bar_label(bh, fmt="%.2f", padding=2, fontsize=8)
    ax2.bar_label(bf, fmt="%.2f", padding=2, fontsize=8)
    ax2.axhline(Q_h[0], color="#1f77b4", ls=":", lw=1)   # melange non traite
    ax2.axhline(Q_f[0], color="#e377c2", ls=":", lw=1)
    ax2.set_ylabel("PESQ (MOS)")
    ax2.set_ylim(1, max(Q_h + Q_f) + 0.5)
    ax2.grid(axis="y", alpha=0.4)
    ax2.legend(loc="upper right")
    ax2.set_xticks(x, ETIQUETTES)

    fig.tight_layout()
    FIGS.mkdir(exist_ok=True)
    sortie = FIGS / "p4_q7_synthese.png"
    fig.savefig(sortie, dpi=150)
    print(f"Figure enregistree : {sortie}")


if __name__ == "__main__":
    main()
