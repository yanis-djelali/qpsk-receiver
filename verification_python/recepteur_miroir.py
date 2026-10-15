"""Réplique Python de recepteur_qpsk.m, étape par étape.

Sert à vérifier la chaîne sans MATLAB. Usage :
    python verification_python/recepteur_miroir.py data/signal_recu.mat
Produit docs/img/verification_python.png.
"""
import sys

import matplotlib
import numpy as np
import scipy.io as sio
from scipy.signal import welch

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def recepteur(r, alpha=0.1, beta=1e-4, T=128):
    N = len(r)

    # 1) dérive de fréquence par blocs
    n_bloc, n_blocs = 1000, 700
    phase = np.array([np.angle(np.sum(r[i * n_bloc:(i + 1) * n_bloc] ** 4)) for i in range(n_blocs)])
    df_blocs = np.diff(np.unwrap(phase)) / (4 * 2 * np.pi * n_bloc)

    # 2) filtre adapté (pwelch 256/128/512, fenêtre de Hamming, spectre centré)
    _, pxx = welch(r[9999:25000], fs=1, window="hamming", nperseg=256, noverlap=128,
                   nfft=512, return_onesided=False, detrend=False)
    pxx = np.fft.fftshift(pxx)
    h = np.real(np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(np.sqrt(pxx)))))
    y = np.convolve(r, np.conj(h[::-1]))[256:256 + N]       # conv(..., 'same') de MATLAB

    # 3) échantillonnage, Ts = 10 Te, offset 9
    symboles = y[9::10]

    # 4) PLL d'ordre 2
    t = np.arange(1, len(df_blocs) + 1) * 100
    delta = (df_blocs[0] + (df_blocs[1] - df_blocs[0]) * (1 - t[0]) / (t[1] - t[0])) * 2 * np.pi
    phi = 0.0
    s_corr = np.empty_like(symboles)
    for k, sym in enumerate(symboles):
        s = sym * np.exp(-1j * phi)
        s_corr[k] = s
        e = np.imag((s / abs(s)) ** 4)
        delta += beta * e
        phi += delta + alpha * e

    # 5) trame utile
    actifs = np.flatnonzero(np.abs(s_corr) > 0.2 * np.abs(s_corr).max())
    utile = s_corr[actifs[0]:actifs[-1] + 1]

    # 6) rotation pi/4, décision par signe, choix de l'ambiguïté
    meilleur = None
    for permute in (0, 1):
        for rot in range(4):
            s = utile * np.exp(1j * np.pi / 4) * np.exp(1j * rot * np.pi / 2)
            b1, b2 = (s.real < 0).astype(int), (s.imag < 0).astype(int)
            if permute:
                b1, b2 = b2, b1
            img = bits_vers_image(np.column_stack([b1, b2]).ravel(), T)
            score = min(correlation_voisins(img, 0), correlation_voisins(img, 1))
            if meilleur is None or score > meilleur[0]:
                meilleur = (score, rot, permute, img)
    return actifs[0] + 1, actifs[-1] + 1, utile, meilleur


def bits_vers_image(bits, T):
    n = T * T
    octets = bits[:n * 8].reshape(8, n).T                  # reshape(..., [], 8) en ordre colonne
    pixels = octets @ (2 ** np.arange(8))
    return pixels.reshape(T, T, order="F")


def correlation_voisins(img, axe):
    x = img.astype(float)
    a, b = (x[:-1, :], x[1:, :]) if axe == 0 else (x[:, :-1], x[:, 1:])
    return np.corrcoef(a.ravel(), b.ravel())[0, 1]


if __name__ == "__main__":
    r = sio.loadmat(sys.argv[1])["signal_recu"].ravel()
    debut, fin, utile, (score, rot, permute, img) = recepteur(r)
    print(f"Trame : symboles {debut} à {fin} ({len(utile)} symboles)")
    print(f"Rotation retenue : {rot} x pi/2, bits permutés : {permute}, régularité {score:.3f}")

    fig, ax = plt.subplots(1, 2, figsize=(8, 4))
    z = utile[::10]
    ax[0].plot(z.real, z.imag, ".", ms=1)
    ax[0].set_aspect("equal")
    ax[0].set_title("Constellation après PLL")
    ax[1].imshow(img, cmap="gray")
    ax[1].set_title("Image décodée")
    ax[1].axis("off")
    plt.tight_layout()
    plt.savefig("docs/img/verification_python.png", dpi=110)
