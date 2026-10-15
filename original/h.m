close all; clear; clc;

load('signal_recu.mat');

r = signal_recu;

% Le bon bloc periodique
debut = 20000;
fin   = 26000;
bloc  = r(debut:fin);
N     = length(bloc);
f     = (-N/2 : N/2-1) / N;

% Spectre du bloc
R   = fftshift(fft(bloc));
psd = abs(R).^2;

figure;
plot(f, 10*log10(psd));
grid on;
xlabel('frequence normalisee');
ylabel('dB');
title('Spectre du bloc  forme de h(t)');
xlim([-0.3, 0.3]);

