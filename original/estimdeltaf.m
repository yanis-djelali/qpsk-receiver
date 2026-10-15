close all; clear; clc;

load('signal_recu.mat');

r = signal_recu(1:650000);

% On eleve a la puissance 4 pour supprimer la modulation QPSK
z = r .^ 4;

% Phase de z = 4 * phi(t) = 4 * (2*pi*delta_f*t + phi0)
phase = unwrap(angle(z));

% delta_f = derivee de la phase / (4 * 2*pi)
delta_f_brut  = diff(phase) / (4 * 2 * pi);
delta_f_lisse = movmean(delta_f_brut, 10000);

figure;
subplot(2,1,1);
plot(phase);
title('Phase de r^4 (unwrapped)  →  droite = delta f constant, courbe = delta f variable');
xlabel('echantillon'); ylabel('phase (rad)');
grid on;

subplot(2,1,2);
plot(delta_f_lisse);
title('delta f estime au cours du temps');
xlabel('echantillon'); ylabel('delta f normalise');
grid on;
yline(0, 'r--', 'delta f = 0');