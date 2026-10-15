close all; clear; clc;

load('signal_recu.mat');

r = signal_recu(20000:22000);

% Spectre lisse avec pwelch
[pxx, f] = pwelch(r, 256, 128, 512, 1, 'centered');

figure;
subplot(3,1,1);
plot(f, 10*log10(pxx));
grid on;
xlabel('frequence normalisee'); ylabel('dB');
title('Spectre pwelch = |H(f)|^2');

% Racine carree pour avoir |H(f)|
H_module = sqrt(pxx);

subplot(3,1,2);
plot(f, H_module);
grid on;
xlabel('frequence normalisee'); ylabel('amplitude');
title('|H(f)| = racine du spectre');

% IFFT pour retrouver h(t)
H_ifft = ifftshift(H_module);
h = real(ifft(H_ifft));
h = fftshift(h);

subplot(3,1,3);
plot(h);
grid on;
xlabel('echantillon'); ylabel('amplitude');
title('h(t) retrouve par IFFT');