debut = 20000;
fin   = 26000;
bloc  = r(debut:fin);
N     = length(bloc);
f     = (-N/2 : N/2-1) / N;

% Spectre de |bloc|² pour trouver Ts
R2  = fftshift(fft(abs(bloc).^2));
psd = abs(R2);

figure;
subplot(2,1,1);
plot(f, 20*log10(psd));
grid on;
xlabel('frequence normalisee'); ylabel('dB');
title('Spectre de |bloc|² raie a f = 1/Ts');
xlim([0, 0.5]);

% Spectre de bloc^4 pour confirmer QPSK
R4 = fftshift(fft(bloc.^4));

subplot(2,1,2);
plot(f, 20*log10(abs(R4)));
grid on;
xlabel('frequence normalisee'); ylabel('dB');
title('Spectre de bloc^4  →  raie a f = 4*df');
xlim([-0.3, 0.3]);