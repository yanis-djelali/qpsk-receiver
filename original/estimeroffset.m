close all; clear; clc;

load('signal_recu.mat');

r = signal_recu;

%% Trouver h depuis le bloc 20000-22000
bloc = r(20000:22000);
[pxx, ~] = pwelch(bloc, 256, 128, 512, 1, 'centered');
H_module  = sqrt(pxx);
H_ifft    = ifftshift(H_module);
h         = real(ifft(H_ifft));
h         = fftshift(h);

%% Filtre adapté = retourné conjugué de h
h_adapte = conj(fliplr(h));

%% Appliquer le filtre adapté sur tout le signal
y = conv(r, h_adapte, 'same');

%% Tester les 10 offsets d'échantillonnage
% On prend 2000 symboles au milieu du signal (zone active)
debut = 30000;   % echantillon de depart

figure;
for offset = 0:9
    % Echantillonner : 1 point tous les 10, avec cet offset
    idx      = debut + offset : 10 : debut + offset + 2000*10;
    symboles = y(idx);
    
    subplot(2, 5, offset+1);
    plot(real(symboles), imag(symboles), '.', 'MarkerSize', 3);
    grid on;
    axis([-0.3 0.3 -0.3 0.3]);
    title(sprintf('offset = %d', offset));
    xlabel('I'); ylabel('Q');
end

sgtitle('Constellation pour chaque offset  →  lequel a 4 clusters nets ?');