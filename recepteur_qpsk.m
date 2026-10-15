%% recepteur_qpsk.m
% Décodage d'une image de 128 x 128 pixels transmise en QPSK, à partir du
% signal reçu en bande de base, sans connaître le filtre d'émission, le
% rythme symbole ni le décalage de fréquence.
%
% Entrée : data/signal_recu.mat (variable signal_recu), fourni par l'école
% et non publié dans ce dépôt (voir data/README.md).

clear; close all; clc;

load(fullfile('data', 'signal_recu.mat'));
r = signal_recu(:).';                       % 893 917 échantillons complexes

%% 1) Dérive de fréquence, estimée par blocs
% r^4 efface la modulation QPSK : il ne reste que exp(j*4*(2*pi*df*t + phi)).
N_bloc  = 1000;
N_blocs = 700;
phase = zeros(1, N_blocs);
for i = 0:N_blocs-1
    phase(i+1) = angle(sum(r(i*N_bloc+1:(i+1)*N_bloc) .^ 4));
end
delta_f_blocs = diff(unwrap(phase)) / (4 * 2*pi * N_bloc);   % cycles par échantillon

%% 2) Filtre adapté estimé à l'aveugle
% Symboles décorrélés : la DSP du signal reçu est proportionnelle à |H(f)|^2.
[pxx, ~] = pwelch(r(10000:25000), 256, 128, 512, 1, 'centered');
h        = real(fftshift(ifft(ifftshift(sqrt(pxx)))));
h_adapte = conj(fliplr(h(:).'));
y        = conv(r, h_adapte, 'same');

%% 3) Échantillonnage
% Raie à f = 0,1 dans le spectre de |r|^2 : Ts = 10 Te. Offset retenu : 9.
symboles = y(10:10:end);
Ns       = length(symboles);

%% 4) Boucle à verrouillage de phase d'ordre 2
% L'ordre 2 suit à la fois la phase et la dérive de fréquence.
t_blocs = (1:length(delta_f_blocs)) * 100;
alpha   = 0.1;
beta    = 1e-4;
phi     = 0;
delta   = interp1(t_blocs, delta_f_blocs, 1, 'linear', 'extrap') * 2*pi;
s_corr  = zeros(1, Ns);
for k = 1:Ns
    s         = symboles(k) * exp(-1j*phi);
    s_corr(k) = s;
    erreur    = imag((s/abs(s))^4);         % détecteur de phase QPSK
    delta     = delta + beta*erreur;
    phi       = phi + delta + alpha*erreur;
end

%% 5) Détection de la trame utile (seuil à 20 % de l'amplitude maximale)
actifs       = find(abs(s_corr) > 0.2*max(abs(s_corr)));
signal_utile = s_corr(actifs(1):actifs(end));
fprintf('Trame : symboles %d à %d (%d symboles)\n', actifs(1), actifs(end), numel(signal_utile));

%% 6) Décision : rotation de pi/4 puis lecture des signes
% La PLL laisse une ambiguïté de k*pi/2 et l'ordre des deux bits est inconnu :
% on garde la combinaison qui donne l'image la plus régulière.
T = 128;
meilleur = -Inf;
for permute_bits = 0:1
    for rot = 0:3
        s  = signal_utile * exp(1j*pi/4) * exp(1j*rot*pi/2);
        b1 = double(real(s) < 0);
        b2 = double(imag(s) < 0);
        if permute_bits
            [b1, b2] = deal(b2, b1);
        end
        img   = bits_vers_image([b1; b2], T);
        score = min(correlation_voisins(img, 1), correlation_voisins(img, 2));
        if score > meilleur
            meilleur     = score;
            image_finale = img;
            choix        = [rot, permute_bits];
        end
    end
end
fprintf('Rotation retenue : %d x pi/2, bits permutés : %d\n', choix(1), choix(2));

figure;
imagesc(image_finale); colormap gray; axis image off;
title('Image décodée');

%% Fonctions locales
function img = bits_vers_image(hatB, T)
    bits   = hatB(:);
    bits   = bits(1:T*T*8);
    octets = reshape(bits, [], 8);          % même convention que le squelette fourni
    pixels = octets * (2.^(0:7)).';         % équivalent de bi2de, bit de poids faible en premier
    img    = reshape(pixels, T, T);
end

function c = correlation_voisins(img, dim)
    x = double(img);
    if dim == 1
        a = x(1:end-1, :);  b = x(2:end, :);
    else
        a = x(:, 1:end-1);  b = x(:, 2:end);
    end
    c = corrcoef(a(:), b(:));
    c = c(1, 2);
end
