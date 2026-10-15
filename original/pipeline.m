close all; clear; clc;

load('signal_recu.mat');

r = signal_recu;

%% 1) Estimation delta_f par blocs
N_bloc  = 1000;
N_blocs = 700;
phase   = zeros(1, N_blocs);
for i = 0:N_blocs-1
    bloc       = r(i*N_bloc+1 : (i+1)*N_bloc);
    phase(i+1) = angle(sum(bloc .^ 4));
end
phase_unwrap  = unwrap(phase);
delta_f_blocs = diff(phase_unwrap) / (4 * 2*pi * N_bloc);

%% 2) Filtre adapté
bloc_h   = r(10000:25000);
[pxx, ~] = pwelch(bloc_h, 256, 128, 512, 1, 'centered');
h        = real(fftshift(ifft(ifftshift(sqrt(pxx)))));
h_adapte = conj(fliplr(h));
y        = conv(r, h_adapte, 'same');

%% 3) Echantillonner offset = 9
symboles = y(10:10:end);
Ns       = length(symboles);

%% 4) PLL ordre 2
t_blocs        = (1:length(delta_f_blocs)) * 100;
delta_f_interp = interp1(t_blocs, delta_f_blocs, 1:Ns, 'linear', 'extrap');

alpha  = 0.1;
beta   = 0.0225;
phi    = 0;
delta  = delta_f_interp(1) * 2*pi;
s_corr = zeros(Ns, 1);

for k = 1:Ns
    s         = symboles(k) * exp(-1j * phi);
    s_corr(k) = s;
    s_norm    = s / abs(s);
    erreur    = imag(s_norm^4);
    delta     = delta + beta  * erreur;
    phi       = phi   + delta + alpha * erreur;
end

%% 5) Détection de trame
amplitude      = abs(s_corr);
seuil          = max(amplitude) * 0.2;
indices_actifs = find(amplitude > seuil);
debut_trame    = indices_actifs(1);
fin_trame      = indices_actifs(end);
signal_utile   = s_corr(debut_trame:fin_trame);

fprintf('debut_trame = %d\n', debut_trame);
fprintf('fin_trame   = %d\n', fin_trame);
fprintf('nb symboles utiles = %d\n', length(signal_utile));

%% 6) T = 128
T         = 128;
nb_pixels = T * T;

%% 7) 8 mappings possibles
figure;
mapping = 0;
for miroir = 0:1
    for rot = 0:3
        mapping = mapping + 1;

        if miroir == 0
            s_m = signal_utile;
        else
            s_m = conj(signal_utile);
        end

        s_rot = s_m .* exp(1j * rot * pi/2);

        bits = zeros(length(s_rot)*2, 1);
        for k = 1:length(s_rot)
            a = angle(s_rot(k));
            if     a >= 0    && a < pi/2,  b = [0 0];
            elseif a >= pi/2 && a <= pi,   b = [0 1];
            elseif a < -pi/2,              b = [1 1];
            else,                          b = [1 0];
            end
            bits(2*k-1:2*k) = b;
        end

        octets = reshape(bits(1:nb_pixels*8), nb_pixels, 8);
        pixels = bi2de(octets);
        matImg = reshape(pixels, T, T);

        subplot(2, 4, mapping);
        imagesc(matImg);
        colormap gray;
        axis image; 
        title(sprintf('rot=%d mir=%d', rot, miroir));
    end
end

sgtitle('T=128 - 8 mappings');