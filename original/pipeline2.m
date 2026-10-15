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
bloc_h   = r(20000:22000);
[pxx, ~] = pwelch(bloc_h, 256, 128, 512, 1, 'centered');
h        = real(fftshift(ifft(ifftshift(sqrt(pxx)))));
h_adapte = conj(fliplr(h));
y        = conv(r, h_adapte, 'same');

%% 3) Echantillonner offset = 9
symboles = y(10:10:end);
Ns       = length(symboles);

t_blocs        = (1:length(delta_f_blocs)) * 100;
delta_f_interp = interp1(t_blocs, delta_f_blocs, 1:Ns, 'linear', 'extrap');

T = 128;

%% Test plusieurs alpha
alpha_list = [0.02, 0.05, 0.1, 0.15];

figure;
for ai = 1:4
    alpha = alpha_list(ai);
    beta  = alpha^2 / 4;
    
    phi   = 0;
    delta = delta_f_interp(1) * 2*pi;
    s_corr = zeros(Ns, 1);
    
    for k = 1:Ns
        s         = symboles(k) * exp(-1j * phi);
        s_corr(k) = s;
        s_norm    = s / abs(s);
        erreur    = imag(s_norm^4);
        delta     = delta + beta  * erreur;
        phi       = phi   + delta + alpha * erreur;
    end
    
    % Detection trame
    amplitude      = abs(s_corr);
    seuil          = max(amplitude) * 0.2;
    indices_actifs = find(amplitude > seuil);
    signal_utile   = s_corr(indices_actifs(1):indices_actifs(end));
    Ns_utile       = length(signal_utile);
    
    % rot=3 mir=1
    s_m  = conj(signal_utile) .* exp(1j * 3 * pi/2);
    hatB = zeros(2, Ns_utile);
    for k = 1:Ns_utile
        a = angle(s_m(k));
        if     a >= 0    && a < pi/2,  label = 0;
        elseif a >= pi/2 && a <= pi,   label = 1;
        elseif a < -pi/2,              label = 3;
        else,                          label = 2;
        end
        hatB(:,k) = de2bi(label, 2)';
    end
    
    hatMatBitImg = reshape(hatB(:), [], 8);
    matImg       = reshape(bi2de(hatMatBitImg), T, T);
    
    subplot(2,2,ai);
    imagesc(matImg);
    colormap gray; axis image;
    title(sprintf('alpha=%.2f  beta=%.4f', alpha, beta));
end

sgtitle('Quel alpha donne le tigre le plus net ?');