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

%% 3) Echantillonner offset = 8
symboles = y(9:10:end);
Ns       = length(symboles);

%% 4) Interpoler delta_f
t_blocs        = (1:length(delta_f_blocs)) * 100;
delta_f_interp = interp1(t_blocs, delta_f_blocs, 1:Ns, 'linear', 'extrap');

%% 5) PLL ordre 2 avec normalisation
figure;
alpha_list = [0.01, 0.05, 0.1, 0.3];

for idx = 1:4
    alpha = alpha_list(idx);
    beta  = alpha^2 / 4;
    
    phi    = 0;
    delta  = delta_f_interp(1) * 2*pi;
    s_corr = zeros(Ns, 1);
    
    for k = 1:Ns
        s         = symboles(k) * exp(-1j * phi);
        s_corr(k) = s;
        
        % Normaliser sur le cercle unité avant l'erreur
        s_norm = s / abs(s);
        erreur = imag(s_norm^4);
        
        delta = delta + beta  * erreur;
        phi   = phi   + delta + alpha * erreur;
    end
    
    subplot(2, 2, idx);
    plot(real(s_corr(1:5000)), imag(s_corr(1:5000)), '.', 'MarkerSize', 2);
    grid on; axis equal;
    title(sprintf('alpha=%.2f  beta=%.4f', alpha, beta));
    xlabel('I'); ylabel('Q');
end

sgtitle('PLL avec normalisation - 4 clusters ?');