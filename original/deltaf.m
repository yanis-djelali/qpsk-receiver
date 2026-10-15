close all; clear; clc;

load('signal_recu.mat');

r = signal_recu;

figure;

for i = 0:700
    % Bloc de 1000 echantillons
    bloc = r(i*1000+1 : (i+1)*1000);
    
    plot(real(bloc), imag(bloc), '.');
    title(sprintf('i = %d  (echantillons %d a %d)', i, i*1000+1, (i+1)*1000));
    xlabel('Partie reelle');
    ylabel('Partie imaginaire');
    grid on;
    axis([-0.2 0.2 -0.2 0.2]);
    
    drawnow limitrate;
end