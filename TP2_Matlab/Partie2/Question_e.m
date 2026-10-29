%% Partie 2 - Question (e) : EQM_min en fonction de l'ordre P
clear; clc; close all;
if ~isfolder('figs_partie2'), mkdir('figs_partie2'); end

% Parametres
N     = 1e4;
n     = (1:N)';
w0    = 0.05*pi;
phi   = pi/4;
sigv2 = 0.1;
Plist = 5:5:50;
nP    = length(Plist);

% Signal propre et signal bruite (meme realisation pour tous les P)
d = sin(n*w0 + phi);
x = d + sqrt(sigv2)*randn(N,1);

EQMmin = zeros(1, nP);
eqmN   = zeros(1, nP);     % eqm(N) mesuree, pour comparaison

for m = 1:nP
    P = Plist(m);

    % Filtre de Wiener d'ordre P
    rx = xcorr(x, P-1, 'biased'); rx = rx(P:end);
    rd = xcorr(d, P-1, 'biased'); rd = rd(P:end);
    wopt = toeplitz(rx) \ rd;
    dhat = filter(wopt, 1, x);

    % EQM_min theorique
    EQMmin(m) = rd(1) - wopt' * rd;

    % eqm(N) mesuree (boucle for)
    S = 0;
    for k = 1:N
        S = S + (d(k) - dhat(k))^2;
    end
    eqmN(m) = S / N;

    fprintf('P = %2d : EQM_min = %.4f | eqm(N) = %.4f\n', P, EQMmin(m), eqmN(m));
end

% Figure 35 : EQM_min en fonction de P
figure;
plot(Plist, EQMmin, 'b-o', 'LineWidth', 1.3, 'MarkerFaceColor', 'b'); hold on;
plot(Plist, eqmN, 'r--s', 'LineWidth', 1.3);
grid on;
xlabel('Ordre P'); ylabel('EQM');
title('EQM_{min} en fonction de P (\sigma_v^2 = 0.1)');
legend('EQM_{min} (formule)', 'eqm(N) (mesuree)', 'Location', 'best');
xticks(Plist);
saveas(gcf, 'figs_partie2/fig35.png');