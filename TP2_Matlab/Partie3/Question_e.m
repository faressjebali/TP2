%% Partie 3 - Question (e) : filtre par trame avec r_d calcule sur la trame 1
clear; clc; close all;
rng(0);
if ~isfolder('figs_partie3'), mkdir('figs_partie3'); end

figDir = 'figs_partie3';

% Signal
[y, Fs] = audioread('monson.wav');
y = y(:,1);
N = 1e5;  debut = round(1*Fs);
d = y(debut+1 : debut+N);
d = d / max(abs(d));
P = 10;  sigv2 = 0.5;

% Decoupage en trames de 3 ms
L = round(3e-3 * Fs);
K = floor(N / L);
M = K*L;                             % echantillons couverts par les trames

% Bruit genere une seule fois : meme signal bruite pour les 3 methodes
v = sqrt(sigv2)*randn(N,1);
x = d + v;

% r_d calcule uniquement sur la trame 1
rd1 = xcorr(d(1:L), P-1, 'biased'); rd1 = rd1(P:end);

gain_i = zeros(K,1);                 % methode (c) : r_d de la trame i
gain_1 = zeros(K,1);                 % methode (e) : r_d de la trame 1
dhat_i = zeros(M,1);
dhat_1 = zeros(M,1);

for i = 1:K
    idx = (i-1)*L+1 : i*L;
    xi = x(idx);  di = d(idx);  vi = v(idx);

    % R_x de la trame i
    rx = xcorr(xi, P-1, 'biased'); rx = rx(P:end);
    Rx = toeplitz(rx);

    % r_d de la trame i
    rdi = xcorr(di, P-1, 'biased'); rdi = rdi(P:end);

    % Les deux filtres
    h_i = Rx \ rdi;
    h_1 = Rx \ rd1;
    deb = max(1, idx(1) - (P-1));
    xe  = x(deb : idx(end));                 % trame + P-1 echantillons avant
    dh_i = filter(h_i, 1, xe);  dh_i = dh_i(end-L+1:end);
    dh_1 = filter(h_1, 1, xe);  dh_1 = dh_1(end-L+1:end);

    % Reconstruction
    dhat_i(idx) = dh_i;
    dhat_1(idx) = dh_1;

    % Gains de la trame (boucle for)
    Sv = 0; Se_i = 0; Se_1 = 0;
    for k = 1:L
        Sv   = Sv   + vi(k)^2;
        Se_i = Se_i + (di(k) - dh_i(k))^2;
        Se_1 = Se_1 + (di(k) - dh_1(k))^2;
    end
    gain_i(i) = 10*log10(Sv / Se_i);
    gain_1(i) = 10*log10(Sv / Se_1);
end

% Methode globale (un seul filtre pour tout le signal)
rx = xcorr(x, P-1, 'biased'); rx = rx(P:end);
rd = xcorr(d, P-1, 'biased'); rd = rd(P:end);
h_g = toeplitz(rx) \ rd;
dhat_g = filter(h_g, 1, x);

% Gain de chaque signal reconstruit, sur les M premiers echantillons
Sv = 0; Se_g = 0; Se_i = 0; Se_1 = 0;
for k = 1:M
    Sv   = Sv   + v(k)^2;
    Se_g = Se_g + (d(k) - dhat_g(k))^2;
    Se_i = Se_i + (d(k) - dhat_i(k))^2;
    Se_1 = Se_1 + (d(k) - dhat_1(k))^2;
end
G_g = 10*log10(Sv / Se_g);
G_i = 10*log10(Sv / Se_i);
G_1 = 10*log10(Sv / Se_1);

fprintf('--- Gain du signal reconstruit (dB) ---\n');
fprintf('global = %.2f | trames r_d^i = %.2f | trames r_d^1 = %.2f\n', G_g, G_i, G_1);
fprintf('--- Gain par trame (dB) ---\n');
fprintf('r_d^i : moyenne = %.2f | mediane = %.2f | ecart-type = %.2f\n', ...
        mean(gain_i), median(gain_i), std(gain_i));
fprintf('r_d^1 : moyenne = %.2f | mediane = %.2f | ecart-type = %.2f\n', ...
        mean(gain_1), median(gain_1), std(gain_1));
fprintf('Trames ou r_d^1 fait mieux que r_d^i : %d / %d\n', sum(gain_1 > gain_i), K);

% Figure 1 : gain par trame avec r_d^1
figure;
plot(1:K, gain_1, 'm.-'); hold on;
yline(mean(gain_1), 'k:', 'LineWidth', 1.5);
yline(G_g, 'r--', 'LineWidth', 1.5);
grid on; xlabel('numero de trame i'); ylabel('gain (dB)');
title(sprintf('Gain par trame avec r_d de la trame 1 (P = %d, \\sigma_v^2 = %.1f)', P, sigv2));
legend('gain par trame (r_d^1)', 'gain moyen par trame', 'gain global', 'Location', 'best');
drawnow;
exportgraphics(gcf, fullfile(figDir, 'p3_qe_gain_trames.png'), 'Resolution', 150);

% Figure 2 : comparaison des deux methodes par trame
figure;
plot(1:K, gain_i, 'b.-'); hold on;
plot(1:K, gain_1, 'm.-');
yline(G_g, 'r--', 'LineWidth', 1.5);
grid on; xlabel('numero de trame i'); ylabel('gain (dB)');
title('Comparaison : r_d de chaque trame et r_d de la trame 1');
legend('r_d^i (question c)', 'r_d^1 (question e)', 'gain global', 'Location', 'best');
drawnow;
exportgraphics(gcf, fullfile(figDir, 'p3_qe_comparaison.png'), 'Resolution', 150);

% Figure 3 : distributions des deux methodes (contours, sans transparence)
figure;
histogram(gain_i, 30, 'DisplayStyle', 'stairs', 'EdgeColor', 'b', 'LineWidth', 1.5); hold on;
histogram(gain_1, 30, 'DisplayStyle', 'stairs', 'EdgeColor', 'm', 'LineWidth', 1.5);
xline(G_g, 'r--', 'LineWidth', 1.5);
grid on; xlabel('gain (dB)'); ylabel('nombre de trames');
title('Distribution du gain par trame');
legend('r_d^i (question c)', 'r_d^1 (question e)', 'gain global', 'Location', 'best');
drawnow;
exportgraphics(gcf, fullfile(figDir, 'p3_qe_hist.png'), 'Resolution', 150);

%% Ecoute (lancer cette section avec "Run Section")
soundsc(dhat_g, Fs);  pause(N/Fs + 0.5);
soundsc(dhat_i, Fs);  pause(N/Fs + 0.5);
soundsc(dhat_1, Fs);