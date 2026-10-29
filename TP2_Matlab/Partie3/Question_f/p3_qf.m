%% Partie 3 - Question (f) : autres scenarios de debruitage
clear; clc; close all;
rng(0);

figDir = '/MATLAB Drive/TP2/Partie3/figs_partie3';

% Signal
[y, Fs] = audioread('/MATLAB Drive/TP2/Partie3/monson.wav');
y = y(:,1);
N = 1e5;  debut = round(1*Fs);
d = y(debut+1 : debut+N);
d = d / max(abs(d));
P = 10;

%% Scenario 1 et 2 : duree des trames et estimation de r_d
sigv2 = 0.5;
v = sqrt(sigv2)*randn(N,1);
x = d + v;

durees  = [3 10 20 50 100];          % en ms
modes   = {'rdi', 'rd1', 'estime'};
nD = length(durees);  nM = length(modes);
G    = zeros(nD, nM);
saut = zeros(nD, nM);

G_glob = gain_global(d, x, v, P);
for k = 1:nD
    L = round(durees(k)*1e-3 * Fs);
    for m = 1:nM
        [G(k,m), saut(k,m)] = gain_trames(d, x, v, P, L, modes{m}, sigv2);
    end
    fprintf('trames %3d ms : r_d^i = %5.2f | r_d^1 = %5.2f | estime = %5.2f dB | global = %5.2f dB\n', ...
            durees(k), G(k,1), G(k,2), G(k,3), G_glob);
end
fprintf('saut aux frontieres / saut ailleurs (r_d^i) : %s\n', mat2str(saut(:,1)', 3));

% Figure 1 : gain en fonction de la duree des trames
figure;
plot(durees, G(:,1), 'b-o', 'LineWidth', 1.3); hold on;
plot(durees, G(:,2), 'm-s', 'LineWidth', 1.3);
plot(durees, G(:,3), 'g-^', 'LineWidth', 1.3);
yline(G_glob, 'r--', 'LineWidth', 1.5);
set(gca, 'XScale', 'log'); xticks(durees); grid on;
xlabel('duree des trames (ms)'); ylabel('gain du signal reconstruit (dB)');
title(sprintf('Gain en fonction de la duree des trames (P = %d, \\sigma_v^2 = %.1f)', P, sigv2));
legend('r_d^i (trame i)', 'r_d^1 (trame 1)', 'r_d estime depuis x', 'global', 'Location', 'best');
drawnow;
exportgraphics(gcf, fullfile(figDir, 'p3_qf_gain_duree.png'), 'Resolution', 150);

% Figure 2 : discontinuites aux frontieres
figure;
plot(durees, saut(:,1), 'b-o', 'LineWidth', 1.3); hold on;
plot(durees, saut(:,2), 'm-s', 'LineWidth', 1.3);
plot(durees, saut(:,3), 'g-^', 'LineWidth', 1.3);
yline(1, 'k:', 'LineWidth', 1.5);
set(gca, 'XScale', 'log'); xticks(durees); grid on;
xlabel('duree des trames (ms)'); ylabel('saut frontieres / saut ailleurs');
title('Discontinuites aux frontieres des trames');
legend('r_d^i (trame i)', 'r_d^1 (trame 1)', 'r_d estime depuis x', 'sans discontinuite', 'Location', 'best');
drawnow;
exportgraphics(gcf, fullfile(figDir, 'p3_qf_saut_duree.png'), 'Resolution', 150);

%% Scenario 3a : niveau de bruit
sigList = [0.05 0.1 0.5 1];
nS = length(sigList);
L20 = round(20e-3 * Fs);             % trames de 20 ms
Gs = zeros(nS, 3);  RSBav = zeros(nS, 1);
for s = 1:nS
    vs = sqrt(sigList(s))*randn(N,1);
    xs = d + vs;
    RSBav(s) = 10*log10(mean(d.^2) / sigList(s));
    Gs(s,1) = gain_global(d, xs, vs, P);
    Gs(s,2) = gain_trames(d, xs, vs, P, L20, 'rdi', sigList(s));
    Gs(s,3) = gain_trames(d, xs, vs, P, L20, 'estime', sigList(s));
    fprintf('sigv2 = %.2f (RSB_av = %5.2f dB) : global = %5.2f | trames r_d^i = %5.2f | estime = %5.2f dB\n', ...
            sigList(s), RSBav(s), Gs(s,1), Gs(s,2), Gs(s,3));
end

figure;
plot(RSBav, Gs(:,1), 'r-o', 'LineWidth', 1.3); hold on;
plot(RSBav, Gs(:,2), 'b-s', 'LineWidth', 1.3);
plot(RSBav, Gs(:,3), 'g-^', 'LineWidth', 1.3);
grid on; xlabel('RSB avant debruitage (dB)'); ylabel('gain (dB)');
title(sprintf('Gain en fonction du niveau de bruit (P = %d, trames de 20 ms)', P));
legend('global', 'trames r_d^i', 'trames r_d estime', 'Location', 'best');
drawnow;
exportgraphics(gcf, fullfile(figDir, 'p3_qf_gain_bruit.png'), 'Resolution', 150);

%% Scenario 3b : bruit blanc ou bruit colore (meme puissance)
sigv2 = 0.5;
vb = sqrt(sigv2)*randn(N,1);              % bruit blanc
vc = filter(1, [1 -0.9], randn(N,1));     % bruit colore passe-bas
vc = vc * sqrt(sigv2 / mean(vc.^2));      % meme puissance que vb
Plist = 5:5:50;
Gb = zeros(size(Plist));  Gc = zeros(size(Plist));
for m = 1:length(Plist)
    Gb(m) = gain_global(d, d + vb, vb, Plist(m));
    Gc(m) = gain_global(d, d + vc, vc, Plist(m));
end
fprintf('P = 10 : bruit blanc = %.2f dB | bruit colore = %.2f dB\n', Gb(2), Gc(2));

figure;
plot(Plist, Gb, 'b-o', 'LineWidth', 1.3); hold on;
plot(Plist, Gc, 'r-s', 'LineWidth', 1.3);
grid on; xticks(Plist);
xlabel('Ordre P'); ylabel('gain (dB)');
title(sprintf('Bruit blanc et bruit colore (\\sigma_v^2 = %.1f, filtre global)', sigv2));
legend('bruit blanc', 'bruit colore (passe-bas)', 'Location', 'best');
drawnow;
exportgraphics(gcf, fullfile(figDir, 'p3_qf_bruit_colore.png'), 'Resolution', 150);

%% Ecoute : trames de 3 ms puis de 20 ms (meme echelle)
L3 = round(3e-3 * Fs);
[~, ~, dh3]  = gain_trames(d, x, v, P, L3,  'rdi', 0.5);
[~, ~, dh20] = gain_trames(d, x, v, P, L20, 'rdi', 0.5);
sc = max(abs([dh3; dh20]));
sound(dh3/sc, Fs);   pause(length(dh3)/Fs + 0.5);
sound(dh20/sc, Fs);

%% ================= Fonctions =================
function G = gain_global(d, x, v, P)
% Gain du filtre de Wiener calcule sur tout le signal
rx = xcorr(x, P-1, 'biased'); rx = rx(P:end);
rd = xcorr(d, P-1, 'biased'); rd = rd(P:end);
h  = toeplitz(rx) \ rd;
dh = filter(h, 1, x);
Sv = 0; Se = 0;
for k = 1:length(d)
    Sv = Sv + v(k)^2;
    Se = Se + (d(k) - dh(k))^2;
end
G = 10*log10(Sv / Se);
end

function [G, saut, dh] = gain_trames(d, x, v, P, L, mode, sigv2)
% Debruitage par trames de L echantillons
% mode : 'rdi' (r_d de la trame i), 'rd1' (r_d de la trame 1),
%        'estime' (r_d = r_x - sigv2*delta, sans le signal propre)
K = floor(length(d) / L);  M = K*L;
dh  = zeros(M,1);
rd1 = xcorr(d(1:L), P-1, 'biased'); rd1 = rd1(P:end);
for i = 1:K
    idx = (i-1)*L+1 : i*L;
    rx = xcorr(x(idx), P-1, 'biased'); rx = rx(P:end);
    switch mode
        case 'rdi'
            rd = xcorr(d(idx), P-1, 'biased'); rd = rd(P:end);
        case 'rd1'
            rd = rd1;
        case 'estime'
            rd = rx;  rd(1) = rd(1) - sigv2;
    end
    h = toeplitz(rx) \ rd;
    deb = max(1, idx(1) - (P-1));            % P-1 echantillons precedents
    ye  = filter(h, 1, x(deb:idx(end)));
    dh(idx) = ye(end-L+1:end);
end
Sv = 0; Se = 0;
for k = 1:M
    Sv = Sv + v(k)^2;
    Se = Se + (d(k) - dh(k))^2;
end
G = 10*log10(Sv / Se);
df = abs(diff(dh));
saut = mean(df(L:L:end)) / mean(df);   % 1 = pas de discontinuite
end