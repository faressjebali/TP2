clear; clc; close all;
dossier = fileparts(mfilename('fullpath'));
rep = fullfile(dossier, 'figs_partie4', filesep);
if ~isfolder(rep), mkdir(rep); end
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');

[d1, d2, x, Fs] = charger_voix(8);
N = length(x);

%% Sans recouvrement / avec recouvrement de 50 % (P = 20)
P   = 20;
Lms = [10 20 40 100];
nL  = length(Lms);
G_sr = zeros(1,nL);  G_ola = zeros(1,nL);      % gains (homme)
S_sr = zeros(1,nL);  S_ola = zeros(1,nL);      % sauts aux frontieres

for l = 1:nL
    L = round(Lms(l)*1e-3*Fs);
    dh_sr  = separation_trames(x, d1, P, L);
    dh_ola = separation_ola(x, d1, P, L);

    E = L+1 : N-L;                   % on exclut les bords (fenetres incompletes)
    G_sr(l)  = 10*log10(mean(d1(E).^2) / mean((d1(E) - dh_sr(E)).^2));
    G_ola(l) = 10*log10(mean(d1(E).^2) / mean((d1(E) - dh_ola(E)).^2));

    % Saut moyen aux debuts de trame / saut moyen partout (1 = pas de discontinuite)
    df = abs(diff(dh_sr));   S_sr(l)  = mean(df(L:L:end))         / mean(df);
    df = abs(diff(dh_ola));  S_ola(l) = mean(df(L/2:L/2:end))     / mean(df);

    fprintf('L = %3d ms | gain : %.2f dB (sans) -> %.2f dB (OLA) | saut : %.2f -> %.2f\n', ...
        Lms(l), G_sr(l), G_ola(l), S_sr(l), S_ola(l));
end

% Meilleur reglage de la section precedente (P = 50, L = 40 ms)
L = round(40e-3*Fs);  E = L+1 : N-L;
dh_sr  = separation_trames(x, d1, 50, L);
dh_ola = separation_ola(x, d1, 50, L);
fprintf('P = 50, L = 40 ms : gain %.2f dB (sans) -> %.2f dB (OLA)\n', ...
    10*log10(mean(d1(E).^2) / mean((d1(E) - dh_sr(E)).^2)), ...
    10*log10(mean(d1(E).^2) / mean((d1(E) - dh_ola(E)).^2)));

%% Figures
figure;
subplot(2,1,1);
semilogx(Lms, G_sr, 'bo-', 'LineWidth', 1.3); hold on;
semilogx(Lms, G_ola, 'rs-', 'LineWidth', 1.3);
grid on; ylabel('gain homme (dB)');
title(sprintf('Sans recouvrement / recouvrement 50 %% + overlap-add (P = %d)', P));
legend('sans recouvrement', 'overlap-add 50 %', 'Location', 'best');
subplot(2,1,2);
semilogx(Lms, S_sr, 'bo-', 'LineWidth', 1.3); hold on;
semilogx(Lms, S_ola, 'rs-', 'LineWidth', 1.3);
yline(1, 'k--');
grid on; xlabel('longueur de trame (ms)'); ylabel('saut frontieres / saut moyen');
legend('sans recouvrement', 'overlap-add 50 %', 'pas de discontinuite', 'Location', 'best');
saveas(gcf, [rep 'p4_q4_ola.png']);

%% Ecoute (L = 10 ms, P = 20) : sans recouvrement, puis overlap-add
L = round(10e-3*Fs);
sound(separation_trames(x, d1, P, L), Fs);  pause(N/Fs + 0.5);
sound(separation_ola(x, d1, P, L), Fs);
