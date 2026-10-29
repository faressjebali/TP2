%% Partie 3 - Question (d) : EQM_min par trame et description statistique
clear; clc; close all;
rng(0);                              % meme bruit qu'en (c)
if ~isfolder('figs_partie3'), mkdir('figs_partie3'); end

figDir = 'figs_partie3';

% Signal (comme en (c))
[y, Fs] = audioread('monson.wav');
y = y(:,1);
N = 1e5;  debut = round(1*Fs);
d = y(debut+1 : debut+N);
d = d / max(abs(d));
P = 10;  sigv2 = 0.5;

% Decoupage en trames de 3 ms
L = round(3e-3 * Fs);
K = floor(N / L);

gains = zeros(K,1);
EQMs  = zeros(K,1);
Pds   = zeros(K,1);                  % puissance de chaque trame

% Debruitage trame par trame
for i = 1:K
    idx = (i-1)*L+1 : i*L;
    [~, ~, gains(i), EQMs(i)] = debruitage_sig(d(idx), P, sigv2);

    % Puissance de la trame (boucle for)
    S = 0;
    for k = idx
        S = S + d(k)^2;
    end
    Pds(i) = S / L;
end

% Debruitage sur tout le signal (reference)
[~, ~, gain_g, EQM_g] = debruitage_sig(d, P, sigv2);

% Description statistique
gs = sort(gains);  es = sort(EQMs);
fprintf('--- Gain (dB) ---\n');
fprintf('moyenne = %.2f | ecart-type = %.2f | mediane = %.2f\n', mean(gains), std(gains), median(gains));
fprintf('min = %.2f | max = %.2f | Q1 = %.2f | Q3 = %.2f | global = %.2f\n', ...
    gs(1), gs(end), gs(round(0.25*K)), gs(round(0.75*K)), gain_g);
fprintf('--- EQM_min ---\n');
fprintf('moyenne = %.4f | ecart-type = %.4f | mediane = %.4f\n', mean(EQMs), std(EQMs), median(EQMs));
fprintf('min = %.4f | max = %.4f | Q1 = %.4f | Q3 = %.4f | global = %.4f\n', ...
    es(1), es(end), es(round(0.25*K)), es(round(0.75*K)), EQM_g);
fprintf('Trames avec EQM_min < 0 : %d\n', sum(EQMs < 0));

% Figure 1 : EQM_min par trame
figure;
plot(1:K, EQMs, 'b.-'); hold on;
yline(EQM_g, 'r--', 'LineWidth', 1.5);
yline(mean(EQMs), 'k:', 'LineWidth', 1.5);
grid on; xlabel('numero de trame i'); ylabel('EQM_{min}');
title(sprintf('EQM_{min} par trame (P = %d, \\sigma_v^2 = %.1f)', P, sigv2));
legend('EQM_{min} par trame', 'EQM_{min} global', 'EQM_{min} moyenne', 'Location', 'best');
saveas(gcf, fullfile(figDir, 'p3_qd_eqm_trames.png'));

% Figure 2 : histogramme du gain
figure;
histogram(gains, 30); hold on;
xline(gain_g, 'r--', 'LineWidth', 1.5);
xline(mean(gains), 'k:', 'LineWidth', 1.5);
grid on; xlabel('gain (dB)'); ylabel('nombre de trames');
title('Distribution du gain par trame');
legend('trames', 'gain global', 'gain moyen', 'Location', 'best');
saveas(gcf, fullfile(figDir, 'p3_qd_hist_gain.png'));

% Figure 3 : histogramme de EQM_min
figure;
histogram(EQMs, 30); hold on;
xline(EQM_g, 'r--', 'LineWidth', 1.5);
xline(mean(EQMs), 'k:', 'LineWidth', 1.5);
grid on; xlabel('EQM_{min}'); ylabel('nombre de trames');
title('Distribution de EQM_{min} par trame');
legend('trames', 'EQM_{min} global', 'EQM_{min} moyenne', 'Location', 'best');
saveas(gcf, fullfile(figDir, 'p3_qd_hist_eqm.png'));

% Figure 4 : gain en fonction de la puissance de la trame
figure;
plot(Pds, gains, 'b.');
grid on; xlabel('puissance de la trame P_d^i'); ylabel('gain (dB)');
title('Gain par trame en fonction de la puissance de la trame');
saveas(gcf, fullfile(figDir, 'p3_qd_gain_vs_puissance.png'));